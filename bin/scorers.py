#!/usr/bin/env python3
"""Score matrices computed in Python, alongside the R steps in call_enrichment.

Each function takes a peptide by sample count matrix and returns a score matrix
of the same shape. None of them writes files. Thresholds and provenance live
with the caller.

Sources for the published parameters:
  arcsinh z   Olin 2026, arcsinh transform then a z score against the
              beads-only samples, called at z > 3.5. The cofactor is not part
              of that paper, see COFACTOR below.
  Aitchison   the centred log ratio, clr(x) = log(x / g(x)) for g the
              geometric mean across peptides, then the same beads-only z.
              Requires positive values, hence the pseudocount.
  true GP     a generalized Poisson with log reference abundance as a linear
              predictor of the mean, scored as -log10 of the upper tail
              probability. This is the same distribution family Larman 2011
              used, not the same construction. See score_true_gp.
"""
import sys

import numpy as np
import pandas as pd

# Cofactor for the arcsinh transform. Five is the convention carried over from
# mass cytometry. Olin 2026 specifies the transform but not a cofactor, so this
# is a choice, not a published value. Changing it changes every score.
COFACTOR = 5.0

# Pseudocount added before the log ratio. The centred log ratio is undefined at
# zero and PhIP-seq count matrices are mostly zero.
PSEUDOCOUNT = 0.5

# Peptides drawn for the generalized Poisson fit. Larger libraries are
# subsampled to keep the fit tractable. A library smaller than this is sampled
# with repeats, which is what the reference implementation does.
GP_FIT_PEPTIDES = 12000


def cpm(counts):
    """Counts per million, per sample."""
    lib = counts.sum(axis=0).replace(0, np.nan)
    return counts.div(lib, axis=1) * 1e6


def _z_against(matrix, reference):
    """Z score every column against the mean and spread of the reference set.

    A peptide with no spread across the reference samples would divide by zero.
    The median non-zero spread is substituted, which keeps that peptide on the
    same scale as the rest instead of sending it to infinity.
    """
    values = matrix[reference].values
    mu = values.mean(axis=1, keepdims=True)
    sd = values.std(axis=1, keepdims=True)
    usable = sd[sd > 1e-9]
    fill = np.median(usable) if usable.size else 1.0
    sd = np.where(sd < 1e-9, fill, sd)
    return pd.DataFrame((matrix.values - mu) / sd,
                        index=matrix.index, columns=matrix.columns)


def score_arcsinh(counts, beads, cofactor=COFACTOR):
    """Olin 2026. Z score of arcsinh transformed CPM against the beads-only set.

    arcsinh is variance stabilising at low counts, where a log transform
    diverges at zero, and preserves the ordering of the high tail.
    """
    transformed = np.arcsinh(cpm(counts).fillna(0) / cofactor)
    return _z_against(transformed, beads)


def score_aitchison(counts, beads, pseudocount=PSEUDOCOUNT):
    """Centred log ratio, then a z score against the beads-only set.

    Subtracting the per sample mean of the log proportions is division by the
    geometric mean, which lifts the data off the simplex so that a rise in one
    peptide no longer forces a fall in the others.
    """
    shifted = counts + pseudocount
    log_p = np.log(shifted / shifted.sum(axis=0))
    clr = log_p - log_p.mean(axis=0)
    return _z_against(clr, beads)


def score_true_gp(counts, reference, n_fit=GP_FIT_PEPTIDES):
    """Generalized Poisson null built from a reference set.

    The mean is modelled as a linear function of log abundance in the reference
    set, with a library size offset, and the score is -log10 of the probability
    of seeing at least the observed count under that null.

    This is not the Larman 2011 construction, which fits a separate generalized
    Poisson at each input abundance level and then regresses both parameters
    against that level, so its dispersion varies with abundance. Here the fit is
    global and the dispersion is one number. Xu 2015 adds zero inflation on top
    of the Larman model, which is also absent. What the two share is the
    distribution family and the use of abundance to set the expected count.

    The reference set supplies both the abundance and the counts the fit is run
    on, so it measures how much a peptide varies between replicates of that set.
    A negative dispersion means the reference varies less than a Poisson would,
    which is a sign the null is too tight for the samples being scored against
    it rather than a property worth using.
    """
    try:
        import statsmodels.api as sm
        from statsmodels.discrete.discrete_model import GeneralizedPoisson
        from statsmodels.distributions.discrete import genpoisson_p
    except ImportError as exc:
        raise SystemExit(
            "the true_gp scorer needs statsmodels, which the other scorers do "
            f"not: {exc}") from exc

    counts = counts.astype(float)
    lib = counts.sum(axis=0)
    mean_lib = float(lib[reference].mean())
    abundance = counts[reference].mean(axis=1)
    log_abundance = np.log1p(abundance.values)

    index = np.linspace(0, len(abundance) - 1, n_fit).astype(int)
    response, predictor, offset = [], [], []
    for column in reference:
        response.append(counts[column].values[index])
        predictor.append(log_abundance[index])
        offset.append(np.full(index.size, np.log(lib[column] / mean_lib)))
    fit = GeneralizedPoisson(np.concatenate(response),
                             sm.add_constant(np.concatenate(predictor)),
                             offset=np.concatenate(offset), p=1
                             ).fit(disp=0, maxiter=100)
    b0, b1, alpha = fit.params[0], fit.params[1], fit.params[-1]
    print(f"[score] generalized Poisson fit: intercept {b0:.3f} slope {b1:.3f} "
          f"dispersion {alpha:.3f}")
    if alpha < 0:
        print("[score] warning: negative dispersion means the reference set "
              "varies less than a Poisson. The null is tighter than the samples "
              "scored against it, so expect the call count to be inflated",
              file=sys.stderr)

    scores = pd.DataFrame(index=counts.index, columns=counts.columns, dtype=float)
    for column in counts.columns:
        mu = np.exp(b0 + b1 * log_abundance + np.log(lib[column] / mean_lib))
        tail = genpoisson_p.sf(np.maximum(counts[column].values - 1, -1),
                               mu, alpha, 1)
        scores[column] = -np.log10(np.clip(tail, 1e-300, 1.0))
    return scores


def apply_replicate_rule(scores, groups):
    """Xu 2015. Replace each replicate score with the minimum across the group.

    Xu called a peptide enriched only when it cleared the threshold in both
    replicates. Taking the minimum makes that a single threshold check, so the
    rule composes with any score where a higher value means more enriched.

    groups maps a group label to the sample columns in it. Groups with fewer
    than two columns present are left alone.
    """
    out = scores.copy()
    for columns in groups.values():
        present = [c for c in columns if c in scores.columns]
        if len(present) < 2:
            continue
        lowest = scores[present].min(axis=1)
        for column in present:
            out[column] = lowest
    return out


# Peptides per bin for the binned generalized Poisson fits. Moment estimates of
# the dispersion are unstable on small samples, so the bin count is derived from
# the library size rather than fixed. phip-stat uses 300 bins on libraries of
# order 100,000 peptides, which is the upper bound here.
GP_PEPTIDES_PER_BIN = 50
GP_MAX_BINS = 300


def _gp_moments(values):
    """Consul-Jain theta and lambda by method of moments.

    For the generalized Poisson the mean is theta / (1 - lambda) and the
    variance is that mean divided by (1 - lambda) squared, so lambda follows
    from the ratio of the two. Underdispersed bins are pinned at lambda zero,
    which reduces the bin to an ordinary Poisson rather than letting a negative
    lambda tighten the null.
    """
    mean = float(np.mean(values))
    var = float(np.var(values, ddof=1)) if len(values) > 1 else 0.0
    if mean <= 0 or var <= 0:
        return None
    lam = 1.0 - np.sqrt(mean / var) if var > mean else 0.0
    lam = float(np.clip(lam, 0.0, 0.99))
    return mean * (1.0 - lam), lam


def _bin_edges(abundance, n_bins):
    """Indices of peptides split into bins of equal size, ordered by abundance."""
    order = np.argsort(np.asarray(abundance), kind="stable")
    return [b for b in np.array_split(order, n_bins) if len(b) > 1]


def _fit_against_abundance(x, y):
    """Straight line through the per bin parameter estimates."""
    if len(x) < 2:
        return 0.0, float(np.mean(y)) if len(y) else 0.0
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept)


def score_binned_gp(counts, reference, zero_inflated=False, n_bins=None):
    """Larman 2011, and with zero_inflated the Xu 2015 variant.

    Peptides are grouped by their abundance in the reference set. Within each
    group a generalized Poisson is fitted to the observed counts of one sample,
    and the two fitted parameters are then regressed against the group's
    abundance to give a null that varies smoothly with it. A peptide's score is
    -log10 of the probability of seeing at least its count under that null.

    This is the construction the published method describes, and it is what
    separates it from score_true_gp: dispersion is a function of abundance
    rather than one number for the whole matrix, and the fit runs on the counts
    being scored rather than on the reference set's own counts.

    With zero_inflated the share of zeros a bin carries beyond what the fitted
    generalized Poisson predicts is treated as a separate zero generating
    process, and the tail probability is scaled by the complement of it. Xu 2015
    added that for libraries where most peptides are never observed.
    """
    try:
        from statsmodels.distributions.discrete import genpoisson_p
    except ImportError as exc:
        raise SystemExit(
            f"the binned generalized Poisson scorers need statsmodels: {exc}") from exc

    counts = counts.astype(float)
    abundance = counts[reference].mean(axis=1).values
    log_abundance = np.log1p(abundance)
    if n_bins is None:
        n_bins = int(np.clip(len(abundance) // GP_PEPTIDES_PER_BIN, 4, GP_MAX_BINS))
    bins = _bin_edges(abundance, n_bins)
    print(f"[score] {len(bins)} abundance bins, "
          f"about {len(abundance) // max(len(bins), 1)} peptides each")

    scores = pd.DataFrame(index=counts.index, columns=counts.columns, dtype=float)
    skipped = 0
    for column in counts.columns:
        observed = counts[column].values
        bin_x, bin_lam, bin_log_theta, bin_zero = [], [], [], []
        for members in bins:
            fitted = _gp_moments(observed[members])
            if fitted is None:
                continue
            theta, lam = fitted
            bin_x.append(float(np.mean(log_abundance[members])))
            bin_lam.append(lam)
            bin_log_theta.append(np.log(max(theta, 1e-9)))
            if zero_inflated:
                seen = float(np.mean(observed[members] == 0))
                predicted = float(np.exp(-theta))
                excess = (seen - predicted) / (1.0 - predicted) if predicted < 1 else 0.0
                bin_zero.append(float(np.clip(excess, 0.0, 0.99)))
        if len(bin_x) < 2:
            skipped += 1
            continue

        lam_slope, lam_intercept = _fit_against_abundance(bin_x, bin_lam)
        theta_slope, theta_intercept = _fit_against_abundance(bin_x, bin_log_theta)
        lam = np.clip(lam_intercept + lam_slope * log_abundance, 0.0, 0.99)
        theta = np.exp(theta_intercept + theta_slope * log_abundance)

        mu = theta / (1.0 - lam)
        alpha = lam / (1.0 - lam)
        tail = genpoisson_p.sf(np.maximum(observed - 1, -1), mu, alpha, 1)
        if zero_inflated:
            zero_slope, zero_intercept = _fit_against_abundance(bin_x, bin_zero)
            inflation = np.clip(zero_intercept + zero_slope * log_abundance, 0.0, 0.99)
            tail = np.where(observed > 0, (1.0 - inflation) * tail, 1.0)
        scores[column] = -np.log10(np.clip(tail, 1e-300, 1.0))

    if skipped:
        print(f"[score] {skipped} sample(s) had too few usable bins to fit",
              file=sys.stderr)
    return scores
