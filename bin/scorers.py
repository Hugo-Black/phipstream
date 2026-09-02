#!/usr/bin/env python3
"""Score matrices computed in Python, alongside the R steps in call_enrichment.

Each function takes a peptide by sample count matrix and returns a score matrix
of the same shape. None of them writes files. Thresholds and provenance live
with the caller.

Both entries here implement a published construction rather than an
approximation of one:

  score_binned_gp        Larman 2011, and with zero_inflated the Xu 2015
                         variant. Its docstring gives the construction.
  apply_replicate_rule   the Xu 2015 reproducibility criterion, which composes
                         with any score where higher means more enriched.
"""
import sys

import numpy as np
import pandas as pd



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

    This is the construction the published method describes. Two properties are
    load bearing and easy to lose in a simpler version: dispersion is a function
    of abundance rather than one number for the whole matrix, and the fit runs
    on the counts being scored rather than on the reference set's own counts.

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
