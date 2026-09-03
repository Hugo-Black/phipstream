<p align="center">
  <img src="phipstream_logo_clean.svg" alt="phipstream" width="640">
</p>

# phipstream

phipstream is a containerized workflow for PhIP-Seq read processing,
oligonucleotide counting, enrichment analysis, and replicate-level peptide
prioritisation. Given sequencing reads and a peptide library table, phipstream
produces per-sample count matrices, enrichment calls, concordance summaries, and
prioritized peptide outputs. External bioinformatics tools are executed through
pinned container images to improve portability and reproducibility across
workstations and high-performance computing environments.

Two execution modes are provided:

1. A direct command line interface composed of independent stages for sample
   sheet generation, adapter trimming, alignment, enrichment scoring, and
   prioritisation.
2. A Nextflow execution mode, based on the bundled phip-flow fork, for managed
   workflow execution and PBS submission.

## Overview

| Command | Function |
|---|---|
| `phipstream qc` | Generate an nf-core/seqinspector sample sheet from the project sample table. |
| `phipstream trim` | Remove amplicon adapter sequences with cutadapt. |
| `phipstream align` | Align R1 reads to library oligos with Bowtie 2 and produce peptide counts. |
| `phipstream score` | Fit the beads-only background with edgeR and call enrichment with BEER. |
| `phipstream prioritise` | Summarise replicate support and produce peptide shortlists. |
| `phipstream images` | List or retrieve the pinned container images. |
| `phipstream nextflow` | Execute the bundled phip-flow based workflow. |
| `phipstream submit` | Submit the bundled workflow as a PBS job. |

The direct command line interface is intended for analyses requiring explicit
stage control and inspectable intermediate files. The Nextflow execution mode is
intended for managed cluster execution and phip-flow compatible deployments.

## Installation and requirements

Required host software:

* Python 3.9 or later.
* pandas 2.0 or later, installed through `requirements.txt`.
* Docker or Apptainer. Docker is appropriate for local workstations. Apptainer is
  appropriate for clusters where Docker is not available to users.

Optional host software:

* Nextflow, required only for `phipstream nextflow` and `phipstream submit`.

Installation:

```bash
git clone https://github.com/Hugo-Black/phipstream.git
cd phipstream
python3 -m pip install -r requirements.txt
make check
```

`make check` reports the selected container engine, Python interpreter, and
optional workflow manager. If both Docker and Apptainer are available, the
container engine is detected automatically. Set `PHIPSTREAM_ENGINE=docker` or
`PHIPSTREAM_ENGINE=apptainer` to select a specific runtime.

## Validation

Retrieve the required images and execute the synthetic test dataset:

```bash
make pull
make test
```

`make pull` retrieves the container images used by the direct command line
interface. `make test` creates a temporary synthetic dataset and executes sample
sheet generation, trimming, alignment, and prioritisation. To include the BEER
enrichment step, run:

```bash
BEER=1 make test
```

Test files are created in a temporary directory and removed after completion.

## Standard analysis workflow

Each stage is an independent command and supports `--help`. Stages should be
executed in the order shown unless an analysis requires only a subset of the
pipeline.

### 1. Generate a read QC sample sheet

```bash
bin/phipstream qc \
    --sample-table  metadata/sample_table.csv \
    --output        results/qc/samplesheet.csv
```

This command writes an nf-core/seqinspector compatible sample sheet.
seqinspector should be executed separately with an appropriate Nextflow
installation.

### 2. Trim adapter sequences

```bash
bin/phipstream trim \
    --sample-table          metadata/sample_table.csv \
    --adapters              metadata/adapters.csv \
    --out-dir               results/trimmed \
    --trimmed-sample-table  metadata/sample_table_trimmed.csv \
    --minimum-length        50
```

The trimming stage writes FASTQ files to `--out-dir` and creates a new sample
table in which read path columns reference the trimmed files. The original sample
table is not modified.

The value of `--minimum-length` should reflect the shortest insert length that is
appropriate for the library and alignment configuration. After trimming, review
the cutadapt JSON and log files in `results/trimmed/cutadapt_logs`. A high
proportion of reads discarded for length can indicate an adapter specification
mismatch or an unsuitable length threshold.

### 3. Align reads and count library oligos

```bash
bin/phipstream align \
    --sample-table   metadata/sample_table_trimmed.csv \
    --peptide-table  metadata/peptide_table.csv \
    --out-dir        results/alignment \
    --preset         local \
    --mode           se_r1
```

Available modes:

* `se_r1`: single-end on R1. Library inserts are represented in a single
  orientation, so Bowtie 2 is restricted to the forward strand with `--norc`.
  This is the default.
* `se_r2`: single-end on R2, restricted to the reverse strand with `--nofw`,
  because R2 reads the reverse complement of the sense insert. Useful when R2
  carries the better base quality for a given run.
* `pe`: paired-end. Each fragment is counted once, from the primary first mate
  of a proper pair. Counting aligned segments would score both mates of a
  fragment, and both mates land on the same oligo.

`se_r2` and `pe` require a populated `fastq_r2_filepath` column and stop with a
named sample if one is missing.

Available presets:

* `local`: uses `--very-sensitive-local`. This is the default and permits soft
  clipping at read ends.
* `end-to-end`: uses `--very-sensitive`. This preset requires full-read
  alignment to an oligo.

Principal outputs:

* `results/alignment/counts/<mode>.<preset>.csv`: peptide by sample count
  matrix. The mode is part of the filename so that two modes written to one
  output directory do not overwrite each other.
* `results/alignment/align_summary.csv`: per-sample read totals, alignment
  rates, and proper pair counts for paired-end runs.

### 4. Call enrichment

```bash
bin/phipstream score \
    --counts         results/alignment/counts/se_r1.local.csv \
    --sample-table   metadata/sample_table_trimmed.csv \
    --peptide-table  metadata/peptide_table.csv \
    --out-dir        results/enrichment \
    --method         beer \
    --min-lib-size   500 \
    --posterior-threshold 0.5
```

The scoring stage prepares the wide count, sample annotation, and peptide
annotation tables consumed by the R scripts. edgeR fits the beads-only
background, and BEER calls enrichment using the resulting PhIPData object.

`--method` selects the scorer. Two run in R, three run in Python.

| Method | Scale | Default cutoff | Writes |
|---|---|---|---|
| `beer` | posterior probability | 0.5 | `beer_posterior.csv.gz`, `beer_hits.csv.gz` |
| `edger` | BH adjusted p value | 0.05 | `edger_logpval.csv.gz`, `edger_hits.csv.gz` |
| `larman_gp` | -log10 p | 2.3 | `larman_gp_mlxp.csv.gz`, `larman_gp_hits.csv.gz` |
| `xu_zigp` | -log10 p | 2.3 | `xu_zigp_mlxp.csv.gz`, `xu_zigp_hits.csv.gz` |

`beer` and `edger` come from the bundled phip-flow R steps. The other two are
computed in this repository and need statsmodels.

Every scorer here is either the published software itself or the published
construction. A method that approximates a paper without reproducing it is worse
than an absent one, because the citation implies a comparability that does not
hold, so nothing of that kind is offered.

#### The generalized Poisson scorers

`larman_gp` is the Larman 2011 construction. Peptides are grouped into bins by
their abundance in a reference set, a generalized Poisson is fitted within each
bin to the counts of the sample being scored, and both fitted parameters are
then regressed against the bin abundance. Two properties follow, and both matter:
the null varies with abundance rather than being one number for the whole
matrix, and the fit runs on the counts under test rather than on the reference
set's own counts. `--gp-bins` sets the bin count, defaulting to one per 50
peptides so the moment estimates have enough observations behind them, capped at
300.

`xu_zigp` adds the Xu 2015 zero inflation. The share of zeros a bin carries
beyond what its fitted generalized Poisson predicts is treated as a separate
zero generating process, and the tail probability is scaled by the complement of
it. That matters for libraries where most peptides are never observed in most
samples. Xu also required a peptide to clear the threshold in both replicates,
which is `--replicate-rule` rather than part of the scorer.

`--gp-null` chooses the reference set. `beads` uses the mock
immunoprecipitations, `input` uses the library reference so expected counts
follow input abundance as the published method describes. On a diluted input
reference the second is the more faithful choice but the less well behaved one,
because such a reference varies far less than the immunoprecipitated samples
being scored against it.

`--replicate-rule` is separate from the choice of scorer. It gives every
replicate in a group the lowest score in that group, so a single threshold check
means every replicate cleared it. That is the reproducibility criterion Xu 2015
applies, and it composes with any of the three Python scorers.

The difference in cost is large. On a run of 96 samples against 1,194 peptides,
edgeR finished in 13 seconds where BEER took 187 minutes, because the BEER
sampler runs serially. That makes `edger` a practical way to confirm a
deployment end to end before committing to a full scoring run.

The two methods do not agree, and are not meant to. They call at different
operating points, so a method chosen for speed is a check on the plumbing rather
than a substitute for the calls BEER makes.

A peptide is called enriched in a replicate when its posterior probability
exceeds `--posterior-threshold`. BEER super-enriched records are retained as hits
when no posterior value is returned.

Samples with fewer than `--min-lib-size` library-mapped reads are excluded before
model fitting and restored as empty output columns after scoring. This preserves
the complete sample layout while preventing shallow samples from destabilising
model fitting. Replicates that are not fitted by the model retain missing
posterior values rather than zero-valued posterior columns.

BEER requires at least two beads-only controls. Four to eight beads-only controls
are recommended. The command issues a warning below four controls and terminates
below two controls.

### 5. Summarise replicate support and prioritise peptides

```bash
bin/phipstream prioritise \
    --posterior      results/enrichment/beer_posterior.csv.gz \
    --hits           results/enrichment/beer_hits.csv.gz \
    --sample-table   metadata/sample_table_trimmed.csv \
    --peptide-table  metadata/peptide_table.csv \
    --out-dir        results/prioritised \
    --role           serum
```

Concordance is calculated within each replicate group with at least
`--min-replicates` fitted replicates. The denominator is the set of peptides
enriched in at least one replicate in the group. The numerator is the subset of
those peptides enriched in at least the requested number of replicates. Groups
without enriched peptides are omitted from the concordance table.

The prioritisation stage begins with peptides enriched in at least
`--min-replicates` replicates of one group. A peptide is retained if at least one
of the following criteria is met:

* Recurrence across at least `--min-groups` replicate groups.
* Enrichment in every fitted replicate of a group.
* Enrichment of a neighbouring tile with at least half overlap in the same group.

These criteria do not use external labels, allowing independent reference data to
be reserved for validation.

Outputs:

* `replicate_concordance.csv`
* `replicate_supported.csv`
* `shortlist.csv`
* `adjacent_pairs.csv`
* `prioritisation_summary.json`

## Input files

### Sample table

The sample table must contain one row per sequenced sample.

Required columns:

| Column | Description |
|---|---|
| `sample_name` | Sample name used as a count matrix column. |
| `fastq_filepath` | Path to the R1 FASTQ file. |
| `technical_replicate_id` | Integer replicate identifier used as the scoring key. |

Optional columns:

| Column | Description |
|---|---|
| `fastq_r2_filepath` | Path to the R2 FASTQ file. Trimmed alongside R1, and required by the `se_r2` and `pe` alignment modes. |
| `control_status` | One of `empirical`, `beads_only`, `library`, or `undetermined`. |
| `sample_role` | Sample class, for example `serum` or `beads`. |
| `participant_ID` | Replicate group used during prioritisation. |
| `in_background_null` | Boolean field used to derive `control_status` when absent. |
| `is_input_reference` | Boolean field used to derive `control_status` when absent. |

`control_status` determines how samples enter the enrichment model. Beads-only
samples define the background distribution, library samples are input references,
empirical samples are scored, and undetermined samples are removed before
scoring. If `control_status` is absent, phipstream derives it from the optional
role and boolean columns and reports the resulting status counts.

### Peptide table

Required columns:

| Column | Description |
|---|---|
| `peptide_id` | Integer peptide identifier. |
| `oligo` | Nucleotide sequence aligned against the reads. |

The prioritisation stage can also use organism, antigen, and tile start columns.
These column names are configurable through command line options.

### Adapter specification

Adapter trimming uses a CSV file with the following columns:

| Column | Description |
|---|---|
| `name` | Label used in cutadapt reports. |
| `read` | `R1` or `R2`. |
| `end` | `5` or `3`. |
| `sequence` | Adapter sequence. |
| `anchored` | `true` when the adapter must occur at the read end. |

Multiple rows may target the same read and end, permitting one trimming command
to process libraries with more than one leader or tail design. The command uses
`--match-read-wildcards`, so N bases in reads can match without consuming the
error allowance. The adapter file should contain canonical adapter sequences
rather than N-masked sequences. An example is provided in
`configs/adapters.example.csv`.

## Containers and pinned software

The direct command line interface invokes the following pinned images:

| Tool | Version | Container image | Used by |
|---|---:|---|---|
| FastQC | 0.12.1 | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` | `stages` |
| cutadapt | 5.2 | `quay.io/biocontainers/cutadapt:5.2--py312hfabe715_2` | `trim` |
| Bowtie 2 | 2.5.4 | `quay.io/biocontainers/bowtie2:2.5.4--he96a11b_7` | `align` |
| SAMtools | 1.19.2 | `quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1` | `align` |
| BEER | 1.2.0 | `quay.io/biocontainers/bioconductor-beer:1.2.0--r42hdfd78af_0` | `score` |
| edgeR | bundled in image | BEER image | `score` |

Container images can be listed or retrieved with:

```bash
bin/phipstream images
bin/phipstream images pull
```

Docker requires approximately 1.8 GB for image layers. Apptainer requires
approximately 650 MB for SIF files. The Apptainer cache path is selected from
`NXF_APPTAINER_CACHEDIR`, then `APPTAINER_CACHEDIR`, then
`~/.apptainer/cache`. On clusters, the cache should be populated from a login
node before job submission when compute nodes lack outbound network access.

For reproducible analyses, retain the pinned image tags and avoid floating tags
such as `latest`.

## Reproducibility

BEER uses the current date as a seed unless a value is supplied. phipstream pins
`--seed` by default and records it in `enrichment_summary.json` along with input
paths, thresholds, and image tags. The Nextflow execution mode uses the same
default through `params.beer_seed`.

With fixed inputs, fixed images, and a fixed seed, trimming and alignment outputs
are deterministic, and BEER returns identical posterior and hit matrices.

## Running every stage from one file

The stages above are separate programs, which suits interactive work but not a
batch scheduler, where a job runs one command. `phipstream stages` reads an INI
file holding the parameters for all of them and runs them in order.

```bash
cp configs/_template.stages.conf configs/mydataset.stages.conf
bin/phipstream stages configs/mydataset.stages.conf --dry-run
bin/phipstream stages configs/mydataset.stages.conf
```

Keys match the stage options with dashes replaced by underscores. Omit the
`adapters` key when the sample table already points at trimmed reads, and the
trimming stage is skipped.

| Key | Default | Purpose |
|---|---|---|
| `sample_table`, `peptide_table`, `out_dir` | required | inputs and destination |
| `adapters` | none | omit to skip trimming |
| `fastqc` | none | set to `true` to write FastQC archives before trimming |
| `mode`, `preset` | `se_r1`, `local` | alignment settings |
| `jobs`, `threads` | `4`, `2` | samples at once, and threads per sample |
| `method` | `beer` | `beer` or `edger`, see the scoring stage above |
| `min_lib_size`, `posterior_threshold`, `edger_fdr`, `seed` | stage defaults | scoring settings |
| `role_column`, `role`, `group_column` | `sample_role`, `serum`, `participant_ID` | which samples are prioritised, and how they group |
| `antigen_columns`, `position_column`, `tile_step` | `virus,antigen`, `tile_start`, `28` | how an adjacent tile is recognised |
| `min_replicates`, `min_groups`, `all_replicates` | `2`, `2`, `3` | prioritisation criteria |
| `walltime`, `ppn`, `mem`, `modules` | `24:00:00`, `8`, `32gb`, `tools apptainer` | batch resources, used with `--submit` |
| `python` | the running interpreter | interpreter with pandas on the cluster |
| `scratch` | none | directory for intermediate BAM files |

Setting `fastqc = true` adds a FastQC stage ahead of trimming, writing one
archive per read file into `<out_dir>/fastqc`. It runs on untrimmed reads,
because trimming removes the leader that overrepresented sequence reporting is
most likely to surface. The workflow's contamination diagnostic reads those
archives rather than the reads themselves, so this stage has to run before it.

`--resume` skips any stage whose output is already present, so a rerun continues
from the stage that failed rather than repeating enrichment calling.

`--submit` writes a PBS job script and submits it when `qsub` is on PATH.

```bash
bin/phipstream stages configs/mydataset.stages.conf \
    --submit --computerome_project <project>
```

Batch resources come from the `walltime`, `ppn`, `mem` and `modules` keys, and
the interpreter from the `python` key. See [docs/computerome.md](docs/computerome.md).

### Driving both routes from the same file

An optional `[workflow]` section holds parameters for the bundled workflow,
which covers the scoring layers the stages do not and adds the library QC and
contamination modules. `--workflow` renders that section into a Nextflow config
and runs it.

```bash
bin/phipstream stages configs/mydataset.stages.conf --workflow --dry-run
bin/phipstream stages configs/mydataset.stages.conf --workflow
bin/phipstream stages configs/mydataset.stages.conf --workflow \
    --submit --computerome_project <project>
```

`sample_table` and `peptide_table` are taken from the `[phipstream]` section, so
a dataset cannot describe itself twice and then disagree with itself. An
`adapters` CSV there is converted into the comma separated lists the workflow
expects, and rows for three prime adapters are reported and dropped because the
workflow has no parameter for them. With `fastqc = true`, `fastqc_dir` is
pointed at the stage output automatically. Anything set in `[workflow]` wins,
and `results` defaults to a `workflow` directory under `out_dir`.

Keys in `[workflow]` keep their case, because they become workflow parameter
names. `run_BEER` is not `run_beer`, and a misspelled parameter is accepted in
silence and then ignored.

The `modules`, `walltime`, `ppn` and `mem` keys are passed to the submitted job
for this route as well, so both routes are described once. `modules` must name
nextflow when `--workflow` is used with `--submit`, and exact versions are
needed where the site publishes no default, since `module load apptainer` on
such a site fails with `Unable to locate a modulefile for apptainer/<nodefault>`
and the job dies on its first line.

The job script, the scheduler output and the Nextflow report and trace are all
written under `<out_dir>/logs`, and Nextflow's scratch space under
`<out_dir>/work`. Scratch holds every intermediate file the run produces and is
the largest thing it writes, so it belongs on the filesystem that holds the
dataset rather than the one holding the code. A standalone
`phipstream submit` has no `out_dir` to work from and reads the config's own
`results` path instead, writing to `<results>/logs`. That path is only used when
it is absolute, since a relative one may hold a Nextflow variable this script
cannot resolve, and in that case output falls back to `logs/` in the repository.
`--log-dir` and `--work-dir` override all of it. A config whose `results` path
is relative gets neither, and Nextflow falls back to `work/` in the launch
directory.

## Workflow manager execution

The bundled workflow under `phip-flow/` supports alignment, count collection,
scoring, and optional modules from the fork. It does not execute the direct
command line prioritisation stage.

Execute the workflow with Nextflow:

```bash
bin/phipstream nextflow configs/dataset.config
```

Submit the workflow to a PBS cluster:

```bash
bin/phipstream submit configs/dataset.config --computerome_project <project>
```

### Scoring steps in the workflow

The workflow provides five scoring steps. Each is an independent switch, so
several can run in one pass, and each adds its own layer to the output dataset.

| Switch | Layer written | Notes |
|---|---|---|
| `run_cpm_enr_workflow` | `enrichment` | fold enrichment over the library control samples |
| `run_zscore_fit_predict` | `zscore` | binned z score against the beads-only samples, threshold `zscore_threshold`, bin floor `zscore_min_peptides_per_bin` |
| `run_gamma_poisson` | `gamma_poisson_mlxp` | gamma-Poisson over size factors, tuned by the `gamma_poisson_*` settings |
| `run_edgeR` | `edgeR_*` | edgeR, and the prior BEER consumes, cutoff `edgeR_threshold` |
| `run_BEER` | `beer_*` | BEER posterior, tuned by the `beer_*` settings |

BEER consumes the prior edgeR fits, so leave `run_edgeR` on whenever `run_BEER`
is on. All five are listed with their settings in `configs/_template.config`.

`run_zscore_fit_predict` measures each peptide against others of similar
abundance rather than against the library as a whole. Peptides are ordered by
their summed abundance across the beads-only samples and merged into bins until
each holds at least `zscore_min_peptides_per_bin`, so the number of bins is
about the library size divided by that floor. The default of 300 comes from
phippery and suits a library of order 100,000 peptides, where it yields
hundreds of bins. A library of one or two thousand peptides reaches the floor
in three or four, at which point a peptide is being compared against others of
quite different abundance and the binning is doing little. Lower the floor on a
small library, keeping in mind that the model trims to the middle 90 percent of
a bin before estimating, so a floor of 100 fits on about 90 peptides.

Each switch writes its layer into a single dataset object rather than into a
file of its own. Three further switches control how that object is written.

| Switch | Default | Output |
|---|---|---|
| `output_wide_csv` | `true` | one table per layer under `results/wide_data/` |
| `output_tall_csv` | `false` | one long-format table under `results/tall_data/` |
| `output_pickle_xarray` | `true` | the binary dataset under `results/pickle_data/`, read with phippery |

Obtaining a z score or CPM matrix as a plain table therefore takes two switches,
the scoring one and `output_wide_csv`.

### Notes on the bundled fork

New dataset configurations should be based on `configs/_template.config`. The
upstream phip-flow image contains bowtie 1, SAMtools, and phippery. phipstream
routes processes requiring Bowtie 2, cutadapt, edgeR, or BEER to the pinned
images listed above. Bowtie 1 remains selectable in the Nextflow execution mode
by setting `aligner = 'bowtie1'`.

Cluster-specific PBS and Apptainer guidance is provided in
`docs/computerome.md`.

## Command options reference

Every option accepted by each command, with its default. The tables are
generated from the programs themselves, so `phipstream <command> --help` always
agrees with what is written here.

### `phipstream qc`

| Option | Default | Purpose |
|---|---|---|
| `--sample-table` | required | table with fastq_filepath and a name column |
| `--output` | required | sample sheet to write |
| `--reads-prefix` | none | prefix joined to relative fastq paths |
| `--tags-column` | `control_status` | column used as the grouping tag |
| `--name-column` | `sample_name` | column used as the sample name |

### `phipstream trim`

| Option | Default | Purpose |
|---|---|---|
| `--sample-table` | required | table with fastq_filepath and optional fastq_r2_filepath |
| `--adapters` | required | adapter spec CSV |
| `--out-dir` | required | destination for trimmed reads |
| `--trimmed-sample-table` | required | table written with fastq paths repointed at the trimmed reads |
| `--minimum-length` | `50` | drop a read once trimming leaves it shorter than this. It is the shortest insert still worth aligning, not a fraction of the... |
| `--jobs` | `4` | samples trimmed at once |
| `--threads` | `2` | cutadapt threads per sample |

### `phipstream align`

| Option | Default | Purpose |
|---|---|---|
| `--sample-table` | required | trimmed sample table with sample_name and fastq_filepath |
| `--peptide-table` | required | peptide table with peptide_id and oligo columns |
| `--out-dir` | required | destination for counts, the reference and the summary |
| `--preset` | `local` | local allows clipped read ends, end-to-end requires the whole read to align. One of `end-to-end`, `local` |
| `--mode` | `se_r1` | reads to align. se_r1 uses R1 on the forward strand, se_r2 uses R2 on the reverse strand, pe uses both and counts each.... One of `pe`, `se_r1`, `se_r2` |
| `--jobs` | `4` | samples aligned at once |
| `--threads` | `2` | threads per sample |
| `--scratch` | none | directory for intermediate BAM files |

### `phipstream score`

| Option | Default | Purpose |
|---|---|---|
| `--counts` | required | peptide by sample counts CSV |
| `--sample-table` | required | sample table carrying control_status or the columns it is derived from |
| `--peptide-table` | required | peptide table with a peptide_id column |
| `--out-dir` | required | destination for the score and hit matrices |
| `--method` | `beer` | scoring method. One of `beer`, `edger`, `larman_gp`, `xu_zigp` |
| `--threshold` | `2.3` | calling threshold for the Python scorers, on the -log10 p scale |
| `--gp-null` | `beads` | set the generalized Poisson scorers build their null from, `beads` or `input` |
| `--gp-bins` | one per 50 peptides | abundance bins for `larman_gp` and `xu_zigp`, capped at 300 |
| `--replicate-rule` | off | give every replicate in a group the lowest score in that group |
| `--replicate-group-column` | `participant_ID` | what counts as one group for the replicate rule |
| `--min-lib-size` | `500` | drop a sample below this many library-mapped reads |
| `--posterior-threshold` | `0.5` | posterior probability above which a peptide is called |
| `--beads-rr` | off | run each beads-only sample against the others |
| `--edger-fdr` | `0.05` | BH cutoff for the edgeR hit matrix |
| `--seed` | `20260101` | sampler seed, pinned so a rerun reproduces the calls. The default value carries no meaning beyond being fixed |
| `--keep-workdir` | none | copy the R working directory here for inspection |

### `phipstream prioritise`

| Option | Default | Purpose |
|---|---|---|
| `--posterior` | required | beer_posterior.csv.gz |
| `--hits` | required | beer_hits.csv.gz |
| `--sample-table` | required | sample table with technical_replicate_id |
| `--peptide-table` | required | peptide table with peptide_id and the grouping columns |
| `--out-dir` | required | destination for the concordance and shortlist tables |
| `--group-column` | `participant_ID` | column that names the replicate group |
| `--role-column` | `sample_role` | column that names the sample class |
| `--role` | `serum` | sample class to prioritise, empty for all |
| `--min-replicates` | `2` | replicates of one group a peptide must be enriched in |
| `--min-groups` | `2` | groups a peptide must recur across to be kept |
| `--all-replicates` | `3` | replicates that count as enrichment in every replicate |
| `--tile-step` | `28` | tile start offset of an adjacent overlapping tile |
| `--antigen-columns` | `virus,antigen` | peptide table columns an adjacent tile must share |
| `--position-column` | `tile_start` | peptide table column holding the tile start |
| `--posterior-threshold` | none | recompute calls at this posterior instead of trusting the hit matrix |

### `phipstream stages`

| Option | Default | Purpose |
|---|---|---|
| `--resume` | off | skip a stage whose output is already present |
| `--workflow` | off | render the `[workflow]` section and run the bundled workflow instead of the stages |
| `--submit` | off | write a PBS job script and submit it instead of running |
| `--computerome_project` | none | PBS account and group list, required with --submit |
| `--dry-run` | off | print the stage commands without running them |

For `--workflow --submit`, the `modules`, `walltime`, `ppn` and `mem` keys are
forwarded to the job script, and its logs are written under `<out_dir>/logs`
alongside the stage route's.

## Repository layout

```text
bin/
  phipstream              command dispatcher
  phipstream-setup        environment check used by make check
  containers.py           image tags and runtime wrapper
  make_qc_samplesheet.py  seqinspector sample sheet writer
  run_fastqc.py           FastQC stage, read for the contamination diagnostic
  trim_reads.py           cutadapt stage
  align_reads.py          Bowtie 2 alignment and count stage
  call_enrichment.py      edgeR and BEER scoring wrapper
  prioritise_peptides.py  replicate support and shortlist builder
  run_stages.py           runs every stage from one config, or submits it
NOTICE                    which license applies to which part of the tree
configs/
  _template.config        dataset config template for the workflow mode
  _template.stages.conf   parameter file for the stage route
  adapters.example.csv    adapter specification example
docs/
  computerome.md          PBS and Apptainer deployment notes
phip-flow/                bundled phip-flow fork
tests/                    synthetic fixtures and end to end checks
```

Generated directories such as `results/`, `work/`, and `logs/` are excluded from
version control.

## Troubleshooting

**No container engine is reported by `make check`.** Install Docker or
Apptainer. On a cluster, load the required module before running `make check`.

**Image retrieval occurs during the first analysis.** Execute `make pull` before
analysis to retrieve the required images explicitly.

**The scoring stage reports a package cache permission error.** The container
requires a writable HOME directory. phipstream sets HOME to the working
directory, so verify write permissions for that directory.

**The scoring stage terminates because beads-only controls are unavailable.**
Add samples with `control_status=beads_only`, or provide the metadata columns
used to derive that status.

**Posterior columns are empty or all peptides are negative.** Review
`align_summary.csv`. Low library-mapped depth or poor alignment can indicate an
adapter mismatch, an incorrect library table, or reads from a different
construct.

## Citations

### Software invoked by the pipeline

Cutadapt 5.2:
Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput
sequencing reads. *EMBnet.journal, 17*(1), 10-12.
https://doi.org/10.14806/ej.17.1.200

Bowtie 2 2.5.4:
Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie
2. *Nature Methods, 9*(4), 357-359. https://doi.org/10.1038/nmeth.1923

SAMtools 1.19.2:
Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.
O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021).
Twelve years of SAMtools and BCFtools. *GigaScience, 10*(2), giab008.
https://doi.org/10.1093/gigascience/giab008

BEER 1.2.0:
Chen, A., Kammers, K., Larman, H. B., Scharpf, R. B., & Ruczinski, I. (2022).
Detecting and quantifying antibody reactivity in PhIP-Seq data with BEER.
*Bioinformatics, 38*(19), 4647-4649.
https://doi.org/10.1093/bioinformatics/btac555

edgeR, bundled in the BEER image:
Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: A Bioconductor
package for differential expression analysis of digital gene expression data.
*Bioinformatics, 26*(1), 139-140. https://doi.org/10.1093/bioinformatics/btp616

pandas 2.0 or later:
McKinney, W. (2010). Data structures for statistical computing in Python. In
*Proceedings of the 9th Python in Science Conference* (pp. 56-61).
https://doi.org/10.25080/Majora-92bf1922-00a

statsmodels, required by the generalized Poisson scorers:
Seabold, S., & Perktold, J. (2010). statsmodels: Econometric and statistical
modeling with Python. In *Proceedings of the 9th Python in Science Conference*
(pp. 92-96). https://doi.org/10.25080/Majora-92bf1922-011

The construction implemented by the larman_gp scorer:
Larman, H. B., Zhao, Z., Laserson, U., Li, M. Z., Ciccia, A., Gakidis, M. A.
M., Church, G. M., Okada, S., Ndung'u, T., Walker, B. D., & Elledge, S. J.
(2011). Autoantigen discovery with a synthetic human peptidome. *Nature
Biotechnology, 29*(6), 535-541. https://doi.org/10.1038/nbt.1856

The zero inflated variant in xu_zigp, and the criterion behind --replicate-rule:
Xu, G. J., Kula, T., Xu, Q., Li, M. Z., Vernon, S. D., Ndung'u, T.,
Ruxrungtham, K., Sanchez, J., Brander, C., Chung, R. T., O'Connor, K. C.,
Walker, B., Larman, H. B., & Elledge, S. J. (2015). Comprehensive serological
profiling of human populations using a synthetic human virome. *Science,
348*(6239), aaa0698. https://doi.org/10.1126/science.aaa0698


Bowtie 1, optional in the Nextflow execution mode:
Langmead, B., Trapnell, C., Pop, M., & Salzberg, S. L. (2009). Ultrafast and
memory-efficient alignment of short DNA sequences to the human genome. *Genome
Biology, 10*(3), R25. https://doi.org/10.1186/gb-2009-10-3-r25

### Infrastructure and workflow components

Nextflow:
Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., &
Notredame, C. (2017). Nextflow enables reproducible computational workflows.
*Nature Biotechnology, 35*(4), 316-319. https://doi.org/10.1038/nbt.3820

Apptainer and Singularity:
Kurtzer, G. M., Sochat, V., & Bauer, M. W. (2017). Singularity: Scientific
containers for mobility of compute. *PLOS ONE, 12*(5), e0177459.
https://doi.org/10.1371/journal.pone.0177459

BioContainers:
da Veiga Leprevost, F., Gruning, B. A., Alves Aflitos, S., et al. (2017).
BioContainers: An open-source and community-driven framework for software
standardization. *Bioinformatics, 33*(16), 2580-2582.
https://doi.org/10.1093/bioinformatics/btx192

nf-core, for seqinspector-compatible QC sheets:
Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A.,
Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for
community-curated bioinformatics pipelines. *Nature Biotechnology, 38*(3),
276-278. https://doi.org/10.1038/s41587-020-0439-x

### Software components without formal publications

* phip-flow: https://github.com/matsengrp/phip-flow. The bundled fork is retained
  under `phip-flow/` with its GPL-3.0 license in `phip-flow/LICENSE`.
* phippery: https://matsen.group/phippery/. The Nextflow execution mode uses
  phippery for the counts object and CPM tables.
* PhIPData: https://bioconductor.org/packages/PhIPData/. PhIPData is the
  Bioconductor container class consumed by BEER and is described in the Chen et
  al. (2022) BEER publication above.
* nf-core/seqinspector: https://github.com/nf-core/seqinspector.
* BiocParallel: https://bioconductor.org/packages/BiocParallel/. The R scoring
  scripts use BiocParallel to enforce serial execution when parallel workers
  would increase memory consumption.

## Authors

James Henderson Pang.

## License

This repository contains code under more than one license. phipstream code
outside `phip-flow/` is licensed under the MIT License. The bundled phip-flow
fork under `phip-flow/` is distributed under the GNU General Public License,
version 3.0. [NOTICE](NOTICE) maps which terms apply where, [LICENSE](LICENSE)
carries the MIT text, and [phip-flow/LICENSE](phip-flow/LICENSE) carries the
GPL-3.0 text.
