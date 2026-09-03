# PBS cluster notes

These notes cover both routes on Computerome and similar PBS or Torque clusters
where Docker is not available on compute nodes.

## Account setup

Load the site modules first, naming exact versions.

```bash
module load tools anaconda3/2024.06-1 apptainer/1.4.5 nextflow/25.10.2
cd /path/to/phipstream
make check
```

**Pin the versions.** A site that publishes no default for a module rejects a
bare name, and `module load apptainer` there fails with `Unable to locate a
modulefile for 'apptainer/<nodefault>'`. Check what exists with `module avail
apptainer`, then record the exact names in the dataset config's `modules` key so
submitted jobs load the same set. Both routes take that key from the config.

Keep the repository somewhere visible from both login and compute nodes. The
default Apptainer cache lives at `~/.apptainer/cache/`, which also needs to be
visible to batch jobs.

Fetch images before submitting work. Run this on a login node, which has
registry access, and only once:

```bash
make pull
```

This matters for more than offline nodes. The commands run the image files that
`make pull` writes. Without them each task resolves the registry reference
itself, and tasks running in parallel race on that conversion and fail.

If your home quota is too small, set `NXF_APPTAINER_CACHEDIR` to a project
filesystem with enough space and rerun the pull.

## Python

The stage commands need an interpreter with pandas, and the generalized Poisson
scorers additionally need statsmodels. The workflow manager route needs neither,
since it runs everything in containers.

Prefer a site module that already carries the scientific stack over building an
environment. On Computerome that is anaconda:

```bash
module load tools anaconda3/2024.06-1
python3 -c "import pandas, statsmodels; print('ok')"
```

Those are the only Python dependencies, so that module is enough and nothing
needs installing.

Do not build a virtual environment or run `pip install` on a login node. Those
nodes cap per-user memory and terminate the process with a bare `Killed` message
part way through, which leaves a partially written package tree behind. If a
site module with pandas is genuinely unavailable, download the wheels on the
login node with `pip download`, then install them from that directory on an
interactive compute node.

Record the interpreter in the dataset config as the `python` key so a batch job
uses it rather than whatever is on the node PATH. A later `module load` can put
an older system Python back in front.

## Dataset setup, stage route

This route runs read quality reporting, trimming, alignment, enrichment calling,
concordance and prioritisation as one job. It is the route that produces a
shortlist.

```bash
cp configs/_template.stages.conf configs/mydataset.stages.conf
```

Edit the paths, the `python` key, and the `modules` line so it loads what this
site actually provides. Check the plan without running anything:

```bash
bin/phipstream stages configs/mydataset.stages.conf --dry-run
```

Then submit it:

```bash
bin/phipstream stages configs/mydataset.stages.conf \
    --submit --computerome_project <project>
```

The project value is the allocation charged by PBS. Without it the job script is
written but not submitted. Job output lands under the `logs/` directory inside
the configured `out_dir`.

Set `fastqc = true` unless the contamination diagnostic is definitely not
wanted. It adds a stage ahead of trimming and writes the archives that module
reads, and nothing else creates them.

Prove the deployment with a fast scorer before committing to BEER. `edger` and
`larman_gp` both finish in seconds where BEER takes hours, and both exercise
every stage including prioritisation. On 96 samples against 1,194 peptides
edgeR scored 3 of 18 donor groups and shortlisted nothing, while `larman_gp`
scored all 18, so the second is the more useful check as well as the faster one.

BEER is the long stage when it is selected. It ran for about three hours on that
dataset and uses one core, because the sampler is forced to serial execution.
Set `walltime` with room to spare, since losing the job at the wall loses that
time.

If a job does die, resubmit it. Every stage is skipped when its output is already
present, so a rerun continues from the stage that failed rather than repeating
the sampler.

## Choosing a route

The stage route submits one job. The workflow manager route submits one job per
task, which for a hundred samples is several hundred of them, each asking for
the profile's default of 4 cores and 16 GB.

That trade is worth taking when a task is real work. Aligning a hundred
thousand peptides is. Counting reads against a library of one thousand is not,
and on such a dataset the stage route finishes the same alignment and counting
in a single job in well under a minute, where the workflow route spends longer
than that waiting in a queue.

Use the workflow route for the scoring layers and analysis modules the stage
route does not provide, and the stage route for everything else, unless the
library is large enough that spreading alignment across nodes pays for the
scheduling overhead.

Submissions are capped at 20 in flight at 10 per minute. Raise `queueSize` in
the profile only after checking what the site tolerates.

## Dataset setup, workflow manager route

This route distributes alignment across nodes but stops after enrichment
calling. It has no concordance or prioritisation step, and its counts are keyed
by row order rather than by sample name, so they are not accepted by the stage
route's scoring step. What it adds is the scoring layers the stages leave out,
CPM, z scores and gamma-Poisson, together with the library QC and contamination
modules.

Prefer driving it from the `[workflow]` section of the stage config, so one file
describes the dataset for both routes:

```bash
bin/phipstream stages configs/mydataset.stages.conf --workflow \
    --submit --computerome_project <project>
```

The contamination module reads FastQC archives rather than reads, and nothing in
the workflow creates them. Set `fastqc = true` in the stage config and run the
stage route first. `fastqc_dir` is then pointed at that output automatically.

A standalone dataset config still works for a run that needs nothing from the
stage route:

```bash
cp configs/_template.config configs/mydataset.config
bin/phipstream submit configs/mydataset.config --computerome_project <project>
```

Before launching a large dataset, run a small one interactively on a login node:

```bash
bin/phipstream nextflow configs/mydataset.config
```

That catches path errors and missing cache files before queue time is spent.

## Resume after failure

Both routes resume. The stage route skips any stage whose output is already
present, so a resubmission continues from the one that failed.

The workflow route keeps finished tasks under `<out_dir>/work` and needs to be
told:

```bash
bin/phipstream stages configs/mydataset.stages.conf --workflow --resume \
    --submit --computerome_project <project>
```

A run stopped part way through can pick up cheaply. One here reported
`Succeeded: 6, Cached: 398` and finished in under three minutes, having reused
the alignment, counting, library QC and contamination it had already done.

Changing an input timestamp or a config value can invalidate part of the cache.
For example, editing the sample table reruns tasks that depend on it.

## Where output goes

Everything a run produces is written under the `out_dir` from the config.

```text
<out_dir>/
  fastqc/        FastQC archives, read by the contamination diagnostic
  trimmed/       trimmed reads
  alignment/     counts matrix and per sample alignment summary
  enrichment/    score and hit matrices for the selected method
  prioritised/   shortlist and replicate concordance
  workflow/      the workflow route's layers, QC and contamination
  logs/          job scripts, scheduler output, Nextflow report and trace
  work/          the workflow route's scratch space
```

`work/` holds every intermediate file the workflow route produces and is the
largest thing written. It can be deleted once a run has finished, at the cost
of not being able to resume it.

A standalone `phipstream submit` has no `out_dir` and instead follows the
config's `results` path, provided that path is absolute. Otherwise it falls
back to `logs/` in the repository. `--log-dir` and `--work-dir` override both.

## Frequent problems

**`Unable to locate a modulefile for '<name>/<nodefault>'`.** The `modules` key
names a module the site publishes no default version for. Run `module avail
<name>` and put an exact version in the config.

**`apptainer: command not found` inside the job.** The `modules` key does not
name a module providing Apptainer.

**The job stays queued.** Run `qstat -f <jobid>` and inspect the scheduler
reason. The usual causes are an invalid project code, a full allocation, or a
request that does not match available resources.

**An R step fails with `cannot create dir '/home/<user>'` or a package failing
to load.** The site is not mounting home into the container, so R cannot
resolve `~` and PhIPData in particular fails on load. Both apptainer scopes
pass `--home $PWD`, which puts the container's home in the task directory. If a
site mounts home but read only, the same option applies.

**`env: 'apptainer': No such file or directory` inside a task.** Nextflow
submits every process as its own job, and those jobs inherit nothing from the
script that submitted them. The `process_modules` parameter is loaded by each
task for exactly this reason, and it must name a module providing Apptainer.
The stage route fills it from the config's `modules` key.

**`Unknown queue MSG=requested queue not found`.** The workflow route submits
each process to the queue named by the `queue` parameter, which defaults to
`batch`. List what the site actually offers with `qstat -Q` and set `queue`
in the dataset config to one of them.

**Image conversion fails while writing the SIF file.** The cache is not writable
or the filesystem handles sparse files poorly. Point `NXF_APPTAINER_CACHEDIR` at
a project filesystem and fetch the images again.
