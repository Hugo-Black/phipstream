# PBS cluster notes

These notes describe the workflow manager path on Computerome and similar PBS or
Torque clusters where Docker is not available on compute nodes.

## Account setup

Load the site modules first. Module names can change, so check `module avail` if
a listed name is not accepted.

```bash
module load tools apptainer nextflow
cd $HOME/phipstream
make check
```

Keep the repository somewhere visible from both login and compute nodes. On
Computerome, `$HOME` works for this. The default Apptainer cache lives at
`~/.apptainer/cache/`, which also needs to be visible to batch jobs.

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

The stage commands need an interpreter with pandas. The workflow manager route
does not.

Prefer a site module that already carries the scientific stack over building an
environment. On Computerome that is anaconda:

```bash
module load tools anaconda3/2024.06-1 apptainer/1.4.5
python3 -c "import pandas; print(pandas.__version__)"
```

pandas is the only Python dependency, so that module is enough and nothing needs
installing.

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

This route runs trimming, alignment, enrichment calling, concordance and
prioritisation as one job. It is the route that produces a shortlist.

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

Before committing to a full run, set `method = edger` in the config and submit
once. That stops after edgeR rather than running the sampler, which took seconds
where BEER took hours on the same data, and it exercises every stage including
prioritisation. Switch back to `method = beer` once the deployment is proven,
and expect a different and much larger set of calls.

Enrichment calling is the long stage. It ran for about three hours on 96 samples
against 1,194 peptides, and it runs on one core because the sampler is forced to
serial execution. Set `walltime` with room to spare, because losing the job at
the wall loses that time.

If a job does die, resubmit it. Every stage is skipped when its output is already
present, so a rerun continues from the stage that failed rather than repeating
the sampler.

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

## Resume after failure, workflow manager route

Nextflow keeps finished tasks under `work/`. After correcting the problem, pass
`--resume`:

```bash
bin/phipstream nextflow configs/mydataset.config --resume
```

Changing an input timestamp or a config value can invalidate part of the cache.
For example, editing the sample table reruns tasks that depend on it.

## Frequent problems

**`apptainer: command not found` inside the job.** The compute job did not load a
module that provides Apptainer. The generated script loads modules itself, so
check that the names in the script match the cluster.

**The job stays queued.** Run `qstat -f <jobid>` and inspect the scheduler
reason. The usual causes are an invalid project code, a full allocation, or a
request that does not match available resources.

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
