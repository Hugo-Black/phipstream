#!/usr/bin/env python3
"""Container runtime helpers for the pipeline stages.

Stage scripts execute tool commands inside pinned Docker or Apptainer images.
The image table in this file is therefore the source of truth for external tool
versions.
"""
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

IMAGES = {
    "fastqc": "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0",
    "cutadapt": "quay.io/biocontainers/cutadapt:5.2--py312hfabe715_2",
    "bowtie2": "quay.io/biocontainers/bowtie2:2.5.4--he96a11b_7",
    "samtools": "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1",
    "beer": "quay.io/biocontainers/bioconductor-beer:1.2.0--r42hdfd78af_0",
}


def engine():
    """Select Docker or Apptainer from the environment or the host PATH."""
    want = os.environ.get("PHIPSTREAM_ENGINE", "").strip().lower()
    if want:
        if want not in ("docker", "apptainer"):
            sys.exit(f"PHIPSTREAM_ENGINE must be docker or apptainer, got {want!r}")
        return want
    if shutil.which("docker") and subprocess.run(
            ["docker", "info"], capture_output=True).returncode == 0:
        return "docker"
    if shutil.which("apptainer") or shutil.which("singularity"):
        return "apptainer"
    sys.exit("no container engine found. Install Docker or Apptainer, "
             "or set PHIPSTREAM_ENGINE")


def cache_dir():
    """Return the Apptainer cache path used by the command wrappers."""
    for name in ("NXF_APPTAINER_CACHEDIR", "APPTAINER_CACHEDIR",
                 "SINGULARITY_CACHEDIR"):
        value = os.environ.get(name)
        if value:
            return Path(value)
    return Path.home() / ".apptainer" / "cache"


def sif_path(key):
    """Where `phipstream images pull` leaves the Apptainer image for a tool."""
    return cache_dir() / (IMAGES[key].replace("/", "-").replace(":", "-") + ".sif")


def image_uri(key, eng):
    """Image reference to execute.

    Apptainer prefers an image file that has already been built. Resolving a
    docker:// reference instead makes every task convert the image again, and
    concurrent tasks racing on that conversion fail. Fetching the images once
    with `phipstream images pull` avoids both.
    """
    uri = IMAGES[key]
    if eng != "apptainer":
        return uri
    local = sif_path(key)
    return str(local) if local.exists() else f"docker://{uri}"


def build(key, command, binds, workdir):
    """Build the argv used to execute a shell command inside a tool image.

    Host directories are mounted at identical paths inside the container, so the
    command can use one path scheme. HOME is pointed at the work directory so R
    packages and other tools have a writable cache location. For Apptainer, the
    same value is passed through --home.
    """
    eng = engine()
    work = str(Path(workdir).resolve())
    mounts = sorted({str(Path(b).resolve()) for b in binds if b} | {work})
    if eng == "docker":
        argv = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
                "-e", f"HOME={work}"]
        for m in mounts:
            argv += ["-v", f"{m}:{m}"]
        argv += ["-w", work, image_uri(key, eng), "sh", "-c", command]
        return argv
    exe = "apptainer" if shutil.which("apptainer") else "singularity"
    argv = [exe, "exec", "--home", work]
    for m in mounts:
        argv += ["-B", f"{m}:{m}"]
    argv += ["--pwd", work, image_uri(key, eng), "sh", "-c", command]
    return argv


def run(key, command, binds, workdir, capture=False, check=True):
    argv = build(key, command, binds, workdir)
    if os.environ.get("PHIPSTREAM_TRACE"):
        print("+ " + " ".join(shlex.quote(a) for a in argv), file=sys.stderr)
    proc = subprocess.run(argv, capture_output=capture, text=True)
    if check and proc.returncode != 0:
        if capture:
            sys.stderr.write(proc.stdout or "")
            sys.stderr.write(proc.stderr or "")
        sys.exit(f"{key} step failed with exit code {proc.returncode}")
    return proc


def pull_all():
    """Pre-fetch all pinned images for offline or scheduled execution.

    Run this on a login node before submitting cluster jobs that will not have
    registry access at runtime.
    """
    eng = engine()
    if eng == "docker":
        for key in IMAGES:
            uri = image_uri(key, eng)
            print(f"[pull] {uri}")
            subprocess.run(["docker", "pull", uri], check=True)
        return
    exe = "apptainer" if shutil.which("apptainer") else "singularity"
    target = cache_dir()
    target.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, APPTAINER_CACHEDIR=str(target),
               SINGULARITY_CACHEDIR=str(target))
    for key in IMAGES:
        # Always fetch from the registry here. image_uri would hand back the
        # local file once one exists, which is right for running and wrong for
        # fetching.
        uri = f"docker://{IMAGES[key]}"
        print(f"[pull] {uri}")
        subprocess.run([exe, "pull", "--force", str(sif_path(key)), uri],
                       check=True, env=env)
    print(f"[pull] images cached under {target}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "pull":
        pull_all()
    else:
        print(f"engine: {engine()}")
        for k, v in IMAGES.items():
            print(f"  {k:9s} {v}")
