# PHIP-FLOW

This directory contains the bundled phip-flow workflow used by phipstream when a
Nextflow deployment is preferred. The workflow is based on the public
matsengrp/phip-flow project for Common Phage Immuno-Precipitation Sequencing
experiments.

Upstream documentation for phippery and phip-flow is available at:

https://matsengrp.github.io/phippery/introduction.html

## Basic upstream run

Install Nextflow:

```bash
curl -s https://get.nextflow.io | bash
```

Install Docker if your system supports it:

```bash
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
```

Run a tagged upstream release:

```bash
nextflow run matsengrp/phip-flow -r V1.12 -profile docker
```

The `-r` option selects the release. Use a version tag for stable runs. The
`main` branch is useful for development testing but should not be treated as a
stable analysis target.

Pipeline parameters are documented here:

https://matsengrp.github.io/phippery/alignments-pipeline.html#parameters

## Container note

The phippery image includes the main phippery environment. edgeR and BEER are
served from separate public images. phipstream pins those images in
`bin/containers.py` and routes its scoring path through them.

[![Docker Repository on Quay](https://quay.io/repository/hdc-workflows/phippery/status "Docker Repository on Quay")](https://quay.io/repository/hdc-workflows/phippery)

## License

This bundled phip-flow fork is distributed under the GNU General Public License,
version 3.0. See `LICENSE` in this directory. The repository root license file
contains the overall license map for the mixed-license distribution.
