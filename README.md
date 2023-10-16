# Virtual Spike-In

This repository contains code for Virtual Spike-In, a hierarchical bayesian model for inference of inter-sample normalization values in sequencing data using a counts table over either an invariant region of the genome or from an exogenous spike-in.
This repository includes scripts used in the initial development of this model, the model itself and an associated NextFlow pipeline to automatically run the model on a selected dataset using a 3' invariant region, as well as scripts used for generating the figures used in the final manuscript. Relevant folders/scripts and their purposes are described below.

## Installation

## Python Dependencies

The python components of this package are presented in a PEP518 compatible `pyproject.toml` file. The project was developed using `poetry`, which is also the easiest way to install the locked versions of the dependencies of the project.

To install the dependencies using poetry, run the following commands:
```
# Clone the repository
git clone <repository>
# Install poetry
cd virtual_spike_in
# <Do whatever virtualenv stuff you want to do>
pip install poetry
# Install the dependencies
poetry install
```

## Pipeline Dependencies

These dependencies should be automatically managed using NextFlow if available on your cluster.
- Nextflow - required for running the pipeline
- Bedtools - used for isoform filtering
- Samtools - used for converting files between different formats


## Optional Dependencies

- R - used for generating figures
	- Tidyverse - used for data manipulation

## Common Issues

- `pymc` is picky about how libraries are linked. If you are using something like `pyenv` to manage your python environments, you might have to install your python version using shared libraries, as follows:
```
env PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.8.0
```
- To squeeze the most possible performance out of your python installation, you can use the following flags to compile python with shared libraries also enabled:
```
env PYTHON_CONFIGURE_OPTS='--enable-shared --enable-optimizations --with-lto' PYTHON_CFLAGS='-march=native -mtune=native' pyenv install 3.8.0
```

## Repository Structure

The repository is structured as a standard python package with some additional folders for project and analysis specific scripts.
Additional definitions are provided for users of the nix package manager which, along with the poetry definitions, define the full development environment used for this project along with locked versions.

- `/` - repository root
  - `pyproject.toml` - python dependency definitions
  - `poetry.lock` - python dependency lockfile
  - `shell.nix` - nix package manager definition for a fully encapsulated development shell for the project
  - `flake.nix` - a wrapper for `shell.nix` to use nix's immutable flake definitions for faster shell access
  - `flake.lock` - locked versions of packages used for the development of this projecta
- `/src` - scripts used for analysis
  - `human` - remnant scripts from initial model development (3' end method)
  - `drosophila` - remnant scripts from initial model development (exogenous spike-in method)
  - `dukler` - scripts for the final paper using data from the Dukler2017Nascent dataset
	- `subsample_crams.sbatch` - remnant script for a proposed analysis based on subsampling "deep" spike-ins in the literature
	- `fastestq_dump.sbatch` - fast script for fetching SRR files from the sequence read archive. Can typically saturate downlink bandwidth.
	- `gen_comparison_plots.r` - script used to generate Figure 3 from the manuscript, in addition to auxillary plots
	- `drosophila_pipeline.sh` - small wrapper to run the Nascent-Flow pipelineon the Dukler2017 dataset for downstream analysis
	- `check_lit_spikein_ratio.r` - script used to generate figures comparing the depth of spikeins in literature samples.
	- `srrs.txt`
	- `all_millionsmapped.txt`
	- `mapped_ratio.txt` -
- `/virtual_spike_in` - model implementation
  - `main.py` - CLI wrapper for running the model
  - `virtual_spike_in.py` - the final implementation of the VSI model
- `/pipeline` - NextFlow pipeline implementation
  - `main.nf` - the main definition of the NextFlow pipeline
  - `nextflow.config` - main configuration file for NextFlow
  - `conf` - additional nextflow configuration files
	- `base.config` - root config file that others are defined from
	- `example.config` - sample experimental configuration
  - `bin` - nextflow auxillary scripts
	- `calc_maximal_isoform.bash` - script for isoform filtering using maximally expressed isoform
	- `merge_counts_with_join.r` - script for merging separate featureCounts counts tables so we can run all counts in parallel on different nodes and accelerate the pipeline.
- `figs` - Inkscape SVGs used for creating figures for the paper

## LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
