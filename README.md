Julia code and processed data to accompany manuscript titled "Parsimonious Creatures: Building Gene Regulatory Networks Robust to Somatic Evolution". All raw data is included, except ~3gb of raw heatmap datasets which can be found on Zenodo.

Within the parent research folder, there are two primary folders: `ode` and `analytic`. To maximize reproducibility, all code used in the paper can be run from one of two docker images which are fully separated into either of these folders.

The `ode` folder contains processed datasets, as well as scripts for plotting data, for running the deterministic model of Parsimonious Creatures, and for running general stochastic simulations of Parsimonious Creatures using `OrdinaryDiffEq` to solve GRN states of any network topology. Finally, a `playground.jl` file is included for getting familiar with running the simulations of Parsimonious Creatures.

The `research` directory contains specialized code for running highly optimized simulations of Parsimonious Creatures using analytic solutions of network topologies which are analytically solvable and thus are much faster to compute. All raw datasets may be reproduced here, but may take days or may require a cluster to complete. This environment is kept separate from `ode` due to dependency incompatibilities.

To best reproduce data and figures, please follow the below steps:
## 1. Download the repository

Open a terminal and navigate to the directory where you want to store the project. Clone the repository and enter the `research/` directory:

```bash
git clone https://github.com/hyrumdiesen/ParsimoniousCreatures1.git
cd ParsimoniousCreatures1/research
```

The remaining Docker commands in this README assume that the terminal is in `research/`, where the Dockerfile is located.

## 2. Build the Docker images

Build the ODE environment:

```bash
docker build --target ode -t parsimonious-ode .
```

Build the analytic environment:

```bash
docker build --target analytic -t parsimonious-analytic .
```

The first build downloads the Julia base image and installs the dependencies. Subsequent builds can reuse Docker's cached layers when the relevant files have not changed.

The two images are independent. Build only the environment needed for your analysis if you do not need both.
## 3. Run code in docker container

To generate all plots used to generate figures:

```
docker run --rm -it ode bash /research/ode/run_all_figures.sh
```

For a more interactive experience of parsimonious creatures, you can explore the script `playground.jl`. To do this, launch Julia within the docker container in the terminal using:

```
docker run --rm -it ode julia --project=/research/ode
```

Or alternatively for running the docker container in VSCode:

```
docker run -d -it --name my_julia_ode_env ode sleep infinity
```

Then use dev containers to attach to running container.
## Citation

If you use this code or build on these analyses, please cite the associated manuscript:

_Parsimonious Creatures: Building Gene Regulatory Networks Robust to Somatic Evolution._

Raw data are available through Zenodo; consult the manuscript or the repository's associated data record for the relevant dataset and citation information.
```