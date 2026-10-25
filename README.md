Code and processed data to accompany manuscript titled "Parsimonious Creatures: Building Gene Regulatory Networks Robust to Somatic Evolution". Raw data can be found on Zenodo.

Within the parent research folder, there are two primary folders: `ode` and `analytic`. To maximize reproducibility, all code used in the paper can be run from one of two docker images which are fully separated into either of these folders.

The `ode` folder contains processed datasets, as well as scripts for plotting data, for running the deterministic model of Parsimonious Creatures, and for running general stochastic simulations of Parsimonious Creatures using `OrdinaryDiffEq` to solve GRN states of any network topology. Finally, a `playground.jl` file is included for getting familiar with running the simulations of Parsimonious Creatures using any GRN topology, as well as some of the GRNs described in the manuscript.

The `research` directory contains specialized code for running highly optimized simulations of Parsimonious Creatures using analytic solutions of network topologies which are analytically solvable and thus are much faster to compute. All raw datasets may be reproduced here, but may take days or may require a cluster to complete. This environment is kept separate from `ode` due to `BayesianOptimization` incompatibilities with `OrdinaryDiffEq`.

To create both docker images please follow the below steps:

Install git and docker if necessary, then use the command line to navigate to a directory where you wish to install the git repo:
```
cd ~/your/path
```
Download the entire git repository:
```
git clone https://github.com//hyrumdiesen/ParsimoniousCreatures1

cd research
```

Build both docker images:

```
docker build --target ode -t ode .

docker build --target analytic -t analytic .
```

If you wish to run the ode environment:
```
docker run -d --name my-simulation -e JULIA_NUM_THREADS=auto -v "${PWD}\data:/research/ode/" ode
```

If you wish to run the analytic environment:
```
docker run -d --name my-simulation -e JULIA_NUM_THREADS=auto -v "${PWD}\data:/research/analytic/" analytic
```
