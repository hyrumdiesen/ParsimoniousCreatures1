#this script consolidates all the plotting scripts in the /plotting folder to produce pdfs for all figures in the paper
#pdfs are sent to the output_figures folder

#!/usr/bin/env bash
set -e

julia --project=/research/ode scripts/plotting/Nvst_plotting.jl
julia --project=/research/ode scripts/plotting/plot_weighted_network.jl
julia --project=/research/ode scripts/plotting/plot_selectedE.jl
julia --project=/research/ode scripts/plotting/plotheatmaps.jl
julia --project=/research/ode scripts/plotting/plothybrids.jl

echo "All figures generated successfully."