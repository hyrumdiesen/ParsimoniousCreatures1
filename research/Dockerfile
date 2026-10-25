FROM julia:1.12.6 AS base
#keep this docker file in the /research directory, I think it will keep different julia
#instances separate?
WORKDIR /research

FROM base AS ode
COPY ode/ /research/

RUN julia --project=/research -e 'using Pkg; Pkg.instantiate()'
CMD ["julia", "--project=/research", "scripts/run_simulations.jl"]

FROM base AS analytic
COPY analytic/ /research/
RUN julia --project=/research -e 'using Pkg; Pkg.instantiate()'
CMD ["julia", "--project=/research", "scripts/Nvst_datagen.jl"]

# Copy environment definitions first
#COPY Project.toml Manifest.toml* ./
#RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Copy all your local .jl scripts into the container
#COPY . .

# Execute the simulation runner when the container starts
#CMD ["julia", "--project=.", "run_selectedE.jl"]

