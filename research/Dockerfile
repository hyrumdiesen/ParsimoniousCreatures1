FROM julia:1.12.6 AS base

ENV JULIA_NUM_THREADS=auto

RUN apt-get update && \
    apt-get install -y ghostscript && \
    rm -rf /var/lib/apt/lists/*

FROM base AS ode

WORKDIR /research/ode
COPY ode/ /research/ode/

RUN julia --project=/research/ode -e 'using Pkg; Pkg.instantiate()'

FROM base AS analytic

WORKDIR /research/analytic
COPY analytic/ /research/analytic/

RUN julia --project=/research/analytic -e 'using Pkg; Pkg.instantiate()'
