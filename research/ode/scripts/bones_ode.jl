using Random
using Distributions
using StatsBase
using Roots
using Plots
using Plots.PlotMeasures
using LinearAlgebra
using OrdinaryDiffEq
using SteadyStateDiffEq
using BlackBoxOptim
using SparseArrays
using InvertedIndices
using HDF5
using DataFrames
using CSV
using Dates
using GLMakie, GraphMakie, Graphs   # GLMakie → interactive; swap for CairoMakie for static/PDF output
using LaTeXStrings

using CairoMakie          # variable edge widths render here; GLMakie collapses them
CairoMakie.activate!() 
using ImageFiltering


NUM_GENOTYPES = 50
B = 1
GRN_TIMESCALE = 1
EPS = 0.25
m = 0.002
mutmag = 0.5
MAX_TIME = 1e7
DELTAI = 0.1
deltaO = 0.001
root = 12

WofN = function(N)
    exp(-0.5 * ((N - root) / 2)^2)
end

solveKCmatrix_signal = function(N, KCmatrix, rs, phis)

    numTFs = trunc(Int, sqrt(length(KCmatrix)))

    p = reshape(KCmatrix, numTFs, numTFs), 1, 1, N, rs
    u0 = ones(numTFs)*0.1
    #u0 = rand(numTFs)
    prob = SteadyStateProblem(GRNode_signal!, u0, p)
    sol = solve(prob, DynamicSS(Rosenbrock23()), abstol=1e-6)
    #sol = solve(prob, DynamicSS(Rodas5P()); abstol=1e-6)
    #sol = solve(prob, DynamicSS(Rodas5P(autodiff = false)); abstol = 1e-6)
    #sol = solve(prob, DynamicSS(Rodas5P(autodiff = ADTypes.AutoFiniteDiff())); abstol = 1e-6)
    #sol = solve(prob, DynamicSS(Tsit5()); abstol = 1e-6)

    #rho = 1/(1+mean(sol .*phis))
    #rho = EPS/(1+mean(sol .* phis))
    rho = EPS/(1+sum(sol .* phis))
    #rho = EPS*mean(sol .*phis)
    return rho
end

solveKCmatrix_anyTF_signal = function(N, KCmatrix, rs, phis, whichTF)

    numTFs = trunc(Int, sqrt(length(KCmatrix)))

    p = reshape(KCmatrix, numTFs, numTFs), 1, 1, N, rs
    #u0 = ones(numTFs)*0.1
    u0 = rand(numTFs)
    prob = SteadyStateProblem(GRNode_signal!, u0, p)
    sol = solve(prob, DynamicSS(Rosenbrock23()), abstol=0.000001)

    return sol[whichTF]

    #rho = 1/(1+mean(sol .*phis))
    #rho = EPS/(1+mean(sol .* phis))
    #rho = EPS/(1+sum(sol .* phis))
    #rho = EPS*mean(sol .*phis)
    #return rho
end

KCmatrixplot = function(KCmatrix, rs, phis)

    lambdalookup = zeros(50)

        #probvect = matrixfunc_celldeath(KCmatrix, rs, phis, deltaO, DELTAI)
        for i in 1:50
            lambdalookup[i] = i*solveKCmatrix_signal(i, KCmatrix, rs, phis)
            #lambdalookup[i] = i*stepfunc(i)
        end
        #println(sum(WofN.(1:50).*probvect)*(1/deltaO))
        Plots.plot(lambdalookup ./ collect(1:1:50), ylim = (0,EPS), linewidth = 5)
        vline!([12])
        hline!([0.1])
end

KCmatrix_anyTF_plot = function(KCmatrix, rs, phis, whichTF)
    sols = zeros(50)
        #probvect = matrixfunc_celldeath(KCmatrix, rs, phis, deltaO, DELTAI)
        for i in 1:50
            #lambdalookup[i] = i*solveKCmatrix_signal(i, KCmatrix, rs, phis)
            sols[i] = solveKCmatrix_anyTF_signal(i, KCmatrix, rs, phis, whichTF)
            #lambdalookup[i] = i*stepfunc(i)
        end
        #println(sum(WofN.(1:50).*probvect)*(1/deltaO))
        #plot(sols ./ collect(1:1:50))
        Plots.plot(sols, ylim = (0,1.0))
        vline!([12])
        #hline!([0.1])
end

foldcopies_weights = function(fold, KCmatrix, rs, phis)
    onefoldidxs = findall(KCmatrix .> 0)
    nfold = zeros(size(KCmatrix).*fold)

    numTFs = size(KCmatrix)[1]

    for idx in onefoldidxs
        for j in 1:(fold)
            nfold[idx[1]+(j-1)*numTFs, idx[2]+(j-1)*numTFs] = KCmatrix[idx]
        end
    end

    return nfold, repeat(rs, fold), repeat(phis, fold)./fold

end

checkerboard7 = function(fold)
    size = fold*2
    KCmatrix = zeros(size,size)
    #for i in 1:(size^2)
    #    KCmatrix[i] = 1 - i%2
    #end
    @inbounds for i in 1:size, j in 1:size
        if i%2 == 0
            mult = 1#3.04339#3.0506#10#2
        else
            mult = 0.3286#1
        end
        KCmatrix[i, j] = isodd(i + j) ? 1*mult : 0
    end
    #KCmatrix[1,size] = 1
    #KCmatrix[size, 1] = 1

    rs = zeros(size)
    phis = zeros(size)
    numsignals = 0
    for i in 1:size
        rs[i] = i%2
        phis[i] = 1 - i%2
        numsignals += 1 - i%2
    end
    rs .*= 0.1744#0.17522#9/12#1/12
    phis .*= 23135.7/(numsignals) 

    KCmatrix .*= 1e8/numsignals
    return KCmatrix, rs, phis
end

separatecircuits7 = function(fold)
    KCmatrix1,rs1,phis1 = checkerboard7(1)
    KCmatrix,rs,phis = foldcopies_weights(fold, KCmatrix1, rs1, phis1)
    return KCmatrix, rs, phis
end

feedforwardavgvar1 = function(nmods)
    KCmatrix = [0 0; exp(18.91928) 0]
    rs = [exp(7.01765), 0]
    phis = [0,exp(9.776)]

    KCmatrix,rs,phis = foldcopies_weights(nmods, KCmatrix,rs,phis)
    return KCmatrix,rs,phis
end

function GRNode_signal!(du, u, p, t)
    GRNmatrix, B, GRNts, N, rs = p
    numTFs = length(u)
    
    # Optimization: Avoid bounds checking with @inbounds
    @inbounds for i in 1:numTFs
        s = 0.0
        # Optimization: Enable SIMD for vectorization
        @simd for j in 1:numTFs
            s += GRNmatrix[i, j] * u[j]
        end

        if rs[i] == 0
            du[i] = GRNts * ((1/(1 + s)) - B * u[i])
        else
            #du[i] = GRNts * ((1/(1 + s + N*rs[i])) - B * u[i])
            du[i] = GRNts * ((1/(1 + s)) - B *(1+N*rs[i])* u[i])
        end
        #du[i] = GRNts * ((N*rs[i])/(1 + s + N*rs[i]) - B * u[i])
        #du[i] = GRNts * (1/(1 + s + N*rs[i]) - B * u[i])

    end
    
    return nothing # Optimization: Return nothing instead of du
end

#deterministic model of parsimonious creatures
function matrixfunc_celldeath(lambdas; deltaO = deltaO)
    #create a lookup table

    lambdalookup = lambdas .*(1:length(lambdas))

    #lambdalookup = zeros(50)
    #for i in 1:50
    #    lambdalookup[i] = i*solveKCmatrix_signal(i, KCmatrix, rs, phis)
    #end

    maxsize = 50
    mat = Matrix{Float64}(undef, maxsize, maxsize)

    # first row
    firstrow = vcat(-lambdalookup[1] - deltaO, 2 * DELTAI, zeros(maxsize - 2))
    mat[1, :] .= firstrow

    # middle rows
    for q in 2:(maxsize - 1)
        row = vcat(
            zeros(q - 2),
            lambdalookup[q-1],
            -lambdalookup[q] - deltaO - q * DELTAI,
            (q + 1) * DELTAI,
            zeros(maxsize - q - 1)
        )
        mat[q, :] .= row
    end

    # last row
    #lastrow = vcat(zeros(maxsize - 1), -deltaO - maxsize * DELTAI)
    #lastrow = vcat(zeros(maxsize - 1), -deltaO)
    lastrow = vcat(zeros(maxsize - 2), lambdalookup[49], -deltaO)
    mat[maxsize, :] .= lastrow
    mat[maxsize-1,maxsize] = 0.0

    #black hole
    #blackholerow = 
    #return mat
    b = vcat(-1.0, zeros(50 - 1))
    probvect = mat \ b
    probvect = probvect / sum(probvect)
    return probvect
    #return mat
end

determ_fitness = function(x; deltaO = deltaO)
    sum((x .*WofN.(1:length(x)))/deltaO)
end

function Nvst_sim(KCmatrix, initialrs, initialphis; MAX_TIME = 1000, deltaO = 0.0, m = m, mutmag = mutmag, NUM_GENOTYPES = 25)

    numTFs   = trunc(Int, sqrt(length(KCmatrix)))
    numphis  = sum(initialphis .> 0)
    numrs    = sum(initialrs   .> 0)
    KCmatrix = reshape(KCmatrix, numTFs, numTFs)

    # ---- Per-genotype arrays ------------------------------------------------
    phis     = zeros(Float64, NUM_GENOTYPES, numTFs)
    rs       = zeros(Float64, NUM_GENOTYPES, numTFs)
    numcells = zeros(Int,     NUM_GENOTYPES)
    ks       = zeros(Float64, NUM_GENOTYPES, numTFs, numTFs)
    states   = fill(0.1, NUM_GENOTYPES, numTFs)

    # ---- Division-rate cache ------------------------------------------------
    rhos            = zeros(Float64, NUM_GENOTYPES)
    rho_per_cell    = zeros(Float64, NUM_GENOTYPES, NUM_GENOTYPES)
    rho_cache_valid = falses(NUM_GENOTYPES, NUM_GENOTYPES)

    # ---- Build the base SteadyStateProblem once (remake is cheap) -----------
    u0_init  = fill(0.1, numTFs)
    p0       = (KCmatrix, B, GRN_TIMESCALE, 1, initialrs)
    base_prob = SteadyStateProblem(GRNode_signal!, u0_init, p0)
    solver    = DynamicSS(Rodas5P())

    # ---- Initialise first genotype ------------------------------------------
    numcells[1]  = 1
    ks[1, :, :] .= KCmatrix
    phis[1, :]  .= initialphis
    rs[1, :]    .= initialrs

    sol0 = solve(base_prob, solver; abstol=1e-6)
    states[1, :] .= sol0.u

    # ---- Static bookkeeping -------------------------------------------------
    nonzero_k_idx   = findall(!iszero, KCmatrix)
    nonzero_r_idx   = findall(!iszero, initialrs)
    nonzero_phi_idx = findall(!iszero, initialphis)
    numnonzeroks    = length(nonzero_k_idx)

    effectivem  = m * (numnonzeroks + numrs)
    mut_weights = Weights(Float64[numphis * 0, numrs, numnonzeroks])
    mut_types   = ["MUT_PHI", "MUT_R", "MUT_K"]

    time       = 0.0
    totalcells = 1
    fitness    = 0.0

    # ---- Discrete-time population snapshot ----------------------------------
    # Record totalcells at each integer t = 1, 2, ..., floor(MAX_TIME).
    # If the simulation ends before a checkpoint, the remaining entries stay 0.
    max_record_steps  = trunc(Int, MAX_TIME)
    totalcells_history = ones(Int, max_record_steps)
    next_record_time   = 1.0   # next integer checkpoint

    # =========================================================================
    # GILLESPIE LOOP
    # =========================================================================

    while time < MAX_TIME && totalcells > 0

        # ---- Division rates (cache-aware) -----------------------------------
        @inbounds for i in 1:NUM_GENOTYPES
            nc = numcells[i]
            if nc == 0
                rhos[i] = 0.0
                continue
            end

            if !rho_cache_valid[i, totalcells]
                p_new    = (ks[i, :, :], B, GRN_TIMESCALE, totalcells, rs[i, :])
                prob_new = remake(base_prob; u0=states[i, :], p=p_new)
                sol_new  = solve(prob_new, solver; abstol=1e-6)

                states[i, :] .= sol_new.u

                phidot = 0.0
                @simd for k in 1:numTFs
                    phidot += sol_new.u[k] * phis[i, k]
                end
                rpc = EPS / (1.0 + phidot)
                rho_per_cell[i, totalcells]    = rpc
                rho_cache_valid[i, totalcells] = true
                rhos[i] = nc * rpc
            else
                rhos[i] = nc * rho_per_cell[i, totalcells]
            end
        end

        # ---- Event rates ----------------------------------------------------
        division_rate  = sum(rhos)
        celldeath_rate = DELTAI * totalcells
        mutation_rate  = effectivem * totalcells
        total_rate     = division_rate + celldeath_rate + deltaO + mutation_rate

        total_rate <= 0 && break

        # ---- Advance time and record snapshots ------------------------------
        dt   = randexp() / total_rate
        t_end = time + dt

        # Fill every integer checkpoint the step [time, t_end) straddles.
        # totalcells hasn't changed yet, so this correctly reflects the
        # population *just before* the upcoming event fires.
        while next_record_time <= t_end && next_record_time <= MAX_TIME
            idx = trunc(Int, next_record_time)
            if 1 <= idx <= max_record_steps
                totalcells_history[idx] = totalcells
            end
            next_record_time += 1.0
        end

        time    = t_end
        fitness += WofN(totalcells) * dt

        # ---- Sample event ---------------------------------------------------
        r = rand() * total_rate
        event = if r < division_rate
            1
        elseif r < division_rate + celldeath_rate
            2
        elseif r < division_rate + celldeath_rate + deltaO
            3
        else
            4
        end

        # ---- Execute event --------------------------------------------------
        if event == 1   # Cell division
            dividingcell = sample(1:NUM_GENOTYPES, Weights(rhos))
            numcells[dividingcell] += 1
            totalcells += 1
            if totalcells > NUM_GENOTYPES
                idx = trunc(Int, next_record_time)
                totalcells_history[idx:end] .= totalcells
                break
            end

        elseif event == 2   # Cell death
            if totalcells > 1
                deadcell = sample(1:NUM_GENOTYPES, Weights(numcells))
                numcells[deadcell] -= 1
                totalcells -= 1
            end

        elseif event == 3   # Organism death
            idx = trunc(Int, next_record_time)
            totalcells_history[idx:end] .= totalcells
            break

        elseif event == 4   # Mutation
            mutatingcell = sample(1:NUM_GENOTYPES, Weights(numcells))
            numcells[mutatingcell] -= 1

            mutatedcell = findfirst(==(0), numcells)
            isnothing(mutatedcell) && break

            ks[mutatedcell, :, :]  .= ks[mutatingcell, :, :]
            phis[mutatedcell, :]   .= phis[mutatingcell, :]
            rs[mutatedcell, :]     .= rs[mutatingcell, :]
            states[mutatedcell, :] .= states[mutatingcell, :]

            whichmutationtype = sample(mut_types, mut_weights)

            if whichmutationtype == "MUT_PHI"
                idx = rand(nonzero_phi_idx)
                phis[mutatedcell, idx] = clamp(
                    phis[mutatedcell, idx] * rand(LogNormal(0, mutmag)), 1e-30, 1e30)

            elseif whichmutationtype == "MUT_K"
                idx = rand(nonzero_k_idx)
                ks[mutatedcell, idx] = clamp(
                    ks[mutatedcell, idx] * rand(LogNormal(0, mutmag)), 1e-30, 1e30)

            elseif whichmutationtype == "MUT_R"
                idx = rand(nonzero_r_idx)
                rs[mutatedcell, idx] = clamp(
                    rs[mutatedcell, idx] * rand(LogNormal(0, mutmag)), 1e-30, 1e30)
            end

            rho_cache_valid[mutatedcell, :] .= false
            numcells[mutatedcell] = 1
        end
    end

    return fitness, totalcells_history
end

# ---- Data -------------------------------------------------------------------
Nvst_data = function(KCmatrix, rs, phis;
                     num_sims = 100, MAX_TIME = 1000, NUM_GENOTYPES = 50,
                     time_bin = 10, log_floor = -3.0)
 
    max_time = trunc(Int, MAX_TIME)
 
    # ---- Run simulations ----------------------------------------------------
    results = Vector{Vector{Int}}(undef, num_sims)
    Threads.@threads for i in 1:num_sims
        _, history = Nvst_sim(KCmatrix, rs, phis;
                              m = m, mutmag = mutmag,       # globals, as before
                              MAX_TIME = MAX_TIME, NUM_GENOTYPES)
        results[i] = history
    end
 
    # ---- Binned cell-count distribution -------------------------------------
    cell_bins = 0:(NUM_GENOTYPES + 2)
    n_tbins   = cld(max_time, time_bin)
    time_axis = ((1:n_tbins) .- 0.5) .* time_bin            # bin centers
 
    freq = zeros(Int, length(cell_bins), n_tbins)
    for history in results, (t, n) in enumerate(history)
        freq[n + 1, cld(t, time_bin)] += 1
    end
 
    prob = freq ./ max.(sum(freq, dims = 1), 1)
 
    # Zero-probability cells become NaN so the heatmap leaves them transparent
    log_prob = ifelse.(prob .> 0, clamp.(log10.(prob), log_floor, 0.0), NaN)
    log_prob = vcat(log_prob, fill(NaN, 1, n_tbins))        # padding row (y = 0:53)
 
    # ---- Mean cells, fitness change, cumulative fitness (full resolution) ---
    cells_mat = reduce(hcat, results)
    dts_mat   = reduce(hcat, [WofN.(r) for r in results])
 
    mean_cells       = vec(mean(cells_mat, dims = 2))
    dtmeanvst        = vec(mean(dts_mat, dims = 2))
    meancumfitnesses = vec(mean(cumsum(dts_mat, dims = 1), dims = 2))
 
    return (; results, cell_bins, max_time, time_bin, time_axis, log_prob,
              log_floor, mean_cells, dtmeanvst, meancumfitnesses)
end

bin_mean(x, b) = [mean(@view x[i:min(i + b - 1, end)]) for i in 1:b:length(x)]

# Trajectory lines colored by how many other trajectories sit nearby.
# Each *segment* (between two consecutive time bins) gets its own color from the
# local density at its endpoints, combined with `segment_agg`. The default, `min`,
# means a segment only looks hot if BOTH ends are in a crowded region, so a line
# leaving (or entering) a hot zone no longer gets a hot connecting segment.
Nvst_hot_trajectories = function(data;
                                 time_bin    = data.time_bin,  # 1 = no binning
                                 bandwidth   = 0.75,           # "nearby" in cell-count units
                                 log_scale   = false,           # log10 of the local count
                                 segment_agg = min,            # min | (a, b) -> (a + b) / 2 | max
                                 nlevels     = 64,             # number of discrete colors
                                 linewidth   = 1.0,
                                 alpha       = 1.0,
                                 colorbar    = true,
                                 bg          = :black,
                                 xlim        = (0, 5000),
                                 ylim        = (0, NUM_GENOTYPES+2),
                                 linewidthbottomplot = 3,
                                 title       = :none)

    results = data.results
    n = length(results)

    Y  = reduce(hcat, [bin_mean(r, time_bin) for r in results])  # time bins x sims
    nb = size(Y, 1)
    x  = ((1:nb) .- 0.5) .* time_bin

    # ---- Local overlap at each (time bin, trajectory) -----------------------
    D = zeros(nb, n)
    for k in 1:nb
        yk = @view Y[k, :]
        for i in 1:n
            D[k, i] = sum(exp.(-0.5 .* ((yk .- yk[i]) ./ bandwidth) .^ 2))  # includes self (=1)
        end
    end
    log_scale && (D = log10.(D))

    # ---- One value per segment, then map to a color level -------------------
    S  = segment_agg.(D[1:end-1, :], D[2:end, :])               # (nb-1) x sims
    hi = max(maximum(D), 1e-6)
    lo = -0.05 * hi                                             # keeps coolest lines visible
    lvl = clamp.(round.(Int, 1 .+ (S .- lo) ./ (hi - lo) .* (nlevels - 1)), 1, nlevels)

    cg = Plots.cgrad(:viridis)
    p = Plots.plot(;
        background_color_inside = bg,
        xlim = xlim, ylim = ylim,
        xlabel = "Time", ylabel = "Total cells",
        legend = false,
        colorbar = colorbar,
        colorbar_title = log_scale ? "log10(nearby trajectories)" : "Nearby trajectories",
        title = title, titlefontsize = 7,
        extra_plot_kwargs = Any[:tight_layout => true],
    )

    # ---- Draw coolest level first, hottest last (hot segments on top) -------
    for l in 1:nlevels
        xs = Float64[]; ys = Float64[]
        for i in 1:n, k in 1:(nb - 1)
            lvl[k, i] == l || continue
            push!(xs, x[k], x[k+1], NaN)
            push!(ys, Y[k, i], Y[k+1, i], NaN)
        end
        isempty(xs) && continue
        Plots.plot!(p, xs, ys;
                    color = cg[(l - 1) / (nlevels - 1)],
                    linewidth = linewidth, alpha = alpha, label = "")
    end

    # ---- Off-screen point whose only job is to supply the colorbar ----------
    if colorbar
        Plots.scatter!(p, [xlim[1] - 1e6], [ylim[1] - 1e6];
                       marker_z = [hi], clims = (lo, hi), color = :viridis,
                       markersize = 0, markerstrokewidth = 0, label = "")
    end

    p2 = Plots.plot(data.time_axis, data.dtmeanvst;
        linewidth = linewidthbottomplot, label = "Mean", legend = false, color = :black,
        ylim = (0, 1.0), xlim = xlim, ylabel = "Reproductive rate", xlabel = "Time", framestyle = :box, extra_plot_kwargs = Any[:tight_layout => true])

    my_layout = @layout [
        a{0.6h}
        b{0.4h}
    ]
    p5 = Plots.plot(p, p2; layout = my_layout, link = :x,
                    size = (400, 400), right_margin = 5Plots.mm)

    return p5
end

# HDF5.jl renamed `attributes` to `attrs` in newer versions; support both.
_attrs(x) = isdefined(HDF5, :attrs) ? HDF5.attrs(x) : HDF5.attributes(x)

# Set an attribute, skipping empty strings (zero-length strings are awkward in HDF5)
function _setattr(obj, key, val)
    (val isa AbstractString && isempty(val)) && return nothing
    _attrs(obj)[string(key)] = val
    return nothing
end

# Save everything needed to replot AND document the file for reviewers.
# `data`   : NamedTuple returned by Nvst_data
# `params` : simulation inputs worth archiving, e.g.
#            (; KCmatrix, rs, phis, m, mutmag, num_sims, MAX_TIME, NUM_GENOTYPES)
function Nvst_save(path, data;
                   params       = (;),
                   seed         = nothing,
                   title        = "",
                   description  = "",
                   figure       = "",
                   code_version = "",
                   contact      = "",
                   license      = "",
                   time_units   = "simulation time steps")

    n_sims      = length(data.results)
    results_mat = permutedims(reduce(hcat, data.results))      # sims x time
    n_steps     = size(results_mat, 2)

    # (name, value, description) for every dataset
    datasets = [
        ("results", results_mat,
         "Total number of cells at each time step for each simulation replicate. " *
         "Axes as written from Julia: (replicate, time_step). NOTE: h5py/numpy report " *
         "the reversed shape (time_step, replicate); MATLAB follows Julia's order."),
        ("mean_cells", data.mean_cells,
         "Mean total cell count across replicates at each time step. Length = number of time steps."),
        ("dtmeanvst", data.dtmeanvst,
         "Mean across replicates of WofN(total cells), the per-time-step change in fitness, " *
         "at each time step (plotted as dW/dt). Length = number of time steps."),
        ("meancumfitnesses", data.meancumfitnesses,
         "Mean across replicates of the cumulative sum over time of WofN(total cells), " *
         "i.e. cumulative fitness. Length = number of time steps."),
        ("time_bin", data.time_bin,
         "Default bin width (in time steps) used when plotting. Does not affect the raw data."),
        ("max_time", data.max_time,
         "Number of time units simulated per replicate."),
        ("NUM_GENOTYPES", last(data.cell_bins) - 2,
         "NUM_GENOTYPES setting passed to the simulation."),
        ("log_floor", get(data, :log_floor, -3.0),
         "Lower clamp on log10(probability) used for the plotted cell-count distribution."),
    ]

    h5open(path, "w") do f
        for (name, val, desc) in datasets
            f[name] = val
            _setattr(f[name], "description", desc)
        end
        _setattr(f["results"], "units", "cells")
        _setattr(f["results"], "time_units", time_units)

        # ---- Simulation inputs ----------------------------------------------
        g = create_group(f, "parameters")
        _setattr(g, "description",
                 "Inputs passed to the simulation that generated /results.")
        for (k, v) in pairs(params)
            v = v isa Bool          ? Int(v)    :
                v isa AbstractRange ? collect(v) : v
            try
                g[string(k)] = v
            catch err
                @warn "Could not store parameter $k of type $(typeof(v)); skipped" err
            end
        end

        # ---- File-level documentation ---------------------------------------
        _setattr(f, "title", title)
        _setattr(f, "description", description)
        _setattr(f, "figure", figure)
        _setattr(f, "code_version", code_version)
        _setattr(f, "contact", contact)
        _setattr(f, "license", license)
        _setattr(f, "format_version", "1.0")
        _setattr(f, "created", string(Dates.now()))
        _setattr(f, "julia_version", string(VERSION))
        _setattr(f, "n_threads", Threads.nthreads())
        _setattr(f, "n_replicates", n_sims)
        _setattr(f, "n_time_steps", n_steps)
        _setattr(f, "time_units", time_units)
        _setattr(f, "random_seed",
                 isnothing(seed) ?
                 "not set; replicates ran in parallel threads, so exact reruns are not guaranteed" :
                 string(seed))
    end
    return path
end

# Load a saved file back into the same NamedTuple Nvst_data returns, so it works
# with Nvst_plot, Nvst_hot_trajectories, etc. Pass time_bin (or log_floor) to
# re-bin differently without rerunning the simulations.
function Nvst_load(path; time_bin = nothing, log_floor = nothing)
    h = h5open(path, "r") do f
        (results          = read(f, "results"),
         mean_cells       = read(f, "mean_cells"),
         dtmeanvst        = read(f, "dtmeanvst"),
         meancumfitnesses = read(f, "meancumfitnesses"),
         time_bin         = read(f, "time_bin"),
         max_time         = read(f, "max_time"),
         NUM_GENOTYPES    = read(f, "NUM_GENOTYPES"),
         log_floor        = read(f, "log_floor"))
    end

    tb = isnothing(time_bin)  ? h.time_bin  : time_bin
    lf = isnothing(log_floor) ? h.log_floor : log_floor

    results  = [Vector{Int}(h.results[i, :]) for i in 1:size(h.results, 1)]
    max_time = h.max_time

    # ---- Rebuild the binned distribution (same as Nvst_data) ----------------
    cell_bins = 0:(h.NUM_GENOTYPES + 2)
    n_tbins   = cld(max_time, tb)
    time_axis = ((1:n_tbins) .- 0.5) .* tb

    freq = zeros(Int, length(cell_bins), n_tbins)
    for history in results, (t, n) in enumerate(history)
        freq[n + 1, cld(t, tb)] += 1
    end
    prob     = freq ./ max.(sum(freq, dims = 1), 1)
    log_prob = ifelse.(prob .> 0, clamp.(log10.(prob), lf, 0.0), NaN)
    log_prob = vcat(log_prob, fill(NaN, 1, n_tbins))

    return (; results, cell_bins, max_time, time_bin = tb, time_axis, log_prob,
              log_floor = lf, mean_cells = h.mean_cells,
              dtmeanvst = h.dtmeanvst, meancumfitnesses = h.meancumfitnesses)
end

function plot_weighted_network_makie(A, win, wout; figwidth = 100, rowwidth = 0.1, edgecurve = 0.05, curves = true, interior_node_labeling = false)
    #A = transpose(A)
    n = size(A, 1)
    input_node  = n + 1
    output_node = n + 2

    # ---- Build expanded graph ----
    g = DiGraph(n + 2)
    for i in 1:n, j in 1:n
        A[i, j] != 0 && add_edge!(g, i, j)
    end
    for i in 1:n
        win[i]  != 0 && add_edge!(g, input_node, i)
        wout[i] != 0 && add_edge!(g, i, output_node)
    end

    # ---- Layered layout ----
    has_in  = [win[i]  != 0 for i in 1:n]
    has_out = [wout[i] != 0 for i in 1:n]

    input_only  = [i for i in 1:n if  has_in[i] && !has_out[i]]
    output_only = [i for i in 1:n if !has_in[i] &&  has_out[i]]
    middle      = [i for i in 1:n if (has_in[i] && has_out[i]) || (!has_in[i] && !has_out[i])]

    #function row_positions(nodes, y)
    #    m = length(nodes)
    #    m == 0 && return Point2f[]
    #    xs = m == 1 ? [0.0] : collect(range(-1.0, 1.0; length=m))
    #    [Point2f(xs[k], y) for k in 1:m]
    #end

    function row_positions(nodes, y)
        m = length(nodes)
        m == 0 && return Point2f[]
        xs = m == 1 ? [0.0] : collect(range(-rowwidth, rowwidth; length=m))   # was -1.0, 1.0
        [Point2f(xs[k], y) for k in 1:m]
    end

    y_input_node  =  0.5
    y_input_row   =  0.2    # was 1.0 — top interior layer pulled toward center
    y_middle_row  =  0.0
    y_output_row  = -0.2    # was -1.0 — bottom interior layer pulled toward center
    y_output_node = -0.5

    positions = Vector{Point2f}(undef, n + 2)
    for (nodes, y) in [(input_only, y_input_row), (middle, y_middle_row), (output_only, y_output_row)]
        for (k, i) in enumerate(nodes)
            positions[i] = row_positions(nodes, y)[k]
        end
    end
    positions[input_node]  = Point2f(0.0, y_input_node)
    positions[output_node] = Point2f(0.0, y_output_node)

    # ---- Node styling ----
    nodecolors = Vector{RGBf}([RGBf(1,1,1) for _ in 1:n])
    #push!(nodecolors, RGBf(0.35, 0.75, 0.35))
    #push!(nodecolors, RGBf(0.90, 0.50, 0.20))
    push!(nodecolors, RGBf(1,1,1))
    push!(nodecolors, RGBf(1,1,1))
    nodesizes  = Vector{Float64}([40.0 for _ in 1:n])
    push!(nodesizes, 60.0); push!(nodesizes, 40.0)

    # No outline on the in/out nodes; keep it on the white interior nodes.
    node_strokewidths = Vector{Float64}([1.0 for _ in 1:n])
    push!(node_strokewidths, 0.0); push!(node_strokewidths, 0.0)

    # Labels: only the input and output terminal nodes (interior nodes stay unlabeled).
    node_labels = fill("", n + 2)
    nlabels_fontsizes = fill(25, n + 2)
    node_labels[input_node]  = "N"
    node_labels[output_node] = "λ"

    if interior_node_labeling == true
        counter = 1
        for labelindex in 1:length(node_labels)
            if node_labels[labelindex] == ""
                node_labels[labelindex] = "T"*string(counter)
                #node_labels[labelindex] = latexstring("Variable_{$counter}")
                #counterstring = string(counter)
                #node_labels[labelindex] = L"Variable_{%$(counterstring)}"
                nlabels_fontsizes[labelindex] = 18
                counter += 1
            end
        end
    end



    # ---- Edge styling ----
    # ---- Edge styling ----
    # ---- Edge styling ----
    edge_list = collect(edges(g))
    nedges    = length(edge_list)

    interior_ew = [log10(A[src(e), dst(e)]) for e in edge_list
                   if src(e) <= n && dst(e) <= n]
    maxew   = isempty(interior_ew)    ? 1.0 : maximum(abs, interior_ew)
    maxwin  = maximum(abs, win)  == 0 ? 1.0 : maximum(abs, win)
    maxwout = maximum(abs, wout) == 0 ? 1.0 : maximum(abs, wout)

    nw_vals     = Vector{Float64}(undef, nedges)
    grp_vals    = Vector{Symbol}(undef, nedges)
    edge_colors = Vector{RGBf}(undef, nedges)
    for (k, e) in enumerate(edge_list)
        s, d = src(e), dst(e)
        if s <= n && d <= n
            nw = abs(log2(A[s, d])) / maxew
            grp_vals[k] = :interior          # interior grey set in the loop below
        elseif s == input_node
            nw = abs(win[d]) / maxwin
            edge_colors[k] = RGBf(0.2, 0.55 + 0.35*nw, 0.2)
            grp_vals[k] = :in
        else
            nw = abs(wout[s]) / maxwout
            edge_colors[k] = RGBf(0.8, 0.3 + 0.3*nw, 0.1)
            grp_vals[k] = :out
        end
        nw_vals[k] = nw
    end

    # Width (all groups) + interior grey, both from a per-group min–max stretch.
    min_width, max_width = 3.0, 6.0
    edge_widths = Vector{Float64}(undef, nedges)
    for grp in (:interior, :in, :out)
        idx = findall(==(grp), grp_vals)
        isempty(idx) && continue
        lo, hi = extrema(nw_vals[idx])
        for k in idx
            t = hi > lo ? (nw_vals[k] - lo) / (hi - lo) : 0.5
            #edge_widths[k] = min_width + (max_width - min_width) * t
            edge_widths[k] = min_width + (max_width - min_width) * (1 - t)
            if grp == :interior
                # lower weight ⇒ higher nw ⇒ higher t ⇒ lighter grey;
                # higher weight ⇒ toward black. (0.75 keeps the weakest visible.)
                v = 0.75 * t
                edge_colors[k] = RGBf(v, v, v)
            end
        end
    end

    # Curve interior edges only; keep the input/output edges straight.
    edge_curves = [grp_vals[k] == :interior ? edgecurve : 0.0 for k in 1:nedges]

    # ---- Plot ----
    fig = Figure(size = (figwidth, 400), backgroundcolor = :transparent)   # was Figure() — first number is the width
    #ax = Axis(fig[1, 1]; aspect = DataAspect())
    # ...and in row_positions, narrow the spread:
    #xs = m == 1 ? [0.0] : collect(range(-0.6, 0.6; length=m))   # was -1.0, 1.0
    
    ax  = Axis(fig[1, 1], backgroundcolor = :transparent)
    hidedecorations!(ax)
    hidespines!(ax)

    CairoMakie.xlims!(ax, -1.0, 1.0)
    CairoMakie.ylims!(ax, -0.7, 0.7)

    pos_ref   = Ref(positions)
    layout_fn = (_) -> pos_ref[]

    GraphMakie.graphplot!(
        ax, g;
        layout           = layout_fn,
        node_color       = nodecolors,
        node_size        = nodesizes,
        node_strokecolor = :black,
        node_strokewidth = node_strokewidths,   # per-node vector (0 = no outline)
        edge_color       = edge_colors,
        edge_width       = edge_widths,          # back to the varying vector
        arrow_marker     = :vline,               # blunt T-bar head
        arrow_size       = 20,
        arrow_shift      = :end,
        curve_distance       = edge_curves,
        curve_distance_usage = curves,
        nlabels          = node_labels,
        nlabels_align    = (:center, :center),   # centered inside each terminal node
        nlabels_color    = :black,
        #nlabels_fontsize = 25,
        nlabels_fontsize = nlabels_fontsizes,
    )

    @show extrema(edge_widths)

    fig
end

checkerboard8 = function(fold)
    size = fold*2
    KCmatrix = zeros(size,size)
    #for i in 1:(size^2)
    #    KCmatrix[i] = 1 - i%2
    #end
    @inbounds for i in 1:size, j in 1:size
        if i%2 == 0
            mult = 1#3.04339#3.0506#10#2
        else
            mult = 0.3286#1
        end
        KCmatrix[i, j] = isodd(i + j) ? 1*mult : 0
    end
    #KCmatrix[1,size] = 1
    #KCmatrix[size, 1] = 1

    rs = zeros(size)
    phis = zeros(size)
    numsignals = 0
    for i in 1:size
        rs[i] = i%2
        phis[i] = 1 - i%2
        numsignals += 1 - i%2
    end
    rs .*= 0.1744#0.17522#9/12#1/12
    phis .*= 23135.7#/(numsignals) 

    KCmatrix .*= 1e8/numsignals
    return KCmatrix, rs, phis
end

hybridnnodesnedges = function(hybridvect)
    #formula for the number of edges in an ssm averaging network is 2*(i^2) + i)
    #for a min network is 3*i
    #where i is the number of ssm modules
    totalnnodes = 0
    totalnnedges = 0
    for i in 1:length(hybridvect)
        totalnnodes += hybridvect[i]*i*2
        totalnnedges += hybridvect[i]*(2*(i^2) + i)
    end

    return totalnnodes, totalnnedges
end

function duplicate_nodes(A, rs, phis, nodes_to_duplicate; connect_duplicates = false)
    N = size(A, 1)
    @assert size(A, 2) == N "Adjacency matrix must be square."
    @assert all(1 .<= nodes_to_duplicate .<= N) "Node indices out of bounds."
    
    num_dups = length(nodes_to_duplicate)
    new_N = N + num_dups
    
    # Pre-allocate the new larger adjacency matrix with the same type
    new_A = zeros(eltype(A), new_N, new_N)
    
    # 1. Copy the original graph into the top-left quadrant
    new_A[1:N, 1:N] = A
    
    # 2. Duplicate rows and columns for each target node
    for (i, original_node) in enumerate(nodes_to_duplicate)
        duplicate_node = N + i
        
        # Copy outgoing edges (Row)
        new_A[duplicate_node, 1:N] = A[original_node, :]
        
        # Copy incoming edges (Column)
        new_A[1:N, duplicate_node] = A[:, original_node]
        
        # Copy cross-connections between already duplicated nodes if applicable
        for (j, other_original) in enumerate(nodes_to_duplicate[1:i-1])
            other_duplicate = N + j
            new_A[duplicate_node, other_duplicate] = A[original_node, other_original]
            new_A[other_duplicate, duplicate_node] = A[other_original, original_node]
        end
        
        # 3. Optional: Connect the original node and its twin
        if connect_duplicates
            new_A[original_node, duplicate_node] = one(eltype(A))
            new_A[duplicate_node, original_node] = one(eltype(A)) # For undirected consistency
        end
    end

    rs_new = zeros(length(rs) + length(nodes_to_duplicate))
    rs_new[1:length(rs)] .= rs
    phis_new = zeros(length(rs_new))
    phis_new[1:length(phis)] .= phis
    

    for i in 1:length(nodes_to_duplicate)
        rs_new[length(rs)+i] = rs[nodes_to_duplicate[i]]
        phis_new[length(rs)+i] = phis[nodes_to_duplicate[i]]
    end
    
    return new_A, rs_new,phis_new
end

slaptogether = function(KCmatrix1, rs1, phis1, KCmatrix2, rs2, phis2)

    KCmatrix = zeros(length(rs1) + length(rs2), length(rs1) + length(rs2))

    KCmatrix[1:length(rs1),1:length(rs1)] .= KCmatrix1
    KCmatrix[(length(rs1)+1):end,(length(rs1)+1):end] .= KCmatrix2

    return KCmatrix, vcat(rs1, rs2), vcat(phis1, phis2)
end

makespecifichybrid = function(vect)
    #the vector describes the number of averaging motifs which average according to it's index.
    #where [N,0,0,0,0,0] is a pure min network and [0,0,0,0,0,1] is a pure averaging network.

    #we can construct this by using checkerboard7 and writing a function to slap networks together.

    KCmatrix = zeros(1,1)
    rs = zeros(1)
    phis = zeros(1)
    for k in 1:length(vect)

        if vect[k] > 0

            for i in 1:vect[k]
                subKCmatrix,subrs,subphis = checkerboard8(k)
                KCmatrix,rs,phis = slaptogether(KCmatrix,rs,phis, subKCmatrix, subrs, subphis)
            end
        end
    end

    KCmatrix = KCmatrix[2:end, 2:end]
    rs = rs[2:end]
    phis = phis[2:end]

    numsignals = sum(rs .> 0)
    phis = phis/numsignals

    return KCmatrix,rs,phis


end

plotmanymutants_withTFs = function(KCmatrix, rs, phis, poissonl, nsamples; mutmag = mutmag)

    allparams = vcat(reshape(rs, 1, length(phis)), KCmatrix)
    nonzeroparams = findall(!iszero, allparams)

    higherorloweratroot = zeros(nsamples)

    #first we just do a single mutation
    lambdas, sols = solvegrowthfunction_withTFs(50, KCmatrix,rs,phis)
    m1lambdas, m1sols = solvegrowthfunction_withTFs(50, KCmatrix, rs/2.0, phis)

    p1 = Plots.plot(collect(1:50), sols[:,1], size = (400, 200), linewidth = 5, label = "before mutation", color = :blue, linestyle = :dot, ylabel = "T1 concentration", xlabel = "Number of cells (N)")
    p1 = Plots.plot!(p1, collect(1:50), m1sols[:,1], linewidth = 3, linestyle = :solid, color = :red, label = "after mutation")

    p2 = Plots.plot(collect(1:50), sols[:,2], size = (400, 200),ylabel = "T2 concentration", xlabel = "Number of cells (N)", linewidth = 5, label = "before mutation", color = :blue, linestyle = :dot)
    p2 = Plots.plot!(p2, collect(1:50), m1sols[:,2], linewidth = 3, linestyle = :solid, color = :red, label = "after mutation")

    for i in 1:nsamples
        allparams = vcat(reshape(rs, 1, length(phis)), KCmatrix)

        nmutations = rand(Poisson(poissonl))
        for mut in 1:nmutations
            mutparam = rand(nonzeroparams)
            allparams[mutparam] *= rand(LogNormal(0, mutmag))
        end
        rsmut = allparams[1,:]
        KCmatrixmut = allparams[2:end,:]
        mlambdas = solvegrowthfunction_withTFs(50, KCmatrixmut, rsmut, phis)[1]
        if mlambdas[12] == lambdas[12]
            higherorloweratroot[i] = 0.5
        else
            higherorloweratroot[i] = mlambdas[12] > lambdas[12]
        end
        if i == 1
            plot = Plots.plot(mlambdas, size = (400, 200), linewidth = 0.8, linestyle = :solid, color = :grey, label = "")
        else
            #Plots.plot!(plot, collect(1:50) .+ (rand(50)*0.5 .* rand([-1,1])), mlambdas, linewidth = 1, linestyle = :solid, color = i+1, label = false)
            #Plots.plot!(plot, collect(1:50), mlambdas, linewidth = 0.8, linestyle = :solid, color = i+1, label = false)
            Plots.plot!(plot, collect(1:50), mlambdas, linewidth = 0.8, linestyle = :solid, color = :grey, label = false)
 
        end

    end
    plot = Plots.plot!(plot, collect(1:50), m1lambdas, linewidth = 3, linestyle = :solid, color = :red, label = "")
 

    plot = Plots.plot!(plot, lambdas, ylim = (-0.0025,EPS+ 0.0025), xlabel = "Number of cells (N)", ylabel = "Cell replication rate (λ)", linewidth = 5, label = "", color = :blue, linestyle = :dot, legend = true)
    #vline!([12])
    #hline!([0.1])

    plot = Plots.plot!(plot, [NaN], [NaN], linewidth = 5, label = "before mutation", color = :blue, linestyle = :dot, legend = true)
 
    plot = Plots.plot!(plot,[NaN], [NaN], linewidth = 3, linestyle = :solid, color = :red, label = "after mutation")
 
    plot = Plots.plot!([NaN], [NaN], size = (400, 200), linewidth = 0.8, linestyle = :solid, color = :grey, label = "other mutants")


    p5 = Plots.plot(p1,p2,plot, layout = (3,1), size = (400,600), left_margin = 12Plots.mm, link = :x)
    #display(p5)
    print(mean(higherorloweratroot))
    return p5
end

solvegrowthfunction_withTFs = function(maxN, KCmatrix, rs, phis)

    rhos = zeros(maxN)

    numTFs = trunc(Int, sqrt(length(KCmatrix)))

    N = 1

    p = reshape(KCmatrix, numTFs, numTFs), 1, 1, N, rs
    u0 = ones(numTFs)*0.1
    prob = SteadyStateProblem(GRNode_signal!, u0, p)
    sol = solve(prob, DynamicSS(Rosenbrock23()), abstol=1e-10)
    rhos[1] = EPS/(1+sum(sol .* phis))
    #rho = EPS*mean(sol .*phis)

    u0 = sol

    sols = zeros(maxN, numTFs)
    sols[1,:] .= sol

    for N in 2:maxN
        p = reshape(KCmatrix, numTFs, numTFs), 1, 1, N, rs
        #prob = SteadyStateProblem(GRNode_signal!, u0, p)
        prob2 = remake(prob, u0 = u0, p=p)
        #sol = solve(prob, DynamicSS(Rosenbrock23()), abstol=1e-10)
        sol = solve(prob2, DynamicSS(Rodas5P()); abstol=1e-10)
        rhos[N] = EPS/(1+sum(sol .* phis))
        #rho = EPS*mean(sol .*phis)
        #u0 = sol .+ 0.001
        sols[N,:] .= sol
    end

    return rhos, sols
end

plotmanymutants = function(KCmatrix, rs, phis, poissonl, nsamples; mutmag = mutmag, NUM_GENOTYPES = 50, size = (400,300))
    lambdas = solvegrowthfunction_withTFs(NUM_GENOTYPES, KCmatrix, rs, phis)[1]

    allparams = vcat(reshape(rs, 1, length(phis)), KCmatrix)
    nonzeroparams = findall(!iszero, allparams)

    higherorloweratroot = zeros(nsamples)

    for i in 1:nsamples
        allparams = vcat(reshape(rs, 1, length(phis)), KCmatrix)

        nmutations = rand(Poisson(poissonl))
        for mut in 1:nmutations
            mutparam = rand(nonzeroparams)
            allparams[mutparam] *= rand(LogNormal(0, mutmag))
        end
        rsmut = allparams[1,:]
        KCmatrixmut = allparams[2:end,:]
        mlambdas = solvegrowthfunction_withTFs(NUM_GENOTYPES, KCmatrixmut, rsmut, phis)[1]
        if mlambdas[12] == lambdas[12]
            higherorloweratroot[i] = 0.5
        else
            higherorloweratroot[i] = mlambdas[12] > lambdas[12]
        end
        if i == 1
            plot = Plots.plot(mlambdas, size = size, linewidth = 0.5, linestyle = :solid, color = :black, label = "")
        else
            #Plots.plot!(plot, collect(1:NUM_GENOTYPES) .+ (rand(NUM_GENOTYPES)*0.5 .* rand([-1,1])), mlambdas, linewidth = 0.8, linestyle = :solid, color = i+1, label = false)
            Plots.plot!(plot, collect(1:NUM_GENOTYPES) .+ (rand(NUM_GENOTYPES)*0.5 .* rand([-1,1])), mlambdas, linewidth = 0.5, linestyle = :solid, color = :black, label = false)
        end

    end

    plot = Plots.plot!(lambdas, ylim = (-0.0025,EPS+ 0.0025), ylabelfontsize = 8, xlabel = "Number of cells (N)", ylabel = "Cell replication rate (λ)", linewidth = 5, label = "", color = :blue, linestyle = :dot, legend = true)
    #vline!([12])
    #hline!([0.1])

    plot = Plots.plot!(plot, [NaN], [NaN], linewidth = 5, label = "before mutation", color = :blue, linestyle = :dot, legend = true)
 
    plot = Plots.plot!([NaN], [NaN], size = size, linewidth = 1.0, linestyle = :solid, color = :black, label = "sampled mutants")


    #display(plot)
    print(mean(higherorloweratroot))
    return plot
end

##plotting analytic data

function color_simplex(n::Int=240)
    v1 = (0.5, sqrt(3)/2)   # strategy 1
    v2 = (1.0, 0.0)         # strategy 2
    v3 = (0.0, 0.0)         # strategy 3
    denom = (v2[2]-v3[2])*(v1[1]-v3[1]) + (v3[1]-v2[1])*(v1[2]-v3[2])

    xs = range(0, 1, length=n)
    ys = range(0, sqrt(3)/2, length=n)
    img = fill(RGBA(1.0, 1.0, 1.0, 0.0), n, n)   # transparent background

    for j in 1:n, i in 1:n
        x, y = xs[j], ys[i]
        p1 = ((v2[2]-v3[2])*(x-v3[1]) + (v3[1]-v2[1])*(y-v3[2])) / denom
        p2 = ((v3[2]-v1[2])*(x-v3[1]) + (v1[1]-v3[1])*(y-v3[2])) / denom
        p3 = 1 - p1 - p2
        if p1 >= 0 && p2 >= 0 && p3 >= 0
            c = ryb_blend(p1, p2, p3)
            img[i,j] = RGBA(red(c), green(c), blue(c), 1.0)
        end
    end
    xs, ys, img
end

function plot_simplex_legend(; labels=("Strategy 1", "Strategy 2", "Strategy 3"),
                               n::Int=240)
    xs, ys, img = color_simplex(n)
    p = Plots.plot(xs, ys, img;
             aspect_ratio=:equal, framestyle=:none, legend=false,
             ticks=nothing, xlims=(-0.25, 1.25), ylims=(-0.15, 1.05))
    # Vertex labels
    annotate!(p, 0.5,         sqrt(3)/2 + 0.1, Plots.text(labels[1], 12, :center, :bottom))
    annotate!(p, 1.0 + 0.02, -0.04,              Plots.text(labels[2], 12, :left,   :top))
    annotate!(p, 0.0 - 0.02, -0.04,              Plots.text(labels[3], 12, :right,  :top))
    # Edge / mix labels (optional, comment out if too busy)
    
    #annotate!(p, 0.98, sqrt(3)/4 + 0.04, Plots.text("FF Min+SSM Avg", 9, :center, :grey20))
    annotate!(p, 0.98, sqrt(3)/4 + 0.04, Plots.text(labels[1]*"+"*labels[2], 9, :center, :grey20))
    #annotate!(p, 0.00, sqrt(3)/4 + 0.04, Plots.text("FF Min+SSM Min", 9, :center, :grey20))
    annotate!(p, 0.00, sqrt(3)/4 + 0.04, Plots.text(labels[1]*"+"*labels[3], 9, :center, :grey20))
    #annotate!(p, 0.50, -0.06,             Plots.text("SSM Avg+SSM Min", 9, :center, :grey20))
    annotate!(p, 0.50, -0.06, Plots.text(labels[2]*"+"*labels[3], 9, :center, :grey20))
    annotate!(p, 0.50,  sqrt(3)/6+0.05,        Plots.text("3-way tie", 11, :center, :black))
    return p
end

# RYB blend: p1,p2,p3 are probabilities (sum to 1) that
# strategies 1,2,3 are the fittest at this point.
function ryb_blend(p1, p2, p3; alpha = 2.4)
    # Pure-strategy anchors
    R = (0.90, 0.10, 0.10)   # strategy 1
    Y = (0.97, 0.90, 0.10)   # strategy 2
    B = (0.10, 0.25, 0.85)   # strategy 3
    # Pairwise-tie anchors (the whole point of using RYB)
    O = (0.98, 0.55, 0.05)   # 1+2 -> orange
    G = (0.15, 0.65, 0.20)   # 2+3 -> green
    P = (0.55, 0.10, 0.60)   # 1+3 -> purple (not requested, but needed)

    wR, wY, wB = p1^alpha, p2^alpha, p3^alpha
    wO, wG, wP = 2p1*p2, 2p2*p3, 2p1*p3   # sum of all weights == 1

    r = wR*R[1] + wY*Y[1] + wB*B[1] + wO*O[1] + wG*G[1] + wP*P[1]
    g = wR*R[2] + wY*Y[2] + wB*B[2] + wO*O[2] + wG*G[2] + wP*P[2]
    b = wR*R[3] + wY*Y[3] + wB*B[3] + wO*O[3] + wG*G[3] + wP*P[3]
    RGB(clamp(r,0,1), clamp(g,0,1), clamp(b,0,1))
end

# P1, P2, P3 :: AbstractMatrix, each entry the prob that strategy i wins
function dominance_image(P1, P2, P3)
    @assert size(P1) == size(P2) == size(P3)
    [ryb_blend(P1[i,j], P2[i,j], P3[i,j]) for i in axes(P1,1), j in axes(P1,2)]
end


function plot_dominance(title, P1, P2, P3; x=nothing, y=nothing, kwargs...)
    img = dominance_image(P1, P2, P3)
    #img = reverse(img, dims = 1)
    xs = x === nothing ? (1:size(P1,2)) : x
    ys = y === nothing ? (1:size(P1,1)) : y
    #Plots.plot(xs, ys, img; yflip=false, aspect_ratio=:equal,
    #     xlabel="x", ylabel="y", framestyle=:none, kwargs...)

    xs = collect(-13.0:0.075:-3.5)
    ys = collect(-13.0:0.075:-3.5)
    xt = vcat(collect(0.00001:0.00001:0.0001), collect(0.0001:0.0001:0.001), collect(0.001:0.001:0.01), collect(0.01:0.01:0.03))
    xlab = vcat(["0.00001"], repeat([""],9), ["0.0001"], repeat([""],9), ["0.001"], repeat([""],9), ["0.01"], repeat([""], 3))
    xt = (log.(xt), xlab)
    Plots.heatmap(xs, ys, img;
        xticks = xt, yticks = xt,
        xlims  = (-13, -3.5), ylims = (-13, -3.5),
        xrotation = 45,
        yrotation = 45,
        yflip = false,
        xlabel = "Organism Death Rate",
        ylabel = "Mutation Rate (/cell/edge)",
        xlabelfontsize = 8,
        ylabelfontsize = 8,
        xtickfont = (6, "sans-serif"),
        ytickfont = (6, "sans-serif"),
        #title = "N* = 12, sigma = 1.0",
        title = title,
        titlefontsize = 10,
        bottom_margin = 8mm,
        tick_direction = :out)
end