#plotting all ssm hybrids for supplemental
include("bones_ode.jl")

const TABLE = NTuple{6,Int}[(0,0,0,0,0,1), (0,0,0,0,1,0), (0,0,0,1,0,0), (0,0,0,2,0,0), (0,0,1,0,0,0),
    (0,0,1,0,1,0), (0,0,1,1,0,0), (0,0,2,0,0,0), (0,0,2,1,0,0), (0,0,3,0,0,0), (0,1,0,0,0,0),
    (0,1,0,0,1,0), (0,1,0,1,0,0), (0,1,1,0,0,0), (0,1,1,1,0,0), (0,1,2,0,0,0), (0,1,3,0,0,0),
    (0,2,0,0,0,0), (0,2,0,0,1,0), (0,2,0,1,0,0), (0,2,1,0,0,0), (0,2,1,1,0,0), (0,2,2,0,0,0),
    (0,3,0,0,0,0), (0,3,0,1,0,0), (0,3,1,0,0,0), (0,3,2,0,0,0), (0,4,0,0,0,0), (0,4,0,1,0,0),
    (0,4,1,0,0,0), (0,5,0,0,0,0), (0,5,1,0,0,0), (0,6,0,0,0,0), (0,7,0,0,0,0), (1,0,0,0,0,0),
    (1,0,0,0,1,0), (1,0,0,1,0,0), (1,0,0,2,0,0), (1,0,1,0,0,0), (1,0,1,0,1,0), (1,0,1,1,0,0),
    (1,0,2,0,0,0), (1,0,3,0,0,0), (1,1,0,0,0,0), (1,1,0,0,1,0), (1,1,0,1,0,0), (1,1,1,0,0,0),
    (1,1,1,1,0,0), (1,1,2,0,0,0), (1,1,3,0,0,0), (1,2,0,0,0,0), (1,2,0,0,1,0), (1,2,0,1,0,0),
    (1,2,1,0,0,0), (1,2,2,0,0,0), (1,3,0,0,0,0), (1,3,0,1,0,0), (1,3,1,0,0,0), (1,3,2,0,0,0),
    (1,4,0,0,0,0), (1,4,0,1,0,0), (1,4,1,0,0,0), (1,5,0,0,0,0), (1,5,1,0,0,0), (1,6,0,0,0,0),
    (1,7,0,0,0,0), (2,0,0,0,0,0), (2,0,0,0,1,0), (2,0,0,1,0,0), (2,0,0,2,0,0), (2,0,1,0,0,0),
    (2,0,1,1,0,0), (2,0,2,0,0,0), (2,0,3,0,0,0), (2,1,0,0,0,0), (2,1,0,0,1,0), (2,1,0,1,0,0),
    (2,1,1,0,0,0), (2,1,1,1,0,0), (2,1,2,0,0,0), (2,1,3,0,0,0), (2,2,0,0,0,0), (2,2,0,1,0,0),
    (2,2,1,0,0,0), (2,2,2,0,0,0), (2,3,0,0,0,0), (2,3,0,1,0,0), (2,3,1,0,0,0), (2,3,2,0,0,0),
    (2,4,0,0,0,0), (2,4,1,0,0,0), (2,5,0,0,0,0), (2,5,1,0,0,0), (2,6,0,0,0,0), (2,7,0,0,0,0),
    (3,0,0,0,0,0), (3,0,0,0,1,0), (3,0,0,1,0,0), (3,0,1,0,0,0), (3,0,1,1,0,0), (3,0,2,0,0,0),
    (3,0,3,0,0,0), (3,1,0,0,0,0), (3,1,0,0,1,0), (3,1,0,1,0,0), (3,1,1,0,0,0), (3,1,1,1,0,0),
    (3,1,2,0,0,0), (3,2,0,0,0,0), (3,2,0,1,0,0), (3,2,1,0,0,0), (3,2,2,0,0,0), (3,3,0,0,0,0),
    (3,3,0,1,0,0), (3,3,1,0,0,0), (3,4,0,0,0,0), (3,4,1,0,0,0), (3,5,0,0,0,0), (3,6,0,0,0,0),
    (3,7,0,0,0,0), (4,0,0,0,0,0), (4,0,0,0,1,0), (4,0,0,1,0,0), (4,0,1,0,0,0), (4,0,1,1,0,0),
    (4,0,2,0,0,0), (4,0,3,0,0,0), (4,1,0,0,0,0), (4,1,0,0,1,0), (4,1,0,1,0,0), (4,1,1,0,0,0),
    (4,1,1,1,0,0), (4,1,2,0,0,0), (4,2,0,0,0,0), (4,2,0,1,0,0), (4,2,1,0,0,0), (4,2,2,0,0,0),
    (4,3,0,0,0,0), (4,3,0,1,0,0), (4,3,1,0,0,0), (4,4,0,0,0,0), (4,4,1,0,0,0), (4,5,0,0,0,0),
    (4,6,0,0,0,0), (5,0,0,0,0,0), (5,0,0,0,1,0), (5,0,0,1,0,0), (5,0,1,0,0,0), (5,0,1,1,0,0),
    (5,0,2,0,0,0), (5,0,3,0,0,0), (5,1,0,0,0,0), (5,1,0,1,0,0), (5,1,1,0,0,0), (5,1,2,0,0,0),
    (5,2,0,0,0,0), (5,2,0,1,0,0), (5,2,1,0,0,0), (5,2,2,0,0,0), (5,3,0,0,0,0), (5,3,1,0,0,0),
    (5,4,0,0,0,0), (5,4,1,0,0,0), (5,5,0,0,0,0), (5,6,0,0,0,0), (6,0,0,0,0,0), (6,0,0,0,1,0),
    (6,0,0,1,0,0), (6,0,1,0,0,0), (6,0,1,1,0,0), (6,0,2,0,0,0), (6,1,0,0,0,0), (6,1,0,1,0,0),
    (6,1,1,0,0,0), (6,1,2,0,0,0), (6,2,0,0,0,0), (6,2,0,1,0,0), (6,2,1,0,0,0), (6,3,0,0,0,0),
    (6,3,1,0,0,0), (6,4,0,0,0,0), (6,4,1,0,0,0), (6,5,0,0,0,0), (6,6,0,0,0,0), (7,0,0,0,0,0),
    (7,0,0,0,1,0), (7,0,0,1,0,0), (7,0,1,0,0,0), (7,0,1,1,0,0), (7,0,2,0,0,0), (7,1,0,0,0,0),
    (7,1,0,1,0,0), (7,1,1,0,0,0), (7,1,2,0,0,0), (7,2,0,0,0,0), (7,2,0,1,0,0), (7,2,1,0,0,0),
    (7,3,0,0,0,0), (7,3,1,0,0,0), (7,4,0,0,0,0), (7,5,0,0,0,0), (8,0,0,0,0,0), (8,0,0,0,1,0),
    (8,0,0,1,0,0), (8,0,1,0,0,0), (8,0,2,0,0,0), (8,1,0,0,0,0), (8,1,0,1,0,0), (8,1,1,0,0,0),
    (8,1,2,0,0,0), (8,2,0,0,0,0), (8,2,1,0,0,0), (8,3,0,0,0,0), (8,3,1,0,0,0), (8,4,0,0,0,0),
    (8,5,0,0,0,0), (9,0,0,0,0,0), (9,0,0,1,0,0), (9,0,1,0,0,0), (9,0,2,0,0,0), (9,1,0,0,0,0),
    (9,1,0,1,0,0), (9,1,1,0,0,0), (9,1,2,0,0,0), (9,2,0,0,0,0), (9,2,1,0,0,0), (9,3,0,0,0,0),
    (9,3,1,0,0,0), (9,4,0,0,0,0), (9,5,0,0,0,0), (10,0,0,0,0,0), (10,0,0,1,0,0), (10,0,1,0,0,0),
    (10,0,2,0,0,0), (10,1,0,0,0,0), (10,1,0,1,0,0), (10,1,1,0,0,0), (10,2,0,0,0,0), (10,2,1,0,0,0),
    (10,3,0,0,0,0), (10,4,0,0,0,0), (11,0,0,0,0,0), (11,0,0,1,0,0), (11,0,1,0,0,0), (11,0,2,0,0,0),
    (11,1,0,0,0,0), (11,1,0,1,0,0), (11,1,1,0,0,0), (11,2,0,0,0,0), (11,2,1,0,0,0), (11,3,0,0,0,0),
    (11,4,0,0,0,0), (12,0,0,0,0,0), (12,0,0,1,0,0), (12,0,1,0,0,0), (12,0,2,0,0,0), (12,1,0,0,0,0),
    (12,1,1,0,0,0), (12,2,0,0,0,0), (12,2,1,0,0,0), (12,3,0,0,0,0), (12,4,0,0,0,0), (13,0,0,0,0,0),
    (13,0,0,1,0,0), (13,0,1,0,0,0), (13,1,0,0,0,0), (13,1,1,0,0,0), (13,2,0,0,0,0), (13,3,0,0,0,0),
    (13,4,0,0,0,0), (14,0,0,0,0,0), (14,0,0,1,0,0), (14,0,1,0,0,0), (14,1,0,0,0,0), (14,1,1,0,0,0),
    (14,2,0,0,0,0), (14,3,0,0,0,0), (15,0,0,0,0,0), (15,0,1,0,0,0), (15,1,0,0,0,0), (15,1,1,0,0,0),
    (15,2,0,0,0,0), (15,3,0,0,0,0), (16,0,0,0,0,0), (16,0,1,0,0,0), (16,1,0,0,0,0), (16,1,1,0,0,0),
    (16,2,0,0,0,0), (16,3,0,0,0,0), (17,0,0,0,0,0), (17,0,1,0,0,0), (17,1,0,0,0,0), (17,2,0,0,0,0),
    (18,0,0,0,0,0), (18,0,1,0,0,0), (18,1,0,0,0,0), (18,2,0,0,0,0), (19,0,0,0,0,0), (19,0,1,0,0,0),
    (19,1,0,0,0,0), (19,2,0,0,0,0), (20,0,0,0,0,0), (20,1,0,0,0,0), (21,0,0,0,0,0), (21,1,0,0,0,0),
    (22,0,0,0,0,0), (22,1,0,0,0,0), (23,0,0,0,0,0), (23,1,0,0,0,0), (24,0,0,0,0,0), (25,0,0,0,0,0),
    (26,0,0,0,0,0),
]

TABLE2 = collect.(TABLE)

#####

E1_BASELINE_HYBRIDS = CSV.read("results/hybriddata/E1_BASELINE_HYBRIDS.csv", DataFrame)
E2_BASELINE_HYBRIDS = CSV.read("results/hybriddata/E2_BASELINE_HYBRIDS.csv", DataFrame)
E3_BASELINE_HYBRIDS = CSV.read("results/hybriddata/E3_BASELINE_HYBRIDS.csv", DataFrame)

E1_deltasigma_HYBRIDS = CSV.read("results/hybriddata/E1_deltasigma_HYBRIDS.csv", DataFrame)
E2_deltasigma_HYBRIDS = CSV.read("results/hybriddata/E2_deltasigma_HYBRIDS.csv", DataFrame)
E3_deltasigma_HYBRIDS = CSV.read("results/hybriddata/E3_deltasigma_HYBRIDS.csv", DataFrame)

E1_deltasigma2_HYBRIDS = CSV.read("results/hybriddata/E1_deltasigma2_HYBRIDS.csv", DataFrame)
E2_deltasigma2_HYBRIDS = CSV.read("results/hybriddata/E2_deltasigma2_HYBRIDS.csv", DataFrame)
E3_deltasigma2_HYBRIDS = CSV.read("results/hybriddata/E3_deltasigma2_HYBRIDS.csv", DataFrame)

E1_deltasigma3_HYBRIDS = CSV.read("results/hybriddata/E1_deltasigma3_HYBRIDS.csv", DataFrame)
E2_deltasigma3_HYBRIDS = CSV.read("results/hybriddata/E2_deltasigma3_HYBRIDS.csv", DataFrame)
E3_deltasigma3_HYBRIDS = CSV.read("results/hybriddata/E3_deltasigma3_HYBRIDS.csv", DataFrame)


ismin = zeros(312)
for i in 1:312
    hybrid = TABLE[i]

    if all(hybrid[2:end] .== 0)
        ismin[i] = 1
    end
    if sum(hybrid) == 1
        ismin[i] = 2
    end
end

minning = sum.(TABLE2)
data = stack(hybridnnodesnedges.(TABLE2))
hybrid_nnodes = data[1,:]
hybrid_nedges = data[2,:]

#now we write a function to automate plotting
plothybrids = function(title, hybriddata)
    y = hybriddata.meanW
    x = minning

    plt = Plots.plot(title = title, legend = false, xlabel = "Number of min() arguments", ylabel = "Expected fitness", yrotation = 45, framestyle = :box)

    Plots.scatter!(plt, x, y; mc = :black, ms = 3)
    Plots.scatter!(plt, x[ismin .== 1], y[ismin .== 1], color = :blue, ms = 7)
    Plots.scatter!(plt, x[ismin .== 2], y[ismin .== 2], color = :red, ms = 7)

    return plt
end

p1 = plothybrids("E1 HYBRIDS, sigma = 1.0", E1_BASELINE_HYBRIDS)
#savefig("E1_BASELINE_HYBRIDS.pdf", p1)
p2 = plothybrids("E2 HYBRIDS, sigma = 1.0", E2_BASELINE_HYBRIDS)
#savefig("E2_BASELINE_HYBRIDS.pdf", p2)
p3 = plothybrids("E3 HYBRIDS, sigma = 1.0", E3_BASELINE_HYBRIDS)
#savefig("E3_BASELINE_HYBRIDS.pdf", p3)

p4 = plothybrids("E1 HYBRIDS, sigma = 0.5", E1_deltasigma3_HYBRIDS)
#savefig("E1_deltasigma3_HYBRIDS.pdf", p4)
p5 = plothybrids("E2 HYBRIDS, sigma = 0.5", E2_deltasigma3_HYBRIDS)
#savefig("E1_deltasigma3_HYBRIDS.pdf", p4)
p6 = plothybrids("E3 HYBRIDS, sigma = 0.5", E3_deltasigma3_HYBRIDS)
#savefig("E1_deltasigma3_HYBRIDS.pdf", p4)

p7 = plothybrids("E1 HYBRIDS, sigma = √0.1", E1_deltasigma2_HYBRIDS)
#savefig("E1_deltasigma2_HYBRIDS.pdf", p7)
p8 = plothybrids("E2 HYBRIDS, sigma = √0.1", E2_deltasigma2_HYBRIDS)
#savefig("E2_deltasigma2_HYBRIDS.pdf", p8)
p9 = plothybrids("E3 HYBRIDS, sigma = √0.1", E3_deltasigma2_HYBRIDS)
#savefig("E3_deltasigma2_HYBRIDS.pdf", p9)

p10 = plothybrids("E1 HYBRIDS, sigma = 0.1", E1_deltasigma_HYBRIDS)
#savefig("E1_deltasigma_HYBRIDS.pdf", p10)
p11 = plothybrids("E2 HYBRIDS, sigma = 0.1", E2_deltasigma_HYBRIDS)
#savefig("E2_deltasigma_HYBRIDS.pdf", p11)
p12 = plothybrids("E3 HYBRIDS, sigma = 0.1", E3_deltasigma_HYBRIDS)
#savefig("E3_deltasigma_HYBRIDS.pdf", p12)

#alternatively we can stitch them all together and save them
my_layout = @layout [
    a  b  c  d ;   # Row 1
    e  f  g  h ;   # Row 2
    i  j  k  l     # Row 3
]

# 2. Feed the actual plot objects in the exact sequential order 
# that matches those layout slots (a through l)
combined_plot = Plots.plot(
    p1, p4, p7, p10,  # Fills Row 1 slots (a, b, c, d)
    p2, p5, p8, p11,  # Fills Row 2 slots (e, f, g, h)
    p3, p6, p9, p12,  # Fills Row 3 slots (i, j, k, l)
    layout = my_layout,
    size = (2400, 1200),
    margin = 10mm,
    left_margin = 12mm,
    title = "",
    )

savefig(combined_plot, "results/plots/selectedE_allhybridsplot.pdf")

#okay, that's the easy part, now we still need to subset just min of 5 arguments in E2,
#then label points

#labels = string.(TABLE2[minning.==5])
labels = string.(TABLE2)

#now we repeat for each sigma, but without the labeling

plot_wrtminargs! = function(slot, hybriddata, subset; title = "", with_labels = false)
    y = hybriddata.meanW
    x = hybrid_nedges
    y2 = y[minning .== subset]
    x2 = x[minning .== subset]
    labels2 = labels[minning .== subset]
    #uppers = map(t -> t[6][2], hybriddata)
    #uppers = uppers[minning .== subset]
    #lowers = map(t -> t[6][3], hybriddata)
    #lowers = lowers[minning .== subset]
    #higherrors = uppers .- y2
    #lowerrors = y2 .- lowers
    
    # Create the Axis directly in the assigned layout slot
    ax = Axis(slot,
        title = title,
        xlabelsize = 12,
        xlabel = "Network size (number of mutating parameters)",
        ylabel = "Expected fitness"
    )

    Makie.scatter!(ax, x2, y2, color=:black)
    #errorbars!(ax, x2, y2, lowerrors, higherrors)
    puremin5 = findfirst(==([subset,0,0,0,0,0]), TABLE2)
    Makie.scatter!(ax, x[puremin5], y[puremin5], color=:blue, markersize = 12)

    if with_labels
        annotation!(ax, x2, y2, 
        text=labels2, 
        fontsize=12,
        algorithm=Makie.automatic,
        )

    end
    return ax
end

# 2. Create your Master Figure first
# Note: Since you are tiling 4 plots horizontally, changed size to (1600, 400) 
# instead of (400, 1600) so they aren't compressed into columns.
master_fig = Figure(size = (1200, 300))

# 3. Direct each function call to its specific master grid cell
#ax1 = plot_wrtminargs!(master_fig[1, 1], E2_BASELINE_HYBRIDS, 5; title = "E2, sigma = 1.0, 5 min() arguments")
#ax2 = plot_wrtminargs!(master_fig[1, 2], E2_deltasigma3_HYBRIDS, 5; title = "E2, sigma = 0.5, 5 min() arguments")
#ax3 = plot_wrtminargs!(master_fig[1, 3], E2_deltasigma2_HYBRIDS, 5; title = "E2, sigma = √0.1, 5 min() arguments")
#ax4 = plot_wrtminargs!(master_fig[1, 4], E2_deltasigma_HYBRIDS, 5; title = "E2, sigma = 0.1, 5 min() arguments")

ax1 = plot_wrtminargs!(master_fig[1, 1], E2_BASELINE_HYBRIDS, 5; title = "")
ax2 = plot_wrtminargs!(master_fig[1, 2], E2_deltasigma3_HYBRIDS, 5; title = "")
ax3 = plot_wrtminargs!(master_fig[1, 3], E2_deltasigma2_HYBRIDS, 5; title = "")
ax4 = plot_wrtminargs!(master_fig[1, 4], E2_deltasigma_HYBRIDS, 5; title = "")


# 4. Clean layout padding and draw
CairoMakie.trim!(master_fig.layout)
master_fig

save("results/plots/HYBRID_E2_5minargs.pdf", master_fig)

labeledfig = Figure(size = (500,400))
ax = plot_wrtminargs!(labeledfig[1,1], E2_BASELINE_HYBRIDS, 5; title = "", with_labels = true)
labeledfig

save("results/plots/HYBRID_E2_5minargs_labeled.pdf", labeledfig)
