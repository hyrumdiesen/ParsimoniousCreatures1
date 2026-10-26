using Plots
using Plots.PlotMeasures
using ImageFiltering
using JLD2

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
    p
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
        yflip = false,
        xlabel = "Organism Death Rate",
        ylabel = "Mutation Rate (/cell/edge)",
        #title = "N* = 12, sigma = 1.0",
        title = title,
        bottom_margin = 8mm,
        tick_direction = :out)
end

#plot simplex legend
plot_simplex_legend(; labels=("FF var2","SSM Avg","SSM Min"))


#format is [1,2,3] from small to large network cutoff
@load "/data/plotting_data/heatmaps/processed_r12_s1p0_paper1.jld2" r12_s1p0_paper1

p1 = plot_dominance("Small Networks (nparams < 9) N* = 12, σ = 1.0", r12_s1p0_paper1[1][1][:,:,3],r12_s1p0_paper1[1][1][:,:,1],r12_s1p0_paper1[1][1][:,:,2])
p2 = plot_dominance("Medium Networks (nparams < 30) N* = 12, σ = 1.0", r12_s1p0_paper1[1][2][:,:,3],r12_s1p0_paper1[1][2][:,:,1],r12_s1p0_paper1[1][2][:,:,2])
p3 = plot_dominance("Large Networks (nparams < 80) N* = 12, σ = 1.0", r12_s1p0_paper1[1][3][:,:,3],r12_s1p0_paper1[1][3][:,:,1],r12_s1p0_paper1[1][3][:,:,2])

savefig(p1, "r12_s1p0_smallnetworks_heatmap.pdf")
savefig(p2, "r12_s1p0_mediumnetworks_heatmap.pdf")
savefig(p3, "r12_s1p0_largenetworks_heatmap.pdf")


@load "/data/plotting_data/heatmaps/processed_r25_s1p0_paper1.jld2" r25_s1p0_paper1

p4 = plot_dominance("Small Networks (nparams < 9) N* = 25, σ = 1.0", r25_s1p0_paper1[1][1][:,:,3],r25_s1p0_paper1[1][1][:,:,1],r25_s1p0_paper1[1][1][:,:,2])
p5 = plot_dominance("Medium Networks (nparams < 30) N* = 25, σ = 1.0", r25_s1p0_paper1[1][2][:,:,3],r25_s1p0_paper1[1][2][:,:,1],r25_s1p0_paper1[1][2][:,:,2])
p6 = plot_dominance("Large Networks (nparams < 80) N* = 25, σ = 1.0", r25_s1p0_paper1[1][3][:,:,3],r25_s1p0_paper1[1][3][:,:,1],r25_s1p0_paper1[1][3][:,:,2])

savefig(p4, "r25_s1p0_smallnetworks_heatmap.pdf")
savefig(p5, "r25_s1p0_mediumnetworks_heatmap.pdf")
savefig(p6, "r25_s1p0_largenetworks_heatmap.pdf")


@load "/data/plotting_data/heatmaps/processed_r12_s0p1_paper1.jld2" r12_s0p1_paper1

plot_dominance("Small Networks (nparams < 9) N* = 12, σ = 0.1", r12_s0p1_paper1[1][1][:,:,3],r12_s0p1_paper1[1][1][:,:,1],r12_s0p1_paper1[1][1][:,:,2])
plot_dominance("Medium Networks (nparams < 30) N* = 12, σ = 0.1", r12_s0p1_paper1[1][2][:,:,3],r12_s0p1_paper1[1][2][:,:,1],r12_s0p1_paper1[1][2][:,:,2])
plot_dominance("Large Networks (nparams < 80) N* = 12, σ = 0.1", r12_s0p1_paper1[1][3][:,:,3],r12_s0p1_paper1[1][3][:,:,1],r12_s0p1_paper1[1][3][:,:,2])

savefig(p7, "r12_s0p1_smallnetworks_heatmap.pdf")
savefig(p8, "r12_s0p1_mediumnetworks_heatmap.pdf")
savefig(p9, "r12_s0p1_largenetworks_heatmap.pdf")