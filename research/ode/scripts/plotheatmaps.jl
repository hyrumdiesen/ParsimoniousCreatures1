include("bones_ode.jl")
#plot simplex legend
p = plot_simplex_legend(; labels=("FF var2","SSM Avg","SSM Min"))
savefig(p, "results/plots/heatmap_legend.pdf")

#format is [1,2,3] from small to large network cutoff for the first dimention of the tensor

data = h5open("results/heatmap_data/heatmaps_processed/processed_r12_s1p0_paper1.h5", "r") do file
    read(file, "my_tensor")
end

p1 = plot_dominance("Small Networks (nparams < 9)\n N* = 12, σ = 1.0", data[1,1,:,:,3],data[1,1,:,:,1],data[1,1,:,:,2])
p2 = plot_dominance("Medium Networks (nparams < 30)\n N* = 12, σ = 1.0", data[2,1,:,:,3],data[2,1,:,:,1],data[2,1,:,:,2])
#p3 = plot_dominance("Large Networks (nparams < 80)\n N* = 12, σ = 1.0", data[3,1,:,:,3],data[3,1,:,:,1],data[3,1,:,:,2])
p3 = plot_dominance("", data[3,1,:,:,3],data[3,1,:,:,1],data[3,1,:,:,2])

#savefig(p1, "r12_s1p0_smallnetworks_heatmap.pdf")
#savefig(p2, "r12_s1p0_mediumnetworks_heatmap.pdf")
savefig(p3, "temp.pdf")
run(`gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.3 -dNOPAUSE -dBATCH -sOutputFile=results/plots/r12_s1p0_largenetworks_heatmap.pdf temp.pdf`)


data = h5open("results/heatmap_data/heatmaps_processed/processed_r25_s1p0_paper1.h5", "r") do file
    read(file, "my_tensor")
end

p4 = plot_dominance("Small Networks (nparams < 9)\n N* = 25, σ = 1.0", data[1,1,:,:,3],data[1,1,:,:,1],data[1,1,:,:,2])
p5 = plot_dominance("Medium Networks (nparams < 30)\n N* = 25, σ = 1.0", data[2,1,:,:,3],data[2,1,:,:,1],data[2,1,:,:,2])
#p6 = plot_dominance("Large Networks (nparams < 80)\n N* = 25, σ = 1.0", data[3,1,:,:,3],data[3,1,:,:,1],data[3,1,:,:,2])
p6 = plot_dominance("", data[3,1,:,:,3],data[3,1,:,:,1],data[3,1,:,:,2])


#savefig(p4, "results/plots/r25_s1p0_smallnetworks_heatmap.pdf")
#savefig(p5, "results/plots/r25_s1p0_mediumnetworks_heatmap.pdf")
#savefig(p6, "results/plots/r25_s1p0_largenetworks_heatmap.pdf")
savefig(p6, "temp.pdf")
run(`gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.3 -dNOPAUSE -dBATCH -sOutputFile=results/plots/r25_s1p0_largenetworks_heatmap.pdf temp.pdf`)


data = h5open("results/heatmap_data/heatmaps_processed/processed_r12_s0p1_paper1.h5", "r") do file
    read(file, "my_tensor")
end

p7 = plot_dominance("Small Networks (nparams < 9)\n N* = 12, σ = 0.1", data[1,1,:,:,3],data[1,1,:,:,1],data[1,1,:,:,2])
p8 = plot_dominance("Medium Networks (nparams < 30)\n N* = 12, σ = 0.1", data[2,1,:,:,3],data[2,1,:,:,1],data[2,1,:,:,2])
#p9 = plot_dominance("Large Networks (nparams < 80)\n N* = 12, σ = 0.1", data[3,1,:,:,3],data[3,1,:,:,1],data[3,1,:,:,2])
p9 = plot_dominance("", data[3,1,:,:,3],data[3,1,:,:,1],data[3,1,:,:,2])

#savefig(p7, "results/plots/r12_s0p1_smallnetworks_heatmap.pdf")
#savefig(p8, "results/plots/r12_s0p1_mediumnetworks_heatmap.pdf")
#savefig(p9, "results/plots/r12_s0p1_largenetworks_heatmap.pdf")
savefig(p9, "temp.pdf")
run(`gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.3 -dNOPAUSE -dBATCH -sOutputFile=results/plots/r12_s0p1_largenetworks_heatmap.pdf temp.pdf`)


#now combine all heatmaps into one figure

#alternatively we can stitch them all together and save them
my_layout = @layout [
    a  b  c ;   # Row 1
    d  e  f ;   # Row 2
    g  h  i     # Row 3
]

# 2. Feed the actual plot objects in the exact sequential order 
# that matches those layout slots (a through l)
combined_plot = Plots.plot(
    p1, p4, p7,  # Fills Row 1 slots (a, b, c, d)
    p2, p5, p8,  # Fills Row 2 slots (e, f, g, h)
    p3, p6, p9,  # Fills Row 3 slots (i, j, k, l)
    layout = my_layout, 
    size = (1200, 1200),
    bottommargin = -8mm,
    topmargin = -8mm,
    leftmargin = 5mm,
    title = ""
)

#savefig(combined_plot, "results/plots/allheatmaps.pdf")
savefig(combined_plot, "temp.pdf")
run(`gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.3 -dNOPAUSE -dBATCH -sOutputFile=results/plots/allheatmaps.pdf temp.pdf`)
