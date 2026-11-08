include(joinpath(pwd(), "scripts/bones_ode.jl"))

#we could just generate the data, but it would be a bit too slow for reviewers, instead we will load pregenerated data from .h5 files
PCsimtraces_ffm_m0 = Nvst_load("datasets/Nvstdata/PCsimtraces_ffm_m0.h5")
PCsimtraces_ffm = Nvst_load("datasets/Nvstdata/PCsimtraces_ffm.h5")
PCsimtraces_ssm_m0 = Nvst_load("datasets/Nvstdata/PCsimtraces_ssm_m0.h5")
PCsimtraces_ssm = Nvst_load("datasets/Nvstdata/PCsimtraces_ssm.h5")
PCsimtraces_min_ssm = Nvst_load("datasets/Nvstdata/PCsimtraces_min_ssm.h5")
PCsimtraces_avg_ssm_real = Nvst_load("datasets/Nvstdata/PCsimtraces_avg_ssm_real.h5")
PCsimtraces_avg_ssm_real_3 = Nvst_load("datasets/Nvstdata/PCsimtraces_avg_ssm_real_3.h5")
PCsimtraces_avg_ssm_fake_3 = Nvst_load("datasets/Nvstdata/PCsimtraces_avg_ssm_fake_3.h5")

PCsimtraces_min_ssm_1 = Nvst_load("datasets/Nvstdata/PCsimtraces_min_ssm_1.h5")
PCsimtraces_min_ssm_5 = Nvst_load("datasets/Nvstdata/PCsimtraces_min_ssm_5.h5")
PCsimtraces_min_ssm_10 = Nvst_load("datasets/Nvstdata/PCsimtraces_min_ssm_10.h5")
PCsimtraces_min_ssm_20 = Nvst_load("datasets/Nvstdata/PCsimtraces_min_ssm_20.h5")


#figure 1 simulations of parsimonious creatures
NUM_GENOTYPES = 50
xlim = 1300
positions = LinRange(100, 1200, 5)

#feedforward motif simulation traces with and without mutation
p = Nvst_figure(PCsimtraces_ffm_m0; positions = positions, xlim = (0,xlim))
savefig(p, "output_figures/plots/PCsimtraces_ffm_m0.pdf")

p = Nvst_figure(PCsimtraces_ffm; positions = positions, xlim = (0,xlim))
savefig(p, "output_figures/plots/PCsimtraces_ffm.pdf")
#simple step motif simulation traces with and without mutation

p = Nvst_figure(PCsimtraces_ssm_m0; positions = positions, xlim = (0,xlim))
savefig(p, "output_figures/plots/PCsimtraces_ssm_m0.pdf")

p = Nvst_figure(PCsimtraces_ssm; positions = positions, xlim = (0,xlim))
savefig(p, "output_figures/plots/PCsimtraces_ssm.pdf")

#figure 2 simulations of parsimonious creatures
NUM_GENOTYPES = 25

#ssm min redundancy
p = Nvst_figure(PCsimtraces_min_ssm; positions = positions, xlim = (0,xlim))
savefig(p, "output_figures/plots/PCsimtraces_min_ssm.pdf")

#ssm averaging redundancy
p = Nvst_figure(PCsimtraces_avg_ssm_real; positions = positions, xlim = (0,xlim))
savefig(p, "output_figures/plots/PCsimtraces_avg_ssm_real.pdf")



#bureaucratic paranoia figure simulations of parsimonious creatures
NUM_GENOTYPES = 25
xlim = 1000
positions = LinRange(100,900,5)

p1 = Nvst_figure(PCsimtraces_min_ssm_1; positions = positions, xlim = (0,xlim))
p2 = Nvst_figure(PCsimtraces_min_ssm_5; positions = positions, xlim = (0,xlim))
p3 = Nvst_figure(PCsimtraces_min_ssm_10; positions = positions, xlim = (0,xlim))
p4 = Nvst_figure(PCsimtraces_min_ssm_20; positions = positions, xlim = (0,xlim))

savefig(p1, "output_figures/plots/PCsimtraces_bureaucraticparanoia1.pdf")
savefig(p2, "output_figures/plots/PCsimtraces_bureaucraticparanoia2.pdf")
savefig(p3, "output_figures/plots/PCsimtraces_bureaucraticparanoia3.pdf")
savefig(p4, "output_figures/plots/PCsimtraces_bureaucraticparanoia4.pdf")


#now combine them together
#p5 = Plots.plot(p1, p2, p3, p4, layout = (1,4), size = (1700,500), margin = 1mm, bottom_margin = 3mm, extra_plot_kwargs = Any[:tight_layout => true])

#p5 = Plots.plot!(p5, subplot = 1)
#savefig(p5, "output_figures/plots/PCsimtraces_bureaucraticparanoia.pdf")


#supplemental figure averaging redundancy approximation

NUM_GENOTYPES = 25
xlim = 600
positions = LinRange(50,550,5)

p = Nvst_figure(PCsimtraces_avg_ssm_real_3; positions = positions, xlim = (0,xlim))
savefig(p, "output_figures/plots/PCsimtraces_avg_ssm_real_3.pdf")

p = Nvst_figure(PCsimtraces_avg_ssm_fake_3; positions = positions, xlim = (0,xlim))
savefig(p, "output_figures/plots/PCsimtraces_avg_ssm_fake_3.pdf")


