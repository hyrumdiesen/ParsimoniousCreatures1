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
#feedforward motif simulation traces with and without mutation
p = Nvst_hot_trajectories(PCsimtraces_ffm_m0; xlim = (0, 1300), time_bin = 10, log_scale = true, colorbar = false, title = "Simulations of Parsimonious Creatures:\nfeedforward motif, no mutation")
savefig(p, "output_figures/plots/PCsimtraces_ffm_m0.pdf")
p = Nvst_hot_trajectories(PCsimtraces_ffm; xlim = (0, 1300), time_bin = 10, log_scale = true, colorbar = false, title = "Simulations of Parsimonious Creatures:\nfeedforward motif, yes mutation")
savefig(p, "output_figures/plots/PCsimtraces_ffm.pdf")
#simple step motif simulation traces with and without mutation
p = Nvst_hot_trajectories(PCsimtraces_ssm_m0; xlim = (0, 1300), time_bin = 10, log_scale = true, colorbar = false, title = "Simulations of Parsimonious Creatures:\nsimple step motif, no mutation")
savefig(p, "output_figures/plots/PCsimtraces_ssm_m0.pdf")
p = Nvst_hot_trajectories(PCsimtraces_ssm; xlim = (0, 1300), time_bin = 10, log_scale = true, colorbar = false, title = "Simulations of Parsimonious Creatures:\nsimple step motif, yes mutation")
savefig(p, "output_figures/plots/PCsimtraces_ssm.pdf")

#figure 2 simulations of parsimonious creatures
NUM_GENOTYPES = 25

#ssm min redundancy
p = Nvst_hot_trajectories(PCsimtraces_min_ssm; xlim = (0, 1300), time_bin = 10, log_scale = true, colorbar = false, title = "Simulations of Parsimonious Creatures:\nSSM min redundancy")
savefig(p, "output_figures/plots/PCsimtraces_min_ssm.pdf")
#ssm averaging redundancy
p = Nvst_hot_trajectories(PCsimtraces_avg_ssm_real; xlim = (0, 1300), time_bin = 10, log_scale = true, colorbar = false, title = "Simulations of Parsimonious Creatures:\nSSM averaging redundancy")
savefig(p, "output_figures/plots/PCsimtraces_avg_ssm_real.pdf")



#bureaucratic paranoia figure simulations of parsimonious creatures
p1 = Nvst_hot_trajectories(PCsimtraces_min_ssm_1; linewidthbottomplot = 3, xlim = (0, 1000), time_bin = 10, log_scale = true, colorbar = false, title = "")
p2 = Nvst_hot_trajectories(PCsimtraces_min_ssm_5; linewidthbottomplot = 3, xlim = (0, 1000), time_bin = 10, log_scale = true, colorbar = false, title = "")
p3 = Nvst_hot_trajectories(PCsimtraces_min_ssm_10; linewidthbottomplot = 3, xlim = (0, 1000), time_bin = 10, log_scale = true, colorbar = false, title = "")
p4 = Nvst_hot_trajectories(PCsimtraces_min_ssm_20; linewidthbottomplot = 3, xlim = (0, 1000), time_bin = 10, log_scale = true, colorbar = false, title = "")

#now combine them together
p5 = Plots.plot(p1, p2, p3, p4, layout = (1,4), size = (1700,500), margin = 1mm, bottom_margin = 10mm, extra_plot_kwargs = Any[:tight_layout => true])

p5 = Plots.plot!(p5, left_margin = 12mm, subplot = 1)
savefig(p5, "output_figures/plots/PCsimtraces_bureaucraticparanoia.pdf")


#supplemental figure averaging redundancy approximation

p = Nvst_hot_trajectories(PCsimtraces_avg_ssm_real_3; xlim = (0, 600), time_bin = 10, log_scale = true, colorbar = false, title = "Simulations of Parsimonious Creatures:\nSSM averaging redundancy, real")
savefig(p, "output_figures/plots/PCsimtraces_avg_ssm_real_3.pdf")

p = Nvst_hot_trajectories(PCsimtraces_avg_ssm_fake_3; xlim = (0, 600), time_bin = 10, log_scale = true, colorbar = false, title = "Simulations of Parsimonious Creatures:\nSSM averaging redundancy, approximated")
savefig(p, "output_figures/plots/PCsimtraces_avg_ssm_fake_3.pdf")
