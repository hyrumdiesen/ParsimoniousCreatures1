#Nvst testing a single feedforward motif
include(joinpath(pwd(), "scripts/bones_ode.jl"))

KCmatrix,rs,phis = feedforwardavgvar1(1)
NUM_GENOTYPES = 50
mutmag = 0.5
m = 0.0
PCsimtraces_ffm_m0 = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 10000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_ffm_m0.h5", PCsimtraces_ffm_m0;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_ffm_m0")

m = 0.0005
PCsimtraces_ffm = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 10000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_ffm.h5", PCsimtraces_ffm;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_ffm")

#now we compare single ssm motif

KCmatrix,rs,phis = separatecircuits7(1)

m = 0.0
PCsimtraces_ssm_m0 = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 10000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_ssm_m0.h5", PCsimtraces_ssm_m0;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_ssm_m0")

m = 0.0005
PCsimtraces_ssm = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 10000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_ssm.h5", PCsimtraces_ssm;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_ssm")

#testing a min network
NUM_GENOTYPES = 25
m = 0.0005
KCmatrix,rs,phis = separatecircuits7(4)
PCsimtraces_min_ssm = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 3000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_min_ssm.h5", PCsimtraces_min_ssm;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_min_ssm")

#in order to make ssm averaging redundancy competitive in this environment at the same timescales for visualization purposes,
#we are going to test an absurdly huge averaging network of over 500 parameters. This will take a few minutes.
KCmatrix,rs,phis = checkerboard7(16)
PCsimtraces_avg_ssm_real = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 3000, time_bin = 1)
p = Nvst_hot_trajectories(PCsimtraces_avg_ssm_real; xlim = (0, 1300), time_bin = 10, log_scale = true, colorbar = false, title = "Simulations of Parsimonious Creatures:\nfeedforward motif, no mutation")

Nvst_save("datasets/Nvstdata/PCsimtraces_avg_ssm_real.h5", PCsimtraces_avg_ssm_real;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_avg_ssm_real")


#now we test bureaucratic paranoia
NUM_GENOTYPES = 25
#we need to jack up the mutation rate to make bureaucratic paranoia visible with relatively small
#min redundancy networks.
m = 0.003
mutmag = 1.0
KCmatrix,rs,phis = separatecircuits7(1)
PCsimtraces_min_ssm_1 = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 3000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_min_ssm_1.h5", PCsimtraces_min_ssm_1;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_min_ssm_1")

KCmatrix,rs,phis = separatecircuits7(5)
PCsimtraces_min_ssm_5 = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 3000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_min_ssm_5.h5", PCsimtraces_min_ssm_5;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_min_ssm_5")

KCmatrix,rs,phis = separatecircuits7(10)
PCsimtraces_min_ssm_10 = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 3000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_min_ssm_10.h5", PCsimtraces_min_ssm_10;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_min_ssm_10")

KCmatrix,rs,phis = separatecircuits7(20)
PCsimtraces_min_ssm_20 = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 3000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_min_ssm_20.h5", PCsimtraces_min_ssm_20;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_min_ssm_20")

#now we test a 3 fold averaging redundancy network against a single simple step motif where the mutation magnitude and mutation
#rate are normalized by the factor of extra parameters in the averaging redundancy network
m = 0.0005
mutmag = 0.5
KCmatrix,rs,phis = checkerboard7(3)
PCsimtraces_avg_ssm_real_3 = Nvst_data(KCmatrix, rs, phis; num_sims = 200, MAX_TIME = 3000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_avg_ssm_real_3.h5", PCsimtraces_avg_ssm_real_3;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_avg_ssm_real_3")

#now for the approximation of an averaging network with a single simple step motif
KCmatrix2,rs2,phis2 = checkerboard7(1)

avg_approximation_factor = (sum(KCmatrix .> 0) + sum(rs .> 0)) / (sum(KCmatrix2 .> 0) + sum(rs2 .> 0))
m = m*avg_approximation_factor
mutmag = mutmag/avg_approximation_factor
PCsimtraces_avg_ssm_fake_3 = Nvst_data(KCmatrix2, rs2, phis2; num_sims = 200, MAX_TIME = 3000, time_bin = 1)
Nvst_save("datasets/Nvstdata/PCsimtraces_avg_ssm_fake_3.h5", PCsimtraces_avg_ssm_fake_3;
    params = (; KCmatrix, rs, phis, m, mutmag,
                num_sims = 200, MAX_TIME = 5000, NUM_GENOTYPES = NUM_GENOTYPES),
    figure = "PCsimtraces_avg_ssm_fake_3")
