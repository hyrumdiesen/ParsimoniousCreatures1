#code for generating network topologies and replication functions, as well as samplings of mutant replication functions
include(joinpath(pwd(), "scripts/bones_ode.jl"))

#single feedforward motif
KCmatrix,rs,phis = feedforwardavgvar1(1)
p =plot_weighted_network_makie(transpose(KCmatrix), rs, phis; figwidth = 400, rowwidth = 0.4, edgecurve = 0.1, curves = false)
save("output_figures/topology_diagrams/topo_ffmotif.pdf", p)
p =plot_weighted_network_makie(transpose(KCmatrix), rs, phis; figwidth = 400, rowwidth = 0.4, edgecurve = 0.1, curves = false, interior_node_labeling = true)
save("output_figures/topology_diagrams/topo_ffmotif_labeled.pdf", p)
#sample mutant replication functions
p = plotmanymutants_withTFs(KCmatrix,rs,phis, 1.0, 20)
savefig(p, "output_figures/plots/repfunc_ffmotif.pdf")

#averaging redundancy feedforward motif var 1
KCmatrix2,rs2,phis2 = duplicate_nodes(KCmatrix*1.0, rs*1.0,phis*1.0, [1,2,1,2])
p = plot_weighted_network_makie(transpose(KCmatrix2), rs2, phis2; figwidth = 400, rowwidth = 0.4, edgecurve = 0.1, curves = false)
save("output_figures/topology_diagrams/topo_ffavgredvar1.pdf", p)

#averaging redundancy feedforward motif var 2
KCmatrix,rs,phis = feedforwardavgvar1(3)
#KCmatrix2,rs2,phis2 = duplicate_nodes(KCmatrix*1.0, rs*1.0,phis*1.0, [1,2,1,2])
p = plot_weighted_network_makie(transpose(KCmatrix), rs, phis; figwidth = 400, rowwidth = 0.4, edgecurve = 0.1, curves = false)
save("output_figures/topology_diagrams/topo_ffavgredvar2.pdf", p)

#single simple step motif
KCmatrix,rs,phis = separatecircuits7(1)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 400, rowwidth = 0.4, edgecurve = 0.15)
save("output_figures/topology_diagrams/topo_ssmotif.pdf", p)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 400, rowwidth = 0.4, edgecurve = 0.15, interior_node_labeling = true)

save("output_figures/topology_diagrams/topo_ssmotif_labeled.pdf", p)

#sample mutant replication functions
p = plotmanymutants_withTFs(KCmatrix,rs,phis, 1.0, 40)
savefig(p, "output_figures/plots/repfunc_ssmotif.pdf")


#averaging redundancy simple step motif
KCmatrix,rs,phis = checkerboard7(3)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 400, rowwidth = 0.4, edgecurve = 0.1)
save("output_figures/topology_diagrams/topo_ssmavgred.pdf", p)

#sample mutant replication functions
p = plotmanymutants(KCmatrix,rs,phis, 2.0, 30)
savefig(p, "output_figures/plots/repfunc_ssmavgred.pdf")


#min redundancy simple step motif
KCmatrix,rs,phis = separatecircuits7(3)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 400, rowwidth = 0.4, edgecurve = 0.1)
save("output_figures/topology_diagrams/topo_ssmminred_3.pdf", p)

p = plotmanymutants(KCmatrix,rs,phis, 4.0, 30)

savefig(p, "output_figures/plots/repfunc_ssmminred.pdf")


#compare averaging redundancy real and approximated
#sample mutant replication functions
KCmatrix,rs,phis = checkerboard7(3)
p = plotmanymutants(KCmatrix,rs,phis, 2.0, 30, NUM_GENOTYPES = 25)
savefig(p, "output_figures/plots/repfunc_ssmavgred_real.pdf")

KCmatrix,rs,phis = separatecircuits7(1)
p = plotmanymutants(KCmatrix,rs,phis, 2.0, 30; mutmag = mutmag/7, NUM_GENOTYPES = 25)
savefig(p, "output_figures/plots/repfunc_ssmavgred_fake.pdf")

####plot bureaucratic paranoia min redundancy topologies
KCmatrix,rs,phis = separatecircuits7(5)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 600, rowwidth = 0.4, edgecurve = 0.08)
save("output_figures/topology_diagrams/topo_ssmminred_5.pdf", p)

KCmatrix,rs,phis = separatecircuits7(10)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 1000, rowwidth = 0.4, edgecurve = 0.03)
save("output_figures/topology_diagrams/topo_ssmminred_10.pdf", p)

KCmatrix,rs,phis = separatecircuits7(20)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 1100, rowwidth = 0.4, edgecurve = 0.03)
save("output_figures/topology_diagrams/topo_ssmminred_20.pdf", p)



#plot cycles of duplication and subfunctionalization
hybridvect = [0,1,0,0,0,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 300, rowwidth = 0.4, edgecurve = 0.15)
save("output_figures/topology_diagrams/topo_hybrid[0,1,0,0,0,0].pdf", p)

hybridvect = [0,0,1,0,0,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 500, rowwidth = 0.3, edgecurve = 0.1)
save("output_figures/topology_diagrams/topo_hybrid[0,0,1,0,0,0].pdf", p)

hybridvect = [0,0,0,1,0,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 600, rowwidth = 0.3, edgecurve = 0.07)
save("output_figures/topology_diagrams/topo_hybrid[0,0,0,1,0,0].pdf", p)

hybridvect = [2,0,0,0,0,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 300, rowwidth = 0.4, edgecurve = 0.18)
save("output_figures/topology_diagrams/topo_hybrid[2,0,0,0,0,0].pdf", p)

hybridvect = [1,1,0,0,0,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 500, rowwidth = 0.4, edgecurve = 0.1)
save("output_figures/topology_diagrams/topo_hybrid[1,1,0,0,0,0].pdf", p)

hybridvect = [3,0,0,0,0,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 500, rowwidth = 0.4, edgecurve = 0.1)
save("output_figures/topology_diagrams/topo_hybrid[3,0,0,0,0,0].pdf", p)


######plot hybrid network topologies
hybridvect = [5,0,0,0,0,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 600, rowwidth = 0.4, edgecurve = 0.06)
save("output_figures/topology_diagrams/topo_hybrid[5,0,0,0,0,0].pdf", p)


hybridvect = [0,5,0,0,0,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 1200, rowwidth = 0.4, edgecurve = 0.03)
save("output_figures/topology_diagrams/topo_hybrid[0,5,0,0,0,0].pdf", p)

hybridvect = [4,1,0,0,0,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 800, rowwidth = 0.4, edgecurve = 0.05)
save("output_figures/topology_diagrams/topo_hybrid[4,1,0,0,0,0].pdf", p)

hybridvect = [4,0,0,0,1,0]
KCmatrix,rs,phis = makespecifichybrid(hybridvect)
p = plot_weighted_network_makie(KCmatrix, rs, phis; figwidth = 1200, rowwidth = 0.4, edgecurve = 0.03)
save("output_figures/topology_diagrams/topo_hybrid[4,0,0,0,1,0].pdf", p)