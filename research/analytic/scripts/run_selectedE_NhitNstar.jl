#run bones file to initialize parameters and functions
include("bones_selectedE.jl")

println("Threads available: ", Threads.nthreads())

mutmag = 1.0
root = 12
NUM_GENOTYPES = 50
FLOOR = 0.1
germlinem = 0.0
gammam = 0.0
Ninit = 1

#initialize analytic_sim functions

test = constructffmodmatrix(1)
ffavg_analytic_sim_nodemandgermm(test)
test = constructffvar2modmatrix(1)
ffvar2_analytic_sim_nodemandgermm(test)
anyssmnet_analytic_sim_nodemandgermm([1,0,0,0,0,0])

NhitNstartoggle = 1

testdf = CSV.read("data/Nstar12and25_baseline.csv", DataFrame)

testdf2 = testallnetworks_nooptim("data/NhitNstar12and25_poo.csv", testdf[1:2,:])