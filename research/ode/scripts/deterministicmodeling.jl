include(joinpath(pwd(), "scripts/bones_ode.jl"))

#running a deterministic model of parsimonious creatures without mutation. Then optimizing to naiively discover the 
#feedforward motif and the simple step motif.

KCmatrix,rs,phis = checkerboard7(1)
test = solvegrowthfunction_withTFs(50, KCmatrix,rs,phis)
test2 = matrixfunc_celldeath(test[1])

Plots.bar(test2)

#that's on average 84% of maximum possible fitness, where 100% maximum fitness is
#an organism that is spontaneously born with 12 cells and never deviates from 12 cells until it dies.

#now we can make an objective function
objective_naiive2nodes = function(x; returnqvect = false)
    x = exp.(x)
    KCmatrix = [x[1] x[2]; x[3] x[4]]
    rs = [x[5], x[6]]
    phis = [x[7], x[8]]

    lambdas = solvegrowthfunction_withTFs(50, KCmatrix,rs,phis)[1]
    qvect = matrixfunc_celldeath(lambdas)

    fitness = determ_fitness(qvect)*deltaO

    if returnqvect
        return qvect, KCmatrix, rs, phis
    else
        return -fitness
    end
end

objective_naiive2nodes_nofeedback = function(x; returnqvect = false)
    x = exp.(x)
    KCmatrix = [0.0 0.0; x[1] 0.0]
    rs = [x[2], x[3]]
    phis = [x[4], x[5]]

    lambdas = solvegrowthfunction_withTFs(50, KCmatrix,rs,phis)[1]
    qvect = matrixfunc_celldeath(lambdas)

    fitness = determ_fitness(qvect)*deltaO

    if returnqvect
        return qvect, KCmatrix, rs, phis
    else
        return -fitness
    end
end

#setting a strong upper bound interaction weight of 1e8
ffm_found = bboptimize(objective_naiive2nodes_nofeedback, SearchRange = (-5.0, log(1e8)), NumDimensions = 5, MaxSteps = 10000)

ffm_found_data = objective_naiive2nodes_nofeedback(best_candidate(ffm_found); returnqvect=true)
Plots.bar(ffm_found_data[1], xlabel = "Number of cells (N)", ylabel = "Probability", title = "Deterministic Model:\n Discovered Feedforward Motif")
#where a bar plot shows the average probability of an organism being at a particular number of cells
#at any point in it's lifetime given an infinite population of organisms.

#compare this to the feedforward motif used in the data

KCmatrix,rs,phis = feedforwardavgvar1(1)
ffm_used_data = matrixfunc_celldeath(solvegrowthfunction_withTFs(50, KCmatrix,rs,phis)[1] )
Plots.bar(ffm_used_data, xlabel = "Number of cells (N)", ylabel = "Probability", title = "Deterministic Model:\n Feedforward Motif Used in Analysis")

#discovering the simple step motif is a bit finicky because it may stuck on local optima 
#and essentially form a feedforward motif with a fitness of ~0.37
#if this happens just try again,
#if the optimization is successful at finding the sime step motif, it will go to a fitness of ~0.86
ssm_found = bboptimize(objective_naiive2nodes, SearchRange = (-5.0, log(1e8)), NumDimensions = 8, MaxSteps = 10000)

ssm_found_data = objective_naiive2nodes(best_candidate(ssm_found); returnqvect=true)
Plots.bar(ssm_found_data[1], xlabel = "Number of cells (N)", ylabel = "Probability", title = "Deterministic Model:\n Discovered Simple Step Motif")

KCmatrix,rs,phis = separatecircuits7(1)
ssm_used_data = matrixfunc_celldeath(solvegrowthfunction_withTFs(50, KCmatrix,rs,phis)[1] )
Plots.bar(ssm_used_data, xlabel = "Number of cells (N)", ylabel = "Probability", title = "Deterministic Model:\n Simple Step Motif Used in Analysis")

#just to show how close the optimized values are numerically
println("DISCOVERED SIMPLE STEP MOTIF:")
println("rs")
println(ssm_found_data[3])
println("KCmatrix:")
println(ssm_found_data[2])
println("phis")
println(ssm_found_data[4])
println("SIMPLE STEP MOTIF USED IN DATA")
println("rs")
println(separatecircuits7(1)[2])
println("KCmatrix:")
println(separatecircuits7(1)[1])
println("phis")
println(separatecircuits7(1)[3])