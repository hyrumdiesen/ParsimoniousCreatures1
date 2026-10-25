#plotting base selected environments
include("bones_ode.jl")
df = DataFrame(CSV.File("results/selectedEdata/selectedEdata_combined.csv"))

plotselectedE = function(df, E, delta, ylimW)

    whichcolor = ["ff_avg_var1", "ff_avg_var2", "ssm_avg", "ssm_min"]
    #colors = [7, :lime, :red, :green, :blue]
    colors = [:pink, :red, :green, :blue]

    groupedEdelta = groupby(df, [:delta, :E])
    #for each plot

    networks = groupby(groupedEdelta[(delta = delta, E = E)], :network_type)

    color = colors[findfirst(whichcolor .== networks[1].network_type[1])]
    #p = Plots.plot(networks[1].nedges, networks[1].meanW, yerr = (networks[1].meanWupper, networks[1].meanWlower))
    #p = Plots.plot(legend = :none, size = (300,300), ylim = (0, maxWs[i]))
    #p = Plots.plot(size = (300,300), ylim = (0, ylimWs[i]))
    p = Plots.plot(size = (300,300))
    #p = Plots.plot(size = (300,300), xlabel = "Network size (mutating parameters)", ylabel = "Mean fitness")
    j = 0

    for net in networks #this iterates over each network type
        j+= 1

        #we only want to bother doing the maximum of two runs thing if there is a
        #rep 2 for that network type.
        
        networktype = networks[j].network_type[1]

        means = zeros(size(networks[j])[1])
        uppers = zeros(size(networks[j])[1])
        lowers = zeros(size(networks[j])[1])

        for k in 1:length(networks[j].network_type) #this iterates over each size of a particular network

            networktype = networks[j].network_type[k]

            means[k] = networks[j].meanW[k]
            #intentionally switching these since they were mislabeled
            uppers[k] = networks[j].meanWlower[k] .- networks[j].meanW[k]  
            lowers[k] = networks[j].meanW[k] .- networks[j].meanWupper[k]

        end
        color = colors[findfirst(whichcolor .== networks[j].network_type[1])]
        xs = networks[j].nedges
        xlim = (0,80)
        #if delta == "nodem_adjusted"
        if contains(delta, "nodem")
            xs = networks[j].nnodes
            xlim = (0,20.3)
        end

        if networktype == "ff_avg_var1" || networktype == "ff_avg_var2"
            style = :dot
        else
            style = :solid
        end
        if networktype == "ff_avg"
            leg = "feedforward averaging"
        elseif networktype == "ff_avg_var1"
            leg = "feedforward averaging (var1)"
        elseif networktype == "ff_avg_var2"
            leg = "feedforward averaging (var2)"
        elseif networktype == "ssm_avg"
            leg = "ssm averaging"
        elseif networktype == "ssm_min"
            leg = "ssm min"
        else
            leg = :none
        end
        
        #Plots.plot!(p, xs, means, yerr = (uppers, lowers), size = (300,300), linewidth = 3, label = leg, color = color, markerstrokecolor=color, xlim = xlim, ylim = (0,ylimW), linestyle = style)
        Plots.plot!(p, xs, means, yerr = (lowers, uppers), size = (300,300), linewidth = 3, label = :none, color = color, markerstrokecolor=color, xlim = xlim, ylim = (0,ylimW), linestyle = style, framestyle = :box, yrotation = 45)
        #Plots.plot!(p, xs, means, yerr = (uppers, lowers), size = (300,300), linewidth = 3, label = leg, color = color, markerstrokecolor=color, xlim = xlim, linestyle = style)
    end

    #Plots.hline!(p, [maxWs[i]], linewidth = 1, color = :black, linestyle = :dash, label = "max fitness")
    return p

end

plotcolumn = function(df, delta, ylimWs)

    #interate over each E
    ps = []
    for E in 1:3
        p = plotselectedE(df, E, delta, ylimWs[E])
        push!(ps, p)
    end

    fig = Plots.plot(ps[1], ps[2], ps[3], layout=(3, 1), size = (400,800) ,left_margin = 10Plots.mm, xguidefontsize = 9)
    fig
    return ps, fig
end

p_base = plotcolumn(df, "baseline", [70, 400, 10000])[2]
savefig(p_base, "results/plots/selectedE_baseline.pdf")

ps = plotcolumn(df, "Nstar", [70, 400, 10000])[2]
savefig(ps, "results/plots/selectedE_Nstar.pdf")

ps = plotcolumn(df, "mutmag", [250, 4000, 10000])[2]
savefig(ps, "results/plots/selectedE_mutmag.pdf")

#now plot supplemental selected E
#first node mutation

ps = plotcolumn(df, "nodem_adjusted", [70, 400, 10000])[2]
savefig(ps, "results/plots/selectedE_nodem.pdf")

#now germline mutation

p1 = plotcolumn(df, "germm0p001", [70, 400, 10000])[2]
p2 = plotcolumn(df, "germm0p01", [70, 400, 10000])[2]
p3 = plotcolumn(df, "germm0p1", [70, 400, 10000])[2]

my_layout = @layout [
    a  b  c d;   # Row 1
]
p_combined = Plots.plot(p_base, p1, p2, p3, layout = my_layout, size = (2400,1200), xlabel = "Network size (number of mutating parameters)", labelfontsize = 10, ylabel = "Expected Fitness", margin = 10mm, leftmargin = 12mm)

savefig(p_combined, "results/plots/selectedE_germlinem.pdf")
