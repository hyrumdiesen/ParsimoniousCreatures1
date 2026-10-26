using JLD2
using HDF5
using Dates
using JSON3


#I think it's easiest to just convert .jld2 to named dictionaries in .json format for the selectedE data
#for the raw heatmap data, best to just make them into a giant tensor and save as .h5

#load in a .jld2 file 
@load "/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_avg_r25_s1p0.jld2" 

data_dict = load("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_avg_r25_s1p0.jld2")
data_dict.keys

first_object = jldopen("filename.jld2", "r") do file
    # Get the very first key string from the keys iterator
    first_key = first(keys(file)) 
    
    println("Loading the first key: ", first_key)
    return file[first_key]
end

alldata[26][2]

data = first(load("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_avg_r25_s1p0.jld2"))[2]

#so it's a 26 long vector of two tensors, I think we can turn that into one big tensor

sizes_ssm_avg = string(collect(1:26), "that's annoying. This isn't number of modules, but multiples of 3 parameters")
sizes_ssm_min = string(collect(1:26))
sizes_ffvar1 = string(collect(1:8))
sizes_ffvar2 = string([1,2,3,4,5,6,8,10,12,15,20,25,30,40])

converttotensor = function(filename, nmods_tested)



    filename2 = filename*".jld2"

    #@load filename2 alldata

    @load filename2 allthedata_min_r25_s1p0
    alldata = allthedata_min_r25_s1p0

    #@load filename2 test_allsizes_r12_s1p0_avg
    #alldata = test_allsizes_r12_s1p0_avg

    nsizes = length(alldata)

    tensordata = Array{Float64}(undef, nsizes, 2, 127,127,50)

    filename3 = split(filename, "/")[end]
    for i in 1:nsizes
        tensordata[i, 1, :,:,:] .= alldata[i][1]
        tensordata[i, 2, :,:,:] .= alldata[i][2][:,:,2:end]
    end
    #return tensordata

    dim_labels = ["network size", "first element is optimized r, second element is fitness", "organism death rate", "mutation rate", "generation"]

    h5open(filename*".h5", "w") do file
        # Save the tensor
        file["my_tensor"] = tensordata

        #create_dataset(file, filename3, tensordata, chunk=(nsizes, 2, 127, 127, 50))

        
        # Attach a "dim_labels" attribute to the dataset
        attributes(file["my_tensor"])["dim_labels"] = dim_labels
        attributes(file["my_tensor"])["number of modules tested"] = nmods_tested
        attributes(file["my_tensor"])["other metadata"] = filename3
    end

end

converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_avg_r12_s0p1", sizes_ssm_avg)
converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_avg_r12_s1p0", sizes_ssm_avg)
converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_avg_r25_s1p0", sizes_ssm_avg)

converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_min_r12_s0p1", sizes_ssm_avg)
converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_min_r12_s1p0", sizes_ssm_avg)
converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_min_r25_s1p0", sizes_ssm_avg)



converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_ffvar1_r12_s0p1", sizes_ffvar1)
converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_ffvar1_r12_s1p0", sizes_ffvar1)
converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_ffvar1_r25_s1p0", sizes_ffvar1)
converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_ffvar2_r12_s0p1", sizes_ffvar2)
converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_ffvar2_r12_s1p0", sizes_ffvar2)
converttotensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_ffvar2_r25_s1p0", sizes_ffvar2)

#I guess we should have tested first if this worked

#try to load .h5 file

loaded_data = h5open("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/heatmaps_raw/alldata_ffvar2_r25_s1p0.h5", "r") do file
    read(file["my_tensor"])
end

loaded_data


#now let's do the same thing with the processed data, turn it into a bigger tensor rather than a vector of tensors

#first pring in a processed data .jld2 file
@load "/Users/hyrumd/Documents/consolidatedcodepaperone/research/ode/results/heatmap_data/heatmaps_processed/processed_r25_s1p0_paper1.jld2"

using Plots
Plots.heatmap(r25_s1p0_paper1[1][1][:,:,1])
#first dimension is small, medium large
#second dimension is 1: bootstrap probability of best network or 2. ideal sizes
#third dimension is each network type

r25_s1p0_paper1[2]


convert_processed_totensor = function(filename)

    filename2 = filename*".jld2"

    #@load filename2 alldata

    @load filename2 r12_s0p1_paper1
    alldata = r12_s0p1_paper1

    tensordata = Array{Float64}(undef, 3, 2, 127,127,3)

    filename3 = split(filename, "/")[end]
    for i in 1:3
        tensordata[i, 1, :,:,:] .= alldata[i][1]
        tensordata[i, 2, :,:,:] .= alldata[i][2]
    end
    #return tensordata

    dim_labels = ["first dimension is 1. small networks, 2. medium networks, or 3. large networks", "second dimension is 1: bootstrap probability of best network or 2. ideal sizes", "organism death rate", "mutation rate", "network type"]

    h5open(filename*".h5", "w") do file
        # Save the tensor
        file["my_tensor"] = tensordata

        #create_dataset(file, filename3, tensordata, chunk=(nsizes, 2, 127, 127, 50))

        
        # Attach a "dim_labels" attribute to the dataset
        attributes(file["my_tensor"])["dim_labels"] = dim_labels
        attributes(file["my_tensor"])["other metadata"] = filename3
    end

end

convert_processed_totensor("/Users/hyrumd/Documents/consolidatedcodepaperone/research/ode/results/heatmap_data/heatmaps_processed/processed_r12_s0p1_paper1")

@load "/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/plotting_data/hybriddata/E1_BASELINE_HYBRIDS.jld2"

E1_BASELINE_HYBRIDS[2]

#the solution is to create a dataframe, then we will flatten the data and save as .csv
#but that's alot of work, I think instead we will just put the fitnesses as columns on the dataframes and call it good.
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

using DataFrames
convert_hybrids = function(filename)

    alldata = first(load(filename))[2]

    filename = split(filename, "/")[end]
    E = parse(Int, split(filename, "_")[1][2])

    delta = split(filename, "_")[2]

    if E == 1
        m = 0.01
        deltaO = 0.003
    elseif E == 2
        m = 0.003
        deltaO = 0.0001
    elseif E == 3
        m = 0.00005
        deltaO = 0.0001
    end

    if delta == "BASELINE"
        mutmag = 1.0
    elseif delta == "deltasigma"
        mutmag = 0.1
    elseif delta == "deltasigma2"
        mutmag = sqrt(0.1)
    elseif delta == "deltasigma3"
        mutmag = 0.5
    end

    dfrows = []

    for i in 1:length(TABLE2)

        data = alldata[i]

        hybridvect = TABLE2[i]
        nnodes,nedges = hybridnnodesnedges(hybridvect)

        dfrow = DataFrame(
                network_type = [string(hybridvect)],
                delta = [delta],
                E = [E],
                nedges = [nedges],
                nnodes = [nnodes],
                Nstar = [12],
                Ninit = [1],
                dO = [deltaO],
                m = [m],
                mutmag = [mutmag],
                gammam = [0.0],
                optr = [data[2][data[3]]],
                germlinem = [0.0],
                objreplicates = [30],
                fitnessreplicates = [length(data[4])],
                meanlogW = [data[5][1]],
                meanlogWupper = [data[5][3]],
                meanlogWlower = [data[5][2]],
                meanW = [data[6][1]],
                meanWupper = [data[6][3]],
                meanWlower = [data[6][2]],
                logWs = [data[4]]
            )

        push!(dfrows, dfrow)

    end

    df = reduce(vcat, dfrows) 

    return df

end

convert_hybrids = function(filename)

    alldata = first(load(filename))[2]

    filename = split(filename, "/")[end]
    E = parse(Int, split(filename, "_")[1][2])

    delta = split(filename, "_")[2]

    if E == 1
        m = 0.01
        deltaO = 0.003
    elseif E == 2
        m = 0.003
        deltaO = 0.0001
    elseif E == 3
        m = 0.00005
        deltaO = 0.0001
    end

    if delta == "BASELINE"
        mutmag = 1.0
    elseif delta == "deltasigma"
        mutmag = 0.1
    elseif delta == "deltasigma2"
        mutmag = sqrt(0.1)
    elseif delta == "deltasigma3"
        mutmag = 0.5
    end

    dfrows = []

    for i in 1:length(TABLE2)

        data = alldata[i]

        hybridvect = TABLE2[i]
        nnodes,nedges = hybridnnodesnedges(hybridvect)

        dfrow = DataFrame(
                network_type = [string(hybridvect)],
                delta = [delta],
                E = [E],
                nedges = [nedges],
                nnodes = [nnodes],
                Nstar = [12],
                Ninit = [1],
                dO = [deltaO],
                m = [m],
                mutmag = [mutmag],
                gammam = [0.0],
                optr = [data[2][data[3]]],
                germlinem = [0.0],
                objreplicates = [30],
                fitnessreplicates = [length(data[4])],
                meanlogW = [data[5][1]],
                meanlogWupper = [data[5][3]],
                meanlogWlower = [data[5][2]],
                meanW = [data[6][1]],
                meanWupper = [data[6][3]],
                meanWlower = [data[6][2]],
                best_y = [data[1]],
                best_x = [data[2]],
                chosen = [data[3]],
                logWs = [data[4]]                

            )

        push!(dfrows, dfrow)

    end

    df = reduce(vcat, dfrows) 

    return df

end

test = convert_hybrids("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/plotting_data/hybriddata/E1_BASELINE_HYBRIDS.jld2")

using CSV
CSV.write("test.csv", test)
#now let's see if we can flatten the dataframe.


#now do it for every file in the folder, then combine them together into one mega dataframe
folder_path = "/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/plotting_data/hybriddata"
csv_files = filter(f -> endswith(f, ".jld2"), readdir(folder_path, join=true))

subdfs = []
for i in 1:length(csv_files)
    subdf = convert_hybrids(csv_files[i])
    push!(subdfs, subdf)
end

masterdf = reduce(vcat, subdfs)

CSV.write("allhybrids_allenvironments.csv", masterdf)



for i in 1:length(csv_files)
    df = convert_hybrids(csv_files[i])
    filename = csv_files[i][1:(end-5)]
    filename = filename*".csv"
    CSV.write(filename, df)
end

#now edit selectedE .csv files

function find_target_folders(root_path::String)
    target_names = Set(["E1", "E2", "E3"])
    matching_paths = String[]

    # walkdir returns (current_root, directories, files) for every subfolder
    for (root, dirs, _) in walkdir(root_path)
        for dir in dirs
            if dir in target_names
                # Combine the current root with the folder name to get the full path
                push!(matching_paths, joinpath(root, dir))
            end
        end
    end

    return matching_paths
end

addnewdata = function(subdf, jld2var)

    #make a new dataframe that has all the extra data we want, then hcat it together
    #with the old dataframe

    dfrows = []
    for i in 1:length(jld2var)
        data = jld2var[i]
        #fs = jld2var[i][4]

        dfrow = DataFrame(best_y = [data[1]],
                best_x = [data[2]],
                chosen = [data[3]],
                logWs = [data[4]])
        push!(dfrows, dfrow)


    end
    newdf = reduce(vcat, dfrows)
    newdf = hcat(subdf, newdf)

    return newdf
end

addnewdata_all = function(folder_path)

    jld2files = filter(f -> endswith(f, ".jld2"), readdir(folder_path, join=true))
    csvfile = filter(f -> endswith(f, ".csv"), readdir(folder_path, join=true))

    df = DataFrame(CSV.File(csvfile))
    #return folder_path*csvfile[1]
    df2 = groupby(df, :network_type)

    df3 = DataFrame()
    #return df2
    for i in 1:length(jld2files)
        if contains(jld2files[i], "ssm_avg")
            @load jld2files[i] ssm_avg_data
            df3 = vcat(df3, addnewdata(df2[1], ssm_avg_data))
        elseif contains(jld2files[i], "ssm_min")
            @load jld2files[i] ssm_min_data
            df3 = vcat(df3, addnewdata(df2[2], ssm_min_data))
        elseif contains(jld2files[i], "ff_var1")
            @load jld2files[i] ff_var1_data
            df3 = vcat(df3, addnewdata(df2[3], ff_var1_data))
        elseif contains(jld2files[i], "ff_var2")
            @load jld2files[i] ff_var2_data
            df3 = vcat(df3, addnewdata(df2[4], ff_var2_data))
        end
    end

    CSV.write(csvfile[1], DataFrame(df3))

    #return df3
end

addnewdata_all()

folderpaths = find_target_folders("/Users/hyrumd/Documents/consolidatedcodepaperone/research/analytic/data/raw_data/selectedE_raw")

test = addnewdata_all(folderpaths[1])

for i in folderpaths
    addnewdata_all(i)
end

#now we should combine these dataframes together.

subdfs = []
for i in folderpaths
    csvfile = filter(f -> endswith(f, ".csv"), readdir(i, join=true))

    subdf = DataFrame(CSV.File(csvfile))
    push!(subdfs, subdf)
end

masterdf = reduce(vcat, subdfs)

CSV.write("allE_allnetworks_alldeltas.csv", masterdf)