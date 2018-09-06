@everywhere begin
using PseudoSpin
Point3 = PseudoSpin.Point3
myargs = remotecall_fetch(() -> ARGS, 1)
println("I got $myargs on $(myid())")

kwargs = Dict{Symbol, Any}()

argfile = open(myargs[1], "r")

for argline in eachline(argfile)
    isempty(argline) && continue
    args = split(chomp(argline), "\t")
    length(args) == 1 && push!(args, "")
    if startswith(args[1], "#")
        continue
    elseif startswith(args[1], "neighbors") || startswith(args[1], "neighbours")
        push!(kwargs, :neighbor_search_depth => parse(Int64, args[2]))
    elseif startswith(args[1], "do_paths") || startswith(args[1], "do paths")
        push!(kwargs, :do_paths => parse(Bool, args[2]))
    elseif startswith(args[1], "L")
        push!(kwargs, :L => parse(Int64, args[2]))
    elseif startswith(args[1], "path")
        push!(kwargs, :path => String(args[2]))  # Paths should not have \t in them
    elseif startswith(args[1], "folder")
        push!(kwargs, :folder => String(args[2]))
    elseif startswith(args[1], "filename") || startswith(args[1], "name")
        push!(kwargs, :filename => String(args[2]))
    elseif (args[1] == "T") || (args[1] == "temperature")
        push!(kwargs, :T => parse(Float64, args[2]))
    elseif (args[1] == "Ts") || (args[1] == "temperatures")
        push!(kwargs, :Ts => map(x -> parse(Float64, x), args[2:end]))
    elseif args[1] == "Js" # How do I do this?
        if length(args) < 9; throw(ErrorException("Failed to read Js")) end
        # push!(kwargs, :Js => [
        #     (parse(Float64, args[2]), parse(Float64, args[3])),
        #     (parse(Float64, args[4]), parse(Float64, args[5])),
        #     (parse(Float64, args[6]), parse(Float64, args[7])),
        #     (parse(Float64, args[8]), parse(Float64, args[9]))
        # ])
        push!(kwargs, :J1s => (parse(Float64, args[2]), parse(Float64, args[3])))
        push!(kwargs, :J2s => (parse(Float64, args[4]), parse(Float64, args[5])))
        push!(kwargs, :K => parse(Float64, args[6]) * parse(Float64, args[9]))
        @assert parse(Float64, args[7]) == 0.0 "General 4-spin term not supported"
        @assert parse(Float64, args[8]) == 0.0 "General 4-spin term not supported"
    elseif args[1] == "J1"
        push!(kwargs, :J1 => parse(Float64, args[2]))
    elseif args[1] == "J2"
        push!(kwargs, :J2 => parse(Float64, args[2]))
    elseif args[1] == "J3"
        push!(kwargs, :J3 => parse(Float64, args[2]))
    elseif args[1] == "J1s"
        push!(kwargs, :J1s => (parse(Float64, args[2]), parse(Float64, args[3])))
    elseif args[1] == "J2s"
        push!(kwargs, :J2s => (parse(Float64, args[2]), parse(Float64, args[3])))
    elseif args[1] == "J3s"
        push!(kwargs, :J3s => (parse(Float64, args[2]), parse(Float64, args[3])))
    elseif args[1] == "K"
        push!(kwargs, :K => parse(Float64, args[2]))
    elseif args[1] == "lambda"
        push!(kwargs, :lambda => parse(Float64, args[2]))
    elseif startswith(args[1], "freeze_temp") || startswith(args[1], "freeze temp")
        push!(kwargs, :Freeze_temperature => parse(Float64, args[2]))
    elseif startswith(args[1], "N_switch")
        push!(kwargs, :N_switch => parse(Int64, args[2]))
    elseif startswith(args[1], "TH") || startswith(args[1], "thermalization")
        push!(kwargs, :TH_sweeps => parse(Int64, args[2]))
    elseif startswith(args[1], "ME") || startswith(args[1], "measurement")
        push!(kwargs, :ME_sweeps => parse(Int64, args[2]))
    elseif startswith(args[1], "h")
        push!(kwargs, :h => Point3(map(x -> parse(Float64, x), args[2:end])))
    elseif args[1] == "g"
        push!(kwargs, :g => parse(Float64, args[2]))
    elseif args[1] == "zeta"
        push!(kwargs, :zeta => parse(Float64, args[2]))

    elseif args[1] in ["TGen_method", "TH_method"]
        push!(kwargs, :TGen_method => args[2] |> Meta.parse |> eval)
    elseif args[1] == "thermalizer_method"
        push!(kwargs, :thermalizer_method => args[2] |> Meta.parse |> eval)
    elseif args[1] == "batch_size"
        push!(kwargs, :batch_size => parse(Int64, args[2]))
    elseif args[1] == "adaptive_sample_size"
        push!(kwargs, :adaptive_sample_size => parse(Int64, args[2]))
    elseif args[1] == "skip"
        push!(kwargs, :skip => parse(Int64, args[2]))

    elseif startswith(args[1], "spins") || startswith(args[1],"spin")
        println("TODO: custom initial spin vectors")
        exit(-1)
        # push!(spins, Point3{Float64}(map(x -> parse(Float64, x), args[2:end])))
    elseif isempty(args[1])
        continue
    else
        println("Did not recognise \"", argline, "\"")
    end
end

close(argfile)

simulate!(; kwargs...)
end
