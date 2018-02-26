tic()
include("Lattice.jl")

# TODO
# unify simulate! to work with single temp input


################################################################################
#### Defaults
################################################################################

neighbor_search_depth = 2
do_paths = true
L = 6

path = ""           # should only end with "/" if folder is ""
folder = ""         # should not end with "/"
filename = ""

Ts = 1.0            # Float64 or Vector{Float64}
# Js = [(-1.0, 0.0), (0.120, 0.0), (-0.1, 0.0), (0.0, 1.0)]
J1 = -1.0
J2 = 0.12
K = -0.1
lambda = 0.0
make_Js = true
Js = [(0., 0.), (0., 0.), (0., 0.), (0., 0.)]
h = Point3(0.)

TH_sweeps = Int64(2e6)
ME_sweeps = Int64(5e6)
Freeze_temperature = 1.5 * maximum(Ts)
N_switch = div(TH_sweeps, 2)

spins = Point3{Float64}[]


################################################################################
#### Argument Parsing
################################################################################
argfile = open(ARGS[1], "r")

for argline in eachline(argfile)
    args = split(chomp(argline), "\t")
    if startswith(args[1], "#")
        continue
    elseif startswith(args[1], "neighbors") || startswith(args[1], "neighbours")
        neighbor_search_depth = parse(Int64, args[2])
    elseif startswith(args[1], "do_paths") || startswith(args[1], "do paths")
        do_paths = parse(Bool, args[2])
    elseif startswith(args[1], "L")
        L = parse(Int64, args[2])
    elseif startswith(args[1], "path")
        path = String(args[2])  # Paths should not have \t in them
    elseif startswith(args[1], "folder")
        folder = String(args[2])
    elseif startswith(args[1], "filename") || startswith(args[1], "name")
        filename = String(args[2])
    elseif (args[1] == "T") || (args[1] == "temperature")
        Ts = parse(Float64, args[2])
    elseif (args[1] == "Ts") || (args[1] == "temperatures")
        Ts = map(x -> parse(Float64, x), args[2:end])
    elseif args[1] == "Js" # How do I do this?
        if length(args) < 9; throw(ErrorException("Failed to read Js")) end
        make_Js = false
        Js = [
            (parse(Float64, args[2]), parse(Float64, args[3])),
            (parse(Float64, args[4]), parse(Float64, args[5])),
            (parse(Float64, args[6]), parse(Float64, args[7])),
            (parse(Float64, args[8]), parse(Float64, args[9]))
        ]
    elseif args[1] == "J1"
        J1 = parse(Float64, args[2])
    elseif args[1] == "J2"
        J2 = parse(Float64, args[2])
    elseif args[1] == "K"
        K = parse(Float64, args[2])
    elseif args[1] == "lambda"
        lambda = parse(Float64, args[2])
    elseif startswith(args[1], "freeze_temp") || startswith(args[1], "freeze temp")
        Freeze_temperature = parse(Float64, args[2])
    elseif startswith(args[1], "N_switch")
        N_switch = parse(Int64, args[2])
    elseif startswith(args[1], "TH") || startswith(args[1], "thermalization")
        TH_sweeps = parse(Int64, args[2])
    elseif startswith(args[1], "ME") || startswith(args[1], "measurement")
        ME_sweeps = parse(Int64, args[2])
    elseif startswith(args[1], "h")
        h = Point3(map(x -> parse(Float64, x), args[2:end]))
    elseif startswith(args[1], "spins") || startswith(args[1],"spin")
        push!(spins, Point3{Float64}(map(x -> parse(Float64, x), args[2:end])))
    end
end

close(argfile)

# safety and additional parsing stuff
if folder == ""
    if !endswith(path, "/")
        path = path * "/"
    end
else
    if !endswith(folder, "/")
        folder = folder * "/"
    end
    if startswith(folder, "/")
        folder = folder[2:end]
    end
end

if make_Js
    Js = [
        (J1, J1*lambda),
        (J2, J2*lambda),
        (K, 0.0),
        (0.0, 1.0)
    ]
end


################################################################################
#### Start simulation
################################################################################

# TODO: fully implement the input of a specifc spin configuration

g = RGraph(diamond("A"), neighbor_search_depth)
if do_paths; findPaths!(g) end
sim, flat_index, lattice_indices = Basisfill(g, L)
if isempty(spins)
    spins = rand_spin(sim.N_nodes)
else
    @assert length(spins) == sim.N_nodes "Spin configuration does not fit System size!"
    spins = deepcopy(spins) # defrag spins?
end
# println(spins)

simulate!(
    sim,
    spins,
    L,
    path * folder,
    filename,
    Ts,
    Js,
    Freezer(TH_sweeps, Freeze_temperature, N_switch=N_switch),
    ME_sweeps,
    h
)
println("Runtime: ", toq())
