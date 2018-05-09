include("Lattice.jl")

# TODO
# unify simulate! to work with single temp input

# defaults
neighbor_search_depth = 2
do_paths = true
L = 6

path = ""           # should only end with "/" if folder is ""
folder = "test/"         # should not end with "/"
filename = "testfile"

Ts = 1.0            # Flaot64 or Vector{Float64}
# Js = [(-1.0, 0.0), (0.120, 0.0), (-0.1, 0.0), (0.0, 1.0)]
J1 = -1.0
J2 = 0.12
K = -0.1
lambda = 0.0
Freeze_temperature = 1.5 * maximum(Ts)

TH_sweeps = 25_000
ME_sweeps = 50_000


# argument parsing
i = 1
while i <= length(ARGS)
    if ARGS[i] == "-neighbors"
        neighbor_search_depth = parse(Int64, ARGS[i+1])
        i += 2
    elseif ARGS[i] == "-do_paths"
        do_paths = parse(Bool, ARGS[i+1])
        i += 2
    elseif ARGS[i] == "-L"
        L = parse(Int64, ARGS[i+1])
        i += 2

    elseif ARGS[i] == "-p"
        path = String(ARGS[i+1])
        i += 2
    elseif ARGS[i] == "-f"
        folder = String(ARGS[i+1])
        i += 1
    elseif ARGS[i] == "-n"
        filename = String(ARGS[i+1])
        i += 2

    elseif ARGS[i] == "-T"
        Ts = parse(Float64, ARGS[i+1])
        i += 2
    elseif ARGS[i] == "-Ts"
        Ts = Float64[]
        j = i+1
        while !startswith(ARGS[j], "-")
            push!(Ts, parse(Float64, ARGS[j]))
            j += 1
        end
        i = j
    elseif ARGS[i] == "-Js" # How do I do this?
        throw(ErrorException("Js doesn't work"))
    elseif ARGS[i] == "-J1"
        J1 = parse(Float64, ARGS[i+1])
        i += 2
    elseif ARGS[i] == "-J2"
        J2 = parse(Float64, ARGS[i+1])
        i += 2
    elseif ARGS[i] == "-K"
        K = parse(Float64, ARGS[i+1])
        i += 2
    elseif ARGS[i] == "-lambda"
        lambda = parse(Float64, ARGS[i+1])
        i += 2
    elseif ARGS[i] == "-Freeze_temp"
        Freeze_temperature = parse(Float64, ARGS[i+1])
        i += 2

    elseif ARGS[i] == "-TH"
        TH_sweeps = parse(Int64, ARGS[i+1])
        i += 2
    elseif ARGS[i] == "-ME"
        ME_sweeps = parse(Int64, ARGS[i+1])
        i += 2
    end
end


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

Js = [
    (J1, J1*lambda),
    (J2, J2*lambda),
    (K, 0.0),
    (0.0, 1.0)
]


# run simulation
g = RGraph(diamond("A"), neighbor_search_depth)
if do_paths; findPaths!(g) end
sim, flat_index, lattice_indices = Basisfill(g, L)
spins = rand_spin(sim.N_nodes)

simulate!(
    sim,
    spins,
    L,
    path * folder,
    filename,
    Ts,
    Js,
    Freezer(TH_sweeps, Freeze_temperature),
    ME_sweeps
)
