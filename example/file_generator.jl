################################################################################
# Number of parameter files/simulations.
# If N = -1, the script will figure this out...
# N = -1

Ts = [
	0.001, 
	collect(0.01:0.01:0.3)..., 
	collect(0.32:0.02:0.7)..., 
	collect(0.73:0.03:1.0)..., 
	collect(1.05:0.05:2.0)...
]
NT = length(Ts)
N = 6NT

input_path = "/home/ffreyer/simulation/param/batch1/"
param_filename = ["$i.param" for i in 0:N-1]
bash_filename = "run_sim1.sh"


# Set all the none-default things
param = Dict{Symbol, Any}(
    # NOTE output_path should be static
    # NOTE should end with /
    :path => "/home/ffreyer/simulation/output/g_phase_diagram/batch1/",
    :folder => [@sprintf(
        "J1_%0.2f_g_%0.2f_L8", x, y
    ) for _ in 1:NT, x in [-0.4, 0.0, 0.4], y in [-0.4, 0.0]][:],

    :L => 8,
    :T => [x for x in Ts, _ in 1:6][:],
    :J1 => [x for _ in 1:NT, x in [-0.4, 0.0, 0.4], __ in 1:2][:],
    #:J2 => 0.,
    #:K => 0.,
    #:lambda => 0.,
    # :h => (0., 0., 0.),
    :g => [x for _ in 1:NT, __ in 1:3, x in [-0.4, 0.0]][:],

    # :TH => 2_000_000,
    :freeze_temp => 2.5,
    # :ME => 5_000_000,

    :filename => ["T$i" for i in 1:NT, _ in 1:6][:],


    :neighbors => 1,
    :do_paths => false,
    # :spins => rand_spin(2^3),
    # :Ts => [T],           # NOTE: You probably don't want this
    # :Js => ((J1, lambda*J1), (J2, lambda*J2), (K, 0.), (0., 1.)),
    # :N_switch => div(TH_sweeps, 2),
)


#=
# Reminder:

> [x for _ in 1:2, x in 1:3][:] creates
[1, 1, 2, 2, 3, 3]

#_lambda_%0.2f   , lambda
replace(@sprintf("J2_%0.2f_K_%0.2f_L%i/", J2, K, L), ".", "_")
=#
################################################################################


# Create paths
!isdir(input_path) && mkdir(input_path)
!isdir(param[:path]) && mkdir(param[:path])
for folder in param[:folder]
    if !isdir(param[:path] * folder)
        mkdir(param[:path] * folder)
    end
end


# Checking if array lengths match
for (k, v) in param
    if typeof(v) == Vector
        if N == -1
            N = length(v)
        else
            if N != length(v)
                println("Length of entry '$k' does not match previous lengths. Abort.")
                exit(-1)
            end
        end
    end
end


# write that splats Tuples correctly (hopefully)
function _write(file, var)
    if typeof(var) <: Tuple{4, Tuple} # Js
        write(file, mapreduce(
            inner -> map(string, (x, y) -> "$x\t$y", inner),
            (x, y) -> "$x\t$y",
            var
        ))
    elseif typeof(var) <: Tuple
        write(file, map(string, (x, y) -> "$x\t$y", var))
    else
        write(file, string(var))
    end
    nothing
end


# generate param files
for i in 1:N
    open(input_path * param_filename[i], "w") do file
        for (k, v) in param
            write(file, string(k) * "\t")
            if typeof(v) <: Vector
                _write(file, v[i])
            else
                _write(file, v)
            end
            write(file, "\n")
        end
    end
end





# "2048mb"  4096 6144
# NOTE -----------------------------------
memory = "4096mb" #"2048mb"
time = "12:00:00" #6
# NOTE -----------------------------------

open(bash_filename, "w") do file		# <--- CHANGE
    write(file, "#!/bin/bash -l\n")
    write(file, "#SBATCH --mem=" * memory * "\n")
    write(file, "#SBATCH --time=" * time * "\n")
    write(file, "#SBATCH --account=AG-Trebst\n")
    write(file, "#SBATCH --array=0-" * string(N-1) * "\n")
    write(file, "\n")

    write(file, "source /projects/ag‐trebst/julia/0.6/load\n")
    write(file, "\n")

    write(file, "input_path=\"" * string(input_path) * "\"\n")
    write(file, "output_path=\"" * string(param[:path]) * "\"\n")
    write(file, "\n")

    write(file, "input=\"\${input_path}\${SLURM_ARRAY_TASK_ID}.param\"\n")
    write(file, "output=\"\${output_path}\${SLURM_ARRAY_TASK_ID}.out\"\n")
    write(file, "\n")

    write(file, "julia --optimize=3 mainf.jl \${input} > \${output}\n")
end

exit()



# DEFAULTS for reference
#=
begin
    path = ""
    folder = ""
    filename = ""

    neighbor_search_depth = 2
    do_paths = true
    L = 6
    spins = rand_spin(2^3)

    T = 1.0
    Ts = [T]
    J1 = 0.
    J2 = 0.
    K = 0.
    lambda = 0.
    Js = ((J1, lambda*J1), (J2, lambda*J2), (K, 0.), (0., 1.))
    h = (0., 0., 0.)
    g = 0.

    TH_sweeps = 2_000_000
    N_switch = div(TH_sweeps, 2)
    Freeze_temperature = 1.5 * maximum(Ts)
    ME_sweeps = 5_000_000
end
=#
