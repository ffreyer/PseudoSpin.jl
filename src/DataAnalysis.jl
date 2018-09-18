################################################################################
#### Binning Anaylsis
################################################################################


# This is a Node in the Binning Analysis tree, that averages two values. There
# is one of these for each binning level. Since when two values should be
# compressed, this is done immediately, so that only one value needs to be saved.
# switch indicates whether value should be written to or averaging should happen.
mutable struct BinningCompressor
    value::Float64
    switch::UInt8
end


mutable struct BinnerA
    # list of Compressors, one per level
    compressors::Vector{BinningCompressor}

    # currently binned values
    output::Vector{Float64}
    # output array has a constant size. i indicates the first index that can
    # be (over-)written
    i::UInt32
    # max length of output
    max_output_size::UInt32

    # sum(x) for all values on a given lvl (excluding the upper most level)
    x_sum::Vector{Float64}
    # sum(x.^2) for all values on a given lvl (excluding the upper most level)
    x2_sum::Vector{Float64}
    # number of values that are summed on a given lvl (excluding the upper most level)
    count::Vector{Int64}
end


"""
    BinnerA(min_output_size)

Creates a new Binning Analysis which keeps min_output_size values around after
binning. Returns a Binning Analysis object. Use push! to add values.
"""
function BinnerA(min_output_size::Integer)
    BinnerA(
        BinningCompressor[],
        Array{Float64}(2 * min_output_size),
        UInt32(1),
        UInt32(2 * min_output_size),
        Int64[],
        Float64[],
        Float64[]
    )
end


"""
    push!(BinningAnalysis, value)

Pushes a new value into the Binning Analysis.
"""
function push!(B::BinnerA, value::Float64)
    # first set of values
    if isempty(B.compressors)
        # fill output if possible
        if B.i < B.max_output_size
            B.output[B.i] = value
            B.i += UInt32(1)
            return nothing

        # once output is full, add compressor and drop one level
        else
            push!(B.compressors, BinningCompressor(0., UInt8(0)))
            push!(B.x_sum, 0.)
            push!(B.x2_sum, 0.)
            push!(B.count, 0)

            B.output[B.i] = value
            B.i = 1
            for x in B.output
                push!(B, 1, x)
            end

            return nothing
        end

    else
        # push value in compressor stack
        if B.i <= B.max_output_size
            push!(B, 1, value)
            return nothing

        # if output is full, add new level and drop down
        else
            push!(B.compressors, BinningCompressor(0., UInt8(0)))
            push!(B.x_sum, 0.)
            push!(B.x2_sum, 0.)
            push!(B.count, 0)

            B.i = 1
            for x in B.output
                push!(B, length(B.compressors), x)
            end
            push!(B, 1, value)

            return nothing
        end
    end

    return nothing
end


# recursion, back-end function
function push!(B::BinnerA, lvl::Int64, value::Float64)
    C = B.compressors[lvl]

    # any value propagating through this function is new to lvl. Therefore we
    # add it to the sums. Note that values pushed to the output arrays are not
    # added here until the array drops to the next level. (New compressors are
    # added)
    B.x_sum[lvl] += value
    B.x2_sum[lvl] += value^2
    B.count[lvl] += 1

    if C.switch == 0
        # Compressor has space -> save value
        C.value = value
        C.switch = 1
        return nothing
    else
        # Do averaging
        if lvl == length(B.compressors)
            # No more propagation possible -> push to output array
            B.output[B.i] = 0.5 * (C.value + value)
            B.i += 1
            C.switch = 0
            return nothing
        else
            # propagate to next lvl
            C.switch = 0
            push!(B, lvl+1, 0.5 * (C.value + value))
            return nothing
        end
    end

    println("How did you get here? recursive push!")
    return nothing
end


# NOTE:
# This is actually Var(x)/N
# This definition is more practical for the binning analysis and error
# calculation since we need the standard error sqrt(Var(x) / N)
"""
    var(BinningAnalysis[, lvl])

Calculates the **variance/N** of a given level in the Binning Analysis.
"""
function var(B::BinnerA, lvl::Int64=1_000_000)
    # The upper most lvl needs to be calculated explicitly
    if lvl < length(B.count)
        return ((B.x2_sum[lvl+1] / B.count[lvl+1]) - (B.x_sum[lvl+1] / B.count[lvl+1])^2) / (B.count[lvl+1] - 1)
    else
        return ((sum(B.output[1:B.i-1].^2) / (B.i-1)) - (sum(B.output[1:B.i-1]) / (B.i-1))^2) / (B.i - 2)
    end
end


"""
    mean(BinningAnalysis[, lvl])

Calculates the mean for a given level in the Binning Analysis
"""
function mean(B::BinnerA, lvl::Int64=1_000_000)
    if lvl < length(B.count)
        return B.x_sum[lvl] / B.count[lvl]
    else
        return mean(B.output[1:B.i-1])
    end
end


"""
    tau(BinningAnalysis, lvl)

Calculates tau...
"""
function tau(B::BinnerA, lvl::Int64)
    var_0 = var(B, 0)
    var_l = var(B, lvl)
    0.5 * (var_l / var_0 - 1)
end


################################################################################
#### Histogram binning
################################################################################

# This is kept very simplistic. The width of bins is fixed and bins are created
# when needed. bin width should be somewhere in the range of 0.001 - 0.0001.

mutable struct BinnerH
    data::Dict{Int64, Int64}
    bin_width::Float64

    """
        BinnerH(bin_width)

    Creates a new Histogram Binner with a given bin width.
    """
    BinnerH(bin_width::Float64) = new(Dict{Int64, UInt32}(), bin_width)
end

"""
    push!(HistogramBinner, x)

Pushes x into the Histogram Binner.
"""
function push!(B::BinnerH, x::Float64)
    id = floor(Int64, x / B.bin_width)

    if haskey(B.data, id)
        B.data[id] += 1
    else
        push!(B.data, id => 1)
    end
end


################################################################################
#### Jackknife
################################################################################


# This requires the full xs arrays :(
"""
    jackknife(function, args...)

Performs a jackknife analysis given a function and its arguments. Returns the
mean and standard error.

Example:
    xs = rand(1000)
    ys = rand(1000)
    value, error = jackknife((x, y) -> x*y, xs, ys)
"""
function jackknife(f::Function, args...)
    N = length(args[1])
    x0s = map(sum, args) # (xs, ...) -> (sum xs, ...)
    x_avs_s = map((xs, x0) -> [(x0 - xs[i]) / (N-1) for i in 1:N], args, x0s)
    x0s = map(x -> x / N, x0s)

    y0 = map(f, x0s...)
    y_avs = map(f, x_avs_s...)

    # outputs mean and standard deviation of the distribution of means (which
    # is the standard error of ys)
    y0, sqrt((N-1) / N * sum((y_avs .- y0).^2))
end


################################################################################
#### Thermalization time-averaging
################################################################################

# Take some data and compute local averages
#

mutable struct Compressor
    compression::Int64
    invN::Float64
    count::Int64
    sum::Float64
    values::Vector{Float64}
end


function Compressor(compression::Int64 = 1_000)
    Compressor(
        compression,
        1.0 / compression,
        0,
        0.0,
        Float64[]
    )
end

function push!(c::Compressor, value::Float64)
    if c.count == c.compression
        push!(c.values, c.sum * c.invN)
        c.sum = value
        c.count = 1
    else
        c.sum += value
        c.count += 1
    end
    nothing
end

function output(c::Compressor)
    if c.count != 0
        push!(c.values, c.sum / c.count)
        c.count = 0
        c.sum = 0.0
    end
    c.values
end
