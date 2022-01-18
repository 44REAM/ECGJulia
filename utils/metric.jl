
using Statistics
include("dsp.jl")

function sdnn(nn_interval::Vector{<:Real})

    sdnn = std(nn_interval)
    return sdnn
end


function sdann_asdnn(nn_interval::Vector{<:Real}, index_array::Vector{<:Integer},
    duration::Integer, min_number_of_interval::Integer = 3)

    if length(nn_interval) < min_number_of_interval
        return nothing, nothing
    end

    segment_index = get_segment_index(index_array, duration)

    average_array::Vector{Float64} = []
    std_array::Vector{Float64} = []

    start_idx = 1
    for end_idx in segment_index
        landmark = @view nn_interval[start_idx:end_idx]

        if length(landmark) >= min_number_of_interval
            push!(average_array, mean(landmark))
            push!(std_array, std(landmark))
        end
        start_idx = end_idx
    end

    if length(average_array) >= min_number_of_interval
        
        sdann = std(average_array)
        asdnn = mean(std_array)
        return sdann, asdnn
    end

    return nothing, nothing
end