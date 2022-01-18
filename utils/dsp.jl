
include("gen_utils.jl")

using DSP

function mask_segment(signals::Vector{<:Real}, threshold::Real, edge::Integer)

    clip_binary = diff(abs.(signals) .>= threshold)

    #TODO slow on findall
    clipping_edges = findall(!=(0), clip_binary)

    question_segment = Vector{Tuple{Int,Int}}(undef, 0)

    n_edge::Int = length(clipping_edges)
    n_signals::Int = length(signals)
    if n_edge == 0
        return question_segment
    end

    start = 1

    if clip_binary[clipping_edges[1]] == -1
        push!(question_segment, (1, min(clipping_edges[1] + edge, n_signals)))
        start += 1
    end

    while start < n_edge
        push!(question_segment, (max(clipping_edges[start] + 1 - edge, 1), min(clipping_edges[start+1] + edge, n_signals)))
        start += 2
    end

    if clip_binary[clipping_edges[end]] == 1
        push!(question_segment, (max(clipping_edges[end] + 1 - edge, 1), n_signals))
    end

    return question_segment
end

function connect_segment(segments::Vector{<:Tuple{T,T}}) where T
    index = 1
    n_segment = length(segments)

    new_question_segments = Vector{Tuple{T,T}}(undef, 0)

    while index <= n_segment
        if index == n_segment
            push!(new_question_segments, segments[index])
            break
        end
        next_index = index + 1
        endpoint = segments[index][2]
        while endpoint >= segments[next_index][1] - 1
            endpoint = segments[next_index][2]
            next_index += 1
            if next_index > n_segment
                break
            end
        end

        push!(new_question_segments, (segments[index][1], segments[next_index-1][2]))
        index = next_index
    end
    return new_question_segments
end

function replace_zeros!(signals::Vector{<:Real}, segments::Vector{<:Tuple{<:Integer,<:Integer}})

    replace_value::Float64 = 0
    n_signals = length(signals)

    for segment in segments
        start_idx = segment[1]
        end_idx = segment[2]

        if start_idx-1 < 1 && end_idx+1 > n_signals
            replace_value = 0
        elseif start_idx-1 < 1
            replace_value = signals[end_idx+1]
        elseif end_idx+1 > n_signals
            replace_value = signals[start_idx-1]
        else
            replace_value = (signals[end_idx+1]+signals[start_idx-1])/2
        end

        for i in start_idx:end_idx
            signals[i] = replace_value
        end
        # signals[start_idx:end_idx] .= 0
    end
end

function ma(signals::Vector{<:Real}, windows::Integer, algo::String = "cumsum")
    return ma_cumsum(signals, windows)
end

function pad_end(signals::Vector{<:Real}, n_pad::Integer, pad::Union{Real, Nothing} = nothing)
    append_back = Vector{Float64}(undef, n_pad)

    if pad === nothing
        fill!(append_back, signals[end])
        return [signals; append_back]
    end

    fill!(append_back, pad)
    return [signals; append_back]
end

function ma_cumsum(signals::Vector{<:Real}, windows::Integer)
    back = convert(Int64, floor(windows / 2))

    signals = pad_end(signals, back)

    signals = cumsum(signals)
    signals[windows-1:end] = (@view signals[windows-1:end]) - (@view signals[1:(end-windows+2)])

    return (@view signals[back:end-1]) ./ windows
end

function butter_lowpass(signals::Vector{<:Real}, low::Real, fs::Real ,order::Integer = 2)
    pass = Lowpass(low; fs)
    coef = digitalfilter(pass, Butterworth(order))
    signals = filtfilt(coef, signals)
    return signals
end

function butter_bandpass(signals::Vector{<:Real}, low::Real, high::Real, fs::Real ,order::Integer = 2)
    pass = Bandpass(low, high; fs)
    coef = digitalfilter(pass, Butterworth(order))
    signals = filtfilt(coef, signals)
    return signals
end


function union_array_of_section(sections_arrays::Vector{<:Vector{<:Tuple{T,T}}}) where T

    all_sections = Vector{Tuple{T,T}}(undef, 0)
    for sections_array in sections_arrays
        append!(all_sections, sections_array)
    end
    sort!(all_sections)
    new_sections = Vector{Vector{T}}(undef, 0)

    for (start_idx, end_idx) in all_sections
        if !isempty(new_sections) && new_sections[end][2] >= start_idx - 1
            new_sections[end][2] = max(new_sections[end][2], end_idx)
        else

            push!(new_sections, [start_idx, end_idx])

        end
    end

    return new_sections
end

function get_segment_index(index_array::Vector{T}, window::Real) where T<: Integer

    segment_index = Vector{Int64}(undef, 0)

    segment_count = convert(Int64, ceil(index_array[1]/window))
    segment_count == 0 ? segment_count = 1 : segment_count

    push!(segment_index, 1)

    for (i, idx) in enumerate(index_array)
        if idx - index_array[segment_index[end]] >= window
            push!(segment_index, i)
        end
    end

    return segment_index

end

function remove_outlier(array::Vector{<:Real}, upper_percentile::Real, lower_percentile::Real)

    index = get_outlier_index(array, upper_percentile, lower_percentile)

    return array[index]
end

function get_outlier_index(array::Vector{<:Real}, upper_percentile::Real, lower_percentile::Real)
    upper, lower = percentile(array, upper_percentile, lower_percentile)
    return findall(x -> (x>lower && x < upper), array)
end

function percentile(array::Vector{<:Real}, upper_percentile::Real, lower_percentile::Real)
    upper_index = convert(Int64, floor(length(array)*upper_percentile/100))
    lower_index = convert(Int64, ceil(length(array)*lower_percentile/100))
    sort_array = sort(array)
    
    upper_index = max(upper_index, 1)
    lower_index = min(lower_index, length(array))
    return sort_array[upper_index], sort_array[lower_index]
end