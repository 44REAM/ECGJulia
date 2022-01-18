
module peakobject

using Statistics

include("utils/utils.jl")

using .utils

using Base

abstract type Peak{T} end

export QRS, Event, Question, PeakArray, Peak
export UNKNOWN, NORMAL, BRADY, TACHY, PAC, PVC, AF, VTVF, PAUSE, QUESTION, EVENT

@enum ECGLabel begin
    UNKNOWN = 0
    NORMAL = 1
    BRADY = 2
    TACHY = 3
    PAC = 4
    PVC = 5
    AF = 6
    VTVF = 7
    PAUSE = 8
    QUESTION = 9
    EVENT = 10
end

mutable struct QRS{T} <: Peak{T}
    idx::T
    diag::ECGLabel
    normal::Bool

    t::Union{T, Nothing}
    t_end::Union{T, Nothing}
    onset::Union{T, Nothing}
    offset::Union{T, Nothing}

    interval::Union{T, Nothing}
    next_interval::Union{T, Nothing}

    hr::Union{Float64, Nothing}
    QRS(idx::T) where T = new{T}(idx, UNKNOWN, true)
end

mutable struct Question{T} <: Peak{T}
    idx::T
    Question(idx::T) where T = new{T}(idx)
end

mutable struct Event{T} <: Peak{T}
    idx::T
    Event(idx::T) where T = new{T}(idx)
end


function assigned_hr!(peak::QRS, hr::Real, tachy::Real, brady::Real)
    peak.hr = hr

    if hr > tachy add_diagnosis!(peak, TACHY) 
    elseif hr > brady add_diagnosis!(peak, BRADY) 
    end
end

function assigned_interval!(peak::QRS, interval::Real, pause_threshold::Real)
    peak.interval = interval
    if interval > pause_threshold add_diagnosis!(peak, PAUSE) end
end

function add_diagnosis!(peak::QRS, diag::ECGLabel)
    if diag>peak.diag peak.diag = diag end
    if (diag == PAC) || (diag == PVC) || (diag == VTVF)
        peak.normal = false
    end
end

function shift!(peak::Peak, index::Integer)
    peak.idx += index
end

function resample!(peak::Peak, ratio::Real)
    peak.idx *= ratio
end

mutable struct PeakArray{T}
    array::Vector{Peak{T}}
    start_ts::UInt64
    fs::Union{UInt16, Nothing}
    PeakArray(array::Vector{<:Peak{T}}) where{T} = new{T}(array, 0)
end

function assigned_peak_by_section!(peakarray::PeakArray, section_array::Vector{<:Vector{<:Integer}}, func::Function)
    section_idx = 1
    n_section = length(section_array)

    if n_section == 0 return end

    for (idx, peak) in enumerate(peakarray.array)
        while peak.idx > section_array[section_idx][2]
            section_idx+=1
            if section_idx>= n_section
                return
            end
        end

        if peak.idx >= section_array[section_idx][1]
            peakarray[idx] = func(peak)
        end
    end
end

function assigned_question!(peakarray::PeakArray, question_section::Vector{T}) where T <:Vector{<:Integer}
    f(peak::Peak) = Question(peak.idx)
    assigned_peak_by_section!(peakarray, question_section, f)
    
    for section::T in question_section
        tmp = max(  convert(  Int,floor((section[1] + section[2])/2)), 1)
        push!(peakarray, Event( tmp))
    end
    sort!(peakarray)
end

function assigned_unknown!(peak::QRS)
    if peak.diag == UNKNOWN
        peak.normal = false
        return
    end

    if isnothing(peak.interval) || isnothing(peak.next_interval)
        peak.diag = UNKNOWN
        peak.normal = false
    end

end

function assigned_unknown!(peakarray::PeakArray)
    for peak in peakarray.array
        if !(peak isa QRS) continue end
        
        if peak.diag == UNKNOWN  assigned_unknown!(peak) end
    end
end


function assigned_interval!(peakarray::PeakArray, pause_threshold::Real)
    n_peak = length(peakarray)

    for i in 2:n_peak
        previous_peak = peakarray.array[i-1]
        this_peak = peakarray.array[i]
        if ( previous_peak isa QRS) && (this_peak isa QRS)
            interval = this_peak.idx - previous_peak.idx
            previous_peak.next_interval = interval

            assigned_interval!(this_peak, interval, pause_threshold)
        end
    end
end


function get_index(peakarray::PeakArray)
    return map(peak -> peak.idx, peakarray.array)
end

function get_interval(peakarray::PeakArray{T}, nn::Bool, to_timestamp::Bool = false, rm_outlier::Bool = false,
    upper_percentile::Real = 95.0, lower_percentile::Real = 5.0) where T <: Integer

    function condition(nn::Bool)
        if nn
            f(peak::Peak) = ((peak isa QRS) && !(peak.interval === nothing)) && peak.normal
            return f
        else
            g(peak::Peak) = ((peak isa QRS) && !(peak.interval === nothing))
            return g
        end
    end 

    f::Function = condition(nn)

    interval::Vector{T} = [peak.interval for peak in peakarray.array if f(peak)]
    index::Vector{T} = [peak.idx for peak in peakarray.array if f(peak)]

    if length(interval) == 0
        return interval, index
    end
    if rm_outlier
        outlier_index = utils.get_outlier_index(interval, upper_percentile, lower_percentile)
        interval = interval[outlier_index]
        index = index[outlier_index]
    end

    if to_timestamp
        map!(sample_to_msec, interval)
        map!(sample_to_msec, index)
    end

    return interval, index
end

function get_different_interval(peakarray::PeakArray{T}, rm_outlier::Bool = false,
    upper_percentile::Real=95.0, lower_percentile::Real=5.0) where T <: Integer

    interval::Vector{T} = Vector{T}([])
    different_interval::Vector{T} = Vector{T}([])

    for peak in peakarray.array
        if (!(peak isa QRS) || isnothing(peak.interval))
            if length(interval) >=10
                if rm_outlier
                    interval = utils.remove_outlier(interval, upper_percentile, lower_percentile)
                end
                append!(different_interval, diff(interval))
            end
            interval = Vector{T}([])
            continue
        end
        push!(interval, peak.interval)
    end

    if length(interval) >=10
        if rm_outlier
            interval = utils.remove_outlier(interval, upper_percentile, lower_percentile)
        end
        append!(different_interval, diff(interval))
    end


    if length(different_interval) <= 10 return Vector{T}(undef,0) end
    if rm_outlier
        different_interval = utils.remove_outlier(different_interval, upper_percentile, lower_percentile)
    end

    return different_interval

end


function get_segment_index(peakarray::PeakArray, window::Real)
    utils.get_segment_index(map(peak -> peak.idx, peakarray.array), window)
end

function get_qt(peakarray::PeakArray)
    
end

function get_metrics(peakarray::PeakArray, upper_percentile::Real = 98.0, lower_percentile::Real = 2.0)
    to_timestamp::Bool = false
    rm_outlier::Bool = true

    min_number_of_peak = 3

    sdann = nothing
    asdnn = nothing
    sdnn = nothing
    hr = nothing
    rmssd = nothing

    interval, index = get_interval(peakarray, true, to_timestamp, rm_outlier, upper_percentile, lower_percentile)
    if length(interval) > min_number_of_peak
        duration = sec_to_sample(5, peakarray.fs, false)
        sdann, asdnn = utils.sdann_asdnn(interval, index, duration, min_number_of_peak)
        sdnn = utils.sdnn(interval)

        if !isnothing(sdann)
            sdann = sample_to_msec(sdann, peakarray.fs)
            asdnn = sample_to_msec(asdnn, peakarray.fs)
        end
        sdnn = sample_to_msec(sdnn, peakarray.fs)
    end

    interval, _ = get_interval(peakarray, false, to_timestamp, rm_outlier, upper_percentile, lower_percentile)
    if length(interval) > min_number_of_peak
        hr = sample_to_hr(mean(interval), peakarray.fs)
    end

    interval = get_different_interval(peakarray, true, upper_percentile, lower_percentile)
    if length(interval) > min_number_of_peak
        rmssd = utils.sdnn(interval)
        rmssd = sample_to_msec(rmssd, peakarray.fs)
    end

    results = Dict(
        [
            ("hr", hr),
            ("sdnn", sdnn),
            ("asdnn", asdnn),
            ("sdann", sdann),
            ("rmssd", rmssd)
        ]
    )
    return results
end

Base.:+(peak1::Peak, peak2::Peak) = peak1.idx+ peak2.idx
Base.:-(peak1::Peak, peak2::Peak) = peak1.idx- peak2.idx
Base.:(==)(peak1::Peak, peak2::Peak) = peak1.idx == peak2.idx
Base.:>(peak1::Peak, peak2::Peak) = peak1.idx > peak2.idx
Base.:<(peak1::Peak, peak2::Peak) = peak1.idx < peak2.idx
Base.:>=(peak1::Peak, peak2::Peak) = peak1.idx >= peak2.idx
Base.:<=(peak1::Peak, peak2::Peak) = peak1.idx <= peak2.idx
Base.isless(peak1::Peak, peak2::Peak) = isless(peak1.idx, peak2.idx)

Base.:+(peakarray1::PeakArray, peakarray2::PeakArray) = peakarray1.array + peakarray2.array
Base.:-(peakarray1::PeakArray, peakarray2::PeakArray) = peakarray1.array - peakarray2.array
Base.push!(peakarray::PeakArray, peak::Peak) = push!(peakarray.array, peak)
Base.lastindex(peakarray:: PeakArray) = length(peakarray.array)
Base.getindex(peakarray::PeakArray, index::Integer) = peakarray.array[index]
Base.getindex(peakarray::PeakArray, index::UnitRange{<:Integer}) = PeakArray(peakarray.array[index])
Base.length(peakarray::PeakArray) = length(peakarray.array)
Base.sort!(peakarray::PeakArray) = sort!(peakarray.array)
Base.setindex!(peakarray::PeakArray, peak::Peak, index::Integer) = setindex!(peakarray.array, peak, index)
Base.diff(peakarray::PeakArray) = diff(peakarray.array)

end
