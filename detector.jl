
include("signalobject.jl")
include("peakobject.jl")
include("config/ecgconfig.jl")
include("utils/utils.jl")

using Wavelets
using PyCall
using Statistics

using .signalobject
using .peakobject
import .ecgconfig as cfg
using .utils
pywt = pyimport("pywt")


function preprocess(signals::Vector{T}, loss_threshold = 6000,
    noise_threshold = 100) where T<:Real


    edge = sec_to_sample(2, cfg.FS, true)
    question_segments = utils.mask_segment(signals, loss_threshold, edge)
    question_segments = utils.connect_segment(question_segments)
    utils.replace_zeros!(signals, question_segments)

    signals::Vector{T} = signals - utils.ma(signals, msec_to_sample(600, cfg.FS, true))

    question_segments2 = utils.mask_segment(signals, noise_threshold, edge)
    question_segments2 = utils.connect_segment(question_segments2)
    utils.replace_zeros!(signals, question_segments2)

    question_segments = utils.union_array_of_section([question_segments, question_segments2])

    signals = utils.butter_lowpass(signals, 25, cfg.FS, 2)
    return signalobject.SignalObject(signals), question_segments

end

function test(signals_object::SignalObject)
    print("OK")
end

function detect_peak(signals_object::SignalObject{T},
    detect_algo::String = "above_mean_fix", enhance_algo::String = "abs_dwt_enhancer",
    return_enhanced::Bool = false) where T

    enhance_signals::Vector{T} = enhancer(signals_object.signals, enhance_algo)

    index_array = peak_detector(enhance_signals, detect_algo)

    if return_enhanced
        return index_array, enhance_signals
    end

    return index_array
end

function detect_r_peak(signals_object::SignalObject{<:Real}, index_array::Vector{<:Integer}) 

    peakarray = PeakArray([Event(0)])
    peakarray.fs = cfg.FS

    r_spacing = msec_to_sample(cfg.DETECTOR_RPEAK_SPACING_MSEC, peakarray.fs)
    r_search_range = msec_to_sample(50, peakarray.fs)

    n_signals = length(signals_object.signal_mean_absolute)

    for index in index_array
        if index - peakarray[end].idx > r_spacing
            start_idx = max(index - r_search_range, 0)
            end_idx = min(index + r_search_range, n_signals)

            _, max_idx = findmax((@view signals_object.signal_mean_absolute[start_idx:end_idx] ))
            push!(peakarray, QRS(max_idx+ start_idx -1))
        end
    end

    return peakarray
end

# function detect_peak_above_mean_fix2(signals::Vector{T}) where T <:Real
#     n_signals = length(signals)

#     smooth::Vector{T} = utils.ma(signals, msec_to_sample(cfg.DETECTOR_AMF_SMOOTH_MSEC, cfg.FS))

#     stride = 2

#     index_array = Vector{Int}(undef, 0)

#     width_threshold = msec_to_sample( cfg.DETECTOR_AMF_WIDTH_THRESHOLD_MSEC,cfg.FS, false)/stride

#     smooth = abs.(smooth) .+  cfg.DETECTOR_AMF_THRESHOLD

#     starting::Bool = false
#     start::Int = 0

#     for i in 1:stride:n_signals

#         if signals[i] > smooth[i]
#             if !starting
#                 start = i
#                 starting = true
#             end
#             continue
#         elseif i-start > width_threshold && starting
#             window = @view signals[start:i]
#             max_val, max_idx = findmax(window)


#             if max_val - smooth[(start + max_idx -1) ]<0.2*max_val
#                 starting = false
#                 continue
#             end
#             push!(index_array, start + max_idx -1)
#         end
#         starting = false

#     end

#     return index_array
# end

# function detect_peak_above_mean_fix(signals::Vector{T}) where T <:Real
#     n_signals = length(signals)

#     smooth::Vector{T} = utils.ma(signals, msec_to_sample(cfg.DETECTOR_AMF_SMOOTH_MSEC, cfg.FS))

#     stride = 2

#     window = Vector{T}(undef, 0)
#     index_array = Vector{Integer}(undef, 0)

#     width_threshold = msec_to_sample( cfg.DETECTOR_AMF_WIDTH_THRESHOLD_MSEC,cfg.FS, false)/stride

#     smooth = abs.(smooth) .+  cfg.DETECTOR_AMF_THRESHOLD

#     for i in 1:stride:n_signals

#         if signals[i] > smooth[i]
#             push!(window, signals[i])
#             continue
#         elseif length(window) > width_threshold
#             max_val, max_idx = findmax(window)

#             if max_val - smooth[i - length(window) + max_idx]<0.2*max_val
#                 window = Vector{T}(undef, 0)
#                 continue
#             end
#             push!(index_array, i - length(window) + max_idx)
#         end
#         window = Vector{T}(undef, 0)

#     end

#     return index_array
# end


function detect_peak_above_mean_fix(sig::SignalObject{T}) where T <:Real
    signals = sig.signal_enhance
    n_signals = length(signals)

    smooth::Vector{T} = utils.ma(signals, msec_to_sample(cfg.DETECTOR_AMF_SMOOTH_MSEC, cfg.FS))

    stride = 2

    index_array = Vector{Int}(undef, 0)

    width_threshold = msec_to_sample( cfg.DETECTOR_AMF_WIDTH_THRESHOLD_MSEC,cfg.FS, false)/stride

    smooth = abs.(smooth) .+  cfg.DETECTOR_AMF_THRESHOLD

    starting::Bool = false
    start::Int = 0

    for i in 1:stride:n_signals

        if signals[i] > smooth[i]
            if !starting
                start = i
                starting = true
            end
            continue
        elseif i-start > width_threshold && starting
            window = @view signals[start:i]
            max_val, max_idx = findmax(window)
            max_idx *= stride

            if max_val - smooth[(start + max_idx -1) ]<0.2*max_val
                starting = false
                continue
            end
            push!(index_array, start + max_idx -1)
        end
        starting = false

    end

    return index_array
end

function detect_t_second_derivative!(signals::SignalObject{T}, peakarray::PeakArray{PeakT}) where {T <:Real, PeakT <: Integer}

    start_from_r_peak = msec_to_sample(cfg.DETECTOR_T_FROM_R_MSEC ,peakarray.fs)
    n_signal = length(signals)
    interval_ratio = cfg.DETECTOR_T_INTERVAL_RATIO
    signal_mean_absolute::Vector{T} = signals.signal_mean_absolute
    signal_second_derivative::Vector{T} = signals.signal_second_derivative

    for peak in peakarray.array
        if (!(peak isa QRS) || isnothing(peak.next_interval) || !peak.normal) continue end
        end_from_r_peak = convert(PeakT, floor(peak.next_interval*interval_ratio))

        if end_from_r_peak <= start_from_r_peak continue end

        start_idx = min(peak.idx+start_from_r_peak, n_signal)
        end_idx = min(peak.idx+end_from_r_peak, n_signal)

        landmark = @view signal_mean_absolute[start_idx:end_idx]

        if length(landmark) > 2
            t_peak_idx = argmax(landmark) + peak.idx + start_from_r_peak - 2
            peak.t = t_peak_idx
        else
            continue
        end
        t_peak = peak.t
        end_idx = min(convert( PeakT ,floor(0.3*peak.next_interval))+t_peak, end_from_r_peak+peak.idx)

        landmark = @view signal_second_derivative[t_peak: end_idx]
        if length(landmark) >2
            peak.t_end = argmax(landmark) + t_peak -1
        end
    end
end

function detect_qrs_onset_offset_second_derivative!(signals::SignalObject{T}, peakarray::PeakArray{PeakT}) where {T <:Real, PeakT <: Integer}
    search_near_r_long =  utils.msec_to_sample(100, peakarray.fs)
    search_near_r_short =  utils.msec_to_sample(20, peakarray.fs)

    second_dif_signals::Vector{T} =  signals.signal_second_derivative
    n_signals = length(second_dif_signals)

    for peak in peakarray.array
        if (!(peak isa QRS) || isnothing(peak.interval) || isnothing(peak.next_interval)) continue end

        start_idx = max(peak.idx-search_near_r_long, 0)
        end_idx = max(peak.idx-search_near_r_short, 0)

        landmask = @view second_dif_signals[start_idx: end_idx]
        if length(landmask) != 0
            onset_idx = argmax(landmask) + start_idx -1
            peak.onset = onset_idx
        end

        start_idx = min(peak.idx + search_near_r_short, n_signals )
        end_idx = min(peak.idx + search_near_r_long, n_signals)

        landmask = @view second_dif_signals[start_idx: end_idx]
        if length(landmask) != 0
            offset_idx = argmax(landmask) + start_idx-1
            peak.offset = offset_idx
        end

    end
end

# function detect_peak_above_mean_fix_mod(sig::SignalObject{T}) where T <:Real
#     signals = sig.signal_enhance
#     n_signals = length(signals)

#     smooth::Vector{T} = utils.ma(signals, msec_to_sample(cfg.DETECTOR_AMF_SMOOTH_MSEC, cfg.FS))

#     stride = 2

#     window = Vector{T}(undef, 0)
#     index_array = Vector{Integer}(undef, 0)

#     width_threshold = msec_to_sample( cfg.DETECTOR_AMF_WIDTH_THRESHOLD_MSEC,cfg.FS, false)/stride

#     smooth = abs.(smooth) .+  cfg.DETECTOR_AMF_THRESHOLD

#     for i in 1:stride:n_signals

#         if signals[i] > smooth[i]
#             push!(window, signals[i])
#             continue
#         elseif length(window) > width_threshold
#             max_val, max_idx = findmax(window)

#             if max_val - smooth[i - length(window) + max_idx]<0.2*max_val
#                 window = Vector{T}(undef, 0)
#                 continue
#             end
#             push!(index_array, i - length(window) + max_idx)
#         end
#         window = Vector{T}(undef, 0)

#     end

#     return index_array
# end

function abs_swt_enhancer!(signals::SignalObject{<:Real})

    level = 3
    check = 2^level
    left = length(signals.signals)%check
    if left !=0
        signals_en = utils.pad_end(signals.signals, check - left, 0)
    else
        signals_en = signals.signals
    end

    signals_en = abs.(pywt.swt(signals_en, "db3", 3)[1][2])
    signals_en = utils.butter_bandpass(signals_en, 0.01, 8.0, cfg.FS).*20
    signals_en = signals_en - utils.ma(signals_en, msec_to_sample(600, cfg.FS))

    signals.signal_enhance = signals_en

end

function abs_dwt_enhancer!(signals::SignalObject{<:Real})
    level = 3
    check = 2^level
    left = length(signals.signals)%check
    if left !=0
        signals_en = utils.pad_end(signals.signals, check - left, 0)
    else
        signals_en = signals.signals
    end


    signals_en = abs.( dwt(signals_en, wavelet(WT.db3), level))
    signals_en = utils.butter_bandpass(signals_en, 0.01, 8.0, cfg.FS).*20
    signals_en = signals_en - utils.ma(signals_en, msec_to_sample(600, cfg.FS))
    signals.signal_enhance = signals_en
end

function tachy_brady!(peakarray::peakobject.PeakArray,segment_index::Vector{<:Integer} ,algo::String = "time")
    if algo == "time"
        tachy_brady_time!(peakarray, segment_index)
    #TODO
    # elseif algo == "peak"
    #     tachy_brady_pause_peak(peakarray)
    end
end

function tachy_brady_time!(peakarray::peakobject.PeakArray, segment_index::Vector{<:Integer})

    start_idx = 1
    threshold_tachy = cfg.HR_TACHY
    threshold_brady = cfg.HR_BRADY
    fs = peakarray.fs

    for end_idx in segment_index
        peakarray_landmark = peakarray[start_idx: end_idx]
        interval, _ = peakobject.get_interval(peakarray_landmark,false, false, true)
        if length(interval) < 3 continue end

        hr = sample_to_hr(mean(interval), fs)

        for peak in peakarray_landmark.array
            if peak isa QRS peakobject.assigned_hr!(peak, hr, threshold_tachy, threshold_brady) end
        end
        start_idx = end_idx

    end
end

function tachy_brady_peak!(peakarray::peakobject.PeakArray)
end

function ectopic!(peakarray::peakobject.PeakArray, segment_index::Vector{<:Integer}, algo::String = "time")
    if algo == "time" 
        ectopic_time!(peakarray, segment_index )
    else 
        ectopic_time!(peakarray, segment_index ) 
    end

end

function ectopic_time!(peakarray::peakobject.PeakArray,segment_index::Vector{<:Integer} )
    start_idx = 1
    ectopic_ratio = cfg.ECTOPIC_RATIO
    qrs_width_sample = msec_to_sample(cfg.ECTOPIC_WIDTH_MSEC, peakarray.fs)

    for end_idx in segment_index

        peakarray_landmark = peakarray[start_idx: end_idx]
        interval_array, _ = peakobject.get_interval(peakarray_landmark,false, false, false)
        if length(interval_array) == 0 continue end
        m_rr = median(interval_array)

        for peak in peakarray_landmark.array
            if !(peak isa QRS) continue end
            if isnothing(peak.interval) continue end

            if peak.interval < m_rr * ectopic_ratio
                pvc_pac!(peak, qrs_width_sample)
            end
        end
        start_idx = end_idx
    end
end

function pvc_pac!(peak::peakobject.QRS, qrs_width_sample)

    if isnothing(peak.onset) || isnothing(peak.offset)
        peakobject.add_diagnosis!(peak, peakobject.PAC)
        return
    end

    qrs_width = peak.offset - peak.onset
    if qrs_width > qrs_width_sample peakobject.add_diagnosis!(peak,peakobject.PVC) 
    else peakobject.add_diagnosis!(peak, peakobject.PAC) end 
end
