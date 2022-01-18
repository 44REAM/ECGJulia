
module signalobject

include("utils/utils.jl")
include("config/ecgconfig.jl")

import Base

import .ecgconfig as cfg
using .utils

export SignalObject
mutable struct SignalObject{T<:Real}
    signals::Vector{T}
    signal_mean_absolute::Union{Vector{T}, UndefInitializer}
    signal_second_derivative::Union{Vector{T}, UndefInitializer}
    signal_enhance::Union{Vector{T}, UndefInitializer}
    SignalObject(signals::Vector{T}) where T = new{T}(signals)
end

function signal_mean_absolute(signals::Vector{<:Real}, windows_msec::Real = 20)

    windows = msec_to_sample(windows_msec, cfg.FS)
    return  utils.ma(abs.(signals), windows)
end

function signal_second_derivative(signals::Vector{T}) where T
    second_dif_signals = zeros(T, length(signals))
    second_dif_signals[2:end-1] = diff(diff(signals))
    return second_dif_signals
end

function Base.getproperty(sig::SignalObject, s::Symbol)
    if s == :signal_mean_absolute
        if isdefined(sig, s)
            return getfield(sig, s)
        else
            setfield!(sig, s, signal_mean_absolute(sig.signals))
            return getfield(sig, s)
        end
    elseif s == :signal_enhance
        if isdefined(sig, s)
            return getfield(sig, s)
        else
            return sig.signals
        end
    elseif s == :signal_second_derivative
        if isdefined(sig, s)
            return getfield(sig, s)
        else
            setfield!(sig, s, signal_second_derivative(sig.signals))
            return getfield(sig, s)
        end
    else
        return getfield(sig, s)
    end
end

Base.length(sig::SignalObject) = length(sig.signals)

end