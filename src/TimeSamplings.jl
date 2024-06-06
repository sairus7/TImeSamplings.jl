module TimeSamplings

using Dates
using Unitful
using IntervalSets

"""
AbstractSampler is interface for forward and backward translation between base and
"self" time samples:
- backward translation - call operator ()
- forward translation - getindex operator []
"""

abstract type AbstractSampler end

"""
Convert base samples to "self" samples: `self = obj(ibase)`
"""
function (obj::AbstractSampler)(base_sample)
    error("$(typeof(obj)) has no implemententation for $(typeof(base_sample))")
end
function (obj::AbstractSampler)(base_sample::AbstractVector{T}) where T # vector -> vector
    broadcast(base_sample) do s
        obj(s)
    end
end
"""
Convert "self" samples to base samples: `ibase = obj[iself]`
"""
function Base.getindex(obj::AbstractSampler, sample)
    error("$(typeof(obj)) has no implemententation for $(typeof(sample))")
end
function Base.getindex(obj::AbstractSampler, sample::AbstractVector{T}) where T # vector -> vector
    broadcast(sample) do s
        getindex(obj, s)
    end
end

const TimePeriod = Union{Period, Unitful.Time}

# internal time conversion functions
@inline offset2ms(fs::Float64, offset::Real) = offset * 1000 / fs
@inline ms2offset(fs::Float64, ms::Real) = ms * fs / 1000

# from 1-based index - to floating-point milliseconds
@inline index2ms(fs::Float64, index::Real) = offset2ms(fs, index - 1)
@inline ms2index(fs::Float64, ms::Real) = ms2offset(fs, ms) + 1

# time Period to integer milliseconds (may error for nanoseconds)
@inline time2ms(time::Union{Unitful.Time, Period}) = Dates.value(Millisecond(time))
@inline ms2time(msec) = Millisecond(msec) # integer milliseconds to time Period
# @inline ms2time(msec) = DateTime(Dates.UTM(msec)) - transforms to absolute DateTime

"""
Conversion from(to) relative time periods to(from) astronomic datetime
"""
struct TimeSampler <: AbstractSampler
    timestart::DateTime
end
(obj::TimeSampler)(time::DateTime) = time - obj.timestart # scalar -> scalar
(obj::TimeSampler)(time::ClosedInterval{<:DateTime}) = begin # interval -> interval
    t1, t2 = endpoints(time)
    ClosedInterval(obj(t1), obj(t2))
end
(obj::TimeSampler)(time::AbstractRange{<:DateTime}) = begin # range -> range
    t1, tstep, t2 = first(time), step(time), last(time)
    obj(t1):tstep:obj(t2)
end
Base.getindex(obj::TimeSampler, period::TimePeriod) = obj.timestart + ms2time(period) # scalar -> scalar
Base.getindex(obj::TimeSampler, period::ClosedInterval{<:TimePeriod}) = begin # interval -> interval
    p1, p2 = endpoints(period)
    ClosedInterval(obj[p1], obj[p2])
end
Base.getindex(obj::TimeSampler, period::AbstractRange{<:TimePeriod}) = begin # range -> range
    p1, tstep, p2 = first(period), step(period), last(period)
    obj[p1]:ms2time(tstep):obj[p2]
end

"""
Conversion from(to) discrete samples to(from) time periods
"""
struct DiscreteSampler <: AbstractSampler
    fs::Float64 # sampling rate
end
(obj::DiscreteSampler)(period::TimePeriod) = trunc(Int, ms2index(obj.fs, time2ms(period))) # scalar -> scalar
(obj::DiscreteSampler)(period::ClosedInterval{<:TimePeriod}) = begin # interval -> unit range
    p1, p2 = endpoints(period)
    range(obj(p1), obj(p2))
end

Base.getindex(obj::DiscreteSampler, index::Real) = ms2time(floor(Int, index2ms(obj.fs, index))) # scalar -> scalar
Base.getindex(obj::DiscreteSampler, index::UnitRange{<:Real}) = begin # unit range -> interval
    i1, i2 = first(index), last(index)
    ClosedInterval(obj[i1], obj[i2])
end

"""
Alignment modes for downsampling
"""
@enum Alignment Left Center Right

"""
Conversion from(to) decimated discrete samples to(from) dense discrete samples
by given decimation factor and alignment
"""
struct DownSampler <: AbstractSampler
    factor::Float64
    offset::Float64
    function DownSampler(factor, alignment::Alignment = Left)
        offset = if alignment == Left
            0.
        elseif alignment == Center
            factor / 2
        else # if alignment == Right
            factor - 1
        end
        new(factor, offset)
    end
end
(obj::DownSampler)(baseindex::Real) = trunc(Int, (baseindex-1 #=-obj.offset=#)/obj.factor + 1)
(obj::DownSampler)(baseindex::UnitRange{<:Real}) = begin # unit range -> unit range
    i1, i2 = first(baseindex), last(baseindex)
    range(obj(i1), obj(i2))
end
Base.getindex(obj::DownSampler, index::Real) = trunc(Int, (index-1)*obj.factor + 1 + obj.offset)
Base.getindex(obj::DownSampler, index::AbstractRange{<:Real}) = begin # unit range -> unit range
    i1, i2 = first(index), last(index)
    range(obj[i1], obj[i2])
end
"""
Shift samples by constant value
"""
struct ShiftSampler <: AbstractSampler
    shift::Float64
end
(obj::ShiftSampler)(baseindex::Real) = trunc(Int, baseindex - obj.shift) # rounding conflicts with range!
(obj::ShiftSampler)(baseindex::AbstractRange{<:Real}) = begin # unit range -> unit range
    i1, i2 = first(baseindex), last(baseindex)
    range(obj(i1), obj(i2))
end
Base.getindex(obj::ShiftSampler, index::Real) = trunc(Int, index + obj.shift)
Base.getindex(obj::ShiftSampler, index::AbstractRange{<:Real}) = begin # unit range -> unit range
    i1, i2 = first(index), last(index)
    range(obj[i1], obj[i2])
end
"""
Conversion from(to) event index to(from) discrete samples or times
by given sorted events positions / timestamps
"""
struct EventSampler{T} <: AbstractSampler
    indices::Vector{T} # sorted values
end
(obj::EventSampler)(baseindex::Union{Real, UnitRange{<:Real}}) = searchsorted(obj.indices, baseindex, lt=(x,y)->isless(last(x), first(y))) # scalars / range -> range
Base.getindex(obj::EventSampler, index::Union{Real, AbstractVector{<:Real}}) = obj.indices[index]

"""
Chained conversion through multiple samplings
"""
struct SamplingChain{S1 <: AbstractSampler, S2 <: AbstractSampler} <: AbstractSampler
    parent::S1
    self::S2
    function SamplingChain(parent, self, children...) # 2 or more samplings chained
        for s in (self, children...)
            parent = new{typeof(parent), typeof(s)}(parent, s)
        end
        return parent
    end
end
(obj::SamplingChain)(base_sample) = obj.self(obj.parent(base_sample))
(obj::SamplingChain)(base_sample::AbstractVector{T}) where T = obj.self(obj.parent(base_sample)) # solve ambuiguity
Base.getindex(obj::SamplingChain, sample) = obj.parent[obj.self[sample]]
Base.getindex(obj::SamplingChain, sample::AbstractVector{T}) where T = obj.parent[obj.self[sample]] # solve ambuiguity

"""
Make sampler from `fs` and/or `timestart`
returns concrete object of AbstractSampler
"""
Sampler(timestart::DateTime) = TimeSampler(timestart)
Sampler(fs::Float64) = DiscreteSampler(fs)
Sampler(indices::Vector{T}) where {T} = EventSampler(indices)
function Sampler(timestart::DateTime, fs::Float64)
    SamplingChain(
        TimeSampler(timestart),
        DiscreteSampler(fs)
    )
end
function Sampler(fs::Float64, indices::Vector)
    SamplingChain(
        DiscreteSampler(fs),
        EventSampler(indices)
    )
end
function Sampler(timestart::DateTime, fs::Float64, indices::Vector)
    SamplingChain(
        TimeSampler(timestart),
        DiscreteSampler(fs),
        EventSampler(indices)
    )
end


"""
Translate from `index_src` into `index_dst`
"""
function translate_index(
    index_src::AbstractVector{<:Integer},
    timestart_src::DateTime,
    fs_src::Real,
    timestart_dst::DateTime,
    fs_dst::Real)

    samp_src = Sampler(timestart_src, fs_src)
    samp_dst = Sampler(timestart_dst, fs_dst)

    broadcast(index_src) do i
        samp_dst(samp_src[i])
    end
end

"""
Translate from `index_src` into `index_dst`
"""
function translate_index(
    index_src::AbstractVector{<:Integer},
    fs_src::Real, fs_dst::Real)

    samp_src = Sampler(fs_src)
    samp_dst = Sampler(fs_dst)

    broadcast(index_src) do i
        samp_dst(samp_src[i])
    end
end

"""
Translate from `periods_src` into `periods_dst`
"""
function translate_periods(
    periods_src::AbstractVector{<:TimePeriod},
    timestart_src::DateTime,
    timestart_dst::DateTime)

    samp_src = Sampler(timestart_src)
    samp_dst = Sampler(timestart_dst)

    broadcast(periods_src) do p
        samp_dst(samp_src[p])
    end
end

end
