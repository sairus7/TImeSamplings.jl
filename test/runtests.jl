import TimeSamplings as ts
using Test
using Dates
using Unitful
using Unitful: s, ms
using IntervalSets

@testset "TimeSamplings.jl" begin

# ==========
# astronomic time
tstart = parse(DateTime, "2021-02-13T23:34:42")
s1 = ts.TimeSampler(tstart)

delta = Millisecond(1000)
tend = tstart + delta
# scalar
@test s1[delta] == s1[1s] == tend
@test s1(tend) == delta
# vector
@test s1[[0s, delta]] == [tstart, tend]
@test s1([tstart, tend]) == [0s, delta]
# interval
@test s1[0s..delta] == tstart..tend
@test s1(tstart..tend) == 0s..delta
# step range
@test s1[0ms:1ms:1000ms] == tstart:Millisecond(1):tend
@test s1(tstart:Millisecond(1):tend) == 0ms:1ms:1000ms

# ==========
# discrete-time samples
fs = 250
s2 = ts.DiscreteSampler(fs)
# scalar -> scalar
@test s2[1] == 0ms
@test s2[2] == 4ms
@test s2(0ms) == s2(3ms) == 1 # all values between two discrete timestamps are truncated to nearest previous value
@test s2(4ms) == 2
# vector -> vector
@test s2[[1, 2]] == [0ms, 4ms]
@test s2([0ms, 4ms]) == [1, 2]
# range
@test s2[1:3] == 0ms..8ms # unit range -> interval
@test s2[1:1:3] == 0ms:4ms:8ms # step range -> vector
@test s2(0ms:4ms:1000ms) == 1:1:251
# interval
# s2[1..3] == 0ms..8ms # interval - not available
@test s2(0ms..1000ms) == 1:251 # interval -> range

# ==========
# downsampled 10-second trend
fdcm = 1/10
factor = fs / fdcm

s3_l = ts.DownSampler(factor, ts.Left)
s3_c = ts.DownSampler(factor, ts.Center)
s3_r = ts.DownSampler(factor, ts.Right)

# scalar
@test s3_l[1] == 1 # left (first) point of decimation window
@test s3_c[1] == 1 + factor/2 # center of decimation window
@test s3_r[1] == factor # right (last) point of decimation window

@test s3_l[2] == factor + 1
@test s3_c[2] == 1 + factor/2 + factor
@test s3_r[2] == factor * 2

@test s3_l(1) == s3_r(1) == s3_c(1) == 1
@test s3_l(factor) == s3_r(factor) == s3_c(factor) == 1
@test s3_l(factor + 1) == s3_r(factor + 1) == s3_c(factor + 1) == 2
@test s3_l(factor * 2) == s3_r(factor * 2) == s3_c(factor * 2) == 2

# vector
@test s3_l[[1,2,3]] == [1, 1+factor, 1+2factor]
@test s3_l([1, factor/2, 1+factor, 1+2factor]) == [1,1,2,3]
# range
@test s3_l[1:3] == 1:1+2factor
@test s3_l(1:factor:1+2factor) == [1,2,3]
# interval - not available
# s3_l[1..3] == 1..1+2factor
# s3_l(1..1+2factor) == 1..3

# ==========
# shifted by constant
shift = 15
s4 = ts.ShiftSampler(shift)

# scalar
@test s4[1] == 1 + shift
@test s4(1 + shift) == 1
# vector
@test s4[[1, 2]] == [1+shift, 2+shift]
@test s4([1+shift, 2+shift]) == [1, 2]
# range
@test s4[1:2] == 1+shift:2+shift
@test s4(1+shift:2+shift) == 1:2
# interval - not available
# s4[1..2] == 1+shift..2+shift
# s4(1+shift..2+shift) == 1..2

# ==========
# events
events = [123, 321, 555]
se1 = ts.EventSampler(events)
# scalar
@test se1[1] == 123
@test se1(123) == 1:1
@test se1(200) == 2:1 # first after : last before
# vector
@test se1[[1,2]] == [123, 321]
@test se1([100, 123]) == [1:0, 1:1]
# range
@test se1[1:2] == [123, 321] # range --> vector
@test se1(100:110) == 1:0
@test se1(100:350) == 1:2
# interval - not available
# se1[1..2] == 123..321
# se1(100..350) == 1..2

# ==========
# chained
tstart = parse(DateTime, "2021-02-13T23:34:42")
sampt = ts.TimeSampler(tstart)
# signal-1: 1000 Hz
samp1 = ts.DiscreteSampler(1000)
# signal-2: 50 Hz (decimated signal-1 with factor=20)
samp2 = ts.DownSampler(20)
# events positions based on signal-2
samp3 = ts.EventSampler([1, 51, 501]) # 0, 1, 10 seconds

# get events positions on signal-1
to_baseindex = ts.SamplingChain(samp2, samp3)
@test to_baseindex[1:3] == [1, 1001, 10001]
# get event indices within base index range
@test to_baseindex(1:5000) == 1:2

# get events time periods
to_period = ts.SamplingChain(samp1, samp2, samp3)
@test to_period[1:3] == [0s, 1s, 10s]
# get event indices within timeperiod interval
@test to_period(0s..5s) == 1:2

# get events astronomic times
to_datetime = ts.SamplingChain(sampt, samp1, samp2, samp3)
@test to_datetime[1:3] == tstart .+ [Second(0), Second(1), Second(10)]
# get event indices within datetime interval
@test to_datetime(tstart..tstart+Second(5)) == 1:2


# event segments - not is not quite accurate,
# because event ends should be aligned by the right border
samp3s = ts.EventSampler([1:5, 51:55, 501:505]) # 0..0.1, 1..1.1, 10..10.1 seconds

# get events positions on signal-1
to_baseindex = ts.SamplingChain(samp2, samp3s)
@test to_baseindex[1:3] == [1:81, 1001:1081, 10001:10081]

# get event indices within base index range
@test to_baseindex(1:5000) == 1:2
@test to_baseindex(5:1001) == 1:2

# get events time periods
to_period = ts.SamplingChain(samp1, samp2, samp3s)
@test to_period[1:3] == [0ms..80ms, 1000ms..1080ms, 10000ms..10080ms]
# get event indices within timeperiod interval
@test to_period(0s..5s) == 1:2

sampt[[0ms..80ms, 1000ms..1080ms, 10000ms..10080ms]]

# get events astronomic times
to_datetime = ts.SamplingChain(sampt, samp1, samp2, samp3s)
to_datetime[1:3]
# get event indices within datetime interval
@test to_datetime(tstart..tstart+Second(5)) == 1:2

# TODO: add date intervals # == tstart .+ [Second(0), Second(1), Second(10)]

# =============
# test translation

tstart1 = parse(DateTime, "2021-02-13T23:34:42")
tstart2 = tstart1 + Hour(1)

ts1 = ts.Sampler(tstart1)
ts2 = ts.Sampler(tstart2)

# period 2 -> time -> period 1
period2 = Minute(1)
t = ts2[period2]
period1 = ts1(t)
@test period1 == Minute(Hour(1)) + period2

periods1 = ts.translate_periods([period2], tstart2, tstart1)
@test periods1[1] == Minute(Hour(1)) + period2

# samples 2 -> period 2 -> time -> period 1 -> samples 1
fs1 = fs2 = 1000.
samp1 = ts.Sampler(tstart1, fs1)
samp2 = ts.Sampler(tstart2, fs2)

point2 = 1
t = samp2[point2]
point1 = samp1(t)
@test point1 == point2 + 3600 * fs2

points1 = ts.translate_index([point2], tstart2, fs2, tstart1, fs1)
@test points1[1] == point2 + 3600 * fs2

end
