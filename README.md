# TimeSamplings


This is a draft package for forward and backward translation between:
- astronomic time
- time periods
- indices of sampling rate (uniform sampling)
- downsampled indices
- indices with offset
- event indices (non-uniform sampling)

Unlike AxisKeys, Samplers can be composed into bi-directional "pipes" 
(see `SamplingChain`) that call each other recursively to make forward and backward 
translations from indices to time scale, from dense to downsampled indices and so on:
- forward translation - getindex operator []
- backward translation - call operator ()

The following behavior is done:
- scalar translates into scalar
- vector translates into vector
- time intervals translate into unit ranges of indexes, and vice versa
- unit index range translate into unit range of other indices (unit step is not changed)
- step index range is translated into a vector (because of rounding)
- indices are always rounded to integers
- precision is fixed down to milliseconds

IntervalSets, Unitful - are used only for timescales and convenient input of time intervals.
Discrete samplers use only index ranges

Currently there are some problems with static alignment modes, open-closed intervals and
rounding ends when event timestamps are specified as index ranges.