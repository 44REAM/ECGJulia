
export msec_to_sample
export sec_to_sample
export sample_to_sec
export sample_to_msec
export sample_to_hr

msec_to_sample(msec:: Real, fs::Real, to_int::Bool = true) = 
    to_int ? convert(Integer, msec*fs/1000) : msec*fs/1000
sec_to_sample(sec:: Real, fs::Real, to_int::Bool = true) = 
    to_int ? convert(Integer, sec*fs) : sec*fs
sample_to_sec(sample:: Real, fs::Real) = sample/fs
sample_to_msec(sample:: Real, fs::Real) = sample/fs * 1000
hr_to_sample(hr::Real, fs::Real, to_int::Bool = true) = sec_to_sample(60/hr, fs, to_int)
sample_to_hr(sample::Real, fs::Real) = 60/sample_to_sec(sample, fs)

msec_to_sample(500.0, 500)