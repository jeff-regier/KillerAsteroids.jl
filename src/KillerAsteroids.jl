module KillerAsteroids

using FITSIO
using WCSLIB
using Distributions

const B = 4
export B # number of filter band

include("model_types.jl")
include("model_probability.jl")
include("wise_psf.jl")

export Image, AsteroidParams, Prior
export compute_log_probability, compute_log_likelihood
export load_wise_psf, extrapolate_position

end # module
