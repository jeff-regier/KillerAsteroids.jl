module KillerAsteroids

using FITSIO
using WCSLIB
using Distributions

include("model_types.jl")
include("model_probability.jl")
include("wise_psf.jl")

export Image, AsteroidParams, Prior
export compute_log_probability, compute_log_likelihood
export load_wise_psf

end # module
