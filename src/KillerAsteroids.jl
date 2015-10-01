module KillerAsteroids

using FITSIO
using WCSLIB
using Distributions

include("model_types.jl")
include("model_probability.jl")
include("wise_psf.jl")
include("wise_utils.jl")

export Image, AsteroidParams, Prior
export compute_log_probability, compute_log_likelihood
export load_wise_psf
export wise_params

end # module
