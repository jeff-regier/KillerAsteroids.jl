module KillerAsteroids

using FITSIO
using WCSLIB
using Distributions
using Interpolations

include("model_types.jl")
include("model_probability.jl")
include("wise_psf.jl")
include("wise_utils.jl")
include("mask_dilate.jl")

export Image, AsteroidParams, Prior
export compute_log_probability, compute_log_likelihood
export load_wise_psf, extrapolate_position
export wise_band_to_params
export l1b_image_name
export l1b_bool_mask
export B # number of filter band
export mask_dilate

end # module
