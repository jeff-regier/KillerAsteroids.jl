"""An image, taken though a particular filter band"""
type Image
    H::Int64  # The image height.
    W::Int64  # The image width.
    pixels::Matrix{Float64}  # An HxW matrix of pixel intensities.
    epsilon::Float64  # The background noise in nanomaggies
    iota::Float64  # expect number of photons per nanomaggie
    psf::Matrix{Float64}  # the point spread function
    t::Float64  # the time the image was taken (Or is this an integer?)
end

"""The parameters that characterize a single asteroid"""
type AsteroidParams
    r::Float64  # brightness in nanomaggies
    u::Vector{Float64}  # position at time 0
    v::Vector{Float64}  # velocity (constant over time)
end


type Prior
    r::LogNormal
    v::MvNormal
end


function sample_prior()
    r_prior = LogNormal(1500., sqrt(500)) # mode is 1000
    v_prior = MvNormal([3., 5], [4. 1; 3 2]) # pixels / second

    Prior(r_prior, v_prior)
end

