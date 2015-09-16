"""An image, taken though a particular filter band"""
immutable Image
    H::Int64  # The image height.
    W::Int64  # The image width.
    pixels::Matrix{Float64}  # An HxW matrix of pixel flux measurements in DN
    nmgy_per_dn::Float64  # expected digital number (DN) per nanomaggie
    sky_noise_mean::Float64  # background noise rate in nanomaggies
    read_noise_var::Float64  # Gaussian noise (mean 0) from the CCD, in DN
    gain::Float64  # expected number of photoelectrons per DN (?)
    psf::Matrix{Float64}  # the point spread function
    t::Float64  # the time the image was taken (Or is this an integer?)
end


"""The parameters that characterize a single asteroid"""
immutable AsteroidParams
    r::Float64  # brightness in nanomaggies
    u::Vector{Float64}  # position at time 0
    v::Vector{Float64}  # velocity (constant over time)
end


immutable Prior
    r::LogNormal
    v::MvNormal
end


