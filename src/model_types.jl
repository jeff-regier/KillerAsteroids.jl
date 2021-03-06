"""An image, taken though a particular filter band"""
immutable Image
    H::Int64  # The image height.
    W::Int64  # The image width.
    pixels::Matrix{Float64}  # An HxW matrix of pixel flux measurements in DN
    sky_noise_mean::Float64  # background noise rate in DN
    psf::Matrix{Float64}  # the point spread function
    band_id::Int64  # identifies the filter band for this image (1,2,3 or 4)
    t::Float64  # the time the image was taken (Or is this an integer?)
    wcs::WCSLIB.wcsprm  # the image position
end


"""The parameters that characterize a single asteroid"""
immutable AsteroidParams
    #TODO: specify the asteroid brightness in each band
    # (for now only use images in the same band)
    r::Vector{Float64}  # brightness in each band in nanomaggies
    u::Vector{Float64}  # position (ra / dec) at time 0
    v::Vector{Float64}  # velocity (constant over time)
end


immutable Prior
    log_r::MvNormal # log brightness in each of the 4 bands (nanomaggies)
    v::MvNormal  # velocity in units of (ra, dec) / whatever units Image.t is in
    S::Poisson  # number of asteroids in the collection of images
end

