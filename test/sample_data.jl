using Distributions
using WCSLIB

const sample_wcs = wcsprm(2; # naxis
           cdelt = [-0.066667, 0.066667],
           ctype = ["RA---AIR", "DEC--AIR"],
           crpix = [0.75, 0.3393],
           crval = [0., 0],
           pv    = [pvcard(2, 1, 5.0)])

const sample_ast = begin
    # this asteroids moves one pixel up and one pixel over per second
    u0 = wcsp2s(sample_wcs, [8 11.]')[:]
    u1 = wcsp2s(sample_wcs, [9 12.]')[:]
    r = fill(10_000., B) # brightness
    AsteroidParams(r, u0, u1 - u0)
end

const sample_psf = begin
    bvn_for_psf = MvNormal([2.99, 3.01], [0.2 0; 0 0.2])
    ret = Array(Float64, 5, 5)
    for w2 in 1:5, h2 in 1:5
        ret[h2, w2] = pdf(bvn_for_psf, [h2, w2])
    end
    ret /= sum(ret)
end

function generate_sample_image(ast=sample_ast, 
        psf=sample_psf, band_id=2, t=0., wcs=sample_wcs)
    H, W = 30, 30

    # TODO: look these up based on the band_id
    # use the band 2 values:
    nmgy_per_dn = 14.5  # approximate value for W2
    sky_noise_mean = 40  # approx value for W2, in DN
    read_noise_var = 7.78  # approx value for W2, in DN^2
    # Not sure if gain matters...the Poisson noise is in DN, not in number
    # of photoelectrons...
    gain = 4.60

    pixel_var = sky_noise_mean + read_noise_var
    sky_read_rv = Normal(sky_noise_mean, sqrt(pixel_var))
    pixels = rand(sky_read_rv, H, W)

    ast_r_dn = ast.r / nmgy_per_dn

    u_t = extrapolate_position(ast.u, ast.v, t)
    u_t_px_crd = wcss2p(wcs, u_t'')
    u_t_px = round(Int, u_t)  # the psf is constant per pixel for now

    @assert size(psf) == (5, 5)
    for w2 in 1:5, h2 in 1:5
        h = round(Int, u_t_px[1]) + h2 - 3
        w = round(Int, u_t_px[2]) + w2 - 3
        ast_r_pixel_dn =  ast_r_dn * psf[h2, w2]
        pixels[h, w] += rand() * sqrt(ast_r_pixel_dn) + ast_r_pixel_dn
    end

    Image(H, W, pixels, nmgy_per_dn, sky_noise_mean, read_noise_var,
        gain, psf, band_id, t, wcs)
end

const sample_img = generate_sample_image()

const sample_prior = Prior(
    # r_prior: mean = 7332; median = 4447; mode=1635; sd=9611
    MvnNormal(fill(8.4, 5), 1), # per band brightness
    MvNormal([3., 5], [4. 1; 3 2]) # velocity
    Poisson(1.3)) # number of sources

