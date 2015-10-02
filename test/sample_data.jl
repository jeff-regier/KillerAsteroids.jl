using Distributions
using WCSLIB


const sample_wcs = wcsprm(2; # naxis
           cdelt = [-0.066667, 0.066667],
           ctype = ["RA---AIR", "DEC--AIR"],
           crpix = [0.75, 0.3393],
           crval = [0., 0],
           pv    = [pvcard(2, 1, 5.0)])

const sample_ast = begin
    # this asteroids moves one pixel down and one pixel right per second
    u0 = wcsp2s(sample_wcs, [20 12.]')[:]
    u1 = wcsp2s(sample_wcs, [19 13.]')[:]
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

function generate_sample_image(; ast=sample_ast, psf=sample_psf, 
    H=30, W=30, band_id=2, t=0., wcs=sample_wcs)

    sky_noise_mean = 40  # approx value for W2, in DN
    nmgy_per_dn = wise_band_to_params[band_id].nmgy_per_dn
    read_noise_var = wise_band_to_params[band_id].read_noise_var

    pixel_var = sky_noise_mean + read_noise_var
    sky_read_rv = Normal(sky_noise_mean, sqrt(pixel_var))
    pixels = rand(sky_read_rv, H, W)

    ast_r_dn = ast.r[band_id] / nmgy_per_dn

    u_t = extrapolate_position(ast.u, ast.v, t)
    u_t_px_crd = wcss2p(wcs, u_t'')
    u_t_px = round(Int, u_t_px_crd)  # the psf is constant per pixel for now

    @assert size(psf) == (5, 5)
    for w2 in 1:5, h2 in 1:5
        h = round(Int, u_t_px[1]) + h2 - 3
        w = round(Int, u_t_px[2]) + w2 - 3
        ast_r_pixel_dn =  ast_r_dn * psf[h2, w2]
        pixels[h, w] += rand() * sqrt(ast_r_pixel_dn) + ast_r_pixel_dn
    end

    Image(H, W, pixels, sky_noise_mean, psf, band_id, t, wcs)
end

const sample_img = generate_sample_image()

const sample_prior = Prior(
    # r_prior: mean = 7332; median = 4447; mode=1635; sd=9611
    MvNormal(fill(8.4, B), 1), # per band log brightness
    MvNormal([3., 5], [4. 1; 3 2]), # velocity
    Poisson(1.3)) # number of sources

