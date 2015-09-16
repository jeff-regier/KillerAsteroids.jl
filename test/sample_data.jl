using Distributions


function wrapped_poisson(rate::Float64)
    0 < rate ? float(rand(Poisson(rate))) : 0.
end


function generate_sample_image(psf::Matrix{Float64})
    H, W = 30, 30
    nmgy_per_dn = 1 / 112
    sky_noise_mean = 2.2
    read_noise_var = 0.3
    gain = 4.4  # Not sure if gain matters...pixel fluxes in DN are Poisson, 
                #right?, not pixels fluxes nn # of photoelectrons?

    pixel_var = sky_noise_mean + read_noise_var
    sky_read_rv = Normal(sky_noise_mean, sqrt(pixel_var))
    pixels = rand(sky_read_rv, H, W)

    ast = AsteroidParams(10.3, [20, 12.], [3.1, 5.1])
    ast_r_dn = ast.r / nmgy_per_dn
    for w2 in 1:5, h2 in 1:5
        h = round(Int, ast.u[1]) + h2 - 3
        w = round(Int, ast.u[2]) + w2 - 3
        ast_r_pixel_dn =  ast_r_dn * psf[h2, w2]
        pixels[h, w] += rand() * sqrt(ast_r_pixel_dn) + ast_r_pixel_dn
    end

    Image(H, W, pixels, nmgy_per_dn, sky_noise_mean, read_noise_var, gain, psf, 0.)
end


function generate_sample_image()
    bvn_for_psf = MvNormal([2.99, 3.01], [0.2 0; 0 0.2])
    psf = Array(Float64, 5, 5)
    for w2 in 1:5, h2 in 1:5
        psf[h2, w2] = pdf(bvn_for_psf, [h2, w2])
    end
    psf /= sum(psf)

    generate_sample_image(psf)
end


function sample_prior()
    r_prior = LogNormal(15., sqrt(50)) # mode is 1000
    v_prior = MvNormal([3., 5], [4. 1; 3 2]) # pixels / second

    Prior(r_prior, v_prior)
end


