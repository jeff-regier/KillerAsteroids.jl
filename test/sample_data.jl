using Distributions


const good_ast = AsteroidParams(10_000., [20, 12.], [3.1, 5.1])

function generate_sample_image(psf::Matrix{Float64})
    H, W = 30, 30
    # use the  band 2 value : 
    nmgy_per_dn = 14.5  # approximate value for W2
    sky_noise_mean = 40  # approx value for W2, in DN
    read_noise_var = 7.78  # approx value for W2, in DN^2
    # Not sure if gain matters...the Poisson noise is in DN, not in number
    # of photoelectrons...
    gain = 4.60

    pixel_var = sky_noise_mean + read_noise_var
    sky_read_rv = Normal(sky_noise_mean, sqrt(pixel_var))
    pixels = rand(sky_read_rv, H, W)

    ast_r_dn = good_ast.r / nmgy_per_dn
    @assert size(psf) == (5, 5)
    for w2 in 1:5, h2 in 1:5
        h = round(Int, good_ast.u[1]) + h2 - 3
        w = round(Int, good_ast.u[2]) + w2 - 3
        ast_r_pixel_dn =  ast_r_dn * psf[h2, w2]
        pixels[h, w] += rand() * sqrt(ast_r_pixel_dn) + ast_r_pixel_dn
    end

    Image(H, W, pixels, nmgy_per_dn, sky_noise_mean, read_noise_var,
        gain, psf, 2, 0.)
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


