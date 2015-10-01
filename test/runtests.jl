using KillerAsteroids
using Base.Test
using FITSIO

include("sample_data.jl")


function test_truth_most_likely_with_all_synthetic_data()
    test_img = generate_sample_image()

    # it's a pretty bright asteroid, so even with the psf,
    # the pixel that contains the asteroid should be the brightest
    @test_approx_eq test_img.pixels[20, 12] maximum(test_img.pixels)

    prior = sample_prior()

    good_ll = compute_log_probability([good_ast,], [test_img,], prior)

    bad_ast = AsteroidParams(good_ast.r, [19.4, 12.], good_ast.v)
    bad_ll = compute_log_probability([bad_ast,], [test_img,], prior)

    info("$good_ll > $bad_ll")
    @test good_ll > bad_ll
end


function test_truth_most_likely_with_wise_psf()
    band_id = 3
    halfsidelen = 2
    psf = load_wise_psf(band_id, halfsidelen) # sidelength will be 2*2 + 1
    psf /= sum(psf)

    test_img = generate_sample_image(psf)

    # the pixel that contains the asteroid should be the brightest
    @test_approx_eq test_img.pixels[20, 12] maximum(test_img.pixels)

    prior = sample_prior()

    good_ll = compute_log_probability([good_ast,], [test_img,], prior)

    bad_ast = AsteroidParams(good_ast.r, [19.4, 12.], good_ast.v)
    bad_ll = compute_log_probability([bad_ast,], [test_img,], prior)

    info("$good_ll > $bad_ll")
    @test good_ll > bad_ll
end

function test_truth_most_likely_with_real_bright_asteroid(band_id::Int64)
    halfsidelen = 2
    psf = load_wise_psf(band_id, halfsidelen) # sidelength will be 2*2 + 1
    psf /= sum(psf)

    fdir = joinpath(Pkg.dir("KillerAsteroids"), "dat")
    fname = string("stereoskopia_w", string(band_id),".fits")
    fname = joinpath(fdir, fname)

    f = FITS(fname)
    pixels = read(f[1])

    sz = size(pixels)

    H = sz[1]
    W = sz[2]
    nmgy_per_dn = [5.0, 14.5]
    sky_noise_mean = [22, 53]
    read_noise_var = [9.95, 7.78]
    gain = [3.75, 4.60]

    real_img = Image(H, W, pixels, nmgy_per_dn[band_id], 
                     sky_noise_mean[band_id], read_noise_var[band_id], 
                     gain[band_id], psf, band_id, 0.)

    prior = sample_prior()

    ast_flux_nmgy = [10729.9, 47165.5]
    ast = AsteroidParams(ast_flux_nmgy[band_id], [15., 15.], [0., 0.])
    good_ll = compute_log_probability([ast,], [real_img,], prior)

    # try an asteroid with the right brightness but wrong position
    bad_ast = AsteroidParams(ast.r, [20., 12.], ast.v)
    bad_ll = compute_log_probability([bad_ast,], [real_img,], prior)

    info("$good_ll > $bad_ll")
    @test good_ll > bad_ll

    # try an asteroid with the right position but too faint
    bad_ast = AsteroidParams(0.1*ast.r,ast.u, ast.v)
    bad_ll = compute_log_probability([bad_ast,], [real_img,], prior)

    info("$good_ll > $bad_ll")
    @test good_ll > bad_ll

    # try an asteroid with the right position but too bright
    bad_ast = AsteroidParams(10.0*ast.r,ast.u, ast.v)
    bad_ll = compute_log_probability([bad_ast,], [real_img,], prior)

    info("$good_ll > $bad_ll")
    @test good_ll > bad_ll

end

# thats the actual path for asteroid 2005_UT453 has the highest
# probability according to our model
function test_truth_most_likely_with_all_real_data()
    band_id = 2
    halfsidelen = 2
    psf = load_wise_psf(band_id, halfsidelen)
    psf /= sum(psf)
    @test(minimum(psf) > 0)

    # Possible steps to writing this test:
    #
    # 1. images = <load images of 2005_UT453 taken in band 2>
    # 2. good_ast = <load the true band 2 brigthness of the 2005_UT453,
    #           its true position (in pixel coordinates at time 0)
    #           and its true velocity during the image exposures>
    #           Note: time 0 can be time of the first image in the stack
    # 3. good_ll compute_log_prob([good_ast,], images, prior), where
    #            prior is set like it's done in sample_prior(), but with
    #            more realistic values
    # 4. for a variety of bad_ast, some perhaps not all that different
    #    from good_ast, do 
    #             bad_ll = compute_log_prob([bad_ast,], images, prior)
    #    and
    #          @test bad_ll < good_ll
end


function test_variable_numbers_of_asteroids()
    test_img = generate_sample_image()
    prior = sample_prior()
    good_ll = compute_log_probability([good_ast,], [test_img,], prior)

    bad_ll_0 = compute_log_probability(AsteroidParams[], [test_img,], prior)
    info("$good_ll > $bad_ll_0")
    @test good_ll > bad_ll_0

    bad_ll_2 = compute_log_probability([good_ast, good_ast], [test_img,], prior)
    info("$good_ll > $bad_ll_2")
    @test good_ll > bad_ll_2
end


test_truth_most_likely_with_all_synthetic_data()
test_truth_most_likely_with_wise_psf()
test_truth_most_likely_with_real_bright_asteroid(2)
test_truth_most_likely_with_all_real_data()
test_variable_numbers_of_asteroids()

