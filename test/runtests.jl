using KillerAsteroids
using Base.Test
using FITSIO

include("sample_data.jl")

const dat_dir = joinpath(Pkg.dir("KillerAsteroids"), "dat")


function test_truth_most_likely_with_all_synthetic_data()
    # it's a pretty bright asteroid, so even with the psf,
    # the pixel that contains the asteroid should be the brightest
    @test_approx_eq sample_img.pixels[20, 12] maximum(sample_img.pixels)

    good_lp = compute_log_probability([sample_ast,], [sample_img,], sample_prior)

    bad_u_px_crd = [19.4, 12.]
    bad_u = wcsp2s(sample_img.wcs, bad_u_px_crd'')[:]
    bad_ast = AsteroidParams(sample_ast.r, bad_u, sample_ast.v)
    bad_lp = compute_log_probability([bad_ast,], [sample_img,], sample_prior)

    info("$good_lp > $bad_lp")
    @test good_lp > bad_lp
end


function test_truth_most_likely_with_wise_psf()
    band_id = 3
    halfsidelen = 2
    psf = load_wise_psf(band_id, halfsidelen) # sidelength will be 2*2 + 1
    psf /= sum(psf)

    test_img = generate_sample_image(psf=psf)

    # the pixel that contains the asteroid should be the brightest
    @test_approx_eq test_img.pixels[20, 12] maximum(test_img.pixels)

    good_lp = compute_log_probability([sample_ast,], [test_img,], sample_prior)

    bad_u_px_crd = [19.4, 12.]
    bad_u = wcsp2s(test_img.wcs, bad_u_px_crd'')[:]
    bad_ast = AsteroidParams(sample_ast.r, bad_u, sample_ast.v)
    bad_lp = compute_log_probability([bad_ast,], [test_img,], sample_prior)

    info("$good_lp > $bad_lp")
    @test good_lp > bad_lp
end


function test_truth_most_likely_with_real_bright_asteroid(band_id::Int64)
    halfsidelen = 2
    psf = load_wise_psf(band_id, halfsidelen) # sidelength will be 2*2 + 1
    psf /= sum(psf)

    fname = joinpath(dat_dir, "stereoskopia_w$band_id.fits")
    f = FITS(fname)
    pixels = read(f[1])
    header_str = read_header(f[1], ASCIIString)
    close(f)

    H, W = size(pixels)
    par = wise_band_to_params[band_id]
    sky_noise_mean = [22., 53.]
    wcs = wcspih(header_str)[1][1]

    real_img = Image(H, W, pixels, sky_noise_mean[band_id], psf, 
                        band_id, 0., wcs)

    # TODO: fill in bands 3 and 4 with their actual real values
    ast_flux_nmgy = [10729.9, 47165.5, 33333., 33333.]
    ast_u_px_crd = [15., 15.]
    ast_u = wcsp2s(wcs, ast_u_px_crd'')[:]
    ast = AsteroidParams(ast_flux_nmgy, ast_u, [0., 0.])
    good_lp = compute_log_probability([ast,], [real_img,], sample_prior)

    # try an asteroid with the right brightness but wrong position
    bad_u_px_crd = [20, 12.]
    bad_u = wcsp2s(wcs, bad_u_px_crd'')[:]
    bad_ast = AsteroidParams(ast.r, bad_u, ast.v)
    bad_lp = compute_log_probability([bad_ast,], [real_img,], sample_prior)

    info("$good_lp > $bad_lp")
    @test good_lp > bad_lp

    # try an asteroid with the right position but too faint
    bad_ast = AsteroidParams(0.1 * ast.r, ast.u, ast.v)
    bad_lp = compute_log_probability([bad_ast,], [real_img,], sample_prior)

    info("$good_lp > $bad_lp")
    @test good_lp > bad_lp

    # try an asteroid with the right position but too bright
    bad_ast = AsteroidParams(10.0 * ast.r, ast.u, ast.v)
    bad_lp = compute_log_probability([bad_ast,], [real_img,], sample_prior)

    info("$good_lp > $bad_lp")
    @test good_lp > bad_lp

end


function test_variable_numbers_of_asteroids()
    good_lp = compute_log_probability([sample_ast,], [sample_img,], sample_prior)

    bad_lp_0 = compute_log_probability(AsteroidParams[], [sample_img,], sample_prior)
    info("$good_lp > $bad_lp_0")
    @test good_lp > bad_lp_0

    bad_lp_2 = compute_log_probability([sample_ast, sample_ast], 
                    [sample_img,], sample_prior)
    info("$good_lp > $bad_lp_2")
    @test good_lp > bad_lp_2
end


function test_asteroid_partially_off_edge()
    band_id = 1
    halfsidelen = 10
    psf = load_wise_psf(band_id, halfsidelen)
    psf /= sum(psf)

    fname = joinpath(dat_dir, "stereoskopia_w$band_id.fits")
    f = FITS(fname)
    pixels = read(f[1])
    header_str = read_header(f[1], ASCIIString)
    close(f)

    H, W = size(pixels)
    par = wise_band_to_params[band_id]
    sky_noise_mean = 22.
    wcs = wcspih(header_str)[1][1]

    real_img = Image(H, W, pixels, sky_noise_mean, psf, band_id, 0., wcs)

    edge_u_px_crd = [25., 15.]
    edge_u = wcsp2s(wcs, edge_u_px_crd'')[:]
    sample_r = fill(10729.9, B)
    edge_ast = AsteroidParams(sample_r, edge_u, [0., 0.])
    edge_lp = compute_log_probability([edge_ast,], [real_img,], sample_prior)
end


function test_truth_most_likely_with_multiple_images()
    times = [0., 1, 3, 3.5]
    images = [generate_sample_image(t=times[i], band_id=(1+i%2)) for i in 1:4]
    correct_ast_lp = compute_log_probability([sample_ast,], [sample_img,],
                        sample_prior)

    wrong_u_px_crd = [19.4, 12.]
    wrong_u = wcsp2s(sample_img.wcs, wrong_u_px_crd'')[:]
    wrong_u_ast = AsteroidParams(sample_ast.r, wrong_u, sample_ast.v)
    wrong_u_lp = compute_log_probability([wrong_u_ast,], images, sample_prior)
    info("$correct_ast_lp > $wrong_u_lp")
    @test correct_ast_lp > wrong_u_lp

    wrong_v_ast = AsteroidParams(sample_ast.r, sample_ast.u, [0., 0.])
    wrong_v_lp = compute_log_probability([wrong_v_ast,], images, sample_prior)
    info("$correct_ast_lp > $wrong_v_lp")
    @test correct_ast_lp > wrong_v_lp

    wrong_v_ast2 = AsteroidParams(sample_ast.r, sample_ast.u, -1 * sample_ast.v)
    wrong_v_lp2 = compute_log_probability([wrong_v_ast2,], images, sample_prior)
    info("$correct_ast_lp > $wrong_v_lp2")
    @test correct_ast_lp > wrong_v_lp2
end

function test_image_file_name()

  band_id = 1
  scan_id = "00918a"
  frame_num = 148

  fname = l1b_image_name(band_id, scan_id, frame_num, dat_dir)

  @test isfile(fname)

  f = FITS(fname)
  l1b_image = read(f[1])

  par = wise_band_to_params[band_id]

  @test size(l1b_image) == (par.l1b_sidelength, par.l1b_sidelength)

  mskname = l1b_image_name(band_id, scan_id, frame_num, dat_dir; im_type="msk")

  println(mskname)

  @test isfile(mskname)

  f = FITS(mskname)
  l1b_msk = read(f[1])

  @test size(l1b_msk) == (par.l1b_sidelength, par.l1b_sidelength)

end

# thats the actual path for asteroid 2005_UT453 has the highest
# probability according to our model
function test_truth_most_likely_with_multiple_real_images()
    band_id = 2
    halfsidelen = 2
    psf = load_wise_psf(band_id, halfsidelen)
    psf /= sum(psf)
    @test(minimum(psf) > 0)

    # Possible steps to writing this test:
    #
    # 1. images = <load images of 2005_UT453 taken in band 2>
    # 2. sample_ast = <load the true band 2 brigthness of the 2005_UT453,
    #           its true position (in pixel coordinates at time 0)
    #           and its true velocity during the image exposures>
    #           Note: time 0 can be time of the first image in the stack
    # 3. good_lp compute_log_prob([sample_ast,], images, prior), where
    #            prior is set like it's done in sample_prior, but with
    #            more realistic values
    # 4. for a variety of bad_ast, some perhaps not all that different
    #    from sample_ast, do 
    #             bad_lp = compute_log_prob([bad_ast,], images, prior)
    #    and
    #          @test bad_lp < good_lp
end


test_truth_most_likely_with_all_synthetic_data()
test_truth_most_likely_with_wise_psf()
test_truth_most_likely_with_real_bright_asteroid(1)
test_truth_most_likely_with_real_bright_asteroid(2)
test_variable_numbers_of_asteroids()
test_asteroid_partially_off_edge()
test_truth_most_likely_with_multiple_images()
test_truth_most_likely_with_multiple_real_images()
test_image_file_name()
