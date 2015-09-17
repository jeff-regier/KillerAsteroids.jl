using KillerAsteroids
using Base.Test

include("sample_data.jl")


function test_truth_most_likely_with_all_synthetic_data()
    test_img = generate_sample_image()

    # it's a pretty bright asteroid, so even with the psf,
    # the pixel that contains the asteroid should be the brightest
    @test_approx_eq test_img.pixels[20, 12] maximum(test_img.pixels)

    prior = sample_prior()

    good_ast = AsteroidParams(10., [20, 12.], [3.1, 5.1])
    good_ll = compute_log_probability(good_ast, test_img, prior)

    bad_ast = AsteroidParams(10., [19.4, 12.], [3.1, 5.1])
    bad_ll = compute_log_probability(bad_ast, test_img, prior)

    info("$good_ll > $bad_ll")
    @test good_ll > bad_ll
end


function test_truth_most_likely_with_wise_psf()
    band_id = 3
    halfsidelen = 2
    psf = load_wise_psf(band_id, halfsidelen) # sidelength will be 2*2 + 1
    psf /= sum(psf)

    test_img = generate_sample_image(psf)
    prior = sample_prior()

    good_ast = AsteroidParams(10., [20, 12.], [3.1, 5.1])
    good_ll = compute_log_probability(good_ast, test_img, prior)

    bad_ast = AsteroidParams(10., [19.4, 12.], [3.1, 5.1])
    bad_ll = compute_log_probability(bad_ast, test_img, prior)

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
    # 1. image_stack = <load images of 2005_UT453 taken in band 2>
    # 2. good_ast = <load the true band 2 brigthness of the 2005_UT453,
    #           its true position (in pixel coordinates at time 0)
    #           and its true velocity during the image exposures>
    #           Note: time 0 can be time of the first image in the stack
    # 3. good_ll compute_log_prob(good_ast, image_stack, prior), where
    #            prior is set like it's done in sample_prior(), but with
    #            more realistic values
    # 4. for a variety of bad_ast, some perhaps not all that different
    #    from good_ast, do 
    #             bad_ll = compute_log_prob(bad_ast, image_stack, prior)
    #    and
    #          @test bad_ll < good_ll
end


test_truth_most_likely_with_all_synthetic_data()
test_truth_most_likely_with_wise_psf()
test_truth_most_likely_with_all_real_data()
