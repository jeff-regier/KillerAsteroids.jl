using KillerAsteroids
using Base.Test


function test_truth_most_likely_with_all_synthetic_data()
    const test_img = Synthetic.generate_sample_image()
    @test test_img.pixels[20, 12] > 300

    prior = ModelInit.sample_prior()

    const good_ast = AsteroidParams(1000., [20, 12.], [3.1, 5.1])
    good_ll = Probability.compute_log_probability(good_ast, test_img, prior)

    const bad_ast = AsteroidParams(1000., [19.4, 12.], [3.1, 5.1])
    bad_ll = Probability.compute_log_probability(bad_ast, test_img, prior)

    info("$good_ll > $bad_ll")
    @test good_ll > bad_ll
end


function test_truth_most_likely_with_wise_psf()
    band_id = 3
    halfsidelen = 2
    psf = WisePSF.load_psf(band_id, halfsidelen) # sidelength will be 2*2 + 1
    psf /= sum(psf)

    const test_img = Synthetic.generate_sample_image(psf)
    @test test_img.pixels[20, 12] > 300

    prior = ModelInit.sample_prior()

    const good_ast = AsteroidParams(1000., [20, 12.], [3.1, 5.1])
    good_ll = Probability.compute_log_probability(good_ast, test_img, prior)

    const bad_ast = AsteroidParams(1000., [19.4, 12.], [3.1, 5.1])
    bad_ll = Probability.compute_log_probability(bad_ast, test_img, prior)

    info("$good_ll > $bad_ll")
    @test good_ll > bad_ll
end


function test_truth_most_likely_with_all_real_data()
    # TODO: test that the actual path for asteroid 2005_UT453 has the highest
    # probability according to our model
end


test_truth_most_likely_with_all_synthetic_data()
#TODO: comment this test in once data from WISE_ETC is in our `dat` directory
#test_truth_most_likely_with_wise_psf()
test_truth_most_likely_with_all_real_data()
