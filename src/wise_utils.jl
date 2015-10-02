immutable WiseBandParams

    # TODO: generalize to handle W3 and W4
    # TODO: add more parameters as they become relevant

    nmgy_per_dn::Float64
    read_noise_var::Float64
    gain::Float64

end

const wise_band_to_params = [
   WiseBandParams(5.0, 9.95, 3.75),
   WiseBandParams(14.5, 7.78, 4.6)
]

# total number of bands in WISE
const B = 4

