immutable WiseBandParams

    # TODO: generalize to handle W3 and W4
    # TODO: add more parameters as they become relevant

    band_id::Int64
    nmgy_per_dn::Float64
    read_noise_var::Float64
    gain::Float64

end

const wise_band_to_params = [
   WiseBandParams(1, 5.0, 9.95, 3.75),
   WiseBandParams(2, 14.5, 7.78, 4.6)
]
