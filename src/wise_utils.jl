function wise_params(band_id::Int64)

    # TODO: generalize to handle W3 and W4
    # TODO: add more parameters as they become relevant

    nmgy_per_dn = [5.0, 14.5]
    read_noise_var = [9.95, 7.78]
    gain = [3.75, 4.60]

    params = ["nmgy_per_dn" => nmgy_per_dn[band_id],
              "read_noise_var" => read_noise_var[band_id], 
              "gain" => gain[band_id]]

end