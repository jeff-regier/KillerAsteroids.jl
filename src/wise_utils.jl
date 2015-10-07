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

function l1b_image_name(band_id::Int64, scan_id::String, frame_num::Int64, 
                      l1b_dir::String; im_type::String="int")

    @assert (im_type == "int") || (im_type == "msk") || (im_type == "unc")

    fname = string(scan_id, lpad(frame_num, 3, "0"), "-w", string(band_id), "-", im_type, "-1b.fits")
    if im_type != "int"
        fname = string(fname, ".gz")
    end
    fname = joinpath(l1b_dir, fname)
end

function l1b_image_dir(scan_id::String, frame_num::Int64, basedir::String)

    subdir = lpad(frame_num, 3, "0")
    dir = joinpath(basedir, scan_id[5:6], scan_id, subdir)

end
