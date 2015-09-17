const dat_dir = joinpath(Pkg.dir("KillerAsteroids"), "dat")


function load_wise_psf(band_id::Int64)
    # TODO: sensible normalization
    # TODO: checks on the input values(?)

    psf_fnames = ["psf_coeff-w1.v1.fits",
                  "psf_coeff-w2.v1.fits",
                  "psf_coeff-big.fits",
                  "psf_coeff-taper-w4.fits" ]

    psf_fname = psf_fnames[band_id]
    psf_fname = joinpath(dat_dir, psf_fname)
    
    f = FITS(psf_fname)
    psf_model = read(f[1])
    close(f)

    # for now just return the zeroth order term, add higher orders later
    convert(Array{Float64}, psf_model[:,:,1])
end


"""
Returns a PSF cutout returned with dimensions 
  (2*halfsidelen+1, 2*halfsielen+1).
"""
function load_wise_psf(band_id::Int64, halfsidelen::Int64)
    psf_model = load_wise_psf(band_id)

    sz = size(psf_model, 1)
    @assert sz == size(psf_model, 2)
    # will break if 2*halfsidelen+1 is larger than PSF model sidelength
    @assert sz > 2halfsidelen

    pixcen = ceil(Int, sz / 2)

    pixrange = (pixcen-halfsidelen):(pixcen+halfsidelen)
    psf_model[pixrange, pixrange]
end

