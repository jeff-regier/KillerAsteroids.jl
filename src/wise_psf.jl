using FITSIO

function wise_psf(band::Int64, halfsidelen::Int64)

  return wise_psf_sidelen(band, halfsidelen)
end

function wise_psf(band::Int64)

  # if no halfsidelen value has been specified, assume user wants full PSF

  sidelens = [325, 325, 499, 285]
  halfsidelen = floor(Int, sidelens[band]/2)

  return wise_psf_sidelen(band, halfsidelen)
end

function wise_psf_sidelen(band::Int64, halfsidelen::Int64)

  # PSF cutout returned will have dimensions (2*halfsidelen+1, 2*halfsielen+1)

  # don't believe any sensible normalization is in place, need to fix that ...

  # could potentially add checks on the input values ...

  # figure out where the PSF model FITS file is located
  psf_dir = ENV["WISE_ETC"]

  psf_fnames = ["psf_coeff-w1.v1.fits",
                "psf_coeff-w2.v1.fits",
                "psf_coeff-big.fits",
                "psf_coeff-taper-w4.fits" ]

  psf_fname = psf_fnames[band]
  psf_fname = joinpath(psf_dir, psf_fname)
  
  f = FITS(psf_fname)
  psf_model = read(f[1])
  close(f)

  # for now just return the zeroth order term, add higher orders later
  psf_model = psf_model[:,:,1]

  sz = size(psf_model, 1)
  @assert 2*halfsidelen < sz
  pixcen = ceil(Int, sz/2)
  # will break if 2*halfsidelen+1 is larger than PSF model sidelength
  psf_model = psf_model[(pixcen-halfsidelen):(pixcen+halfsidelen),
                        (pixcen-halfsidelen):(pixcen+halfsidelen)]

  return psf_model
end