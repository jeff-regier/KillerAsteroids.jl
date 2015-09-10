using FITSIO

function wise_psf(band; sidelen=-1)

  # don't believe any sensible normalization is in place, need to fix that ...

  # could potentially add checks on the input values ...

  # figure out where the PSF model FITS file is located
  psf_dir = ENV["WISE_ETC"]

  psf_fname_dict = [1=>"psf_coeff-w1.v1.fits",
                    2=>"psf_coeff-w2.v1.fits",
                    3=>"psf_coeff-big.fits",
                    4=>"psf_coeff-taper-w4.fits" ]

  psf_fname = psf_fname_dict[band]
  psf_fname = string(psf_dir, "/", psf_fname)
  
  f = FITS(psf_fname)
  psf_model = read(f[1])

  # for now just return the zeroth order term, add higher orders later
  psf_model = psf_model[:,:,1]

  # even sidelengths won't really make sense in my opinion as
  # middle of central pixel will not be the center of the PSF in that case
  if (sidelen != -1) & ((sidelen % 2) != 0)
      # 1-indexing ...
      sz = size(psf_model, 1)
      pixcen = int(ceil(sz/2))
      half = int(floor(sidelen/2))
      # will break if sidelen is greater than PSF model size ...
      psf_model = psf_model[(pixcen-half):(pixcen+half),
                            (pixcen-half):(pixcen+half)]
  end

  return psf_model
end