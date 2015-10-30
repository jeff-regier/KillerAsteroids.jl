function mask_dilate(mask::Array{Bool, 2}, kern::Array{Bool, 2})

  sz_kern = size(kern)

  # dilation kernel needs to be square
  @assert sz_kern[1] == sz_kern[2]

  # dilation kernel needs to have odd sidelength
  @assert (sz_kern[1] % 2) == 1

  # dilation kernel needs to have nonzero total 'weight'
  @assert sum(kern) != 0

  w = find(mask)
  nw = length(w)

  # if the mask is empty then do nothing
  if nw == 0
      return mask
  end

  sidelen_kern = sz_kern[1]

  sz_mask = size(mask)
  nx = sz_mask[1]
  ny = sz_mask[2]

  x_mask = ((w-1) % nx) + 1
  y_mask = ceil(Int64, w / ny)

  wkern = find(kern)
  nwkern = length(wkern)
  if nwkern == 0
      return 0*mask
  end

  dx_kern = ((wkern-1) % sidelen_kern) - floor(Int64, sidelen_kern/2)
  dy_kern = floor(Int64, (wkern-1)/ sidelen_kern) - floor(Int64, sidelen_kern/2)

  x_add_mask = zeros(Int64, nw, nwkern)
  y_add_mask = zeros(Int64, nw, nwkern)

  for i=1:nw
      x_add_mask[i, :] = x_mask[i] + dx_kern
      y_add_mask[i, :] = y_mask[i] + dy_kern
  end

  x_add_mask = reshape(x_add_mask, length(x_add_mask))
  y_add_mask = reshape(y_add_mask, length(y_add_mask))

  wgood = find((x_add_mask .> 0) & (x_add_mask .<= nx) &
               (y_add_mask .> 0) & (y_add_mask .<= ny))

  nwgood = length(wgood)
  if nwgood > 0
     out = copy(mask)
     [(out[x_add_mask[wgood[i]], y_add_mask[wgood[i]]] = true) for i=1:nwgood]
  end
  return out

end