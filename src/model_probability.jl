
function compute_log_prior(asteroids::Vector{AsteroidParams}, prior::Prior)
    lp = logpdf(prior.S, length(asteroids))
    for ast in asteroids
        lp += logpdf(prior.log_r, log(ast.r))
        lp += logpdf(prior.v, ast.v)
    end
    lp
end

"""
Returns the position of the asteroid at time t.
Note: Can make this a nonlinear function here in the future, once the time
of integration is months or years rather than hours.
args
  - u0: position of the asteroid at time 0.
  - v: velocity of the asteroid
  - t: time
"""
function extrapolate_position(u0::Vector{Float64}, v::Vector{Float64}, t::Float64)
    u0 + v * t 
end

function compute_log_likelihood(asteroids::Vector{AsteroidParams},
        images::Vector{Image})
    ll = 0.

    for img in images
        psf_dims = size(img.psf)
        psf_center = [round(Int, (dim + 1) / 2) for dim in size(img.psf)]
        psf_itp = interpolate(img.psf, BSpline(Linear()), OnGrid())
        expected_dn = similar(img.pixels)
        fill!(expected_dn, img.sky_noise_mean)

        read_noise_var = wise_band_to_params[img.band_id].read_noise_var
        nmgy_per_dn = wise_band_to_params[img.band_id].nmgy_per_dn

        for ast in asteroids
            u_t = extrapolate_position(ast.u, ast.v, img.t)
            u_t_px_crd = wcss2p(img.wcs, u_t'')
            u_t_px = round(Int, u_t_px_crd)
            offset = u_t_px_crd - u_t_px

            ast_r_dn = ast.r[img.band_id] / nmgy_per_dn
            for w2 in 2:(psf_dims[2]-1), h2 in 2:(psf_dims[1]-1)
                h = u_t_px[1] + h2 - psf_center[1]
                w = u_t_px[2] + w2 - psf_center[2]
                if (h > img.H) || (w > img.W) || (h < 1) || (w < 1)  
                    continue
                end
                h3 = h2 + offset[1]
                w3 = w2 + offset[2]
                expected_ast_dn = ast_r_dn * psf_itp[h3, w3]
                if (psf_itp[h3, w3] < 0)
                    println(h3, "  ", w3)
                    println(img.psf)
                end
                @assert(psf_itp[h3, w3] >= 0)
                expected_dn[h, w] += expected_ast_dn
            end
        end

        for w in 1:img.W, h in 1:img.H
            pixel_dn_var = expected_dn[h, w] + read_noise_var
            pixel_dist = Normal(expected_dn[h, w], sqrt(pixel_dn_var))
            ll += logpdf(pixel_dist, img.pixels[h, w])
        end
    end

    ll
end



"""
Computes the unnormalized log probability for a particular
candidate asteroid (higher is better)

arguments:
  ast: parameters for a candidate asteroid
  img: an astronomical image
"""
function compute_log_probability(asteroids::Vector{AsteroidParams},
        images::Vector{Image}, prior::Prior)
    lp = compute_log_prior(asteroids, prior)
    ll = compute_log_likelihood(asteroids, images)

    lp + ll
end

