
function compute_log_prior(asteroids::Vector{AsteroidParams}, prior::Prior)
    lp = logpdf(prior.S, length(asteroids))
    for ast in asteroids
        lp += logpdf(prior.r, ast.r)
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
        expected_dn = similar(img.pixels)
        fill!(expected_dn, img.sky_noise_mean)

        for ast in asteroids
            u_t = extrapolate_position(ast.u, ast.v, img.t)
            u_t_px = round(Int, u_t)
            ast_r_dn = ast.r / img.nmgy_per_dn
            for w2 in 1:psf_dims[2], h2 in 1:psf_dims[1]
                h = u_t_px[1] + h2 - psf_center[1]
                w = u_t_px[2] + w2 - psf_center[2]
                expected_ast_dn = ast_r_dn * img.psf[h2, w2]
                expected_dn[h, w] += expected_ast_dn
            end
        end

        for w in 1:img.W, h in 1:img.H
            pixel_dn_var = expected_dn[h, w] + img.read_noise_var
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

