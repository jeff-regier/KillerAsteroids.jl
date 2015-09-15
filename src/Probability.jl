module Probability

VERSION < v"0.4.0-dev" && using Docile

using Distributions
using ModelInit


@doc """
Computes the unnormalized log likelihood for a particular
candidate asteroid (higher is better)

arguments:
  ast: parameters for a candidate asteroid
  img: an astronomical image
""" ->
function compute_log_probability(ast::AsteroidParams, img::Image, prior::Prior)
    log_prior = logpdf(prior.r, ast.r) * logpdf(prior.v, ast.v)

    log_like = 0.
    for w2 in 1:5, h2 in 1:5
        h = round(Int, ast.u[1]) + h2 - 3
        w = round(Int, ast.u[2]) + w2 - 3
        expected_x = img.iota * (img.epsilon + ast.r * img.psf[h2, w2])
        log_like += logpdf(Poisson(expected_x), round(Int, img.pixels[h, w]))
    end

    log_prior + log_like
end

end
