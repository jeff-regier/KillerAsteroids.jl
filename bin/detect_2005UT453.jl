#!/usr/bin/env julia

# Detects the asteroids in a syntheic image now and eventually
# the asteroid 2005 UT453 too.

VERSION < v"0.4.0-dev" && using Docile

import FITSIO
import WCSLIB
using Distributions


const r_prior = LogNormal(1500., sqrt(500)) # mode is 1000
const v_prior = MvNormal([3., 5], [4. 1; 3 2]) # pixels / second


@doc """An image, taken though a particular filter band""" ->
type Image
    H::Int64  # The image height.
    W::Int64  # The image width.
    pixels::Matrix{Float64}  # An HxW matrix of pixel intensities.
    epsilon::Float64  # The background noise in nanomaggies
    iota::Float64  # expect number of photons per nanomaggie
    psf::Matrix{Float64}  # the point spread function
end


function wrapped_poisson(rate::Float64)
    0 < rate ? float(rand(Poisson(rate))) : 0.
end


type AsteroidParams
    r::Float64  # brightness in nanomaggies
    u::Vector{Float64}  # position at time 0
    v::Vector{Float64}  # velocity (constant over time)
end


function generate_sample_synthetic_image()
    H, W = 30, 30
    epsilon = 2.2
    iota = 42.0

    pixels = Array(Float64, H, W)
    for w in 1:W, h in 1:H 
        pixels[h, w] = wrapped_poisson(epsilon * iota)
    end

    bvn_for_psf = MvNormal([2.99, 3.01], [0.2 0; 0 0.2])
    psf = Array(Float64, 5, 5)
    for w2 in 1:5, h2 in 1:5
        psf[h2, w2] = pdf(bvn_for_psf, [h2, w2])
    end
    psf /= sum(psf)

    ast = AsteroidParams(1000., [20, 12.], [3.1, 5.1])
    for w2 in 1:5, h2 in 1:5
        h = int(ast.u[1]) + h2 - 3
        w = int(ast.u[2]) + w2 - 3
        pixels[h, w] += wrapped_poisson(ast.r * iota) * psf[h2, w2]
    end

    Image(H, W, pixels, epsilon, iota, psf)
end


@doc """
arguments:
  ast: parameters for a candidate asteroid
  img: an astronomical image
""" ->
function compute_log_likelihood(ast::AsteroidParams, img::Image)
    log_prior = logpdf(r_prior, ast.r) * logpdf(v_prior, ast.v)

    log_like = 0.
    for w2 in 1:5, h2 in 1:5
        h = int(ast.u[1]) + h2 - 3
        w = int(ast.u[2]) + w2 - 3
        expected_x = img.iota * (img.epsilon + ast.r * img.psf[h2, w2])
        log_like += logpdf(Poisson(expected_x), int(img.pixels[h, w]))
    end

    log_prior + log_like
end


const test_img = generate_sample_synthetic_image()
@assert test_img.pixels[20, 12] > 300


const good_ast = AsteroidParams(1000., [20, 12.], [3.1, 5.1])
good_ll = compute_log_likelihood(good_ast, test_img)

const bad_ast = AsteroidParams(1000., [19.4, 12.], [3.1, 5.1])
bad_ll = compute_log_likelihood(bad_ast, test_img)

println("$good_ll > $bad_ll")
@assert good_ll > bad_ll



