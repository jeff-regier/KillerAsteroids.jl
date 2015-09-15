module Synthetic

using Distributions
using ModelInit


function wrapped_poisson(rate::Float64)
    0 < rate ? float(rand(Poisson(rate))) : 0.
end


function generate_sample_image(psf::Matrix{Float64})
    H, W = 30, 30
    epsilon = 2.2
    iota = 42.0

    pixels = Array(Float64, H, W)
    for w in 1:W, h in 1:H
        pixels[h, w] = wrapped_poisson(epsilon * iota)
    end

    ast = AsteroidParams(1000., [20, 12.], [3.1, 5.1])
    for w2 in 1:5, h2 in 1:5
        h = round(Int, ast.u[1]) + h2 - 3
        w = round(Int, ast.u[2]) + w2 - 3
        pixels[h, w] += wrapped_poisson(ast.r * iota) * psf[h2, w2]
    end

    Image(H, W, pixels, epsilon, iota, psf, 0.)
end


function generate_sample_image()
    bvn_for_psf = MvNormal([2.99, 3.01], [0.2 0; 0 0.2])
    psf = Array(Float64, 5, 5)
    for w2 in 1:5, h2 in 1:5
        psf[h2, w2] = pdf(bvn_for_psf, [h2, w2])
    end
    psf /= sum(psf)

    generate_sample_image(psf)
end


end
