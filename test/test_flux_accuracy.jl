using KillerAsteroids

include("sample_data.jl")

function flux_accuracy_all_synthetic_data()

    fac_min = 0.75 # 0.75 times as bright
    fac_max = 1.25 # 1.25 times as bright
    nfac = 51
    fac_vals = linspace(fac_min, fac_max, nfac)
    lp_vals = zeros(nfac)

    for i=1:nfac
         test_ast = AsteroidParams(fac_vals[i]*sample_ast.r, sample_ast.u, sample_ast.v)
         test_lp = compute_log_probability([test_ast,], [sample_img,], sample_prior)
         lp_vals[i] = test_lp
    end

    lp_max, indbest = findmax(lp_vals)
    return fac_vals[indbest]
end
