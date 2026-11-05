# Independent checks of ROOT's fitted parameters, covariance, and chi2 against
# closed-form weighted least squares (pol1) and the analytic Jacobian (gaus).
# The reference numbers are computed here, not taken from ROOT.

"""Inverse of a small nonsingular matrix by Gauss-Jordan elimination with pivoting (test only)."""
function small_inverse(A::AbstractMatrix)
    n = size(A, 1)
    M = hcat(Float64.(A), [Float64(i == j) for i in 1:n, j in 1:n])
    for c in 1:n
        p = c - 1 + argmax(abs.(M[c:n, c]))
        M[[c, p], :] = M[[p, c], :]
        M[c, :] ./= M[c, c]
        for r in 1:n
            r == c || (M[r, :] .-= M[r, c] .* M[c, :])
        end
    end
    return M[:, n+1:end]
end

"""Weighted linear least squares: parameters, covariance inv(X'WX), chi2, residual sum of squares."""
function weighted_least_squares(X, y, sigma)
    w = 1 ./ sigma .^ 2
    XtWX = [sum(w .* X[:, i] .* X[:, j]) for i in axes(X, 2), j in axes(X, 2)]
    cov = small_inverse(XtWX)
    p = cov * [sum(w .* X[:, i] .* y) for i in axes(X, 2)]
    r = y .- X * p
    return (parameters=p, covariance=cov, chi2=sum(w .* r .^ 2), rss=sum(r .^ 2))
end

"""Jacobian of A exp(-(x-mu)^2/(2 s^2)) with respect to (A, mu, s)."""
gaus_jacobian(x, (A, mu, s)) = hcat(exp.(-(x .- mu) .^ 2 ./ (2s^2)),
                                    A .* exp.(-(x .- mu) .^ 2 ./ (2s^2)) .* (x .- mu) ./ s^2,
                                    A .* exp.(-(x .- mu) .^ 2 ./ (2s^2)) .* (x .- mu) .^ 2 ./ s^3)

@testset "fit errors: closed-form references (no ROOT)" begin
    x, y = collect(0.0:4.0), [1.1, 2.9, 5.2, 6.8, 9.1]
    X = hcat(ones(5), x)
    ref = weighted_least_squares(X, y, fill(0.2, 5))
    @test ref.parameters ≈ [1.04, 1.99] && ref.chi2 ≈ 2.675 && ref.rss ≈ 0.107
    @test ref.covariance ≈ [0.024 -0.008; -0.008 0.004]
    @test small_inverse([5.0 10; 10 30]) ≈ [0.6 -0.2; -0.2 0.1]
end

@testset "fit errors: ROOT pol1 and gaus covariance" begin
    if TEST_NATIVE_ROOT
        x, y = collect(0.0:4.0), [1.1, 2.9, 5.2, 6.8, 9.1]
        X = hcat(ones(5), x)
        # Tolerances: a pol1 chi2 is quadratic, so Migrad's minimum and its Hessian are exact up
        # to the convergence tolerance (EDM <= 1e-3 relative to the chi2 scale): 1e-6 relative.
        tol = 1e-6
        # 1. Known errors (sigma = 0.2): covariance inv(X'WX), errors not normalized.
        weighted = fit_root(GraphSeries(x, y; yerrors=fill(0.2, 5)), :pol1; initial=[0.0, 1.0])
        @test isapprox(weighted.parameters, [1.04, 1.99]; rtol=tol)
        @test isapprox(weighted.covariance, [0.024 -0.008; -0.008 0.004]; rtol=tol, atol=1e-12)
        @test isapprox(weighted.errors, sqrt.([0.024, 0.004]); rtol=tol)
        @test isapprox(weighted.chi2, 2.675; rtol=tol) && weighted.ndf == 3
        @test !weighted.normalized_errors && weighted.valid
        # 2. No errors (plain TGraph): ROOT fits with unit errors and scales the covariance by
        # chi2/ndf = RSS/3, reporting normalized errors.
        plain = fit_root(GraphSeries(x, y), :pol1; initial=[0.0, 1.0])
        @test isapprox(plain.parameters, [1.04, 1.99]; rtol=tol)
        @test isapprox(plain.covariance, (0.107 / 3) .* small_inverse(X' * X); rtol=tol, atol=1e-12)
        @test isapprox(plain.chi2, 0.107; rtol=tol) && plain.ndf == 3
        @test plain.normalized_errors
        @test_throws ArgumentError chi2_pvalue(plain)       # RSS with unit errors is no goodness of fit
        @test 0.4 < chi2_pvalue(weighted) < 0.5               # P(χ²₃ > 2.675) ≈ 0.444
        # 3. Nonuniform errors.
        sigma = [0.1, 0.3, 0.2, 0.5, 0.15]
        ref = weighted_least_squares(X, y, sigma)
        nonuniform = fit_root(GraphSeries(x, y; yerrors=sigma), :pol1; initial=[0.0, 1.0])
        @test isapprox(nonuniform.parameters, ref.parameters; rtol=tol)
        @test isapprox(nonuniform.covariance, ref.covariance; rtol=tol, atol=1e-12)
        @test isapprox(nonuniform.chi2, ref.chi2; rtol=tol) && !nonuniform.normalized_errors
        # 4. Hist1D: bin centres, errors sqrt(sumw2); bins with zero error are excluded by ROOT's
        # chi2 fit (an empty bin and a nonzero bin without variance).
        h = Hist1D(-0.5:1.0:6.5, [y; 0.0; 100.0], [sigma .^ 2; 0.0; 0.0])
        hfit = fit_root(h, :pol1; initial=[0.0, 1.0])
        @test isapprox(hfit.parameters, ref.parameters; rtol=tol)
        @test isapprox(hfit.covariance, ref.covariance; rtol=tol, atol=1e-12)
        @test isapprox(hfit.chi2, ref.chi2; rtol=tol) && hfit.ndf == 3 && hfit.npoints == 5
        # 5. Graph x errors: :ignore (ROOT option EX0) is the y-only fit; :effective uses
        # ROOT's effective variance sigma_y^2 + (f'(x) sigma_x)^2, here a much smaller chi2.
        withx = GraphSeries(x, y; xerrors=fill(0.5, 5), yerrors=fill(0.2, 5))
        ignored = fit_root(withx, :pol1; initial=[0.0, 1.0], xerrors=:ignore)
        @test isapprox(ignored.parameters, weighted.parameters; rtol=tol)
        @test isapprox(ignored.covariance, weighted.covariance; rtol=tol, atol=1e-12)
        @test isapprox(ignored.chi2, 2.675; rtol=tol)
        effective = fit_root(withx, :pol1; initial=[0.0, 1.0])
        slope = effective.parameters[2]
        eff_sigma = sqrt.(0.2^2 .+ (slope * 0.5)^2)
        @test isapprox(effective.chi2, sum((y .- X * effective.parameters) .^ 2 ./ eff_sigma .^ 2); rtol=1e-6)
        @test effective.chi2 < 0.1 * ignored.chi2
        # 6. Noiseless Gaussian (A = 10, mu = 0.3, sigma = 1.2), errors 0.1. At a zero-residual
        # minimum the chi2 Hessian is exactly 2 J'WJ, so the covariance is inv(J'WJ). Minuit
        # estimates it by finite differences (HESSE); 1e-3 relative allows for that step error
        # while still failing for a wrong error definition (factor 2 or 4) or normalization.
        truth = [10.0, 0.3, 1.2]
        xs = collect(-3.0:0.25:3.0)
        model(x) = truth[1] * exp(-(x - truth[2])^2 / (2 * truth[3]^2))
        J = gaus_jacobian(xs, truth)
        expected = small_inverse((J' * J) ./ 0.1^2)
        g = fit_root(GraphSeries(xs, model.(xs); yerrors=fill(0.1, length(xs))), :gaus; initial=[8.0, 0.0, 1.0])
        hg = fit_root(Hist1D(-3.125:0.25:3.125, model.(xs), fill(0.01, length(xs))), :gaus; initial=[8.0, 0.0, 1.0])
        for fit in (g, hg)
            @test isapprox(fit.parameters, truth; rtol=1e-5)
            @test fit.chi2 < 1e-6 && fit.ndf == length(xs) - 3 && !fit.normalized_errors
            @test isapprox(fit.covariance, expected; rtol=1e-3, atol=1e-3 * maximum(abs, expected))
            @test isapprox(fit.errors, sqrt.([expected[i, i] for i in 1:3]); rtol=1e-3)
        end
    else
        @test_skip false
    end
end
