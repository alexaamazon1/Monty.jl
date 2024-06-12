using Test
using Monty
using Unitful
using Distributions
using LinearAlgebra: diag
using Meshes: Point, PointSet, PolyArea, Triangle
using GeoStatsFunctions
using Random

using Monty.Checks

const NRAND = 2_000_000

const origin = Point(0.0, 0.0)

@testset "Checks" begin
    @testset "ispositive" begin
        @test !ispositive(-1)
        @test !ispositive(0)
        @test ispositive(1)
    end

    @testset "isnonnegative" begin
        @test !isnonnegative(-1)
        @test isnonnegative(0)
        @test isnonnegative(1)
    end

    @testset "isfractional" begin
        @test !isfractional(-1)
        @test isfractional(0)
        @test isfractional(0.5)
        @test isfractional(1)
        @test !isfractional(2)
    end

    @testset "islessthanone" begin
        @test !islessthanone(2)
        @test !islessthanone(1)
        @test islessthanone(0)
    end
end

@testset "Mixing Model" begin
    @testset "Core Functions" begin
        @testset "Feedstock Fraction" begin
            @test feedstockfraction(Exponential(), 1) ≈ 1 - exp(-1)
            @test feedstockfraction(Uniform(0.5, 1), 2) ≈ 1
            @test feedstockfraction(Beta(2, 2), 0.5u"m") ≈ 0.5
            @test_throws AssertionError feedstockfraction(Normal(), 1)
            @test_throws AssertionError feedstockfraction(Exponential(), -1)
        end

        @testset "Feedstock Mass" begin
            @test feedstockmass(1, 1, 1) == 0
            @test feedstockmass(1, 1, 0) == 1
            @test 0.5 ≈ Monty.feedstockmass(5, 0.2, 0.5)
            @test begin
                0.5u"kg/m^2" ≈
                Monty.feedstockmass(5u"kg/m^2", 0.2u"kg/kg", 0.5u"kg/kg")
            end
            @test_throws AssertionError feedstockmass(-1, 1, 1)
            @test_throws AssertionError feedstockmass(1, 2, 1)
            @test_throws AssertionError feedstockmass(1, 1, 2)
        end

        @testset "Feedstock Element Mass" begin
            @test feedstockelementmass(1, 1, 1, 1) == 0
            @test feedstockelementmass(1, 1, 0, 1) == 1
            @test feedstockelementmass(1, 1, 0, 0) == 0
            @test 5e-5 ≈ feedstockelementmass(5, 0.2, 0.5, 100e-6)
            @test 5e-5u"kg/m^2" ≈ feedstockelementmass(
                5u"kg/m^2",
                0.2u"kg/kg",
                0.5u"kg/kg",
                100e-6u"kg/kg",
            )
            @test_throws AssertionError feedstockelementmass(-1, 1, 1, 0.1)
            @test_throws AssertionError feedstockelementmass(1, -1, 1, 0.1)
            @test_throws AssertionError feedstockelementmass(1, 1, -1, 0.1)
            @test_throws AssertionError feedstockelementmass(1, 1, 1, 2)
        end

        @testset "Feedstock Thickness" begin
            @test feedstockthickness(1, 1, Inf, 0) ≈ 0
            @test feedstockthickness(1, 0.5, 2, 0) ≈ 0.25
            @test feedstockthickness(1, 0.5, 2, 1) ≈ 0
            @test 6e-4 ≈ feedstockthickness(5, 0.3, 2.5e3, 0)
            @test 6e-4u"m" ≈ feedstockthickness(
                5u"kg/m^2",
                0.3u"kg/kg",
                2.5e3u"kg/m^3",
                0u"kg/kg",
            )
            @test_throws AssertionError feedstockthickness(-1, 1, 1, 0)
            @test_throws AssertionError feedstockthickness(1, 2, 1, 0)
            @test_throws AssertionError feedstockthickness(1, 1, 0, 0)
        end

        @testset "Soil Mass" begin
            @test soilmass(1, 0) == 0
            @test 100 ≈ soilmass(1e3, 0.1)
            @test 100u"kg/m^2" ≈ soilmass(1e3u"kg/m^3", 0.1u"m")
            @test_throws AssertionError soilmass(0, 1)
            @test_throws AssertionError soilmass(1, -1)
        end

        @testset "Soil Element Mass" begin
            @test soilelementmass(1, 0, 1) == 0
            @test soilelementmass(1, 1, 0) == 0
            @test 0.015 ≈ soilelementmass(1e3, 0.15, 100e-6)
            @test 0.015u"kg/m^2" ≈
                  soilelementmass(1e3u"kg/m^3", 0.15u"m", 100e-6u"kg/kg")
            @test_throws AssertionError soilelementmass(0, 1, 0.1)
            @test_throws AssertionError soilelementmass(1, -1, 0.1)
            @test_throws AssertionError soilelementmass(1, 1, -0.1)
        end
    end

    @testset "Interface" begin
        core = singlecore((x=0.1, y=0.2), 1.0)
        @test core[1] ≈ 0.1
        @test core[2] ≈ 0.2
        @test core[:x] ≈ 0.1
        @test core[:y] ≈ 0.2
        @test core.mass ≈ 1
        @test nanalyte(core) == 2
        @test analytes(core) == (:x, :y)
        z = zero(core)
        @test z[1] ≈ 0
        @test z[2] ≈ 0
        @test z.mass ≈ 0
    end

    @testset "Model" begin

        # no dissolution, infinite density feedstock
        γ = 1.0
        d = 0.1
        a = 1.0
        Q = 1.0
        ρf = Inf
        cf = 1e-3
        ρs = 1e3
        cs = 1e-4
        𝓁 = 0.0
        ℒ = 0.0
        α = Q / (ρs * d + Q)
        c₁ = cf * α + (1 - α) * cs
        c₂ = mixing(γ, d, a, Q, ρf, cf, ρs, cs, 𝓁, ℒ)[:analyte]
        @test c₁ ≈ c₂

        # equal masses, infinite density feedstock, no dissolution
        Γ = Uniform(0, 0.005)
        d = 0.005
        a = 1.0
        Q = 5.0
        ρf = Inf
        cf = 1e-3
        ρs = 1e3
        cs = 2e-3
        𝓁 = 0.0
        ℒ = 0.0
        c₁ = 0.0015
        c₂ = mixing(Γ, d, a, Q, ρf, cf, ρs, cs, 𝓁, ℒ)[:analyte]
        @test c₁ ≈ c₂

        # equal masses, finite density feedstock, no dissolution
        γ = 1.0
        d = 0.01
        a = 1.0
        Q = 5.0
        ρf = 1e3
        cf = 1e-3
        ρs = 1e3
        cs = 2e-3
        𝓁 = 0.0
        ℒ = 0.0
        c₁ = 0.0015
        c₂ = mixing(γ, d, a, Q, ρf, cf, ρs, cs, 𝓁, ℒ)[:analyte]
        @test c₁ ≈ c₂

        #infinite density feedstock, half sampling fraction, no dissolution
        γ = 0.5
        d = 0.05
        a = 1.0
        Q = 5.0
        ρf = Inf
        cf = 1e-3
        ρs = 1e3
        cs = 1e-4
        𝓁 = 0.0
        ℒ = 0.0
        α = (Q / 2) / (ρs * d + Q / 2)
        c₁ = α * cf + (1 - α) * cs
        c₂ = mixing(γ, d, a, Q, ρf, cf, ρs, cs, 𝓁, ℒ)[:analyte]
        @test c₁ ≈ c₂

        #infinite density feedstock, 40 % dissolution with 40 % mass loss
        γ = 1.0
        d = 0.1
        a = 1.0
        Q = 5.0
        ρf = Inf
        cf = 1e-3
        ρs = 1e3
        cs = 1e-4
        𝓁 = 0.4
        ℒ = 0.4
        α = 0.6Q / (ρs * d + 0.6Q)
        c₁ = α * cf + (1 - α) * cs
        c₂ = mixing(γ, d, a, Q, ρf, cf, ρs, cs, 𝓁, ℒ)[:analyte]
        @test c₁ ≈ c₂

        #infinite density feedstock, 70 % dissolution with no mass loss
        γ = 1.0
        d = 0.1
        a = 1.0
        Q = 5.0
        ρf = Inf
        cf = 1e-3
        ρs = 1e3
        cs = 1e-4
        𝓁 = 0.7
        ℒ = 0.0
        α = Q / (ρs * d + Q)
        c₁ = α * 0.3cf + (1 - α) * cs
        c₂ = mixing(γ, d, a, Q, ρf, cf, ρs, cs, 𝓁, ℒ)[:analyte]
        @test c₁ ≈ c₂

        #arbitrary case to catch unintended changes
        γ = 0.9
        d = 0.144
        a = 0.025
        Q = 3.3
        ρf = 2563.0
        cf = 1230e-6
        ρs = 885.0
        cs = 287e-6
        𝓁 = 0.55
        ℒ = 0.48
        c₁ = 0.00029634715236715744
        c₂ = mixing(γ, d, a, Q, ρf, cf, ρs, cs, 𝓁, ℒ)[:analyte]
        @test c₁ ≈ c₂
    end

    @testset "Compositing" begin
        core₁ = singlecore((x=0.1, y=0.2), 1.0)
        core₂ = singlecore((x=0.2, y=0.1), 2.0)
        core₃ = singlecore((x=0.3, y=0.5), 3.0)

        composite = core₁ + core₂
        @test composite[:x] ≈ (0.1 + 2 * 0.2) / (1 + 2)
        @test composite[:y] ≈ (0.2 + 2 * 0.1) / (1 + 2)
        @test composite.mass ≈ 3
        @test composite.cores == 2

        composite = core₂ + core₃
        @test composite[:x] ≈ (2 * 0.2 + 3 * 0.3) / (2 + 3)
        @test composite[:y] ≈ (2 * 0.1 + 3 * 0.5) / (2 + 3)
        @test composite.mass ≈ 5
        @test composite.cores == 2

        composite = core₁ + core₂ + core₃
        @test composite[:x] ≈ (0.1 + 2 * 0.2 + 3 * 0.3) / (1 + 2 + 3)
        @test composite[:y] ≈ (0.2 + 2 * 0.1 + 3 * 0.5) / (1 + 2 + 3)
        @test composite.mass ≈ 6
        @test composite.cores == 3
    end
end

@testset "Measurement" begin
    @testset "Noise" begin
        rng = Xoshiro(1)
        x = 2.0
        σ = 0.01
        y = map(_ -> noise(rng, x, σ), 1:NRAND)
        @test isapprox(y |> std, σ * x, rtol=2)
    end

    @testset "Measurement" begin
        rng = Xoshiro(1)

        core = singlecore((x=1000.0, y=200.0), 100.0)
        M = map(_ -> measure(rng, core, (x=0.01, y=0.01), 0.01), 1:NRAND)
        @test isapprox(std(getindex.(M, :x)), 10.0, rtol=1e-2)
        @test isapprox(std(getindex.(M, :y)), 2.0, rtol=1e-2)
        @test isapprox(std(getfield.(M, :mass)), 1.0, rtol=1e-2)
    end
end

@testset "Leaching" begin

    # an arbitrary leached fraction asymptote
    C = 0.89
    # some times to test monotonicity
    t = LinRange(0.0, 10.0, 101)

    𝓌 = ExponentialLeaching(C=C)
    @test 𝓌(0.0) ≈ 0
    @test 𝓌(-0.1) ≈ 1
    @test 𝓌(2.0) ≈ C - C * exp(-2)
    @test 𝓌(1e9) ≈ C
    @test all(diff(𝓌.(t)) .>= 0.0)

    𝓌 = MultiExponentialLeaching(λ=(0.1, 1.0, 10.0), C=C)
    @test 𝓌(0.0) ≈ 0
    @test 𝓌(-0.1) ≈ 1
    @test 𝓌(2.0) ≈ C - C * (exp(-0.2) + exp(-2) + exp(-20)) / 3
    @test 𝓌(1e9) ≈ C
    @test all(diff(𝓌.(t)) .>= 0.0)

    𝓌 = SeasonalLeaching(λ=2, C=C, floor=0.01, power=2, phase=π / 2)
    @test 𝓌(0.0) ≈ 0
    @test 𝓌(-0.1) ≈ 1
    @test 𝓌(1e9) ≈ C
    @test all(diff(𝓌.(t)) .>= 0.0)
end

@testset "Jitters" begin
    rng = Xoshiro(1)
    points = fill(Point(1.0, 2.0), NRAND)
    jitterpoints!(rng, points, 0.5)
    @test isapprox(std(map(p -> p.coords[1], points)), 0.5, rtol=1e-2)
    @test isapprox(std(map(p -> p.coords[2], points)), 0.5, rtol=1e-2)

    point = Point(1.0, 2.0)
    J = NoJitter()
    @test J(point) == point

    J = Jitter(1.0, seed=1)
    points = map(J, fill(point, NRAND))
    @test isapprox(std(map(p -> p.coords[1], points)), 1.0, rtol=1e-2)
    @test isapprox(std(map(p -> p.coords[2], points)), 1.0, rtol=1e-2)

    J = Jitter(Uniform(-1, 1))
    points = map(J, fill(point, NRAND))
    @test isapprox(
        std(map(p -> p.coords[1], points)),
        0.5773502691896257,
        rtol=1e-2,
    )
    @test isapprox(
        std(map(p -> p.coords[2], points)),
        0.5773502691896257,
        rtol=1e-2,
    )

    @test_throws AssertionError Jitter(Uniform(0, 1))

    J = GridCentroidJitter(1.0, 1.0)
    points = map(J, fill(origin, NRAND))
    x = map(p -> p.coords[1], points)
    y = map(p -> p.coords[2], points)
    @test all(-0.5 .<= x .<= 0.5)
    @test isapprox(x |> mean, 0.0, atol=1e-2)
    @test isapprox(x |> std, 0.28867513459481287, rtol=1e-2)
    @test all(-0.5 .<= y .<= 0.5)
    @test isapprox(y |> mean, 0.0, atol=1e-2)
    @test isapprox(y |> std, 0.28867513459481287, rtol=1e-2)
end

@testset "Stencils" begin
    @testset "Single Core" begin
        S = SingleCoreStencil(Float64)
        @test ncore(S) == 1.0
        @test S(Point(1.0, 1.0)) == Point(1.0, 1.0)
    end

    @testset "Random" begin
        S = RandomStencil(5, 1.0, seed=1)
        @test ncore(S) == 5
        points = map(S, fill(Point(0.0, 0.0), NRAND))
        @test isapprox(map(p -> p.coords[1], points) |> std, 1.0, rtol=1e-2)
        @test isapprox(map(p -> p.coords[2], points) |> std, 1.0, rtol=1e-2)
    end

    @testset "Circle" begin
        S = CircleStencil(4, 1.0)
        @test ncore(S) == 4
        @test_throws AssertionError S(origin, 0)
        @test isapprox(S(origin, 1), Point(1.0, 0.0), atol=1e-6)
        @test isapprox(S(origin, 2), Point(0.0, 1.0), atol=1e-6)
        @test isapprox(S(origin, 3), Point(-1.0, 0.0), atol=1e-6)
        @test isapprox(S(origin, 4), Point(0.0, -1.0), atol=1e-6)
        @test_throws AssertionError S(origin, 5)
    end

    @testset "HubSpoke" begin
        S = HubSpokeStencil(5, 1.0)
        @test ncore(S) == 5
        @test_throws AssertionError S(origin, 0)
        @test isapprox(S(origin, 1), Point(0.0, 0.0), atol=1e-6)
        @test isapprox(S(origin, 2), Point(1.0, 0.0), atol=1e-6)
        @test isapprox(S(origin, 3), Point(0.0, 1.0), atol=1e-6)
        @test isapprox(S(origin, 4), Point(-1.0, 0.0), atol=1e-6)
        @test isapprox(S(origin, 5), Point(0.0, -1.0), atol=1e-6)
        @test_throws AssertionError S(origin, 6)
    end

    @testset "Line" begin
        S = LineStencil(3, 1.0, π / 4)
        @test ncore(S) == 3
        @test_throws AssertionError S(origin, 0)
        @test isapprox(S(origin, 1), Point(-1 / √2, -1 / √2))
        @test isapprox(S(origin, 2), origin, atol=1e-6)
        @test isapprox(S(origin, 3), Point(1 / √2, 1 / √2))
        @test_throws AssertionError S(origin, 4)
    end
end

@testset "Sampling" begin
    treatment = PolyArea([(0, 60), (0, 160), (50, 160), (50, 50)])
    control = Triangle((80, 60), (60, 90), (90, 100))
    times = Float64[-0.03, 0.01, 1]
    rng = Xoshiro(1)
    N = 10
    points = Monty.Sampling.randin(rng, treatment, N)

    @test all(points .∈ treatment)
    @test all(gridpointoverlay(rng, treatment, (10, 10)) .∈ treatment)

    @testset "Paired Plans" begin
        plan = pairedsampleplan(points, times)
        @test nsample(plan) == N * length(times)
        plan = pairedsampleplan(points, control, times)
        @test nsample(plan) == N * length(times)
        plan = pairedsampleplan(
            Monty.Sampling.randin(rng, treatment, N),
            Monty.Sampling.randin(rng, control, N),
            times,
        )
        @test nsample(plan) == 2N * length(times)
        plan = pairedsampleplan(rng, treatment, N, times)
        @test nsample(plan) == N * length(times)
        plan = pairedsampleplan(rng, treatment, N, control, N, times)
        @test nsample(plan) == 2N * length(times)
    end

    @testset "Independent Plans" begin
        plan = randomsampleplan(
            rng,
            treatment,
            [2, 3, 4],
            control,
            [1, 2, 3],
            times,
        )
        @test nsample(plan) == 9 + 6
        plan = randomsampleplan(rng, treatment, 8, control, 5, times)
        @test nsample(plan) == length(times) * (8 + 5)
        plan = randomsampleplan(rng, treatment, N, times)
        @test nsample(plan) == length(times) * N
    end
end

@testset "Gaussian" begin
    @testset "MvNormal" begin
        μ = randn(10)
        σ = GaussianCovariance(nugget=0.1, sill=4.0, range=1.0)
        points = randn(2, 10) |> PointSet |> collect
        X = MvNormal(μ, σ, points)
        rng = Xoshiro(0)
        r = rand(rng, X, NRAND)
        @test all(isapprox.(mean(r, dims=2), μ, atol=1e-2))
        @test all(isapprox.(std(r, dims=2), 2.0, atol=1e-2))
    end

    @testset "Simulator" begin
        rng = Xoshiro(1)
        N = 10
        points = randn(rng, 2, N) |> PointSet |> collect
        nugget = abs(randn(rng))
        sill = abs(randn(rng))
        sill, nugget = max(sill, nugget), min(sill, nugget)
        range = abs(randn(rng))
        c = ExponentialCovariance(nugget=nugget^2, sill=sill^2, range=range)
        μ = randn(rng)
        gs = GaussianSimulator(points, μ, c)
        y = Matrix{Float64}(undef, N, NRAND)
        for j ∈ 1:NRAND
            rand!(rng, gs, view(y, :, j))
        end
        @test all(isapprox.(mean(y, dims=2), μ, atol=2e-2))
        @test all(isapprox.(std(y, dims=2), sill, atol=2e-2))
        @test all(isapprox.(std(y, dims=2), sill, atol=2e-2))
    end

    @testset "Cosimulator" begin
        rng = Xoshiro(2)
        N = 10
        points = randn(rng, 2, N) |> PointSet |> collect
        for ρ ∈ Float64[0.0, 0.1, 0.5, 0.9, 0.99, 0.999]
            nugget = abs(randn(rng))
            sill = abs(randn(rng))
            sill, nugget = max(sill, nugget), min(sill, nugget)
            range = abs(randn(rng))
            c = SphericalCovariance(nugget=nugget^2, sill=sill^2, range=range)
            μ = (a=randn(rng), b=randn(rng))
            gc = GaussianCosimulator(points, μ, c, ρ)
            y = (
                a=Matrix{Float64}(undef, N, NRAND),
                b=Matrix{Float64}(undef, N, NRAND),
            )
            for j ∈ 1:NRAND
                rand!(rng, gc, (a=view(y[:a], :, j), b=view(y[:b], :, j)))
            end
            @test all(isapprox.(mean(y[:a], dims=2), μ[:a], atol=2e-2))
            @test all(isapprox.(mean(y[:b], dims=2), μ[:b], atol=2e-2))
            @test all(isapprox.(std(y[:a], dims=2), sill, atol=2e-2))
            @test all(isapprox.(std(y[:b], dims=2), sill, atol=2e-2))
            @test all(isapprox.(diag(cor(y[:a], y[:b], dims=2)), ρ, atol=2e-2))
        end
    end
end

@testset "Simulation" begin
    @testset "Moisture" begin
        @test moisturefraction(0.5) ≈ 1 / 3
        @test moistureratio(0.5) ≈ 1
        @test 0.5 |> moisturefraction |> moistureratio ≈ 0.5
        @test 0.5 |> moistureratio |> moisturefraction ≈ 0.5
    end
end

@testset "CDRPotential" begin
    @testset "cation to CO2" begin
        @test cationCO2(:Ca) ≈ 2.1962173761165724
        @test cationCO2(:Mg) ≈ 3.621477062332853
        @test cationCO2(:Na) ≈ 1.914329781390481
        @test cationCO2(:K) ≈ 1.1256243877611045
        @test cationCO2("Ca") ≈ 2.1962173761165724
        @test cationCO2("Mg") ≈ 3.621477062332853
        @test cationCO2("Na") ≈ 1.914329781390481
        @test cationCO2("K") ≈ 1.1256243877611045
        @test_throws AssertionError cationCO2("Fe")
        @test_throws AssertionError cationCO2(:Ti)
    end
    @testset "CDR potential" begin
        @test cdrpotential(:Ca, 0.1) ≈ 0.21962173761165724
        @test cdrpotential(:Mg, 0.1) ≈ 0.3621477062332853
        @test cdrpotential(:Na, 0.1) ≈ 0.1914329781390481
        @test cdrpotential(:K, 0.1) ≈ 0.11256243877611045
        @test cdrpotential(
            Dict(:Ca => 0.07, :Mg => 0.05, :Na => 0.01, :K => 0.01),
        ) ≈ 0.36520861113631864
        @test cdrpotential((Ca=0.07, Mg=0.05, Na=0.01, K=0.01)) ≈
              0.36520861113631864
    end
end
