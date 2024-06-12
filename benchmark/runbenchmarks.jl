using BenchmarkTools
using Monty
using Unitful
using Meshes
using GeoStatsFunctions
using Distributions
using Random: Xoshiro, rand!

##

const core = Monty.singlecore((Ca=1e4u"ppm", Mg=1e3u"ppm", T=50.0u"ppm"), 1u"g")

##

println("setting up benchmarks...")
suite = BenchmarkGroup(["Monty"])

## mixing model evaluation

suite["mixing model"]["sample"]["single analyte"] = @benchmarkable begin
    mixing(Exponential(0.025), 0.01, 0.005, 5.0, 3e3, 0.05, 1e3, 0.001, 0.4, 0.3)
end

suite["mixing model"]["sample"]["no units"] = @benchmarkable begin
    mixing(
        Exponential(0.025),
        0.01,
        0.005,
        5.0,
        3e3,
        (Ca=0.05, Mg=0.04, T=1e-4),
        1e3,
        (Ca=0.001, Mg=0.0005, T=1e-5),
        (Ca=0.4, Mg=0.6, T=0.0),
        0.3,
    )
end

suite["mixing model"]["sample"]["units"] = @benchmarkable begin
    mixing(
        0.5u"kg/kg",
        0.01u"m",
        0.005u"m^2",
        5.0u"kg/m^2",
        3e3u"kg/m^3",
        (Ca=0.05u"kg/kg", Mg=0.04u"kg/kg", T=1e-4u"kg/kg"),
        1e3u"kg/m^3",
        (Ca=0.001u"kg/kg", Mg=0.0005u"kg/kg", T=1e-5u"kg/kg"),
        (Ca=0.4u"kg/kg", Mg=0.6u"kg/kg", T=0.0u"kg/kg"),
        0.3u"kg/kg",
    )
end

## Sample struct operations

suite["mixing model"]["Sample"]["+"] = @benchmarkable $core + $core

suite["mixing model"]["Sample"]["sum"] = @benchmarkable sum($[core, core, core])

## noise application and measurement

suite["measurement"] = @benchmarkable begin
    Monty.measure(rng, c, (Ca=0.03, Mg=0.04), 0.0)
end setup = (c = singlecore((Ca=0.07, Mg=0.05), 0.1); rng = Xoshiro())

## evaluation of leaching models over time

suite["leaching"]["exponential"] = @benchmarkable begin
    𝓌(1.0)
end setup = (𝓌 = ExponentialLeaching(C=0.05))

suite["leaching"]["multi-exponential"] = @benchmarkable begin
    𝓌(2.0)
end setup = (𝓌 = MultiExponentialLeaching(λ=(0.1, 1.0, 10.0), C=0.1))

suite["leaching"]["seasonal"] = @benchmarkable begin
    𝓌(3.0)
end setup = (𝓌 = SeasonalLeaching(C=0.05, floor=0.01, power=2, phase=π / 4))

## jitter and stencil

suite["sampling"]["jitterpoints!"] = @benchmarkable begin
    jitterpoints!(rng, points, 1.0)
end setup = (points = map(x -> Point(x...), randn(2, 2000) |> eachcol);
rng = Xoshiro())

suite["sampling"]["Jitter"]["normal"] = @benchmarkable begin
    J(point)
end setup = (J = Jitter(1.0, seed=1); point = Point(randn(2)...))

suite["sampling"]["Jitter"]["NoJitter"] = @benchmarkable begin
    J(point)
end setup = (J = NoJitter(Float32); point = Point(randn(Float32, 2)...))

suite["sampling"]["stencils"]["random"] = @benchmarkable begin
    S(point, 2)
end setup = (S = RandomStencil(3, 1.0, seed=1); point = Point(randn(2)...))

suite["sampling"]["stencils"]["circular"] = @benchmarkable begin
    S(point, 2)
end setup = (S = CircleStencil(10, 1.0); point = Point(randn(2)...))

suite["sampling"]["stencils"]["hubspoke"] = @benchmarkable begin
    S(point, 2)
end setup = (S = HubSpokeStencil(10, 1.0); point = Point(randn(2)...))

suite["sampling"]["stencils"]["line"] = @benchmarkable begin
    S(point, 2)
end setup = (S = LineStencil(10, 1.0, π / 4); point = Point(randn(2)...))

## sample plans

suite["sampling"]["SamplePlan"]["paired"]["points"]["controlled"] =
    @benchmarkable begin
        pairedsampleplan(points, control, times)
    end setup = (
        points =
            CartesianGrid((0.0, 0.0), (100.0, 100.0), dims=(5, 5)) .|>
            centroid |>
            collect;
        control = Box((0.0, 0.0), (20.0, 20.0));
        times = Float64[0, 1]
    )

suite["sampling"]["SamplePlan"]["paired"]["points"]["uncontrolled"] =
    @benchmarkable begin
        pairedsampleplan(points, times)
    end setup = (
        points =
            CartesianGrid((0.0, 0.0), (100.0, 100.0), dims=(5, 5)) .|>
            centroid |>
            collect;
        times = Float64[0, 1]
    )

suite["sampling"]["SamplePlan"]["paired"]["geometry"]["controlled"] =
    @benchmarkable begin
        pairedsampleplan(rng, treatment, 25, control, 10, times)
    end setup = (
        treatment = PolyArea([
            (0, 0),
            (0, 100),
            (50, 100),
            (50, 50),
            (100, 50),
            (100, 0),
        ]);
        control = Ball((100, 100), 25);
        times = Float64[0, 1];
        rng = Xoshiro(1)
    )

suite["sampling"]["SamplePlan"]["paired"]["geometry"]["uncontrolled"] =
    @benchmarkable begin
        pairedsampleplan(rng, field, 25, times)
    end setup = (
        field = PolyArea([
            (0, 0),
            (0, 100),
            (50, 100),
            (50, 50),
            (100, 50),
            (100, 0),
        ]);
        times = Float64[0, 1];
        rng = Xoshiro(1)
    )

suite["sampling"]["SamplePlan"]["random"]["controlled"] = @benchmarkable begin
    randomsampleplan(rng, treatment, 25, control, 10, times)
end setup = (
    treatment = PolyArea([
        (0, 0),
        (0, 100),
        (50, 100),
        (50, 50),
        (100, 50),
        (100, 0),
    ]);
    control = Ball((100, 100), 25);
    times = Float64[0, 1];
    rng = Xoshiro(1)
)

suite["sampling"]["SamplePlan"]["random"]["uncontrolled"] = @benchmarkable begin
    randomsampleplan(rng, field, 25, times)
end setup = (
    field = PolyArea([
        (0, 0),
        (0, 100),
        (50, 100),
        (50, 50),
        (100, 50),
        (100, 0),
    ]);
    times = Float64[0, 1];
    rng = Xoshiro(1)
)

suite["simulation"]["simulator"] = @benchmarkable begin
    rand!(rng, gc, y)
end setup = (
    rng = Xoshiro();
    gc = GaussianSimulator(
        randn(2, 10) |> PointSet |> collect,
        0.0,
        SphericalCovariance(nugget=2.0, sill=5.0, range=10.0),
    );
    y = zeros(10)
)

suite["simulation"]["cosimulator"] = @benchmarkable begin
    rand!(rng, gc, y)
end setup = (
    rng = Xoshiro();
    gc = GaussianCosimulator(
        randn(2, 10) |> PointSet |> collect,
        (a=1.0, b=2.0),
        SphericalCovariance(nugget=2.0, sill=5.0, range=10.0),
        0.5,
    );
    y = (a=zeros(10), b=zeros(10))
)

suite["simulation"]["no spatial structure"] = @benchmarkable begin
    executeplan!(
        samp,
        plan,
        stencil=stencil,
        samplerjitter=samplerjitter,
        corejitter=corejitter,
    )
    spreading!(rng, sim, plan, Q)
    unmixed!(rng, sim, depth=depth)
    sim.ρf .= 3e3
    rand!(rng, ρs, sim.ρs)
    feedstockconcentration!(rng, sim, (Ca=0.07, Mg=0.05), 0.05)
    rand!(rng, cs[:Ca], sim.cs[:Ca])
    rand!(rng, cs[:Mg], sim.cs[:Mg])
    clamp!(sim.cs[:Ca], 0.0, Inf)
    clamp!(sim.cs[:Mg], 0.0, Inf)
    leaching!(sim, :Ca, 𝓁Ca, plan)
    leaching!(sim, :Mg, 𝓁Mg, plan)
    massloss!(sim, plan, x -> (x[:Ca] / 2 + x[:Mg] / 2))
    analyze!(rng, sim, (Ca=0.01, Mg=0.005), 0.005)
end setup = (
    rng = Xoshiro();
    plan = pairedsampleplan(
        rand(2, 50) |> PointSet |> collect,
        rand(2, 50) |> PointSet |> collect,
        [0.0, 1.0],
    );
    samp = CoreSet(Float64, 200, 5);
    sim = Simulation((:Ca, :Mg), 200, 5);
    stencil = CircleStencil(5, 1.5);
    samplerjitter = Jitter(2.0);
    corejitter = Jitter(0.1);
    Q = truncated(Normal(4.0, 1.0), 0.0, Inf);
    depth = TriangularDist(0.05, 0.15);
    ρs = Normal(1e3, 50);
    cs = (Ca=Normal(0.002, 0.0002), Mg=Normal(0.001, 0.0001));
    𝓁Ca = ExponentialLeaching(λ=1.0);
    𝓁Mg = ExponentialLeaching(λ=1.5)
)

suite["simulation"]["spatial structure"] = @benchmarkable begin
    executeplan!(
        samp,
        plan,
        stencil=stencil,
        samplerjitter=samplerjitter,
        corejitter=corejitter,
    )
    spreading!(rng, sim, plan, Q)
    unmixed!(rng, sim, depth=depth)
    sim.ρf .= 3e3
    rand!(rng, ρs, sim.ρs)
    feedstockconcentration!(rng, sim, (Ca=0.07, Mg=0.05), 0.05)
    updategaussian!(cs, samp.points, cs_cov)
    soilconcentration!(rng, sim, cs)
    clamp!(sim.cs[:Ca], 0.0, Inf)
    clamp!(sim.cs[:Mg], 0.0, Inf)
    leaching!(sim, :Ca, 𝓁Ca, plan)
    leaching!(sim, :Mg, 𝓁Mg, plan)
    massloss!(sim, plan, x -> (x[:Ca] / 2 + x[:Mg] / 2))
    analyze!(rng, sim, (Ca=0.01, Mg=0.005), 0.005)
end setup = (
    rng = Xoshiro();
    plan = pairedsampleplan(
        rand(2, 50) |> PointSet |> collect,
        rand(2, 50) |> PointSet |> collect,
        [0.0, 1.0],
    );
    samp = CoreSet(Float64, 200, 5);
    sim = Simulation((:Ca, :Mg), 200, 5);
    stencil = CircleStencil(5, 1.5);
    samplerjitter = Jitter(2.0);
    corejitter = Jitter(0.1);
    Q = truncated(Normal(4.0, 1.0), 0.0, Inf);
    depth = TriangularDist(0.05, 0.15);
    ρs = Normal(1e3, 50);
    cs = GaussianCosimulator(
        1_000,
        (Ca=fill(0.002, 1_000), Mg=fill(0.001, 1_000)),
        0.5,
    );
    cs_cov = SphericalCovariance(nugget=5e-5^2, sill=1e-4^2, range=20.0);
    𝓁Ca = ExponentialLeaching(λ=1.0);
    𝓁Mg = ExponentialLeaching(λ=1.5)
)

##

println("tuning benchmarks...")
tune!(suite)

##

println("running benchmarks...")
results = run(suite)

##

show(results)

BenchmarkTools.save("benchmarks.json", results)

##
