import express from 'express'
import cors from 'cors'
import { spawn } from 'child_process'
import { fileURLToPath } from 'url'
import { dirname, join } from 'path'
import fs from 'fs/promises'

const __filename = fileURLToPath(import.meta.url)
const __dirname = dirname(__filename)

const app = express()
const PORT = process.env.PORT || 3001

app.use(cors())
app.use(express.json())

// Store simulation results temporarily
let simulationCache = new Map()

// Check if Julia is available
async function checkJulia() {
  return new Promise((resolve) => {
    const julia = spawn('julia', ['--version'])
    julia.on('close', (code) => {
      resolve(code === 0)
    })
    julia.on('error', () => {
      resolve(false)
    })
  })
}

// Generate Julia simulation script
function generateJuliaScript(params) {
  // Parse time points
  const timePoints = params.sampling?.timePoints 
    ? params.sampling.timePoints.split(',').map(t => parseFloat(t.trim())).filter(t => !isNaN(t))
    : [-0.003, 0.0, 1.0];

  return `
using DrWatson
@quickactivate "cases"
using Monty
using Random
using Distributions
using Meshes
using GeoStatsFunctions
using GeoStatsProcesses
using JSON

# Set random seed for reproducibility
const rng = Xoshiro(${params.execution?.randomSeed || 1})

# Grid configuration
grid = CartesianGrid((${params.grid?.xSize || 8}, ${params.grid?.ySize || 8}), (0.0, 0.0), (${params.grid?.cellWidth || 10.0}, ${params.grid?.cellHeight || 10.0}))

# Split treatment and control cells, alternating
treatment = reduce(vcat, [grid[i, :] for i ∈ 1:2:${params.grid?.xSize || 8}])
control = reduce(vcat, [grid[i, :] for i ∈ 2:2:${params.grid?.xSize || 8}])

# Take cell centroids for each group
treatment_points = treatment .|> centroid
control_points = control .|> centroid

# Sampling plan with specified time points
time_points = [${timePoints.join(', ')}]
plan = pairedsampleplan(grid .|> centroid, control, time_points)

# Core sampling configuration
samp = CoreSet(plan, ${params.sampling?.coresPerSample || 5})

# Allocate space for the simulation
sim = Simulation((:Ca, :Mg), samp)

# Feedstock density
sim.ρf .= ${params.feedstock?.density || 1000}

# Spreading patterns with spatial correlation
Q_μ = ${params.feedstock?.applicationRate || 3.5}
Q = GaussianSimulator(samp, Q_μ)
Q_cov = GaussianCovariance(
    MetricBall((${params.covariance?.appRangeX || 5.0}, ${params.covariance?.appRangeY || 500.0})),
    nugget=(3 / 10)^2,
    sill=(3 / 6)^2,
)

# Cross-correlated baseline soil concentrations
cs = GaussianCosimulator(samp, (Ca=${params.soil?.caConcentration || 0.002}, Mg=${params.soil?.mgConcentration || 0.001}), ${params.soil?.crossCorrelation || 0.75})
cs_cov = SphericalCovariance(nugget=(0.001 / 10)^2, sill=(0.001 / 6)^2, range=${params.covariance?.soilRange || 20.0})

# Soil spatial structure
ρs = GaussianSimulator(samp, ${params.soil?.density || 1000.0})
ρs_cov = SphericalCovariance(nugget=20^2, sill=100^2, range=${params.covariance?.densityRange || 30.0})

# Stencil configuration
${getStencilCode(params.sampling?.stencilType, params.sampling?.stencilSize, params.sampling?.coresPerSample)}

# Jitter configuration
samplerjitter = Jitter(${params.jitter?.sampler || 0.75})
corejitter = Jitter(${params.jitter?.core || 0.1})
plannedjitter = GridCentroidJitter(${params.grid?.cellWidth || 10.0}, ${params.jitter?.planned || 5.0}, seed=1)

# Leaching models for each element
Ca_leaching = ExponentialLeaching(λ=${params.leaching?.caRate || 0.4})
Mg_leaching = ExponentialLeaching(λ=${params.leaching?.mgRate || 0.8})

# Number of simulations
nsim = ${params.execution?.numRealizations || 1000}

# Run simulation stack
sims = simulationstack(nsim, sim, samp, plan) do
    executeplan!(
        samp,
        plan,
        stencil=stencil,
        plannedjitter=plannedjitter,
        samplerjitter=samplerjitter,
        corejitter=corejitter,
    )

    updategaussian!(Q, samp.points, Q_cov)
    spreading!(rng, sim, plan, Q)

    ${getMixingCode(params.mixing?.type, params.mixing?.minDepth, params.mixing?.maxDepth)}

    updategaussian!(ρs, samp.points, ρs_cov)
    rand!(rng, ρs, view(sim.ρs, :))

    feedstockconcentration!(rng, sim, (Ca=${params.feedstock?.caConcentration || 0.07}, Mg=${params.feedstock?.mgConcentration || 0.05}), ${params.feedstock?.concentrationVariability || 0.03})

    updategaussian!(cs, samp.points, cs_cov)
    soilconcentration!(rng, sim, cs)

    leaching!(sim, :Ca, Ca_leaching, plan)
    leaching!(sim, :Mg, Mg_leaching, plan)

    massloss!(sim, plan, x -> (x[:Ca] + x[:Mg]) / 2)

    analyze!(rng, sim, (Ca=${params.execution?.caMeasurementError || 0.03}, Mg=${params.execution?.mgMeasurementError || 0.03}), ${params.execution?.massMeasurementError || 0.005})
end

# Convert simulation results to JSON format
all_results = []
for realization in 1:nsim
    for sample in 1:size(sims[:data], 3)
        push!(all_results, Dict(
            "realization" => realization,
            "sample" => sample,
            "time" => sims[:time][sample],
            "Ca" => sims[:data][realization, 1, sample],
            "Mg" => sims[:data][realization, 2, sample], 
            "mass" => sims[:data][realization, 3, sample],
            "control" => sims[:control][sample],
            "x" => sims[:x][realization, sample],
            "y" => sims[:y][realization, sample]
        ))
    end
end

# Calculate summary statistics
ca_values = [r["Ca"] for r in all_results if !isnan(r["Ca"])]
mg_values = [r["Mg"] for r in all_results if !isnan(r["Mg"])]

summary_stats = Dict(
    "Ca" => Dict(
        "mean" => isempty(ca_values) ? NaN : mean(ca_values),
        "std" => isempty(ca_values) ? NaN : std(ca_values),
        "min" => isempty(ca_values) ? NaN : minimum(ca_values),
        "max" => isempty(ca_values) ? NaN : maximum(ca_values)
    ),
    "Mg" => Dict(
        "mean" => isempty(mg_values) ? NaN : mean(mg_values),
        "std" => isempty(mg_values) ? NaN : std(mg_values),
        "min" => isempty(mg_values) ? NaN : minimum(mg_values),
        "max" => isempty(mg_values) ? NaN : maximum(mg_values)
    )
)

# Prepare output
output = Dict(
    "data" => all_results,
    "summary" => summary_stats,
    "parameters" => $(JSON.stringify(params))
)

# Write results to JSON file
open("simulation_results.json", "w") do f
    JSON.print(f, output)
end

println("Simulation completed successfully!")
`
}

// Helper function to generate stencil code
function getStencilCode(stencilType, stencilSize, coresPerSample) {
  const size = stencilSize || 2.0;
  const cores = coresPerSample || 5;
  
  switch (stencilType) {
    case 'circle':
      return `stencil = CircleStencil(${cores}, ${size})`;
    case 'hubspoke':
      return `stencil = HubSpokeStencil(${cores}, ${size})`;
    case 'line':
      return `stencil = LineStencil(${cores}, ${size})`;
    case 'random':
      return `stencil = RandomStencil(${cores}, ${size})`;
    default:
      return `stencil = CircleStencil(${cores}, ${size})`;
  }
}

// Helper function to generate mixing code
function getMixingCode(mixingType, minDepth, maxDepth) {
  const min = minDepth || 0.05;
  const max = maxDepth || 0.15;
  
  switch (mixingType) {
    case 'triangular':
      return `triangularmixing!(rng, sim, depth=TriangularDist(${min}, ${max}), upper=Uniform(0.0, 0.5))`;
    case 'uniform':
      return `uniformmixing!(rng, sim, depth=TriangularDist(${min}, ${max}), upper=Uniform(0.0, 0.5))`;
    case 'exponential':
      return `exponentialmixing!(rng, sim, depth=TriangularDist(${min}, ${max}), scale=Uniform(0.01, 0.1))`;
    case 'unmixed':
    default:
      return `unmixed!(rng, sim, depth=TriangularDist(${min}, ${max}))`;
  }
}

// API Routes
app.get('/api/status', async (req, res) => {
  const juliaAvailable = await checkJulia()
  res.json({
    julia: juliaAvailable,
    server: 'running',
    timestamp: new Date().toISOString()
  })
})

app.post('/api/simulate', async (req, res) => {
  try {
    const params = req.body
    
    // Validate parameters
    if (!params || typeof params !== 'object') {
      return res.status(400).json({ error: 'Invalid parameters' })
    }
    
    // Check if Julia is available
    const juliaAvailable = await checkJulia()
    if (!juliaAvailable) {
      return res.status(500).json({ 
        error: 'Julia is not available. Please install Julia and ensure it is in your PATH.' 
      })
    }
    
    console.log('Starting Monty simulation with parameters:', params)
    
    // Generate Julia script
    const script = generateJuliaScript(params)
    const scriptPath = join(__dirname, 'temp_simulation.jl')
    
    // Write script to temporary file
    await fs.writeFile(scriptPath, script)
    
    // Run Julia simulation
    const julia = spawn('julia', [scriptPath], {
      cwd: join(__dirname, '..'), // Run from project root
      stdio: ['pipe', 'pipe', 'pipe']
    })
    
    let stdout = ''
    let stderr = ''
    
    julia.stdout.on('data', (data) => {
      stdout += data.toString()
      console.log('Julia stdout:', data.toString())
    })
    
    julia.stderr.on('data', (data) => {
      stderr += data.toString()
      console.error('Julia stderr:', data.toString())
    })
    
    julia.on('close', async (code) => {
      try {
        // Clean up script file
        await fs.unlink(scriptPath).catch(() => {})
        
        if (code !== 0) {
          console.error('Julia process failed with code:', code)
          console.error('stderr:', stderr)
          return res.status(500).json({ 
            error: `Simulation failed with exit code ${code}`,
            details: stderr
          })
        }
        
        // Read results
        try {
          const resultsPath = join(__dirname, '..', 'simulation_results.json')
          const resultsData = await fs.readFile(resultsPath, 'utf8')
          const results = JSON.parse(resultsData)
          
          // Clean up results file
          await fs.unlink(resultsPath).catch(() => {})
          
          console.log('Simulation completed successfully')
          res.json(results)
        } catch (fileError) {
          console.error('Error reading results:', fileError)
          res.status(500).json({ 
            error: 'Failed to read simulation results',
            details: fileError.message
          })
        }
      } catch (error) {
        console.error('Error in simulation cleanup:', error)
        res.status(500).json({ 
          error: 'Internal server error',
          details: error.message
        })
      }
    })
    
    julia.on('error', (error) => {
      console.error('Failed to start Julia process:', error)
      res.status(500).json({ 
        error: 'Failed to start Julia process',
        details: error.message
      })
    })
    
  } catch (error) {
    console.error('Simulation error:', error)
    res.status(500).json({ 
      error: 'Internal server error',
      details: error.message
    })
  }
})

app.listen(PORT, () => {
  console.log(`Monty simulation server running on port ${PORT}`)
  console.log(`API available at http://localhost:${PORT}/api`)
})