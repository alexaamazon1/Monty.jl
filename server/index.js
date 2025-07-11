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
  return `
using Pkg
Pkg.activate(".")

using Monty
using Random
using Distributions
using Meshes
using JSON

# Set random seed for reproducibility
const rng = Xoshiro(42)

# Simulation parameters
const application_rate = ${params.applicationRate}
const feedstock_ca = ${params.feedstockCa}
const feedstock_mg = ${params.feedstockMg}
const soil_ca = ${params.soilCa}
const soil_mg = ${params.soilMg}
const leaching_rate_ca = ${params.leachingRateCa}
const leaching_rate_mg = ${params.leachingRateMg}
const sampling_depth = ${params.samplingDepth}
const num_samples = ${params.numSamples}
const num_realizations = ${params.numRealizations}
const time_points = [${params.timePoints.join(', ')}]

# Create a simple field geometry
field = Ball((0.0, 0.0), 50.0)

# Create sample plan
plan = randomsampleplan(rng, field, num_samples, time_points)
samp = CoreSet(plan, 5)  # 5 cores per sample
sim = Simulation((:Ca, :Mg), samp)

# Run multiple realizations
all_results = []

for realization in 1:num_realizations
    # Execute plan with some jitter
    executeplan!(
        samp,
        plan,
        stencil=CircleStencil(5, 1.0),
        samplerjitter=Jitter(1.0),
        corejitter=Jitter(0.1)
    )
    
    # Set up simulation parameters
    unmixed!(rng, sim, depth=TriangularDist(sampling_depth * 0.8, sampling_depth * 1.2))
    
    # Application rate with some variability
    rand!(rng, Normal(application_rate, application_rate * 0.1), sim.Q)
    
    # Feedstock properties
    sim.ρf .= 2e3  # feedstock density
    feedstockconcentration!(rng, sim, (Ca=feedstock_ca, Mg=feedstock_mg), 0.05)
    
    # Soil properties
    rand!(rng, Normal(1e3, 100), sim.ρs)  # soil density
    rand!(rng, Normal(soil_ca, soil_ca * 0.1), sim.cs[:Ca])
    rand!(rng, Normal(soil_mg, soil_mg * 0.1), sim.cs[:Mg])
    
    # Leaching models
    leaching!(sim, :Ca, ExponentialLeaching(λ=leaching_rate_ca), plan)
    leaching!(sim, :Mg, ExponentialLeaching(λ=leaching_rate_mg), plan)
    
    # Mass loss (simple average of elemental losses)
    massloss!(sim, plan, x -> (x[:Ca] + x[:Mg]) / 2)
    
    # Analyze samples with measurement error
    analyze!(rng, sim, (Ca=0.03, Mg=0.03), 0.005)
    
    # Extract results for this realization
    for i in 1:length(sim.measurements)
        push!(all_results, Dict(
            "realization" => realization,
            "sample" => i,
            "time" => plan.time[i],
            "Ca" => sim.measurements[i][:Ca],
            "Mg" => sim.measurements[i][:Mg],
            "mass" => sim.measurements[i].mass,
            "control" => plan.control[i]
        ))
    end
end

# Calculate summary statistics
ca_values = [r["Ca"] for r in all_results]
mg_values = [r["Mg"] for r in all_results]

summary_stats = Dict(
    "Ca" => Dict(
        "mean" => mean(ca_values),
        "std" => std(ca_values),
        "min" => minimum(ca_values),
        "max" => maximum(ca_values)
    ),
    "Mg" => Dict(
        "mean" => mean(mg_values),
        "std" => std(mg_values),
        "min" => minimum(mg_values),
        "max" => maximum(mg_values)
    )
)

# Prepare output
output = Dict(
    "data" => all_results,
    "summary" => summary_stats,
    "parameters" => Dict(
        "applicationRate" => application_rate,
        "feedstockCa" => feedstock_ca,
        "feedstockMg" => feedstock_mg,
        "soilCa" => soil_ca,
        "soilMg" => soil_mg,
        "leachingRateCa" => leaching_rate_ca,
        "leachingRateMg" => leaching_rate_mg,
        "samplingDepth" => sampling_depth,
        "numSamples" => num_samples,
        "numRealizations" => num_realizations,
        "timePoints" => time_points
    )
)

# Write results to JSON file
open("simulation_results.json", "w") do f
    JSON.print(f, output)
end

println("Simulation completed successfully!")
`
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