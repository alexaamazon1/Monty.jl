import React, { useState } from 'react'
import { Play, Settings, BarChart3, Download, Loader2 } from 'lucide-react'
import SimulationForm from './components/SimulationForm'
import ResultsVisualization from './components/ResultsVisualization'
import { runSimulation } from './services/api'

function App() {
  const [simulationParams, setSimulationParams] = useState({
    // Grid configuration (from main.jl: 8x8 grid, 10x10m cells)
    grid: {
      xSize: 8,
      ySize: 8,
      cellWidth: 10.0,
      cellHeight: 10.0
    },
    
    // Sampling configuration (from main.jl)
    sampling: {
      coresPerSample: 5,
      timePoints: "-0.003, 0.0, 1.0",
      stencilType: "circle",
      stencilSize: 2.0
    },
    
    // Feedstock parameters (from main.jl)
    feedstock: {
      applicationRate: 3.5,           // Q_μ
      density: 1000,                  // sim.ρf
      caConcentration: 0.07,          // feedstock Ca
      mgConcentration: 0.05,          // feedstock Mg
      concentrationVariability: 0.03 // RSD for concentrations
    },
    
    // Soil parameters (from main.jl)
    soil: {
      density: 1000,                  // soil density mean
      caConcentration: 0.002,         // background Ca
      mgConcentration: 0.001,         // background Mg
      crossCorrelation: 0.75          // cross-correlation for Ca/Mg
    },
    
    // Mixing model (from main.jl: unmixed with triangular depth)
    mixing: {
      type: "unmixed",
      minDepth: 0.05,
      maxDepth: 0.15
    },
    
    // Leaching models (from main.jl)
    leaching: {
      caRate: 0.4,                    // Ca_leaching λ
      mgRate: 0.8,                    // Mg_leaching λ
      modelType: "exponential"
    },
    
    // Jitter parameters (from main.jl)
    jitter: {
      sampler: 0.75,                  // samplerjitter
      core: 0.1,                      // corejitter
      planned: 5.0                    // plannedjitter GridCentroidJitter
    },
    
    // Covariance parameters (from main.jl)
    covariance: {
      applicationType: "gaussian",
      appRangeX: 5.0,                 // Q_cov MetricBall
      appRangeY: 500.0,               // Q_cov MetricBall
      soilRange: 20.0,                // cs_cov SphericalCovariance range
      densityRange: 30.0              // ρs_cov SphericalCovariance range
    },
    
    // Execution parameters (from main.jl)
    execution: {
      numRealizations: 10000,         // nsim
      randomSeed: 1,
      caMeasurementError: 0.03,       // analyze! Ca error
      mgMeasurementError: 0.03,       // analyze! Mg error
      massMeasurementError: 0.005     // analyze! mass error
    }
  })
  
  const [results, setResults] = useState(null)
  const [isRunning, setIsRunning] = useState(false)
  const [error, setError] = useState(null)

  const handleRunSimulation = async () => {
    setIsRunning(true)
    setError(null)
    
    try {
      const simulationResults = await runSimulation(simulationParams)
      setResults(simulationResults)
    } catch (err) {
      setError(err.message || 'Failed to run simulation')
      console.error('Simulation error:', err)
    } finally {
      setIsRunning(false)
    }
  }

  const handleDownloadResults = () => {
    if (!results) return
    
    const dataStr = JSON.stringify(results, null, 2)
    const dataBlob = new Blob([dataStr], { type: 'application/json' })
    const url = URL.createObjectURL(dataBlob)
    const link = document.createElement('a')
    link.href = url
    link.download = `monty_simulation_${new Date().toISOString().split('T')[0]}.json`
    link.click()
    URL.revokeObjectURL(url)
  }

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100">
      <div className="container mx-auto px-4 py-8">
        {/* Header */}
        <div className="text-center mb-8">
          <h1 className="text-4xl font-bold text-gray-900 mb-2">
            Monty Simulation Platform
          </h1>
          <p className="text-lg text-gray-600 max-w-2xl mx-auto">
            Simulate geochemical data for enhanced rock weathering (ERW) field trials 
            and commercial deployments with advanced mixing and leaching models.
          </p>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {/* Simulation Parameters */}
          <div className="lg:col-span-1">
            <div className="card">
              <div className="flex items-center gap-2 mb-4">
                <Settings className="w-5 h-5 text-primary-600" />
                <h2 className="text-xl font-semibold">Simulation Parameters</h2>
              </div>
              
              <SimulationForm 
                params={simulationParams}
                onChange={setSimulationParams}
              />
              
              <div className="mt-6 space-y-3">
                <button
                  onClick={handleRunSimulation}
                  disabled={isRunning}
                  className="btn-primary w-full flex items-center justify-center gap-2"
                >
                  {isRunning ? (
                    <>
                      <Loader2 className="w-4 h-4 animate-spin" />
                      Running Simulation...
                    </>
                  ) : (
                    <>
                      <Play className="w-4 h-4" />
                      Run Simulation
                    </>
                  )}
                </button>
                
                {results && (
                  <button
                    onClick={handleDownloadResults}
                    className="btn-secondary w-full flex items-center justify-center gap-2"
                  >
                    <Download className="w-4 h-4" />
                    Download Results
                  </button>
                )}
              </div>
              
              {error && (
                <div className="mt-4 p-3 bg-red-50 border border-red-200 rounded-lg">
                  <p className="text-red-700 text-sm">{error}</p>
                </div>
              )}
            </div>
          </div>

          {/* Results Visualization */}
          <div className="lg:col-span-2">
            <div className="card">
              <div className="flex items-center gap-2 mb-4">
                <BarChart3 className="w-5 h-5 text-primary-600" />
                <h2 className="text-xl font-semibold">Simulation Results</h2>
              </div>
              
              {isRunning ? (
                <div className="flex items-center justify-center h-96">
                  <div className="text-center">
                    <Loader2 className="w-8 h-8 animate-spin text-primary-600 mx-auto mb-4" />
                    <p className="text-gray-600">Running Monty simulation...</p>
                    <p className="text-sm text-gray-500 mt-2">
                      This may take a few moments depending on the number of realizations
                    </p>
                  </div>
                </div>
              ) : results ? (
                <ResultsVisualization results={results} />
              ) : (
                <div className="flex items-center justify-center h-96">
                  <div className="text-center text-gray-500">
                    <BarChart3 className="w-16 h-16 mx-auto mb-4 opacity-50" />
                    <p className="text-lg">No simulation results yet</p>
                    <p className="text-sm">Configure parameters and run a simulation to see results</p>
                  </div>
                </div>
              )}
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}

export default App