import React, { useState } from 'react'
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, ScatterChart, Scatter } from 'recharts'

const ResultsVisualization = ({ results }) => {
  const [activeTab, setActiveTab] = useState('concentrations')

  if (!results || !results.data) {
    return (
      <div className="text-center text-gray-500 py-8">
        <p>No results data available</p>
      </div>
    )
  }

  const { data, summary, parameters } = results

  // Process data for visualization
  const processConcentrationData = () => {
    const timePoints = [...new Set(data.map(d => d.time))].sort((a, b) => a - b)
    
    return timePoints.map(time => {
      const timeData = data.filter(d => d.time === time)
      const caValues = timeData.map(d => d.Ca * 1e6) // Convert to ppm
      const mgValues = timeData.map(d => d.Mg * 1e6) // Convert to ppm
      
      return {
        time,
        Ca_mean: caValues.reduce((a, b) => a + b, 0) / caValues.length,
        Ca_std: Math.sqrt(caValues.reduce((a, b) => a + (b - caValues.reduce((x, y) => x + y, 0) / caValues.length) ** 2, 0) / caValues.length),
        Mg_mean: mgValues.reduce((a, b) => a + b, 0) / mgValues.length,
        Mg_std: Math.sqrt(mgValues.reduce((a, b) => a + (b - mgValues.reduce((x, y) => x + y, 0) / mgValues.length) ** 2, 0) / mgValues.length),
        Ca_min: Math.min(...caValues),
        Ca_max: Math.max(...caValues),
        Mg_min: Math.min(...mgValues),
        Mg_max: Math.max(...mgValues)
      }
    })
  }

  const concentrationData = processConcentrationData()

  const tabs = [
    { id: 'concentrations', label: 'Concentrations Over Time' },
    { id: 'summary', label: 'Summary Statistics' },
    { id: 'parameters', label: 'Simulation Parameters' }
  ]

  return (
    <div className="space-y-4">
      {/* Tab Navigation */}
      <div className="flex space-x-1 bg-gray-100 p-1 rounded-lg">
        {tabs.map(tab => (
          <button
            key={tab.id}
            onClick={() => setActiveTab(tab.id)}
            className={`px-3 py-2 text-sm font-medium rounded-md transition-colors ${
              activeTab === tab.id
                ? 'bg-white text-primary-600 shadow-sm'
                : 'text-gray-600 hover:text-gray-900'
            }`}
          >
            {tab.label}
          </button>
        ))}
      </div>

      {/* Tab Content */}
      <div className="min-h-96">
        {activeTab === 'concentrations' && (
          <div className="space-y-6">
            <div>
              <h3 className="text-lg font-medium mb-4">Calcium Concentrations</h3>
              <ResponsiveContainer width="100%" height={300}>
                <LineChart data={concentrationData}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis 
                    dataKey="time" 
                    label={{ value: 'Time (years)', position: 'insideBottom', offset: -5 }}
                  />
                  <YAxis 
                    label={{ value: 'Concentration (ppm)', angle: -90, position: 'insideLeft' }}
                  />
                  <Tooltip 
                    formatter={(value, name) => [value.toFixed(2), name]}
                    labelFormatter={(value) => `Time: ${value} years`}
                  />
                  <Legend />
                  <Line 
                    type="monotone" 
                    dataKey="Ca_mean" 
                    stroke="#2563eb" 
                    strokeWidth={2}
                    name="Mean Ca"
                    dot={{ fill: '#2563eb', strokeWidth: 2, r: 4 }}
                  />
                  <Line 
                    type="monotone" 
                    dataKey="Ca_min" 
                    stroke="#93c5fd" 
                    strokeDasharray="5 5"
                    name="Min Ca"
                    dot={false}
                  />
                  <Line 
                    type="monotone" 
                    dataKey="Ca_max" 
                    stroke="#93c5fd" 
                    strokeDasharray="5 5"
                    name="Max Ca"
                    dot={false}
                  />
                </LineChart>
              </ResponsiveContainer>
            </div>

            <div>
              <h3 className="text-lg font-medium mb-4">Magnesium Concentrations</h3>
              <ResponsiveContainer width="100%" height={300}>
                <LineChart data={concentrationData}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis 
                    dataKey="time" 
                    label={{ value: 'Time (years)', position: 'insideBottom', offset: -5 }}
                  />
                  <YAxis 
                    label={{ value: 'Concentration (ppm)', angle: -90, position: 'insideLeft' }}
                  />
                  <Tooltip 
                    formatter={(value, name) => [value.toFixed(2), name]}
                    labelFormatter={(value) => `Time: ${value} years`}
                  />
                  <Legend />
                  <Line 
                    type="monotone" 
                    dataKey="Mg_mean" 
                    stroke="#dc2626" 
                    strokeWidth={2}
                    name="Mean Mg"
                    dot={{ fill: '#dc2626', strokeWidth: 2, r: 4 }}
                  />
                  <Line 
                    type="monotone" 
                    dataKey="Mg_min" 
                    stroke="#fca5a5" 
                    strokeDasharray="5 5"
                    name="Min Mg"
                    dot={false}
                  />
                  <Line 
                    type="monotone" 
                    dataKey="Mg_max" 
                    stroke="#fca5a5" 
                    strokeDasharray="5 5"
                    name="Max Mg"
                    dot={false}
                  />
                </LineChart>
              </ResponsiveContainer>
            </div>
          </div>
        )}

        {activeTab === 'summary' && summary && (
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div className="space-y-4">
              <h3 className="text-lg font-medium">Calcium Statistics</h3>
              <div className="bg-gray-50 p-4 rounded-lg space-y-2">
                <div className="flex justify-between">
                  <span className="text-gray-600">Mean:</span>
                  <span className="font-medium">{(summary.Ca.mean * 1e6).toFixed(2)} ppm</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Std Dev:</span>
                  <span className="font-medium">{(summary.Ca.std * 1e6).toFixed(2)} ppm</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Min:</span>
                  <span className="font-medium">{(summary.Ca.min * 1e6).toFixed(2)} ppm</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Max:</span>
                  <span className="font-medium">{(summary.Ca.max * 1e6).toFixed(2)} ppm</span>
                </div>
              </div>
            </div>

            <div className="space-y-4">
              <h3 className="text-lg font-medium">Magnesium Statistics</h3>
              <div className="bg-gray-50 p-4 rounded-lg space-y-2">
                <div className="flex justify-between">
                  <span className="text-gray-600">Mean:</span>
                  <span className="font-medium">{(summary.Mg.mean * 1e6).toFixed(2)} ppm</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Std Dev:</span>
                  <span className="font-medium">{(summary.Mg.std * 1e6).toFixed(2)} ppm</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Min:</span>
                  <span className="font-medium">{(summary.Mg.min * 1e6).toFixed(2)} ppm</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Max:</span>
                  <span className="font-medium">{(summary.Mg.max * 1e6).toFixed(2)} ppm</span>
                </div>
              </div>
            </div>
          </div>
        )}

        {activeTab === 'parameters' && parameters && (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            <div className="space-y-4">
              <h3 className="text-lg font-medium">Feedstock</h3>
              <div className="bg-gray-50 p-4 rounded-lg space-y-2 text-sm">
                <div className="flex justify-between">
                  <span className="text-gray-600">Application Rate:</span>
                  <span className="font-medium">{parameters.applicationRate} kg/m²</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Ca Concentration:</span>
                  <span className="font-medium">{parameters.feedstockCa}</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Mg Concentration:</span>
                  <span className="font-medium">{parameters.feedstockMg}</span>
                </div>
              </div>
            </div>

            <div className="space-y-4">
              <h3 className="text-lg font-medium">Soil</h3>
              <div className="bg-gray-50 p-4 rounded-lg space-y-2 text-sm">
                <div className="flex justify-between">
                  <span className="text-gray-600">Ca Concentration:</span>
                  <span className="font-medium">{parameters.soilCa} kg/kg</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Mg Concentration:</span>
                  <span className="font-medium">{parameters.soilMg} kg/kg</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Sampling Depth:</span>
                  <span className="font-medium">{parameters.samplingDepth} m</span>
                </div>
              </div>
            </div>

            <div className="space-y-4">
              <h3 className="text-lg font-medium">Simulation</h3>
              <div className="bg-gray-50 p-4 rounded-lg space-y-2 text-sm">
                <div className="flex justify-between">
                  <span className="text-gray-600">Samples:</span>
                  <span className="font-medium">{parameters.numSamples}</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Realizations:</span>
                  <span className="font-medium">{parameters.numRealizations}</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-gray-600">Time Points:</span>
                  <span className="font-medium">{parameters.timePoints.join(', ')} yr</span>
                </div>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  )
}

export default ResultsVisualization