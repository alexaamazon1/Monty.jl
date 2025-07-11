import React from 'react'

const SimulationForm = ({ params, onChange }) => {
  const handleChange = (field, value) => {
    onChange(prev => ({
      ...prev,
      [field]: value
    }))
  }

  const handleTimePointsChange = (value) => {
    try {
      const timePoints = value.split(',').map(t => parseFloat(t.trim())).filter(t => !isNaN(t))
      onChange(prev => ({
        ...prev,
        timePoints
      }))
    } catch (err) {
      // Invalid input, ignore
    }
  }

  return (
    <div className="space-y-4">
      {/* Feedstock Parameters */}
      <div>
        <h3 className="text-sm font-medium text-gray-700 mb-3">Feedstock Parameters</h3>
        <div className="space-y-3">
          <div>
            <label className="block text-xs text-gray-600 mb-1">
              Application Rate (kg/m²)
            </label>
            <input
              type="number"
              step="0.1"
              value={params.applicationRate}
              onChange={(e) => handleChange('applicationRate', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
          
          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="block text-xs text-gray-600 mb-1">
                Ca Concentration
              </label>
              <input
                type="number"
                step="0.001"
                value={params.feedstockCa}
                onChange={(e) => handleChange('feedstockCa', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">
                Mg Concentration
              </label>
              <input
                type="number"
                step="0.001"
                value={params.feedstockMg}
                onChange={(e) => handleChange('feedstockMg', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
          </div>
        </div>
      </div>

      {/* Soil Parameters */}
      <div>
        <h3 className="text-sm font-medium text-gray-700 mb-3">Soil Parameters</h3>
        <div className="grid grid-cols-2 gap-2">
          <div>
            <label className="block text-xs text-gray-600 mb-1">
              Soil Ca (kg/kg)
            </label>
            <input
              type="number"
              step="0.0001"
              value={params.soilCa}
              onChange={(e) => handleChange('soilCa', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">
              Soil Mg (kg/kg)
            </label>
            <input
              type="number"
              step="0.0001"
              value={params.soilMg}
              onChange={(e) => handleChange('soilMg', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
        </div>
      </div>

      {/* Leaching Parameters */}
      <div>
        <h3 className="text-sm font-medium text-gray-700 mb-3">Leaching Rates</h3>
        <div className="grid grid-cols-2 gap-2">
          <div>
            <label className="block text-xs text-gray-600 mb-1">
              Ca Leaching Rate
            </label>
            <input
              type="number"
              step="0.1"
              value={params.leachingRateCa}
              onChange={(e) => handleChange('leachingRateCa', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">
              Mg Leaching Rate
            </label>
            <input
              type="number"
              step="0.1"
              value={params.leachingRateMg}
              onChange={(e) => handleChange('leachingRateMg', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
        </div>
      </div>

      {/* Sampling Parameters */}
      <div>
        <h3 className="text-sm font-medium text-gray-700 mb-3">Sampling Parameters</h3>
        <div className="space-y-3">
          <div>
            <label className="block text-xs text-gray-600 mb-1">
              Sampling Depth (m)
            </label>
            <input
              type="number"
              step="0.01"
              value={params.samplingDepth}
              onChange={(e) => handleChange('samplingDepth', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
          
          <div>
            <label className="block text-xs text-gray-600 mb-1">
              Number of Samples
            </label>
            <input
              type="number"
              value={params.numSamples}
              onChange={(e) => handleChange('numSamples', parseInt(e.target.value))}
              className="input-field"
            />
          </div>
          
          <div>
            <label className="block text-xs text-gray-600 mb-1">
              Number of Realizations
            </label>
            <input
              type="number"
              value={params.numRealizations}
              onChange={(e) => handleChange('numRealizations', parseInt(e.target.value))}
              className="input-field"
            />
          </div>
          
          <div>
            <label className="block text-xs text-gray-600 mb-1">
              Time Points (years, comma-separated)
            </label>
            <input
              type="text"
              value={params.timePoints.join(', ')}
              onChange={(e) => handleTimePointsChange(e.target.value)}
              className="input-field"
              placeholder="0, 0.5, 1.0"
            />
          </div>
        </div>
      </div>
    </div>
  )
}

export default SimulationForm