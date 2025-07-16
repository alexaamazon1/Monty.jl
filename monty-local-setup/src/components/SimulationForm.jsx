import React, { useState } from 'react'
import { ChevronDown, ChevronRight, Info } from 'lucide-react'

const SimulationForm = ({ params, onChange }) => {
  const [expandedSections, setExpandedSections] = useState({
    grid: true,
    sampling: true,
    feedstock: true,
    soil: false,
    mixing: false,
    leaching: false,
    jitter: false,
    covariance: false,
    execution: false
  })

  const toggleSection = (section) => {
    setExpandedSections(prev => ({
      ...prev,
      [section]: !prev[section]
    }))
  }

  const handleChange = (field, value) => {
    onChange(prev => ({
      ...prev,
      [field]: value
    }))
  }

  const handleNestedChange = (section, field, value) => {
    onChange(prev => ({
      ...prev,
      [section]: {
        ...prev[section],
        [field]: value
      }
    }))
  }

  const Section = ({ title, sectionKey, children, tooltip }) => (
    <div className="border border-gray-200 rounded-lg mb-4">
      <button
        type="button"
        onClick={() => toggleSection(sectionKey)}
        className="w-full px-4 py-3 flex items-center justify-between bg-gray-50 hover:bg-gray-100 rounded-t-lg"
      >
        <div className="flex items-center gap-2">
          {expandedSections[sectionKey] ? (
            <ChevronDown className="w-4 h-4" />
          ) : (
            <ChevronRight className="w-4 h-4" />
          )}
          <span className="font-medium text-gray-900">{title}</span>
          {tooltip && <Info className="w-4 h-4 text-gray-500" title={tooltip} />}
        </div>
      </button>
      {expandedSections[sectionKey] && (
        <div className="p-4 space-y-4">
          {children}
        </div>
      )}
    </div>
  )

  return (
    <div className="space-y-0">
      {/* Grid Configuration */}
      <Section 
        title="Grid Configuration" 
        sectionKey="grid"
        tooltip="Configure the spatial grid for sampling locations"
      >
        <div className="grid grid-cols-2 gap-3">
          <div>
            <label className="block text-xs text-gray-600 mb-1">Grid X Size</label>
            <input
              type="number"
              value={params.grid?.xSize || 8}
              onChange={(e) => handleNestedChange('grid', 'xSize', parseInt(e.target.value))}
              className="input-field"
              min="1"
              max="20"
            />
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Grid Y Size</label>
            <input
              type="number"
              value={params.grid?.ySize || 8}
              onChange={(e) => handleNestedChange('grid', 'ySize', parseInt(e.target.value))}
              className="input-field"
              min="1"
              max="20"
            />
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Cell Width (m)</label>
            <input
              type="number"
              step="0.1"
              value={params.grid?.cellWidth || 10}
              onChange={(e) => handleNestedChange('grid', 'cellWidth', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Cell Height (m)</label>
            <input
              type="number"
              step="0.1"
              value={params.grid?.cellHeight || 10}
              onChange={(e) => handleNestedChange('grid', 'cellHeight', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
        </div>
      </Section>

      {/* Sampling Configuration */}
      <Section 
        title="Sampling Configuration" 
        sectionKey="sampling"
        tooltip="Configure sampling design and timing"
      >
        <div className="space-y-3">
          <div>
            <label className="block text-xs text-gray-600 mb-1">Cores per Sample</label>
            <input
              type="number"
              value={params.sampling?.coresPerSample || 5}
              onChange={(e) => handleNestedChange('sampling', 'coresPerSample', parseInt(e.target.value))}
              className="input-field"
              min="1"
              max="10"
            />
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Sampling Times (years, comma-separated)</label>
            <input
              type="text"
              value={params.sampling?.timePoints || "-0.003, 0.0, 1.0"}
              onChange={(e) => handleNestedChange('sampling', 'timePoints', e.target.value)}
              className="input-field"
              placeholder="-0.003, 0.0, 1.0"
            />
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Core Stencil Type</label>
            <select
              value={params.sampling?.stencilType || "circle"}
              onChange={(e) => handleNestedChange('sampling', 'stencilType', e.target.value)}
              className="input-field"
            >
              <option value="circle">Circle Stencil</option>
              <option value="hubspoke">Hub-Spoke Stencil</option>
              <option value="line">Line Stencil</option>
              <option value="random">Random Stencil</option>
            </select>
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Stencil Radius/Size (m)</label>
            <input
              type="number"
              step="0.1"
              value={params.sampling?.stencilSize || 2.0}
              onChange={(e) => handleNestedChange('sampling', 'stencilSize', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
        </div>
      </Section>

      {/* Feedstock Parameters */}
      <Section 
        title="Feedstock Parameters" 
        sectionKey="feedstock"
        tooltip="Configure feedstock properties and application"
      >
        <div className="space-y-3">
          <div>
            <label className="block text-xs text-gray-600 mb-1">Application Rate Mean (kg/m²)</label>
            <input
              type="number"
              step="0.1"
              value={params.feedstock?.applicationRate || 3.5}
              onChange={(e) => handleNestedChange('feedstock', 'applicationRate', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Feedstock Density (kg/m³)</label>
            <input
              type="number"
              value={params.feedstock?.density || 1000}
              onChange={(e) => handleNestedChange('feedstock', 'density', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="block text-xs text-gray-600 mb-1">Ca Concentration</label>
              <input
                type="number"
                step="0.001"
                value={params.feedstock?.caConcentration || 0.07}
                onChange={(e) => handleNestedChange('feedstock', 'caConcentration', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">Mg Concentration</label>
              <input
                type="number"
                step="0.001"
                value={params.feedstock?.mgConcentration || 0.05}
                onChange={(e) => handleNestedChange('feedstock', 'mgConcentration', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Concentration Variability (RSD)</label>
            <input
              type="number"
              step="0.001"
              value={params.feedstock?.concentrationVariability || 0.03}
              onChange={(e) => handleNestedChange('feedstock', 'concentrationVariability', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
        </div>
      </Section>

      {/* Soil Parameters */}
      <Section 
        title="Soil Parameters" 
        sectionKey="soil"
        tooltip="Configure background soil properties"
      >
        <div className="space-y-3">
          <div>
            <label className="block text-xs text-gray-600 mb-1">Soil Density Mean (kg/m³)</label>
            <input
              type="number"
              value={params.soil?.density || 1000}
              onChange={(e) => handleNestedChange('soil', 'density', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="block text-xs text-gray-600 mb-1">Background Ca (kg/kg)</label>
              <input
                type="number"
                step="0.0001"
                value={params.soil?.caConcentration || 0.002}
                onChange={(e) => handleNestedChange('soil', 'caConcentration', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">Background Mg (kg/kg)</label>
              <input
                type="number"
                step="0.0001"
                value={params.soil?.mgConcentration || 0.001}
                onChange={(e) => handleNestedChange('soil', 'mgConcentration', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Cross-correlation</label>
            <input
              type="number"
              step="0.1"
              min="0"
              max="1"
              value={params.soil?.crossCorrelation || 0.75}
              onChange={(e) => handleNestedChange('soil', 'crossCorrelation', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
        </div>
      </Section>

      {/* Mixing Model */}
      <Section 
        title="Mixing Model" 
        sectionKey="mixing"
        tooltip="Configure how feedstock mixes with soil"
      >
        <div className="space-y-3">
          <div>
            <label className="block text-xs text-gray-600 mb-1">Mixing Type</label>
            <select
              value={params.mixing?.type || "unmixed"}
              onChange={(e) => handleNestedChange('mixing', 'type', e.target.value)}
              className="input-field"
            >
              <option value="unmixed">Unmixed (Layered)</option>
              <option value="triangular">Triangular Mixing</option>
              <option value="uniform">Uniform Mixing</option>
              <option value="exponential">Exponential Mixing</option>
            </select>
          </div>
          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="block text-xs text-gray-600 mb-1">Min Depth (m)</label>
              <input
                type="number"
                step="0.01"
                value={params.mixing?.minDepth || 0.05}
                onChange={(e) => handleNestedChange('mixing', 'minDepth', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">Max Depth (m)</label>
              <input
                type="number"
                step="0.01"
                value={params.mixing?.maxDepth || 0.15}
                onChange={(e) => handleNestedChange('mixing', 'maxDepth', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
          </div>
        </div>
      </Section>

      {/* Leaching Models */}
      <Section 
        title="Leaching Models" 
        sectionKey="leaching"
        tooltip="Configure element leaching behavior"
      >
        <div className="space-y-3">
          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="block text-xs text-gray-600 mb-1">Ca Leaching Rate (λ)</label>
              <input
                type="number"
                step="0.1"
                value={params.leaching?.caRate || 0.4}
                onChange={(e) => handleNestedChange('leaching', 'caRate', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">Mg Leaching Rate (λ)</label>
              <input
                type="number"
                step="0.1"
                value={params.leaching?.mgRate || 0.8}
                onChange={(e) => handleNestedChange('leaching', 'mgRate', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Leaching Model Type</label>
            <select
              value={params.leaching?.modelType || "exponential"}
              onChange={(e) => handleNestedChange('leaching', 'modelType', e.target.value)}
              className="input-field"
            >
              <option value="exponential">Exponential Leaching</option>
              <option value="linear">Linear Leaching</option>
              <option value="power">Power Law Leaching</option>
            </select>
          </div>
        </div>
      </Section>

      {/* Jitter Parameters */}
      <Section 
        title="Jitter Parameters" 
        sectionKey="jitter"
        tooltip="Configure spatial uncertainty in sampling"
      >
        <div className="space-y-3">
          <div className="grid grid-cols-3 gap-2">
            <div>
              <label className="block text-xs text-gray-600 mb-1">Sampler Jitter (m)</label>
              <input
                type="number"
                step="0.1"
                value={params.jitter?.sampler || 0.75}
                onChange={(e) => handleNestedChange('jitter', 'sampler', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">Core Jitter (m)</label>
              <input
                type="number"
                step="0.01"
                value={params.jitter?.core || 0.1}
                onChange={(e) => handleNestedChange('jitter', 'core', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">Planned Jitter (m)</label>
              <input
                type="number"
                step="0.1"
                value={params.jitter?.planned || 5.0}
                onChange={(e) => handleNestedChange('jitter', 'planned', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
          </div>
        </div>
      </Section>

      {/* Covariance Models */}
      <Section 
        title="Spatial Covariance" 
        sectionKey="covariance"
        tooltip="Configure spatial correlation structure"
      >
        <div className="space-y-3">
          <div>
            <label className="block text-xs text-gray-600 mb-1">Application Covariance Type</label>
            <select
              value={params.covariance?.applicationType || "gaussian"}
              onChange={(e) => handleNestedChange('covariance', 'applicationType', e.target.value)}
              className="input-field"
            >
              <option value="gaussian">Gaussian</option>
              <option value="spherical">Spherical</option>
              <option value="exponential">Exponential</option>
            </select>
          </div>
          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="block text-xs text-gray-600 mb-1">App Range X (m)</label>
              <input
                type="number"
                value={params.covariance?.appRangeX || 5.0}
                onChange={(e) => handleNestedChange('covariance', 'appRangeX', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">App Range Y (m)</label>
              <input
                type="number"
                value={params.covariance?.appRangeY || 500.0}
                onChange={(e) => handleNestedChange('covariance', 'appRangeY', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
          </div>
          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="block text-xs text-gray-600 mb-1">Soil Range (m)</label>
              <input
                type="number"
                value={params.covariance?.soilRange || 20.0}
                onChange={(e) => handleNestedChange('covariance', 'soilRange', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">Density Range (m)</label>
              <input
                type="number"
                value={params.covariance?.densityRange || 30.0}
                onChange={(e) => handleNestedChange('covariance', 'densityRange', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
          </div>
        </div>
      </Section>

      {/* Execution Parameters */}
      <Section 
        title="Execution Parameters" 
        sectionKey="execution"
        tooltip="Configure simulation execution settings"
      >
        <div className="space-y-3">
          <div>
            <label className="block text-xs text-gray-600 mb-1">Number of Realizations</label>
            <input
              type="number"
              value={params.execution?.numRealizations || 10000}
              onChange={(e) => handleNestedChange('execution', 'numRealizations', parseInt(e.target.value))}
              className="input-field"
              min="1"
              max="50000"
            />
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Random Seed</label>
            <input
              type="number"
              value={params.execution?.randomSeed || 1}
              onChange={(e) => handleNestedChange('execution', 'randomSeed', parseInt(e.target.value))}
              className="input-field"
            />
          </div>
          <div className="grid grid-cols-2 gap-2">
            <div>
              <label className="block text-xs text-gray-600 mb-1">Ca Measurement Error</label>
              <input
                type="number"
                step="0.001"
                value={params.execution?.caMeasurementError || 0.03}
                onChange={(e) => handleNestedChange('execution', 'caMeasurementError', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
            <div>
              <label className="block text-xs text-gray-600 mb-1">Mg Measurement Error</label>
              <input
                type="number"
                step="0.001"
                value={params.execution?.mgMeasurementError || 0.03}
                onChange={(e) => handleNestedChange('execution', 'mgMeasurementError', parseFloat(e.target.value))}
                className="input-field"
              />
            </div>
          </div>
          <div>
            <label className="block text-xs text-gray-600 mb-1">Mass Measurement Error</label>
            <input
              type="number"
              step="0.001"
              value={params.execution?.massMeasurementError || 0.005}
              onChange={(e) => handleNestedChange('execution', 'massMeasurementError', parseFloat(e.target.value))}
              className="input-field"
            />
          </div>
        </div>
      </Section>
    </div>
  )
}

export default SimulationForm