# Monty Simulation Web Interface Guide

## Overview

This web interface provides a comprehensive front-end for the Monty geochemical simulation package. It allows you to configure all the sophisticated parameters from the main simulation script through an intuitive web UI.

## Getting Started

1. Start the development server:
   ```bash
   npm run dev
   ```

2. Open your browser to `http://localhost:5173`

3. Configure your simulation parameters in the collapsible sections

4. Click "Run Simulation" to execute

## Parameter Sections

### Grid Configuration
- **Grid X/Y Size**: Number of cells in each direction (default: 8x8)
- **Cell Width/Height**: Physical size of each grid cell in meters (default: 10m x 10m)

### Sampling Configuration
- **Cores per Sample**: Number of cores collected per composite sample (default: 5)
- **Sampling Times**: Comma-separated list of sampling times in years (default: -0.003, 0.0, 1.0)
- **Core Stencil Type**: Spatial arrangement of cores (Circle, Hub-Spoke, Line, Random)
- **Stencil Radius/Size**: Physical size parameter for core arrangement (default: 2.0m)

### Feedstock Parameters
- **Application Rate**: Mean application rate in kg/m² (default: 3.5)
- **Feedstock Density**: Bulk density of feedstock in kg/m³ (default: 1000)
- **Ca/Mg Concentration**: Element concentrations in the feedstock (default: 0.07, 0.05)
- **Concentration Variability**: Relative standard deviation for concentrations (default: 0.03)

### Soil Parameters
- **Soil Density**: Mean soil bulk density in kg/m³ (default: 1000)
- **Background Ca/Mg**: Background soil concentrations in kg/kg (default: 0.002, 0.001)
- **Cross-correlation**: Spatial correlation between Ca and Mg (default: 0.75)

### Mixing Model
- **Mixing Type**: How feedstock integrates with soil
  - Unmixed: Layered application (default)
  - Triangular: Wedge-shaped mixing profile
  - Uniform: Uniform mixing over depth
  - Exponential: Exponentially decaying mixing
- **Min/Max Depth**: Range of mixing depths in meters (default: 0.05-0.15m)

### Leaching Models
- **Ca/Mg Leaching Rate**: Exponential decay rates (λ) for each element (default: 0.4, 0.8)
- **Leaching Model Type**: Mathematical model for leaching behavior (default: Exponential)

### Jitter Parameters
- **Sampler Jitter**: Location uncertainty for sampler positioning (default: 0.75m)
- **Core Jitter**: Location uncertainty for individual cores (default: 0.1m)
- **Planned Jitter**: Intentional offset within grid cells (default: 5.0m)

### Spatial Covariance
- **Application Covariance Type**: Spatial correlation model for application rates
- **App Range X/Y**: Correlation ranges for anisotropic application patterns (default: 5m, 500m)
- **Soil Range**: Correlation range for soil properties (default: 20m)
- **Density Range**: Correlation range for soil density (default: 30m)

### Execution Parameters
- **Number of Realizations**: Monte Carlo simulation count (default: 10,000)
- **Random Seed**: For reproducible results (default: 1)
- **Measurement Errors**: Analytical uncertainty for Ca, Mg, and mass measurements

## Output

The simulation generates:
- Time series data for Ca and Mg concentrations
- Mass measurements with uncertainty
- Spatial coordinates for each sample
- Summary statistics across all realizations
- Downloadable JSON results

## Requirements

- Julia with Monty.jl package installed
- Node.js and npm
- All dependencies listed in package.json

## Tips

1. Start with fewer realizations (100-1000) for testing
2. Expand parameter sections as needed - some are collapsed by default
3. Use the default values as a starting point - they match the reference simulation
4. Download results for further analysis in your preferred tools

## Troubleshooting

- **Julia not found**: Ensure Julia is installed and in your system PATH
- **Simulation fails**: Check parameter ranges and ensure positive values where required
- **Long run times**: Reduce the number of realizations for faster testing