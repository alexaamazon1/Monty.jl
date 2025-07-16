# 🚀 Run Monty Simulation Locally

## Quick Start (3 Steps)

### Option 1: Simple Static Version (No Server Required)
1. Open `index.html` directly in your browser
2. This gives you the full UI to explore parameters
3. Note: Simulations won't run without the backend server

### Option 2: Full Setup with Simulation Support

#### Step 1: Install Node.js
- Download from: https://nodejs.org/
- Choose LTS version

#### Step 2: Install Dependencies
```bash
npm install
```

#### Step 3: Start the Application
```bash
# Start both frontend and backend
npm run dev
```

Then open: http://localhost:5173

## 🎯 What You Get

### Complete Monty Simulation Interface
- **Grid Configuration**: 8×8 field layout
- **Sampling Design**: Core arrangements and timing
- **Feedstock Parameters**: Application rates and chemistry
- **Soil Properties**: Background concentrations and density
- **Mixing Models**: 4 different mixing approaches
- **Leaching Behavior**: Element-specific decay rates
- **Spatial Correlation**: Advanced covariance models
- **Uncertainty Modeling**: Jitter and measurement errors

### Default Parameters (from main.jl)
All defaults match your reference simulation:
- Application rate: 3.5 kg/m²
- Grid: 8×8 cells, 10m each
- 5 cores per sample in 2m circle
- Ca leaching λ=0.4, Mg λ=0.8
- 10,000 Monte Carlo realizations

## 🔧 For Full Simulations

### Install Julia
1. Download from: https://julialang.org/downloads/
2. Add to system PATH
3. Install Monty package:
   ```julia
   using Pkg
   Pkg.add("Monty")
   ```

### Verify Setup
- Backend API: http://localhost:3001/api/status
- Should show: `{"julia":true,"server":"running"}`

## 📊 Using the Interface

1. **Explore Sections**: Click to expand parameter groups
2. **Adjust Parameters**: Modify values or use proven defaults
3. **Run Simulation**: Click "Run Simulation" button
4. **Download Results**: Get JSON data for analysis

## 🛠️ Troubleshooting

### Port Already in Use
```bash
# Kill existing processes
pkill -f "vite|node"
npm run dev
```

### Julia Not Found
- Ensure Julia is installed and in PATH
- Test: `julia --version`

### Simulation Errors
- Start with fewer realizations (100-1000)
- Check parameter ranges (positive values)
- Verify Julia package installation

## 📁 Files Included

- `index.html` - Static UI (no server needed)
- `src/` - React source code
- `server/` - Node.js backend
- `package.json` - Dependencies
- This README

## 🎉 Ready to Simulate!

The interface provides the full power of your sophisticated Monty simulation with an intuitive web UI. All the complex parameters from your Julia script are now accessible through organized, collapsible sections.

Perfect for:
- Testing different scenarios
- Educational demonstrations  
- Rapid parameter exploration
- Production simulations