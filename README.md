# Monty Simulation Web App

A modern web application for running Monty geochemical simulations with an intuitive interface.

## Features

- **Interactive Parameter Configuration**: Easy-to-use form for setting simulation parameters
- **Real-time Visualization**: Charts showing concentration changes over time
- **Multiple Realizations**: Support for Monte Carlo simulations with statistical analysis
- **Results Export**: Download simulation results as JSON files
- **Responsive Design**: Works on desktop and mobile devices

## Prerequisites

Before running the web app, make sure you have:

1. **Node.js** (v16 or higher)
2. **Julia** (v1.10 or higher) installed and in your PATH
3. **Monty.jl package** installed in Julia

### Installing Julia and Monty

1. Install Julia from [https://julialang.org/downloads/](https://julialang.org/downloads/)
2. Add Julia to your system PATH
3. Install the Monty package:
   ```julia
   using Pkg
   Pkg.add("https://github.com/LithosCarbon/Monty.jl")
   ```

## Installation

1. Install dependencies:
   ```bash
   npm install
   ```

2. Install server dependencies:
   ```bash
   cd server
   npm install
   cd ..
   ```

## Running the Application

Start both the client and server:

```bash
npm run dev
```

This will start:
- Frontend development server on `http://localhost:3000`
- Backend API server on `http://localhost:3001`

## Usage

1. **Configure Parameters**: Use the form on the left to set simulation parameters:
   - Feedstock properties (application rate, Ca/Mg concentrations)
   - Soil properties (baseline concentrations, sampling depth)
   - Leaching rates for different elements
   - Sampling parameters (number of samples, realizations, time points)

2. **Run Simulation**: Click "Run Simulation" to execute the Monty simulation

3. **View Results**: Explore the results in different tabs:
   - **Concentrations Over Time**: Line charts showing how Ca and Mg concentrations change
   - **Summary Statistics**: Statistical summary of all simulation results
   - **Simulation Parameters**: Review the parameters used for the simulation

4. **Export Results**: Download simulation results as JSON files for further analysis

## Simulation Details

The web app runs a simplified version of the Monty simulation that:

- Creates a circular field geometry for random sampling
- Uses exponential leaching models for Ca and Mg
- Applies measurement noise to simulate real analytical uncertainty
- Runs multiple realizations to capture variability
- Provides statistical summaries of results

## API Endpoints

- `GET /api/status` - Check server and Julia availability
- `POST /api/simulate` - Run a Monty simulation with provided parameters

## Troubleshooting

### "Julia is not available" Error
- Ensure Julia is installed and in your system PATH
- Test by running `julia --version` in your terminal

### Simulation Fails
- Check that the Monty.jl package is properly installed
- Verify all parameter values are reasonable (positive numbers, etc.)
- Check the browser console and server logs for detailed error messages

### Long Simulation Times
- Reduce the number of realizations for faster results
- Decrease the number of samples if needed
- Consider the complexity of your parameter settings

## Development

The application is built with:
- **Frontend**: React + Vite + Tailwind CSS
- **Backend**: Node.js + Express
- **Visualization**: Recharts
- **Simulation Engine**: Julia + Monty.jl

To modify the simulation logic, edit the Julia script generation in `server/index.js`.