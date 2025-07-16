import axios from 'axios'

const API_BASE_URL = '/api'

export const runSimulation = async (params) => {
  try {
    const response = await axios.post(`${API_BASE_URL}/simulate`, params, {
      timeout: 300000, // 5 minutes timeout for long simulations
    })
    return response.data
  } catch (error) {
    if (error.response) {
      throw new Error(error.response.data.error || 'Simulation failed')
    } else if (error.request) {
      throw new Error('No response from server. Please check if Julia and Monty are installed.')
    } else {
      throw new Error('Request failed: ' + error.message)
    }
  }
}

export const getSimulationStatus = async () => {
  try {
    const response = await axios.get(`${API_BASE_URL}/status`)
    return response.data
  } catch (error) {
    throw new Error('Failed to get simulation status')
  }
}