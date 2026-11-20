"""Archived exploratory residual-analysis example; not used in the final workflow.

This script expects a matching statepoint and named mesh tallies. Its example
cross-section constants are illustrative and were not used for reported results.
"""
import numpy as np
import openmc

# 1. Load your simulation results
sp = openmc.StatePoint('statepoint.35.h5')

# 2. Extract your physical tallies (Assumes you set up a Mesh Tally for Flux and Fission)
flux_tally = sp.get_tally(name='flux_mesh')
fission_tally = sp.get_tally(name='fission_mesh')

# Get the raw mean matrices (reshaped to fit your 3D mesh dimensions)
# Example assuming a 10x10x10 spatial mesh
nx, ny, nz = 10, 10, 10
flux = flux_tally.get_values(value='mean').reshape((nx, ny, nz))
fission = fission_tally.get_values(value='mean').reshape((nx, ny, nz))

# 3. Illustrative macroscopic cross-section constants for this archived example.
# They are not material data for the final design or used in reported results.
Sigma_t = 0.35   # Total Macroscopic Cross Section (cm^-1)
nu_Sigma_f = 0.08 # Fission Production Cross Section (cm^-1)
k_eff = sp.k_combined.nominal_value # The converged multiplication factor

# 4. Calculate the Residual Matrix (r = Production / k - Destruction)
# This implements the standard eigenvalue residual: (1/k)*F*phi - A*phi
production = (nu_Sigma_f * fission) / k_eff
destruction = Sigma_t * flux

# The raw numerical residual for every single mesh cell
residuals = production - destruction

# 5. Export the raw numbers to a CSV / text file to hold them "in hand"
flattened_residuals = residuals.flatten()
np.savetxt("openmc_residuals.csv", flattened_residuals, delimiter=",", header="Cell_Residual")

print(f" Residual calculation complete. Max residual value: {np.max(np.abs(flattened_residuals)):.6e}")
print("Raw data exported to 'openmc_residuals.csv'")
