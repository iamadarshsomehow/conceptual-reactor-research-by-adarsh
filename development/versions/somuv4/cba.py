import openmc
import numpy as np

# Load statepoint file
sp = openmc.StatePoint('statepoint.35.h5')

# Extract axial power tally (1 x 1 x 50 mesh)
tally = sp.get_tally(name='axial_profile')
mesh = tally.find_filter(openmc.MeshFilter).mesh

# Axial bounds and mesh geometry
z_min, z_max = mesh.lower_left[2], mesh.upper_right[2]
n_slices = mesh.dimension[2]  # 50
dz = (z_max - z_min) / n_slices  # cm

# Extract raw kappa-fission score (eV per source neutron)
power_raw = tally.get_slice(scores=['kappa-fission']).mean.flatten()

# Calculate Z-centroids for each slice (in meters for OpenFOAM)
z_centroids_m = np.linspace(z_min + dz/2, z_max - dz/2, n_slices) / 100.0

# Normalize to relative power shape (Normalized Peaking Factor F_z)
# Mean of F_z across the height will equal 1.0
f_z = power_raw / np.mean(power_raw)

# Save to CSV for plotting / records
data = np.column_stack((z_centroids_m, power_raw, f_z))
np.savetxt(
    'axial_power_profile.csv', 
    data, 
    delimiter=',', 
    header='z_center_m,kappa_fission_eV,normalized_Fz', 
    comments=''
)

# Generate OpenFOAM NonUniform List format for fvOptions
of_list = "(\n" + "\n".join([f"    ({z:.5f} {val:.5f})" for z, val in zip(z_centroids_m, f_z)]) + "\n);"

with open('openfoam_power_input.txt', 'w') as f:
    f.write(of_list)

print("========================================")
print("     AXIAL PROFILE EXPORT COMPLETE      ")
print("========================================")
print(f"Extracted {n_slices} axial slices across Z = [{z_min:.2f}, {z_max:.2f}] cm")
print("Saved outputs:")
print("  1. axial_power_profile.csv  (For Excel / Python plotting)")
print("  2. openfoam_power_input.txt (For OpenFOAM fvOptions)")
print("========================================")

