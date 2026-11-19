import openmc

sp = openmc.StatePoint('statepoint.35.h5')
tally = sp.get_tally(name='flux_3d')
mesh = tally.find_filter(openmc.MeshFilter).mesh

# Extract the raw 1D flux and power data
flux_data = tally.get_slice(scores=['flux']).mean.flatten()
power_data = sp.get_tally(name='power_3d').get_slice(scores=['kappa-fission']).mean.flatten()

# Export straight to VTK
mesh.write_data_to_vtk(
    filename='somuv4_results.vtk',
    datasets={'flux': flux_data, 'power': power_data}
)
print("Exported to 'somuv4_results.vtk'. Open this file in ParaView!")
