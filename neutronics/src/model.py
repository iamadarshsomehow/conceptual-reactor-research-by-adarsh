import os
from pathlib import Path

import numpy as np
import openmc

H5M = os.environ.get(
    "SOMUV4_H5M", str(Path(__file__).resolve().parents[2] / "private" / "somuv4.h5m")
)

# ==============================================================================
# 1. MATERIALS DEFINITION
# ==============================================================================
air = openmc.Material(name="air")
air.set_density("g/cm3", 0.00016)
air.add_nuclide("He4", 1.0)
air.temperature = 900

carrier = openmc.Material(name="carrier")
carrier.set_density("g/cm3", 2.0)
carrier.add_element("C", 1.0)
carrier.add_s_alpha_beta("c_Graphite")
carrier.temperature = 294

clad = openmc.Material(name="clad")
clad.set_density("g/cm3", 4.0)
clad.add_element("Zr", 1)
clad.temperature = 600

d20 = openmc.Material(name="d20")   # same name, different contents
d20.add_nuclide("N14", 0.78)        # approximate atom fractions
d20.add_nuclide("O16", 0.21)
d20.add_nuclide("Ar40", 0.01)
d20.set_density("g/cm3", 0.0012)
d20.temperature = 294

fuel = openmc.Material(name="fuel")
fuel.set_density("g/cm3", 20)
fuel.add_nuclide("O16", 2.0)
fuel.add_nuclide("U235", 0.0073)
fuel.add_nuclide("U238", 0.9927)
fuel.temperature = 900

graphite = openmc.Material(name="graphite")
graphite.set_density("g/cm3", 2.0)
graphite.add_element("C", 1.0)
graphite.add_s_alpha_beta("c_Graphite")
graphite.temperature = 294

h20 = openmc.Material(name="h20")
h20.set_density("g/cm3", 0.6)
h20.add_nuclide("H1", 2.0)
h20.add_nuclide("O16", 1.0)
h20.add_s_alpha_beta("c_H_in_H2O")
h20.temperature = 294

pressure = openmc.Material(name="pressure")
pressure.set_density("g/cm3", 4.0)
pressure.add_element("Zr", 1)
pressure.temperature = 600

materials = openmc.Materials([
    air, carrier, clad, d20, fuel, graphite, h20, pressure
])

# ==============================================================================
# 2. DAGMC GEOMETRY
# ==============================================================================
dag = openmc.DAGMCUniverse(H5M, auto_geom_ids=True)

print("\nDAGMC material tags:")
for name in dag.material_names:
    print(" ", name)

expected = {
    "air", "carrier", "clad", "d20",
    "fuel", "graphite", "h20", "pressure"
}

missing = expected - set(dag.material_names)
if missing:
    raise RuntimeError(
        f"Missing expected DAGMC tags: {sorted(missing)}"
    )

lower_left, upper_right = dag.bounding_box
center = 0.5 * (lower_left + upper_right)

print("\nBounding box:")
print(" Lower :", lower_left)
print(" Upper :", upper_right)
print(" Source:", center)

# ==============================================================================
# 3. REFLECTIVE BOUNDING REGION
# ==============================================================================
min_x = openmc.XPlane(x0=lower_left[0], boundary_type="reflective")
max_x = openmc.XPlane(x0=upper_right[0], boundary_type="reflective")
min_y = openmc.YPlane(y0=lower_left[1], boundary_type="reflective")
max_y = openmc.YPlane(y0=upper_right[1], boundary_type="reflective")
min_z = openmc.ZPlane(z0=lower_left[2], boundary_type="reflective")
max_z = openmc.ZPlane(z0=upper_right[2], boundary_type="reflective")

bounding_region = (
    +min_x & -max_x &
    +min_y & -max_y &
    +min_z & -max_z
)

root_cell = openmc.Cell(
    cell_id=99999,
    fill=dag,
    region=bounding_region
)

geometry = openmc.Geometry([root_cell])

# ==============================================================================
# 4. SIMULATION SETTINGS
# ==============================================================================
settings = openmc.Settings()
settings.run_mode = "eigenvalue"
settings.batches = 35
settings.particles = 50000
settings.inactive = 5
settings.temperature = {"method": "nearest"}

source = openmc.IndependentSource()
source.space = openmc.stats.Point(center)
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([2.0e6], [1.0])

settings.source = source

# ==============================================================================
# 5. COMPREHENSIVE TALLIES SETUP
# ==============================================================================
# A. 3D Mesh for Spatial Flux & Power Distributions
mesh_3d = openmc.RegularMesh()
mesh_3d.dimension = (50, 50, 50)
mesh_3d.lower_left = lower_left
mesh_3d.upper_right = upper_right
mesh_3d_filter = openmc.MeshFilter(mesh_3d)

# 1. 3D Neutron Flux Distribution Tally
tally_flux_3d = openmc.Tally(name="flux_3d")
tally_flux_3d.filters = [mesh_3d_filter]
tally_flux_3d.scores = ["flux"]

# 2. 3D Fission & Power Distribution Tally (kappa-fission = energy deposition in eV)
tally_power_3d = openmc.Tally(name="power_3d")
tally_power_3d.filters = [mesh_3d_filter]
tally_power_3d.scores = ["fission", "kappa-fission"]

# 3. Axial Power Profile (1D Mesh along Z)
mesh_axial = openmc.RegularMesh()
mesh_axial.dimension = (1, 1, 50)
mesh_axial.lower_left = lower_left
mesh_axial.upper_right = upper_right
mesh_axial_filter = openmc.MeshFilter(mesh_axial)

tally_axial = openmc.Tally(name="axial_profile")
tally_axial.filters = [mesh_axial_filter]
tally_axial.scores = ["kappa-fission", "flux"]

# 4. Radial Power Profile (1D Mesh along X)
mesh_radial = openmc.RegularMesh()
mesh_radial.dimension = (50, 1, 1)
mesh_radial.lower_left = lower_left
mesh_radial.upper_right = upper_right
mesh_radial_filter = openmc.MeshFilter(mesh_radial)

tally_radial = openmc.Tally(name="radial_profile")
tally_radial.filters = [mesh_radial_filter]
tally_radial.scores = ["kappa-fission", "flux"]

# 5. Energy-Dependent Flux Spectrum Tally (200 bins: 1 meV to 10 MeV)
energy_bins = np.logspace(-3, 7, 201)
energy_filter = openmc.EnergyFilter(energy_bins)

tally_spectrum = openmc.Tally(name="energy_spectrum")
tally_spectrum.filters = [energy_filter]
tally_spectrum.scores = ["flux"]

# 6. Global Reaction Rates
tally_rx_rates = openmc.Tally(name="reaction_rates")
tally_rx_rates.scores = ["fission", "nu-fission", "absorption", "(n,gamma)"]

# 7. Surface Current / Boundary Leakage Tally
boundary_surfaces = [min_x, max_x, min_y, max_y, min_z, max_z]
surface_filter = openmc.SurfaceFilter(boundary_surfaces)

tally_current = openmc.Tally(name="boundary_current")
tally_current.filters = [surface_filter]
tally_current.scores = ["current"]

# Register all tallies
tallies = openmc.Tallies([
    tally_flux_3d,
    tally_power_3d,
    tally_axial,
    tally_radial,
    tally_spectrum,
    tally_rx_rates,
    tally_current
])

# ==============================================================================
# 6. MODEL EXECUTION
# ==============================================================================
model = openmc.model.Model(
    geometry=geometry,
    materials=materials,
    settings=settings,
    tallies=tallies
)

print("\n========================================")
print("          SOMUV4 TRANSPORT TEST")
print("========================================")
print("H5M        :", H5M)
print("Run mode   :", settings.run_mode)
print("Batches    :", settings.batches)
print("Particles  :", settings.particles)
print("Mesh       :", mesh_3d.dimension)
print("========================================\n")

model.run()
