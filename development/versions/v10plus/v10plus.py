
import openmc

H5M = "/root/v10plus/somup.h5m"

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
clad.add_element("Zr", 0.99)
clad.add_element("Nb", 0.01)
clad.temperature = 600

d20 = openmc.Material(name="d20")
d20.set_density("g/cm3", 1.1)
d20.add_nuclide("H2", 2.0)
d20.add_nuclide("O16", 1.0)
d20.add_s_alpha_beta("c_D_in_D2O")
d20.temperature = 294

# Fixed-source surrogate material used for transport testing
fuel = openmc.Material(name="fuel")
fuel.set_density("g/cm3", 20.0)
fuel.add_nuclide("O16", 2.0)
fuel.add_nuclide("U235",0.0073)
fuel.add_nuclide("U238",0.9927)
fuel.temperature = 1200

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
pressure.add_element("Zr", 0.975)
pressure.add_element("Nb", 0.025)
pressure.temperature = 600

# Confirmed: ref = graphite
ref = openmc.Material(name="ref")
ref.set_density("g/cm3", 2.0)
ref.add_element("C", 1.0)
ref.add_s_alpha_beta("c_Graphite")
ref.temperature = 294

materials = openmc.Materials([
    air,
    carrier,
    clad,
    d20,
    fuel,
    graphite,
    h20,
    pressure,
    ref
])

# ==============================================================================
# 2. DAGMC GEOMETRY
# ==============================================================================

dag = openmc.DAGMCUniverse(
    H5M,
    auto_geom_ids=True
)

print("\nDAGMC material tags:")
for name in dag.material_names:
    print(" ", name)

expected = {
    "air",
    "carrier",
    "clad",
    "d20",
    "fuel",
    "graphite",
    "h20",
    "pressure",
    "ref"
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
# 3. REFLECTIVE BOUNDARY SETUP
# ==============================================================================

min_x = openmc.XPlane(
    x0=lower_left[0],
    boundary_type="reflective"
)

max_x = openmc.XPlane(
    x0=upper_right[0],
    boundary_type="reflective"
)

min_y = openmc.YPlane(
    y0=lower_left[1],
    boundary_type="reflective"
)

max_y = openmc.YPlane(
    y0=upper_right[1],
    boundary_type="reflective"
)

min_z = openmc.ZPlane(
    z0=lower_left[2],
    boundary_type="reflective"
)

max_z = openmc.ZPlane(
    z0=upper_right[2],
    boundary_type="reflective"
)

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
settings.batches = 30
settings.particles = 50000
settings.inactive = 0
settings.temperature = {"method": "nearest"}

source = openmc.IndependentSource()
source.space = openmc.stats.Point(center)
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete(
    [2.0e6],
    [1.0]
)

settings.source = source

# ==============================================================================
# 5. FLUX TALLY
# ==============================================================================

mesh = openmc.RegularMesh()
mesh.dimension = (20, 20, 20)
mesh.lower_left = lower_left
mesh.upper_right = upper_right

mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name="V10plus_flux")
tally.filters = [mesh_filter]
tally.scores = ["flux"]

tallies = openmc.Tallies([tally])

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
print("          V10PLUS TRANSPORT TEST")
print("========================================")
print("H5M        :", H5M)
print("Run mode   :", settings.run_mode)
print("Batches    :", settings.batches)
print("Inactive   :", settings.inactive)
print("Particles  :", settings.particles)
print("Mesh       :", mesh.dimension)
print("========================================\n")

model.run()
