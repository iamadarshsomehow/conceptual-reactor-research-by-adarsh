import openmc

H5M = "/root/reactor_v9/autism.h5m"

# ============================================================
# MATERIALS
# ============================================================

air = openmc.Material(name="air")
air.set_density("g/cm3", 0.00016)
air.add_nuclide("He4", 1.0)
air.temperature = 900

carrier = openmc.Material(name="carrier")
carrier.set_density("g/cm3", 2.0)
carrier.add_element("C", 1.0)
carrier.temperature = 294

clad = openmc.Material(name="clad")
clad.set_density("g/cm3", 4.0)
clad.add_element("Zr", 0.99)
clad.add_element("Nb", 0.01)
clad.temperature = 600

d20 = openmc.Material(name="d20")
d20.set_density("g/cm3", 1.0)
d20.add_nuclide("H2", 2.0)
d20.add_nuclide("O16", 1.0)
d20.temperature = 294

# Safe O2 surrogate for transport testing
fuel = openmc.Material(name="fuel")
fuel.set_density("g/cm3", 10.5)
fuel.add_nuclide("O16", 2.0)
fuel.add_nuclide("U235",0.0073)
fuel.add_nuclide("U238",0.9927)
fuel.temperature = 1200

graphite = openmc.Material(name="graphite")
graphite.set_density("g/cm3", 2.0)
graphite.add_element("C", 1.0)
graphite.temperature = 294

h20 = openmc.Material(name="h20")
h20.set_density("g/cm3", 0.6)
h20.add_nuclide("H1", 2.0)
h20.add_nuclide("O16", 1.0)
h20.temperature = 294

pressure = openmc.Material(name="pressure")
pressure.set_density("g/cm3", 4.0)
pressure.add_element("Zr", 0.975)
pressure.add_element("Nb", 0.025)
pressure.temperature = 600

materials = openmc.Materials([
    air,
    carrier,
    clad,
    d20,
    fuel,
    graphite,
    h20,
    pressure
])

# ============================================================
# DAGMC GEOMETRY
# ============================================================

dag = openmc.DAGMCUniverse(H5M)

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
    "pressure"
}

missing = expected - set(dag.material_names)

if missing:
    raise RuntimeError(
        f"Missing expected DAGMC tags: {sorted(missing)}"
    )

universe = dag.bounded_universe()
geometry = openmc.Geometry(universe)

lower_left = dag.bounding_box[0]
upper_right = dag.bounding_box[1]
center = 0.5 * (lower_left + upper_right)

print("\nBounding box:")
print(" Lower:", lower_left)
print(" Upper:", upper_right)
print(" Source:", center)

# ============================================================
# FIXED-SOURCE SETTINGS
# ============================================================

settings = openmc.Settings()
settings.run_mode = "eigenvalue"
settings.batches = 30
settings.particles = 50000
settings.inactive = 10

settings.temperature = {
    "method": "nearest"
}

source = openmc.Source()
source.space = openmc.stats.Point(center)
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete(
    [1.0e6],
    [1.0]
)

settings.source = source

# ============================================================
# 100 x 100 x 100 FLUX MESH
# ============================================================

mesh = openmc.RegularMesh()
mesh.dimension = (100, 100, 100)
mesh.lower_left = lower_left
mesh.upper_right = upper_right

mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name="V9_flux")
tally.filters = [mesh_filter]
tally.scores = ["flux"]

tallies = openmc.Tallies([tally])

# ============================================================
# MODEL
# ============================================================

model = openmc.model.Model(
    geometry=geometry,
    materials=materials,
    settings=settings,
    tallies=tallies
)

print("\n========================================")
print("          V9 TRANSPORT TEST")
print("========================================")
print("H5M       :", H5M)
print("Run mode  :", settings.run_mode)
print("Batches   :", settings.batches)
print("Particles :", settings.particles)
print("Mesh      :", mesh.dimension)
print("========================================\n")

model.run()
