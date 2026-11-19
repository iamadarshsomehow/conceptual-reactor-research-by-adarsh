import openmc

H5M = "/root/reactor_v8_5/dom.h5m"

# ============================================================
# MATERIALS
# ============================================================

air = openmc.Material(name="air")
air.set_density("g/cm3", 0.00016)
air.add_nuclide("He4", 1.0)
air.temperature = 900

# carrier = GRAPHITE
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

# O2 surrogate — NOT uranium
fuel = openmc.Material(name="fuel")
fuel.set_density("g/cm3", 10.5)
fuel.add_nuclide("O16", 2.0)
fuel.add_nuclide("U235",0.0073)
fuel.add_nuclide("U238",0.9927)
fuel.temperature = 1200

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
    carrier,
    clad,
    d20,
    fuel,
    h20,
    pressure
])

# ============================================================
# DAGMC
# ============================================================

dag = openmc.DAGMCUniverse(H5M)

print("\nDAGMC material tags:")
for name in dag.material_names:
    print(" ", name)

expected = {
    "carrier",
    "clad",
    "d20",
    "fuel",
    "h20",
    "pressure"
}

missing = expected - set(dag.material_names)

if missing:
    raise RuntimeError(
        f"Missing expected DAGMC tags: {sorted(missing)}"
    )

# ============================================================
# GEOMETRY
# ============================================================

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
# SETTINGS
# ============================================================

settings = openmc.Settings()

settings.run_mode = "eigenvalue"
settings.batches = 20
settings.particles = 50000

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
# FLUX MESH
# ============================================================

mesh = openmc.RegularMesh()
mesh.dimension = (100, 100, 100)
mesh.lower_left = lower_left
mesh.upper_right = upper_right

mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name="V85_flux")
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
print("       V8.5 TRANSPORT TEST")
print("========================================")
print("H5M       :", H5M)
print("Run mode  :", settings.run_mode)
print("Batches   :", settings.batches)
print("Particles :", settings.particles)
print("Mesh      :", mesh.dimension)
print("========================================\n")

model.run()
