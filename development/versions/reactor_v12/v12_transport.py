import openmc

H5M = "/root/reactor_v12/hai.h5m"

# =========================
# V12 MATERIAL DEFINITIONS
# =========================

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
d20.set_density("g/cm3", 1.11)
d20.add_nuclide("H2", 2.0)
d20.add_nuclide("O16", 1.0)
d20.temperature = 294

# Non-fissile surrogate for transport testing
fuel = openmc.Material(name="fuel")
fuel.set_density("g/cm3", 20.0)
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

# =========================
# DAGMC GEOMETRY
# =========================

dag = openmc.DAGMCUniverse(H5M)

expected_tags = {
    "air",
    "carrier",
    "clad",
    "d20",
    "fuel",
    "graphite",
    "h20",
    "pressure"
}

actual_tags = set(dag.material_names)

missing = expected_tags - actual_tags

if missing:
    raise RuntimeError(
        f"Missing expected DAGMC tags: {sorted(missing)}"
    )

universe = dag.bounded_universe()
geometry = openmc.Geometry(universe)

# =========================
# V12 FIXED-SOURCE SETTINGS
# =========================

settings = openmc.Settings()

settings.run_mode = "eigenvalue"

settings.batches = 40
settings.inactive = 10
settings.particles = 50000

settings.temperature = {
    "method": "nearest"
}

# =========================
# SOURCE
# =========================

lower_left = dag.bounding_box[0]
upper_right = dag.bounding_box[1]

center = 0.5 * (lower_left + upper_right)

source = openmc.Source()
source.space = openmc.stats.Point(center)
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete(
    [1.0e6],
    [1.0]
)

settings.source = source

# =========================
# 100 x 100 x 100 MESH
# =========================

mesh = openmc.RegularMesh()

mesh.dimension = (100, 100, 100)
mesh.lower_left = lower_left
mesh.upper_right = upper_right

mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name="V12_flux")

tally.filters = [
    mesh_filter
]

tally.scores = [
    "flux"
]

tallies = openmc.Tallies([
    tally
])

# =========================
# BUILD AND RUN
# =========================

model = openmc.model.Model(
    geometry=geometry,
    materials=materials,
    settings=settings,
    tallies=tallies
)

print("================================")
print("V12 TRANSPORT TEST")
print("================================")
print("D2O density       : 1.11 g/cm3")
print("O2 surrogate      : 20.0 g/cm3")
print("Inactive batches  : 10")
print("Active batches    : 30")
print("Particles/batch   : 50000")
print("Mesh              : 100 x 100 x 100")
print("================================")

model.run()

print("V12 transport finished.")
