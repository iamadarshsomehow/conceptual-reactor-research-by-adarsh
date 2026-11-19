import openmc
import numpy as np

# ============================================================
# V8 — BKL
# Fixed-source transport test
# ============================================================

H5M = "/root/reactor_v8/bkl.h5m"

# ============================================================
# MATERIALS
# ============================================================

# air
air = openmc.Material(name="air")
air.set_density("g/cm3", 0.00016)
air.add_nuclide("He4", 1.0)
air.temperature = 900

# mat:carrier -> 99% Zr + 1% Nb
carrier = openmc.Material(name="carrier")
carrier.set_density("g/cm3", 2)
carrier.add_element("C", 1.0)
carrier.temperature = 294

# mat:clad -> 99% Zr + 1% Nb
clad = openmc.Material(name="clad")
clad.set_density("g/cm3", 4)
clad.add_element("Zr", 0.99)
clad.add_element("Nb", 0.01)
clad.temperature = 600

# mat:d20 -> PURE D2O
d20 = openmc.Material(name="d20")
d20.set_density("g/cm3", 1.10)
d20.add_nuclide("H2", 2.0)
d20.add_nuclide("O16", 1.0)
d20.temperature = 294

# mat:fuel -> O2 SURROGATE
fuel = openmc.Material(name="fuel")
fuel.set_density("g/cm3", 10.5)
fuel.add_nuclide("O16", 2.0)
fuel.add_nuclide("U238",0.9927)
fuel.add_nuclide("U235",0.0073)
fuel.temperature = 1200

# mat:graphite -> PURE CARBON
graphite = openmc.Material(name="graphite")
graphite.set_density("g/cm3", 2)
graphite.add_element("C", 1.0)
graphite.temperature = 294

# mat:h20 -> PURE H2O
h20 = openmc.Material(name="h20")
h20.set_density("g/cm3", 0.5)
h20.add_nuclide("H1", 2.0)
h20.add_nuclide("O16", 1.0)
h20.temperature = 600

# mat:pressure -> 97.5% Zr + 2.5% Nb
pressure = openmc.Material(name="pressure")
pressure.set_density("g/cm3", 4)
pressure.add_element("Zr", 0.975)
pressure.add_element("Nb", 0.025)
pressure.temperature = 600

# ref -> BeO reflector
reflector = openmc.Material(name="ref")
reflector.set_density("g/cm3", 2)
reflector.add_element("C", 1.0)
reflector.temperature = 294

# ============================================================
# MATERIAL COLLECTION
# ============================================================

materials = openmc.Materials([
    air,
    carrier,
    clad,
    d20,
    fuel,
    graphite,
    h20,
    pressure,
    reflector,
])

# ============================================================
# DAGMC GEOMETRY
# ============================================================

dag = openmc.DAGMCUniverse(H5M)

print("\nMaterial tags found in DAGMC:")
for name in dag.material_names:
    print(" ", name)

# ============================================================
# MATERIAL TAG -> OPENMC MATERIAL
# ============================================================

material_map = {
    "air": air,
    "carrier": carrier,
    "clad": clad,
    "d20": d20,
    "fuel": fuel,
    "graphite": graphite,
    "h20": h20,
    "pressure": pressure,
    "ref":reflector
}

# Check every DAGMC material tag
for name in dag.material_names:
    if name not in material_map:
        raise RuntimeError(
            f"No OpenMC material mapping found for DAGMC tag: {name}"
        )

print("\nMaterial mapping:")
for name in dag.material_names:
    print(f"  {name} -> {material_map[name].name}")

# ============================================================
# GEOMETRY
# ============================================================

universe = dag.bounded_universe()
geometry = openmc.Geometry(universe)

# ============================================================
# SETTINGS
# ============================================================

settings = openmc.Settings()

settings.run_mode = "eigenvalue"

settings.batches = 50
settings.particles = 50000
settings.inactive = 25
settings.temperature = {
    "method": "nearest"
}

# Source at center of DAGMC bounding box
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

# ============================================================
# 100 x 100 x 100 FLUX MESH
# ============================================================

mesh = openmc.RegularMesh()
mesh.dimension = (100, 100, 100)
mesh.lower_left = lower_left
mesh.upper_right = upper_right

mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name="V8_flux")
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

# ============================================================
# RUN
# ============================================================

print("\n========================================")
print("        V8 — BKL TRANSPORT TEST")
print("========================================")
print("H5M:", H5M)
print("Source:", center)
print("Particles:", settings.particles)
print("Batches:", settings.batches)
print("Mesh:", mesh.dimension)
print("Run mode:", settings.run_mode)
print("========================================\n")

model.run()
