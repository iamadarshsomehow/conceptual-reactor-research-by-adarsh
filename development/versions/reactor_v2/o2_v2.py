import openmc

H5M = "/root/reactor_v2/vo.h5m"

# ---------------- MATERIALS ----------------

air = openmc.Material(name="air")
air.set_density("g/cm3", 0.00016)
air.add_nuclide("He4", 1.0)
air.temperature = 900.0

clad = openmc.Material(name="clad")
clad.set_density("g/cm3", 6.5)
clad.add_element("Zr", 0.99)
clad.add_element("Nb", 0.01)
clad.temperature = 600.0

pressure = openmc.Material(name="pressure")
pressure.set_density("g/cm3", 6.5)
pressure.add_element("Zr", 0.975)
pressure.add_element("Nb", 0.025)
pressure.temperature = 600.0

carrier = openmc.Material(name="Part")
carrier.set_density("g/cm3", 6.5)
carrier.add_element("Zr", 0.99)
carrier.add_element("Nb", 0.01)
carrier.temperature = 600.0

d20 = openmc.Material(name="d20")
d20.set_density("g/cm3", 0.86)
d20.add_nuclide("H2", 2.0)
d20.add_nuclide("O16", 1.0)
d20.temperature = 600.0

h20 = openmc.Material(name="h20")
h20.set_density("g/cm3", 0.74)
h20.add_nuclide("H1", 2.0)
h20.add_nuclide("O16", 1.0)
h20.temperature = 600.0

graphite = openmc.Material(name="graphite")
graphite.set_density("g/cm3", 1.7)
graphite.add_element("C", 1.0)
graphite.temperature = 600.0

# O2 TEST MATERIAL
fuel = openmc.Material(name="fuel")
fuel.set_density("g/cm3", 10.5)
fuel.add_nuclide("O16", 2.0)
fuel.add_nuclide("U235",0.0073)
fuel.add_nuclide("U238",0.9927)
fuel.temperature = 1200.0

materials = openmc.Materials([
    carrier,
    air,
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

dag = openmc.DAGMCUniverse("vo.h5m")

print("DAGMC materials:")
for name in dag.material_names:
    print("  ", name)

geometry = openmc.Geometry(dag.bounded_universe())

# ============================================================
# FIXED-SOURCE SETTINGS
# ============================================================

settings = openmc.Settings()
settings.run_mode = "eigenvalue"
settings.batches = 100
settings.particles = 20000
settings.inactive = 30

settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((0.0, 0.0, 30.0)),
    angle=openmc.stats.Isotropic(),
    energy=openmc.stats.Discrete(
        [1.0e6],
        [1.0]
    )
)

# ============================================================
# TALLY
# ============================================================

mesh = openmc.RegularMesh()
mesh.dimension = (10, 10, 20)
mesh.lower_left = (-17.5, -20.207260131835938, 0.0)
mesh.upper_right = (17.5, 20.207260131835938, 60.0)

mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name="flux")
tally.filters = [mesh_filter]
tally.scores = ["flux"]

tallies = openmc.Tallies([tally])

# ============================================================
# MODEL
# ============================================================

model = openmc.Model(
    geometry=geometry,
    materials=materials,
    settings=settings,
    tallies=tallies
)

# ============================================================
# SYNCHRONIZE DAGMC MATERIALS
# ============================================================

model.init_lib()
model.sync_dagmc_universes()

print("\nDAGMC cells after synchronization:")
for cell_id, cell in dag.cells.items():
    print(
        cell_id,
        "->",
        cell.name,
        "| fill =",
        cell.fill
    )

model.finalize_lib()

# ============================================================
# RUN
# ============================================================

print("\nStarting fixed-source O2 transport calculation...")

model.run()

print("\nSimulation finished.")
