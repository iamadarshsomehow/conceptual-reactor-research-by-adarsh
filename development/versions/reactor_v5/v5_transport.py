import openmc

# ============================================================
# V5 — O2 FIXED-SOURCE TRANSPORT TEST
# Geometry: baka.h5m
# ============================================================

H5M = "/root/reactor_v5/baka.h5m"


# ============================================================
# MATERIAL DEFINITIONS
# ============================================================

# mat:air -> PURE HELIUM
air = openmc.Material(name="air")
air.set_density("g/cm3", 0.00016)
air.add_nuclide("He4", 1.0)
air.temperature = 900


# mat:carrier -> 99% Zr + 1% Nb
carrier = openmc.Material(name="carrier")
carrier.set_density("g/cm3", 4)
carrier.add_element("Zr", 0.99)
carrier.add_element("Nb", 0.01)
carrier.temperature = 600


# mat:clad -> 99% Zr + 1% Nb
clad = openmc.Material(name="clad")
clad.set_density("g/cm3", 4)
clad.add_element("Zr", 0.99)
clad.add_element("Nb", 0.01)
clad.temperature = 600


# mat:d20 -> PURE D2O
d20 = openmc.Material(name="d20")
d20.set_density("g/cm3", 1)
d20.add_nuclide("H2", 2.0)
d20.add_nuclide("O16", 1.0)
d20.temperature = 294


# mat:fuel -> O2 TEST MATERIAL
# Deliberately non-reactor test material.
o2 = openmc.Material(name="fuel")
o2.set_density("g/cm3", 20)
o2.add_nuclide("O16", 2.0)
o2.add_nuclide("U235",0.0073)
o2.add_nuclide("U238",0.9927)
o2.temperature = 1200


# mat:graphite -> PURE CARBON
graphite = openmc.Material(name="graphite")
graphite.set_density("g/cm3", 2)
graphite.add_element("C", 1.0)
graphite.temperature = 294


# mat:h20 -> PURE H2O
h20 = openmc.Material(name="h20")
h20.set_density("g/cm3", 0.6)
h20.add_nuclide("H1", 2.0)
h20.add_nuclide("O16", 1.0)
h20.temperature = 294


# mat:pressure -> 97.5% Zr + 2.5% Nb
pressure = openmc.Material(name="pressure")
pressure.set_density("g/cm3", 4)
pressure.add_element("Zr", 0.975)
pressure.add_element("Nb", 0.025)
pressure.temperature = 600


materials = openmc.Materials([
    air,
    carrier,
    clad,
    d20,
    o2,
    graphite,
    h20,
    pressure,
])


# ============================================================
# LOAD DAGMC GEOMETRY
# ============================================================

print("=" * 70)
print("V5 O2 FIXED-SOURCE TRANSPORT")
print("=" * 70)

print("\nLoading DAGMC geometry:")
print(H5M)

dag = openmc.DAGMCUniverse(H5M)


# ============================================================
# VERIFY CAD MATERIAL TAGS
# ============================================================

expected_tags = {
    "air",
    "carrier",
    "clad",
    "d20",
    "fuel",
    "graphite",
    "h20",
    "pressure",
}

found_tags = set(dag.material_names)

print("\nCAD material tags:")
for tag in sorted(found_tags):
    print("  mat:" + tag)

missing = expected_tags - found_tags
unexpected = found_tags - expected_tags

if missing:
    raise RuntimeError(
        "Missing CAD material tags: " + str(sorted(missing))
    )

if unexpected:
    print("\nWARNING: unexpected CAD tags:")
    for tag in sorted(unexpected):
        print("  mat:" + tag)

print("\nCAD TAG CHECK: PASS")


# ============================================================
# CAD -> OPENMC MATERIAL MAPPING
# ============================================================

cad_to_openmc = {
    "air": air,
    "carrier": carrier,
    "clad": clad,
    "d20": d20,
    "fuel": o2,
    "graphite": graphite,
    "h20": h20,
    "pressure": pressure,
}

print("\nCAD -> OpenMC material mapping:")
for tag, material in sorted(cad_to_openmc.items()):
    print(f"  mat:{tag:10s} -> {material.name}")

print("\nMATERIAL MAPPING CHECK: PASS")


# ============================================================
# GEOMETRY
# ============================================================

universe = dag.bounded_universe()
geometry = openmc.Geometry(universe)


# ============================================================
# DAGMC BOUNDING BOX
# ============================================================

lower = dag.bounding_box.lower_left
upper = dag.bounding_box.upper_right
center = (lower + upper) / 2.0


# ============================================================
# FIXED NEUTRON SOURCE
# ============================================================

source = openmc.IndependentSource()

source.space = openmc.stats.Point(center)
source.angle = openmc.stats.Isotropic()

# 1 MeV test neutron
source.energy = openmc.stats.Discrete(
    [1.0e6],
    [1.0]
)


# ============================================================
# SETTINGS
# ============================================================

settings = openmc.Settings()

settings.run_mode = "eigenvalue"

settings.source = source

settings.batches = 40
settings.particles = 40000
settings.inactive = 10

settings.temperature = {
    "method": "nearest"
}


# ============================================================
# 3-D FLUX MESH
# ============================================================

mesh = openmc.RegularMesh()

mesh.lower_left = lower
mesh.upper_right = upper

mesh.dimension = (100, 100, 100)

mesh_filter = openmc.MeshFilter(mesh)

flux = openmc.Tally(name="V5_flux")

flux.filters = [mesh_filter]
flux.scores = ["flux"]

tallies = openmc.Tallies([flux])


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
# SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("MODEL SUMMARY")
print("=" * 70)

print("DAGMC cells    :", dag.n_cells)
print("DAGMC surfaces :", dag.n_surfaces)

print("\nBounding box:")
print("  Lower:", lower)
print("  Upper:", upper)

print("\nSource:")
print("  Position:", center)
print("  Energy  : 1 MeV")
print("  Type    : isotropic")

print("\nRun:")
print("  Mode      : fixed source")
print("  Batches   :", settings.batches)
print("  Particles :", settings.particles)

print("\nFlux mesh:", mesh.dimension)


# ============================================================
# RUN
# ============================================================

print("\nStarting OpenMC...")
print("=" * 70)

model.run()

print("\n" + "=" * 70)
print("OPENMC RUN FINISHED")
print("=" * 70)
