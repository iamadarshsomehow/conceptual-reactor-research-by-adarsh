"""
RBMK benchmark cell (DAGMC): k-infinity versus coolant void fraction.
OpenMC, ENDF/B-VII.1. Same method and settings as sweep_all.py.

Usage (nuclear-env venv; CoolProp installed):
  python rbmk_benchmark_sweep.py --quick                     # short test
  python rbmk_benchmark_sweep.py --pressure 7  --enrich 2.0  # RBMK as built
  python rbmk_benchmark_sweep.py --pressure 10 --enrich 2.0  # RBMK at your pressure
  python rbmk_benchmark_sweep.py --pressure 7  --enrich 2.4  # post-Chernobyl enrichment
  add --box-source if the default source gives an error (see SOURCE below)

Long run in the background:
  nohup python rbmk_benchmark_sweep.py --pressure 7 > rbmk7.log 2>&1 &

Output: one folder per case and rbmk_<P>MPa_<enrich>pct.csv (rewritten after each case).

ASSUMPTIONS (state them in the paper)
  * 'carrier' is the central carrier rod: zirconium alloy (same as clad).
  * 'air' is helium gas (pellet hole and gaps).
  * The cell is a 25 x 25 cm graphite block fully filled by the DAGMC volumes,
    so the reflective bounding box is a true lattice boundary.
  * Fuel 900 K (shared with the other design), graphite 700 K (set the same
    value in sweep_all.py), beta assumed 650 pcm.
  * Fresh fuel; no control rods, absorbers or burnup.
"""
import argparse
import csv
import os
from math import sqrt

import openmc
from CoolProp.CoolProp import PropsSI

H5M = "/root/rbmkagain/benchmark.h5m"

ap = argparse.ArgumentParser()
ap.add_argument("--pressure", type=float, default=7.0, help="coolant pressure, MPa")
ap.add_argument("--enrich", type=float, default=2.0, help="U-235 weight percent")
ap.add_argument("--quick", action="store_true", help="short test run")
ap.add_argument("--box-source", action="store_true",
                help="use a fixed box around the cell centre instead of 'only_fissionable'")
args = ap.parse_args()

BATCHES, INACTIVE, PARTICLES = (60, 15, 20000) if args.quick else (100, 30, 75000)
VOID_FRACTIONS = [0.0, 0.2,]
BETA_PCM = 650.0  # assumed delayed neutron fraction (about 0.0065)

# ---------------- temperatures (K) and densities (g/cm3) ----------------
FUEL_T = 900
CLAD_T = 600
PRESS_T = 600
GRAPHITE_T = 700
AIR_T = 900

UO2_DENSITY = 10.5
ZR_DENSITY = 6.5
GRAPHITE_DENSITY = 1.65
HE_DENSITY = 0.00016

# ---------------- coolant properties ----------------
P = args.pressure * 1e6
RHO_WATER = PropsSI("D", "P", P, "Q", 0, "Water") / 1000.0
RHO_STEAM = PropsSI("D", "P", P, "Q", 1, "Water") / 1000.0
T_SAT = PropsSI("T", "P", P, "Q", 0, "Water")


def zr_alloy(name, nb_wo, temperature):
    m = openmc.Material(name=name)
    m.add_element("Zr", 1.0 - nb_wo, "wo")
    m.add_element("Nb", nb_wo, "wo")
    m.set_density("g/cm3", ZR_DENSITY)
    m.temperature = temperature
    return m


def build_model(alpha):
    rho_cool = alpha * RHO_STEAM + (1.0 - alpha) * RHO_WATER

    # ---------- 1. MATERIALS ----------
    air = openmc.Material(name="air")
    air.set_density("g/cm3", HE_DENSITY)
    air.add_nuclide("He4", 1.0)
    air.temperature = AIR_T

    carrier = zr_alloy("carrier", 0.01, CLAD_T)  # central carrier rod, Zr-1%Nb
    clad = zr_alloy("clad", 0.01, CLAD_T)  # fuel cladding, Zr-1%Nb
    press = zr_alloy("press", 0.025, PRESS_T)  # pressure tube, Zr-2.5%Nb

    fuel = openmc.Material(name="fuel")
    fuel.add_element("U", 1.0, enrichment=args.enrich)
    fuel.add_nuclide("O16", 2.0)
    fuel.set_density("g/cm3", UO2_DENSITY)
    fuel.temperature = FUEL_T

    graphite = openmc.Material(name="graphite")
    graphite.set_density("g/cm3", GRAPHITE_DENSITY)
    graphite.add_element("C", 1.0)
    graphite.add_s_alpha_beta("c_Graphite")
    graphite.temperature = GRAPHITE_T

    h20 = openmc.Material(name="h20")
    h20.set_density("g/cm3", rho_cool)
    h20.add_nuclide("H1", 2.0)
    h20.add_nuclide("O16", 1.0)
    h20.add_s_alpha_beta("c_H_in_H2O")
    h20.temperature = T_SAT

    materials = openmc.Materials([air, carrier, clad, fuel, graphite, h20, press])

    # ---------- 2. DAGMC GEOMETRY ----------
    dag = openmc.DAGMCUniverse(H5M, auto_geom_ids=True)

    expected = {"air", "carrier", "clad", "fuel", "graphite", "h20", "press"}
    missing = expected - set(dag.material_names)
    if missing:
        raise RuntimeError(f"Missing expected DAGMC tags: {sorted(missing)}")

    lower_left, upper_right = dag.bounding_box
    center = 0.5 * (lower_left + upper_right)

    # ---------- 3. REFLECTIVE BOUNDING REGION ----------
    min_x = openmc.XPlane(x0=lower_left[0], boundary_type="reflective")
    max_x = openmc.XPlane(x0=upper_right[0], boundary_type="reflective")
    min_y = openmc.YPlane(y0=lower_left[1], boundary_type="reflective")
    max_y = openmc.YPlane(y0=upper_right[1], boundary_type="reflective")
    min_z = openmc.ZPlane(z0=lower_left[2], boundary_type="reflective")
    max_z = openmc.ZPlane(z0=upper_right[2], boundary_type="reflective")

    bounding_region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z
    root_cell = openmc.Cell(cell_id=99999, fill=dag, region=bounding_region)
    geometry = openmc.Geometry([root_cell])

    # ---------- 4. SETTINGS (eigenvalue) ----------
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.batches = BATCHES
    settings.inactive = INACTIVE
    settings.particles = PARTICLES
    settings.temperature = {"method": "interpolation"}
    settings.source_rejection_fraction = 0.01

    # SOURCE: the cell centre is the carrier rod, not fuel, so a plain point source is poor.
    if args.box_source:
        # fixed box of +/- 4 cm (pressure-tube radius) around the centre, full height
        ll = (center[0] - 4.0, center[1] - 4.0, lower_left[2])
        ur = (center[0] + 4.0, center[1] + 4.0, upper_right[2])
        space = openmc.stats.Box(ll, ur)
    else:
        # fuel is only ~3% of the whole cell, so sample inside the pressure-tube region only
        ll = (center[0] - 4.0, center[1] - 4.0, lower_left[2])
        ur = (center[0] + 4.0, center[1] + 4.0, upper_right[2])
        space = openmc.stats.Box(ll, ur, only_fissionable=True)
    settings.source = openmc.IndependentSource(
        space=space,
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([2.0e6], [1.0]),
    )

    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.dimension = (8, 8, 8)
    entropy_mesh.lower_left = lower_left
    entropy_mesh.upper_right = upper_right
    settings.entropy_mesh = entropy_mesh

    # ---------- 5. MODEL (no tallies in the sweep) ----------
    return openmc.model.Model(geometry=geometry, materials=materials, settings=settings)


def write_csv(path, results):
    k0, s0 = results[0][1], results[0][2]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["void_fraction", "k_inf", "k_std", "delta_rho_pcm_vs_0pct",
                    "delta_rho_err_pcm", "delta_rho_in_beta"])
        for a, k, s in results:
            d = (1.0 / k0 - 1.0 / k) * 1e5  # positive = reactivity up
            e = 1e5 * sqrt((s / k**2) ** 2 + (s0 / k0**2) ** 2) if a > 0 else 0.0
            w.writerow([a, k, s, d, e, d / BETA_PCM])


def main():
    tag = f"rbmk_{args.pressure:g}MPa_{args.enrich:g}pct"
    print(f"{tag}: quick={args.quick} batches={BATCHES} inactive={INACTIVE} particles={PARTICLES}")
    print(f"liquid {RHO_WATER:.4f} g/cm3, steam {RHO_STEAM:.4f} g/cm3, T_sat {T_SAT:.1f} K\n")
    results = []
    for alpha in VOID_FRACTIONS:
        folder = f"{tag}_void{int(alpha * 100):03d}"
        os.makedirs(folder, exist_ok=True)
        print(f"--- running void fraction {alpha:.2f} ---", flush=True)
        sp_path = build_model(alpha).run(cwd=folder)
        with openmc.StatePoint(sp_path) as sp:
            k, s = sp.keff.nominal_value, sp.keff.std_dev
        results.append((alpha, k, s))
        print(f"void {alpha:.2f}: k = {k:.5f} +/- {s:.5f}", flush=True)
        write_csv(f"{tag}.csv", results)
    print(f"\nDone. See {tag}.csv")


if __name__ == "__main__":
    main()
