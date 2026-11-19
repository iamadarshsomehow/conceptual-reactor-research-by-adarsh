"""
Sweeps for the somuv4 block (OpenMC, ENDF/B-VII.1, DAGMC).

Cases
  baseline            coolant = saturated liquid at 10 MPa
  void_020 ... 100    coolant steam-water mixture at 10 MPa (void fraction)
  d2o_075/050/025pct  D2O density scaled down (partial drain)
  d2o_drained_air     D2O replaced by air (full drain)
  d2o_replaced_by_h2o D2O replaced by light water (contamination case)

Usage (inside your nuclear-env venv, after: pip install CoolProp matplotlib pandas)
  python sweep_all.py --quick       # short test run (few minutes), check it works first
  python sweep_all.py               # all cases, full statistics
  python sweep_all.py void          # baseline + void sweep only
  python sweep_all.py drain         # baseline + drain cases only

Long run in the background:
  nohup python sweep_all.py > sweep.log 2>&1 &
  tail -f sweep.log

Output: sweep_results.csv (rewritten after every case), void_results.png,
        baseline/power_map.csv (flux + kappa-fission mesh tally of the baseline)
"""
import csv
import os
import sys
from math import sqrt

import openmc
from CoolProp.CoolProp import PropsSI

H5M = "/root/somuv4/somuv4.h5m"

QUICK = "--quick" in sys.argv
GROUP = next((a for a in sys.argv[1:] if not a.startswith("--")), "all")
BATCHES, INACTIVE, PARTICLES = (60, 15, 20000) if QUICK else (100, 30, 50000)

# ---------------- coolant properties at 10 MPa ----------------
P = 10e6
RHO_WATER = PropsSI("D", "P", P, "Q", 0, "Water") / 1000.0  # g/cm3 saturated liquid
RHO_STEAM = PropsSI("D", "P", P, "Q", 1, "Water") / 1000.0  # g/cm3 saturated steam
T_SAT = PropsSI("T", "P", P, "Q", 0, "Water")  # K

VOID_FRACTIONS = [0.0, 0.2]

# ---------------- material assumptions (values from your script) ----------------
# These are kept exactly as you had them. Check each one is intended
# (see the notes in chat), then change here.
FUEL_DENSITY = 20.0
CLAD_DENSITY = 4.0
PRESSURE_TUBE_DENSITY = 4.0
GRAPHITE_DENSITY = 2.0
CARRIER_DENSITY = 2.0


def coolant_density(alpha):
    """Homogeneous steam-water mixture density in g/cm3 for void fraction alpha."""
    return alpha * RHO_STEAM + (1.0 - alpha) * RHO_WATER


def build_model(alpha=0.0, d2o_fraction=1.0, d2o_mode="d2o", tallies=False):
    # ---------- materials ----------
    air = openmc.Material(name="air")
    air.set_density("g/cm3", 0.00016)
    air.add_nuclide("He4", 1.0)
    air.temperature = 900

    carrier = openmc.Material(name="carrier")
    carrier.set_density("g/cm3", CARRIER_DENSITY)
    carrier.add_element("C", 1.0)
    carrier.add_s_alpha_beta("c_Graphite")
    carrier.temperature = 294

    clad = openmc.Material(name="clad")
    clad.set_density("g/cm3", CLAD_DENSITY)
    clad.add_element("Zr", 1)
    clad.temperature = 600

    # D2O cavity material (tag "d20" must stay the same name)
    d20 = openmc.Material(name="d20")
    if d2o_mode == "d2o":
        d20.set_density("g/cm3", 1.1 * d2o_fraction)
        d20.add_nuclide("H2", 2.0)
        d20.add_nuclide("O16", 1.0)
        d20.add_s_alpha_beta("c_D_in_D2O")
    elif d2o_mode == "air":
        d20.set_density("g/cm3", 0.0012)
        d20.add_nuclide("N14", 0.78)
        d20.add_nuclide("O16", 0.21)
        d20.add_nuclide("Ar40", 0.01)
    elif d2o_mode == "h2o":
        d20.set_density("g/cm3", 1.0)
        d20.add_nuclide("H1", 2.0)
        d20.add_nuclide("O16", 1.0)
        d20.add_s_alpha_beta("c_H_in_H2O")
    d20.temperature = 294

    fuel = openmc.Material(name="fuel")
    fuel.set_density("g/cm3", FUEL_DENSITY)
    fuel.add_nuclide("U235", 0.0073)
    fuel.add_nuclide("U238", 0.9927)
    fuel.temperature = 900

    graphite = openmc.Material(name="graphite")
    graphite.set_density("g/cm3", GRAPHITE_DENSITY)
    graphite.add_element("C", 1.0)
    graphite.add_s_alpha_beta("c_Graphite")
    graphite.temperature = 294

    # coolant: saturated water/steam mixture at 10 MPa, at saturation temperature
    h20 = openmc.Material(name="h20")
    h20.set_density("g/cm3", coolant_density(alpha))
    h20.add_nuclide("H1", 2.0)
    h20.add_nuclide("O16", 1.0)
    h20.add_s_alpha_beta("c_H_in_H2O")
    h20.temperature = T_SAT

    pressure = openmc.Material(name="pressure")
    pressure.set_density("g/cm3", PRESSURE_TUBE_DENSITY)
    pressure.add_element("Zr", 1)
    pressure.temperature = 600

    materials = openmc.Materials([air, carrier, clad, d20, fuel, graphite, h20, pressure])

    # ---------- DAGMC geometry ----------
    dag = openmc.DAGMCUniverse(H5M, auto_geom_ids=True)
    expected = {"air", "carrier", "clad", "d20", "fuel", "graphite", "h20", "pressure"}
    missing = expected - set(dag.material_names)
    if missing:
        raise RuntimeError(f"Missing expected DAGMC tags: {sorted(missing)}")

    lower_left, upper_right = dag.bounding_box
    center = 0.5 * (lower_left + upper_right)

    min_x = openmc.XPlane(x0=lower_left[0], boundary_type="reflective")
    max_x = openmc.XPlane(x0=upper_right[0], boundary_type="reflective")
    min_y = openmc.YPlane(y0=lower_left[1], boundary_type="reflective")
    max_y = openmc.YPlane(y0=upper_right[1], boundary_type="reflective")
    min_z = openmc.ZPlane(z0=lower_left[2], boundary_type="reflective")
    max_z = openmc.ZPlane(z0=upper_right[2], boundary_type="reflective")
    region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z

    root_cell = openmc.Cell(cell_id=99999, fill=dag, region=region)
    geometry = openmc.Geometry([root_cell])

    # ---------- settings ----------
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.batches = BATCHES
    settings.inactive = INACTIVE
    settings.particles = PARTICLES
    # interpolation is needed because the coolant sits at ~584 K, between library temperatures
    settings.temperature = {"method": "interpolation"}

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point(center)
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([2.0e6], [1.0])
    settings.source = source

    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.dimension = (8, 8, 8)
    entropy_mesh.lower_left = lower_left
    entropy_mesh.upper_right = upper_right
    settings.entropy_mesh = entropy_mesh

    model = openmc.model.Model(geometry=geometry, materials=materials, settings=settings)

    # ---------- tallies (baseline only) ----------
    if tallies:
        mesh = openmc.RegularMesh()
        mesh.dimension = (20, 20, 20)
        mesh.lower_left = lower_left
        mesh.upper_right = upper_right
        tally = openmc.Tally(name="power_map")
        tally.filters = [openmc.MeshFilter(mesh)]
        tally.scores = ["flux", "kappa-fission"]
        model.tallies = openmc.Tallies([tally])

    return model


def case_list(group):
    cases = [("baseline", dict(alpha=0.0, tallies=True))]
    if group in ("void", "all"):
        cases += [(f"void_{int(a * 100):03d}", dict(alpha=a)) for a in VOID_FRACTIONS[1:]]
    if group in ("drain", "all"):
        cases += [(f"d2o_{int(f * 100):03d}pct", dict(d2o_fraction=f)) for f in (0.75, 0.5, 0.25)]
        cases += [
            ("d2o_drained_air", dict(d2o_mode="air")),
            ("d2o_replaced_by_h2o", dict(d2o_mode="h2o")),
        ]
    return cases


def run_case(name, kwargs):
    os.makedirs(name, exist_ok=True)
    model = build_model(**kwargs)
    sp_path = model.run(cwd=name)
    with openmc.StatePoint(sp_path) as sp:
        k = sp.keff
        if kwargs.get("tallies"):
            try:
                df = sp.get_tally(name="power_map").get_pandas_dataframe()
                df.to_csv(os.path.join(name, "power_map.csv"), index=False)
            except Exception as exc:  # keep the sweep going if post-processing fails
                print("could not export power map:", exc, flush=True)
    return k.nominal_value, k.std_dev


def write_results(results):
    k0, s0 = results[0][1], results[0][2]
    rows = []
    for name, k, s in results:
        d_rho = (1.0 / k0 - 1.0 / k) * 1e5  # pcm vs baseline; positive = reactivity up
        err = 1e5 * sqrt((s / k**2) ** 2 + (s0 / k0**2) ** 2) if name != "baseline" else 0.0
        rows.append((name, k, s, d_rho, err))
    with open("sweep_results.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "k_inf", "k_std", "delta_rho_pcm_vs_baseline", "delta_rho_err_pcm"])
        w.writerows(rows)
    return rows


def plot_void(rows):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    pts = [(0.0, r[3], r[4]) for r in rows if r[0] == "baseline"]
    pts += [(int(r[0][5:]) / 100.0, r[3], r[4]) for r in rows if r[0].startswith("void_")]
    if len(pts) < 2:
        return
    plt.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[p[2] for p in pts],
                 marker="o", capsize=3)
    plt.axhline(0, color="gray", lw=0.8)
    plt.xlabel("Void fraction")
    plt.ylabel("Reactivity change vs 0% void (pcm)")
    plt.title("Void effect at 10 MPa")
    plt.savefig("void_results.png", dpi=150, bbox_inches="tight")


def main():
    print(f"group={GROUP} quick={QUICK} batches={BATCHES} inactive={INACTIVE} particles={PARTICLES}")
    print(f"10 MPa: liquid {RHO_WATER:.4f} g/cm3, steam {RHO_STEAM:.4f} g/cm3, T_sat {T_SAT:.1f} K\n")
    results = []
    rows = []
    for name, kwargs in case_list(GROUP):
        print(f"--- running {name} ---", flush=True)
        k, s = run_case(name, kwargs)
        results.append((name, k, s))
        print(f"{name}: k = {k:.5f} +/- {s:.5f}", flush=True)
        rows = write_results(results)
    plot_void(rows)
    print("\nDone. See sweep_results.csv")


if __name__ == "__main__":
    main()
