"""
Baseline / D2O-drain / LOCA runs for the hexagonal-block design (OpenMC, ENDF/B-VII.1, DAGMC).

MODES
  python drain_runs.py                 baseline only (1 run)
  python drain_runs.py drain           D2O drained 0, 25, 50, 75, 100 %      (5 runs)
  python drain_runs.py loca            D2O AND light water drained together,
                                       0, 25, 50, 75, 100 %                   (5 runs)
  python drain_runs.py plot            rebuild tables and graphs from saved CSVs
                                       (no OpenMC needed; overlays drain vs loca if both exist)

OPTIONS
  --quick              short test run (results get a _quick suffix, nothing is overwritten)
  --rect               old rectangular reflective box instead of the hexagonal prism
  --hex-orientation    auto (default), x or y  (only if the hexagon does not line up with your CAD)

Long run in the background:
  nohup python drain_runs.py drain > drain.log 2>&1 &
  tail -f drain.log

OUTPUT (per mode)
  results_<mode>.csv          raw numbers
  results_<mode>.md           presentable table:  k-inf = 1.04745 +/- 0.00069
  results_<mode>_table.png    the same table as an image for slides or the paper
  results_<mode>_plot.png     k-inf and reactivity change versus % drained
  results_compare_plot.png    drain vs loca on one graph (plot mode, when both exist)
  power_map.csv               (baseline run) flux + kappa-fission on a 20x20x20 mesh
  power_map_plot.png          (baseline run / plot mode) radial power and flux maps + axial profile
  entropy_<mode>.csv          Shannon entropy per batch for every case
  entropy_<mode>_plot.png     entropy vs batch, one line per case, dashed line where tallies start;
                              each line is labelled "settled" or "STILL DRIFTING" (rough guide)

HOW A DRAIN IS MODELLED
  Uniform density scaling: drained fraction f -> density x (1 - f), thermal-scattering table kept.
  At 100 % the material is replaced by air. 0 % is the undisturbed baseline.

BEFORE THE FIRST RUN
  Set the values in the MATERIAL SETTINGS block to match the model you have been using.
  The 0 % case should then reproduce your earlier baseline k within its error bar.
"""
import argparse
import csv
import os
from pathlib import Path
from math import sqrt

ap = argparse.ArgumentParser(description="baseline / D2O drain / LOCA runs")
ap.add_argument("mode", nargs="?", default="baseline", choices=["baseline", "drain", "loca", "plot"])
ap.add_argument("--quick", action="store_true", help="short test run")
ap.add_argument("--rect", action="store_true", help="rectangular reflective box instead of hex prism")
ap.add_argument("--hex-orientation", default="auto", choices=["auto", "x", "y"])
args = ap.parse_args()

SUF = "_quick" if args.quick else ""
BATCHES, INACTIVE, PARTICLES = (60, 15, 20000) if args.quick else (200, 30, 100000)
LEVELS = [0, 25, 50, 75, 100]  # percent drained
LEAKAGE = 0.03                 # assumed neutron leakage of a large RBMK-sized core (ASSUMPTION)
BASELINE_NOTE = ("k-eff = k-inf \u00d7 (1 \u2212 0.03) = k-inf \u00d7 0.97. This assumes 3 % neutron leakage "
                 "for a large RBMK-sized core. It is an assumption, not a calculated value; "
                 "the real leakage depends on core size and reflector.")

# ======================= MATERIAL SETTINGS (match your current model) =======================
H5M = os.environ.get(
    "SOMUV4_H5M", str(Path(__file__).resolve().parents[2] / "private" / "somuv4.h5m")
)
PRESSURE_MPA = 10.0

FUEL_DENSITY = 19.05           # g/cm3, natural uranium metal (textbook value); set YOUR value
CLAD_DENSITY = 4.0             # g/cm3  (as in sweep_all.py)
PRESSURE_TUBE_DENSITY = 4.0    # g/cm3  (as in sweep_all.py)
GRAPHITE_DENSITY = 2.0         # g/cm3  (as in sweep_all.py)
CARRIER_DENSITY = 2.0          # g/cm3  (as in sweep_all.py)
D2O_DENSITY = 1.1              # g/cm3

FUEL_T, CLAD_T, PRESSURE_T = 900, 600, 600
GRAPHITE_T = 700               # K, use the same value in every model you compare
D2O_T = 294
AIR_T = 900
# ============================================================================================

CSV_HEADER = ["case", "drained_percent", "k_inf", "k_std", "delta_rho_pcm", "delta_rho_err_pcm"]
LABELS = {"drain": "D2O drained", "loca": "D2O + water drained (LOCA)", "baseline": "baseline"}


# ------------------------------------------------------------------ model
def build_model(d2o_drain=0.0, h2o_drain=0.0, tallies=False):
    import openmc

    def make_air(name, temperature):
        m = openmc.Material(name=name)
        m.add_nuclide("N14", 0.78)
        m.add_nuclide("O16", 0.21)
        m.add_nuclide("Ar40", 0.01)
        m.set_density("g/cm3", 0.0012)
        m.temperature = temperature
        return m

    air = openmc.Material(name="air")  # helium gas tag
    air.set_density("g/cm3", 0.00016)
    air.add_nuclide("He4", 1.0)
    air.temperature = AIR_T

    carrier = openmc.Material(name="carrier")
    carrier.set_density("g/cm3", CARRIER_DENSITY)
    carrier.add_element("C", 1.0)
    carrier.add_s_alpha_beta("c_Graphite")
    carrier.temperature = GRAPHITE_T

    clad = openmc.Material(name="clad")
    clad.set_density("g/cm3", CLAD_DENSITY)
    clad.add_element("Zr", 1)
    clad.temperature = CLAD_T

    if d2o_drain >= 1.0:
        d20 = make_air("d20", D2O_T)
    else:
        d20 = openmc.Material(name="d20")
        d20.set_density("g/cm3", D2O_DENSITY * (1.0 - d2o_drain))
        d20.add_nuclide("H2", 2.0)
        d20.add_nuclide("O16", 1.0)
        d20.add_s_alpha_beta("c_D_in_D2O")
        d20.temperature = D2O_T

    fuel = openmc.Material(name="fuel")  # natural uranium metal
    fuel.set_density("g/cm3", FUEL_DENSITY)
    fuel.add_nuclide("U235", 0.0073)
    fuel.add_nuclide("U238", 0.9927)
    fuel.temperature = FUEL_T

    graphite = openmc.Material(name="graphite")
    graphite.set_density("g/cm3", GRAPHITE_DENSITY)
    graphite.add_element("C", 1.0)
    graphite.add_s_alpha_beta("c_Graphite")
    graphite.temperature = GRAPHITE_T

    if h2o_drain >= 1.0:
        h20 = make_air("h20", T_SAT)
    else:
        h20 = openmc.Material(name="h20")  # saturated liquid water at the chosen pressure
        h20.set_density("g/cm3", RHO_WATER * (1.0 - h2o_drain))
        h20.add_nuclide("H1", 2.0)
        h20.add_nuclide("O16", 1.0)
        h20.add_s_alpha_beta("c_H_in_H2O")
        h20.temperature = T_SAT

    pressure = openmc.Material(name="pressure")
    pressure.set_density("g/cm3", PRESSURE_TUBE_DENSITY)
    pressure.add_element("Zr", 1)
    pressure.temperature = PRESSURE_T

    materials = openmc.Materials([air, carrier, clad, d20, fuel, graphite, h20, pressure])

    dag = openmc.DAGMCUniverse(H5M, auto_geom_ids=True)
    expected = {"air", "carrier", "clad", "d20", "fuel", "graphite", "h20", "pressure"}
    missing = expected - set(dag.material_names)
    if missing:
        raise RuntimeError(f"Missing expected DAGMC tags: {sorted(missing)}")

    lower_left, upper_right = dag.bounding_box
    center = 0.5 * (lower_left + upper_right)

    min_z = openmc.ZPlane(z0=lower_left[2], boundary_type="reflective")
    max_z = openmc.ZPlane(z0=upper_right[2], boundary_type="reflective")

    if args.rect:
        min_x = openmc.XPlane(x0=lower_left[0], boundary_type="reflective")
        max_x = openmc.XPlane(x0=upper_right[0], boundary_type="reflective")
        min_y = openmc.YPlane(y0=lower_left[1], boundary_type="reflective")
        max_y = openmc.YPlane(y0=upper_right[1], boundary_type="reflective")
        region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z
    else:
        dx = upper_right[0] - lower_left[0]
        dy = upper_right[1] - lower_left[1]
        # flat-to-flat is the SMALLER extent of the bounding box, edge = flat-to-flat / sqrt(3)
        edge = min(dx, dy) / sqrt(3.0)
        orient = args.hex_orientation
        if orient == "auto":
            orient = "y" if dx < dy else "x"
        prism = openmc.model.HexagonalPrism(edge_length=edge, orientation=orient,
                                            origin=(center[0], center[1]),
                                            boundary_type="reflective")
        region = -prism & +min_z & -max_z

    root_cell = openmc.Cell(cell_id=99999, fill=dag, region=region)
    geometry = openmc.Geometry([root_cell])

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.batches = BATCHES
    settings.inactive = INACTIVE
    settings.particles = PARTICLES
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

    if tallies:  # flux + fission power mesh map (baseline run only)
        mesh = openmc.RegularMesh()
        mesh.dimension = (50,50, 50)
        mesh.lower_left = lower_left
        mesh.upper_right = upper_right
        tally = openmc.Tally(name="power_map")
        tally.filters = [openmc.MeshFilter(mesh)]
        tally.scores = ["flux", "kappa-fission"]
        model.tallies = openmc.Tallies([tally])
    return model


# ------------------------------------------------------------------ results handling
def csv_path(mode):
    return f"results_{mode}{SUF}.csv"


def save_csv(mode, rows):
    with open(csv_path(mode), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(CSV_HEADER)
        for r in rows:
            w.writerow([r["case"], r["pct"], r["k"], r["k_std"], r["d_rho"], r["d_rho_err"]])


def load_csv(mode):
    path = csv_path(mode)
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        rows = [dict(case=r["case"], pct=float(r["drained_percent"]), k=float(r["k_inf"]),
                     k_std=float(r["k_std"]), d_rho=0.0, d_rho_err=0.0)
                for r in csv.DictReader(f)]
    return add_reactivity(rows) if rows else None


def add_reactivity(rows):
    """Reactivity change versus the first row (0 % drained), in pcm, with propagated error."""
    k0, s0 = rows[0]["k"], rows[0]["k_std"]
    for r in rows:
        k, s = r["k"], r["k_std"]
        r["d_rho"] = (1.0 / k0 - 1.0 / k) * 1e5
        r["d_rho_err"] = 0.0 if r is rows[0] else 1e5 * sqrt((s / k**2) ** 2 + (s0 / k0**2) ** 2)
    return rows


def fmt_k(r):
    return f"{r['k']:.5f} \u00b1 {r['k_std']:.5f}"


def table_cells(mode, rows):
    if mode == "baseline":
        head = ["Case", "k-inf (\u00b1 1\u03c3)", "k-eff, 3 % leakage*"]
        body = [[r["case"], fmt_k(r),
                 f"{r['k'] * (1 - LEAKAGE):.5f} \u00b1 {r['k_std'] * (1 - LEAKAGE):.5f}"] for r in rows]
        return head, body
    head = ["Drained (%)", "k-inf (\u00b1 1\u03c3)", "\u0394\u03c1 vs 0% (pcm)"]
    body = []
    for r in rows:
        d = "-" if r is rows[0] else f"{r['d_rho']:+.0f} \u00b1 {r['d_rho_err']:.0f}"
        body.append([f"{r['pct']:.0f}", fmt_k(r), d])
    return head, body


def write_table(mode, rows):
    head, body = table_cells(mode, rows)
    print("\n" + "  |  ".join(head))
    for b in body:
        print("  |  ".join(b))
    if mode == "baseline":
        print("* " + BASELINE_NOTE)
    with open(f"results_{mode}{SUF}.md", "w", encoding="utf-8") as f:
        f.write(f"**{LABELS[mode]}**: k-infinity versus drained fraction\n\n")
        f.write("| " + " | ".join(head) + " |\n")
        f.write("|" + "|".join(["---"] * len(head)) + "|\n")
        for b in body:
            f.write("| " + " | ".join(b) + " |\n")
        if mode == "baseline":
            f.write("\n\\* " + BASELINE_NOTE + "\n")
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7.5, 0.7 + 0.42 * (len(body) + 1)))
        ax.axis("off")
        tbl = ax.table(cellText=body, colLabels=head, loc="center", cellLoc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(11)
        tbl.scale(1, 1.5)
        for j in range(len(head)):
            tbl[0, j].set_text_props(weight="bold")
        ax.set_title(LABELS[mode], fontsize=12, pad=8)
        if mode == "baseline":
            import textwrap

            ax.text(0.5, 0.0, "* " + "\n".join(textwrap.wrap(BASELINE_NOTE, 80)), transform=ax.transAxes,
                    ha="center", va="top", fontsize=8.5, color="#333")
        fig.savefig(f"results_{mode}{SUF}_table.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
    except ImportError:
        print("matplotlib not installed: skipped the table image")


def draw(series, path, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    for label, rows in series:
        x = [r["pct"] for r in rows]
        ax1.errorbar(x, [r["k"] for r in rows], yerr=[r["k_std"] for r in rows],
                     marker="o", capsize=3, label=label)
        ax2.errorbar(x, [r["d_rho"] for r in rows], yerr=[r["d_rho_err"] for r in rows],
                     marker="o", capsize=3, label=label)
    ax1.axhline(1.0, color="gray", ls="--", lw=1, label="k = 1")
    ax2.axhline(0.0, color="gray", ls="--", lw=1)
    ax1.set_xlabel("Drained (%)")
    ax1.set_ylabel("k-infinity")
    ax1.set_title("k-infinity vs drained fraction")
    ax2.set_xlabel("Drained (%)")
    ax2.set_ylabel("Reactivity change vs 0% (pcm)")
    ax2.set_title("Reactivity change")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3)
    ax1.legend()
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_graph(mode, rows):
    try:
        draw([(LABELS[mode], rows)], f"results_{mode}{SUF}_plot.png", LABELS[mode])
    except ImportError:
        print("matplotlib not installed: skipped the graph")


# ------------------------------------------------------------------ Shannon entropy
def entropy_csv_path(mode):
    return f"entropy_{mode}{SUF}.csv"


def save_entropy(mode, entropies):
    names = list(entropies)
    n = max(len(v) for v in entropies.values())
    with open(entropy_csv_path(mode), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["batch"] + names)
        for i in range(n):
            w.writerow([i + 1] + [f"{entropies[c][i]:.6f}" if i < len(entropies[c]) else "" for c in names])


def load_entropy(mode):
    path = entropy_csv_path(mode)
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        names = [c for c in rd.fieldnames if c != "batch"]
        data = {c: [] for c in names}
        for r in rd:
            for c in names:
                if r[c] != "":
                    data[c].append(float(r[c]))
    return data or None


def entropy_verdicts(entropies, inactive):
    """Rough numeric guide: is the end of the inactive phase at the same level as the active phase?"""
    import statistics as st

    out = {}
    for name, e in entropies.items():
        act = e[inactive:]
        tail = e[max(0, inactive - 10):inactive]
        if len(act) < 5 or len(tail) < 3:
            out[name] = ("n/a", 0.0)
            continue
        sd = st.pstdev(act) or 1e-12
        z = abs(st.mean(tail) - st.mean(act)) / sd
        out[name] = ("settled" if z < 2.0 else "STILL DRIFTING", z)
    return out


def draw_entropy(mode, entropies):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed: skipped the entropy graph")
        return
    verdicts = entropy_verdicts(entropies, INACTIVE)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for name, e in entropies.items():
        ax.plot(range(1, len(e) + 1), e, lw=1.2, label=f"{name}: {verdicts[name][0]}")
    ax.axvline(INACTIVE + 0.5, color="gray", ls="--", lw=1)
    ax.text(INACTIVE + 1, ax.get_ylim()[0], " tallies start", va="bottom", fontsize=8, color="gray")
    ax.set_xlabel("Batch")
    ax.set_ylabel("Shannon entropy (bits)")
    ax.set_title(f"Source convergence: {LABELS[mode]}")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.savefig(f"entropy_{mode}{SUF}_plot.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Entropy check (rough guide; judge by the graph):")
    for name, (v, z) in verdicts.items():
        print(f"  {name}: {v} (offset {z:.1f} sigma)")


# ------------------------------------------------------------------ flux / power map
def power_csv_path():
    return f"power_map{SUF}.csv"


def save_power(sp):
    """Write the baseline mesh tally (flux, kappa-fission) to power_map.csv: ix,iy,iz,score,mean,std_dev."""
    df = sp.get_tally(name="power_map").get_pandas_dataframe()
    df.columns = ["_".join(str(x) for x in c if str(x)) if isinstance(c, tuple) else str(c)
                  for c in df.columns]
    xc = next(c for c in df.columns if c.endswith("_x"))
    yc = next(c for c in df.columns if c.endswith("_y"))
    zc = next(c for c in df.columns if c.endswith("_z"))
    sdc = next(c for c in df.columns if c.startswith("std"))
    with open(power_csv_path(), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ix", "iy", "iz", "score", "mean", "std_dev"])
        for _, r in df.iterrows():
            w.writerow([int(r[xc]), int(r[yc]), int(r[zc]), r["score"], r["mean"], r[sdc]])


def draw_power():
    path = power_csv_path()
    if not os.path.exists(path):
        return
    try:
        import numpy as np
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("numpy/matplotlib not installed: skipped the power map")
        return
    data = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            data.setdefault(r["score"], []).append((int(r["ix"]), int(r["iy"]), int(r["iz"]), float(r["mean"])))
    arrs = {}
    for score, pts in data.items():
        n = max(max(p[0], p[1], p[2]) for p in pts)
        a = np.zeros((n, n, n))
        for ix, iy, iz, m in pts:
            a[ix - 1, iy - 1, iz - 1] = m
        arrs[score] = a
    pw, fl = arrs["kappa-fission"], arrs["flux"]

    def radial(a):
        img = a.sum(axis=2).T  # x horizontal, y vertical, summed over height
        img = img / img.max()
        return np.ma.masked_where(img <= 0, img)

    def axial(a):
        prof = a.sum(axis=(0, 1))
        return prof / prof.mean()

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    for ax, a, title in ((axs[0], pw, "Fission power (kappa-fission)\nsummed over height, relative to peak"),
                         (axs[1], fl, "Neutron flux\nsummed over height, relative to peak")):
        im = ax.imshow(radial(a), origin="lower", cmap="viridis", extent=(0, 100, 0, 100))
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("x (% of bounding box)")
        ax.set_ylabel("y (% of bounding box)")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    z = (np.arange(pw.shape[2]) + 0.5) / pw.shape[2] * 100
    axs[2].plot(z, axial(pw), marker="o", ms=3, label="fission power")
    axs[2].plot(z, axial(fl), marker="s", ms=3, label="flux")
    axs[2].axhline(1.0, color="gray", ls="--", lw=1)
    axs[2].set_xlabel("Height (% of cell)")
    axs[2].set_ylabel("Relative to axial average")
    axs[2].set_title("Axial profile", fontsize=10)
    axs[2].grid(alpha=0.3)
    axs[2].legend()
    fig.suptitle("Baseline flux and power map (20 x 20 x 20 mesh, per source particle, shape only)")
    fig.tight_layout()
    fig.savefig(f"power_map{SUF}_plot.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Axial peak / average: power {axial(pw).max():.2f}, flux {axial(fl).max():.2f} "
          f"(coarse mesh: peaking values depend on mesh size)")


# ------------------------------------------------------------------ runs
def run_mode(mode):
    import openmc

    if mode == "baseline":
        cases = [("baseline", 0, 0.0, 0.0)]
    elif mode == "drain":
        cases = [(f"d2o_drain_{p:03d}", p, p / 100.0, 0.0) for p in LEVELS]
    else:  # loca
        cases = [(f"loca_{p:03d}", p, p / 100.0, p / 100.0) for p in LEVELS]

    rows = []
    entropies = {}
    for name, pct, d2o, h2o in cases:
        folder = name + SUF
        os.makedirs(folder, exist_ok=True)
        print(f"--- running {name} (D2O drained {d2o * 100:.0f}%, water drained {h2o * 100:.0f}%) ---",
              flush=True)
        sp_path = build_model(d2o, h2o, tallies=(mode == "baseline")).run(cwd=folder)
        with openmc.StatePoint(sp_path) as sp:
            k, s = sp.keff.nominal_value, sp.keff.std_dev
            entropies[name] = [float(e) for e in sp.entropy]
            if mode == "baseline":
                save_power(sp)
        rows.append(dict(case=name, pct=pct, k=k, k_std=s, d_rho=0.0, d_rho_err=0.0))
        print(f"{name}: k = {k:.5f} +/- {s:.5f}", flush=True)
        save_csv(mode, add_reactivity(rows))  # saved after every case
        save_entropy(mode, entropies)
    draw_entropy(mode, entropies)
    if mode == "baseline":
        draw_power()
    return rows


def report_all():
    found = {}
    for mode in ("baseline", "drain", "loca"):
        rows = load_csv(mode)
        if rows:
            found[mode] = rows
            write_table(mode, rows)
            if mode != "baseline":
                make_graph(mode, rows)
            ent = load_entropy(mode)
            if ent:
                draw_entropy(mode, ent)
    draw_power()
    if not found:
        raise SystemExit(f"No result CSVs found (looked for results_<mode>{SUF}.csv in this folder)")
    series = [(LABELS[m], found[m]) for m in ("drain", "loca") if m in found]
    if len(series) == 2:
        draw(series, f"results_compare{SUF}_plot.png", "D2O drain vs LOCA (D2O + water)")
        print(f"\nWrote results_compare{SUF}_plot.png")


if args.mode == "plot":
    report_all()
else:
    import openmc  # noqa: F401
    from CoolProp.CoolProp import PropsSI

    P = PRESSURE_MPA * 1e6
    RHO_WATER = PropsSI("D", "P", P, "Q", 0, "Water") / 1000.0  # g/cm3, saturated liquid
    T_SAT = PropsSI("T", "P", P, "Q", 0, "Water")  # K
    print(f"mode={args.mode} quick={args.quick} batches={BATCHES} inactive={INACTIVE} "
          f"particles={PARTICLES}")
    print(f"{PRESSURE_MPA:g} MPa saturated liquid {RHO_WATER:.4f} g/cm3 at {T_SAT:.1f} K\n")
    result_rows = run_mode(args.mode)
    write_table(args.mode, result_rows)
    if args.mode != "baseline":
        make_graph(args.mode, result_rows)
    print("\nDone.")
