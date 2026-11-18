"""power_scan.py - how much channel power can the concept carry?

Runs channel.py for a list of powers in several cases, then finds by linear
interpolation the power at which each limit is reached:
  fuel centre 650 C, fuel centre 600 C, interface (clad inner) 550 C,
  derated CHF margin = limit.
Run from the folder with channel.py, config.yaml, props.py, power_map.csv:
    python power_scan.py
Writes scan_all.csv, allowable_power.csv, power_scan.png.  About 3-6 minutes.
"""
import copy, json, os, subprocess, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "config.yaml")) as fh:
    BASE = yaml.safe_load(fh)
BASE["power_shape"]["power_map_csv"] = os.path.join(HERE, BASE["power_shape"]["power_map_csv"])
L = BASE["limits"]
POWERS = [round(float(p), 3) for p in np.arange(0.10, 0.401, 0.025)]    # MW

def variant(**ch):
    c = copy.deepcopy(BASE)
    for k, v in ch.items():
        a, b = k.split("__"); c[a][b] = v
    return c

CASES = {
    "base (openmc axial, gap 15000)": variant(),
    "gap 10000": variant(geometry__gap_conductance_W_m2K=10000),
    "gap 5000": variant(geometry__gap_conductance_W_m2K=5000),
    "hot channel x1.2": variant(power_shape__radial_peaking=1.2),
    "cosine axial": variant(power_shape__axial_shape="cosine"),
    "cosine + hot x1.2": variant(power_shape__axial_shape="cosine", power_shape__radial_peaking=1.2),
}

def run(cfg, tag):
    folder = os.path.join(HERE, "scan_runs", tag)
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, "config.yaml")
    with open(path, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    env = {**os.environ, "CHANNEL_NOPLOT": "1",
           "PYTHONPATH": HERE + os.pathsep + os.environ.get("PYTHONPATH", "")}
    r = subprocess.run([sys.executable, os.path.join(HERE, "channel.py"), path],
                       cwd=folder, capture_output=True, text=True, env=env)
    try:
        with open(os.path.join(folder, "summary.json")) as fh:
            return json.load(fh)
    except Exception:
        print(f"  FAILED {tag}: {r.stderr[-300:]}"); return None

rows = []
for cname, cfg in CASES.items():
    print(f"case: {cname}")
    for p in POWERS:
        c = copy.deepcopy(cfg); c["operating_point"]["channel_power_MW"] = p
        s = run(c, f"{cname.split(' ')[0]}_{abs(hash(cname))%10000}_{p}")
        if s is None:
            continue
        rows.append(dict(case=cname, power_MW=p, T_fuel_centre_C=s["T_fuel_centre_max_C"],
                         T_clad_inner_C=s["T_clad_inner_max_C"],
                         chf_derated=s["chf_margin_derated_min"], chf_raw=s["chf_margin_raw_min"],
                         exit_quality=s["exit_quality"], exit_void=s["exit_void"],
                         mdot=s["mdot_kg_s"]))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(HERE, "scan_all.csv"), index=False)

def cross(x, y, target, rising=True):
    """x where y first crosses target (linear interpolation); None if never."""
    x = np.asarray(x); y = np.asarray(y)
    for i in range(len(x) - 1):
        a, b = y[i], y[i + 1]
        if (rising and a <= target < b) or ((not rising) and a >= target > b):
            return x[i] + (target - a) * (x[i + 1] - x[i]) / (b - a)
    if (rising and y[0] > target) or ((not rising) and y[0] < target):
        return float(x[0])          # already beyond the limit at the lowest power
    return None                      # limit never reached within the scan

out = []
for cname, d in df.groupby("case", sort=False):
    d = d.sort_values("power_MW")
    out.append(dict(case=cname,
        P_fuel_650C=cross(d.power_MW, d.T_fuel_centre_C, L["fuel_temp_limit_C"]),
        P_fuel_600C=cross(d.power_MW, d.T_fuel_centre_C, L["fuel_temp_limit_strict_C"]),
        P_interface=cross(d.power_MW, d.T_clad_inner_C, L["interface_temp_limit_C"]),
        P_chf=cross(d.power_MW, d.chf_derated, L["min_chf_margin"], rising=False)))
res = pd.DataFrame(out).round(3)
res["allowable_MW_(650C)"] = res[["P_fuel_650C", "P_interface", "P_chf"]].min(axis=1, skipna=True)
res["allowable_MW_(600C)"] = res[["P_fuel_600C", "P_interface", "P_chf"]].min(axis=1, skipna=True)
res.to_csv(os.path.join(HERE, "allowable_power.csv"), index=False)
print("\n=== POWER AT WHICH EACH LIMIT IS REACHED (MW, channel; None = not reached up to "
      f"{POWERS[-1]} MW; a value equal to {POWERS[0]} means already exceeded at the lowest power) ===")
print(res.to_string(index=False))

fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))
for cname, d in df.groupby("case", sort=False):
    d = d.sort_values("power_MW")
    ax[0].plot(d.power_MW, d.T_fuel_centre_C, marker="o", ms=3, label=cname)
    ax[1].plot(d.power_MW, d.chf_derated, marker="o", ms=3, label=cname)
ax[0].axhline(L["fuel_temp_limit_C"], color="k", ls="--", label="650 C")
ax[0].axhline(L["fuel_temp_limit_strict_C"], color="r", ls=":", label="600 C")
ax[0].axhline(L["interface_temp_limit_C"], color="m", ls="-.", label="interface 550 C (clad inner)")
ax[0].set_xlabel("channel power (MW)"); ax[0].set_ylabel("max fuel centreline T (C)")
ax[0].set_title("Fuel temperature versus channel power"); ax[0].legend(fontsize=7)
ax[1].axhline(L["min_chf_margin"], color="k", ls="--", label="CHF margin limit")
ax[1].set_xlabel("channel power (MW)"); ax[1].set_ylabel("minimum derated CHF margin")
ax[1].set_title("CHF margin versus channel power"); ax[1].legend(fontsize=7)
for a in ax: a.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "power_scan.png"), dpi=130)
print("\nSaved scan_all.csv, allowable_power.csv, power_scan.png")
