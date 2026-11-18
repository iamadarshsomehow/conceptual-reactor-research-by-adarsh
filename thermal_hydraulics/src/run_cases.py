"""run_cases.py - run the planned cases through channel.py and print one table.

Run from the folder that holds channel.py, config.yaml, props.py, power_map.csv:
    python run_cases.py
Each case runs in cases/<name>/ (its own CSV, plot and summary.json).
Writes cases_summary.csv with the table.
"""
import copy
import json
import os
import subprocess
import sys

import pandas as pd
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "config.yaml")) as fh:
    BASE = yaml.safe_load(fh)

# Make the power-map path absolute so every case folder can find it.
BASE["power_shape"]["power_map_csv"] = os.path.join(HERE, BASE["power_shape"]["power_map_csv"])


def variant(**changes):
    cfg = copy.deepcopy(BASE)
    for key, val in changes.items():
        section, name = key.split("__")
        cfg[section][name] = val
    return cfg


CASES = {
    "base":         variant(),
    "nodes_100":    variant(numerics__n_axial_nodes=100),
    "hot_1.2":      variant(power_shape__radial_peaking=1.2),
    "gap_5000":     variant(geometry__gap_conductance_W_m2K=5000),
    "gap_10000":    variant(geometry__gap_conductance_W_m2K=10000),
    "cosine_axial": variant(power_shape__axial_shape="cosine"),
    "cosine_hot1.2": variant(power_shape__axial_shape="cosine", power_shape__radial_peaking=1.2),
}

rows = []
for name, cfg in CASES.items():
    folder = os.path.join(HERE, "cases", name)
    os.makedirs(folder, exist_ok=True)
    cfg_path = os.path.join(folder, "config.yaml")
    with open(cfg_path, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
    res = subprocess.run([sys.executable, os.path.join(HERE, "channel.py"), cfg_path],
                         cwd=folder, capture_output=True, text=True,
                         env={**os.environ, "PYTHONPATH": HERE + os.pathsep + os.environ.get("PYTHONPATH", "")})
    with open(os.path.join(folder, "log.txt"), "w") as fh:
        fh.write(res.stdout + res.stderr)
    try:
        with open(os.path.join(folder, "summary.json")) as fh:
            s = json.load(fh)
    except Exception:
        print(f"{name}: FAILED (see cases/{name}/log.txt)")
        print(res.stderr[-800:])
        continue
    rows.append({
        "case": name,
        "exit_x": round(s["exit_quality"], 4),
        "exit_void": round(s["exit_void"], 3),
        "dp_kPa": round(s["dp_total_kPa"], 2),
        "T_wall_max_C": round(s["T_wall_max_C"], 1),
        "T_fuel_centre_max_C": round(s["T_fuel_centre_max_C"], 1),
        "margin_650_K": round(s["fuel_margin_K"], 1),
        "margin_600_K": round(s["fuel_margin_strict_K"], 1),
        "T_clad_in_max_C": round(s["T_clad_inner_max_C"], 1),
        "interface_margin_K": round(s["interface_margin_K"], 1),
        "CHF_margin_raw": round(s["chf_margin_raw_min"], 2),
        "CHF_margin_derated": round(s["chf_margin_derated_min"], 2),
        "checks_passed": s["all_checks_passed"],
    })

df = pd.DataFrame(rows)
df.to_csv(os.path.join(HERE, "cases_summary.csv"), index=False)
print("\n=== CASE COMPARISON ===")
print(df.to_string(index=False))
print(f"\nLimits: fuel {BASE['limits']['fuel_temp_limit_C']} C (strict {BASE['limits']['fuel_temp_limit_strict_C']} C), interface {BASE['limits']['interface_temp_limit_C']} C, "
      f"CHF margin {BASE['limits']['min_chf_margin']} (derating {BASE['limits'].get('chf_derating_factor', 1.0)}).")
print("Negative fuel_margin_K = limit exceeded.  Saved cases_summary.csv")

# ---- comparison graph: fuel centre and clad-inner temperature per case
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 5.5))
xs = range(len(df))
ax.bar([i - 0.2 for i in xs], df["T_fuel_centre_max_C"], width=0.4, label="max fuel centreline T")
ax.bar([i + 0.2 for i in xs], df["T_clad_in_max_C"], width=0.4, label="max clad inner T")
ax.axhline(BASE["limits"]["fuel_temp_limit_C"], color="k", ls="--", label="fuel limit 650 C")
ax.axhline(BASE["limits"]["fuel_temp_limit_strict_C"], color="r", ls=":", label="fuel limit 600 C")
ax.axhline(BASE["limits"]["interface_temp_limit_C"], color="m", ls="-.", label="interface limit 550 C")
ax.set_xticks(list(xs)); ax.set_xticklabels(df["case"], rotation=30, ha="right")
ax.set_ylabel("temperature (C)"); ax.set_title("Temperature limits by case (channel power from config.yaml)")
ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "cases_comparison.png"), dpi=130)
print("Saved cases_comparison.png")
