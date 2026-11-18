"""channel.py - 1D channel model, LAYERS 1 + 2 + 3 + 4.

Layer 1: axial energy balance (enthalpy, water temperature, steam quality).
Layer 2: void fraction (drift-flux) and pressure drop.
Layer 3: wall, clad and fuel temperatures.
Layer 4: critical heat flux (CHF) margin, Biasi correlation.

Run:  python channel.py            (uses config.yaml in this folder)
      python channel.py other.yaml
Outputs: channel_results.csv, channel_energy_balance.png

Modelling choices and their limits (state these in the paper):
  * pressure_MPa in config.yaml is the channel INLET pressure.
  * Void fraction comes from equilibrium quality only; void from subcooled
    boiling is neglected (understates void low in the channel).
  * Void: Zuber-Findlay drift-flux with Ishii C0 and Vgj (round-tube
    correlations, used here for a rod bundle in a pressure tube).
  * Friction: homogeneous two-phase model, Blasius factor, McAdams viscosity;
    no spacers, inlet/exit losses or rod roughness.
  * Wall heat transfer: wall temperature is the LOWER of (a) single-phase
    Dittus-Boelter convection and (b) the Thom nucleate-boiling curve
    (valid roughly 5 to 14 MPa). Dittus-Boelter is formally valid for Re > 1e4;
    the code counts nodes below that.
  * Conduction: constant conductivities, steady, uniform heat generation in a
    solid fuel cylinder, no radial power depression, all heat generated in the
    fuel (clad and coolant direct heating neglected). Gap = one conductance.
  * CHF: Biasi correlation (tube data), local-conditions form, evaluated at the
    local equilibrium quality; for subcooled nodes quality is clamped to 0
    (conservative). Constants: Biasi et al., Energia Nucleare 14 (1967) 530-536;
    original constants take G in g/cm2s (converted inside); the SI form was
    cross-checked against a secondary source, the original paper was not read.
    Ranges: 2.7-14 MPa, 100-6000 kg/m2s, D 3-37.5 mm, L 0.2-6 m. Rod-bundle effects are covered only by the user-set
    limits.chf_derating_factor multiplier (an assumption, not a result).
  * Hot channel: radial_peaking multiplies the channel power while the mass
    flow stays at its average-channel value (common-header assumption).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

import props

K0 = 273.15
G_ACC = 9.80665


# ---------------------------------------------------------------- power shape
def axial_layers_from_openmc(csv_path):
    """Sum kappa-fission over x,y for each z layer; keep the layers with power."""
    d = pd.read_csv(csv_path, skiprows=[1])
    d.columns = ["x", "y", "z", "nuclide", "score", "mean", "std"]
    k = d[d.score == "kappa-fission"]
    per_z = k.groupby("z")["mean"].sum().sort_index()
    active = per_z[per_z > 0]
    return active.values / active.values.mean()   # relative power, mean = 1


def build_profile(cfg, n_nodes):
    """Relative axial power f(z) at node centres, normalised so mean(f) = 1."""
    shape = cfg["power_shape"]["axial_shape"]
    zc = (np.arange(n_nodes) + 0.5) / n_nodes     # 0..1 along the heated length
    if shape == "cosine":
        f = np.cos(np.pi * (zc - 0.5))            # peak mid-height
        f = f / f.mean()
        print(f"  cosine shape: peak/mean = {f.max():.3f}")
    else:
        base_dir = os.path.dirname(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "config.yaml"))
        csv = os.path.join(base_dir, cfg["power_shape"]["power_map_csv"])
        layers = axial_layers_from_openmc(csv)
        lc = (np.arange(len(layers)) + 0.5) / len(layers)
        f = np.interp(zc, lc, layers)             # flat beyond first/last layer centre
        print(f"  OpenMC: {len(layers)} active z-layers, layer peak/mean = {layers.max():.3f}")
    return f / f.mean()


# ------------------------------------------------------- two-phase hydraulics
def darcy_f(Re):
    """Darcy friction factor (smooth tube)."""
    if Re < 2300.0:
        return 64.0 / Re
    if Re < 2.0e4:
        return 0.316 * Re ** -0.25
    return 0.184 * Re ** -0.2


def local_state(P, h, G):
    """Fluid state at pressure P [MPa], enthalpy h [kJ/kg], mass flux G [kg/m2s].

    Returns dict with x, alpha, rho_m (void-based mixture density), rho_h
    (homogeneous density), mu (friction viscosity), M (momentum flux, Pa),
    T_C (bulk temperature), Tsat_C, and liquid-phase k_l, cp_l, mu_l used for
    the single-phase heat transfer coefficient.
    """
    s = props.sat(P)
    if h < s["hf"]:                                # subcooled liquid
        T = props.T_from_h_liquid(P, h)
        w = props.liquid_from_T(P, T)
        return dict(x=0.0, alpha=0.0, rho_m=w["rho"], rho_h=w["rho"], mu=w["mu"],
                    M=G ** 2 / w["rho"], T_C=T - K0, Tsat_C=s["Tsat"] - K0,
                    k_l=w["k"], cp_l=w["cp"], mu_l=w["mu"])
    x = (h - s["hf"]) / s["hfg"]
    rf, rg = s["rho_f"], s["rho_g"]
    j = G * (x / rg + (1.0 - x) / rf)              # total volumetric flux
    jg = G * x / rg
    C0 = 1.2 - 0.2 * np.sqrt(rg / rf)              # Ishii distribution parameter
    Vgj = 1.41 * (s["sigma"] * G_ACC * (rf - rg) / rf ** 2) ** 0.25
    alpha = jg / (C0 * j + Vgj)
    alpha = min(max(alpha, 0.0), 0.999)
    rho_m = alpha * rg + (1.0 - alpha) * rf
    rho_h = 1.0 / (x / rg + (1.0 - x) / rf)
    mu = 1.0 / (x / s["mu_g"] + (1.0 - x) / s["mu_f"])        # McAdams
    M = G ** 2 * (x ** 2 / (alpha * rg) + (1.0 - x) ** 2 / ((1.0 - alpha) * rf))
    return dict(x=x, alpha=alpha, rho_m=rho_m, rho_h=rho_h, mu=mu, M=M,
                T_C=s["Tsat"] - K0, Tsat_C=s["Tsat"] - K0,
                k_l=s["k_f"], cp_l=s["cp_f"], mu_l=s["mu_f"])


def wall_temperature(st, P, G, Dh, q_wall):
    """Clad outer-surface temperature [C] for wall heat flux q_wall [W/m2].

    Returns (T_wall_C, h_conv, regime) with regime 'liquid' or 'boiling'.
    """
    Re = G * Dh / st["mu_l"]
    Pr = st["cp_l"] * 1000.0 * st["mu_l"] / st["k_l"]
    Nu = max(0.023 * Re ** 0.8 * Pr ** 0.4, 4.36)          # Dittus-Boelter, laminar floor
    h_conv = Nu * st["k_l"] / Dh
    T_sp = st["T_C"] + q_wall / h_conv                     # single-phase wall temperature
    # Thom: dT_sat = 22.65 * q''[MW/m2]^0.5 * exp(-p[bar]/87)
    T_nb = st["Tsat_C"] + 22.65 * np.sqrt(q_wall / 1.0e6) * np.exp(-(P * 10.0) / 87.0)
    if T_nb < T_sp:
        return T_nb, h_conv, "boiling", Re
    return T_sp, h_conv, "liquid", Re


def biasi_chf(P_MPa, G, Dh, x_e):
    """Biasi CHF [W/m2]. P in MPa, G in kg/m2s, Dh in m, x_e = equilibrium quality.

    q1, q2 in W/cm2, D in cm, G given in kg/m2s (converted to g/cm2s inside), p in bar;
    the larger of the two is taken. Subcooled quality is clamped to 0.
    """
    p = P_MPa * 10.0
    D = Dh * 100.0
    n = 0.4 if D >= 1.0 else 0.6
    x = min(max(x_e, 0.0), 1.0)
    fP = 0.7249 + 0.099 * p * np.exp(-0.032 * p)
    hP = -1.159 + 0.149 * p * np.exp(-0.019 * p) + 8.99 * p / (10.0 + p ** 2)
    # Biasi's constants (1883, 3780) are for G in g/cm2s; 1 kg/m2s = 0.1 g/cm2s.
    # (Equivalent SI form: 2.764e7 and 15.048e7 W/m2 with G in kg/m2s.)
    Gc = G / 10.0
    q1 = 1883.0 / (D ** n * Gc ** (1.0 / 6.0)) * (fP / Gc ** (1.0 / 6.0) - x)
    q2 = 3780.0 * hP / (D ** n * Gc ** 0.6) * (1.0 - x)
    return max(q1, q2) * 1.0e4


# ------------------------------------------------------------------- main
def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)

    op, geo, ps = cfg["operating_point"], cfg["geometry"], cfg["power_shape"]
    mat, lim = cfg["materials"], cfg["limits"]
    P = float(op["pressure_MPa"])
    T_in = float(op["inlet_temp_C"]) + K0
    n_ch = int(op["n_channels"])
    Q_ch_kW = float(op["channel_power_MW"]) * 1000.0 / n_ch
    frac = ps["fraction_to_coolant"]
    if frac is None:
        raise SystemExit("Set power_shape.fraction_to_coolant in config.yaml")
    # frac = share of the heat generated in the fuel that stays in the coolant.
    # The rest (assumed to go from the water into the D2O) does NOT change the
    # fuel/clad heat flux, only the coolant enthalpy rise.
    rp = ps.get("radial_peaking")
    rp = 1.0 if rp is None else float(rp)          # hot-channel power multiplier
    Q_avg_kW = Q_ch_kW * float(frac)               # average-channel NET heat kept by coolant
    Q_fuel_kW = Q_ch_kW * rp                       # heat generated in the fuel in this run
    Q_cool_kW = Q_avg_kW * rp                      # net heat kept by coolant in this run
    L = float(geo["heated_length_m"])
    Dh = float(geo["hydraulic_diameter_m"])
    n_rods = int(geo["n_rods"])
    r_fo = float(geo["fuel_outer_radius_m"])
    r_ci = float(geo["clad_inner_radius_m"])
    r_co = float(geo["clad_outer_radius_m"])
    h_gap = float(geo["gap_conductance_W_m2K"])
    k_clad = float(mat["clad_conductivity_W_mK"])
    k_fuel = float(mat["fuel_conductivity_W_mK"])
    T_lim = float(lim["fuel_temp_limit_C"])
    T_lim_strict = float(lim.get("fuel_temp_limit_strict_C", 600.0))
    T_int_lim = float(lim.get("interface_temp_limit_C", 550.0))
    chf_limit = float(lim["min_chf_margin"])
    chf_derate = float(lim.get("chf_derating_factor", 1.0))
    N = int(cfg["numerics"]["n_axial_nodes"])
    dz = L / N
    x_target = float(op["target_exit_quality"])

    print("== Inputs ==")
    f = build_profile(cfg, N)
    perim_calc = n_rods * 2.0 * np.pi * r_co
    if abs(perim_calc - float(geo["heated_perimeter_m"])) > 0.01 * perim_calc:
        print(f"  WARNING: heated_perimeter_m ({geo['heated_perimeter_m']}) differs from "
              f"n_rods*2*pi*r_clad_out ({perim_calc:.4f})")
    if rp != 1.0:
        print(f"  HOT-CHANNEL RUN: power x {rp:.3f}, mass flow kept at the average-channel value")

    # --- water states at inlet and target exit (at the inlet pressure)
    s = props.sat(P)
    w_in = props.liquid_from_T(P, T_in)
    h_in = w_in["h"]
    h_out = s["hf"] + x_target * s["hfg"]
    subcool_K = s["Tsat"] - T_in

    # --- mass flow from the AVERAGE-channel energy balance (design choice via target quality)
    mdot = Q_avg_kW / (h_out - h_in)              # kg/s  (kW / (kJ/kg))
    G = mdot / float(geo["flow_area_m2"])         # kg/m2s
    q_flux = Q_fuel_kW * 1000.0 / (float(geo["heated_perimeter_m"]) * L)  # W/m2 (fuel heat crosses the clad)

    # --- layer 1: march enthalpy up the channel
    z = (np.arange(N) + 0.5) * dz
    dQ = Q_cool_kW * f * dz / L                   # kW kept by the coolant in each slice
    dQf = Q_fuel_kW * f * dz / L                  # kW generated in the fuel in each slice
    h = np.empty(N + 1)
    h[0] = h_in
    for i in range(N):
        h[i + 1] = h[i] + dQ[i] / mdot
    h_node = 0.5 * (h[:-1] + h[1:])               # enthalpy at slice centre

    # --- layer 2: void fraction and pressure drop (explicit pressure march)
    Pn = P
    M_prev = G ** 2 / w_in["rho"]
    T = np.empty(N); xq = np.empty(N); alpha = np.empty(N); rho_m = np.empty(N)
    P_end = np.empty(N); P_start_arr = np.empty(N)
    dp_f = np.empty(N); dp_g = np.empty(N); dp_a = np.empty(N)
    states = []
    for i in range(N):
        P_start = Pn
        P_start_arr[i] = P_start
        st = local_state(P_start, h_node[i], G)
        states.append(st)
        Re = G * Dh / st["mu"]
        dp_f[i] = darcy_f(Re) / Dh * G ** 2 / (2.0 * st["rho_h"]) * dz
        dp_g[i] = st["rho_m"] * G_ACC * dz
        dp_a[i] = st["M"] - M_prev
        M_prev = st["M"]
        Pn -= (dp_f[i] + dp_g[i] + dp_a[i]) / 1.0e6
        P_end[i] = Pn
        T[i], alpha[i], rho_m[i] = st["T_C"], st["alpha"], st["rho_m"]
        xq[i] = props.quality_from_h(P_start, h_node[i])   # negative = subcooled
    # Close the last half-slice so exit void and acceleration refer to the true exit.
    st_exit = local_state(Pn, h[-1], G)
    corr = st_exit["M"] - M_prev
    dp_a[-1] += corr
    Pn -= corr / 1.0e6
    P_end[-1] = Pn
    alpha_exit = st_exit["alpha"]
    x_exit_nominal = props.quality_from_h(P, h[-1])
    x_exit_actual = props.quality_from_h(Pn, h[-1])

    idx = np.where(h >= s["hf"])[0]
    z_sat = (idx[0] * dz) if len(idx) else float("nan")

    # --- layer 3: wall, clad and fuel temperatures
    q_lin = np.empty(N); q_wall = np.empty(N); h_conv = np.empty(N); Re_arr = np.empty(N)
    T_wall = np.empty(N); T_ci = np.empty(N); T_fs = np.empty(N); T_cl = np.empty(N)
    d_clad = np.empty(N); d_gap = np.empty(N); d_fuel = np.empty(N)
    regime = []
    for i in range(N):
        q_lin[i] = dQf[i] * 1000.0 / (dz * n_rods)              # W/m per rod
        q_wall[i] = q_lin[i] / (2.0 * np.pi * r_co)            # W/m2 at clad outer surface
        q_fuel_surf = q_lin[i] / (2.0 * np.pi * r_fo)          # W/m2 at fuel surface
        T_wall[i], h_conv[i], reg, Re_arr[i] = wall_temperature(states[i], P_start_arr[i], G, Dh, q_wall[i])
        regime.append(reg)
        d_clad[i] = q_lin[i] * np.log(r_co / r_ci) / (2.0 * np.pi * k_clad)
        d_gap[i] = q_fuel_surf / h_gap
        d_fuel[i] = q_lin[i] / (4.0 * np.pi * k_fuel)
        T_ci[i] = T_wall[i] + d_clad[i]
        T_fs[i] = T_ci[i] + d_gap[i]
        T_cl[i] = T_fs[i] + d_fuel[i]
    # --- layer 4: CHF margin
    chf = np.array([biasi_chf(P_start_arr[i], G, Dh, xq[i]) for i in range(N)])
    mchf_raw = chf / q_wall
    mchf = chf_derate * mchf_raw
    i_c = int(np.argmin(mchf))
    i_pk = int(np.argmax(T_cl))
    boil_idx = [i for i, r in enumerate(regime) if r == "boiling"]
    n_lowRe = int(np.sum(Re_arr < 1.0e4))

    # --- checks and report
    heat_sum = dQ.sum()
    heat_abs = mdot * (h[-1] - h[0])
    print("\n== Layer 1 results: energy balance ==")
    print(f"  channel power to coolant : {Q_cool_kW:9.2f} kW "
          f"({frac*100:.0f}% of {Q_ch_kW*rp:.1f} kW fuel heat; fraction is an ASSUMPTION)")
    print(f"  heat generated in fuel   : {Q_fuel_kW:9.2f} kW (sets fuel/clad temperatures)")
    print(f"  inlet T / Tsat / subcool : {T_in-K0:6.1f} C / {s['Tsat']-K0:6.1f} C / {subcool_K:6.1f} K")
    print(f"  mass flow                : {mdot:9.4f} kg/s")
    print(f"  mass flux G              : {G:9.1f} kg/m2s")
    print(f"  mean wall heat flux      : {q_flux/1000:9.1f} kW/m2")
    print(f"  bulk saturation reached  : z = {z_sat*1000:6.1f} mm of {L*1000:.0f} mm")
    print(f"  exit quality (design P)  : {x_exit_nominal:9.4f} (target {x_target}"
          f"{'' if rp == 1.0 else ', hot channel is expected to exceed it'})")

    dp_tot = dp_f.sum() + dp_g.sum() + dp_a.sum()
    print("\n== Layer 2 results: void fraction and pressure drop ==")
    print(f"  exit void fraction       : {alpha_exit:9.3f}")
    print(f"  exit quality (exit P)    : {x_exit_actual:9.4f}")
    print(f"  pressure drop total      : {dp_tot/1000:9.3f} kPa")
    print(f"    friction               : {dp_f.sum()/1000:9.3f} kPa")
    print(f"    gravity                : {dp_g.sum()/1000:9.3f} kPa")
    print(f"    acceleration           : {dp_a.sum()/1000:9.3f} kPa")
    print(f"  exit pressure            : {Pn:9.4f} MPa (inlet {P} MPa)")

    print("\n== Layer 3 results: wall, clad and fuel temperatures ==")
    print(f"  peak linear power        : {q_lin.max()/1000:9.2f} kW/m per rod "
          f"(mean {q_lin.mean()/1000:.2f}); peak wall flux {q_wall.max()/1000:.0f} kW/m2")
    if boil_idx:
        print(f"  wall boiling (Thom) from : z = {z[boil_idx[0]]*1000:6.1f} mm "
              f"({len(boil_idx)} of {N} nodes boiling)")
    else:
        print("  wall boiling             : none (single-phase convection everywhere)")
    print(f"  max clad outer (wall) T  : {T_wall.max():9.1f} C")
    print(f"  max clad inner T         : {T_ci.max():9.1f} C")
    print(f"  max fuel surface T       : {T_fs.max():9.1f} C")
    print(f"  max fuel centerline T    : {T_cl.max():9.1f} C at z = {z[i_pk]*1000:.0f} mm")
    print(f"  drops at that height     : clad {d_clad[i_pk]:.0f} K, gap {d_gap[i_pk]:.0f} K, "
          f"fuel {d_fuel[i_pk]:.0f} K, wall-to-saturation {T_wall[i_pk]-states[i_pk]['Tsat_C']:.1f} K")
    margin = T_lim - T_cl.max()
    print(f"  fuel temperature limit   : {T_lim:.0f} C -> margin {margin:+.1f} K "
          f"[{'OK' if margin >= 0 else 'LIMIT EXCEEDED'}]")
    margin_strict = T_lim_strict - T_cl.max()
    print(f"  strict fuel limit        : {T_lim_strict:.0f} C -> margin {margin_strict:+.1f} K "
          f"[{'OK' if margin_strict >= 0 else 'LIMIT EXCEEDED'}]")
    margin_int = T_int_lim - T_ci.max()
    print(f"  fuel-clad interface limit: {T_int_lim:.0f} C vs max clad inner {T_ci.max():.1f} C -> margin {margin_int:+.1f} K "
          f"[{'OK' if margin_int >= 0 else 'LIMIT EXCEEDED'}] (interface taken as clad inner surface)")
    print(f"  Dittus-Boelter validity  : {n_lowRe} of {N} nodes have Re < 1e4 "
          f"(min Re {Re_arr.min():.0f}); treat the liquid-region wall temperature as approximate")

    print("\n== Layer 4 results: CHF margin (Biasi) ==")
    print(f"  CHF range along channel  : {chf.min()/1e6:6.2f} to {chf.max()/1e6:.2f} MW/m2")
    print(f"  minimum margin, raw Biasi: {mchf_raw.min():9.2f} at z = {z[int(np.argmin(mchf_raw))]*1000:.0f} mm")
    print(f"  with derating factor {chf_derate:.2f}: {mchf.min():9.2f} at z = {z[i_c]*1000:.0f} mm")
    print(f"  CHF margin limit         : {chf_limit:.2f} -> "
          f"[{'OK' if mchf.min() >= chf_limit else 'BELOW LIMIT'}] (judged on the derated value)")
    rng = [("pressure 2.7-14 MPa", 2.7 <= P_start_arr.min() and P_start_arr.max() <= 14.0),
           ("mass flux 100-6000 kg/m2s", 100.0 <= G <= 6000.0),
           ("diameter 3-37.5 mm", 3.0 <= Dh * 1000 <= 37.5),
           ("length 0.2-6 m", 0.2 <= L <= 6.0)]
    print("  Biasi validity: " + "; ".join(f"{n}: {'inside' if ok else 'OUTSIDE'}" for n, ok in rng))
    x_edge = 1.0 / (1.0 + s["rho_f"] / s["rho_g"])
    print(f"  Biasi quality range starts near x = {x_edge:.3f}; exit quality here = {x_exit_actual:.3f}"
          f"{' (below the range edge: extrapolation)' if x_exit_actual < x_edge else ''}")
    if Dh * 1000 < 5.0:
        print("  note: hydraulic diameter is near the lower edge of the Biasi range")

    print("\n== Checks ==")
    print(f"  sum of slice heat = {heat_sum:.3f} kW, channel heat to coolant = {Q_cool_kW:.3f} kW")
    print(f"  mdot*(h_out-h_in) = {heat_abs:.3f} kW")
    ok1 = abs(heat_sum - Q_cool_kW) < 1e-6 * Q_cool_kW and abs(heat_abs - Q_cool_kW) < 1e-6 * Q_cool_kW
    if rp == 1.0:
        ok1 = ok1 and abs(x_exit_nominal - x_target) < 1e-6
    print("  ENERGY BALANCE CHECK:", "PASSED" if ok1 else "FAILED - send me this output")
    ok2 = bool(np.all(alpha >= 0) and np.all(alpha <= 1) and np.all(np.diff(alpha) >= -1e-9)
               and np.isfinite(dp_tot) and np.all(P_end > 0))
    print("  VOID/PRESSURE SANITY CHECK (0<=alpha<=1, alpha rising, finite dp):",
          "PASSED" if ok2 else "FAILED - send me this output")
    ok3 = bool(np.all(T_cl >= T_fs) and np.all(T_fs >= T_ci) and np.all(T_ci >= T_wall)
               and np.all(T_wall >= T - 1e-9) and abs(q_lin.sum() * dz * n_rods / 1000.0 - Q_fuel_kW) < 1e-6 * Q_fuel_kW)
    print("  TEMPERATURE CHAIN CHECK (centre >= fuel surf >= clad in >= wall >= bulk, rod power sums to channel power):",
          "PASSED" if ok3 else "FAILED - send me this output")

    ok4 = bool(np.all(np.isfinite(chf)) and np.all(chf > 0))
    print("  CHF SANITY CHECK (finite, positive):", "PASSED" if ok4 else "FAILED - send me this output")

    # --- save
    summary = dict(
        mdot_kg_s=float(mdot), G_kg_m2s=float(G), exit_void=float(alpha_exit),
        exit_quality=float(x_exit_actual), dp_total_kPa=float(dp_tot / 1000.0),
        T_wall_max_C=float(T_wall.max()), T_fuel_centre_max_C=float(T_cl.max()),
        z_T_fuel_max_mm=float(z[i_pk] * 1000.0), fuel_limit_C=T_lim,
        fuel_margin_K=float(margin), fuel_limit_strict_C=T_lim_strict,
        fuel_margin_strict_K=float(margin_strict), T_clad_inner_max_C=float(T_ci.max()),
        interface_limit_C=T_int_lim, interface_margin_K=float(margin_int),
        Q_fuel_kW=float(Q_fuel_kW), Q_to_coolant_kW=float(Q_cool_kW), chf_margin_raw_min=float(mchf_raw.min()),
        chf_margin_derated_min=float(mchf.min()), chf_limit=chf_limit,
        chf_derating_factor=chf_derate, z_chf_min_mm=float(z[i_c] * 1000.0),
        all_checks_passed=bool(ok1 and ok2 and ok3 and ok4))
    import json
    with open("summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    pd.DataFrame({"z_m": z, "rel_power": f, "slice_kW": dQ, "h_kJ_kg": h_node,
                  "T_water_C": T, "quality": xq, "void_fraction": alpha,
                  "P_MPa_end_of_node": P_end, "rho_mix_kg_m3": rho_m,
                  "dp_friction_Pa": dp_f, "dp_gravity_Pa": dp_g, "dp_accel_Pa": dp_a,
                  "q_lin_W_m_per_rod": q_lin, "q_wall_W_m2": q_wall, "h_conv_W_m2K": h_conv,
                  "Re_liquid": Re_arr, "wall_regime": regime,
                  "T_wall_C": T_wall, "T_clad_inner_C": T_ci,
                  "T_fuel_surface_C": T_fs, "T_fuel_centre_C": T_cl,
                  "CHF_W_m2": chf, "CHF_margin_raw": mchf_raw,
                  "CHF_margin_derated": mchf}).to_csv("channel_results.csv", index=False)

    zmm = z * 1000
    if os.environ.get("CHANNEL_NOPLOT"):
        print("\nSaved channel_results.csv (plot skipped)")
        return
    fig, ax = plt.subplots(4, 3, figsize=(13, 16))
    ax[0, 0].plot(f, zmm); ax[0, 0].set_xlabel("relative power"); ax[0, 0].set_title("Axial power shape")
    ax[0, 1].plot(T, zmm); ax[0, 1].set_xlabel("water temperature (C)"); ax[0, 1].set_title("Bulk water temperature")
    ax[0, 2].plot(np.clip(xq, 0, None) * 100, zmm, label="quality (%)")
    ax[0, 2].plot(alpha * 100, zmm, label="void fraction (%)")
    ax[0, 2].set_xlabel("percent"); ax[0, 2].set_title("Quality and void fraction"); ax[0, 2].legend()
    ax[1, 0].plot((P - P_end) * 1000, zmm); ax[1, 0].set_xlabel("pressure drop from inlet (kPa)")
    ax[1, 0].set_title("Pressure drop along channel")
    ax[1, 1].plot(np.cumsum(dp_f) / 1000, zmm, label="friction")
    ax[1, 1].plot(np.cumsum(dp_g) / 1000, zmm, label="gravity")
    ax[1, 1].plot(np.cumsum(dp_a) / 1000, zmm, label="acceleration")
    ax[1, 1].set_xlabel("cumulative (kPa)"); ax[1, 1].set_title("Pressure drop components"); ax[1, 1].legend()
    ax[1, 2].plot(rho_m, zmm); ax[1, 2].set_xlabel("mixture density (kg/m3)"); ax[1, 2].set_title("Coolant density")
    ax[2, 0].plot(T, zmm, label="bulk water")
    ax[2, 0].plot(T_wall, zmm, label="clad outer")
    ax[2, 0].plot(T_ci, zmm, label="clad inner")
    ax[2, 0].plot(T_fs, zmm, label="fuel surface")
    ax[2, 0].plot(T_cl, zmm, label="fuel centre")
    ax[2, 0].axvline(T_lim, color="k", ls="--", label=f"fuel limit {T_lim:.0f} C")
    ax[2, 0].axvline(T_lim_strict, color="r", ls=":", label=f"strict limit {T_lim_strict:.0f} C")
    ax[2, 0].axvline(T_int_lim, color="m", ls="-.", label=f"interface limit {T_int_lim:.0f} C")
    ax[2, 0].set_xlabel("temperature (C)"); ax[2, 0].set_title("Temperatures through the rod"); ax[2, 0].legend(fontsize=7)
    ax[2, 1].plot(d_clad, zmm, label="across clad")
    ax[2, 1].plot(d_gap, zmm, label="across gap")
    ax[2, 1].plot(d_fuel, zmm, label="through fuel")
    ax[2, 1].set_xlabel("temperature rise (K)"); ax[2, 1].set_title("Where the temperature rise happens"); ax[2, 1].legend()
    ax[2, 2].plot(q_wall / 1000, zmm); ax[2, 2].set_xlabel("wall heat flux (kW/m2)"); ax[2, 2].set_title("Clad surface heat flux")
    ax[3, 0].plot(chf / 1000, zmm, label="CHF (Biasi)")
    ax[3, 0].plot(q_wall / 1000, zmm, label="wall heat flux")
    ax[3, 0].set_xlabel("kW/m2"); ax[3, 0].set_title("CHF versus actual heat flux"); ax[3, 0].legend()
    ax[3, 1].plot(mchf_raw, zmm, label="raw Biasi")
    ax[3, 1].plot(mchf, zmm, label=f"derated x{chf_derate:.2f}")
    ax[3, 1].axvline(chf_limit, color="k", ls="--", label=f"limit {chf_limit:.2f}")
    ax[3, 1].set_xlabel("CHF margin"); ax[3, 1].set_title("CHF margin"); ax[3, 1].legend(fontsize=8)
    ax[3, 2].axis("off")
    for a in ax.ravel():
        if a.axison:
            a.grid(alpha=0.3); a.set_ylabel("height (mm)")
    fig.tight_layout(); fig.savefig("channel_energy_balance.png", dpi=130)
    print("\nSaved channel_results.csv and channel_energy_balance.png")


if __name__ == "__main__":
    main()
