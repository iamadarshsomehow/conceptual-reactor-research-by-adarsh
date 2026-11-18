"""props.py - water/steam properties for the 1D channel model.

Uses iapws (IAPWS-97). Install once in your venv:
    pip install iapws pyyaml numpy scipy matplotlib
Units used everywhere in this project:
    P [MPa], T [K], h [kJ/kg], rho [kg/m3], mu [Pa.s], k [W/m.K], cp [kJ/kg.K], sigma [N/m]
"""
from iapws import IAPWS97


def sat(P):
    """Saturated liquid and vapour properties at pressure P [MPa]."""
    f = IAPWS97(P=P, x=0)
    g = IAPWS97(P=P, x=1)
    return {
        "Tsat": f.T,
        "hf": f.h, "hg": g.h, "hfg": g.h - f.h,
        "rho_f": f.rho, "rho_g": g.rho,
        "mu_f": f.mu, "mu_g": g.mu,
        "k_f": f.k, "k_g": g.k,
        "cp_f": f.cp, "cp_g": g.cp,
        "sigma": f.sigma,
    }


def liquid_from_T(P, T):
    """Subcooled liquid properties at P [MPa], T [K]."""
    w = IAPWS97(P=P, T=T)
    return {"h": w.h, "rho": w.rho, "mu": w.mu, "k": w.k, "cp": w.cp}


def T_from_h_liquid(P, h):
    """Liquid temperature [K] from enthalpy h [kJ/kg] (valid below saturation)."""
    return IAPWS97(P=P, h=h).T


def quality_from_h(P, h):
    """Equilibrium quality x from enthalpy. Negative = subcooled, >1 = superheated."""
    s = sat(P)
    return (h - s["hf"]) / s["hfg"]


if __name__ == "__main__":
    # Self-test against steam-table values at 10 MPa (rounded table values).
    s = sat(10.0)
    expected = {            # value, tolerance
        "Tsat":  (584.15, 0.3),    # 311.0 C
        "hf":    (1407.9, 3.0),
        "hg":    (2725.5, 3.0),
        "rho_f": (688.4, 2.0),
        "rho_g": (55.5, 0.5),
    }
    ok = True
    for key, (val, tol) in expected.items():
        diff = s[key] - val
        flag = "OK " if abs(diff) <= tol else "BAD"
        if flag == "BAD":
            ok = False
        print(f"{flag} {key:6s} computed {s[key]:10.3f}  table {val:10.3f}  diff {diff:+.3f}")
    print("Tsat in C:", round(s["Tsat"] - 273.15, 2))
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED - send me this output")
