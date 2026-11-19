"""Left: void reactivity effect (0 -> 20 % void). Right: how much of the heated channel length boils.
Fill in the two BOILING numbers from your 1D thermal-hydraulics run, then: python make_void_chart.py"""
from math import sqrt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- void results: (k at 0 %, std, k at 20 %, std) ----------
D = {
    "Your design (10 MPa)": (1.04745, 0.00069, 1.05163, 0.00072),
    "RBMK cell (7 MPa)":    (1.31174, 0.00087, 1.31472, 0.00095),
}
COL = {"Your design (10 MPa)": "#2a6fbb", "RBMK cell (7 MPa)": "#d9782d"}

# ---------- boiling length from the T-H run (FILL THESE IN) ----------
HEATED_LENGTH_M = 0.40
BOILING_LENGTH_M = {            # length of channel where boiling occurs, metres
    "Average channel": 0.035,    # e.g. 0.18
    "Hot channel": 0.093,        # e.g. 0.25
}
BOIL_COL = {"Average channel": "#2a6fbb", "Hot channel": "#a23b3b"}
BOIL_CASE = "10 MPa, 0.2 MW, gap 15000 W/m2K"   # label the case you used

if any(v is None for v in BOILING_LENGTH_M.values()):
    raise SystemExit("Fill in BOILING_LENGTH_M (values in metres) from your T-H run first.")

fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.5, 4.8))

for j, (name, (k0, s0, k1, s1)) in enumerate(D.items()):
    d = (1 / k0 - 1 / k1) * 1e5
    e = 1e5 * sqrt((s0 / k0**2) ** 2 + (s1 / k1**2) ** 2)
    ax.bar(j, d, width=0.55, color=COL[name], zorder=2)
    ax.errorbar(j, d, yerr=e, color="#222", capsize=5, lw=1.5, zorder=3)
    ax.annotate(f"{d:,.0f} ± {e:,.0f}", (j, d + e), xytext=(0, 6), textcoords="offset points",
                ha="center", va="bottom", fontsize=10, color="#222")
ax.set_xticks(range(len(D)))
ax.set_xticklabels(list(D), fontsize=10)
ax.set_ylabel("pcm")
ax.set_title("Reactivity change, 0 → 20 % void", fontsize=11)
ax.set_ylim(0, ax.get_ylim()[1] * 1.15)

for j, (name, L) in enumerate(BOILING_LENGTH_M.items()):
    pct = 100.0 * L / HEATED_LENGTH_M
    bx.bar(j, pct, width=0.55, color=BOIL_COL[name], zorder=2)
    bx.annotate(f"{pct:.0f} %  ({L:.2f} of {HEATED_LENGTH_M:.2f} m)", (j, pct), xytext=(0, 6),
                textcoords="offset points", ha="center", va="bottom", fontsize=10, color="#222")
bx.set_xticks(range(len(BOILING_LENGTH_M)))
bx.set_xticklabels(list(BOILING_LENGTH_M), fontsize=10)
bx.set_ylabel("% of heated length")
bx.set_ylim(0, 115)
bx.set_title("Share of channel length with boiling\n" + BOIL_CASE, fontsize=11)

for a in (ax, bx):
    a.axhline(0, color="#888", lw=0.8)
    a.grid(axis="y", alpha=0.25, zorder=0)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)

fig.suptitle("Void reactivity effect and extent of boiling", fontsize=12)
fig.text(0.5, 0.005, "Provisional: void results are a single fresh-fuel cell, no control rods, 1σ error bars from k only; "
         "boiling extent is from the 1D T-H model.", ha="center", fontsize=8, color="#555")
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig("void_effect_chart.png", dpi=200, bbox_inches="tight")
print("saved void_effect_chart.png")
