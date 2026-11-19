import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ------------------------------------------------------------
# Load statepoint
# ------------------------------------------------------------
sp = openmc.StatePoint("statepoint.5.h5")
tally = sp.get_tally(name="V3_flux")

flux = tally.get_reshaped_data(value="mean")[:, 0, 0]
flux = flux.reshape((100, 100, 100))

# ------------------------------------------------------------
# Geometry bounds from V3
# ------------------------------------------------------------
lower = np.array([-20.0, -23.0940094, 0.0])
upper = np.array([20.0, 23.0940094, 40.0])

# Cell-center coordinates
x = np.linspace(lower[0], upper[0], 100)
y = np.linspace(lower[1], upper[1], 100)
z = np.linspace(lower[2], upper[2], 100)

# ------------------------------------------------------------
# Central slices
# ------------------------------------------------------------
ix = 50
iy = 50
iz = 50

xy = flux[:, :, iz]
xz = flux[:, iy, :]
yz = flux[ix, :, :]

# Avoid log(0)
positive = flux[flux > 0]
vmin = positive.min()
vmax = positive.max()

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# XY
im1 = axes[0].imshow(
    xy.T,
    origin="lower",
    extent=[lower[0], upper[0], lower[1], upper[1]],
    aspect="equal",
    norm=LogNorm(vmin=vmin, vmax=vmax)
)
axes[0].set_title("XY slice — central Z")
axes[0].set_xlabel("X")
axes[0].set_ylabel("Y")

# XZ
im2 = axes[1].imshow(
    xz.T,
    origin="lower",
    extent=[lower[0], upper[0], lower[2], upper[2]],
    aspect="equal",
    norm=LogNorm(vmin=vmin, vmax=vmax)
)
axes[1].set_title("XZ slice — central Y")
axes[1].set_xlabel("X")
axes[1].set_ylabel("Z")

# YZ
im3 = axes[2].imshow(
    yz.T,
    origin="lower",
    extent=[lower[1], upper[1], lower[2], upper[2]],
    aspect="equal",
    norm=LogNorm(vmin=vmin, vmax=vmax)
)
axes[2].set_title("YZ slice — central X")
axes[2].set_xlabel("Y")
axes[2].set_ylabel("Z")

fig.colorbar(im3, ax=axes, label="Flux")

plt.tight_layout()
plt.savefig("V3_flux_highres.png", dpi=300, bbox_inches="tight")
plt.show()
