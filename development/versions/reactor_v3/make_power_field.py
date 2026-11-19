import openmc
import numpy as np

# ------------------------------------------------------------
# Load OpenMC result
# ------------------------------------------------------------
sp = openmc.StatePoint("statepoint.5.h5")
tally = sp.get_tally(name="V3_flux")

flux = tally.get_reshaped_data(value="mean")[:, 0, 0]
flux = flux.reshape((100, 100, 100))

# ------------------------------------------------------------
# Convert flux distribution into a normalized heat-source field
# ------------------------------------------------------------
# This is a transport-to-CFD demonstration.
# It is NOT an absolute reactor-power calculation.

mean_flux = np.mean(flux[flux > 0])

q_rel = flux / mean_flux

# ------------------------------------------------------------
# Save the field
# ------------------------------------------------------------
np.save("V3_power_field.npy", q_rel)

# Also save a simple text version for inspection
np.savetxt(
    "V3_power_field_flat.txt",
    q_rel.ravel()
)

print("==============================================")
print("OPENMC -> POWER FIELD")
print("==============================================")
print("Flux shape       :", flux.shape)
print("Mean positive flux:", mean_flux)
print("Relative q''' min :", q_rel.min())
print("Relative q''' max :", q_rel.max())
print("Relative q''' mean:", q_rel.mean())
print()
print("Saved:")
print("  V3_power_field.npy")
print("  V3_power_field_flat.txt")
