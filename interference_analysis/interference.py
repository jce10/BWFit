import numpy as np
import matplotlib.pyplot as plt

# Representative trapped/superradiant pair
M_T = 7.688   # trapped centroid, MeV
G_T = 0.08   # trapped width, MeV
f_T = 1.0    # trapped strength scale

M_S = 8.220  # superradiant centroid, MeV
G_S = 1.20   # superradiant width, MeV
f_S = 2.0    # superradiant strength scale

E = np.linspace(4.8, 10.6, 1200)

def interference(E, M_T, G_T, f_T, M_S, G_S, f_S, delta_phys):
    X = (E - M_T) * (E - M_S) + (G_T * G_S) / 4.0
    Y = (G_T / 2.0) * (E - M_S) - (G_S / 2.0) * (E - M_T)
    prefactor = (2.0 * np.sqrt(f_T * f_S * G_T * G_S)) / (2.0 * np.pi)
    return prefactor * (X * np.cos(delta_phys) + Y * np.sin(delta_phys)) / (X**2 + Y**2)

delta_fit_values = [0, 0.5, 1.0, 1.5, 2.0]

plt.figure(figsize=(9, 5.5))

for dfit in delta_fit_values:
    dphys = dfit * np.pi / 2.0
    y = interference(E, M_T, G_T, f_T, M_S, G_S, f_S, dphys)
    plt.plot(E, y, label=fr"$\delta_{{fit}}={dfit}$")

plt.axhline(0, linestyle="--", linewidth=1)
plt.axvline(M_T, linestyle=":", linewidth=1, label="trapped centroid")
plt.axvline(M_S, linestyle="-.", linewidth=1, label="superradiant centroid")

plt.xlabel("Excitation energy E [MeV]")
plt.ylabel("Interference contribution [arb. units]")
plt.title(r"Trapped–Superradiant Interference Term vs Phase, $\delta=\delta_{fit}\pi/2$")
plt.legend(fontsize=9)
plt.tight_layout()
plt.show()
