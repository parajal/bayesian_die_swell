import numpy as np
import matplotlib.pyplot as plt

# Read data, skipping comment lines starting with '#'
data = np.loadtxt("curve5_0010.txt", comments="#")

x        = data[:, 0]
pressure = data[:, 2]
c_xy     = data[:, 3]
gammadot = data[:, 4]
N1       = data[:, 7]

plt.figure(figsize=(9, 6))
plt.plot(x, pressure, label="pressure",  marker=".")
plt.plot(x, c_xy,     label="c_xy",      marker=".")
plt.plot(x, gammadot, label="gammadot",  marker=".")
plt.plot(x, N1,       label="N1",         marker=".")

plt.xlabel("x")
plt.ylabel("value")
# plt.title("x vs pressure, c_xy, gammadot, N1")
plt.legend()
plt.grid(True, alpha=0.3)
# plt.tight_layout()
plt.savefig("curve5_plot.png", dpi=150)
plt.show()
