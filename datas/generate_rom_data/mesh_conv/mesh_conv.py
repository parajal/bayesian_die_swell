import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

folder = Path(__file__).parent

for i in range(4):
    x = np.loadtxt(folder / f"curve4_x_m{i}.txt")
    y = np.loadtxt(folder / f"curve4_m{i}.txt")

    plt.plot(x, y, label=f"m{i}")

plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.tight_layout()
plt.show()