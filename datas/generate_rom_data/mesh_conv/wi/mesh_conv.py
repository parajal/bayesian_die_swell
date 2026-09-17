import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

folder = Path(__file__).parent

wi = [0.5, 1, 2, 4]
s_r_tan = [1.133, 1.142,1.177,1.277]
s_r_fem = [1.136, 1.185, 1.342, 1.695]
plt.plot(wi, s_r_tan, marker = 'o', label = "tanner")
plt.plot(wi, s_r_fem, marker = 'o', label = "FEM")

plt.xlabel("Wi")
plt.ylabel("Swell_ratio")
plt.legend()
plt.show()