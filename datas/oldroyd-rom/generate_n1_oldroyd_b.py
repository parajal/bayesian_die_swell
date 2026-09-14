"""Generate N1.txt from parameters.txt using the Oldroyd-B model.

For each row of parameters.txt (columns: lambda, beta, alpha) this computes
N1 = tau_xx - tau_yy at the axisymmetric wall shear rate for U_avg / radius,
and writes one N1 value per line to N1.txt (single column).

The third column (alpha) is ignored: Oldroyd-B is the alpha = 0 special case.

Run:
    python generate_n1_oldroyd_b.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# The rheology module lives in the sibling swell_inference-2 repo. Add it to the
# path only if it is not already importable from the current environment.
try:
    import rheology  # noqa: F401
except ModuleNotFoundError:
    fallback = Path(r"C:/Users/paraj/Documents/swell_inference-2")
    if (fallback / "rheology.py").exists():
        sys.path.insert(0, str(fallback))

from compute_n1_oldroyd_b_function import compute_n1_oldroyd_b

# Flow / material settings (match compute_n1_uncertainty.py conventions).
U_AVG = 0.1
RADIUS = 1.0
ETA0 = 1.0

SCRIPT_DIR = Path(__file__).resolve().parent
PARAMETERS_PATH = SCRIPT_DIR / "parameters.txt"
OUTPUT_PATH = SCRIPT_DIR / "N1.txt"


def main() -> None:
    params = np.loadtxt(PARAMETERS_PATH)
    params = np.atleast_2d(params)

    n1_values = []
    for lam, beta, _alpha in params[:, :3]:
        result = compute_n1_oldroyd_b(
            U_avg=U_AVG,
            radius=RADIUS,
            lam=float(lam),
            beta=float(beta),
            eta0=ETA0,
        )
        n1_values.append(float(result["N1"][0]))

    np.savetxt(OUTPUT_PATH, np.array(n1_values), fmt="%.12g")
    print(f"Wrote {len(n1_values)} N1 values to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
