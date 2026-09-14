"""Generate tau_xy.txt from pressure_all.txt.

Each row of pressure_all.txt holds five centerline pressure values sampled at
x = [-3, -2.5, -2, -1.5, -1]. For each row we fit a line p(x) through those
five points (least squares), take the slope dp/dx, and compute the wall shear
stress

    tau_xy = -(1/2) * dp/dx        (R/2 factor with R = 1)

The slope is negated because pressure decreases with x (dp/dx < 0), so tau_xy
comes out positive, matching the analytic Oldroyd-B wall shear stress.

One tau_xy value per row is written to tau_xy.txt (single column).

Run:
    python generate_tau_xy.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

X_POINTS = np.array([-3.0, -2.5, -2.0, -1.5, -1.0])

SCRIPT_DIR = Path(__file__).resolve().parent
PRESSURE_PATH = SCRIPT_DIR / "pressure_all.txt"
OUTPUT_PATH = SCRIPT_DIR / "tau_xy.txt"


def main() -> None:
    pressures = np.loadtxt(PRESSURE_PATH)
    pressures = np.atleast_2d(pressures)

    if pressures.shape[1] != X_POINTS.size:
        raise ValueError(
            f"expected {X_POINTS.size} columns in pressure_all.txt, "
            f"found {pressures.shape[1]}"
        )

    tau_xy_values = []
    for row in pressures:
        slope = np.polyfit(X_POINTS, row, 1)[0]
        tau_xy_values.append(-0.5 * slope)

    np.savetxt(OUTPUT_PATH, np.array(tau_xy_values), fmt="%.12g")
    print(f"Wrote {len(tau_xy_values)} tau_xy values to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
