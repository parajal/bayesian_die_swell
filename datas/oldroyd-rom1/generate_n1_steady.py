"""Generate N1.txt: steady-state (fully-developed) N1 per case folder.

Each case folder holds ``curve5_0004.txt`` (the last time snapshot, i.e. the
steady state) sampled along the die wall. Its columns are

    x  y  pressure  c_xy  gammadot  c_xx  c_yy  N1

N1 = tau_xx - tau_yy is essentially constant in the fully-developed region
upstream of the die and then spikes near the die exit (x -> 0). To get the
steady-state value we average N1 over the upstream window x in [-3, -1]
(the same region used for the tau_xy pressure slope), which excludes the exit.

One N1 value per folder is written to N1.txt (single column), in sorted
folder order (case001 .. case125), row-aligned with parameters.txt,
pressure_all.txt, curve4_y_all.txt and tau_xy.txt.

Run:
    python generate_n1_steady.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

X_LO, X_HI = -3.0, -1.0  # upstream fully-developed window (excludes die exit)
N1_COL = 7               # zero-based column index of N1
X_COL = 0

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_PATH = SCRIPT_DIR / "N1.txt"
CURVE5_NAME = "curve5_0004.txt"


def steady_n1(path: Path) -> float:
    data = np.loadtxt(path, comments="#")
    x = data[:, X_COL]
    n1 = data[:, N1_COL]
    mask = (x >= X_LO) & (x <= X_HI)
    if not np.any(mask):
        raise ValueError(f"no points in x=[{X_LO},{X_HI}] for {path}")
    return float(np.mean(n1[mask]))


def main() -> None:
    folders = sorted(p for p in SCRIPT_DIR.glob("case*") if p.is_dir())
    if not folders:
        raise FileNotFoundError(f"no case* folders found in {SCRIPT_DIR}")

    n1_values = [steady_n1(folder / CURVE5_NAME) for folder in folders]

    np.savetxt(OUTPUT_PATH, np.array(n1_values), fmt="%.12g")
    print(f"Wrote {len(n1_values)} steady-state N1 values to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
