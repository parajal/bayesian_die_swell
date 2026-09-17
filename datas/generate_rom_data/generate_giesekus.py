"""
Generate final curve-4 surface rows for the Giesekus cluster solver.

Edit the constants below, then run:

    python3 generate_giesekus.py

The default grid is 9 x 10 x 5 = 450 runs:
    lambda in [1, 10]
    betav  in [0.1, 0.9]
    alpha  in [0.1, 0.5]
"""

from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path


# =============================================================================
# Easy-to-change settings
# =============================================================================

# Set USE_DOCKER = True only when running this script on the host machine.
# Inside the TFEM container, use the local executable in this directory.
USE_DOCKER = False
BUILD_SOLVER = True
DOCKER_IMAGE = os.environ.get(
    "ROMLAB_DOCKER_IMAGE",
    "fedora-tfem-gfortran-gmsh:latest",
)
DOCKER_WORKDIR = "/shared"

OUTPUT_DIR = "lambda-beta-alpha-curve4"
EXECUTABLE = "./extrudate_swell2d_c_cluster"

LAMBDA_MIN = 1.0
LAMBDA_MAX = 10.0
N_LAMBDA = 6         # 9 gives 450 rows; use 10 for 500 rows.

BETA_MIN = 0.1
BETA_MAX = 0.9
N_BETA = 5

ALPHA_MIN = 0.1
ALPHA_MAX = 0.5
N_ALPHA = 4

MODEL = 3            # Giesekus
PLANAR = False        # False = axisymmetric, True = planar
ETA_0 = 1.0
U_AVG = 0.025
DELTAT = 0.8
NUMTIMESTEPS = 700
DX_BOX = 0.1
DX_WALL = 0.1
DX_INLET = 0.1
VTKEVERY = 0

STOP_WHEN_STEADY = True
STEADY_HEIGHT_TOL = 1.0e-5
STEADY_END_MAX_TOL = 1.0e-2
STEADY_HEIGHT_STEPS = 10
STEADY_MIN_STEPS = 80


# =============================================================================
# Script internals
# =============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / OUTPUT_DIR
LOG_DIR = OUT_DIR / "run_logs"


def linspace(start: float, stop: float, count: int) -> list[float]:
    if count == 1:
        return [(start + stop) / 2.0]
    step = (stop - start) / (count - 1)
    return [start + i * step for i in range(count)]


def parameter_grid() -> list[tuple[float, float, float]]:
    lambda_values = linspace(LAMBDA_MIN, LAMBDA_MAX, N_LAMBDA)
    beta_values = linspace(BETA_MIN, BETA_MAX, N_BETA)
    alpha_values = linspace(ALPHA_MIN, ALPHA_MAX, N_ALPHA)

    return [
        (lam, beta, alpha)
        for alpha in alpha_values
        for beta in beta_values
        for lam in lambda_values
    ]


def solver_input(lam: float, beta: float, alpha: float) -> str:
    planar = ".true." if PLANAR else ".false."
    stop_when_steady = ".true." if STOP_WHEN_STEADY else ".false."
    return (
        "&comppar\n"
        f"  model        = {MODEL}\n"
        f"  lambda       = {lam:.16f}\n"
        f"  mobility     = {alpha:.16f}\n"
        f"  betav        = {beta:.16f}\n"
        f"  eta_0        = {ETA_0:.16f}\n"
        f"  U_avg        = {U_AVG:.16f}\n"
        f"  planar       = {planar}\n"
        f"  deltat       = {DELTAT:.16f}\n"
        f"  numtimesteps = {NUMTIMESTEPS}\n"
        f"  dx_box       = {DX_BOX:.16f}\n"
        f"  dx_wall      = {DX_WALL:.16f}\n"
        f"  dx_inlet     = {DX_INLET:.16f}\n"
        f"  vtkevery     = {VTKEVERY}\n"
        f"  stop_when_steady  = {stop_when_steady}\n"
        f"  steady_height_tol = {STEADY_HEIGHT_TOL:.16e}\n"
        f"  steady_end_max_tol = {STEADY_END_MAX_TOL:.16e}\n"
        f"  steady_height_steps = {STEADY_HEIGHT_STEPS}\n"
        f"  steady_min_steps = {STEADY_MIN_STEPS}\n"
        "/\n"
    )


def docker_command(shell_command: str) -> list[str]:
    return [
        "docker",
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={SCRIPT_DIR},target={DOCKER_WORKDIR}",
        "-w",
        DOCKER_WORKDIR,
        DOCKER_IMAGE,
        "bash",
        "-lc",
        shell_command,
    ]


def run_command(
    command: list[str],
    log_file: Path,
    input_text: str | None = None,
) -> None:
    result = subprocess.run(
        command,
        cwd=SCRIPT_DIR,
        input=input_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    log_file.write_text(
        "COMMAND:\n"
        + " ".join(command)
        + "\n\nSTDOUT:\n"
        + result.stdout
        + "\n\nSTDERR:\n"
        + result.stderr,
        encoding="utf-8",
    )

    if result.returncode != 0:
        tail = (result.stderr or result.stdout)[-2000:]
        raise RuntimeError(f"command failed; see {log_file}\n{tail}")


def prepare_output_files() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    LOG_DIR.mkdir(exist_ok=True)

    for name in [
        "parameters.txt",
        "parameters_completed.txt",
        "curve4_x.txt",
        "curve4_y.txt",
        "curve4_x_mismatch.log",
        "pressure_drop.txt",
    ]:
        path = OUT_DIR / name
        if path.exists():
            path.unlink()


def clear_solver_outputs() -> None:
    for name in [
        "input.txt",
        "curve4_x.txt",
        "curve4_y.txt",
        "curve4_x.out",
        "curve4_y.out",
        "pressure_drop.txt",
        "mesh.geo",
        "mesh.msh",
        "mesh.vtk",
        "mesh_inlet.geo",
        "mesh_inlet.msh",
        "mesh_inlet.vtk",
        "outputmesh.out",
        "outputmesh_inlet.out",
        "cval.out",
        "timings.txt",
    ]:
        path = SCRIPT_DIR / name
        if path.exists():
            path.unlink()

    for vtk_file in SCRIPT_DIR.glob("*.vtk"):
        vtk_file.unlink()


def read_last_row(path: Path) -> str:
    rows = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not rows:
        raise RuntimeError(f"{path} is empty")
    return " ".join(rows[-1].split())


def build_solver() -> None:
    if not BUILD_SOLVER:
        return

    if USE_DOCKER:
        command = docker_command("make extrudate_swell2d_c_cluster")
    else:
        command = ["make", "extrudate_swell2d_c_cluster"]

    run_command(command, LOG_DIR / "build.log")


def run_solver(input_text: str, log_file: Path) -> None:
    input_path = SCRIPT_DIR / "input.txt"
    input_path.write_text(input_text, encoding="utf-8")

    if USE_DOCKER:
        command = docker_command(f"{EXECUTABLE} < input.txt")
        run_command(command, log_file)
    else:
        command = [str(SCRIPT_DIR / "extrudate_swell2d_c_cluster")]
        run_command(command, log_file, input_text=input_text)


def write_requested_parameters(grid: list[tuple[float, float, float]]) -> None:
    with (OUT_DIR / "parameters.txt").open("w", encoding="utf-8") as handle:
        handle.write("# row lambda betav alpha\n")
        for row, (lam, beta, alpha) in enumerate(grid, start=1):
            handle.write(f"{row} {lam:.16f} {beta:.16f} {alpha:.16f}\n")


def append_completed_parameters(
    row: int,
    lam: float,
    beta: float,
    alpha: float,
) -> None:
    with (OUT_DIR / "parameters_completed.txt").open("a", encoding="utf-8") as handle:
        handle.write(f"{row} {lam:.16f} {beta:.16f} {alpha:.16f}\n")


def main() -> None:
    grid = parameter_grid()

    print(f"Output directory: {OUT_DIR}")
    print(f"Grid size: {N_LAMBDA} x {N_BETA} x {N_ALPHA} = {len(grid)} runs")
    print("curve4_y.txt will get one row per run.")
    print("curve4_x.txt will be copied once from the first run.")

    prepare_output_files()
    write_requested_parameters(grid)
    build_solver()

    curve4_x_reference: str | None = None
    start_time = time.time()

    for row, (lam, beta, alpha) in enumerate(grid, start=1):
        print(
            f"[{row:04d}/{len(grid):04d}] "
            f"lambda={lam:.8f} beta={beta:.8f} alpha={alpha:.8f}"
        )

        clear_solver_outputs()
        input_text = solver_input(lam, beta, alpha)
        run_solver(input_text, LOG_DIR / f"run_{row:04d}.log")

        curve4_x_row = read_last_row(SCRIPT_DIR / "curve4_x.txt")
        curve4_y_row = read_last_row(SCRIPT_DIR / "curve4_y.txt")

        if curve4_x_reference is None:
            curve4_x_reference = curve4_x_row
            (OUT_DIR / "curve4_x.txt").write_text(
                curve4_x_row + "\n",
                encoding="utf-8",
            )
        elif curve4_x_row != curve4_x_reference:
            with (OUT_DIR / "curve4_x_mismatch.log").open("a", encoding="utf-8") as handle:
                handle.write(
                    f"row {row}: curve4_x changed for "
                    f"lambda={lam:.16f} betav={beta:.16f} alpha={alpha:.16f}\n"
                )

        with (OUT_DIR / "curve4_y.txt").open("a", encoding="utf-8") as handle:
            handle.write(curve4_y_row + "\n")

        # pressure drop of this run: step time p_in p_out dp_pressure dp_total area
        pressure_drop_path = SCRIPT_DIR / "pressure_drop.txt"
        pressure_drop_row = (
            read_last_row(pressure_drop_path) if pressure_drop_path.exists() else "nan"
        )
        with (OUT_DIR / "pressure_drop.txt").open("a", encoding="utf-8") as handle:
            handle.write(f"{row} {pressure_drop_row}\n")

        append_completed_parameters(row, lam, beta, alpha)

    elapsed = time.time() - start_time
    print(f"Finished {len(grid)} runs in {elapsed:.2f} seconds.")
    print(f"Wrote {OUT_DIR / 'curve4_y.txt'}")
    print(f"Wrote {OUT_DIR / 'curve4_x.txt'}")


if __name__ == "__main__":
    main()
