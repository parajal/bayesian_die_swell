#!/bin/bash
#SBATCH --job-name=swell_sweep
#SBATCH --output=%j.out
#SBATCH --partition=mech-cem.cpu.q
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

set -euo pipefail

module purge
module load Umbrella/2024
module load foss/2025a
module load gmsh/4.15.0-foss-2025a

which gmsh
gmsh --version

export LD_LIBRARY_PATH="$EBROOTFLEXIBLAS/lib:$LD_LIBRARY_PATH"

ROOTDIR="$PWD"

mkdir -p sweep_results

lambda_values=(1.0 3.0 5.0 7.0 9.0)
betav_values=(0.1 0.3 0.5 0.7 0.9)
alpha_values=(0.0)
eta0_values=(0.1 0.5 1.0 2.0 3.0)

row=0

for lambda_value in "${lambda_values[@]}"; do
  for betav_value in "${betav_values[@]}"; do
    for alpha_value in "${alpha_values[@]}"; do
      for eta0_value in "${eta0_values[@]}"; do

        row=$((row+1))

        run_id=$(printf "case%03d_lambda_%s_betav_%s_alpha_%s_eta0_%s" \
          "$row" "$lambda_value" "$betav_value" "$alpha_value" "$eta0_value")

        run_id=${run_id//./p}

        run_dir="$ROOTDIR/sweep_results/$run_id"

        # Skip already completed cases
        if [ -f "$run_dir/curve4_y.txt" ]; then
          echo "Skipping completed case $row"
          continue
        fi

        mkdir -p "$run_dir"

        cp "$ROOTDIR/mesh_2D_bc.igo" "$run_dir/"
        cp "$ROOTDIR/mesh_inlet_2D_bc.igo" "$run_dir/"

        cat > "$run_dir/input.txt" << EOF
&comppar
model        = 2
lambda       = $lambda_value
mobility     = $alpha_value
betav        = $betav_value
eta_0        = $eta0_value
U_avg        = 0.1
planar       = .false.
deltat       = 0.4
numtimesteps = 1000
dx_box       = 0.1
dx_wall      = 0.1
dx_inlet     = 0.1
vtkevery     = 100
stop_when_steady = .true.
steady_height_tol = 1.0e-05
steady_end_max_tol = 2.0e-04
steady_height_steps = 1
steady_min_steps = 80
steady_round_decimals = 3
/
EOF

        echo "================================="
        echo "CASE $row"
        echo "lambda = $lambda_value"
        echo "betav  = $betav_value"
        echo "alpha  = $alpha_value"
        echo "eta_0  = $eta0_value"
        echo "dir    = $run_dir"
        echo "================================="

        (
          cd "$run_dir"

          "$ROOTDIR/extrudate_swell" < input.txt \
            > simulation.out 2>&1
        )

        if [ ! -f "$run_dir/curve4_y.txt" ]; then
          echo "ERROR: case failed: $run_id"
          tail -50 "$run_dir/simulation.out"
          exit 1
        fi

        echo "Finished case $row"

      done
    done
  done
done

echo "Sweep finished"
