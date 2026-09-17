#!/bin/bash
#SBATCH --job-name=swell2d
#SBATCH --output=%j.out
#SBATCH --partition=mech-cem.cpu.q
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --mem=20G

module load intel

set -euo pipefail

mkdir -p input_files sweep_results

lambda_values=(1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0)
betav_values=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
alpha_values=(0.1 0.2 0.3 0.4 0.45)

aggregate_tmp="sweep_results/curve4_y_sweep.tmp"
params_file="sweep_results/curve4_params.out"

: > "$aggregate_tmp"
printf "%s\n" "# row lambda betav alpha" > "$params_file"

row=0

for lambda_value in "${lambda_values[@]}"; do
  for betav_value in "${betav_values[@]}"; do
    for alpha_value in "${alpha_values[@]}"; do
      row=$((row + 1))
      run_id=$(printf "row%03d_lambda_%s_betav_%s_alpha_%s" \
        "$row" "$lambda_value" "$betav_value" "$alpha_value")
      run_id=${run_id//./p}
      input_file="input_files/input_${run_id}.txt"

      cat > "$input_file" <<EOF
&comppar
  model        = 3
  lambda       = $lambda_value
  mobility     = $alpha_value  ! Giesekus alpha; PTT epsilon for model=5/6
  betav        = $betav_value
  eta_0        = 1.0
  U_avg        = 0.1
  planar       = .false.
  deltat       = 0.4
  numtimesteps = 600
  dx_box       = 0.1
  dx_wall      = 0.1
  dx_inlet     = 0.1
  vtkevery     = 100
/
EOF

      echo "Starting simulation $row: lambda=$lambda_value betav=$betav_value alpha=$alpha_value"

      ./extrudate_swell2d_c_cluster < "$input_file"

      cp curve4_y.out "sweep_results/curve4_y_${run_id}.out"
      cp curve4_x.out "sweep_results/curve4_x_${run_id}.out"
      cat curve4_y.out >> "$aggregate_tmp"
      printf "%d %s %s %s\n" "$row" "$lambda_value" "$betav_value" "$alpha_value" >> "$params_file"

      echo "Finished simulation $row"
    done
  done
done

cp "$aggregate_tmp" curve4_y.out

echo "Sweep finished: wrote curve4_y.out and $params_file"
