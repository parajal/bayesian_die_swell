#!/bin/bash
#SBATCH --job-name=swell2d
#SBATCH --output=%j.out
#SBATCH --partition=mech-cem.cpu.q
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=20G

set -e

module purge
module load Umbrella/2024
module load intel/2024a
module load gmsh/4.14.0-foss-2024a

mkdir -p input_files
input_file="input_files/input_swell.txt"

cat > "$input_file" <<EOF
&comppar
  model = 2
  lambda = 2.5
  betav = 0.6
  mobility = 0.2
  eta_0 = 1.0
  U_avg = 0.05
  planar = .false.
  deltat = 0.2
  numtimesteps = 700
  dx_box = 0.4
  dx_wall = 0.4
  dx_inlet = 0.4
  vtkevery = 100
  stop_when_steady = .true.
  steady_height_tol = 1.01e-05
  steady_end_max_tol = 1.0e-04
  steady_height_steps = 10
  steady_min_steps = 80
/
EOF

./extrudate_swell_cluster < "$input_file"