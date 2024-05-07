#!/bin/bash
#
#SBATCH --job-name="compress"
#SBATCH --partition=thin
#SBATCH --time=2-00:00:00
#SBATCH -n 588
#SBATCH -o stdout-batch/slurm-compress.out
#SBATCH -e stdout-batch/slurm-compress.out

tar cvzf /gpfs/scratch1/nodespecific/int4/colomeso/tests/PerforatedCylinder.jl/data/VTKs/serial_results.tgz /gpfs/scratch1/nodespecific/int4/colomeso/tests/PerforatedCylinder.jl/data/VTKs/serial_*

