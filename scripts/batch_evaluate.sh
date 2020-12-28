#!/bin/bash
#SBATCH -J  GFOLD
#SBATCH -o GFOLD-%j.out
#SBATCH --partition Lewis,hpc5,hpc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

#SBATCH --mem-per-cpu=10G
#SBATCH --time 48:00:00


python P1_run_MC_parallel.py /storage/hpc/data/esdft/modified_pdbs /storage/hpc/data/esdft/Res_Con /storage/hpc/data/esdft/results_docking_MC /storage/hpc/data/esdft/talaris2013.wts 10 /storage/hpc/data/esdft/targets.txt

