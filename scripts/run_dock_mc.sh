#!/bin/bash
#SBATCH -J  GFOLD
#SBATCH -o GFOLD-%j.out
#SBATCH --partition Lewis,hpc5,hpc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

#SBATCH --mem-per-cpu=10G
#SBATCH --time 04:00:00

module load miniconda3
source activate /storage/hpc/data/esdft/envs/pyrosetta_env

targetid=$1 
initial_start=$2 
restraints=$3 
outputfolder=$4 
weight_file=$5



curr_dir="/storage/hpc/data/esdft/results_docking_MC/${targetid}/"

echo "$curr_dir"

mkdir -p "$curr_dir"
cd $curr_dir


python /storage/hpc/data/esdft/docking_mc_parallel.py $initial_start  $restraints  $outputfolder  $weight_file
