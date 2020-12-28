
import os
import os.path
import time
import glob
import sys

if len(sys.argv) != 7:
    print("Wrong input parameters\n\n")
    exit()

initial_starts = sys.argv[1]
restraints_dir = sys.argv[2]
output_dir = sys.argv[3]
weight_file = sys.argv[4]
proc_num = int(sys.argv[5])
target_path = sys.argv[6]

if not os.path.exists(output_dir):
    os.system("mkdir -p " + output_dir)

if not os.path.exists(initial_starts):
    print("Failed to find " + initial_starts)
    exit()

if not os.path.exists(restraints_dir):
    print("Failed to find " + restraints_dir)
    exit()

if not os.path.exists(weight_file):
    print("Failed to find " + weight_file)
    exit()

if not os.path.exists(target_path):
    print("Failed to find " + target_path)
    exit()

# reads all the file in a folder

# fileNames = glob.glob(fasta_dir+"/*fasta")

fileNames = []
with open(target_path, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            clos = line.strip().split()
            fileNames.append(clos[0])

shell_dir = output_dir + '/shell_files'

if not os.path.exists(shell_dir):
    os.system("mkdir -p " + shell_dir)

for filepath in fileNames:
    filename = os.path.basename(filepath)
    target_id = filename.split('.')[0]  # T0949.fasta
    initial_start = initial_starts + '/' + target_id + '_modified.pdb'
    restraint_file = restraints_dir + '/' + target_id + '_.rr'

    if not os.path.isfile(initial_start):
        continue

    if not os.path.isfile(restraint_file):
        continue

    print("Generating shell script for " + target_id + "\n")
    work_dir = output_dir + '/' + target_id + '_out'
    if not os.path.exists(work_dir):
        os.system("mkdir -p " + work_dir)

    run_file = shell_dir + '/' + target_id + ".sh"
    os.system("touch " + run_file + ".queued")
    with open(run_file, "w") as f:
        f.write("#!/bin/bash -l\n")
        f.write("#SBATCH -J  " + str(target_id) + "\n")
        f.write("#SBATCH -o " + target_id + "-\%j.out\n")
        f.write("#SBATCH --partition Lewis,hpc4,hpc5\n")
        f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --ntasks=1\n")
        f.write("#SBATCH --cpus-per-task=1\n")
        f.write("#SBATCH --mem-per-cpu=10G\n")
        f.write("#SBATCH --time 2-00:00\n")
        f.write("mv " + run_file + ".queued" + " " + run_file + ".running" + "\n")
        cmd = "sh /storage/hpc/data/esdft/run_dock_mc.sh " + target_id + " " + initial_start + " " + restraint_file + " " + work_dir + " " + weight_file
        f.write("printf \"" + cmd + "\"\n")
        f.write(
            "sh /storage/hpc/data/esdft/run_dock_mc.sh " + target_id + " " + initial_start + " " + restraint_file + " " + work_dir + " " + weight_file + "\n")
        f.write("mv " + run_file + ".running" + " " + run_file + ".done" + "\n")


starttime = time.time()
for filepath in fileNames:
    filename = os.path.basename(filepath)
    target_id = filename.split('.')[0]  # T0949.fasta
    run_file = shell_dir + '/' + target_id + ".sh"
    print("Running " + run_file + "\n")
    ### check the running jobs
    while True:
        running_jobs_num = 1
        done_jobs_num = 0
        total_jobs_num = 0
        files_list = glob.glob(shell_dir + "/*")
        for item in files_list:
            if item[-3:] == '.sh':
                total_jobs_num += 1
            if item[-5:] == '.done':
                done_jobs_num += 1
            if item[-8:] == '.running':
                running_jobs_num += 1

        time.sleep(2)
        if running_jobs_num <= proc_num:
            break
        else:
            continue

    ### submit jobs
    print("sh " + run_file + " &> " + run_file + ".log" + " &\n")
    os.system("sh " + run_file + " &> " + run_file + ".log" + " &")  # run it in background
    print("Running(" + str(running_jobs_num) + ")" + " Done(" + str(done_jobs_num) + ")" + " Total(" + str(
        total_jobs_num) + ") " + "\n")
    time.sleep(1)

print(fileNames)

#### wait until all jobs are done
while True:
    running_jobs_num = 1
    done_jobs_num = 0
    total_jobs_num = 0
    files_list = glob.glob(shell_dir + "/*")
    for item in files_list:
        if item[-8:] == '.running':
            running_jobs_num += 1

    time.sleep(2)
    if running_jobs_num == 0:
        break
    else:
        continue

print("Execution Time: " + str(time.time() - starttime))

#### run in parallel with 10 proteins one time
