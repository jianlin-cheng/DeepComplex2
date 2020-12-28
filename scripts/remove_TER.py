
import os
import os.path
import time
import glob
import sys

if len(sys.argv) != 2:
    print("Wrong input parameters\n\n")
    exit()

working_dir = sys.argv[1]

if not os.path.exists(working_dir):
    os.system("mkdir -p " + working_dir)


def remove_ter(pdb_file):
    new_pdb = []
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

        for line in lines:
            if not line.startswith('TER'):
                new_pdb.append(line)

    with open(pdb_file, 'w') as f:
        for new_line in new_pdb:
            f.write(new_line)


fileNames = glob.glob(working_dir+"/*pdb")

for pdb_file in fileNames:
    remove_ter(pdb_file)

