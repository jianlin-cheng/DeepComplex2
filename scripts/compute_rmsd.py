

import os
import os.path
import glob
import sys

from pyrosetta import *
from rosetta import *
from rosetta.protocols.rigid import *
from rosetta.core.scoring import *
init()

if len(sys.argv) != 4:
    print('Wrong input parameters\n\n')
    print(len(sys.argv))
    exit()


native_pdb = sys.argv[1]
predicted_pdb = sys.argv[2]
OUT = sys.argv[3]

file_name = os.path.basename(native_pdb)
target_id = file_name.split('.')[0]

true_pose = pyrosetta.pose_from_pdb(native_pdb)
predicted_pose = pyrosetta.pose_from_pdb(predicted_pdb)


rmsd = CA_rmsd(true_pose, predicted_pose)
print(rmsd)

file_path = OUT + '/' + 'rmsd_MC.txt'

with open(file_path, 'a') as f:
    f.write(target_id)
    f.write(' ')
    f.write(str(rmsd))
    f.write('\n')
