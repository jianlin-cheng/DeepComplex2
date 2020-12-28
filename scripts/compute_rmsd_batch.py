

import os
import os.path
import glob
import sys

if len(sys.argv) != 4:
    print('Wrong input parameters\n\n')
    print(len(sys.argv))
    exit()


native_dir = sys.argv[1]
predicted_dir = sys.argv[2]
OUT = sys.argv[3]


fileNames = glob.glob(native_dir+'/*pdb')

for file_path in fileNames:
    file_name = os.path.basename(file_path)
    target_id = file_name.split('.')[0]

    native_file = native_dir + '/' + target_id + '.pdb'
    print(native_file)
    predicted_file = predicted_dir + '/' + target_id + '_MC.pdb'
    print(predicted_file)

    if os.path.isfile(predicted_file):

        cmd_n = 'python' + ' ' + '/storage/hpc/data/esdft/compute_rmsd.py' + ' ' + native_file + ' ' + predicted_file + ' ' + OUT

        os.system(cmd_n)

