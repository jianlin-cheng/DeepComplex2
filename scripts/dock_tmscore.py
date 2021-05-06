import os
import os.path
import glob
import sys

if len(sys.argv) != 3:
    print('Wrong input parameters\n\n')
    print(len(sys.argv))
    exit()


native_dir = sys.argv[1]
predicted_dir = sys.argv[2]

fileNames = glob.glob(native_dir+'/*pdb')


for file_path in fileNames:
    file_name = os.path.basename(file_path)
    target_id = file_name.split('.')[0]

    native_file = native_dir + '/' + target_id + '.pdb'
    
    predicted_file = predicted_dir + '/' + target_id + '_MC.pdb'
    print(native_file)
    print(predicted_file)


    if os.path.isfile(predicted_file):
        path = "/home/esdft/data/results_dock/{}.txt".format(target_id)

        with open(path, "w") as f:

            cmd_n = '/storage/hpc/data/esdft/TMalign' + ' ' + native_file + ' ' + predicted_file + ' ' + '>' + path
            os.system(cmd_n)


