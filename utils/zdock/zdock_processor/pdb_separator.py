import copy
import os


def file_reader(input_dir):
    _input_dir=copy.deepcopy(input_dir)
    contents = ""
    f = open(_input_dir, "r")
    if f.mode == 'r':
        contents = f.read()
        f.close()
    return contents.split("END")


def write2File(_filename, _cont):
    with open(_filename, "w") as f:
        f.writelines(_cont)
        if _cont[len(_cont) - 1].strip() != "END":
            f.write("END")
    return


def specific_dir_reader(_input_dir):
    file_names = []
    i = 0
    for root, directories, files in os.walk(_input_dir):
        i = i + 1
        file_names.append(directories)

    return file_names[0]


_input_pdb = "/home/rajroy/q3_het30/pdbs/"
list_of_pdb_dir = specific_dir_reader(_input_pdb)
input_file = "/home/rajroy/complex.1.pdb"
output_file = "/home/rajroy/q3_het30/sep_pdbs/"
for pdb in list_of_pdb_dir:
    pdb_complex_1 = _input_pdb + "/" + pdb + '/' + "complex.1.pdb"
    print(pdb_complex_1)
    names = pdb.split("_")
    model_a = file_reader(pdb_complex_1)
    file_name_a = output_file + "/" + names[0]+".atom"
    file_name_b = output_file + "/" + names[1]+".atom"
    write2File(file_name_a, model_a[0])
    write2File(file_name_b, model_a[1])
