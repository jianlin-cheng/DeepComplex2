import os


def loadFastaDictionary(dict_file):
    fasta_dict = {}
    with open(dict_file, "r") as f:
        for line in f:
            fasta_dict[line.strip().split(":")[0].strip()] = line.strip().split(":")[1].strip()
    return fasta_dict


def binding_only_file_reader(_seq_file):
    binding_fasta = loadFastaDictionary("/exports/store1/raj/future_task_new_db/binding_only_modified_fasta_dict.txt")
    with open(_seq_file) as infile:
        number_counter = 0
        data_counter = 0
        number_counter = 0
        data_counter_2220 = 0
        data_counter_1110 = 0
        file_object = open('/exports/store1/raj/future_task_new_db/new_binding_file.txt', 'a')
        for line in infile:
            # print(line)
            # getting the binding file

            temp_array = line.split("\t")
            fasta_name_a = temp_array[0]
            fasta_name_b = temp_array[1]
            binding_fasta_a = binding_fasta.get(fasta_name_a)
            binding_fasta_b = binding_fasta.get(fasta_name_b)
            total_len=len(binding_fasta_a)+len(binding_fasta_b)
            # if total_len <=2220:

            # Append 'hello' at the end of file


            if total_len <=1110:
                file_object.write(line)
                data_counter_1110=data_counter_1110+1
                # os.system('echo '+line +' >> '+'/exports/store1/raj/future_task_new_db/new_binding_file.txt')
            # name_order_a=fasta_name_a+"_"+fasta_name_b
            # fasta_order_a = binding_fasta_a+binding_fasta_b
            # name_order_b=fasta_name_b+"_"+fasta_name_a
            # fasta_order_b=binding_fasta_b+binding_fasta_a
            # print(name_order_a)
            # print(total_len)
            data_counter_2220 = data_counter_2220 + 1
            if number_counter%1000000==0:
                print("data_counter_1110 " + str(data_counter_1110) + "\n")
                # print("data_counter_1440 "+str(data_counter_1440)+"\n")
                # print("data_counter_2220 " + str(data_counter_2220) + "\n")
                # print("number_counter " + str(number_counter) + "\n")
            number_counter = number_counter + 1

# Close the file
    file_object.close()

binding_only_file_reader('/exports/store1/raj/future_task_new_db/binding_file.txt')