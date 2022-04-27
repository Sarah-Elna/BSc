## This code has the purpouse of extracting the first four unique cifers of my data file names, add them to a list, and print them to a textfile names names.txt

import os
from os import listdir
from os.path import isfile, join

# set path to directory
my_path1 = os.path.abspath("C:/Users/Sarah/Documents/AU/6/Bachelor/GenomeDK/Filer_fra_wolf/trimmed/")

# set path to names.txt
my_path2 = os.path.abspath("C:/Users/Sarah/Documents/AU/6/Bachelor/GitHub/BSc")

# get the file names in a list
file_names = [f for f in listdir(my_path1) if isfile(join(my_path1, f))]

# get the first four unique cifers of the files
name_list=[]
for file_name in file_names:
    if file_name.endswith("12-single.fastq"):
        name = file_name[0:4]
        if name not in name_list:
            name_list.append(name)

# create the empty names.txt file
n = open('{}\\numbers.txt'.format(my_path2), 'w')
n.write("")
n.close()

# adds the unique four cifers to the names.txt file
for i in range(len(name_list)):
    n = open('{}\\numbers.txt'.format(my_path2), 'a')
    n.write("{}\n".format(name_list[i]))
    n.close()