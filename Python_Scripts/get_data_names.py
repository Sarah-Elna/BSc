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
first_four_list=[]
for four in file_names:
    first_four = four[0:4]
    if first_four not in first_four_list:
        first_four_list.append(first_four)

# create the empty names.txt file
n = open('{}\\names.txt'.format(my_path2), 'w')
n.write("")
n.close()

# adds the unique four cifers to the names.txt file
for i in range(len(first_four_list)):
    n = open('{}\\names.txt'.format(my_path2), 'a')
    n.write("{}\n".format(first_four_list[i]))
    n.close()