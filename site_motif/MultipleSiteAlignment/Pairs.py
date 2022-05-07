import os
import sys

if len(sys.argv) == 3:
    site_folder = sys.argv[1]
    output_folder = sys.argv[2]
else:
    print('Pairs.py <Site-Folder> <Output-Folder>')
    sys.exit()

dire = os.getcwd()
output_file = os.path.join(output_folder, 'PairList.txt')

with open(output_file, 'w') as out:
    for i in os.listdir(dire+'/'+ site_folder):
        for j in os.listdir(dire+'/'+ site_folder):
            out.write(i + '\t' + j + '\n')
