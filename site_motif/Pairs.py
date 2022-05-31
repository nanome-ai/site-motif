import os
import sys

if len(sys.argv) == 3:
    site_folder = sys.argv[1]
    output_folder = sys.argv[2]
else:
    print('Pairs.py <Site-Folder> <Output-Folder>')
    sys.exit()

output_file = os.path.join(output_folder, 'PairList.txt')

with open(output_file, 'w') as out:
    for i in os.listdir(site_folder):
        for j in os.listdir(site_folder):
            out.write(i + '\t' + j + '\n')
