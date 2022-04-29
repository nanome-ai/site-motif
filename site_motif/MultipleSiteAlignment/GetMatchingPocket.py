"""
Find the best match for a given pocket in provided alignment output file.
"""
import sys
import logging
import operator


logging.basicConfig(level=logging.DEBUG, format="[%(levelname)s] %(message)s")


if len(sys.argv) == 3:
    align_output = sys.argv[1]
    fixed_file_name = sys.argv[2]
else:
    print('GetMatchingPocket.py <align_output.txt> <fixed_file_name>')
    sys.exit()


align_output_data = open(align_output, 'r').readlines()
arr = []
best_match = ''
best_match_score = -1

for line in align_output_data:
    line = line.strip()
    line_split = line.split('\t')
    # We're only concerned with alignments involving the provided file.
    if not operator.xor(
        line_split[0] == fixed_file_name,
            line_split[1] == fixed_file_name):
        continue

    if len(line_split) == 4 and line_split[2] != 'None':
        pdb_name = line_split[0]
        scores = line_split[2].split(' ')
        alignment = line_split[3]
        if len(scores) == 5:
            x1, res_match, mapp_score, x2, x3 = scores
            if float(mapp_score) > best_match_score:
                best_match_score = float(mapp_score)
                best_match = line
    
print(best_match)
