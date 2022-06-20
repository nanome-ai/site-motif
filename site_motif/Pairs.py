import os
import sys

__all__ = ['write_pairs']

def write_pairs(site_dir, output_dir):
    """Writes a file enumerating all pairs of pdb files to check."""
    output_file = os.path.join(output_dir, 'PairList.txt')
    with open(output_file, 'w') as out:
        for i in os.listdir(site_dir):
            for j in os.listdir(site_dir):
                out.write(i + '\t' + j + '\n')
    return output_file


if __name__ == '__main__':
    if len(sys.argv) == 3:
        site_folder = sys.argv[1]
        output_folder = sys.argv[2]
    else:
        print('Pairs.py <Site-Folder> <Output-Folder>')
        sys.exit()
    output_file = write_pairs(site_folder, output_folder)
    print(output_file)
