import os
import sys

__all__ = ['write_pairs']

def write_pairs(site_dir: str, output_dir: str, reference_pdb=''):
    """Writes a file enumerating all pairs of pdb files to check.
    
    If reference pdb provided, only write pairs comparing to that pdb.
    """
    output_file = os.path.join(output_dir, 'PairList.txt')
    lines = 0
    reference_pdb_filename = os.path.basename(reference_pdb) if reference_pdb else ''
    with open(output_file, 'w') as out:
        for i in os.listdir(site_dir):
            for j in os.listdir(site_dir):
                if not reference_pdb_filename or reference_pdb_filename in [i, j]:
                    out.write(i + '\t' + j + '\n')
                    lines += 1
    return output_file


if __name__ == '__main__':
    if len(sys.argv) >= 3:
        site_folder = sys.argv[1]
        output_folder = sys.argv[2]
        reference_pdb = None
        if len(sys.argv) == 4:
            reference_pdb = sys.argv[3]
    else:
        print('Pairs.py <Site-Folder> <Output-Folder> <Reference-PDB (optional)>')
        sys.exit()
    output_file = write_pairs(site_folder, output_folder, reference_pdb=reference_pdb)
