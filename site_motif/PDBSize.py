import os
import sys

__all__ = ['write_pdb_size']


def write_pdb_size(site_dir, output_dir):
	"""Writes a file containing the number of alpha carbons in each pdb file."""
	output_path = os.path.join(output_dir, "PDBSize.txt")
	with open(output_path,'w') as out:
		for pdb_file in os.listdir(site_dir):
			pdb_path = os.path.join(site_dir, pdb_file)
			with open(pdb_path, 'r') as f:
				pdb_lines = f.readlines()
			alpha_carbon_count = 0
			for line in pdb_lines:
				line_split = line.strip().split()
				if line_split[0] == "ATOM":
					if line_split[2] == "CA":
						alpha_carbon_count += 1
			out.write(pdb_file + "\t" + str(alpha_carbon_count) + "\n")

if __name__ == '__main__':
	if len(sys.argv) == 3:
		site_dir = sys.argv[1]
		output_dir = sys.argv[2]
	else:
		print("PDBSize.py <sites_dir> <output_dir>")
		sys.exit()
	write_pdb_size(site_dir, output_dir)
