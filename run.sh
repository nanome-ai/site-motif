data_folder=$1
pairs_list="$(pwd)/PairList.txt"
pdb_sizes="$(pwd)/PDBSize.txt"
python $(pwd)/site_motif/MultipleSiteAlignment/Pairs.py $data_folder
python $(pwd)/site_motif/MultipleSiteAlignment/PDBSize.py $data_folder
mpiexec -n 5 python $(pwd)/site_motif/MultipleSiteAlignment/pocket_matrix_mpi7.py $data_folder $pairs_list $pdb_sizes
