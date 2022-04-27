data_folder=$1
pairs_list="$(pwd)/PairList.txt"
pdb_sizes="$(pwd)/PDBSize.txt"
echo "Running on $1"
python $(pwd)/MultipleSiteAlignment/Pairs.py $data_folder
python $(pwd)/MultipleSiteAlignment/PDBSize.py $data_folder
mpiexec -n 4 python $(pwd)/MultipleSiteAlignment/pocket_matrix_mpi7.py $data_folder $pairs_list $pdb_sizes
