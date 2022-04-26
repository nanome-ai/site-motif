data_folder=$1
pairs_list="PairList.txt"
pdb_sizes="PDBSize.txt"
echo "Running on $data_folder"
python $(pwd)/MultipleSiteAlignment/Pairs.py $data_folder
python $(pwd)/MultipleSiteAlignment/PDBSize.py $data_folder
mpiexec -n 2 python $(pwd)/MultipleSiteAlignment/pocket_matrix_mpi7.py $data_folder $pairs_list $pdb_sizes