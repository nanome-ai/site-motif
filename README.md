# SiteMotif is a graph based method for aligning protein binding sites in sequence-order independent fashion

This repo is forked from https://github.com/santhoshgits/MAPP-3D. We converted it to run on python 3, removed the PairwiseAlignment portion, and modified the directory structure to suit our needs

## Usage

```markdown
./bin/site-motif <sites_dir> <output_dir>

 * sites_dir = folder containing pdb files you wish to align
 * output_dir = folder to store site-motif results
```

---

<span style="color:blue">INSTALLATION INSTRUCTION</span><br>
<p>An MPI library such as MPICH is required</p>
<p>Please refer to https://www.mpich.org/downloads/ for install instructions.</p>
python requirements

<p> After MPICH installed on your host, install python requirements</p>

```markdown
pip install -r requirements.txt 
```

---
<h3><span style="color:red">Handling ERRORS</span></h3>

1) Both version of SiteMotif requires binding site coordinates in .pdb format and chain identifier have to be present for all the residues.

2) While running MPI version of our code. If a progam gets terminated in between due to some unknown reason. Please don't delete 'align_output.txt', pocket_matrix7_mpi.py read  file 'align_output.txt' as a checkpoint and runs only those pairs that are not compared before.

File align_output.txt is tab delimited

---

# Steps for Running MPI version of MAPP

For the purpose of this tutorial, we have added a total of 100 ATP binding sites to a folder 'ATP'. The goal is to run SiteMotif(MPI) for all ATP pairs and find if any motif could be inferred.

1. **Run Pairs.py script -** This will read all sites in the given folder and create a tab separated paired entries stored in 'Pairs.txt'.  
User can provide their custom pair wise entries to avoid pairs that you dont want to compare.
Usage: python Pairs.py ATP <br> 
Output: PairList.txt
PairList.txt is tab separated


2. **Run pdb_res.py -** This is sort all binding sites based on the number of residue present in them.<br> Usage: python PDBSize.py ATP <br>
Output: PDBSize.txt

3. **Running SiteMotif**<br>
USAGE: mpirun -n 4 python pocket_matrix_mpi7.py arg1 arg2 arg3<br>
arg1 - ATP site folder<br>
arg2 - output of Pairs.py<br>
arg3 - output of PDBSize.py

Example) mpirun -n 4 python pocket_matrix_mpi7.py ATP PairList.txt PDBSize.txt <br>
align_output.txt is the output file generated after running the step 3.

4. **Analysing SiteMotif Result**-
File 'align_output.txt' is a tab separated data containing MAPP scores and residue-residue correspondances for all combinations of site pairs as specified in file 'PairList.txt' for which the coordinate is present in the folder 'ATP'


><span style="color:red"> To find representative for the ATP binding site, any clustering algorithm can be used. Here we pick one representative based on number connections with the other ATP sites
after imposing a cutoff of M-dist-min > 0.6 and M-dist-max > 0.4.</span>

NOTE: The nature of binding site varies from one ligand to the next, as well as from one binding sites to another. As a result, it is recommended to do an initial analysis of the site's network and select an acceptable representative.

Run the below provided script to find the representative site<br>
python Analyse.py align_output.txt .4 7

4. **Generating Motif based on the chosen representative (For undermining purposes only)**<br>
Usage: Motif.py <align_output.txt> < representative-site> <No. of residue match><br>
Example: python Motif.py align_output.txt 1B0U_ATP_A_301.pdb 4

Output: [WW]-x-x-[IHF]-x-[VA]-x(18)-[PSA]-[TS]-G-[SA]-G-K-[ST]-T-x(22)-[EQ]-x(78)<br>
The sequence [TS]-G-[SA]-G-K-[ST]-T is a well characterized Walker motif associated with nucleotide binding was identified correctly by SiteMotif. 

## Citations
Sankar S, Chandra N (2022) SiteMotif: A graph-based algorithm for deriving structural motifs in Protein Ligand binding sites. PLoS Comput Biol 18(2): e1009901. https://doi.org/10.1371/journal.pcbi.1009901
