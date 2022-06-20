import copy
import logging
import math
import numpy as np
import re
import sys
import time
from collections import defaultdict, Counter
from mpi4py import MPI

logging.basicConfig(level=logging.DEBUG, format="[%(levelname)s] %(message)s")

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size() - 1
status = MPI.Status()
logging.debug(f"rank={rank}")
logging.debug(f"size={size}")
'''
 mpirun -mca btl ^openib -n 2 python pocket_matrix_mpi7.py sam_atp_ip sam_atp_ip pdb_res_list 
'''


def PairWise(res, coord):
    logging.debug("running PairWise")
    arr = []
    for i in range(len(res)):
        x_ca, y_ca, z_ca = coord[i][0]
        x_cb, y_cb, z_cb = coord[i][1]
        x_cn, y_cn, z_cn = coord[i][2]

        for j in range(len(res)):
            if i != j:
                x1_ca, y1_ca, z1_ca = coord[j][0]
                x1_cb, y1_cb, z1_cb = coord[j][1]
                x1_cn, y1_cn, z1_cn = coord[j][2]
                ans_ca = math.sqrt(
                    pow((x_ca - x1_ca), 2) + pow((y_ca - y1_ca), 2) + pow((z_ca - z1_ca), 2))
                ans_cb = math.sqrt(
                    pow((x_cb - x1_cb), 2) + pow((y_cb - y1_cb), 2) + pow((z_cb - z1_cb), 2))
                ans_cn = math.sqrt(
                    pow((x_cn - x1_cn), 2) + pow((y_cn - y1_cn), 2) + pow((z_cn - z1_cn), 2))
                arr.append([res[i], res[j], ans_ca, ans_cb, ans_cn])
    return arr


def center_residues(arr, coord):
    coord -= np.mean(coord, axis=0)
    dic = {1: "   ", 2: "  ", 3: " ", 4: ""}
    brr = []
    for i in range(len(arr)):
        j = ["%.3f" % j for j in coord[i]]
        val = ''.join([dic[len(k.split(".")[0])] + k for k in j])
        brr.append(arr[i][:30] + val + arr[i][54:])
    return brr


def file_process(arr):
    '''
    remove redundancy and h_bond
    store coordinate
    1) From whole coordinate -> to centroid
    2) CB -> Again Simple. If Glycine, then Ca
    3) CA
    '''
    logging.debug("running file_process")
    whole_dic = {}
    h_dic = {"H": 0}
    brr, het_arr, coord = [], [], []

    for i in arr:
        i = i.strip()

        if i[:4] == "ATOM":
            if i[74:78].strip() not in h_dic:
                var = i[13:16].strip() + " " + i[17:26].strip()
                #var = i[13:16].strip()+" "+i[17:26].strip()
                # time.sleep(.1)
                if var not in whole_dic:
                    whole_dic[var] = 0
                    brr.append(i)
                    coord.append(
                        [i[28:38].strip(), i[38:46].strip(), i[46:54].strip()])

        if i[:6] == "HETATM":
            brr.append(i)
            coord.append(
                [i[28:38].strip(), i[38:46].strip(), i[46:54].strip()])

    coord = np.asarray(coord, dtype='float')
    brr = center_residues(brr, coord)

    dic1, dic2 = defaultdict(list), defaultdict(list)

    for line in brr:
        if line[:4] == "ATOM":
            val = line[17:20] + '-' + line[21:22] + '-' + line[22:26].strip()
            dic1[val].append(
                [line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
            dic2[val].append(line[13:16].strip())

    res, coord1 = [], []
    for i, j in dic2.items():
        coord = np.asarray(dic1[i], dtype='float')
        cn = np.mean(coord, axis=0)
        if i[:3] == "GLY":
            for j1 in range(len(j)):
                if j[j1] == 'CA':
                    cb = coord[j1]
                if j[j1] == 'CA':
                    ca = coord[j1]

        else:
            for j1 in range(len(j)):
                if j[j1] == 'CB':
                    cb = coord[j1]
                if j[j1] == 'CA':
                    ca = coord[j1]

        res.append(i)
        if len(coord1) == 0:
            coord1 = [[ca, cb, cn]]
        else:
            coord1 = np.append(coord1, [[ca, cb, cn]], axis=0)

    arr = PairWise(res, coord1)
    arr = sorted(arr, key=lambda x: float(x[3]))
    dic, dic1_pairs = defaultdict(list), defaultdict(list)
    arr1 = []

    for i in arr:
        dic[i[0] + ' ' + i[1]].append(i[2])
        dic[i[0] + ' ' + i[1]].append(i[3])
        dic[i[0] + ' ' + i[1]].append(i[4])
        arr1.append(i[0] + ' ' + i[1])
        dic1_pairs[i[0]].append(i[1])
    return arr1, dic, dic1_pairs, brr, het_arr


def compare_matrices(mat1, mat2):

    if (mat2[0] - 1.1) < mat1[0] < (mat2[0] + 1.1):
        if (mat2[1] - 1.1) < mat1[1] < (mat2[1] + 1.1):
            if (mat2[2] - 1.51) < mat1[2] < (mat2[2] + 1.51):
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def CheckDistance(S1, S2, v1, v2):
    # S1, S2 = S1[1:], S2[1:]
    if compare_matrices(res_dic1[v1], res_dic2[v2]):
        v1a = v1.split(' ')[1]
        v2a = v2.split(' ')[1]

        for i in zip(S1, S2):
            p1, p2 = i
            p1 = p1.split(' ')[1]
            p2 = p2.split(' ')[1]

            if p1 == v1a or p2 == v2a:
                return False
            if not compare_matrices(res_dic1[p1 + ' ' + v1a], res_dic2[p2 + ' ' + v2a]):
                return False

        return True
    else:
        return False


def Recursion(res_pair1, sequence):

    if sequence == 'start':
        for i in res_arr2:
            # if res_dic1[res_pair1] == res_dic2[i]:

            if compare_matrices(res_dic1[res_pair1], res_dic2[i]):
                return res_pair1, i, 'First'
        return None, None, 'First'

    else:
        res1a, res1b = res_pair1.split(' ')
        res2a, res2b = sequence.split(' ')
        for i in res_pairs_dic1[res1b]:
            if i != res1a and i not in dic_single1:
                for j in res_pairs_dic2[res2b]:
                    if j != res2a and j not in dic_single2:
                        Check = CheckDistance(
                            SequenceArrays1, SequenceArrays2, res1b + ' ' + i, res2b + ' ' + j)
                        #Check = CheckDistance(SequenceArrays1, SequenceArrays2, [[res1b],[i]], [[res2b],[j]])
                        if Check:
                            if i + '\t' + j not in dic_pair_captch:
                                seq1 = ' '.join(
                                    SequenceArrays1) + ' ' + res1b + ' ' + i
                                seq2 = ' '.join(
                                    SequenceArrays2) + ' ' + res2b + ' ' + j
                                if seq1 + '\t' + seq2 not in dic_whole:
                                    return res1b + ' ' + i, res2b + ' ' + j, 'Next'
        return None, None, 'Next'


def SortedArr():
    return True
    arr = [i[0].split(' ')[1] + ' ' + i[1].split(' ')[1]
           for i in zip(SequenceArrays1, SequenceArrays2)]
    arr = sorted(arr)
    if len(arr) < 100:  # return False for sensitive cases
        return True
    arr = '_'.join(arr)

    # time.sleep(11)
    if arr not in SortedArrDic:
        SortedArrDic[arr] = 0
        return True
    return False


def PairNext(S1, S2):
    dic1, dic2 = {}, {}
    if len(S1) <= 1:
        return None, None, None, None
    for i in zip(S1[:-1], S2[:-1]):

        dic1[i[0].split(' ')[0]] = 0
        dic1[i[0].split(' ')[1]] = 0
        dic2[i[1].split(' ')[0]] = 0
        dic2[i[1].split(' ')[1]] = 0

    S1, S2 = S1[:-1], S2[:-1]
    p1 = S1[-1].split(' ')[1]
    p2 = S2[-1].split(' ')[1]
    for i in res_pairs_dic1[p1]:
        # if i+' '+p1 not in dic1 and p1+' '+i not in dic1:
        if i not in dic1:
            for j in res_pairs_dic2[p2]:
                if j not in dic2:
                    # if j+' '+p2 not in dic2 and p2+' '+j not in dic2:  # res1b+' '+i, res2b+' '+j
                    Check = CheckDistance(S1, S2, p1 + ' ' + i, p2 + ' ' + j)
                    #Check = CheckDistance(SequenceArrays1, SequenceArrays2, p1+' '+i, p2+' '+j)
                    if Check:
                        seq1 = ' '.join(S1) + ' ' + p1 + ' ' + i
                        seq2 = ' '.join(S2) + ' ' + p2 + ' ' + j

                        if seq1 + '\t' + seq2 not in dic_whole:
                            return(S1, S2, p1 + ' ' + i, p2 + ' ' + j)

    return PairNext(S1, S2)


def run():
    global dic_single1, dic_single2, SequenceArrays1, SequenceArrays2, dic_pair_captch, dic_whole
    Final1, Final2 = [], []
    BreakLoop = False

    dic_whole = {}
    for i in res_arr1:
        if BreakLoop:
            break

        ans = 'start'
        SequenceArrays1, SequenceArrays2 = [], []
        dic_single1, dic_single2 = {}, {}

        # add two nested while loops
        InitiateFirstBreak = False
        ObtainedCount = 0

        while True:
            dic_pair_captch = {}
            if ans != 'start':
                if not SequenceArrays1:
                    break
            if InitiateFirstBreak:
                if not SequenceArrays1:
                    break
                if ObtainedCount <= 1:
                    break

                if len(SequenceArrays1) <= 2:
                    seq1 = ' '.join(SequenceArrays1)
                    seq2 = ' '.join(SequenceArrays2)
                    dic_whole[seq1 + '\t' + seq2] = 0
                    break
                ObtainedCount = 0
                SequenceArrays1, SequenceArrays2, i1, ans = PairNext(
                    copy.deepcopy(SequenceArrays1), copy.deepcopy(SequenceArrays2))
                if not i1:
                    break
                SequenceArrays1.append(i1)
                SequenceArrays2.append(ans)

                for j in zip(SequenceArrays1, SequenceArrays2):
                    dic_pair_captch[j[0].split(
                        ' ')[0] + '\t' + j[1].split(' ')[0]] = 0

                i = i1
            while True:
                if not ans:
                    break
                i, ans, CheckPoint = Recursion(i, ans)
                ObtainedCount += 1
                if not ans:
                    break
                if i + '\t' + ans in dic_whole:
                    break

                InitiateFirstBreak = True
                dic_single1[i.split(' ')[0]] = 0
                dic_single2[ans.split(' ')[0]] = 0
                # dic_pair_captch[i.split(' ')[1]+'\t'+ans.split(' ')[1]] = 0
                dic_pair_captch[i.split(' ')[0] + '\t' + ans.split(' ')[0]] = 0
                SequenceArrays1.append(i)
                SequenceArrays2.append(ans)

            if not SortedArr():
                BreakLoop = True

            seq1 = ' '.join(SequenceArrays1)
            seq2 = ' '.join(SequenceArrays2)

            # time.sleep(.1)
            dic_whole[seq1 + '\t' + seq2] = 0

            Final1.append(SequenceArrays1)
            Final2.append(SequenceArrays2)
    logging.debug("Run completed")
    return Final1, Final2


def process_hits(Final1, Final2):
    logging.debug("Starting process_hits()")
    arr = []
    for i in zip(Final1, Final2):
        arr.append([i[0], i[1], len(i[0])])

    arr = sorted(arr, key=lambda x: int(x[2]), reverse=True)
    NewArr = []
    for i in arr:
        if int(i[2]) > 3:
            val1 = [j.split(' ')[1] for j in i[0][:-1]]
            val2 = [j.split(' ')[1] for j in i[1][:-1]]
            NewArr.append([val1, val2, len(val1)])
    '''
	if len(NewArr) < 10:
		for i in arr:
			if int(i[2]) == 3:
				val1 = [ j.split(' ')[1] for j in i[0][:-1] ]
				val2 = [ j.split(' ')[1] for j in i[1][:-1] ]
				if res_dic1[' '.join(val1)][2] < 4:
					
					NewArr.append([val1, val2, len(val1)])	
	'''
    NewArray = []
    if not NewArr:
        return None
    NewArray.append(NewArr[0])  # base condition
    for i in NewArr:
        check = True
        for j in NewArray:
            dic = {j[0][k] + ' ' + j[1][k]: 0 for k in range(len(j[0]))}

            dic_count = sum([1 for k in range(len(i[0]))
                             if i[0][k] + ' ' + i[1][k] in dic])

            if dic_count == len(j[0]):
                check = False
                break
        if check:
            NewArray.append(i)
    logging.debug("finished process_hits()")
    return NewArray

# END OF MAIN MAPP CODE
# START OF ALIGNMENT CODE


def rmsd(V, W):
    D = len(V[0])
    N = len(V)
    result = 0.0
    for v, w in zip(V, W):
        result += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(result / N)


def kabsch_rmsd(P, Q, translate=False):
    P = kabsch_rotate(P, Q)
    return rmsd(P, Q)


def kabsch_rotate(P, Q):
    U = kabsch(P, Q)
    P = np.dot(P, U)
    return P


def kabsch(P, Q, check):
    if check:
        P -= np.mean(P, axis=0)
        Q -= np.mean(Q, axis=0)
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return U


def centroid(X):
    centroid = X.mean(axis=0)
    return centroid


def blosum():
    global residue_pairs_dictionary
    residue_dict_single = {'G': 'GLY', 'A': 'ALA', 'V': 'VAL', 'L': 'LEU', 'I': 'ILE', 'T': 'THR', 'S': 'SER', 'Y': 'TYR', 'W': 'TRP',
                           'P': 'PRO', 'F': 'PHE', 'N': 'ASN', 'Q': 'GLN', 'D': 'ASP', 'E': 'GLU', 'R': 'ARG', 'K': 'LYS', 'H': 'HIS',
                           'C': 'CYS', 'M': 'MET'
                           }
    aline = []
    aline.append(
        "A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *")
    aline.append(
        "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4")
    aline.append(
        "R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4")
    aline.append(
        "N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4")
    aline.append(
        "D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4")
    aline.append(
        "C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4")
    aline.append(
        "Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4")
    aline.append(
        "E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4")
    aline.append(
        "G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4")
    aline.append(
        "H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4")
    aline.append(
        "I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4")
    aline.append(
        "L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4")
    aline.append(
        "K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4")
    aline.append(
        "M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4")
    aline.append(
        "F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4")
    aline.append(
        "P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4")
    aline.append(
        "S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4")
    aline.append(
        "T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4")
    aline.append(
        "W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4")
    aline.append(
        "Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4")
    aline.append(
        "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4")
    res_info = aline[0].strip()
    res_info = re.sub(" {1,}", " ", res_info)
    res_info_split = res_info.split(" ")[:-4]
    ans = []
    for line in aline[1:]:
        line = line.strip()
        line = re.sub(" {1,}", " ", line)
        residue = line.split(" ")[0]
        min_max = []
        for i in line.split(" ")[1:-4]:
            min_max.append(int(i))
        minimum = min(min_max)
        if minimum < 0:
            new_min = abs(minimum)
            new_min_max = []
            for i in min_max:
                new_min_max.append(i + new_min)

        for j, k in zip(res_info_split, new_min_max):
            ans.append(residue_dict_single[residue] +
                       " " + residue_dict_single[j] + " " + str(k))

    dic_temp = {}
    residue_pairs_dictionary = {}
    for i in ans:
        i0 = i.split(" ")[0]
        i1 = i.split(" ")[1]
        val = int(i.split(" ")[2])
        if i0 not in dic_temp:
            dic_temp[i0] = 0
            residue_pairs_dictionary[i0] = {}
        residue_pairs_dictionary[i0][i1] = val


def dihedral1(aa1, aa2):
    arr = ["_CA", "_CN", "_N"]
    arr3 = []
    # [('_CA', '_CN', '_N'), ('_CA', '_N', '_CN'), ('_CN', '_CA', '_N'), ('_CN', '_N', '_CA'), ('_N', '_CA', '_CN'), ('_N', '_CN', '_CA')]
    arr3.append(["_CA", "_CN", "_N"])
    arr1, arr2 = [], []
    ans = []
    for j in arr3:
        arr1, arr2 = [], []
        for i, k in zip(arr, j):
            if aa1 + i in arr1_dihedral_dic and aa2 + k in arr2_dihedral_dic:
                arr1.append(arr1_dihedral_dic[aa1 + i])
                arr2.append(arr2_dihedral_dic[aa2 + k])
        if len(arr1) == 3:
            try:
                arr1 = np.asarray(arr1, dtype=float)
                arr2 = np.asarray(arr2, dtype=float)
                ap1, ap2, ap3 = arr1[0], arr1[1], arr1[2]
                bp1, bp2, bp3 = arr2[0], arr2[1], arr2[2]
                av1 = ap3 - ap1
                av2 = ap2 - ap1
                acp = np.cross(av1, av2)
                bv1 = bp3 - bp1
                bv2 = bp2 - bp1
                bcp = np.cross(bv1, bv2)
                a1, b1, c1, a2, b2, c2 = acp[0], acp[1], acp[2], bcp[0], bcp[1], bcp[2]
                d = (a1 * a2 + b1 * b2 + c1 * c2)
                e1 = math.sqrt(a1 * a1 + b1 * b1 + c1 * c1)
                e2 = math.sqrt(a2 * a2 + b2 * b2 + c2 * c2)
                d = d / (e1 * e2)
                A = math.degrees(math.acos(d))
                ans.append(A)
            except:
                ans.append(10.0)
        else:
            ans.append(10.0)

    for i in range(0, 360, 10):
        if i <= ans[0] <= i + 10:
            return 360.0 - i
    return 360.0


def SiteGen():
    arr = copy.deepcopy(B_all)
    minim = []

    ans_dic = defaultdict(list)
    for i in range(0, len(arr)):
        info1, name1 = pdb1_res_info[i]
        if name1 == "CA":
            # p#rint name1
            x, y, z = arr[i]
            for j in range(0, len(site2_coord)):
                info2, name2 = pdb2_res_info[j]
                if name2 == "CA":
                    x1, y1, z1 = site2_coord[j]
                    x_ans = pow((x - x1), 2)
                    y_ans = pow((y - y1), 2)
                    z_ans = pow((z - z1), 2)
                    ans = math.sqrt(x_ans + y_ans + z_ans)
                    minim.append(ans)
                    if ans < 1.25:

                        ans_dic[info1].append(info2)

    global site_check_dic_len, site_check_res_corr, site_check_sum1, arr1_dihedral_dic, arr2_dihedral_dic
    arr1_dihedral_dic, arr2_dihedral_dic = {}, {}

    dic_cent, dic_ca, dic_n = defaultdict(
        list), defaultdict(list), defaultdict(list)
    for i in range(0, len(arr)):
        info, name = pdb1_res_info[i]
        dic_cent[info].append(arr[i])
        if name == "CA":
            dic_ca[info].append(arr[i])
        if name == "N":

            dic_n[info].append(arr[i])

    for i in dic_cent.items():

        arr1_dihedral_dic[i[0] + "_CN"] = np.mean(i[1], axis=0)
        arr1_dihedral_dic[i[0] + "_CA"] = np.asarray(dic_ca[i[0]])[0]
        arr1_dihedral_dic[i[0] + "_N"] = np.asarray(dic_n[i[0]])[0]

    dic_cent, dic_ca, dic_n = defaultdict(
        list), defaultdict(list), defaultdict(list)
    for i in range(0, len(site2_coord)):
        info, name = pdb2_res_info[i]
        dic_cent[info].append(site2_coord[i])
        if name == "CA":
            dic_ca[info].append(site2_coord[i])
        if name == "N":
            dic_n[info].append(site2_coord[i])

    for i in dic_cent.items():
        arr2_dihedral_dic[i[0] + "_CN"] = np.mean(i[1], axis=0)
        arr2_dihedral_dic[i[0] + "_CA"] = np.asarray(dic_ca[i[0]])[0]
        arr2_dihedral_dic[i[0] + "_N"] = np.asarray(dic_n[i[0]])[0]

    arr, dihed_factor, site_check_sum1 = [], [], []
    for i in ans_dic.items():
        dihed_val = dihedral1(i[0], Counter(i[1]).most_common(1)[0][0])
        dihed_factor.append(dihed_val)
        if dihed_val > 245:
            arr.append([i[0] + " " + Counter(i[1]).most_common(1)[0][0]])

        site_check_sum1.append(
            residue_pairs_dictionary[i[0][:3]][Counter(i[1]).most_common(1)[0][0][:3]])

    site_check_sum1 = sum(site_check_sum1) * sum(dihed_factor)
    return site_check_sum1, arr


def SiteGen1(arr):

    arr1 = []
    dic = {1: "   ", 2: "  ", 3: " ", 4: ""}
    for i in arr:
        var = ""
        for j in i:
            j1 = "%.3f" % j
            var += dic[len(j1.split(".")[0])] + j1
        arr1.append(var)

    out = open("frag.pdb", 'w')
    for i1 in range(len(arr1)):
        i = arr1[i1]
        j = site1a[i1]
        out.write(j[:30] + i + j[54:] + "\n")
    out.close()

    out = open('fixed.pdb', 'w')
    for i in site2a:
        out.write(i + "\n")
    out.close()


def site_gen_het(site_gen_het):
    dic1, dic2 = {}, {}
    out = open("align", 'w')

    for i in site_gen_het:
        i1 = i[0].split(" ")
        out.write(i1[1] + " " + i1[0] + "\n")
        dic1[i1[0]] = 0
        dic2[i1[1]] = 0
    out.close()
    arr = B_all
    arr1 = []
    dic = {1: "   ", 2: "  ", 3: " ", 4: ""}
    for i in arr:
        var = ""
        for j in i:
            j1 = "%.3f" % j
            var += dic[len(j1.split(".")[0])] + j1
        arr1.append(var)

    out = open("site1.pdb", 'w')
    for i1 in range(0, len(arr1)):
        i = arr1[i1]
        j = site1a[i1]
        if j[:4] == "ATOM":
            res_info = j[17:20] + "-" + j[21:22] + "-" + j[22:26].strip()

            if res_info in dic1:

                out.write(j[:30] + i + j[54:] + "\n")
        else:
            out.write(j[:30] + i + j[54:] + "\n")
    out.close()

    out = open("site2.pdb", 'w')
    for i in site2a:
        if i[:4] == "ATOM":
            res_info = i[17:20] + "-" + i[21:22] + "-" + i[22:26].strip()

            if res_info in dic2:

                out.write(i + "\n")
        else:
            out.write(i + "\n")
    out.close()


def print_scores(arr):
    res1, res2, summer = [], [], []
    control1, control2 = [], []
    min_max, min_max2 = 0, []
    fin_val = ""
    for i in arr:
        i1, i2 = i[0].split(' ')
        res1.append(i2)
        res2.append(i1)
        i1, i2 = i1[:3], i2[:3]
        summer.append(residue_pairs_dictionary[i1][i2])
        control1.append(residue_pairs_dictionary[i1][i1])
        control2.append(residue_pairs_dictionary[i2][i2])
    res1 = sorted(set(res1))
    res2 = sorted(set(res2))
    min_max = len(arr)
    min_max2.append(pdb2_ln)
    min_max2.append(pdb1_ln)
    summer = float(sum(summer))
    if len(control1) == 0 or len(control2) == 0:
        return 'No Score'

    blosum_score = summer / ((sum(control1) / 2) + (sum(control2) / 2))

    fin_val = str(len(res1)) + "/" + str(pdb1_ln) + " " + str(len(res2)) + "/" + str(pdb2_ln) + " " + str(float(min_max) / min(min_max2)) + " " + str(float(min_max) / max(min_max2)) +\
        " " + str(blosum_score)
    return fin_val


def MainCode(aline, bline):
    logging.debug("Starting MainCode")
    aline = aline.split('__')  # __
    bline = bline.split('__')
    global res_dic1, res_dic2, res_arr1, res_arr2, res_pairs_dic1, res_pair2_dic2
    global dic_loop1, dic_loop2, dic_whole, SortedArrDic
    global dic_single1, dic_single2, site1a, site2a
    global pdb1_ln, pdb2_ln

    global res_arr1, res_arr2, res_dic1, res_dic2, res_pairs_dic1, res_pairs_dic2, B_all
    global pdb1_res_info, pdb2_res_info, site1_coord, site2_coord

    res_arr1, res_dic1, res_pairs_dic1, pdb1_lines, pdb1_het_lines = file_process(
        aline)
    res_arr2, res_dic2, res_pairs_dic2, pdb2_lines, pdb2_het_lines = file_process(
        bline)
    logging.debug("Starting run()")
    Final1, Final2 = run()
    logging.debug("Finished run()")
    NewArray = process_hits(Final1, Final2)
    if not NewArray:
        return 'None\tNone'
        # sys.exit()

    # Note bring coordinate to center before doing any opertation
    pdb1_trans_coord, pdb1_res_info, pdb1_generated_coord, pdb1_ca_dic = [
    ], [], [], defaultdict(list)
    pdb2_trans_coord, pdb2_res_info, pdb2_generated_coord, pdb2_ca_dic = [
    ], [], [], defaultdict(list)
    pdb1_ln, pdb2_ln = 0, 0
    logging.debug("Processing pdb1_lines")
    for line in pdb1_lines:
        if line[:4] == "ATOM":
            res1 = line[17:20] + "-" + line[21:22] + "-" + line[22:26].strip()
            pdb1_res_info.append([res1, line[13:16].strip()])
            pdb1_generated_coord.append(line)
            # if res1 == 'THR-A-70':
            #	print(res1, line)
            pdb1_trans_coord.append(
                [line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
            if line[13:16].strip() == "CA":
                pdb1_ln += 1
                pdb1_ca_dic[res1].append(line[28:38].strip())
                pdb1_ca_dic[res1].append(line[38:46].strip())
                pdb1_ca_dic[res1].append(line[46:54].strip())
    logging.debug("Finished pdb1_lines")
    logging.debug("Processing pdb2_lines")
    for line in pdb2_lines:
        if line[:4] == "ATOM":
            res1 = line[17:20] + "-" + line[21:22] + "-" + line[22:26].strip()
            pdb2_res_info.append(
                [line[17:20] + "-" + line[21:22] + "-" + line[22:26].strip(), line[13:16].strip()])
            pdb2_generated_coord.append(line)
            pdb2_trans_coord.append(
                [line[28:38].strip(), line[38:46].strip(), line[46:54].strip()])
            if line[13:16].strip() == "CA":
                pdb2_ln += 1
                pdb2_ca_dic[res1].append(line[28:38].strip())
                pdb2_ca_dic[res1].append(line[38:46].strip())
                pdb2_ca_dic[res1].append(line[46:54].strip())
    pdb1_trans_coord = np.asarray(pdb1_trans_coord, dtype='float')
    pdb2_trans_coord = np.asarray(pdb2_trans_coord, dtype='float')

    ResLists = []
    maxi, index_ln = 0, 0

    NewCount = 0
    for i in NewArray:
        site1_arr, site2_arr = [], []
        site1_coord, site2_coord = copy.deepcopy(
            pdb1_trans_coord), copy.deepcopy(pdb2_trans_coord)
        for j in range(len(i[0])):
            site1_arr.append(pdb1_ca_dic[i[0][j]])
            site2_arr.append(pdb2_ca_dic[i[1][j]])

        site1_arr, site2_arr = np.asarray(
            site1_arr, dtype=float), np.asarray(site2_arr, dtype=float)
        U = kabsch(copy.deepcopy(site1_arr), copy.deepcopy(site2_arr), True)
        B_all = copy.deepcopy(site1_coord)
        B_all -= site1_arr.mean(axis=0)
        B_all = np.dot(B_all, U)
        B_all += site2_arr.mean(axis=0)

        site1a = pdb1_generated_coord
        site2a = pdb2_generated_coord

        score, new_res_list = SiteGen()
        if score > maxi:
            maxi = score
            index_ln = len(i[0])
        else:
            if len(i[0]) < index_ln:
                NewCount += 1

        ResLists.append([new_res_list, score])
    ResLists = sorted(ResLists, key=lambda x: float(x[1]), reverse=True)
    line1 = ResLists[0][0]
    site1_arr, site2_arr = [], []

    MAPP_scores = print_scores(line1)
    MAPP_seqs = '_'.join([i[0] for i in line1])
    if not MAPP_seqs:
        MAPP_seqs = 'No Match'

    # return MAPP_scores+'\t'+MAPP_seqs
    logging.debug("finishing MainCode()")
    return MAPP_scores + '\t' + MAPP_seqs

# MPI CODE START


def pdb_res():
    res_dic = {}
    aline = open(pdb_size_file, 'r').readlines()
    for line in aline:
        line = line.strip()
        l = line.split("\t")

        res_dic[l[0]] = l[1]
    return res_dic


def s2():
    completed_alignment_dict = {}
    out = open(align_output_file, 'a+')
    out.close()
    aline = open(align_output_file, 'r')
    for line in aline:
        line = line.strip()
        completed_alignment_dict[line.split("\t")[0] + " " + line.split("\t")[1]] = 0
    return completed_alignment_dict


def chunk_mem_mpi(arr1, runnable_rank):
    new_arr1 = []
    for i in arr1:
        if int(i.split(" ")[2]) == runnable_rank:
            new_arr1.append(i)
    new_arr2 = []
    for i in new_arr1:
        aline = open(sites_folder + "/" + i.split(" ")[0], 'r').readlines()
        bline = open(sites_folder + "/" + i.split(" ")[1], 'r').readlines()
        ac = ""
        bc = ""
        for al in aline:
            al = al.strip()
            ac += al
            ac += "__"
        for bl in bline:
            bl = bl.strip()
            bc += bl
            bc += "__"
        new_arr2.append(i.replace(" ", "ttt") + "ttt" + ac[:-2] + "ttt" + bc[:-2])

    return new_arr2


def s1(completed_alignment_dict, res_dic):
    aline = open(paired_list, 'r').readlines()
    paired_arr = []
    for line in aline:
        line = line.strip()
        l = line.split("\t")
        paired_arr.append(line + "\t" + str(res_dic[l[0]] + res_dic[l[1]]))

    paired_arr = sorted(paired_arr, key=lambda x: int(
        x.split("\t")[2]), reverse=True)

    list_of_file_pairs = []
    file_pair_count = 0
    step_size = size * 10  # how many pairs to batch into each process
    for pdb_file_tsv_line in paired_arr:
        pdb_file1, pdb_file2, _ = pdb_file_tsv_line.split("\t")
        list_key = ' '.join([pdb_file1, pdb_file2])
        if list_key not in completed_alignment_dict:
            list_of_file_pairs.append(list_key)
            file_pair_count = len(list_of_file_pairs)

            if file_pair_count == step_size or file_pair_count == len(paired_arr):
                if rank == 0:
                    logging.info(f"Batch run has started")
                    count = 0
                    arr1 = []
                    for count_i in list_of_file_pairs:
                        logging.info(count_i)
                        count += 1
                        arr1.append(count_i + " " + str(count))
                        if count >= size:
                            count = 0
                    for runnable_rank in range(1, size + 1):
                        new_arr1 = []
                        new_arr1 = chunk_mem_mpi(arr1, runnable_rank)
                        comm.send(new_arr1, dest=runnable_rank)

                    out = open(align_output_file, 'a+')
                    for gettable_rank in range(1, size + 1):
                        ans = comm.recv(source=MPI.ANY_SOURCE)
                        logging.info(f"Batch run {gettable_rank} has finished")
                        for l2 in ans:
                            out.write(l2 + "\n")
                    time.sleep(.1)
                    out.close()

                elif rank != 0:
                    else_ans = []
                    while True:
                        else_arr = comm.recv(source=0)
                        else_ans = []
                        for line in else_arr:
                            ls = line.split("ttt")
                            start_time = time.time()

                            else_ans1 = 'SiteError'
                            try:
                                else_ans1 = MainCode(ls[3], ls[4])
                            except:
                                logging.error("MainCode() failed", exc_info=True)
                                MPI.COMM_WORLD.Abort()
                            else_ans.append(ls[0] + "\t" + ls[1] + "\t" + else_ans1)

                            end_time = time.time()
                            l_child = ls[:3]

                            msg = (
                                f"Pair {l_child[0]} {l_child[1]} was completed by processor --> "
                                f"{l_child[2]} in time {round((end_time-start_time), 2)} sec"
                            )
                            logging.info(msg)
                        comm.send(else_ans, dest=0)
                        # break

                        if len(else_arr) < 9:
                            break

            if file_pair_count >= step_size:
                list_of_file_pairs = []
                file_pair_count = 0
    if rank == 0:
        if len(list_of_file_pairs) < 1:
            print("\nYour job is completed!!!\n")
            MPI.COMM_WORLD.Abort()
        logging.debug('exit stage')
        arr11 = []
        count_1 = 0
        for count_i1 in list_of_file_pairs:
            count_1 += 1
            arr11.append(count_i1 + " " + str(count_1))
            if count_1 >= size:
                count_1 = 0

        for runnable_rank in range(1, size + 1):
            new_arr1 = []
            new_arr1 = chunk_mem_mpi(arr11, runnable_rank)
            comm.send(new_arr1, dest=runnable_rank)

        out = open(f"${output_dir}/align_output.txt", 'a+')
        for gettable_rank in range(1, size + 1):
            ans = comm.recv(source=MPI.ANY_SOURCE)
            logging.debug(ans)
            for l2 in ans:
                out.write(l2 + "\n")
        time.sleep(.1)
        out.close()
        logging.info("Alignment completed.")
        MPI.COMM_WORLD.Abort()

    elif rank != 0:
        else_ans = []
        while True:
            else_arr = comm.recv(source=0)
            else_ans = []
            for line in else_arr:
                ls = line.split("ttt")
                start_time = time.time()
                else_ans1 = 'SiteError'
                try:
                    else_ans1 = MainCode(ls[3], ls[4])
                except:
                    pass
                else_ans.append(ls[0] + "\t" + ls[1] + "\t" + else_ans1)

                end_time = time.time()
                l_child = ls[:3]

                logging.info((
                    f"Pair {l_child[0]} {l_child[1]} is completed from the processor --> "
                    f"{l_child[2]} in time {round((end_time-start_time), 2)} sec)"
                ))

            comm.send(else_ans, dest=0)
            if len(else_arr) < 9:
                break


if __name__ == "__main__":
    if len(sys.argv) == 5:
        sites_folder = sys.argv[1]
        paired_list = sys.argv[2]
        pdb_size_file = sys.argv[3]
        output_dir = sys.argv[4]
    else:
        print("python pocket_matrix_mpi.py <SiteFolder> <PairList.txt> <PDBSize.txt> ")
        sys.exit()

    align_output_file = f'{output_dir}/align_output.txt'

    blosum()
    global res_dic
    res_dic = pdb_res()

    global completed_alignment_dict
    completed_alignment_dict = s2()

    s1(completed_alignment_dict, res_dic)
