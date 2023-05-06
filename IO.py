# Copyright (C) 2018-2023, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This section contains the input and output functions.
# Author: Xiao Wang, 2018-12; Modified by Jie Li and Bin-Guang Ma, 2023-03. 


import os
import sys
import numpy as np


def CheckFile(intra_chr_path, inter_chr_path, chr_num):
    '''Check and load input files.'''
    if os.path.exists(intra_chr_path):
        all_intra_chr = os.listdir(intra_chr_path)
        intra_bin_num = [0,] * chr_num  
        intra_bin_num_zero = [0,] * chr_num 
        intra_dict = {}
        for each in all_intra_chr:
            
            # each: chr1.txt chr2.txt
            # other: chr1_circle.txt
            file_path = intra_chr_path + "/" + each
            each = each.split('.')[0].split('_')[0]
            intra_bin_num[int(each[3:])-1], intra_dict[int(each[3:])], intra_bin_num_zero[int(each[3:])-1]  = ReadfileIntra(file_path)
    else:
        print ("error: no such file path of intrachromosomal contact matrix: %s" % intra_chr_path)
        sys.exit()
    if chr_num == 1:
        inter_dict = {}
        pass
    else:
        if os.path.exists(inter_chr_path):
            all_inter_chr = os.listdir(inter_chr_path)
            inter_dict = {}
            for each in all_inter_chr:
                # each: chr1_chr2.txt chr2_chr3.txt 
                # other: chr1_chr2_fork.txt
                file_path = inter_chr_path + "/" + each
                each = each.split('.')[0]
                if len(each)>=11:
                    chr1 = each.split('_')[0]
                    chr2 = each.split('_')[1]
                    each = chr1 + "_" + chr2
                    
                    
                inter_dict[each] = ReadfileInter(each, file_path, intra_bin_num, intra_bin_num_zero)
        else:
            print ("error: no such file path of interchromosomal contact matrix: %s" % inter_chr_path)
            sys.exit()

    return intra_bin_num, intra_bin_num_zero, intra_dict, inter_dict


def ReadfileIntra(file_path):
    '''
    Read interaction frequency matrix file of intra-chromosomes.
    Matrix file must be a text file. The format can be:
    ```
    97.40402117	15.07634967	6.004951709	3.802337111	3.338865904	2.989131799	 ...  2.449228619
    15.07634967	88.14832493	16.24036544	6.241608089	5.616626943	4.121033363	 ...  4.094041487
    6.004951709	16.24036544	74.80135268	19.14015601	6.81180746	4.659733536	 ...  4.540593717
    3.802337111	6.241608089	19.14015601	81.84660783	16.192136	6.67701795	 ...  5.221994188
    3.338865904	5.616626943	6.81180746	16.192136	77.08564297	15.67730956	 ...  6.461680865
    2.989131799	4.121033363	4.659733536	6.67701795	15.67730956	78.79137425	 ...  18.20757172
    ...			...			...			...			...			...			 ...  ...
    2.449228619	4.094041487	4.540593717	5.221994188	6.461680865	18.20757172	 ...  75.30837229
    ```
    or be:
    ```
    1	1	995.119511949
    1	2	528.471583395
    1	3	593.973763883
    1	4	221.270598604
    1	5	94.9033081901
    1	6	80.9735274218
    1	7	194.371570516
    1	8	42.404366
    ...
    100	100	940.06730994
    ```
    '''
    f_list=[]
    f = open(file_path ,"r")
    for eachline in f:
        if eachline: 
            pass 
        else:
            continue
        f_list.append(eachline.split())
    f_array = np.array(f_list, dtype = np.float32)
    
    f_shape = np.shape(f_array)
    if f_shape[0] == f_shape[1]:
        f_matrix = f_array 
        bin_num = f_shape[0]
        bin_num_min = 1
    else:
        bin_num_max = int(max(np.max(f_array[:,0]), np.max(f_array[:,1])) )
        bin_num_min = int(min(np.min(f_array[:,0]), np.min(f_array[:,1])) )
        bin_num = bin_num_max - bin_num_min + 1
        
        f_matrix = np.zeros((bin_num_max+1, bin_num_max+1), dtype=np.float32)
        for each in f_array:
            f_matrix[int(each[0]), int(each[1])] = each[2]
            f_matrix[int(each[1]), int(each[0])] = each[2]
            if int(each[0]) == int(each[1]):
                each[2] = 0.0    
        
        f_matrix =  f_matrix[bin_num_min:, bin_num_min:]
    
    return bin_num, f_matrix, bin_num_min-1
             
 
def ReadfileInter(each, file_path, intra_bin_num, intra_bin_num_zero):
    '''
    Read interaction frequency matrix file of inter-chromosomes.
    Matrix file must be a text file.
    '''
    f_list=[]
    f = open(file_path, "r")

    for eachline in f:
        if eachline:
            pass
        else:
            continue
        f_list.append(eachline.split())
    f.close()
    f_array = np.array(f_list, dtype=np.float32)
    
    f_shape = np.shape(f_array)
    chr1 = each.split('_')[0][3:]
    chr2 = each.split('_')[1][3:]
    
    bin_num1 = intra_bin_num[int(chr1)-1]
    bin_num2 = intra_bin_num[int(chr2)-1]
    
    bin_num1_zero = intra_bin_num_zero[int(chr1)-1]
    bin_num2_zero = intra_bin_num_zero[int(chr2)-1]
    
    f_matrix = np.zeros((int(bin_num1 + bin_num1_zero), int(bin_num2 + bin_num2_zero)), dtype=np.float32)
    for each in f_array:
        
        f_matrix[int(each[0]) - 1, int(each[1]) - 1] = each[2]
    
    f_matrix = f_matrix[bin_num1_zero:, bin_num2_zero:]
    
    return f_matrix


def OutputChromosomes(positions, o_file, intra_bin_num_zero):
    '''Save the reconstruction result as a .pdb file for each chromosome and a .pdb file for all chromosomes together.'''
    pdb_allchrom = o_file + "/chr_all.pdb"
    oall_file = open(pdb_allchrom, "w")
    oall_file.write("\n")    
    atom_num = 0
    
    for key in positions:
  
        pdb_file = o_file + "/chr" + str(key+1) + ".pdb"
        oo_file = open(pdb_file, "w")
        oo_file.write("\n")
        col4 = str(key+1)
        
        WritePDB(positions[key], oo_file, col4, atom_num, intra_bin_num_zero[key])
        oo_file.write("END")
        oo_file.close()
        
        WritePDB(positions[key], oall_file, col4, atom_num, intra_bin_num_zero[key])
        
        atom_num += len(positions[key])
        
    oall_file.write("END")
    oall_file.close()


def WritePDB(positions, oo_file, col4, atom_num, intra_bin_num_zero):    
    '''Save the reconstruction result as a .pdf file.'''
    
    col1 = "ATOM"
    col3 = "BIN  C" + col4

    bin_num = len(positions)
    chrtype = 1
    
    
    for i in range(1, bin_num+1):
        col2 = str(atom_num + i)
        col5 = str(i + intra_bin_num_zero)
        col6 = "%.3f" % positions[i-1][0]
        col7 = "%.3f" % positions[i-1][1]
        col8 = "%.3f" % positions[i-1][2]
        col2 = " "*(5 - len(col2)) + col2
        col5 = " " * (4 - len(col5)) + col5
        col6 = " " * (8 - len(col6)) + col6
        col7 = " " * (8 - len(col7)) + col7
        col8 = " " * (8 - len(col8)) + col8

        col = (col1, col2, col3, col4, col5, col6, col7, col8)
        line = "%s  %s  %s %s%s   %s%s%s\n" % col
        oo_file.write(line)
    col1 = "CONECT"
    for i in range(1, bin_num+1):
        col2 = str(atom_num + i)
        j = atom_num + i + 1
        if j > atom_num + bin_num:
            if chrtype == 1:     
                continue
            j = atom_num + 1
        col3 = str(j)

        col2 = " " * (5 - len(col2)) + col2
        col3 = " " * (5 - len(col3)) + col3

        line = "%s%s%s\n" % (col1, col2, col3)
        oo_file.write(line)

