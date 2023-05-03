# Copyright (C) 2018-2023, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This section contains the core procedures.
# Author: Xiao Wang, 2018-12; Modified by Jie Li and Bin-Guang Ma, 2023-03. 

import numpy as np
from math import sqrt
from Matrix import *

def  calc_coCC(dis_intra_dict, dis_inter_dict, intra_bin_num, chr_num):
    '''Calculate the co-clustering coefficient matrixes.'''
    # chr1_chr2 chr1_chr3 not: chr2_chr1
    all_matrix = get_all_matrix(dis_intra_dict, dis_inter_dict, intra_bin_num, chr_num)
        
    row_num = np.shape (all_matrix)[0]
    col_num = np.shape (all_matrix)[1]
    cc_ij_matrix = np.zeros((row_num, col_num), np.float32)
    for i in range(row_num):
        for j in range(col_num):
            if all_matrix [i, j] :
                # k_ij
                k_ij_matrix = all_matrix [i, :] + all_matrix [j, :]
                not_zeros =(k_ij_matrix != 0.0)
                
                k_ij_matrix [not_zeros ] = 1.0
                k_ij = np.sum(k_ij_matrix )
                
                k_ij_list = not_zeros.tolist()
                ij_all_matrix = all_matrix [k_ij_list, :][:, k_ij_list]
                # e_ij
                not_zeros =(ij_all_matrix != 0.0)
                ij_all_matrix [not_zeros ] = 1.0
                # if have smooth: e_ij = (np.sum(ij_all_matrix ) - k_ij) / 2.0
                e_ij = (np.sum(ij_all_matrix )) / 2.0
                # cc_ij = 1.0
                cc_ij = 2.0 * e_ij / k_ij / (k_ij - 1 )
                cc_ij_matrix [i, j] = cc_ij
    cc_intra_dict, cc_inter_dict = matrix_to_dict(cc_ij_matrix, intra_bin_num, chr_num)
    return cc_intra_dict, cc_inter_dict


def evrc(position1, position2, dis_matrix, cc_matrix, beta):
    '''
    Calculate the transformation matrix and F value based on
    the co-Clustering-Coefficient matrix and distance matrix using 
    error vector resultant.
    '''
    
    bin_num1, bin_num2 = np.shape(dis_matrix)
    trans_matrix = np.zeros((bin_num1, 3), np.float32)
    
    f_value = 0.0
    for i in range(bin_num1):
        pos_i = position1[i]
        
        divec = np.array([0.0, 0.0, 0.0], dtype=np.float)
        ervec = np.array([0.0, 0.0, 0.0], dtype=np.float)
        
        for j in range(bin_num2):
            
            pos_j = position2[j]
            divec = pos_j - pos_i

            mod = sqrt(divec[0]**2 + divec[1]**2 + divec[2]**2)
            dis = dis_matrix[i, j]
            
            if mod == 0:
                continue
            divec = divec / mod        

            if dis == 0:
                err = 0
            else:
                err = mod - dis
            
            divec = divec * err * cc_matrix [i, j] / (bin_num2**beta)
            ervec = ervec + divec
        
        not_zeros = dis_matrix[i,] != 0.0 
        not_zeros_sum = np.sum(np.int64(not_zeros))
        if not_zeros_sum != 0 :
            ervec /=  not_zeros_sum
            trans_matrix[i] += ervec
    
    return trans_matrix
    
