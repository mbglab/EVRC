# Copyright (C) 2018-2023, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This section contains the core procedures in Cython.
# Author: Xiao Wang, 2018-12; Modified by Jie Li, Wei-Cheng Gu and Bin-Guang Ma, 2023-03. 

import cython
import numpy as np
from Matrix import *

# Cython: boundscheck = False
from cython.parallel import prange
from libc.math cimport sqrt
cimport numpy as np
cimport openmp

@cython.wraparound(False)
@cython.boundscheck(False)
def calc_coCC(dict dis_intra_dict, dict dis_inter_dict, list intra_bin_num, int chr_num, int thread):
    '''Calculate the co-clustering coefficient matrixes.'''
    #chr1_chr2 chr1_chr3 not: chr2_chr1

    all_matrix = get_all_matrix(dis_intra_dict, dis_inter_dict, intra_bin_num, chr_num)
    cc_ij_matrix = np.zeros((all_matrix.shape[0], all_matrix.shape[0]), np.float32)
    not_zeros = np.zeros((all_matrix.shape[0], all_matrix.shape[0]), np.int32)

    cdef :
        int row_num = all_matrix.shape[0]
        int col_num = all_matrix.shape[1]
        float[:,:] all_matrix_view = all_matrix
        float[:,:] cc_ij_matrix_view = cc_ij_matrix
        int[:,:] not_zeros_view = not_zeros

        int i, j, k, m, o, p
        float k_ij, e_ij, cc_ij, kmatrix_element, ematrix_element
        dict cc_intra_dict, cc_inter_dict

    if thread == 0:
        thread = openmp.omp_get_num_procs()

    for i in prange(row_num, nogil=True, num_threads=thread):
        for j in xrange(col_num):
            if all_matrix_view[i, j]:
                m = 0
                k_ij, e_ij, kmatrix_element, ematrix_element = 0.0, 0.0, 0.0, 0.0

                # k_ij
                for k in range(col_num):
                    kmatrix_element = all_matrix_view[i, k] + all_matrix_view[j, k]

                    if kmatrix_element != 0.0:
                        not_zeros_view[i, m] = k
                        m = m + 1

                        k_ij = k_ij + 1.0

                # e_ij
                for o in range(m):
                    for p in range(m):
                        ematrix_element = all_matrix_view[not_zeros_view[i, o], not_zeros_view[i, p]]

                        if ematrix_element != 0.0:
                            e_ij = e_ij + 1.0

                e_ij = e_ij / 2.0

                # cc_ij
                cc_ij = 2.0 * e_ij / k_ij / (k_ij - 1)
                cc_ij_matrix_view[i, j] = cc_ij

    cc_intra_dict, cc_inter_dict = matrix_to_dict(cc_ij_matrix, intra_bin_num, chr_num)

    return cc_intra_dict, cc_inter_dict


def evrc(float[:, :] position1, float[:, :] position2, float[:, :] dis_matrix, int thread, float[:,:] cc_matrix, float beta):
    '''
    Calculate the transformation matrix and F value based on
    the co-Clustering-Coefficient matrix and distance matrix using
    error vector resultant.
    '''
    cdef:
        int bin_num1 = dis_matrix.shape[0]
        int bin_num2 = dis_matrix.shape[1]

        float[:, :] trans_matrix = np.zeros((bin_num1, 3), np.float32)  
        
        float pos_ix, pos_iy, pos_iz
        float divec_x, divec_y, divec_z
        float ervec_x, ervec_y, ervec_z
        float pos_jx, pos_jy, pos_jz
        
        float mod, dis, err
        int i, j
        float num
        
    if thread == 0:
        thread = openmp.omp_get_num_procs()
        
    for i in prange(bin_num1, nogil=True, num_threads=thread):
        pos_ix, pos_iy, pos_iz = position1[i, 0], position1[i, 1], position1[i, 2]
        
        ervec_x, ervec_y, ervec_z = 0.0, 0.0, 0.0
        divec_x, divec_y, divec_z = 0.0, 0.0, 0.0
      
        num = 0.0 
        
        for j in xrange(bin_num2):
           
            pos_jx, pos_jy, pos_jz = position2[j, 0], position2[j, 1], position2[j, 2]
            divec_x = pos_jx - pos_ix
            divec_y = pos_jy - pos_iy
            divec_z = pos_jz - pos_iz
            
            mod = sqrt(divec_x**2 + divec_y**2 + divec_z**2)
            dis = dis_matrix[i, j]
            
            if mod == 0:
                continue
                
            divec_x = divec_x / mod
            divec_y = divec_y / mod
            divec_z = divec_z / mod
            
            if dis == 0.0 :
                err = 0.0
            
            else:
                err = mod - dis             
            
            divec_x = divec_x * err * cc_matrix [i, j] / (bin_num2**beta)
            divec_y = divec_y * err * cc_matrix [i, j] / (bin_num2**beta)
            divec_z = divec_z * err * cc_matrix [i, j] / (bin_num2**beta)
            
            ervec_x = ervec_x + divec_x
            ervec_y = ervec_y + divec_y
            ervec_z = ervec_z + divec_z
            
            if dis != 0.0:
                num = num + 1.0
       
        if num != 0:
            ervec_x = ervec_x / num
            ervec_y = ervec_y / num
            ervec_z = ervec_z / num
       
        trans_matrix[i, 0] += ervec_x
        trans_matrix[i, 1] += ervec_y
        trans_matrix[i, 2] += ervec_z       
     
    return trans_matrix
