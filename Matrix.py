# Copyright (C) 2018-2023, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This section contains the basic classes and functions for data processing.
# Author: Xiao Wang, 2018-12; Modified by Jie Li and Bin-Guang Ma, 2023-03. 

import numpy as np 

class Matrix:
    '''The matrix object involves interaction frequency (IF) matrix processing and transforming into the expected distance matrix.'''
    def __init__(self, if_matrix, bin_num1, bin_num2):
        self.if_matrix = if_matrix
        self.bin_num1 = bin_num1
        self.bin_num2 = bin_num2
        self.dis_matrix = np.zeros([self.bin_num1, self.bin_num2], np.float32)

    def IF2Dis(self, alpha):
        '''Convert the  matrix to distance matrix.'''
        self.dis_matrix = np.zeros([self.bin_num1, self.bin_num2], np.float32)
       
        less_1 = self.if_matrix > 0
        self.dis_matrix[less_1 ] = (1.0 / self.if_matrix[less_1]) ** alpha

    def SizeChange(self, scale1):
        prop = 1000000.0 / scale1 
        self.dis_matrix *= prop
        
    def DisFilterIntra(self):
        '''This function involves the process of handling default values.'''
        dis_neighbor = 0
        count = 0
        for i in range(np.shape(self.dis_matrix)[0]):
            self.dis_matrix[i, i] = 0.0
            if i < np.shape(self.dis_matrix)[0]-1 and self.dis_matrix[i, i+1] > 0:
                dis_neighbor += self.dis_matrix[i, i+1]
                count += 1
        dis_neighbor = dis_neighbor / count
        for i in range(0, np.shape(self.dis_matrix)[0] - 1):
            if self.dis_matrix[i, i+1] == 0.0:
                self.dis_matrix[i, i+1] = dis_neighbor
                self.dis_matrix[i+1, i] = dis_neighbor

def IfDis(intra_bin_num, if_intra_dict, if_inter_dict, \
           smooth, alpha, thread):
    '''Convert the IF matrix to the expected distance matrix.'''
    if thread == "auto":
        thread = 0   
    max_dis_intra = [0,] * len(intra_bin_num) 
    sum_dis = 0
    dis_intra_dict={}
    
    for key in if_intra_dict:
        
        bin_num = intra_bin_num[key-1]
        intra_matrix = Matrix(if_intra_dict[key], bin_num, bin_num)
        intra_matrix.IF2Dis(alpha)   
        max_dis_intra[key-1] = np.max(intra_matrix.dis_matrix)
        dis_intra_dict[key] = intra_matrix.dis_matrix  
        sum_dis += np.sum(intra_matrix.dis_matrix)
 
    max_dis_inter = []
    dis_inter_dict={}
    
    if if_inter_dict == {}:
        max_dis_dis = max(max_dis_intra)
        max_dis_inter = 0
    else:
        for key in if_inter_dict:
            key1 = key.split('_')[0][3:]
            key2 = key.split('_')[1][3:]
            #print "if_inter_dict: ", if_inter_dict[key]
            bin_num1 = intra_bin_num[int(key1)-1]       
            bin_num2 = intra_bin_num[int(key2)-1] 
            inter_matrix = Matrix(if_inter_dict[key], bin_num1, bin_num2)
            inter_matrix.IF2Dis(alpha)
            max_dis_inter.append(np.max(inter_matrix.dis_matrix))
            dis_inter_dict[key] = inter_matrix.dis_matrix
            sum_dis += np.sum(inter_matrix.dis_matrix)

        max_dis_dis = max(max(max_dis_intra), max(max_dis_inter))
        max_dis_inter = max(max_dis_inter )
    
    for key in dis_intra_dict:
        
        bin_num = intra_bin_num[key-1]
        intra_matrix = Matrix(0, bin_num, bin_num)
        intra_matrix.dis_matrix = dis_intra_dict[key]
        intra_matrix.SizeChange(sum_dis)
        intra_matrix.DisFilterIntra()
        dis_intra_dict[key] = intra_matrix.dis_matrix        

    if dis_inter_dict == {}:
        pass
    else: 
        for key in dis_inter_dict:
            key1 = key.split('_')[0][3:]
            key2 = key.split('_')[1][3:]
            bin_num1 = intra_bin_num[int(key1)-1]       
            bin_num2 = intra_bin_num[int(key2)-1] 
            inter_matrix = Matrix(0, bin_num1, bin_num2)
            inter_matrix.dis_matrix = dis_inter_dict[key]
            inter_matrix.SizeChange(sum_dis)
            dis_inter_dict[key] = inter_matrix.dis_matrix     

    return dis_intra_dict, dis_inter_dict, max_dis_inter, thread

def get_all_matrix(dis_intra_dict, dis_inter_dict, intra_bin_num, chr_num):
    '''Merge the separated matrices into a whole matrix.'''
    if chr_num == 1:
        all_matrix = dis_intra_dict[1]
    else:
        matrix_i_dict = {}
        for i in range(1,chr_num+1):
            intra_matrix = dis_intra_dict [i]
            matrix_i = intra_matrix 
            #ahead matrix
            for ahead in reversed(range(1, i)):
                inter_name ='chr' + str(ahead) + '_chr' + str(i) 
                inter_matrix = np.transpose(dis_inter_dict[inter_name])
                matrix_i = np.hstack((inter_matrix, matrix_i))
            #later matrix   
            for later in range(i+1, chr_num+1):
                inter_name = 'chr' + str(i) + '_chr' + str(later) 
                inter_matrix = dis_inter_dict[inter_name]
                matrix_i = np.hstack((matrix_i, inter_matrix))
            matrix_i_dict[i] = matrix_i 
        all_matrix = matrix_i_dict[1]
        for i in range(2, chr_num+1):
            all_matrix = np.vstack((all_matrix, matrix_i_dict[i]))
    return all_matrix

def matrix_to_dict(cc_ij_matrix, intra_bin_num, chr_num):
    '''Assign the co-clustering coefficient matrix to dict.'''
    cc_ij_intra_dict={}
    cc_ij_inter_dict={}
        
    if chr_num == 1:
        cc_ij_intra_dict [1] = cc_ij_matrix 
    else:
        for i in range(1, chr_num+1):
            for j in range(1, chr_num+1):
                if i == j:        
                    index2 = sum(intra_bin_num[0:i])
                    index1 = index2 - intra_bin_num[i-1]
                        
                    cc_ij_intra_dict[i] = cc_ij_matrix[index1:index2, index1:index2]
                
                else:
                    inter_name ='chr' + str(i) + '_chr' + str(j)
                        
                    index1_2 = sum(intra_bin_num[0:i])
                    index1_1 = index1_2 - intra_bin_num[i-1]
                    index2_2 = sum(intra_bin_num[0:j])
                    index2_1 = index2_2 - intra_bin_num[j-1]
                            
                    cc_ij_inter_dict[inter_name] = cc_ij_matrix[index1_1:index1_2, index2_1:index2_2]
    return cc_ij_intra_dict, cc_ij_inter_dict 


