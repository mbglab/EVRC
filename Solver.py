# Copyright (C) 2018-2023, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This section contains the basic data structure of the program.
# Chromosome Class: Defining the coordinates of chromosomes and some basic transformations.
# Author: Xiao Wang, 2018-12; Modified by Jie Li and Bin-Guang Ma, 2023-03. 

import numpy as np

def solver_p(coordinate, dis_intra_dict, dis_inter_dict, iter_num, cc_intra_dict, cc_inter_dict, beta):
    '''
    Use pure python code to optimize the structure.
    Single-thread solution.
    '''
    from Core import evrc
    coord_dict = coordinate.coord_dict
    i = 0 
    f0 = 0
    while True:
        
        f = 0
        for key1 in dis_intra_dict:
            trans1 = evrc(coord_dict[key1-1], coord_dict[key1-1], dis_intra_dict[key1], cc_intra_dict[key1], beta)
            trans = trans1
            
            if dis_inter_dict  == {}:
                pass
            else:
                for key2 in dis_inter_dict:
                    if int(key2[3]) == key1:
                        
                        trans_inter = evrc(coord_dict[key1-1], coord_dict[int(key2[-1])-1], dis_inter_dict[key2], cc_inter_dict[key2], beta)
                        trans += trans_inter

                    if int(key2[-1]) == key1:
                        
                        trans_inter = evrc(coord_dict[key1-1], coord_dict[int(key2[3])-1], np.transpose(dis_inter_dict[key2]), np.transpose(cc_inter_dict[key2]), beta)
                        trans += trans_inter

                trans = (trans - trans1) / (len(dis_intra_dict)-1) + trans1

            for j in range(trans.shape[0]):
                f += (trans[j, 0]**2 + trans[j, 1]**2 + trans[j, 2]**2) ** 0.5
            
            coord_dict[key1-1] += trans
            
        coordinate.coord_dict = coord_dict
        i += 1
        if iter_num != "auto":
            if i > int(iter_num):
                break
        else:
            if abs(f - f0) < 0.000001:
                break
            f0 = f
        print (i, f)  
        

def solver_c(coordinate, dis_intra_dict, dis_inter_dict, iter_num, thread, cc_intra_dict, cc_inter_dict, beta):
    '''
    Optimize the structure with Cython and openmp.
    This function can only run on CPUs.
    This function requires the .pyx file to be compiled in advance.
    '''
    try:
        from CoreC import evrc
    except:
        solver_p(coordinate, dis_intra_dict, dis_inter_dict, iter_num, cc_intra_dict, cc_inter_dict, beta)

        return
    coord_dict = coordinate.coord_dict
    if thread == "auto":
        thread = 0
    ii = 0 
    f0 = 0
    while True:
        f = 0
        for key1 in dis_intra_dict:
            trans1 = evrc(coord_dict[key1-1], coord_dict[key1-1], dis_intra_dict[key1], \
                              thread, cc_intra_dict[key1], beta)
            trans1 = np.matrix(trans1, np.float32)
            trans = trans1
            if dis_inter_dict == {}:
                pass
            else:
                for key2 in dis_inter_dict:
                    key_first = key2.split('_')[0][3:]
                    key_later = key2.split('_')[1][3:]
                    
                    if int(key_first ) == key1:
                        
                        trans_inter = evrc(coord_dict[key1-1], coord_dict[int(key_later )-1], dis_inter_dict[key2], \
                                               thread, cc_inter_dict[key2], beta)
                        trans_inter = np.matrix(trans_inter, np.float32)      
                        trans += trans_inter
                        
                    if int(key_later ) == key1:
                        
                        trans_inter = evrc(coord_dict[key1-1], coord_dict[int(key_first )-1], \
                                               np.transpose(dis_inter_dict[key2]), \
                                               thread, np.transpose(cc_inter_dict[key2]), beta)
                        trans_inter = np.matrix(trans_inter, np.float32) 
                        trans += trans_inter
                        
                trans = (trans - trans1)/(len(dis_intra_dict)-1) + trans1
            for i in range(trans.shape[0]):
                f += (trans[i, 0]**2 + trans[i, 1]**2 + trans[i, 2]**2) ** 0.5
                
            coord_dict[key1-1] += trans   
        coordinate.coord_dict = coord_dict
        ii += 1
        if iter_num != "auto":
            if ii > int(iter_num):
                break 
        else:
            if abs(f - f0) < 0.000001:
                break
            f0 = f
    
        print (ii, f)
        
