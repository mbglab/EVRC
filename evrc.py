#!user/bin/env python3 

# Copyright (C) 2018-2023, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This section is the main program.
# Author: Xiao Wang, 2018-12; Modified by Jie Li, Wei-Cheng Gu and Bin-Guang Ma, 2023-03. 

'''
Usage: python evrc.py -i INPUT-FILE -o OUTPUT-FILE [...]
'''

import Opt
import IO
import Cell
import Matrix
import Core
import Solver

# Get command line parameters
intra_chr_path, inter_chr_path, chr_num, o_file, smooth, alpha, \
iter_num, seed, thread, beta = Opt.GetArg()

# Get interaction frequency matrix
intra_bin_num, intra_bin_num_zero, if_intra_dict, if_inter_dict = IO.CheckFile(intra_chr_path, inter_chr_path, chr_num)
        
# Get distance matrix
dis_intra_dict, dis_inter_dict, max_dis_inter, thread \
    = Matrix.IfDis(intra_bin_num, if_intra_dict, if_inter_dict, \
                       smooth, alpha, thread)

# Get co-clustering coefficient
try:
    from CoreC import calc_coCC

    cc_intra_dict, cc_inter_dict = calc_coCC(dis_intra_dict, dis_inter_dict, intra_bin_num, chr_num, thread)
except:
    print("Dynamic Library CoreC not found. Run the command 'python setup.py build_ext --inplace' to create it.")
    print("Using pure python engine to solve it.")

    from Core import calc_coCC

    cc_intra_dict, cc_inter_dict = calc_coCC(dis_intra_dict, dis_inter_dict, intra_bin_num, chr_num)

# Get initial chromosomes structure coordinates
coordinate = Cell.Chromosomes(intra_bin_num, seed)

# Get optimal chromosomes structure coordinates
Solver.solver_c(coordinate, dis_intra_dict, dis_inter_dict, \
                iter_num, thread, cc_intra_dict, cc_inter_dict, beta)

if smooth != "auto":
    coordinate.Smooth(smooth)
coordinate.GetCenter()

# Get PDB file
IO.OutputChromosomes(coordinate.coord_dict, o_file, \
                 intra_bin_num_zero)
print("Solve Done!")






