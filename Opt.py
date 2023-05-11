# Copyright (C) 2018-2023, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This section contains the functions for processing command line options.
# Author: Xiao Wang, 2018-12; Modified by Jie Li and Bin-Guang Ma, 2023-03. 

import sys
import getopt

def GetArg():
    '''
    Get command line parameters.
    '''
    
    intra_chr_path = ''
    inter_chr_path = ''
    chr_num = ''
    o_file = ''
        
    smooth = "auto"
    alpha = 0.5
    scale = 20.0
    iter_num = "auto"
    seed = "auto"
    thread = "auto"
    beta = 0.1
    
    try:
        options, args = getopt.getopt(sys.argv[1:], "ho:", \
                                      ["intra_chr=", "inter_chr=", "chr_num=",  \
                                       "help", "sf=", "alpha=", \
                                       "iter_num=", "seed=", "thread=", "beta="])

        for opt_name, opt_value in options:
            
            if opt_name in ("--intra_chr"):
                intra_chr_path = opt_value
                
            elif opt_name in ("--inter_chr"):
                inter_chr_path = opt_value
                
            elif opt_name in ("--chr_num"):
                chr_num = int(opt_value)
                if chr_num <= 0:
                    print ("error: wrong parameter in --chr_num")
                    ShowHelp()
                    sys.exit()                           
                
            elif opt_name in ("-o"):
                o_file = opt_value
            
            elif opt_name in ('-h','--help'):
                ShowHelp()
                sys.exit()
                
            elif opt_name in ("--sf"):
                if opt_value != "auto":
                    smooth = float(opt_value)
                
            elif opt_name in ("--alpha"):
                alpha = float(opt_value)
                
            elif opt_name in ("--iter_num"):
                if opt_value != "auto":          
                    iter_num = int(opt_value)    
         
            elif opt_name in ("--seed"):
                if opt_value != "auto":
                    seed = int(opt_value)
                
            elif opt_name in ("--thread"):
                if opt_value != "auto":
                    thread = int(opt_value)

            elif opt_name in ("--beta"):
                beta = float(opt_value)
                    
                                
        if not intra_chr_path:
            print ("error: no file path of intrachromosomal interaction matrix")
            sys.exit() 
        if chr_num != 1:           
            if not inter_chr_path:
                print ("error: no file path of interchromosomal interaction matrix")
                sys.exit()
        if not chr_num:
            print ("error: no the number of chromosomes")
            sys.exit()
        if not o_file:
            print ("error: no output file")
            sys.exit()
            
        return intra_chr_path, inter_chr_path, chr_num, o_file, smooth, alpha, iter_num, seed, thread, beta

    except TypeError as e:
        print (e)
        sys.exit()
        
        

def ShowHelp():
    '''
    Print help information to the terminal.
    '''
    
    help_info = ''' 
Usage: 

python evrc.py --intra_chr FilePath --inter_chr FilePath --chr_num Chromosome_Number -o OutputFilePath \
--sf Smoothing_factor --alpha Transform_Exponent  \
--iter_num Iteration_Number --seed Seed --thread Thread_Number --beta Convergence_Factor

Required parameters: 

--intra_chr   file path of the intrachromosomal contact matrix 
--chr_num     the number of chromosomes 
-o            output file path

Optional prameters: 

--inter_chr     file path of the interchromosomal contact matrix 
--help          show help (also -h) 
--sf            smoothing factor (default: "auto". Do not smooth the structure)
--alpha         the transformation exponent from interaction matrix into distance matrix (default: 0.5)
--iter_num      the number of iterations (default: "auto". The program will terminate when F value is less than 1E-6)
--seed          seed for generating the random initial 3D coordinates (default: "auto")
--thread        the number of threads (default: "auto". All the CPU cores will be used)
--beta          the convergence factor(default: 0.1)
''' 

    print (help_info)
    
