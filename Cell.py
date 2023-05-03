# Copyright (C) 2018-2023, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This section contains the basic data structure of the program.
# Author: Xiao Wang, 2018-12; Modified by Jie Li and Bin-Guang Ma, 2023-03. 

import numpy as np 

la = np.seterr(divide='ignore', invalid='ignore')

from scipy.ndimage import gaussian_filter1d

class Chromosomes:
    '''The chromosomes object contains the coordinates of each bin of the chromosomes structure.'''
    def __init__(self, intra_bin_num, seed, scale = 20.0):
        self.intra_bin_num = intra_bin_num # tuple
        self.seed = seed
        self.scale = scale
        self.coord_dict = {}
        for i in range(len(self.intra_bin_num)):
            self.each = self.intra_bin_num[i]
            self.CoordinateInit()  
            self.coord_dict[i] = self.coor
        self.GetCenter()

    def CoordinateInit(self):
        '''Randomly generates the coordinates of each bin using a random or specified seed.'''
        if self.seed != "auto":
            np.random.seed(seed = self.seed)
        self.coor = np.random.rand(self.each, 3) * self.scale
        self.coor = self.coor.astype(np.float32)
        
    def GetCenter(self):
        '''
        Get the coordinates of the chromosomes structure center, and transform coordinates
        so that the chromosomes center is at the origin of coordinate system.
        '''
        self.center = np.array([0.0, 0.0, 0.0])   
        for key in self.coord_dict:
            chr_coord = self.coord_dict[key]
            
            self.center += np.array([sum(chr_coord[:, 0]), sum(chr_coord[:, 1]), sum(chr_coord[:, 2])])
        
        self.center /= sum(self.intra_bin_num)
        
        for key in self.coord_dict:
            self.coord_dict[key] -= self.center

    def Smooth(self, smooth):
        '''Gaussian smoothing.'''
        for key in self.coord_dict:
            self.coord_dict[key][:, 0] = gaussian_filter1d(self.coord_dict[key][:, 0], smooth, axis = 0)
            self.coord_dict[key][:, 1] = gaussian_filter1d(self.coord_dict[key][:, 1], smooth, axis = 0)
            self.coord_dict[key][:, 2] = gaussian_filter1d(self.coord_dict[key][:, 2], smooth, axis = 0)
                   


