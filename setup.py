# Copyright (C) 2018-2023, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This section is for compiling Cython.
# Author: Xiao Wang, 2018-12; Modified by Jie Li, Wei-Cheng Gu and Bin-Guang Ma, 2023-03. 

'''
Compile the CORE.pyx file into an executable file.
Compile command: python setup.py build_ext --inplace.
Before compiling, you need to install the cython library.
'''

import sys
import numpy as np

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


if "win" in sys.platform:
    ext_modules = [
        Extension(
            "CoreC",
            ["CoreC.pyx"],
            extra_compile_args=['/openmp'],
            extra_link_args=['/openmp']
        )
    ]
else:
    ext_modules = [
        Extension(
            "CoreC",
            ["CoreC.pyx"],
            extra_compile_args=['-fopenmp'],
            extra_link_args=['-fopenmp'],
            libraries=["m"]
        )
    ]


setup(
    name='Core',
    ext_modules=cythonize(ext_modules, annotate = True, language_level = "3"),include_dirs=[np.get_include()])

