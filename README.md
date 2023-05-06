# EVRC
The EVRC algorithm first calculates the co-clustering coefficients between chromatin segments, and then sums up the error vectors of segments to reconstruct the 3D structure of chromatin through continuous iterative optimization.

# Installation
Requirement:
*Python 3.x
*numpy
*scipy

Optinal (Recommended):
*Cython

# Usage
```bash
python evrc.py --intra_chr FilePath --inter_chr FilePath --chr_num Chromosome_Number -o OutputFilePath \
--sf Smooth_factor --alpha Transform_Exponent --scale Chromosomes_Size  \
--iter_num Iteration_Number --seed Seed --thread Thread_Number --beta Convergence_Factor
```

Required parameters: 

--intra_chr     file path of intrachromosomal interaction matrix  
--chr_num       the number of chromosomes  
-o              output file path  

Optional prameters: 

--inter_chr     file path of interchromosomal interaction matrix  
--help          show help (also -h)  
--sf            smooth factor (default: "auto". Do not smooth the structure)  
--alpha         the transformation exponent from interaction matrix into distance matrix (default: 0.5)  
--scale         the size of chromosomes (default: 20.0)  
--iter_num      the number of iterations (default: "auto". The program will terminate when F value is less than 1E-6)  
--seed          seed for generating random initial 3D coordinates (default: "auto")  
--thread        the number of threads (default: "auto". All the CPU cores will be used)  
--beta          convergence factor(default: 0.1; The higher the beta value, the slower the convergence speed. The beta value should be increased if it does not converge)  



# Input file format
1. Matrix format
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

2. Three columns format
The format of the input file can be three columns, separated by tabs:
```
bin1	bin2	frequency
```
Example:
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

For the 3D structure reconstruction of a single chromosome, an input file directory is required, which contains an interaction frequency file of the chromosome, and the prefix of this file name should be 'chr1_'.

For the 3D structure reconstruction of multiple chromosomes, it is necessary to provide a file directory for intra-chromosome interaction and a file directory for inter-chromosome interaction.  
The file directory of intra-chromosome interaction contains the file name prefix of chr1_, chr2_, etc.  
The file directory of inter-chromosome interaction contains the file name prefix of chr1_chr2, chr1_chr3, chr2_chr3, etc.  

# Output file format
The output file is a '.pdb' file that follows the standard PDB format. You can use any PDB viewer program to draw the structure.


# Speed up calculation with Cython
The program uses Cython for acceleration by default. Iterative optimization code is written in Cython syntax (see CoreC.pyx file). 
When running the program for the first time, you need to compile CoreC.pyx into a dynamic link library by this command:
```bash
python setup.py build_ext --inplace
```
Before compiling, you need to install Cython by "pip install cython" and something alike. 
This will produce a CoreC.so file on a Linux system or a CoreC.pyd file on a Windows system.
