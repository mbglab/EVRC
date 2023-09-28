import os
import math
import random
import numpy as np


def write_pdb_double_chr(positions1, positions2, structure, chrtype):
    '''
    Generate simulating structure of double chromosomes.
    '''
    col1 = "ATOM"
    col3 = "BIN  C" + "1"
    col4 = 1
    col9 = str(0) + " " + str(0)
    bin_num1 = len(positions1)
    bin_num2 = len(positions2)
    atom_num = 0

    os.mkdir(structure)
    oo_file = open(f"{structure}/{structure}.pdb", "w")
    for i in range(1, bin_num1 + 1):
        col2 = str(atom_num + i)
        col5 = str(i)
        col6 = "%.1f" % positions1[i - 1][0]
        col7 = "%.1f" % positions1[i - 1][1]
        col8 = "%.1f" % positions1[i - 1][2]
        col2 = " " * (5 - len(col2)) + col2
        col5 = " " * (4 - len(col5)) + col5
        col6 = " " * (8 - len(col6)) + col6
        col7 = " " * (8 - len(col7)) + col7
        col8 = " " * (8 - len(col8)) + col8

        col = (col1, col2, col3, col4, col5, col6, col7, col8, col9)
        line = "%s  %s  %s %s%s   %s%s%s  %s\n" % col
        oo_file.write(line)

    col1 = "CONECT"
    for i in range(1, bin_num1 + 1):
        col2 = str(atom_num + i)
        j = atom_num + i + 1
        if j > atom_num + bin_num1:
            if chrtype == 1:
                continue
            j = atom_num + 1

        col3 = str(j)
        col2 = " " * (5 - len(col2)) + col2
        col3 = " " * (5 - len(col3)) + col3
        line = "%s%s%s\n" % (col1, col2, col3)
        oo_file.write(line)

    atom_num = bin_num1
    col1 = "ATOM"
    col3 = "BIN  C" + "2"
    col4 = 2
    for i in range(1, bin_num2 + 1):
        col2 = str(atom_num + i)
        col5 = str(atom_num + i)
        col6 = "%.1f" % positions2[i - 1][0]
        col7 = "%.1f" % positions2[i - 1][1]
        col8 = "%.1f" % positions2[i - 1][2]
        col2 = " " * (5 - len(col2)) + col2
        col5 = " " * (4 - len(col5)) + col5
        col6 = " " * (8 - len(col6)) + col6
        col7 = " " * (8 - len(col7)) + col7
        col8 = " " * (8 - len(col8)) + col8

        col = (col1, col2, col3, col4, col5, col6, col7, col8, col9)
        line = "%s  %s  %s %s%s   %s%s%s  %s\n" % col
        oo_file.write(line)

    col1 = "CONECT"
    for i in range(1, bin_num2 + 1):
        col2 = str(atom_num + i)
        j = atom_num + i + 1
        if j > atom_num + bin_num2:
            if chrtype == 1:
                continue
            j = atom_num + 1

        col3 = str(j)
        col2 = " " * (5 - len(col2)) + col2
        col3 = " " * (5 - len(col3)) + col3
        line = "%s%s%s\n" % (col1, col2, col3)
        oo_file.write(line)

    oo_file.write("END")
    oo_file.close()

def add_random_noise2(positions1, positions2, structure):
    '''
    Generate the interaction frequency (IF) matrix of double chromosome simulating structure
    with different levels of random noises.
    '''
    for k in range(1, 11): # Repeat 10 times
        bin_num1 = len(positions1)
        bin_num2 = len(positions2)

        for s in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]: # different noise level
            os.makedirs(f"{structure}/{k}/s{s}/intra", exist_ok= True)
            # Chr1
            oo_file = open(f"{structure}/{k}/s{s}/intra/chr1_{structure}_if_s{s}_{k}.txt", "w")
            for i in range(bin_num1):
                for j in range(i):
                    divec = positions1[i, :] - positions1[j, :]
                    dis = math.sqrt(divec[0] ** 2 + divec[1] ** 2 + divec[2] ** 2)
                    r = random.uniform(-1, 1)
                    dis *= (1 + r * s)
                    if dis == 0:
                        line = str(i + 1) + "\t" + str(j + 1) + "\t" + str(dis) + "\n"
                    else:
                        line = str(i + 1) + "\t" + str(j + 1) + "\t" + str(10000.0 / dis) + "\n"

                    oo_file.write(line)
            oo_file.close()
            # Chr2
            oo_file = open(f"{structure}/{k}/s{s}/intra/chr2_{structure}_if_s{s}_{k}.txt", "w")
            for i in range(bin_num2):
                for j in range(i):
                    divec = positions2[i, :] - positions2[j, :]
                    dis = math.sqrt(divec[0] ** 2 + divec[1] ** 2 + divec[2] ** 2)
                    r = random.uniform(-1, 1)
                    dis *= (1 + r * s)
                    if dis == 0:
                        line = str(i + 1) + "\t" + str(j + 1) + "\t" + str(dis) + "\n"
                    else:
                        line = str(i + 1) + "\t" + str(j + 1) + "\t" + str(10000.0 / dis) + "\n"

                    oo_file.write(line)
            oo_file.close()
            # inter_chromosome
            os.makedirs(f"{structure}/{k}/s{s}/inter/", exist_ok=True)
            oo_file = open(f"{structure}/{k}/s{s}/inter/chr1_chr2_{structure}_if_s{s}_{k}.txt", "w")
            for i in range(bin_num1):
                for j in range(bin_num2):
                    divec = positions1[i, :] - positions2[j, :]

                    dis = math.sqrt(divec[0] ** 2 + divec[1] ** 2 + divec[2] ** 2)
                    r = random.uniform(-1, 1)
                    dis *= (1 + r * s)
                    if dis == 0:
                        line = str(i + 1) + "\t" + str(j + 1) + "\t" + str(dis) + "\n"
                    else:
                        line = str(i + 1) + "\t" + str(j + 1) + "\t" + str(10000.0 / dis) + "\n"

                    oo_file.write(line)
            oo_file.close()


### replication fork ###
replication_fork_positions1 = np.zeros((100, 3), np.float32)
replication_fork_positions2 = np.zeros((100, 3), np.float32)
r = 10

for i in range(0, 26):  # Chr 1
    x = i
    y = math.sqrt(350 * (1 - (i ** 2 / 625.00)))
    z = 0

    i1 = i
    i2 = 50 - i1
    replication_fork_positions1[i1, 0] = x
    replication_fork_positions1[i1, 1] = y
    replication_fork_positions1[i2, 0] = x
    replication_fork_positions1[i2, 1] = -y
    replication_fork_positions1[i1, 2] = z
    replication_fork_positions1[i2, 2] = z

    replication_fork_positions2[i1, 0] = x
    replication_fork_positions2[i1, 1] = y
    replication_fork_positions2[i2, 0] = x
    replication_fork_positions2[i2, 1] = -y
    replication_fork_positions2[i1, 2] = z
    replication_fork_positions2[i2, 2] = z

for i in range(51, 76): # Chr 2
    r = 25
    x = 50 - i
    y = math.sqrt(350 * (1 - (x ** 2 / 625.00)))
    angle = 90 / 25 * (-x)
    z = r - r * math.cos(angle * math.pi / 180)

    i1 = i
    i2 = 100 - (i - 50)
    replication_fork_positions1[i1, 0] = x
    replication_fork_positions1[i1, 1] = -y
    replication_fork_positions1[i2, 0] = x
    replication_fork_positions1[i2, 1] = y
    replication_fork_positions1[i1, 2] = z
    replication_fork_positions1[i2, 2] = z

    replication_fork_positions2[i1, 0] = x
    replication_fork_positions2[i1, 1] = -y
    replication_fork_positions2[i2, 0] = x
    replication_fork_positions2[i2, 1] = y
    replication_fork_positions2[i1, 2] = -z
    replication_fork_positions2[i2, 2] = -z

write_pdb_double_chr(replication_fork_positions1, replication_fork_positions2, "fork", 0)
add_random_noise2(replication_fork_positions1, replication_fork_positions2, "fork")


### double helix ###
double_helix_positions1 = np.zeros((100, 3), np.float32)
double_helix_positions2 = np.zeros((100, 3), np.float32)
r = 20

z = 0
for i in range(0, 100):
    angle = 360.0 / 20 * i
    x = r * math.cos(angle * math.pi / 180)
    y = r * math.sin(angle * math.pi / 180)
    z += 1
    double_helix_positions1[i, 0] = x
    double_helix_positions1[i, 1] = y
    double_helix_positions1[i, 2] = z

    double_helix_positions2[i, 0] = x
    double_helix_positions2[i, 1] = y
    double_helix_positions2[i, 2] = z + 10

write_pdb_double_chr(double_helix_positions1, double_helix_positions2, "double_helix", 1)
add_random_noise2(double_helix_positions1, double_helix_positions2, "double_helix")


### double spherical helix ###
double_spherical_helix_positions1 = np.zeros((100, 3), np.float32)
double_spherical_helix_positions2 = np.zeros((100, 3), np.float32)

for i in range(0, 100):
    angle = 360.0 / 10 * i

    if i < 50:
        x = i * math.cos(angle * math.pi / 180)
        y = i * math.sin(angle * math.pi / 180)
        z = i
    else:
        x = (100.0 - i) * math.cos(angle * math.pi / 180)
        y = (100.0 - i) * math.sin(angle * math.pi / 180)
        z = i

    double_spherical_helix_positions1[i, 0] = x  
    double_spherical_helix_positions1[i, 1] = y
    double_spherical_helix_positions1[i, 2] = z

    double_spherical_helix_positions2[i, 0] = x + 60
    double_spherical_helix_positions2[i, 1] = y
    double_spherical_helix_positions2[i, 2] = z + 50

write_pdb_double_chr(double_spherical_helix_positions1, double_spherical_helix_positions2, "double_spherical_helix", 1)
add_random_noise2(double_spherical_helix_positions1, double_spherical_helix_positions2, "double_spherical_helix")
