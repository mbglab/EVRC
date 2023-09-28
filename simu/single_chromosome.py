import os
import math
import random
import numpy as np


def write_pdb_single_chr(positions, structure, chrtype):
    '''
    Generate simulating structure of single chromosome.
    '''
    col1 = "ATOM"
    col3 = "BIN  C" + "1"
    col4 = 1
    col9 = str(0) + " " + str(0)
    bin_num = len(positions)
    atom_num = 0

    os.mkdir(structure)
    oo_file = open(f"{structure}/{structure}.pdb", "w")
    for i in range(1, bin_num + 1):
        col2 = str(atom_num + i)
        col5 = str(i)
        col6 = "%.2f" % positions[i - 1][0]
        col7 = "%.2f" % positions[i - 1][1]
        col8 = "%.2f" % positions[i - 1][2]
        col2 = " " * (5 - len(col2)) + col2
        col5 = " " * (4 - len(col5)) + col5
        col6 = " " * (8 - len(col6)) + col6
        col7 = " " * (8 - len(col7)) + col7
        col8 = " " * (8 - len(col8)) + col8

        col = (col1, col2, col3, col4, col5, col6, col7, col8, col9)
        line = "%s  %s  %s %s%s   %s%s%s  %s\n" % col
        oo_file.write(line)

    col1 = "CONECT"
    for i in range(1, bin_num + 1):
        col2 = str(atom_num + i)
        j = atom_num + i + 1
        if j > atom_num + bin_num:
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

def add_random_noise1(positions, structure):
    '''
    Generate the interaction frequency (IF) matrix of single chromosome simulating structure
    with different levels of random noises.
    '''
    for k in range(1, 11): # Repeat 10 times
        bin_num = len(positions)

        for s in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]: # different noise level

            os.makedirs(structure + "/" + str(k) + "/s" + str(s), exist_ok = True)
            oo_file = open(structure + "/" +str(k) + "/s" + str(s) + f"/chr1_{structure}_if_s" + str(s) + ".txt", "w")

            for i in range(bin_num):
                for j in range(i):
                    divec = positions[i, :] - positions[j, :]
                    dis = math.sqrt(divec[0] ** 2 + divec[1] ** 2 + divec[2] ** 2)
                    r = random.uniform(-1, 1)
                    dis *= (1 + r * s)
                    line = str(i + 1) + "\t" + str(j + 1) + "\t" + str(10000.0 / dis) + "\n"
                    oo_file.write(line)

            oo_file.close()


### circle ###
circle_r = 100
circle_positions = np.zeros ((200, 3), np.float32)

for i in range(0, 200):
    angle = 360.0 / 200 * i
    x = circle_r * math.cos (angle * math.pi / 180)
    y = circle_r * math.sin (angle * math.pi / 180)
    z = 0
    circle_positions [i, 0] = x
    circle_positions [i, 1] = y
    circle_positions [i, 2] = z

write_pdb_single_chr(circle_positions, "circle", 0)
add_random_noise1(circle_positions, "circle")


### spiral ###
spiral_r = 100
spiral_positions = np.zeros ((200, 3), np.float32)

z = 0
for i in range(0,200):
    angle = 360.0 / 50 * i
    x = spiral_r * math.cos (angle * math.pi / 180)
    y = spiral_r * math.sin (angle * math.pi / 180)
    z += 1
    spiral_positions [i, 0] = x
    spiral_positions [i, 1] = y
    spiral_positions [i, 2] = z

write_pdb_single_chr(spiral_positions, "spiral", 1)
add_random_noise1(spiral_positions, "spiral")


### circular_spiral ###
circular_spiral_r = 10
circular_spiral_r1 = 10
circular_spiral_positions = np.zeros((200, 3), np.float32)

z = 0
for i in range(0, 200):
    angle = 360.0 / 25 * (i % 25)
    z = circular_spiral_r * math.cos(angle * math.pi / 180)

    angle1 = 360.0 / 200 * i
    rr = circular_spiral_r1 + circular_spiral_r - circular_spiral_r * math.sin(angle * math.pi / 180)
    x = rr * math.cos(angle1 * math.pi / 180)
    y = rr * math.sin(angle1 * math.pi / 180)
    circular_spiral_positions[i, 0] = x
    circular_spiral_positions[i, 1] = y
    circular_spiral_positions[i, 2] = z

write_pdb_single_chr(circular_spiral_positions, "circular_spiral", 0)
add_random_noise1(circular_spiral_positions, "circular_spiral")



