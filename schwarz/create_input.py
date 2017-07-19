# Should take as input some atomic modeling export, such as a .xyz file, and create at least part of a PythTB
# file by first defining the lattice parameters and then taking the data and putting it in crystal
# coordinates, if it isn't already.

import numpy as np
import copy

def get_lat(f):
    lat = []
    for x in range(3):
        line = f.readline().split()
        lat.append(line)
    return lat

def add_lat(lat, out):
    out.write("lat = [")
    for i in range(3):
        if i < 2:
            out.write("["+lat[i][0]+", "+lat[i][1]+", "+lat[i][2]+"], \\"+"\n\t")
        else:
            out.write("["+lat[i][0]+", "+lat[i][1]+", "+lat[i][2]+"]")
    out.write("]\n\n")

def get_orb(f):
    orb = []
    for line in f:
        line = line.split()
        orb.append(line)
    return orb

def add_orb(orb, out, length): 
    out.write("orb = [")
    for i in range(length):
        if i < length - 1:
            out.write("["+orb[i][0]+", "+orb[i][1]+", "+orb[i][2]+"], \\"+"\n\t")
        else:
            out.write("["+orb[i][0]+", "+orb[i][1]+", "+orb[i][2]+"]")
    out.write("]\n\n")

def set_onsite(length):
    out.write("my_model.set_onsite([")
    for i in range(length):
        if i < length - 1 and (i+1)%5 != 0:
            out.write("delta, ")
        elif i < length - 1 and (i+1)%5 == 0:
            out.write("delta, \\"+"\n\t")
        else:
            out.write("delta])\n\n")

def set_hops(lat, orb, length, maxi):
    out.write("# set hopping parameters for connected orbitals\n")
    out.write("# (amplitude, i, j, [lattice vector to cell containing j])\n")
    m = maxi - .055
    for i in range(0, length):
        if orb[i][0] > m:
            # shift structure by positive lat_vec_a and check
            orb_cp = copy.deepcopy(orb)
            for k in range(length):
                orb_cp[k][0] += 1.0000
            for j in range(length):
                orb_i_real = []
                orb_j_real = []
                for y in range(3):
                    orb_i_y = orb[i][0]*lat[0][y] + orb[i][1]*lat[1][y] + orb[i][2]*lat[2][y]
                    orb_i_real.append(orb_i_y)
                    orb_j_y = orb_cp[j][0]*lat[0][y] + orb_cp[j][1]*lat[1][y] + orb_cp[j][2]*lat[2][y]
                    orb_j_real.append(orb_j_y)
                dist = np.sqrt((orb_i_real[0]-orb_j_real[0])**2 + (orb_i_real[1]-orb_j_real[1])**2 + (orb_i_real[2]-orb_j_real[2])**2)
                if dist < 1.6 and dist > 0.01:
                    out.write("my_model.set_hop(t1, %d, %d, [1, 0, 0])\n" %(i, j))
                elif dist > 1.7 and dist < 2.7:
                    out.write("my_model.set_hop(t2, %d, %d, [1, 0, 0])\n" %(i, j))
                elif dist > 2.7 and dist < 3.8:
                    out.write("my_model.set_hop(t3, %d, %d, [1, 0, 0])\n" %(i, j))
        elif orb[i][1] > m:
            # shift structure by positive lat_vec_b and check
            orb_cp = copy.deepcopy(orb)
            for k in range(length):
                orb_cp[k][1] += 1.0000
            for j in range(length):
                orb_i_real = []
                orb_j_real = []
                for y in range(3):
                    orb_i_y = orb[i][0]*lat[0][y] + orb[i][1]*lat[1][y] + orb[i][2]*lat[2][y]
                    orb_i_real.append(orb_i_y)
                    orb_j_y = orb_cp[j][0]*lat[0][y] + orb_cp[j][1]*lat[1][y] + orb_cp[j][2]*lat[2][y]
                    orb_j_real.append(orb_j_y)
                dist = np.sqrt((orb_i_real[0]-orb_j_real[0])**2 + (orb_i_real[1]-orb_j_real[1])**2 + (orb_i_real[2]-orb_j_real[2])**2)
                if dist < 1.6 and dist > 0.01:
                    out.write("my_model.set_hop(t1, %d, %d, [0, 1, 0])\n" %(i, j))
                elif dist > 1.7 and dist < 2.7:
                    out.write("my_model.set_hop(t2, %d, %d, [0, 1, 0])\n" %(i, j))
                elif dist > 2.7 and dist < 3.8:
                    out.write("my_model.set_hop(t3, %d, %d, [0, 1, 0])\n" %(i, j))
        if orb[i][2] > m:
            # shift structure by positive lat_vec_c and check
            orb_cp = copy.deepcopy(orb)
            for k in range(length):
                orb_cp[k][2] += 1.0000
            for j in range(length):
                orb_i_real = []
                orb_j_real = []
                for y in range(3):
                    orb_i_y = orb[i][0]*lat[0][y] + orb[i][1]*lat[1][y] + orb[i][2]*lat[2][y]
                    orb_i_real.append(orb_i_y)
                    orb_j_y = orb_cp[j][0]*lat[0][y] + orb_cp[j][1]*lat[1][y] + orb_cp[j][2]*lat[2][y]
                    orb_j_real.append(orb_j_y)
                dist = np.sqrt((orb_i_real[0]-orb_j_real[0])**2 + (orb_i_real[1]-orb_j_real[1])**2 + (orb_i_real[2]-orb_j_real[2])**2)
                if dist < 1.6 and dist > 0.01:
                    out.write("my_model.set_hop(t1, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 1.7 and dist < 2.7:
                    out.write("my_model.set_hop(t2, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 2.7 and dist < 3.8:
                    out.write("my_model.set_hop(t3, %d, %d, [0, 0, 1])\n" %(i, j))
        if orb[i][0] > m and orb[i][1] > m:
            #shift by both lattice vectors
            orb_cp = copy.deepcopy(orb)
            for k in range(length):
                orb_cp[k][0] += 1.0000
                orb_cp[k][1] += 1.0000
            for j in range(length):
                orb_i_real = []
                orb_j_real = []
                for y in range(3):
                    orb_i_y = orb[i][0]*lat[0][y] + orb[i][1]*lat[1][y] + orb[i][2]*lat[2][y]
                    orb_i_real.append(orb_i_y)
                    orb_j_y = orb_cp[j][0]*lat[0][y] + orb_cp[j][1]*lat[1][y] + orb_cp[j][2]*lat[2][y]
                    orb_j_real.append(orb_j_y)
                dist = np.sqrt((orb_i_real[0]-orb_j_real[0])**2 + (orb_i_real[1]-orb_j_real[1])**2 + (orb_i_real[2]-orb_j_real[2])**2)
                if dist < 1.6 and dist > 0.01:
                    out.write("my_model.set_hop(t1, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 1.7 and dist < 2.7:
                    out.write("my_model.set_hop(t2, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 2.7 and dist < 3.8:
                    out.write("my_model.set_hop(t3, %d, %d, [0, 0, 1])\n" %(i, j))
        if orb[i][0] > m and orb[i][2] > m:
            #shift by both lattice vectors
            orb_cp = copy.deepcopy(orb)
            for k in range(length):
                orb_cp[k][0] += 1.0000
                orb_cp[k][2] += 1.0000
            for j in range(length):
                orb_i_real = []
                orb_j_real = []
                for y in range(3):
                    orb_i_y = orb[i][0]*lat[0][y] + orb[i][1]*lat[1][y] + orb[i][2]*lat[2][y]
                    orb_i_real.append(orb_i_y)
                    orb_j_y = orb_cp[j][0]*lat[0][y] + orb_cp[j][1]*lat[1][y] + orb_cp[j][2]*lat[2][y]
                    orb_j_real.append(orb_j_y)
                dist = np.sqrt((orb_i_real[0]-orb_j_real[0])**2 + (orb_i_real[1]-orb_j_real[1])**2 + (orb_i_real[2]-orb_j_real[2])**2)
                if dist < 1.6 and dist > 0.01:
                    out.write("my_model.set_hop(t1, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 1.7 and dist < 2.7:
                    out.write("my_model.set_hop(t2, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 2.7 and dist < 3.8:
                    out.write("my_model.set_hop(t3, %d, %d, [0, 0, 1])\n" %(i, j))
        if orb[i][1] > m and orb[i][2] > m:
            #shift by both lattice vectors
            orb_cp = copy.deepcopy(orb)
            for k in range(length):
                orb_cp[k][1] += 1.0000
                orb_cp[k][2] += 1.0000
            for j in range(length):
                orb_i_real = []
                orb_j_real = []
                for y in range(3):
                    orb_i_y = orb[i][0]*lat[0][y] + orb[i][1]*lat[1][y] + orb[i][2]*lat[2][y]
                    orb_i_real.append(orb_i_y)
                    orb_j_y = orb_cp[j][0]*lat[0][y] + orb_cp[j][1]*lat[1][y] + orb_cp[j][2]*lat[2][y]
                    orb_j_real.append(orb_j_y)
                dist = np.sqrt((orb_i_real[0]-orb_j_real[0])**2 + (orb_i_real[1]-orb_j_real[1])**2 + (orb_i_real[2]-orb_j_real[2])**2)
                if dist < 1.6 and dist > 0.01:
                    out.write("my_model.set_hop(t1, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 1.7 and dist < 2.7:
                    out.write("my_model.set_hop(t2, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 2.7 and dist < 3.8:
                    out.write("my_model.set_hop(t3, %d, %d, [0, 0, 1])\n" %(i, j))
        if orb[i][0] > m and orb[i][1] > m and orb[i][2] > m:
            # shift by all lattice vectors
            orb_cp = copy.deepcopy(orb)
            for k in range(length):
                orb_cp[k][0] += 1.0000
                orb_cp[k][1] += 1.0000
                orb_cp[k][2] += 1.0000
            for j in range(length):
                orb_i_real = []
                orb_j_real = []
                for y in range(3):
                    orb_i_y = orb[i][0]*lat[0][y] + orb[i][1]*lat[1][y] + orb[i][2]*lat[2][y]
                    orb_i_real.append(orb_i_y)
                    orb_j_y = orb_cp[j][0]*lat[0][y] + orb_cp[j][1]*lat[1][y] + orb_cp[j][2]*lat[2][y]
                    orb_j_real.append(orb_j_y)
                dist = np.sqrt((orb_i_real[0]-orb_j_real[0])**2 + (orb_i_real[1]-orb_j_real[1])**2 + (orb_i_real[2]-orb_j_real[2])**2)
                if dist < 1.6 and dist > 0.01:
                    out.write("my_model.set_hop(t1, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 1.7 and dist < 2.7:
                    out.write("my_model.set_hop(t2, %d, %d, [0, 0, 1])\n" %(i, j))
                elif dist > 2.7 and dist < 3.8:
                    out.write("my_model.set_hop(t3, %d, %d, [0, 0, 1])\n" %(i, j))
        
        for j in range(i+1, length):
            orb_i_real = []
            orb_j_real = []
            for y in range(3):
                orb_i_y = orb[i][0]*lat[0][y] + orb[i][1]*lat[1][y] + orb[i][2]*lat[2][y]
                orb_i_real.append(orb_i_y)
                orb_j_y = orb[j][0]*lat[0][y] + orb[j][1]*lat[1][y] + orb[j][2]*lat[2][y]
                orb_j_real.append(orb_j_y)
            dist = np.sqrt((orb_i_real[0]-orb_j_real[0])**2 + (orb_i_real[1]-orb_j_real[1])**2 + (orb_i_real[2]-orb_j_real[2])**2)
            if dist < 1.6 and dist > 0.01:
                out.write("my_model.set_hop(t1, %d, %d, [0, 0, 0])\n" %(i, j))
            elif dist > 1.7 and dist < 2.7:
                out.write("my_model.set_hop(t2, %d, %d, [0, 0, 0])\n" %(i, j))
            elif dist > 2.7 and dist < 3.8:
                out.write("my_model.set_hop(t3, %d, %d, [0, 0, 0])\n" %(i, j))



if __name__ == "__main__":
    structure = input('Enter the .vasp file name: ')
    out = open(structure+"_TB.py", "w")
    f = open(structure+".vasp", "r")
    
    out.write("# Greg Stewart 2017\n\n")
    out.write("from __future__ import print_function\n")
    out.write("from pythtb import *\n")
    out.write("import numpy as np\n")
    out.write("import pylab as plt\n\n")

    f.readline()
    f.readline()

    lat = get_lat(f)

    add_lat(lat, out)
    
    f.readline()
    f.readline()
    f.readline()

    orb = get_orb(f)
    # print('orbitals:')
    # print(orb)
    length = len(orb)

    add_orb(orb, out, length)
    
    out.write("my_model = tb_model(3, 3, lat, orb)\n\n# set model parameters\n")
    out.write("delta = 0.0\nt1 = -2.8\nt2 = -0.09\nt3 = -.3\n\n")
    
    set_onsite(length)
    
    maxi = 0.0
    for i in range(length):
        for j in range(3):
            orb[i][j] = float(orb[i][j])
            if float(orb[i][2]) > maxi:
                maxi = orb[i][j]
    print('maximum orbital coordinate: %f' %maxi)

    for i in range(3):
        for j in range(3):
            lat[i][j] = float(lat[i][j])

    set_hops(lat, orb, length, maxi)
    
    out.write("\n"+"my_model.display()\n\n") # print model

    out.write("path = [ [0., 0., 0.], [-0.5, 0.5, 0.5], [0., .5, 0.], [.25, .25, .25], [0., 0., 0.], [0., 0.5, 0.] ]\n")
    out.write("label = (r'$\Gamma $',r'$H$', r'$N$', r'$P$', r'$\Gamma $', r'$N$')\n")
    out.write("nk = 800\n\n")

    # make path
    out.write("(k_vec, k_dist, k_node) = my_model.k_path(path, nk)\n\n")

    out.write("print('-'*20)\n")
    out.write("print('starting calculation')\n")
    out.write("print('-'*20)\n")
    out.write("print('Calculating bands . . .')\n\n")

    # get eigenvalues
    out.write("evals = my_model.solve_all(k_vec)\n\n")

    # get band structure figure

    out.write("fig, ax = plt.subplots(figsize=(3,4))\n\n")
    out.write("# ax.set_xlim([0, k_node[-1]])\n")
    out.write("ax.set_xticks(k_node)\n")
    out.write("ax.set_xticklabels(label)\n\n")
    out.write("for n in range(len(k_node)):\n\t")
    out.write("ax.axvline(x=k_node[n], linewidth=0.3, color='black')\n\n")
    out.write("ax.set_title('"+structure+" band structure')\n")
    out.write("ax.set_xlabel('Path in k-space')\n")
    out.write("ax.set_ylabel('Band energy')\n\n")
    out.write("for xband in range("+str(length)+"):\n")
    out.write("\tax.plot(k_dist, evals[xband], linewidth=.5, color='black')\n")
    out.write("\n# fig.tight_layout()\n")
    out.write("fig.savefig('"+structure+".eps')\n\n")
    out.write("print('Done.')")

    out.close()

    print('\nPython Tight Binding file has been created.')
