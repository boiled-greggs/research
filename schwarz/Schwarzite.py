#!/usr/bin/env python

# Toy graphene model

# Copyright under GNU General Public License 2010, 2012, 2016
# by Sinisa Coh and David Vanderbilt (see gpl-pythtb.txt)


from __future__ import print_function
from pythtb import * # import TB model class
import numpy as np
import pylab as plt

# define lattice vectors
lat=[[7.77, 0.00, 0.00], [0.00, 7.77, 0.00],[0.00, 0.00,7.77]]
# define coordinates of orbitals
orb=[[0.579999983, 0.180000007, 0.180000007],\
    [0.420000017, 0.819999993, 0.180000007],\
    [0.420000017, 0.180000007, 0.819999993],\
    [0.579999983, 0.819999993, 0.819999993],\
    [0.180000007, 0.579999983, 0.180000007],\
    [0.180000007, 0.420000017, 0.819999993],\
    [0.819999993, 0.420000017, 0.180000007],\
    [0.819999993, 0.579999983, 0.819999993],\
    [0.180000007, 0.180000007, 0.579999983],\
    [0.819999993, 0.180000007, 0.420000017],\
    [0.180000007, 0.819999993, 0.420000017],\
    [0.819999993, 0.819999993, 0.579999983],\
    [0.079999983, 0.680000007, 0.680000007],\
    [0.920000017, 0.319999993, 0.680000007],\
    [0.920000017, 0.680000007, 0.319999993],\
    [0.079999983, 0.319999993, 0.319999993],\
    [0.680000007, 0.079999983, 0.680000007],\
    [0.680000007, 0.920000017, 0.319999993],\
    [0.319999993, 0.920000017, 0.680000007],\
    [0.319999993, 0.079999983, 0.319999993],\
    [0.680000007, 0.680000007, 0.079999983],\
    [0.319999993, 0.680000007, 0.920000017],\
    [0.680000007, 0.319999993, 0.920000017],\
    [0.319999993, 0.319999993, 0.079999983],\
    [0.660000026, 0.100000001, 0.340000004],\
    [0.339999974, 0.899999976, 0.340000004],\
    [0.339999974, 0.100000001, 0.659999967],\
    [0.660000026, 0.899999976, 0.659999967],\
    [0.340000004, 0.660000026, 0.100000001],\
    [0.340000004, 0.339999974, 0.899999976],\
    [0.659999967, 0.339999974, 0.100000001],\
    [0.659999967, 0.660000026, 0.899999976],\
    [0.100000001, 0.340000004, 0.660000026],\
    [0.899999976, 0.340000004, 0.339999974],\
    [0.100000001, 0.659999967, 0.339999974],\
    [0.899999976, 0.659999967, 0.660000026],\
    [0.160000026, 0.600000024, 0.840000033],\
    [0.839999974, 0.400000006, 0.840000033],\
    [0.839999974, 0.600000024, 0.159999996],\
    [0.160000026, 0.400000006, 0.159999996],\
    [0.840000033, 0.160000026, 0.600000024],\
    [0.840000033, 0.839999974, 0.400000006],\
    [0.159999996, 0.839999974, 0.600000024],\
    [0.159999996, 0.160000026, 0.400000006],\
    [0.600000024, 0.840000033, 0.160000026],\
    [0.400000006, 0.840000033, 0.839999974],\
    [0.600000024, 0.159999996, 0.839999974],\
    [0.400000006, 0.159999996, 0.160000026],]

# make two dimensional tight-binding graphene model
my_model = tb_model(3,3,lat,orb)

# set model parameters
delta=0.0
t1n=-1
t2n=-0.35
t3n=-0.225

# set on-site energies
my_model.set_onsite([delta, delta, delta, delta, delta, delta, delta, delta, \
        delta, delta, delta, delta, delta, delta, delta, delta, delta, delta, \
        delta, delta, delta, delta, delta, delta, delta, delta, delta, delta, \
        delta, delta, delta, delta, delta, delta, delta, delta, delta, delta, \
        delta, delta, delta, delta, delta, delta, delta, delta, delta, delta])
# set hoppings (one for each connected pair of orbitals)
# (amplitude, i, j, [lattice vector to cell containing j])
# print('Check! ' + str(orb[0][0])+ ' ' + str(orb[0][1])+ ' ' + str(orb[0][2]))
# print('Check! ' + str(orb[47][0])+ ' ' + str(orb[47][1])+ ' ' + str(orb[47][2]))
for x in range(0,47):
    for y in range(x+1,47):
        dn_x=(orb[x][0] - orb[y][0])
        dn_y=(orb[x][1] - orb[y][1])
        dn_z=(orb[x][2] - orb[y][2])
        dist = np.sqrt( dn_x**2 + dn_y**2 + dn_z**2)
        if (5.51*dist < 1.3000000):
            my_model.set_hop(t1n, x, y, [0, 0, 0])
        if (5.51*dist > 1.300000) and (5.51*dist < 2.600000):
            my_model.set_hop(t3n, x, y, [0, 0, 0])
        if (orb[x][0] > 0.9) and (orb[y][0] > 0.5):
            my_model.set_hop(t1n, x, y, [1, 0, 0])
        elif (orb[x][1] > 0.8) and (orb[y][1] > 0.5):
            my_model.set_hop(t1n, x, y, [0, 1, 0])
        elif (orb[x][2] > 0.8) and (orb[y][2] > 0.5):
            my_model.set_hop(t1n, x, y, [0, 0, 1])
"""
my_model.set_hop(t, 1, 0, [ 1, 0, 0])
my_model.set_hop(t, 1, 0, [ 0, 1, 0])
my_model.set_hop(t, 0, 1, [ 1, 0, 0])
my_model.set_hop(t, 0, 1, [ 0, 1, 0])
"""
# print tight-binding model
my_model.display()

# generate list of k-points following a segmented path in the BZ
# list of nodes (high-symmetry points) that will be connected
# DON'T KNOW HOW TO DO THIS
# path=[[0., 0., 0.0], [2./15., 1./15., 0.0], [0.0, 0.10, 0.0], [0., 0., 0.0]]
# label=(r'$\Gamma $',r'$K$', r'$M$', r'$\Gamma $')
path=[[0.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.5, 0.5, 0.5], [0.5, 0.5, 0.0], [0.0, 0.0, 0.0],[0.5, 0.5, 0.5]]
label=(r'$\Gamma $', r'$X$', r'$R$', r'$M$', r'$\Gamma $', r'$R$')
# total number of interpolated k-points along the path
nk=400

# call function k_path to construct the actual path
(k_vec,k_dist,k_node)=my_model.k_path(path,nk)
# inputs:
#   path, nk: see above
#   my_model: the pythtb model
# outputs:
#   k_vec: list of interpolated k-points
#   k_dist: horizontal axis position of each k-point in the list
#   k_node: horizontal axis position of each original node

print('---------------------------------------')
print('starting calculation')
print('---------------------------------------')
print('Calculating bands...')

# obtain eigenvalues to be plotted
evals=my_model.solve_all(k_vec)

# figure for bandstructure

fig, ax = plt.subplots()
# specify horizontal axis details
# set range of horizontal axis
# ax.set_xlim([0,k_node[-1]])
# put tickmarks and labels at node positions
ax.set_xticks(k_node)
ax.set_xticklabels(label)
# add vertical lines at node positions
for n in range(len(k_node)):
  ax.axvline(x=k_node[n],linewidth=0.5, color='k')
# put title
ax.set_ylim(-6, 6)
ax.set_title("Band of P_8-0 Schwarzite")
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")

# plot first and second band
for xband in range(0,47):
    ax.plot(k_dist,evals[xband])
    # ax.plot(k_dist,evals[1])

# make an PDF figure of a plot
# fig.tight_layout()
fig.savefig("P_8-0.pdf")

print('Done.\n')

# def edgeCheck(orbX,orbY,orbZ,edge_pos):
#     thr=0.85
#     rayon=np.sqrt(orbX**2+orbY**2+orbZ**2)
#     rayonXY=np.sqrt(orbX**2+orbY**2)
#     rayonYZ=np.sqrt(orbY**2+orbZ**2)
#     rayonXZ=np.sqrt(orbZ**2+orbX**2)
#     if rayon > 0.93:
#         if (orbX > 0.90 and rayonYZ < thr):
#             edge_pos=[1, 0, 0]
#         if (orbY > 0.90 and rayonXZ < thr):
#             edge_pos=[0, 1, 0]
#         if (orbZ > 0.90 and rayonXY < thr):
#             edge_pos=[0, 0, 1]
#     return
