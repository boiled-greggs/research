# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[10.1660003662, 0.0000000000, 0.0000000000], \
	[0.0000000000, 10.1660003662, 0.0000000000], \
	[0.0000000000, 0.0000000000, 2.4619998932]]

orb = [[0.835424662, 0.500000000, 0.500072479], \
	[0.806210279, 0.636727333, 0.500072479], \
	[0.771290660, 0.697123468, 0.000000000], \
	[0.667417228, 0.790570199, 0.000000000], \
	[0.603676736, 0.818997741, 0.500072479], \
	[0.464588583, 0.833555698, 0.500072479], \
	[0.396323293, 0.818997741, 0.000000000], \
	[0.275334358, 0.749060154, 0.000000000], \
	[0.228709340, 0.697123468, 0.500072479], \
	[0.171854377, 0.569445670, 0.500072479], \
	[0.164575338, 0.500000000, 0.000000000], \
	[0.193789750, 0.363272667, 0.000000000], \
	[0.228709340, 0.302876532, 0.500072479], \
	[0.332582772, 0.209429801, 0.500072479], \
	[0.396323293, 0.181002289, 0.000000000], \
	[0.535411417, 0.166444272, 0.000000000], \
	[0.603676736, 0.181002289, 0.500072479], \
	[0.724665642, 0.250939816, 0.500072479], \
	[0.771290660, 0.302876532, 0.000000000], \
	[0.828145623, 0.430554330, 0.000000000]]

my_model = tb_model(1, 3, lat, orb, per=[2])

# set model parameters
delta = 0.0
t1 = -1.8
t2 = 0.
t3 = 0.

my_model.set_onsite([delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta])

# set hopping parameters for connected orbitals
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(t2, 0, 0, [0, 0, 1])
my_model.set_hop(t3, 0, 1, [0, 0, 1])
my_model.set_hop(t2, 0, 2, [0, 0, 1])
my_model.set_hop(t3, 0, 3, [0, 0, 1])
my_model.set_hop(t3, 0, 17, [0, 0, 1])
my_model.set_hop(t2, 0, 18, [0, 0, 1])
my_model.set_hop(t1, 0, 19, [0, 0, 1])
my_model.set_hop(t1, 0, 1, [0, 0, 0])
my_model.set_hop(t2, 0, 2, [0, 0, 0])
my_model.set_hop(t3, 0, 3, [0, 0, 0])
my_model.set_hop(t3, 0, 17, [0, 0, 0])
my_model.set_hop(t2, 0, 18, [0, 0, 0])
my_model.set_hop(t1, 0, 19, [0, 0, 0])
my_model.set_hop(t3, 1, 0, [0, 0, 1])
my_model.set_hop(t2, 1, 1, [0, 0, 1])
my_model.set_hop(t1, 1, 2, [0, 0, 1])
my_model.set_hop(t2, 1, 3, [0, 0, 1])
my_model.set_hop(t3, 1, 4, [0, 0, 1])
my_model.set_hop(t3, 1, 18, [0, 0, 1])
my_model.set_hop(t2, 1, 19, [0, 0, 1])
my_model.set_hop(t1, 1, 2, [0, 0, 0])
my_model.set_hop(t2, 1, 3, [0, 0, 0])
my_model.set_hop(t3, 1, 4, [0, 0, 0])
my_model.set_hop(t3, 1, 18, [0, 0, 0])
my_model.set_hop(t2, 1, 19, [0, 0, 0])
my_model.set_hop(t1, 2, 3, [0, 0, 0])
my_model.set_hop(t2, 2, 4, [0, 0, 0])
my_model.set_hop(t3, 2, 5, [0, 0, 0])
my_model.set_hop(t3, 2, 19, [0, 0, 0])
my_model.set_hop(t1, 3, 4, [0, 0, 0])
my_model.set_hop(t2, 3, 5, [0, 0, 0])
my_model.set_hop(t3, 3, 6, [0, 0, 0])
my_model.set_hop(t3, 4, 1, [0, 0, 1])
my_model.set_hop(t2, 4, 2, [0, 0, 1])
my_model.set_hop(t1, 4, 3, [0, 0, 1])
my_model.set_hop(t2, 4, 4, [0, 0, 1])
my_model.set_hop(t3, 4, 5, [0, 0, 1])
my_model.set_hop(t2, 4, 6, [0, 0, 1])
my_model.set_hop(t3, 4, 7, [0, 0, 1])
my_model.set_hop(t1, 4, 5, [0, 0, 0])
my_model.set_hop(t2, 4, 6, [0, 0, 0])
my_model.set_hop(t3, 4, 7, [0, 0, 0])
my_model.set_hop(t3, 5, 2, [0, 0, 1])
my_model.set_hop(t2, 5, 3, [0, 0, 1])
my_model.set_hop(t3, 5, 4, [0, 0, 1])
my_model.set_hop(t2, 5, 5, [0, 0, 1])
my_model.set_hop(t1, 5, 6, [0, 0, 1])
my_model.set_hop(t2, 5, 7, [0, 0, 1])
my_model.set_hop(t3, 5, 8, [0, 0, 1])
my_model.set_hop(t1, 5, 6, [0, 0, 0])
my_model.set_hop(t2, 5, 7, [0, 0, 0])
my_model.set_hop(t3, 5, 8, [0, 0, 0])
my_model.set_hop(t1, 6, 7, [0, 0, 0])
my_model.set_hop(t2, 6, 8, [0, 0, 0])
my_model.set_hop(t3, 6, 9, [0, 0, 0])
my_model.set_hop(t1, 7, 8, [0, 0, 0])
my_model.set_hop(t2, 7, 9, [0, 0, 0])
my_model.set_hop(t3, 7, 10, [0, 0, 0])
my_model.set_hop(t3, 8, 5, [0, 0, 1])
my_model.set_hop(t2, 8, 6, [0, 0, 1])
my_model.set_hop(t1, 8, 7, [0, 0, 1])
my_model.set_hop(t2, 8, 8, [0, 0, 1])
my_model.set_hop(t3, 8, 9, [0, 0, 1])
my_model.set_hop(t2, 8, 10, [0, 0, 1])
my_model.set_hop(t3, 8, 11, [0, 0, 1])
my_model.set_hop(t1, 8, 9, [0, 0, 0])
my_model.set_hop(t2, 8, 10, [0, 0, 0])
my_model.set_hop(t3, 8, 11, [0, 0, 0])
my_model.set_hop(t3, 9, 6, [0, 0, 1])
my_model.set_hop(t2, 9, 7, [0, 0, 1])
my_model.set_hop(t3, 9, 8, [0, 0, 1])
my_model.set_hop(t2, 9, 9, [0, 0, 1])
my_model.set_hop(t1, 9, 10, [0, 0, 1])
my_model.set_hop(t2, 9, 11, [0, 0, 1])
my_model.set_hop(t3, 9, 12, [0, 0, 1])
my_model.set_hop(t1, 9, 10, [0, 0, 0])
my_model.set_hop(t2, 9, 11, [0, 0, 0])
my_model.set_hop(t3, 9, 12, [0, 0, 0])
my_model.set_hop(t1, 10, 11, [0, 0, 0])
my_model.set_hop(t2, 10, 12, [0, 0, 0])
my_model.set_hop(t3, 10, 13, [0, 0, 0])
my_model.set_hop(t1, 11, 12, [0, 0, 0])
my_model.set_hop(t2, 11, 13, [0, 0, 0])
my_model.set_hop(t3, 11, 14, [0, 0, 0])
my_model.set_hop(t3, 12, 9, [0, 0, 1])
my_model.set_hop(t2, 12, 10, [0, 0, 1])
my_model.set_hop(t1, 12, 11, [0, 0, 1])
my_model.set_hop(t2, 12, 12, [0, 0, 1])
my_model.set_hop(t3, 12, 13, [0, 0, 1])
my_model.set_hop(t2, 12, 14, [0, 0, 1])
my_model.set_hop(t3, 12, 15, [0, 0, 1])
my_model.set_hop(t1, 12, 13, [0, 0, 0])
my_model.set_hop(t2, 12, 14, [0, 0, 0])
my_model.set_hop(t3, 12, 15, [0, 0, 0])
my_model.set_hop(t3, 13, 10, [0, 0, 1])
my_model.set_hop(t2, 13, 11, [0, 0, 1])
my_model.set_hop(t3, 13, 12, [0, 0, 1])
my_model.set_hop(t2, 13, 13, [0, 0, 1])
my_model.set_hop(t1, 13, 14, [0, 0, 1])
my_model.set_hop(t2, 13, 15, [0, 0, 1])
my_model.set_hop(t3, 13, 16, [0, 0, 1])
my_model.set_hop(t1, 13, 14, [0, 0, 0])
my_model.set_hop(t2, 13, 15, [0, 0, 0])
my_model.set_hop(t3, 13, 16, [0, 0, 0])
my_model.set_hop(t1, 14, 15, [0, 0, 0])
my_model.set_hop(t2, 14, 16, [0, 0, 0])
my_model.set_hop(t3, 14, 17, [0, 0, 0])
my_model.set_hop(t1, 15, 16, [0, 0, 0])
my_model.set_hop(t2, 15, 17, [0, 0, 0])
my_model.set_hop(t3, 15, 18, [0, 0, 0])
my_model.set_hop(t3, 16, 13, [0, 0, 1])
my_model.set_hop(t2, 16, 14, [0, 0, 1])
my_model.set_hop(t1, 16, 15, [0, 0, 1])
my_model.set_hop(t2, 16, 16, [0, 0, 1])
my_model.set_hop(t3, 16, 17, [0, 0, 1])
my_model.set_hop(t2, 16, 18, [0, 0, 1])
my_model.set_hop(t3, 16, 19, [0, 0, 1])
my_model.set_hop(t1, 16, 17, [0, 0, 0])
my_model.set_hop(t2, 16, 18, [0, 0, 0])
my_model.set_hop(t3, 16, 19, [0, 0, 0])
my_model.set_hop(t3, 17, 0, [0, 0, 1])
my_model.set_hop(t3, 17, 14, [0, 0, 1])
my_model.set_hop(t2, 17, 15, [0, 0, 1])
my_model.set_hop(t3, 17, 16, [0, 0, 1])
my_model.set_hop(t2, 17, 17, [0, 0, 1])
my_model.set_hop(t1, 17, 18, [0, 0, 1])
my_model.set_hop(t2, 17, 19, [0, 0, 1])
my_model.set_hop(t1, 17, 18, [0, 0, 0])
my_model.set_hop(t2, 17, 19, [0, 0, 0])
my_model.set_hop(t1, 18, 19, [0, 0, 0])

my_model.display()

path = 'fullc'
label = (r'$\pi/a$',r'$\Gamma $',r'$\pi/a$')
nk = 800

(k_vec, k_dist, k_node) = my_model.k_path(path, nk)

print('-'*20)
print('starting calculation')
print('-'*20)
print('Calculating bands . . .')

evals = my_model.solve_all(k_vec)

fig, ax = plt.subplots(figsize=(3,4))

# ax.set_xlim([0, k_node[-1]])
ax.set_xticks(k_node)
ax.set_xticklabels(label)

for n in range(len(k_node)):
	ax.axvline(x=k_node[n], linewidth=0.3, color='black')

ax.set_title('5-5-tube band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(20):
	ax.plot(k_dist, evals[xband], linewidth=.5, color='black')

# fig.tight_layout()
fig.savefig('5-5-tube.eps')

"""
kmesh=20
kpts=[]
for i in range(kmesh):
    for j in range(kmesh):
        kpts.append([float(i)/float(kmesh)]) # ,float(j)/float(kmesh)])
# solve the model on this mesh
"""
# k_vecs = my_model.k_uniform_mesh([10,20,30])
kmesh = 100
kpts = []
for i in range(kmesh):
    kpts.append(float(i)/float(kmesh))
evals=my_model.solve_all(kpts)
# flatten completely the matrix
# evals=evals.flatten()

# plotting DOS
print('Plotting DOS...')

# now plot density of states
fig, ax = plt.subplots()
# ax.hist(evals,100,range=(-6.,6.))
print(evals)
print(evals[0][50])
ax.set_ylim(0.0,100.0)
# put title
ax.set_title("5,5 tube density of states")
ax.set_xlabel("Band energy")
ax.set_ylabel("Number of states")
# make an PDF figure of a plot
fig.tight_layout()
fig.savefig("5-5-tube_dos.eps")

print('Done.')
