# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[7.7800002098, 0.0000000000, 0.0000000000], \
	[0.0000000000, 7.7800002098, 0.0000000000], \
	[0.0000000000, 0.0000000000, 7.7800002098]]

orb = [[0.310631037, 0.310664266, 0.086759776], \
	[0.810665309, 0.810657680, 0.586753130], \
	[0.689308941, 0.689322472, 0.086759776], \
	[0.189274669, 0.189329118, 0.586753130], \
	[0.688773572, 0.311178386, 0.913226962], \
	[0.188821062, 0.811171770, 0.413233578], \
	[0.311166346, 0.688808382, 0.913226962], \
	[0.811247408, 0.188814968, 0.413233578], \
	[0.086724997, 0.311178386, 0.310921311], \
	[0.586759269, 0.811171770, 0.810914695], \
	[0.086874716, 0.689322472, 0.689065397], \
	[0.586827278, 0.189329118, 0.189072043], \
	[0.913343489, 0.688808382, 0.310921311], \
	[0.413309216, 0.188814968, 0.810914695], \
	[0.913065195, 0.310664266, 0.689065397], \
	[0.413112700, 0.810657680, 0.189072043], \
	[0.311174929, 0.086759776, 0.310921311], \
	[0.811209083, 0.586753130, 0.810914695], \
	[0.689318597, 0.086888306, 0.689065397], \
	[0.189237550, 0.586881697, 0.189072043], \
	[0.310749829, 0.913098454, 0.689065397], \
	[0.810702324, 0.413105071, 0.189072043], \
	[0.688893557, 0.913226962, 0.310921311], \
	[0.188859314, 0.413233578, 0.810914695], \
	[0.310631037, 0.310664266, 0.913226962], \
	[0.810665309, 0.810657680, 0.413233578], \
	[0.689308941, 0.689322472, 0.913226962], \
	[0.189274669, 0.189329118, 0.413233578], \
	[0.311166346, 0.688808382, 0.086759776], \
	[0.811247408, 0.188814968, 0.586881697], \
	[0.688773572, 0.311178386, 0.086759776], \
	[0.188821062, 0.811171770, 0.586881697], \
	[0.311174929, 0.086759776, 0.689065397], \
	[0.811209083, 0.586753130, 0.189072043], \
	[0.689318597, 0.086888306, 0.310921311], \
	[0.189237550, 0.586881697, 0.810914695], \
	[0.688893557, 0.913226962, 0.689065397], \
	[0.188859314, 0.413233578, 0.189072043], \
	[0.310749829, 0.913098454, 0.310921311], \
	[0.810702324, 0.413105071, 0.810914695], \
	[0.086724997, 0.311178386, 0.689065397], \
	[0.586759269, 0.811171770, 0.189072043], \
	[0.086874716, 0.689322472, 0.310921311], \
	[0.586827278, 0.189329118, 0.810914695], \
	[0.913065195, 0.310664266, 0.310921311], \
	[0.413112700, 0.810657680, 0.810914695], \
	[0.913343489, 0.688808382, 0.689065397], \
	[0.413309216, 0.188814968, 0.189072043]]

my_model = tb_model(3, 3, lat, orb)

# set model parameters
delta = 0.0
t1 = -1.0
t2 = -0.30
t3 = -0.20

my_model.set_onsite([delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta])

# set hopping parameters for connected orbitals
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(t2, 0, 8, [0, 0, 0])
my_model.set_hop(t2, 0, 11, [0, 0, 0])
my_model.set_hop(t2, 0, 16, [0, 0, 0])
my_model.set_hop(t2, 0, 19, [0, 0, 0])
my_model.set_hop(t1, 0, 37, [0, 0, 0])
my_model.set_hop(t1, 0, 47, [0, 0, 0])
my_model.set_hop(t2, 1, 9, [0, 0, 0])
my_model.set_hop(t2, 1, 12, [0, 0, 0])
my_model.set_hop(t2, 1, 17, [0, 0, 0])
my_model.set_hop(t2, 1, 22, [0, 0, 0])
my_model.set_hop(t1, 1, 25, [0, 0, 0])
my_model.set_hop(t1, 1, 36, [0, 0, 0])
my_model.set_hop(t1, 1, 46, [0, 0, 0])
my_model.set_hop(t2, 2, 12, [0, 0, 0])
my_model.set_hop(t2, 2, 15, [0, 0, 0])
my_model.set_hop(t2, 2, 21, [0, 0, 0])
my_model.set_hop(t2, 2, 22, [0, 0, 0])
my_model.set_hop(t1, 2, 33, [0, 0, 0])
my_model.set_hop(t1, 2, 41, [0, 0, 0])
my_model.set_hop(t2, 3, 8, [0, 0, 0])
my_model.set_hop(t2, 3, 13, [0, 0, 0])
my_model.set_hop(t2, 3, 16, [0, 0, 0])
my_model.set_hop(t2, 3, 23, [0, 0, 0])
my_model.set_hop(t1, 3, 27, [0, 0, 0])
my_model.set_hop(t1, 3, 32, [0, 0, 0])
my_model.set_hop(t1, 3, 40, [0, 0, 0])
my_model.set_hop(t2, 4, 11, [0, 0, 1])
my_model.set_hop(t2, 4, 21, [0, 0, 1])
my_model.set_hop(t1, 4, 30, [0, 0, 1])
my_model.set_hop(t2, 4, 13, [0, 0, 0])
my_model.set_hop(t2, 4, 14, [0, 0, 0])
my_model.set_hop(t2, 4, 17, [0, 0, 0])
my_model.set_hop(t2, 4, 18, [0, 0, 0])
my_model.set_hop(t1, 4, 39, [0, 0, 0])
my_model.set_hop(t1, 4, 43, [0, 0, 0])
my_model.set_hop(t2, 5, 10, [0, 0, 0])
my_model.set_hop(t2, 5, 15, [0, 0, 0])
my_model.set_hop(t2, 5, 19, [0, 0, 0])
my_model.set_hop(t2, 5, 20, [0, 0, 0])
my_model.set_hop(t1, 5, 31, [0, 0, 0])
my_model.set_hop(t1, 5, 38, [0, 0, 0])
my_model.set_hop(t1, 5, 42, [0, 0, 0])
my_model.set_hop(t2, 6, 15, [0, 0, 1])
my_model.set_hop(t2, 6, 19, [0, 0, 1])
my_model.set_hop(t1, 6, 28, [0, 0, 1])
my_model.set_hop(t2, 6, 9, [0, 0, 0])
my_model.set_hop(t2, 6, 10, [0, 0, 0])
my_model.set_hop(t2, 6, 20, [0, 0, 0])
my_model.set_hop(t2, 6, 23, [0, 0, 0])
my_model.set_hop(t1, 6, 35, [0, 0, 0])
my_model.set_hop(t1, 6, 45, [0, 0, 0])
my_model.set_hop(t2, 7, 11, [0, 0, 0])
my_model.set_hop(t2, 7, 14, [0, 0, 0])
my_model.set_hop(t2, 7, 18, [0, 0, 0])
my_model.set_hop(t2, 7, 21, [0, 0, 0])
my_model.set_hop(t1, 7, 29, [0, 0, 0])
my_model.set_hop(t1, 7, 34, [0, 0, 0])
my_model.set_hop(t1, 7, 44, [0, 0, 0])
my_model.set_hop(t2, 8, 16, [0, 0, 0])
my_model.set_hop(t2, 8, 19, [0, 0, 0])
my_model.set_hop(t1, 8, 27, [0, 0, 0])
my_model.set_hop(t1, 8, 37, [0, 0, 0])
my_model.set_hop(t2, 9, 17, [0, 0, 0])
my_model.set_hop(t2, 9, 20, [0, 0, 0])
my_model.set_hop(t1, 9, 26, [0, 0, 0])
my_model.set_hop(t1, 9, 36, [0, 0, 0])
my_model.set_hop(t1, 9, 45, [0, 0, 0])
my_model.set_hop(t2, 10, 20, [0, 0, 0])
my_model.set_hop(t2, 10, 23, [0, 0, 0])
my_model.set_hop(t1, 10, 31, [0, 0, 0])
my_model.set_hop(t1, 10, 35, [0, 0, 0])
my_model.set_hop(t2, 11, 16, [0, 0, 0])
my_model.set_hop(t2, 11, 21, [0, 0, 0])
my_model.set_hop(t1, 11, 30, [0, 0, 0])
my_model.set_hop(t1, 11, 34, [0, 0, 0])
my_model.set_hop(t1, 11, 47, [0, 0, 0])
my_model.set_hop(t2, 12, 5, [1, 0, 0])
my_model.set_hop(t2, 12, 19, [1, 0, 0])
my_model.set_hop(t1, 12, 42, [1, 0, 0])
my_model.set_hop(t2, 12, 21, [0, 0, 0])
my_model.set_hop(t2, 12, 22, [0, 0, 0])
my_model.set_hop(t1, 12, 25, [0, 0, 0])
my_model.set_hop(t1, 12, 33, [0, 0, 0])
my_model.set_hop(t2, 13, 18, [0, 0, 0])
my_model.set_hop(t2, 13, 23, [0, 0, 0])
my_model.set_hop(t1, 13, 24, [0, 0, 0])
my_model.set_hop(t1, 13, 32, [0, 0, 0])
my_model.set_hop(t1, 13, 43, [0, 0, 0])
my_model.set_hop(t2, 14, 3, [1, 0, 0])
my_model.set_hop(t2, 14, 23, [1, 0, 0])
my_model.set_hop(t1, 14, 40, [1, 0, 0])
my_model.set_hop(t2, 14, 17, [0, 0, 0])
my_model.set_hop(t2, 14, 18, [0, 0, 0])
my_model.set_hop(t1, 14, 29, [0, 0, 0])
my_model.set_hop(t1, 14, 39, [0, 0, 0])
my_model.set_hop(t2, 15, 19, [0, 0, 0])
my_model.set_hop(t2, 15, 22, [0, 0, 0])
my_model.set_hop(t1, 15, 28, [0, 0, 0])
my_model.set_hop(t1, 15, 38, [0, 0, 0])
my_model.set_hop(t1, 15, 41, [0, 0, 0])
my_model.set_hop(t1, 16, 27, [0, 0, 0])
my_model.set_hop(t1, 16, 47, [0, 0, 0])
my_model.set_hop(t1, 17, 26, [0, 0, 0])
my_model.set_hop(t1, 17, 39, [0, 0, 0])
my_model.set_hop(t1, 17, 46, [0, 0, 0])
my_model.set_hop(t1, 18, 29, [0, 0, 0])
my_model.set_hop(t1, 18, 43, [0, 0, 0])
my_model.set_hop(t1, 19, 28, [0, 0, 0])
my_model.set_hop(t1, 19, 37, [0, 0, 0])
my_model.set_hop(t1, 19, 42, [0, 0, 0])
my_model.set_hop(t2, 20, 3, [0, 1, 0])
my_model.set_hop(t2, 20, 13, [0, 1, 0])
my_model.set_hop(t1, 20, 32, [0, 1, 0])
my_model.set_hop(t1, 20, 31, [0, 0, 0])
my_model.set_hop(t1, 20, 45, [0, 0, 0])
my_model.set_hop(t1, 21, 30, [0, 0, 0])
my_model.set_hop(t1, 21, 33, [0, 0, 0])
my_model.set_hop(t1, 21, 44, [0, 0, 0])
my_model.set_hop(t2, 22, 7, [0, 1, 0])
my_model.set_hop(t2, 22, 11, [0, 1, 0])
my_model.set_hop(t1, 22, 34, [0, 1, 0])
my_model.set_hop(t1, 22, 25, [0, 0, 0])
my_model.set_hop(t1, 22, 41, [0, 0, 0])
my_model.set_hop(t1, 23, 24, [0, 0, 0])
my_model.set_hop(t1, 23, 35, [0, 0, 0])
my_model.set_hop(t1, 23, 40, [0, 0, 0])
my_model.set_hop(t1, 24, 0, [0, 0, 1])
my_model.set_hop(t2, 24, 37, [0, 0, 1])
my_model.set_hop(t2, 24, 47, [0, 0, 1])
my_model.set_hop(t2, 24, 32, [0, 0, 0])
my_model.set_hop(t2, 24, 35, [0, 0, 0])
my_model.set_hop(t2, 24, 40, [0, 0, 0])
my_model.set_hop(t2, 24, 43, [0, 0, 0])
my_model.set_hop(t2, 25, 33, [0, 0, 0])
my_model.set_hop(t2, 25, 36, [0, 0, 0])
my_model.set_hop(t2, 25, 41, [0, 0, 0])
my_model.set_hop(t2, 25, 46, [0, 0, 0])
my_model.set_hop(t1, 26, 2, [0, 0, 1])
my_model.set_hop(t2, 26, 33, [0, 0, 1])
my_model.set_hop(t2, 26, 41, [0, 0, 1])
my_model.set_hop(t2, 26, 36, [0, 0, 0])
my_model.set_hop(t2, 26, 39, [0, 0, 0])
my_model.set_hop(t2, 26, 45, [0, 0, 0])
my_model.set_hop(t2, 26, 46, [0, 0, 0])
my_model.set_hop(t2, 27, 32, [0, 0, 0])
my_model.set_hop(t2, 27, 37, [0, 0, 0])
my_model.set_hop(t2, 27, 40, [0, 0, 0])
my_model.set_hop(t2, 27, 47, [0, 0, 0])
my_model.set_hop(t2, 28, 37, [0, 0, 0])
my_model.set_hop(t2, 28, 38, [0, 0, 0])
my_model.set_hop(t2, 28, 41, [0, 0, 0])
my_model.set_hop(t2, 28, 42, [0, 0, 0])
my_model.set_hop(t2, 29, 34, [0, 0, 0])
my_model.set_hop(t2, 29, 39, [0, 0, 0])
my_model.set_hop(t2, 29, 43, [0, 0, 0])
my_model.set_hop(t2, 29, 44, [0, 0, 0])
my_model.set_hop(t2, 30, 33, [0, 0, 0])
my_model.set_hop(t2, 30, 34, [0, 0, 0])
my_model.set_hop(t2, 30, 44, [0, 0, 0])
my_model.set_hop(t2, 30, 47, [0, 0, 0])
my_model.set_hop(t2, 31, 35, [0, 0, 0])
my_model.set_hop(t2, 31, 38, [0, 0, 0])
my_model.set_hop(t2, 31, 42, [0, 0, 0])
my_model.set_hop(t2, 31, 45, [0, 0, 0])
my_model.set_hop(t2, 32, 40, [0, 0, 0])
my_model.set_hop(t2, 32, 43, [0, 0, 0])
my_model.set_hop(t2, 33, 41, [0, 0, 0])
my_model.set_hop(t2, 33, 44, [0, 0, 0])
my_model.set_hop(t2, 34, 44, [0, 0, 0])
my_model.set_hop(t2, 34, 47, [0, 0, 0])
my_model.set_hop(t2, 35, 40, [0, 0, 0])
my_model.set_hop(t2, 35, 45, [0, 0, 0])
my_model.set_hop(t1, 36, 18, [0, 1, 0])
my_model.set_hop(t2, 36, 29, [0, 1, 0])
my_model.set_hop(t2, 36, 43, [0, 1, 0])
my_model.set_hop(t2, 36, 45, [0, 0, 0])
my_model.set_hop(t2, 36, 46, [0, 0, 0])
my_model.set_hop(t2, 37, 42, [0, 0, 0])
my_model.set_hop(t2, 37, 47, [0, 0, 0])
my_model.set_hop(t1, 38, 16, [0, 1, 0])
my_model.set_hop(t2, 38, 27, [0, 1, 0])
my_model.set_hop(t2, 38, 47, [0, 1, 0])
my_model.set_hop(t2, 38, 41, [0, 0, 0])
my_model.set_hop(t2, 38, 42, [0, 0, 0])
my_model.set_hop(t2, 39, 43, [0, 0, 0])
my_model.set_hop(t2, 39, 46, [0, 0, 0])
my_model.set_hop(t1, 44, 8, [1, 0, 0])
my_model.set_hop(t2, 44, 27, [1, 0, 0])
my_model.set_hop(t2, 44, 37, [1, 0, 0])
my_model.set_hop(t1, 46, 10, [1, 0, 0])
my_model.set_hop(t2, 46, 31, [1, 0, 0])
my_model.set_hop(t2, 46, 35, [1, 0, 0])

my_model.display()

path = [ [0., 0., 0.], [-.5, .5, .5], [0., .5, 0.], [.25, .25, .25], [0., 0., 0.], [0., .5, .5] ]
label = (r'$\Gamma $',r'$H$', r'$N$', r'$P$', r'$\Gamma $', r'$N$')

nk = 800

(k_vec, k_dist, k_node) = my_model.k_path(path, nk)

print('-'*20)
print('starting calculation')
print('-'*20)
print('Calculating bands . . .')

evals = my_model.solve_all(k_vec)

fig, ax = plt.subplots(figsize=(3,5))

# ax.set_xlim([0, k_node[-1]])
ax.set_xticks(k_node)
ax.set_xticklabels(label)

for n in range(len(k_node)):
	ax.axvline(x=k_node[n], linewidth=0.5, color='k')

ax.set_title('P8-0 schwarzite band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(48):
	ax.plot(k_dist, evals[xband])

# fig.tight_layout()
fig.savefig('P8-0.pdf')

print('Done.')