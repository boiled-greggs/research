# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[6.754998207, 0.000000000, 0.000000000], \
	[-2.251666212, 6.368673335, 0.000000000], \
	[-2.251666212, -3.184336971, 5.515432721]]

orb = [[0.619104012, 0.395019974, 0.395026097], \
	[0.375860610, 0.773400326, 0.773410639], \
	[0.997483430, 0.221453756, 0.599838797], \
	[0.997485602, 0.599839847, 0.221467194], \
	[0.395081484, 0.619049976, 0.395035391], \
	[0.773451852, 0.375807501, 0.773402973], \
	[0.599886527, 0.997426169, 0.221460862], \
	[0.221516145, 0.997431002, 0.599847280], \
	[0.395079566, 0.395024846, 0.619055400], \
	[0.773454623, 0.773399899, 0.375815100], \
	[0.221510507, 0.599829118, 0.997436898], \
	[0.599893695, 0.221455460, 0.997436500], \
	[0.619101506, 0.221456506, 0.221462402], \
	[0.375869227, 0.599842374, 0.599851478], \
	[0.997486098, 0.773399102, 0.395022869], \
	[0.997483615, 0.395011414, 0.773407006], \
	[0.395067008, 0.773398951, 0.997435091], \
	[0.773456727, 0.395019029, 0.997435522], \
	[0.599897792, 0.599841510, 0.375820141], \
	[0.221510690, 0.221454441, 0.619050282], \
	[0.395073111, 0.997429387, 0.773408452], \
	[0.773453702, 0.997426381, 0.395019117], \
	[0.221511127, 0.619041044, 0.221463314], \
	[0.599898667, 0.375810620, 0.599848241]]

my_model = tb_model(3, 3, lat, orb)

# set model parameters
delta = -0.5
t1 = -2.8
t2 = -0.09
t3 = -.3

my_model.set_onsite([delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta])

# set hopping parameters for connected orbitals
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(0.3716, 0, 4, [0, 0, 0])
my_model.set_hop(0.3648, 0, 5, [0, 0, 0])
my_model.set_hop(0.3716, 0, 8, [0, 0, 0])
my_model.set_hop(0.3648, 0, 9, [0, 0, 0])
my_model.set_hop(3.3526, 0, 12, [0, 0, 0])
my_model.set_hop(0.2338, 0, 13, [0, 0, 0])
my_model.set_hop(0.2203, 0, 14, [0, 0, 0])
my_model.set_hop(0.2203, 0, 15, [0, 0, 0])
my_model.set_hop(2.5463, 0, 18, [0, 0, 0])
my_model.set_hop(2.5462, 0, 23, [0, 0, 0])
my_model.set_hop(0.3648, 1, 4, [0, 0, 0])
my_model.set_hop(0.3648, 1, 7, [0, 0, 0])
my_model.set_hop(0.3648, 1, 8, [0, 0, 0])
my_model.set_hop(0.3648, 1, 10, [0, 0, 0])
my_model.set_hop(3.3526, 1, 13, [0, 0, 0])
my_model.set_hop(2.5464, 1, 16, [0, 0, 0])
my_model.set_hop(2.5463, 1, 20, [0, 0, 0])
my_model.set_hop(0.3648, 2, 8, [1, 0, 0])
my_model.set_hop(0.3648, 2, 10, [1, 0, 0])
my_model.set_hop(0.2203, 2, 13, [1, 0, 0])
my_model.set_hop(2.5464, 2, 19, [1, 0, 0])
my_model.set_hop(0.3648, 2, 5, [0, 0, 0])
my_model.set_hop(0.2203, 2, 12, [0, 0, 0])
my_model.set_hop(3.3526, 2, 15, [0, 0, 0])
my_model.set_hop(0.3648, 3, 4, [1, 0, 0])
my_model.set_hop(0.3648, 3, 7, [1, 0, 0])
my_model.set_hop(0.2203, 3, 13, [1, 0, 0])
my_model.set_hop(2.5463, 3, 22, [1, 0, 0])
my_model.set_hop(0.3648, 3, 9, [0, 0, 0])
my_model.set_hop(0.2203, 3, 12, [0, 0, 0])
my_model.set_hop(3.3529, 3, 14, [0, 0, 0])
my_model.set_hop(0.3717, 4, 8, [0, 0, 0])
my_model.set_hop(0.3648, 4, 9, [0, 0, 0])
my_model.set_hop(2.5464, 4, 13, [0, 0, 0])
my_model.set_hop(2.5463, 4, 18, [0, 0, 0])
my_model.set_hop(0.2203, 4, 20, [0, 0, 0])
my_model.set_hop(0.2203, 4, 21, [0, 0, 0])
my_model.set_hop(3.3524, 4, 22, [0, 0, 0])
my_model.set_hop(0.2338, 4, 23, [0, 0, 0])
my_model.set_hop(0.3648, 5, 8, [0, 0, 0])
my_model.set_hop(0.3648, 5, 11, [0, 0, 0])
my_model.set_hop(2.5462, 5, 15, [0, 0, 0])
my_model.set_hop(2.5462, 5, 17, [0, 0, 0])
my_model.set_hop(3.3530, 5, 23, [0, 0, 0])
my_model.set_hop(0.3648, 6, 0, [0, 1, 0])
my_model.set_hop(0.3648, 6, 2, [0, 1, 0])
my_model.set_hop(2.5463, 6, 12, [0, 1, 0])
my_model.set_hop(0.2203, 6, 23, [0, 1, 0])
my_model.set_hop(0.3648, 6, 9, [0, 0, 0])
my_model.set_hop(3.3526, 6, 21, [0, 0, 0])
my_model.set_hop(0.2203, 6, 22, [0, 0, 0])
my_model.set_hop(0.3648, 7, 8, [0, 1, 0])
my_model.set_hop(0.3648, 7, 11, [0, 1, 0])
my_model.set_hop(2.5464, 7, 19, [0, 1, 0])
my_model.set_hop(0.2203, 7, 23, [0, 1, 0])
my_model.set_hop(3.3528, 7, 20, [0, 0, 0])
my_model.set_hop(0.2203, 7, 22, [0, 0, 0])
my_model.set_hop(2.5464, 8, 13, [0, 0, 0])
my_model.set_hop(0.2203, 8, 16, [0, 0, 0])
my_model.set_hop(0.2203, 8, 17, [0, 0, 0])
my_model.set_hop(0.2338, 8, 18, [0, 0, 0])
my_model.set_hop(3.3524, 8, 19, [0, 0, 0])
my_model.set_hop(2.5463, 8, 23, [0, 0, 0])
my_model.set_hop(2.5462, 9, 14, [0, 0, 0])
my_model.set_hop(3.3528, 9, 18, [0, 0, 0])
my_model.set_hop(2.5463, 9, 21, [0, 0, 0])
my_model.set_hop(0.3648, 10, 4, [0, 0, 1])
my_model.set_hop(0.3648, 10, 6, [0, 0, 1])
my_model.set_hop(0.2203, 10, 18, [0, 0, 1])
my_model.set_hop(2.5464, 10, 22, [0, 0, 1])
my_model.set_hop(3.3525, 10, 16, [0, 0, 0])
my_model.set_hop(0.2203, 10, 19, [0, 0, 0])
my_model.set_hop(0.3648, 11, 0, [0, 0, 1])
my_model.set_hop(0.3648, 11, 3, [0, 0, 1])
my_model.set_hop(2.5464, 11, 12, [0, 0, 1])
my_model.set_hop(0.2203, 11, 18, [0, 0, 1])
my_model.set_hop(3.3525, 11, 17, [0, 0, 0])
my_model.set_hop(0.2203, 11, 19, [0, 0, 0])
my_model.set_hop(0.3648, 12, 18, [0, 0, 0])
my_model.set_hop(0.3648, 12, 23, [0, 0, 0])
my_model.set_hop(0.3648, 13, 16, [0, 0, 0])
my_model.set_hop(0.3716, 13, 18, [0, 0, 0])
my_model.set_hop(0.3648, 13, 19, [0, 0, 0])
my_model.set_hop(0.3648, 13, 20, [0, 0, 0])
my_model.set_hop(0.3648, 13, 22, [0, 0, 0])
my_model.set_hop(0.3716, 13, 23, [0, 0, 0])
my_model.set_hop(0.2203, 14, 1, [1, 0, 0])
my_model.set_hop(2.5462, 14, 7, [1, 0, 0])
my_model.set_hop(0.3648, 14, 20, [1, 0, 0])
my_model.set_hop(0.3648, 14, 22, [1, 0, 0])
my_model.set_hop(0.3648, 14, 18, [0, 0, 0])
my_model.set_hop(0.3716, 14, 21, [0, 0, 0])
my_model.set_hop(0.2203, 15, 1, [1, 0, 0])
my_model.set_hop(2.5463, 15, 10, [1, 0, 0])
my_model.set_hop(0.3648, 15, 16, [1, 0, 0])
my_model.set_hop(0.3648, 15, 19, [1, 0, 0])
my_model.set_hop(0.3716, 15, 17, [0, 0, 0])
my_model.set_hop(0.3648, 15, 23, [0, 0, 0])
my_model.set_hop(2.5464, 16, 6, [0, 0, 1])
my_model.set_hop(0.2203, 16, 9, [0, 0, 1])
my_model.set_hop(0.3648, 16, 21, [0, 0, 1])
my_model.set_hop(0.3648, 16, 22, [0, 0, 1])
my_model.set_hop(0.3716, 16, 20, [0, 0, 0])
my_model.set_hop(2.5462, 17, 3, [0, 0, 1])
my_model.set_hop(0.2203, 17, 9, [0, 0, 1])
my_model.set_hop(0.3648, 17, 12, [0, 0, 1])
my_model.set_hop(0.3648, 17, 14, [0, 0, 1])
my_model.set_hop(0.3648, 17, 23, [0, 0, 0])
my_model.set_hop(0.3648, 18, 21, [0, 0, 0])
my_model.set_hop(0.3648, 18, 22, [0, 0, 0])
my_model.set_hop(0.3716, 18, 23, [0, 0, 0])
my_model.set_hop(0.3648, 19, 23, [0, 0, 0])
my_model.set_hop(0.2203, 20, 5, [0, 1, 0])
my_model.set_hop(2.5463, 20, 11, [0, 1, 0])
my_model.set_hop(0.3648, 20, 17, [0, 1, 0])
my_model.set_hop(0.3648, 20, 19, [0, 1, 0])
my_model.set_hop(2.5463, 21, 2, [0, 1, 0])
my_model.set_hop(0.2203, 21, 5, [0, 1, 0])
my_model.set_hop(0.3648, 21, 12, [0, 1, 0])
my_model.set_hop(0.3648, 21, 15, [0, 1, 0])

my_model.display()

path = [ [0., 0., 0.], [0.5, -0.5, 0.5], [0., 0., 0.5], [.25, .25, .25], [0., 0., 0.], [0., 0., 0.5] ]
label = (r'$\Gamma $',r'$H$', r'$N$', r'$P$', r'$\Gamma $', r'$N$')
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

ax.set_title('bccP80 band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(24):
	ax.plot(k_dist, evals[xband], linewidth=.5, color='black')

# fig.tight_layout()
fig.savefig('bccP80-1.eps')

print('Done.')
