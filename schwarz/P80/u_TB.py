# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[3.8900001049, 3.8900001049, -3.8900001049], \
	[-3.8900001049, 3.8900001049, 3.8900001049], \
	[3.8900001049, -3.8900001049, 3.8900001049]]

orb = [[0.619701982, 0.396404892, 0.396371782], \
	[0.375096351, 0.774092138, 0.774078548], \
	[0.997387648, 0.221265852, 0.597892642], \
	[0.997410417, 0.597927332, 0.221253812], \
	[0.396882981, 0.620504320, 0.396626592], \
	[0.774206698, 0.374853492, 0.773950338], \
	[0.598043561, 0.997165978, 0.221125662], \
	[0.220591605, 0.997166038, 0.598022342], \
	[0.396914274, 0.396661282, 0.620500922], \
	[0.774216413, 0.773963869, 0.374849558], \
	[0.220710203, 0.598055601, 0.997251391], \
	[0.599280596, 0.221009418, 0.997230589], \
	[0.618433774, 0.220753014, 0.220740259], \
	[0.376364619, 0.599728763, 0.599694848], \
	[0.997454107, 0.773643434, 0.395724863], \
	[0.997408092, 0.395693183, 0.773629069], \
	[0.395645946, 0.773835659, 0.997695744], \
	[0.774172723, 0.395500869, 0.997631311], \
	[0.599280596, 0.599472344, 0.375693947], \
	[0.220689818, 0.220881268, 0.618767858], \
	[0.395614713, 0.997678816, 0.773821414], \
	[0.774186373, 0.997678876, 0.395466954], \
	[0.220612064, 0.618702888, 0.220868409], \
	[0.599311769, 0.375629306, 0.599568188]]

my_model = tb_model(3, 3, lat, orb)

# set model parameters
delta = 0.0
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
my_model.set_hop(t3, 0, 2, [0, 0, 0])
my_model.set_hop(t3, 0, 3, [0, 0, 0])
my_model.set_hop(t2, 0, 4, [0, 0, 0])
my_model.set_hop(t2, 0, 5, [0, 0, 0])
my_model.set_hop(t2, 0, 8, [0, 0, 0])
my_model.set_hop(t2, 0, 9, [0, 0, 0])
my_model.set_hop(t1, 0, 12, [0, 0, 0])
my_model.set_hop(t3, 0, 13, [0, 0, 0])
my_model.set_hop(t3, 0, 14, [0, 0, 0])
my_model.set_hop(t3, 0, 15, [0, 0, 0])
my_model.set_hop(t1, 0, 18, [0, 0, 0])
my_model.set_hop(t3, 0, 19, [0, 0, 0])
my_model.set_hop(t3, 0, 22, [0, 0, 0])
my_model.set_hop(t1, 0, 23, [0, 0, 0])
my_model.set_hop(t2, 1, 4, [0, 0, 0])
my_model.set_hop(t2, 1, 7, [0, 0, 0])
my_model.set_hop(t2, 1, 8, [0, 0, 0])
my_model.set_hop(t2, 1, 10, [0, 0, 0])
my_model.set_hop(t1, 1, 13, [0, 0, 0])
my_model.set_hop(t1, 1, 16, [0, 0, 0])
my_model.set_hop(t3, 1, 18, [0, 0, 0])
my_model.set_hop(t3, 1, 19, [0, 0, 0])
my_model.set_hop(t1, 1, 20, [0, 0, 0])
my_model.set_hop(t3, 1, 22, [0, 0, 0])
my_model.set_hop(t3, 1, 23, [0, 0, 0])
my_model.set_hop(t3, 2, 1, [1, 0, 0])
my_model.set_hop(t2, 2, 8, [1, 0, 0])
my_model.set_hop(t2, 2, 10, [1, 0, 0])
my_model.set_hop(t3, 2, 13, [1, 0, 0])
my_model.set_hop(t3, 2, 16, [1, 0, 0])
my_model.set_hop(t1, 2, 19, [1, 0, 0])
my_model.set_hop(t2, 2, 5, [0, 0, 0])
my_model.set_hop(t3, 2, 12, [0, 0, 0])
my_model.set_hop(t1, 2, 15, [0, 0, 0])
my_model.set_hop(t3, 2, 17, [0, 0, 0])
my_model.set_hop(t3, 2, 23, [0, 0, 0])
my_model.set_hop(t3, 3, 1, [1, 0, 0])
my_model.set_hop(t2, 3, 4, [1, 0, 0])
my_model.set_hop(t2, 3, 7, [1, 0, 0])
my_model.set_hop(t3, 3, 13, [1, 0, 0])
my_model.set_hop(t3, 3, 20, [1, 0, 0])
my_model.set_hop(t1, 3, 22, [1, 0, 0])
my_model.set_hop(t2, 3, 9, [0, 0, 0])
my_model.set_hop(t3, 3, 12, [0, 0, 0])
my_model.set_hop(t1, 3, 14, [0, 0, 0])
my_model.set_hop(t3, 3, 18, [0, 0, 0])
my_model.set_hop(t3, 3, 21, [0, 0, 0])
my_model.set_hop(t3, 4, 6, [0, 0, 0])
my_model.set_hop(t3, 4, 7, [0, 0, 0])
my_model.set_hop(t2, 4, 8, [0, 0, 0])
my_model.set_hop(t2, 4, 9, [0, 0, 0])
my_model.set_hop(t3, 4, 12, [0, 0, 0])
my_model.set_hop(t1, 4, 13, [0, 0, 0])
my_model.set_hop(t1, 4, 18, [0, 0, 0])
my_model.set_hop(t3, 4, 19, [0, 0, 0])
my_model.set_hop(t3, 4, 20, [0, 0, 0])
my_model.set_hop(t3, 4, 21, [0, 0, 0])
my_model.set_hop(t1, 4, 22, [0, 0, 0])
my_model.set_hop(t3, 4, 23, [0, 0, 0])
my_model.set_hop(t2, 5, 8, [0, 0, 0])
my_model.set_hop(t2, 5, 11, [0, 0, 0])
my_model.set_hop(t3, 5, 12, [0, 0, 0])
my_model.set_hop(t3, 5, 13, [0, 0, 0])
my_model.set_hop(t1, 5, 15, [0, 0, 0])
my_model.set_hop(t1, 5, 17, [0, 0, 0])
my_model.set_hop(t3, 5, 18, [0, 0, 0])
my_model.set_hop(t3, 5, 19, [0, 0, 0])
my_model.set_hop(t1, 5, 23, [0, 0, 0])
my_model.set_hop(t2, 6, 0, [0, 1, 0])
my_model.set_hop(t2, 6, 2, [0, 1, 0])
my_model.set_hop(t3, 6, 5, [0, 1, 0])
my_model.set_hop(t1, 6, 12, [0, 1, 0])
my_model.set_hop(t3, 6, 15, [0, 1, 0])
my_model.set_hop(t3, 6, 23, [0, 1, 0])
my_model.set_hop(t2, 6, 9, [0, 0, 0])
my_model.set_hop(t3, 6, 14, [0, 0, 0])
my_model.set_hop(t3, 6, 18, [0, 0, 0])
my_model.set_hop(t1, 6, 21, [0, 0, 0])
my_model.set_hop(t3, 6, 22, [0, 0, 0])
my_model.set_hop(t3, 7, 5, [0, 1, 0])
my_model.set_hop(t2, 7, 8, [0, 1, 0])
my_model.set_hop(t2, 7, 11, [0, 1, 0])
my_model.set_hop(t3, 7, 17, [0, 1, 0])
my_model.set_hop(t1, 7, 19, [0, 1, 0])
my_model.set_hop(t3, 7, 23, [0, 1, 0])
my_model.set_hop(t3, 7, 13, [0, 0, 0])
my_model.set_hop(t3, 7, 16, [0, 0, 0])
my_model.set_hop(t1, 7, 20, [0, 0, 0])
my_model.set_hop(t3, 7, 22, [0, 0, 0])
my_model.set_hop(t3, 8, 10, [0, 0, 0])
my_model.set_hop(t3, 8, 11, [0, 0, 0])
my_model.set_hop(t3, 8, 12, [0, 0, 0])
my_model.set_hop(t1, 8, 13, [0, 0, 0])
my_model.set_hop(t3, 8, 16, [0, 0, 0])
my_model.set_hop(t3, 8, 17, [0, 0, 0])
my_model.set_hop(t3, 8, 18, [0, 0, 0])
my_model.set_hop(t1, 8, 19, [0, 0, 0])
my_model.set_hop(t3, 8, 22, [0, 0, 0])
my_model.set_hop(t1, 8, 23, [0, 0, 0])
my_model.set_hop(t3, 9, 12, [0, 0, 0])
my_model.set_hop(t3, 9, 13, [0, 0, 0])
my_model.set_hop(t1, 9, 14, [0, 0, 0])
my_model.set_hop(t1, 9, 18, [0, 0, 0])
my_model.set_hop(t1, 9, 21, [0, 0, 0])
my_model.set_hop(t3, 9, 22, [0, 0, 0])
my_model.set_hop(t3, 9, 23, [0, 0, 0])
my_model.set_hop(t2, 10, 4, [0, 0, 1])
my_model.set_hop(t2, 10, 6, [0, 0, 1])
my_model.set_hop(t3, 10, 9, [0, 0, 1])
my_model.set_hop(t3, 10, 18, [0, 0, 1])
my_model.set_hop(t3, 10, 21, [0, 0, 1])
my_model.set_hop(t1, 10, 22, [0, 0, 1])
my_model.set_hop(t3, 10, 13, [0, 0, 0])
my_model.set_hop(t1, 10, 16, [0, 0, 0])
my_model.set_hop(t3, 10, 19, [0, 0, 0])
my_model.set_hop(t3, 10, 20, [0, 0, 0])
my_model.set_hop(t2, 11, 0, [0, 0, 1])
my_model.set_hop(t2, 11, 3, [0, 0, 1])
my_model.set_hop(t3, 11, 9, [0, 0, 1])
my_model.set_hop(t1, 11, 12, [0, 0, 1])
my_model.set_hop(t3, 11, 14, [0, 0, 1])
my_model.set_hop(t3, 11, 18, [0, 0, 1])
my_model.set_hop(t3, 11, 15, [0, 0, 0])
my_model.set_hop(t1, 11, 17, [0, 0, 0])
my_model.set_hop(t3, 11, 19, [0, 0, 0])
my_model.set_hop(t3, 11, 23, [0, 0, 0])
my_model.set_hop(t3, 12, 14, [0, 0, 0])
my_model.set_hop(t3, 12, 15, [0, 0, 0])
my_model.set_hop(t2, 12, 18, [0, 0, 0])
my_model.set_hop(t2, 12, 23, [0, 0, 0])
my_model.set_hop(t2, 13, 16, [0, 0, 0])
my_model.set_hop(t2, 13, 18, [0, 0, 0])
my_model.set_hop(t2, 13, 19, [0, 0, 0])
my_model.set_hop(t2, 13, 20, [0, 0, 0])
my_model.set_hop(t2, 13, 22, [0, 0, 0])
my_model.set_hop(t2, 13, 23, [0, 0, 0])
my_model.set_hop(t3, 14, 1, [1, 0, 0])
my_model.set_hop(t3, 14, 4, [1, 0, 0])
my_model.set_hop(t1, 14, 7, [1, 0, 0])
my_model.set_hop(t3, 14, 13, [1, 0, 0])
my_model.set_hop(t2, 14, 20, [1, 0, 0])
my_model.set_hop(t2, 14, 22, [1, 0, 0])
my_model.set_hop(t2, 14, 18, [0, 0, 0])
my_model.set_hop(t2, 14, 21, [0, 0, 0])
my_model.set_hop(t3, 15, 1, [1, 0, 0])
my_model.set_hop(t3, 15, 8, [1, 0, 0])
my_model.set_hop(t1, 15, 10, [1, 0, 0])
my_model.set_hop(t3, 15, 13, [1, 0, 0])
my_model.set_hop(t2, 15, 16, [1, 0, 0])
my_model.set_hop(t2, 15, 19, [1, 0, 0])
my_model.set_hop(t2, 15, 17, [0, 0, 0])
my_model.set_hop(t2, 15, 23, [0, 0, 0])
my_model.set_hop(t3, 16, 4, [0, 0, 1])
my_model.set_hop(t1, 16, 6, [0, 0, 1])
my_model.set_hop(t3, 16, 9, [0, 0, 1])
my_model.set_hop(t3, 16, 18, [0, 0, 1])
my_model.set_hop(t2, 16, 21, [0, 0, 1])
my_model.set_hop(t2, 16, 22, [0, 0, 1])
my_model.set_hop(t3, 16, 19, [0, 0, 0])
my_model.set_hop(t2, 16, 20, [0, 0, 0])
my_model.set_hop(t3, 17, 0, [0, 0, 1])
my_model.set_hop(t1, 17, 3, [0, 0, 1])
my_model.set_hop(t3, 17, 9, [0, 0, 1])
my_model.set_hop(t2, 17, 12, [0, 0, 1])
my_model.set_hop(t2, 17, 14, [0, 0, 1])
my_model.set_hop(t3, 17, 18, [0, 0, 1])
my_model.set_hop(t3, 17, 19, [0, 0, 0])
my_model.set_hop(t2, 17, 23, [0, 0, 0])
my_model.set_hop(t2, 18, 21, [0, 0, 0])
my_model.set_hop(t2, 18, 22, [0, 0, 0])
my_model.set_hop(t2, 18, 23, [0, 0, 0])
my_model.set_hop(t2, 19, 23, [0, 0, 0])
my_model.set_hop(t3, 20, 5, [0, 1, 0])
my_model.set_hop(t3, 20, 8, [0, 1, 0])
my_model.set_hop(t1, 20, 11, [0, 1, 0])
my_model.set_hop(t2, 20, 17, [0, 1, 0])
my_model.set_hop(t2, 20, 19, [0, 1, 0])
my_model.set_hop(t3, 20, 23, [0, 1, 0])
my_model.set_hop(t3, 20, 22, [0, 0, 0])
my_model.set_hop(t3, 21, 0, [0, 1, 0])
my_model.set_hop(t1, 21, 2, [0, 1, 0])
my_model.set_hop(t3, 21, 5, [0, 1, 0])
my_model.set_hop(t2, 21, 12, [0, 1, 0])
my_model.set_hop(t2, 21, 15, [0, 1, 0])
my_model.set_hop(t3, 21, 23, [0, 1, 0])
my_model.set_hop(t3, 21, 22, [0, 0, 0])

my_model.display()

path = [ [0., 0., 0.], [0., 0.5, 0.], [0.25, .25, 0.], [.25, .25, .25], [0., 0., 0.], [0.25, 0.25, 0.] ]
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

ax.set_title('u band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(6):
	ax.plot(k_dist, evals[xband+10], linewidth=.5) #, color='black')

# fig.tight_layout()
fig.savefig('u.eps')

print('Done.')
