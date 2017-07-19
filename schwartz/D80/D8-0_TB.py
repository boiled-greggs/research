# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[6.0430002213, 0.0000000000, 0.0000000000], \
	[0.0000000000, 6.0430002213, 0.0000000000], \
	[0.0000000000, 0.0000000000, 6.0430002213]]

orb = [[0.499945283, 0.336606652, 0.663449407], \
	[0.000000000, 0.163504109, 0.663449407], \
	[0.000000000, 0.336606652, 0.836551905], \
	[0.499945283, 0.163504109, 0.836551905], \
	[0.663449407, 0.499945283, 0.336606652], \
	[0.663449407, 0.000000000, 0.163504109], \
	[0.836551905, 0.000000000, 0.336606652], \
	[0.836551905, 0.499945283, 0.163504109], \
	[0.336606652, 0.663449407, 0.499945283], \
	[0.163504109, 0.663449407, 0.000000000], \
	[0.336606652, 0.836551905, 0.000000000], \
	[0.163504109, 0.836551905, 0.499945283], \
	[0.499945283, 0.663449407, 0.336606652], \
	[0.000000000, 0.836551905, 0.336606652], \
	[0.000000000, 0.663449407, 0.163504109], \
	[0.499945283, 0.836551905, 0.163504109], \
	[0.336606652, 0.499945283, 0.663449407], \
	[0.336606652, 0.000000000, 0.836551905], \
	[0.163504109, 0.000000000, 0.663449407], \
	[0.163504109, 0.499945283, 0.836551905], \
	[0.663449407, 0.336606652, 0.499945283], \
	[0.836551905, 0.336606652, 0.000000000], \
	[0.663449407, 0.163504109, 0.000000000], \
	[0.836551905, 0.163504109, 0.499945283]]

my_model = tb_model(3, 3, lat, orb)

# set model parameters
delta = 0.0
t1 = 1.8
t2 = 0.49
t3 = -.4

my_model.set_onsite([delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta])

# set hopping parameters for connected orbitals
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(t3, 0, 1, [0, 0, 0])
my_model.set_hop(t3, 0, 2, [0, 0, 0])
my_model.set_hop(t1, 0, 3, [0, 0, 0])
my_model.set_hop(t2, 0, 4, [0, 0, 0])
my_model.set_hop(t3, 0, 5, [0, 0, 0])
my_model.set_hop(t3, 0, 6, [0, 0, 0])
my_model.set_hop(t3, 0, 7, [0, 0, 0])
my_model.set_hop(t2, 0, 8, [0, 0, 0])
my_model.set_hop(t3, 0, 11, [0, 0, 0])
my_model.set_hop(t3, 0, 12, [0, 0, 0])
my_model.set_hop(t1, 0, 16, [0, 0, 0])
my_model.set_hop(t2, 0, 17, [0, 0, 0])
my_model.set_hop(t3, 0, 18, [0, 0, 0])
my_model.set_hop(t2, 0, 19, [0, 0, 0])
my_model.set_hop(t1, 0, 20, [0, 0, 0])
my_model.set_hop(t2, 0, 23, [0, 0, 0])
my_model.set_hop(t1, 1, 2, [0, 0, 0])
my_model.set_hop(t3, 1, 3, [0, 0, 0])
my_model.set_hop(t3, 1, 8, [0, 0, 0])
my_model.set_hop(t3, 1, 16, [0, 0, 0])
my_model.set_hop(t2, 1, 17, [0, 0, 0])
my_model.set_hop(t1, 1, 18, [0, 0, 0])
my_model.set_hop(t2, 1, 19, [0, 0, 0])
my_model.set_hop(t2, 2, 9, [0, 0, 1])
my_model.set_hop(t3, 2, 10, [0, 0, 1])
my_model.set_hop(t3, 2, 14, [0, 0, 1])
my_model.set_hop(t3, 2, 3, [0, 0, 0])
my_model.set_hop(t3, 2, 8, [0, 0, 0])
my_model.set_hop(t3, 2, 11, [0, 0, 0])
my_model.set_hop(t2, 2, 16, [0, 0, 0])
my_model.set_hop(t3, 2, 17, [0, 0, 0])
my_model.set_hop(t2, 2, 18, [0, 0, 0])
my_model.set_hop(t1, 2, 19, [0, 0, 0])
my_model.set_hop(t3, 3, 4, [0, 0, 1])
my_model.set_hop(t2, 3, 5, [0, 0, 1])
my_model.set_hop(t3, 3, 6, [0, 0, 1])
my_model.set_hop(t3, 3, 7, [0, 0, 1])
my_model.set_hop(t3, 3, 9, [0, 0, 1])
my_model.set_hop(t2, 3, 21, [0, 0, 1])
my_model.set_hop(t1, 3, 22, [0, 0, 1])
my_model.set_hop(t3, 3, 4, [0, 0, 0])
my_model.set_hop(t3, 3, 6, [0, 0, 0])
my_model.set_hop(t3, 3, 8, [0, 0, 0])
my_model.set_hop(t2, 3, 16, [0, 0, 0])
my_model.set_hop(t1, 3, 17, [0, 0, 0])
my_model.set_hop(t2, 3, 18, [0, 0, 0])
my_model.set_hop(t3, 3, 19, [0, 0, 0])
my_model.set_hop(t2, 3, 20, [0, 0, 0])
my_model.set_hop(t3, 3, 23, [0, 0, 0])
my_model.set_hop(t3, 4, 5, [0, 0, 0])
my_model.set_hop(t3, 4, 6, [0, 0, 0])
my_model.set_hop(t1, 4, 7, [0, 0, 0])
my_model.set_hop(t2, 4, 8, [0, 0, 0])
my_model.set_hop(t3, 4, 9, [0, 0, 0])
my_model.set_hop(t3, 4, 10, [0, 0, 0])
my_model.set_hop(t3, 4, 11, [0, 0, 0])
my_model.set_hop(t1, 4, 12, [0, 0, 0])
my_model.set_hop(t2, 4, 15, [0, 0, 0])
my_model.set_hop(t3, 4, 16, [0, 0, 0])
my_model.set_hop(t1, 4, 20, [0, 0, 0])
my_model.set_hop(t2, 4, 21, [0, 0, 0])
my_model.set_hop(t3, 4, 22, [0, 0, 0])
my_model.set_hop(t2, 4, 23, [0, 0, 0])
my_model.set_hop(t1, 5, 6, [0, 0, 0])
my_model.set_hop(t3, 5, 7, [0, 0, 0])
my_model.set_hop(t3, 5, 20, [0, 0, 0])
my_model.set_hop(t2, 5, 21, [0, 0, 0])
my_model.set_hop(t1, 5, 22, [0, 0, 0])
my_model.set_hop(t2, 5, 23, [0, 0, 0])
my_model.set_hop(t2, 6, 1, [1, 0, 0])
my_model.set_hop(t3, 6, 2, [1, 0, 0])
my_model.set_hop(t3, 6, 18, [1, 0, 0])
my_model.set_hop(t3, 6, 7, [0, 0, 0])
my_model.set_hop(t2, 6, 20, [0, 0, 0])
my_model.set_hop(t3, 6, 21, [0, 0, 0])
my_model.set_hop(t2, 6, 22, [0, 0, 0])
my_model.set_hop(t1, 6, 23, [0, 0, 0])
my_model.set_hop(t3, 7, 1, [1, 0, 0])
my_model.set_hop(t3, 7, 8, [1, 0, 0])
my_model.set_hop(t2, 7, 9, [1, 0, 0])
my_model.set_hop(t3, 7, 10, [1, 0, 0])
my_model.set_hop(t3, 7, 11, [1, 0, 0])
my_model.set_hop(t2, 7, 13, [1, 0, 0])
my_model.set_hop(t1, 7, 14, [1, 0, 0])
my_model.set_hop(t3, 7, 8, [0, 0, 0])
my_model.set_hop(t3, 7, 10, [0, 0, 0])
my_model.set_hop(t2, 7, 12, [0, 0, 0])
my_model.set_hop(t3, 7, 15, [0, 0, 0])
my_model.set_hop(t2, 7, 20, [0, 0, 0])
my_model.set_hop(t1, 7, 21, [0, 0, 0])
my_model.set_hop(t2, 7, 22, [0, 0, 0])
my_model.set_hop(t3, 7, 23, [0, 0, 0])
my_model.set_hop(t3, 8, 9, [0, 0, 0])
my_model.set_hop(t3, 8, 10, [0, 0, 0])
my_model.set_hop(t1, 8, 11, [0, 0, 0])
my_model.set_hop(t1, 8, 12, [0, 0, 0])
my_model.set_hop(t2, 8, 13, [0, 0, 0])
my_model.set_hop(t3, 8, 14, [0, 0, 0])
my_model.set_hop(t2, 8, 15, [0, 0, 0])
my_model.set_hop(t1, 8, 16, [0, 0, 0])
my_model.set_hop(t2, 8, 19, [0, 0, 0])
my_model.set_hop(t3, 8, 20, [0, 0, 0])
my_model.set_hop(t1, 9, 10, [0, 0, 0])
my_model.set_hop(t3, 9, 11, [0, 0, 0])
my_model.set_hop(t3, 9, 12, [0, 0, 0])
my_model.set_hop(t2, 9, 13, [0, 0, 0])
my_model.set_hop(t1, 9, 14, [0, 0, 0])
my_model.set_hop(t2, 9, 15, [0, 0, 0])
my_model.set_hop(t2, 10, 5, [0, 1, 0])
my_model.set_hop(t3, 10, 6, [0, 1, 0])
my_model.set_hop(t3, 10, 22, [0, 1, 0])
my_model.set_hop(t3, 10, 11, [0, 0, 0])
my_model.set_hop(t2, 10, 12, [0, 0, 0])
my_model.set_hop(t3, 10, 13, [0, 0, 0])
my_model.set_hop(t2, 10, 14, [0, 0, 0])
my_model.set_hop(t1, 10, 15, [0, 0, 0])
my_model.set_hop(t3, 11, 0, [0, 1, 0])
my_model.set_hop(t2, 11, 1, [0, 1, 0])
my_model.set_hop(t3, 11, 2, [0, 1, 0])
my_model.set_hop(t3, 11, 3, [0, 1, 0])
my_model.set_hop(t3, 11, 5, [0, 1, 0])
my_model.set_hop(t2, 11, 17, [0, 1, 0])
my_model.set_hop(t1, 11, 18, [0, 1, 0])
my_model.set_hop(t2, 11, 12, [0, 0, 0])
my_model.set_hop(t1, 11, 13, [0, 0, 0])
my_model.set_hop(t2, 11, 14, [0, 0, 0])
my_model.set_hop(t3, 11, 15, [0, 0, 0])
my_model.set_hop(t2, 11, 16, [0, 0, 0])
my_model.set_hop(t3, 11, 19, [0, 0, 0])
my_model.set_hop(t3, 12, 13, [0, 0, 0])
my_model.set_hop(t3, 12, 14, [0, 0, 0])
my_model.set_hop(t1, 12, 15, [0, 0, 0])
my_model.set_hop(t2, 12, 16, [0, 0, 0])
my_model.set_hop(t3, 12, 19, [0, 0, 0])
my_model.set_hop(t2, 12, 20, [0, 0, 0])
my_model.set_hop(t3, 12, 21, [0, 0, 0])
my_model.set_hop(t3, 12, 22, [0, 0, 0])
my_model.set_hop(t3, 12, 23, [0, 0, 0])
my_model.set_hop(t3, 13, 1, [0, 1, 0])
my_model.set_hop(t3, 13, 17, [0, 1, 0])
my_model.set_hop(t2, 13, 18, [0, 1, 0])
my_model.set_hop(t1, 13, 14, [0, 0, 0])
my_model.set_hop(t3, 13, 15, [0, 0, 0])
my_model.set_hop(t3, 13, 16, [0, 0, 0])
my_model.set_hop(t3, 13, 19, [0, 0, 0])
my_model.set_hop(t3, 14, 15, [0, 0, 0])
my_model.set_hop(t3, 14, 16, [0, 0, 0])
my_model.set_hop(t1, 15, 5, [0, 1, 0])
my_model.set_hop(t2, 15, 6, [0, 1, 0])
my_model.set_hop(t3, 15, 18, [0, 1, 0])
my_model.set_hop(t3, 15, 20, [0, 1, 0])
my_model.set_hop(t3, 15, 21, [0, 1, 0])
my_model.set_hop(t2, 15, 22, [0, 1, 0])
my_model.set_hop(t3, 15, 23, [0, 1, 0])
my_model.set_hop(t3, 15, 16, [0, 0, 0])
my_model.set_hop(t3, 15, 20, [0, 0, 0])
my_model.set_hop(t3, 15, 21, [0, 0, 0])
my_model.set_hop(t3, 16, 17, [0, 0, 0])
my_model.set_hop(t3, 16, 18, [0, 0, 0])
my_model.set_hop(t1, 16, 19, [0, 0, 0])
my_model.set_hop(t2, 16, 20, [0, 0, 0])
my_model.set_hop(t3, 16, 23, [0, 0, 0])
my_model.set_hop(t3, 17, 5, [0, 0, 1])
my_model.set_hop(t3, 17, 21, [0, 0, 1])
my_model.set_hop(t2, 17, 22, [0, 0, 1])
my_model.set_hop(t1, 17, 18, [0, 0, 0])
my_model.set_hop(t3, 17, 19, [0, 0, 0])
my_model.set_hop(t3, 17, 20, [0, 0, 0])
my_model.set_hop(t3, 17, 23, [0, 0, 0])
my_model.set_hop(t3, 18, 19, [0, 0, 0])
my_model.set_hop(t3, 18, 20, [0, 0, 0])
my_model.set_hop(t1, 19, 9, [0, 0, 1])
my_model.set_hop(t2, 19, 10, [0, 0, 1])
my_model.set_hop(t3, 19, 12, [0, 0, 1])
my_model.set_hop(t3, 19, 13, [0, 0, 1])
my_model.set_hop(t2, 19, 14, [0, 0, 1])
my_model.set_hop(t3, 19, 15, [0, 0, 1])
my_model.set_hop(t3, 19, 22, [0, 0, 1])
my_model.set_hop(t3, 19, 20, [0, 0, 0])
my_model.set_hop(t3, 20, 21, [0, 0, 0])
my_model.set_hop(t3, 20, 22, [0, 0, 0])
my_model.set_hop(t1, 20, 23, [0, 0, 0])
my_model.set_hop(t3, 21, 9, [1, 0, 0])
my_model.set_hop(t3, 21, 13, [1, 0, 0])
my_model.set_hop(t2, 21, 14, [1, 0, 0])
my_model.set_hop(t1, 21, 22, [0, 0, 0])
my_model.set_hop(t3, 21, 23, [0, 0, 0])
my_model.set_hop(t3, 22, 23, [0, 0, 0])
my_model.set_hop(t1, 23, 1, [1, 0, 0])
my_model.set_hop(t2, 23, 2, [1, 0, 0])
my_model.set_hop(t3, 23, 14, [1, 0, 0])
my_model.set_hop(t3, 23, 16, [1, 0, 0])
my_model.set_hop(t3, 23, 17, [1, 0, 0])
my_model.set_hop(t2, 23, 18, [1, 0, 0])
my_model.set_hop(t3, 23, 19, [1, 0, 0])

my_model.display()

path = [ [0., 0., 0.], [0., 0.5, 0.], [0.5, .5, 0.], [0., 0., 0.], [0.5, 0.5, 0.5], [0., 0.5, 0.] ]
label = (r'$\Gamma $',r'$X$', r'$M$', r'$\Gamma $', r'$R$', r'$X$')
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

ax.set_title('D8-0 band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(24):
	ax.plot(k_dist, evals[xband], linewidth=.5, color='black')

# fig.tight_layout()
fig.savefig('D8-0-11.eps')

print('Done.')
