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
my_model.set_hop(0.3784, 0, 4, [0, 0, 0])
my_model.set_hop(0.3662, 0, 5, [0, 0, 0])
my_model.set_hop(0.3784, 0, 8, [0, 0, 0])
my_model.set_hop(0.3663, 0, 9, [0, 0, 0])
my_model.set_hop(3.2936, 0, 12, [0, 0, 0])
my_model.set_hop(0.2373, 0, 13, [0, 0, 0])
my_model.set_hop(0.2224, 0, 14, [0, 0, 0])
my_model.set_hop(0.2224, 0, 15, [0, 0, 0])
my_model.set_hop(2.5973, 0, 18, [0, 0, 0])
my_model.set_hop(2.5913, 0, 23, [0, 0, 0])
my_model.set_hop(0.3665, 1, 4, [0, 0, 0])
my_model.set_hop(0.3662, 1, 7, [0, 0, 0])
my_model.set_hop(0.3665, 1, 8, [0, 0, 0])
my_model.set_hop(0.3661, 1, 10, [0, 0, 0])
my_model.set_hop(3.2935, 1, 13, [0, 0, 0])
my_model.set_hop(2.5935, 1, 16, [0, 0, 0])
my_model.set_hop(2.5944, 1, 20, [0, 0, 0])
my_model.set_hop(0.3669, 2, 8, [1, 0, 0])
my_model.set_hop(0.3671, 2, 10, [1, 0, 0])
my_model.set_hop(0.2224, 2, 13, [1, 0, 0])
my_model.set_hop(2.6072, 2, 19, [1, 0, 0])
my_model.set_hop(0.3673, 2, 5, [0, 0, 0])
my_model.set_hop(0.2224, 2, 12, [0, 0, 0])
my_model.set_hop(3.2898, 2, 15, [0, 0, 0])
my_model.set_hop(0.3670, 3, 4, [1, 0, 0])
my_model.set_hop(0.3673, 3, 7, [1, 0, 0])
my_model.set_hop(0.2224, 3, 13, [1, 0, 0])
my_model.set_hop(2.6101, 3, 22, [1, 0, 0])
my_model.set_hop(0.3673, 3, 9, [0, 0, 0])
my_model.set_hop(0.2224, 3, 12, [0, 0, 0])
my_model.set_hop(3.2898, 3, 14, [0, 0, 0])
my_model.set_hop(0.3760, 4, 8, [0, 0, 0])
my_model.set_hop(0.3669, 4, 9, [0, 0, 0])
my_model.set_hop(2.5944, 4, 13, [0, 0, 0])
my_model.set_hop(2.6109, 4, 18, [0, 0, 0])
my_model.set_hop(0.2224, 4, 20, [0, 0, 0])
my_model.set_hop(0.2224, 4, 21, [0, 0, 0])
my_model.set_hop(3.2872, 4, 22, [0, 0, 0])
my_model.set_hop(0.2363, 4, 23, [0, 0, 0])
my_model.set_hop(0.3670, 5, 8, [0, 0, 0])
my_model.set_hop(0.3688, 5, 11, [0, 0, 0])
my_model.set_hop(2.6117, 5, 15, [0, 0, 0])
my_model.set_hop(2.5953, 5, 17, [0, 0, 0])
my_model.set_hop(3.2925, 5, 23, [0, 0, 0])
my_model.set_hop(0.3663, 6, 0, [0, 1, 0])
my_model.set_hop(0.3674, 6, 2, [0, 1, 0])
my_model.set_hop(2.5912, 6, 12, [0, 1, 0])
my_model.set_hop(0.2224, 6, 23, [0, 1, 0])
my_model.set_hop(0.3668, 6, 9, [0, 0, 0])
my_model.set_hop(3.2947, 6, 21, [0, 0, 0])
my_model.set_hop(0.2224, 6, 22, [0, 0, 0])
my_model.set_hop(0.3665, 7, 8, [0, 1, 0])
my_model.set_hop(0.3645, 7, 11, [0, 1, 0])
my_model.set_hop(2.5968, 7, 19, [0, 1, 0])
my_model.set_hop(0.2224, 7, 23, [0, 1, 0])
my_model.set_hop(3.2850, 7, 20, [0, 0, 0])
my_model.set_hop(0.2224, 7, 22, [0, 0, 0])
my_model.set_hop(2.5935, 8, 13, [0, 0, 0])
my_model.set_hop(0.2224, 8, 16, [0, 0, 0])
my_model.set_hop(0.2224, 8, 17, [0, 0, 0])
my_model.set_hop(0.2365, 8, 18, [0, 0, 0])
my_model.set_hop(3.2859, 8, 19, [0, 0, 0])
my_model.set_hop(2.6109, 8, 23, [0, 0, 0])
my_model.set_hop(2.6107, 9, 14, [0, 0, 0])
my_model.set_hop(3.2861, 9, 18, [0, 0, 0])
my_model.set_hop(2.5937, 9, 21, [0, 0, 0])
my_model.set_hop(0.3668, 10, 4, [0, 0, 1])
my_model.set_hop(0.3669, 10, 6, [0, 0, 1])
my_model.set_hop(0.2224, 10, 18, [0, 0, 1])
my_model.set_hop(2.5968, 10, 22, [0, 0, 1])
my_model.set_hop(3.2861, 10, 16, [0, 0, 0])
my_model.set_hop(0.2224, 10, 19, [0, 0, 0])
my_model.set_hop(0.3646, 11, 0, [0, 0, 1])
my_model.set_hop(0.3695, 11, 3, [0, 0, 1])
my_model.set_hop(2.5840, 11, 12, [0, 0, 1])
my_model.set_hop(0.2215, 11, 18, [0, 0, 1])
my_model.set_hop(3.3240, 11, 17, [0, 0, 0])
my_model.set_hop(0.2215, 11, 19, [0, 0, 0])
my_model.set_hop(0.3670, 12, 18, [0, 0, 0])
my_model.set_hop(0.3667, 12, 23, [0, 0, 0])
my_model.set_hop(0.3668, 13, 16, [0, 0, 0])
my_model.set_hop(0.3784, 13, 18, [0, 0, 0])
my_model.set_hop(0.3665, 13, 19, [0, 0, 0])
my_model.set_hop(0.3669, 13, 20, [0, 0, 0])
my_model.set_hop(0.3667, 13, 22, [0, 0, 0])
my_model.set_hop(0.3780, 13, 23, [0, 0, 0])
my_model.set_hop(0.2224, 14, 1, [1, 0, 0])
my_model.set_hop(2.6133, 14, 7, [1, 0, 0])
my_model.set_hop(0.3675, 14, 20, [1, 0, 0])
my_model.set_hop(0.3675, 14, 22, [1, 0, 0])
my_model.set_hop(0.3674, 14, 18, [0, 0, 0])
my_model.set_hop(0.3772, 14, 21, [0, 0, 0])
my_model.set_hop(0.2224, 15, 1, [1, 0, 0])
my_model.set_hop(2.6088, 15, 10, [1, 0, 0])
my_model.set_hop(0.3673, 15, 16, [1, 0, 0])
my_model.set_hop(0.3673, 15, 19, [1, 0, 0])
my_model.set_hop(0.3774, 15, 17, [0, 0, 0])
my_model.set_hop(0.3677, 15, 23, [0, 0, 0])
my_model.set_hop(2.6109, 16, 6, [0, 0, 1])
my_model.set_hop(0.2224, 16, 9, [0, 0, 1])
my_model.set_hop(0.3675, 16, 21, [0, 0, 1])
my_model.set_hop(0.3668, 16, 22, [0, 0, 1])
my_model.set_hop(0.3760, 16, 20, [0, 0, 0])
my_model.set_hop(2.6091, 17, 3, [0, 0, 1])
my_model.set_hop(0.2224, 17, 9, [0, 0, 1])
my_model.set_hop(0.3665, 17, 12, [0, 0, 1])
my_model.set_hop(0.3673, 17, 14, [0, 0, 1])
my_model.set_hop(0.3673, 17, 23, [0, 0, 0])
my_model.set_hop(0.3666, 18, 21, [0, 0, 0])
my_model.set_hop(0.3670, 18, 22, [0, 0, 0])
my_model.set_hop(0.3760, 18, 23, [0, 0, 0])
my_model.set_hop(0.3672, 19, 23, [0, 0, 0])
my_model.set_hop(0.2224, 20, 5, [0, 1, 0])
my_model.set_hop(2.5987, 20, 11, [0, 1, 0])
my_model.set_hop(0.3673, 20, 17, [0, 1, 0])
my_model.set_hop(0.3670, 20, 19, [0, 1, 0])
my_model.set_hop(2.6102, 21, 2, [0, 1, 0])
my_model.set_hop(0.2224, 21, 5, [0, 1, 0])
my_model.set_hop(0.3667, 21, 12, [0, 1, 0])
my_model.set_hop(0.3674, 21, 15, [0, 1, 0])

my_model.display()

path = [ [0., 0., 0.], [0., 0.5, 0.], [0.25, .25, 0.], [.25, .25, .25], [0., 0., 0.], [0.25, 0.25, 0.] ]
label = (r'$\Gamma $',r'$H$', r'$N$', r'$P$', r'$\Gamma $', r'$N$')
nk = 800

(k_vec, k_dist, k_node) = my_model.k_path(path, nk)

print('-'*20)
print('starting calculation')
print('-'*20)
print('Calculating bands . . .')

out = open(u'_bands.txt', 'w')
nodes = open(u'_nodes.txt', 'w')
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

for xband in range(24):
	for k in range(len(k_dist)):
		out.write(str(k_dist[k]) + '\t' + str(evals[xband][k]) + '\n')
	ax.plot(k_dist, evals[xband], linewidth=.5, color='black')
	out.write('\n')

# fig.tight_layout()
for n in range(len(k_node)):
	nodes.write(label[n] + '\t' + str(k_node[n]) + '\n')
fig.savefig('u.eps')

print('Done.')
