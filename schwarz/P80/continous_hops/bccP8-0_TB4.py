# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[6.7549982071, 0.0000000000, 0.0000000000], \
	[-2.2516662123, 6.3686733348, 0.0000000000], \
	[-2.2516662123, -3.1843369714, 5.5154327209]]

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
t1 = -2.958
t2 = -0.09
t3 = -.3

my_model.set_onsite([delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta])

# set hopping parameters for connected orbitals
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(-0.0654, 0, 2, [0, 0, 0])
my_model.set_hop(-0.0654, 0, 3, [0, 0, 0])
my_model.set_hop(-0.3308, 0, 4, [0, 0, 0])
my_model.set_hop(-0.3154, 0, 5, [0, 0, 0])
my_model.set_hop(-0.3308, 0, 8, [0, 0, 0])
my_model.set_hop(-0.3155, 0, 9, [0, 0, 0])
my_model.set_hop(-3.2505, 0, 12, [0, 0, 0])
my_model.set_hop(-0.1431, 0, 13, [0, 0, 0])
my_model.set_hop(-0.1215, 0, 14, [0, 0, 0])
my_model.set_hop(-0.1215, 0, 15, [0, 0, 0])
my_model.set_hop(-2.6133, 0, 18, [0, 0, 0])
my_model.set_hop(-0.0339, 0, 19, [0, 0, 0])
my_model.set_hop(-0.0339, 0, 22, [0, 0, 0])
my_model.set_hop(-2.6077, 0, 23, [0, 0, 0])
my_model.set_hop(-0.3157, 1, 4, [0, 0, 0])
my_model.set_hop(-0.3154, 1, 7, [0, 0, 0])
my_model.set_hop(-0.3157, 1, 8, [0, 0, 0])
my_model.set_hop(-0.3153, 1, 10, [0, 0, 0])
my_model.set_hop(-3.2505, 1, 13, [0, 0, 0])
my_model.set_hop(-2.6098, 1, 16, [0, 0, 0])
my_model.set_hop(-0.0336, 1, 18, [0, 0, 0])
my_model.set_hop(-0.0725, 1, 19, [0, 0, 0])
my_model.set_hop(-2.6106, 1, 20, [0, 0, 0])
my_model.set_hop(-0.0725, 1, 22, [0, 0, 0])
my_model.set_hop(-0.0335, 1, 23, [0, 0, 0])
my_model.set_hop(-0.3168, 2, 5, [0, 0, 0])
my_model.set_hop(-0.1216, 2, 12, [0, 0, 0])
my_model.set_hop(-3.2471, 2, 15, [0, 0, 0])
my_model.set_hop(-0.0335, 2, 17, [0, 0, 0])
my_model.set_hop(-0.0728, 2, 23, [0, 0, 0])
my_model.set_hop(-0.0654, 2, 1, [1, 0, 0])
my_model.set_hop(-0.3162, 2, 8, [1, 0, 0])
my_model.set_hop(-0.3166, 2, 10, [1, 0, 0])
my_model.set_hop(-0.1215, 2, 13, [1, 0, 0])
my_model.set_hop(-0.0727, 2, 16, [1, 0, 0])
my_model.set_hop(-2.6225, 2, 19, [1, 0, 0])
my_model.set_hop(-0.3168, 3, 9, [0, 0, 0])
my_model.set_hop(-0.1215, 3, 12, [0, 0, 0])
my_model.set_hop(-3.2471, 3, 14, [0, 0, 0])
my_model.set_hop(-0.0727, 3, 18, [0, 0, 0])
my_model.set_hop(-0.0335, 3, 21, [0, 0, 0])
my_model.set_hop(-0.0654, 3, 1, [1, 0, 0])
my_model.set_hop(-0.3164, 3, 4, [1, 0, 0])
my_model.set_hop(-0.3168, 3, 7, [1, 0, 0])
my_model.set_hop(-0.1216, 3, 13, [1, 0, 0])
my_model.set_hop(-0.0727, 3, 20, [1, 0, 0])
my_model.set_hop(-2.6253, 3, 22, [1, 0, 0])
my_model.set_hop(-0.0660, 4, 6, [0, 0, 0])
my_model.set_hop(-0.0653, 4, 7, [0, 0, 0])
my_model.set_hop(-0.3278, 4, 8, [0, 0, 0])
my_model.set_hop(-0.3163, 4, 9, [0, 0, 0])
my_model.set_hop(-0.0339, 4, 12, [0, 0, 0])
my_model.set_hop(-2.6106, 4, 13, [0, 0, 0])
my_model.set_hop(-2.6260, 4, 18, [0, 0, 0])
my_model.set_hop(-0.0337, 4, 19, [0, 0, 0])
my_model.set_hop(-0.1215, 4, 20, [0, 0, 0])
my_model.set_hop(-0.1215, 4, 21, [0, 0, 0])
my_model.set_hop(-3.2448, 4, 22, [0, 0, 0])
my_model.set_hop(-0.1417, 4, 23, [0, 0, 0])
my_model.set_hop(-0.3164, 5, 8, [0, 0, 0])
my_model.set_hop(-0.3186, 5, 11, [0, 0, 0])
my_model.set_hop(-0.0725, 5, 12, [0, 0, 0])
my_model.set_hop(-0.0336, 5, 13, [0, 0, 0])
my_model.set_hop(-2.6268, 5, 15, [0, 0, 0])
my_model.set_hop(-2.6114, 5, 17, [0, 0, 0])
my_model.set_hop(-0.0334, 5, 18, [0, 0, 0])
my_model.set_hop(-0.0721, 5, 19, [0, 0, 0])
my_model.set_hop(-3.2496, 5, 23, [0, 0, 0])
my_model.set_hop(-0.3162, 6, 9, [0, 0, 0])
my_model.set_hop(-0.0334, 6, 14, [0, 0, 0])
my_model.set_hop(-0.0731, 6, 18, [0, 0, 0])
my_model.set_hop(-3.2515, 6, 21, [0, 0, 0])
my_model.set_hop(-0.1215, 6, 22, [0, 0, 0])
my_model.set_hop(-0.3155, 6, 0, [0, 1, 0])
my_model.set_hop(-0.3170, 6, 2, [0, 1, 0])
my_model.set_hop(-0.0654, 6, 5, [0, 1, 0])
my_model.set_hop(-2.6077, 6, 12, [0, 1, 0])
my_model.set_hop(-0.0727, 6, 15, [0, 1, 0])
my_model.set_hop(-0.1215, 6, 23, [0, 1, 0])
my_model.set_hop(-0.0725, 7, 13, [0, 0, 0])
my_model.set_hop(-0.0334, 7, 16, [0, 0, 0])
my_model.set_hop(-3.2428, 7, 20, [0, 0, 0])
my_model.set_hop(-0.1215, 7, 22, [0, 0, 0])
my_model.set_hop(-0.0648, 7, 5, [0, 1, 0])
my_model.set_hop(-0.3157, 7, 8, [0, 1, 0])
my_model.set_hop(-0.3133, 7, 11, [0, 1, 0])
my_model.set_hop(-0.0721, 7, 17, [0, 1, 0])
my_model.set_hop(-2.6129, 7, 19, [0, 1, 0])
my_model.set_hop(-0.1215, 7, 23, [0, 1, 0])
my_model.set_hop(-0.0653, 8, 10, [0, 0, 0])
my_model.set_hop(-0.0655, 8, 11, [0, 0, 0])
my_model.set_hop(-0.0339, 8, 12, [0, 0, 0])
my_model.set_hop(-2.6098, 8, 13, [0, 0, 0])
my_model.set_hop(-0.1215, 8, 16, [0, 0, 0])
my_model.set_hop(-0.1216, 8, 17, [0, 0, 0])
my_model.set_hop(-0.1420, 8, 18, [0, 0, 0])
my_model.set_hop(-3.2436, 8, 19, [0, 0, 0])
my_model.set_hop(-0.0337, 8, 22, [0, 0, 0])
my_model.set_hop(-2.6260, 8, 23, [0, 0, 0])
my_model.set_hop(-0.0725, 9, 12, [0, 0, 0])
my_model.set_hop(-0.0336, 9, 13, [0, 0, 0])
my_model.set_hop(-2.6258, 9, 14, [0, 0, 0])
my_model.set_hop(-3.2438, 9, 18, [0, 0, 0])
my_model.set_hop(-2.6100, 9, 21, [0, 0, 0])
my_model.set_hop(-0.0721, 9, 22, [0, 0, 0])
my_model.set_hop(-0.0334, 9, 23, [0, 0, 0])
my_model.set_hop(-0.0725, 10, 13, [0, 0, 0])
my_model.set_hop(-3.2438, 10, 16, [0, 0, 0])
my_model.set_hop(-0.1215, 10, 19, [0, 0, 0])
my_model.set_hop(-0.0334, 10, 20, [0, 0, 0])
my_model.set_hop(-0.3162, 10, 4, [0, 0, 1])
my_model.set_hop(-0.3163, 10, 6, [0, 0, 1])
my_model.set_hop(-0.0649, 10, 9, [0, 0, 1])
my_model.set_hop(-0.1216, 10, 18, [0, 0, 1])
my_model.set_hop(-0.0722, 10, 21, [0, 0, 1])
my_model.set_hop(-2.6129, 10, 22, [0, 0, 1])
my_model.set_hop(-0.0339, 11, 15, [0, 0, 0])
my_model.set_hop(-3.2779, 11, 17, [0, 0, 0])
my_model.set_hop(-0.1203, 11, 19, [0, 0, 0])
my_model.set_hop(-0.0733, 11, 23, [0, 0, 0])
my_model.set_hop(-0.3133, 11, 0, [0, 0, 1])
my_model.set_hop(-0.3196, 11, 3, [0, 0, 1])
my_model.set_hop(-0.0650, 11, 9, [0, 0, 1])
my_model.set_hop(-2.6009, 11, 12, [0, 0, 1])
my_model.set_hop(-0.0729, 11, 14, [0, 0, 1])
my_model.set_hop(-0.1203, 11, 18, [0, 0, 1])
my_model.set_hop(-0.0647, 12, 14, [0, 0, 0])
my_model.set_hop(-0.0647, 12, 15, [0, 0, 0])
my_model.set_hop(-0.3163, 12, 18, [0, 0, 0])
my_model.set_hop(-0.3160, 12, 23, [0, 0, 0])
my_model.set_hop(-0.3161, 13, 16, [0, 0, 0])
my_model.set_hop(-0.3308, 13, 18, [0, 0, 0])
my_model.set_hop(-0.3157, 13, 19, [0, 0, 0])
my_model.set_hop(-0.3162, 13, 20, [0, 0, 0])
my_model.set_hop(-0.3160, 13, 22, [0, 0, 0])
my_model.set_hop(-0.3303, 13, 23, [0, 0, 0])
my_model.set_hop(-0.3169, 14, 18, [0, 0, 0])
my_model.set_hop(-0.3294, 14, 21, [0, 0, 0])
my_model.set_hop(-0.1216, 14, 1, [1, 0, 0])
my_model.set_hop(-0.0727, 14, 4, [1, 0, 0])
my_model.set_hop(-2.6282, 14, 7, [1, 0, 0])
my_model.set_hop(-0.0648, 14, 13, [1, 0, 0])
my_model.set_hop(-0.3171, 14, 20, [1, 0, 0])
my_model.set_hop(-0.3170, 14, 22, [1, 0, 0])
my_model.set_hop(-0.3296, 15, 17, [0, 0, 0])
my_model.set_hop(-0.3173, 15, 23, [0, 0, 0])
my_model.set_hop(-0.1215, 15, 1, [1, 0, 0])
my_model.set_hop(-0.0726, 15, 8, [1, 0, 0])
my_model.set_hop(-2.6240, 15, 10, [1, 0, 0])
my_model.set_hop(-0.0647, 15, 13, [1, 0, 0])
my_model.set_hop(-0.3168, 15, 16, [1, 0, 0])
my_model.set_hop(-0.3168, 15, 19, [1, 0, 0])
my_model.set_hop(-0.0647, 16, 19, [0, 0, 0])
my_model.set_hop(-0.3278, 16, 20, [0, 0, 0])
my_model.set_hop(-0.0731, 16, 4, [0, 0, 1])
my_model.set_hop(-2.6260, 16, 6, [0, 0, 1])
my_model.set_hop(-0.1216, 16, 9, [0, 0, 1])
my_model.set_hop(-0.0652, 16, 18, [0, 0, 1])
my_model.set_hop(-0.3170, 16, 21, [0, 0, 1])
my_model.set_hop(-0.3161, 16, 22, [0, 0, 1])
my_model.set_hop(-0.0642, 17, 19, [0, 0, 0])
my_model.set_hop(-0.3167, 17, 23, [0, 0, 0])
my_model.set_hop(-0.0725, 17, 0, [0, 0, 1])
my_model.set_hop(-2.6243, 17, 3, [0, 0, 1])
my_model.set_hop(-0.1215, 17, 9, [0, 0, 1])
my_model.set_hop(-0.3158, 17, 12, [0, 0, 1])
my_model.set_hop(-0.3168, 17, 14, [0, 0, 1])
my_model.set_hop(-0.0647, 17, 18, [0, 0, 1])
my_model.set_hop(-0.3159, 18, 21, [0, 0, 0])
my_model.set_hop(-0.3164, 18, 22, [0, 0, 0])
my_model.set_hop(-0.3278, 18, 23, [0, 0, 0])
my_model.set_hop(-0.3167, 19, 23, [0, 0, 0])
my_model.set_hop(-0.0647, 20, 22, [0, 0, 0])
my_model.set_hop(-0.1215, 20, 5, [0, 1, 0])
my_model.set_hop(-0.0731, 20, 8, [0, 1, 0])
my_model.set_hop(-2.6146, 20, 11, [0, 1, 0])
my_model.set_hop(-0.3168, 20, 17, [0, 1, 0])
my_model.set_hop(-0.3165, 20, 19, [0, 1, 0])
my_model.set_hop(-0.0653, 20, 23, [0, 1, 0])
my_model.set_hop(-0.0641, 21, 22, [0, 0, 0])
my_model.set_hop(-0.0725, 21, 0, [0, 1, 0])
my_model.set_hop(-2.6253, 21, 2, [0, 1, 0])
my_model.set_hop(-0.1215, 21, 5, [0, 1, 0])
my_model.set_hop(-0.3160, 21, 12, [0, 1, 0])
my_model.set_hop(-0.3169, 21, 15, [0, 1, 0])
my_model.set_hop(-0.0647, 21, 23, [0, 1, 0])

my_model.display()

path = [ [0., 0., 0.], [0.5, -0.5, 0.5], [0., 0., 0.5], [.25, .25, .25], [0., 0., 0.], [0., 0., 0.5] ]
label = (r'$\Gamma $',r'$H$', r'$N$', r'$P$', r'$\Gamma $', r'$N$')
nk = 800

(k_vec, k_dist, k_node) = my_model.k_path(path, nk)

print('-'*20)
print('starting calculation')
print('-'*20)
print('Calculating bands . . .')

out = open('bccP8-0_bands.txt', 'w')
nodes = open('bccP8-0_nodes.txt', 'w')
evals = my_model.solve_all(k_vec)

fig, ax = plt.subplots(figsize=(3,4))

# ax.set_xlim([0, k_node[-1]])
ax.set_xticks(k_node)
ax.set_xticklabels(label)

for n in range(len(k_node)):
	ax.axvline(x=k_node[n], linewidth=0.3, color='black')

ax.set_title('bccP8-0 band structure')
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
fig.savefig('bccP8-0.eps')

print('Done.')