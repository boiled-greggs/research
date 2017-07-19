# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[8.3147096634, 0.0000000000, 0.0000000000], \
	[-2.7715700642, 7.8391833864, 0.0000000000], \
	[-2.7715700642, -3.9195920674, 6.7889317415]]

orb = [[0.042018708, 0.074983709, 0.876981139], \
	[0.467079133, 0.345995039, 0.031030415], \
	[0.696067750, 0.183256760, 0.639303327], \
	[0.813029945, 0.454268068, 0.485254139], \
	[0.873387754, 0.052069660, 0.068526022], \
	[0.027489064, 0.477182120, 0.339537352], \
	[0.635709941, 0.706170857, 0.176747084], \
	[0.481608719, 0.823080838, 0.447758347], \
	[0.064880721, 0.883438766, 0.045559939], \
	[0.335944146, 0.037540015, 0.470724553], \
	[0.173153803, 0.645813048, 0.699661076], \
	[0.444217145, 0.491711766, 0.816623271], \
	[0.042018767, 0.491711795, 0.176743194], \
	[0.467079133, 0.645757020, 0.447758406], \
	[0.696063936, 0.037540015, 0.068526201], \
	[0.813033879, 0.883494735, 0.339541316], \
	[0.873391747, 0.346051067, 0.816627264], \
	[0.027489064, 0.074983768, 0.699657261], \
	[0.635706067, 0.454268068, 0.470724612], \
	[0.481608778, 0.183200821, 0.045559999], \
	[0.064880721, 0.706114829, 0.031030295], \
	[0.335948020, 0.823136866, 0.876985013], \
	[0.173149809, 0.052069601, 0.485254139], \
	[0.444217145, 0.477182209, 0.639299452], \
	[0.967079103, 0.954268038, 0.139303342], \
	[0.542018592, 0.683256686, 0.985254169], \
	[0.313029975, 0.845995009, 0.376981169], \
	[0.196067825, 0.574983656, 0.531030357], \
	[0.135710001, 0.977182090, 0.947758317], \
	[0.981608689, 0.552069604, 0.676747024], \
	[0.373387724, 0.323080927, 0.839537323], \
	[0.527489007, 0.206170887, 0.568526089], \
	[0.944217145, 0.145813063, 0.970724523], \
	[0.673153698, 0.991711795, 0.545559883], \
	[0.835944116, 0.383438736, 0.316623300], \
	[0.564880669, 0.537540019, 0.199661151], \
	[0.967079103, 0.537540019, 0.839541256], \
	[0.542018652, 0.383494765, 0.568526089], \
	[0.313033879, 0.991711795, 0.947758377], \
	[0.196063980, 0.145757064, 0.676743209], \
	[0.135706052, 0.683200777, 0.199657306], \
	[0.981608748, 0.954268038, 0.316627294], \
	[0.373391777, 0.574983716, 0.545559943], \
	[0.527489066, 0.846051037, 0.970724583], \
	[0.944217145, 0.323136896, 0.985254109], \
	[0.673149765, 0.206114873, 0.139299497], \
	[0.835947990, 0.977182209, 0.531030238], \
	[0.564880610, 0.552069545, 0.376985043]]

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
my_model.set_hop(t1, 0, 17, [0, 0, 0])
my_model.set_hop(t2, 1, 6, [0, 0, 0])
my_model.set_hop(t1, 1, 19, [0, 0, 0])
my_model.set_hop(t1, 1, 35, [0, 0, 0])
my_model.set_hop(t2, 1, 45, [0, 0, 0])
my_model.set_hop(t2, 1, 47, [0, 0, 0])
my_model.set_hop(t2, 2, 9, [0, 0, 0])
my_model.set_hop(t1, 2, 16, [0, 0, 0])
my_model.set_hop(t1, 2, 31, [0, 0, 0])
my_model.set_hop(t2, 2, 44, [0, 0, 0])
my_model.set_hop(t1, 3, 18, [0, 0, 0])
my_model.set_hop(t1, 3, 29, [0, 0, 0])
my_model.set_hop(t1, 3, 34, [0, 0, 0])
my_model.set_hop(t2, 3, 45, [0, 0, 0])
my_model.set_hop(t2, 3, 47, [0, 0, 0])
my_model.set_hop(t1, 4, 14, [0, 0, 0])
my_model.set_hop(t2, 5, 10, [0, 0, 0])
my_model.set_hop(t1, 5, 12, [0, 0, 0])
my_model.set_hop(t1, 5, 27, [0, 0, 0])
my_model.set_hop(t2, 5, 40, [0, 0, 0])
my_model.set_hop(t2, 5, 42, [0, 0, 0])
my_model.set_hop(t1, 6, 15, [0, 0, 0])
my_model.set_hop(t1, 6, 35, [0, 0, 0])
my_model.set_hop(t2, 6, 41, [0, 0, 0])
my_model.set_hop(t1, 7, 13, [0, 0, 0])
my_model.set_hop(t1, 7, 26, [0, 0, 0])
my_model.set_hop(t1, 7, 33, [0, 0, 0])
my_model.set_hop(t2, 7, 40, [0, 0, 0])
my_model.set_hop(t2, 7, 42, [0, 0, 0])
my_model.set_hop(t1, 8, 20, [0, 0, 0])
my_model.set_hop(t1, 9, 22, [0, 0, 0])
my_model.set_hop(t1, 9, 31, [0, 0, 0])
my_model.set_hop(t2, 9, 37, [0, 0, 0])
my_model.set_hop(t2, 9, 39, [0, 0, 0])
my_model.set_hop(t1, 10, 21, [0, 0, 0])
my_model.set_hop(t1, 10, 27, [0, 0, 0])
my_model.set_hop(t2, 10, 38, [0, 0, 0])
my_model.set_hop(t1, 11, 23, [0, 0, 0])
my_model.set_hop(t1, 11, 25, [0, 0, 0])
my_model.set_hop(t1, 11, 30, [0, 0, 0])
my_model.set_hop(t2, 11, 37, [0, 0, 0])
my_model.set_hop(t2, 11, 39, [0, 0, 0])
my_model.set_hop(t2, 12, 20, [0, 0, 0])
my_model.set_hop(t1, 12, 40, [0, 0, 0])
my_model.set_hop(t2, 13, 18, [0, 0, 0])
my_model.set_hop(t2, 13, 23, [0, 0, 0])
my_model.set_hop(t2, 13, 33, [0, 0, 0])
my_model.set_hop(t2, 13, 35, [0, 0, 0])
my_model.set_hop(t1, 13, 42, [0, 0, 0])
my_model.set_hop(t1, 13, 47, [0, 0, 0])
my_model.set_hop(t2, 14, 19, [0, 0, 0])
my_model.set_hop(t2, 14, 34, [0, 0, 0])
my_model.set_hop(t1, 14, 45, [0, 0, 0])
my_model.set_hop(t2, 15, 33, [0, 0, 0])
my_model.set_hop(t2, 15, 35, [0, 0, 0])
my_model.set_hop(t1, 15, 41, [0, 0, 0])
my_model.set_hop(t1, 15, 46, [0, 0, 0])
my_model.set_hop(t2, 16, 29, [0, 0, 0])
my_model.set_hop(t2, 16, 31, [0, 0, 0])
my_model.set_hop(t1, 16, 36, [0, 0, 0])
my_model.set_hop(t1, 16, 44, [0, 0, 0])
my_model.set_hop(t2, 17, 22, [0, 0, 0])
my_model.set_hop(t2, 17, 30, [0, 0, 0])
my_model.set_hop(t1, 17, 39, [0, 0, 0])
my_model.set_hop(t2, 18, 23, [0, 0, 0])
my_model.set_hop(t2, 18, 29, [0, 0, 0])
my_model.set_hop(t2, 18, 31, [0, 0, 0])
my_model.set_hop(t1, 18, 37, [0, 0, 0])
my_model.set_hop(t1, 18, 47, [0, 0, 0])
my_model.set_hop(t1, 19, 45, [0, 0, 0])
my_model.set_hop(t2, 20, 26, [0, 0, 0])
my_model.set_hop(t1, 20, 40, [0, 0, 0])
my_model.set_hop(t2, 21, 25, [0, 0, 0])
my_model.set_hop(t2, 21, 27, [0, 0, 0])
my_model.set_hop(t1, 21, 38, [0, 0, 0])
my_model.set_hop(t1, 21, 43, [0, 0, 0])
my_model.set_hop(t1, 22, 39, [0, 0, 0])
my_model.set_hop(t2, 23, 25, [0, 0, 0])
my_model.set_hop(t2, 23, 27, [0, 0, 0])
my_model.set_hop(t1, 23, 37, [0, 0, 0])
my_model.set_hop(t1, 23, 42, [0, 0, 0])
my_model.set_hop(t1, 24, 8, [1, 0, 0])
my_model.set_hop(t2, 24, 20, [1, 0, 0])
my_model.set_hop(t1, 24, 41, [0, 0, 0])
my_model.set_hop(t1, 25, 6, [0, 0, 1])
my_model.set_hop(t2, 25, 35, [0, 0, 1])
my_model.set_hop(t2, 25, 30, [0, 0, 0])
my_model.set_hop(t1, 25, 43, [0, 0, 0])
my_model.set_hop(t2, 26, 33, [0, 0, 0])
my_model.set_hop(t1, 26, 40, [0, 0, 0])
my_model.set_hop(t1, 27, 42, [0, 0, 0])
my_model.set_hop(t1, 28, 0, [0, 1, 0])
my_model.set_hop(t2, 28, 17, [0, 1, 0])
my_model.set_hop(t1, 28, 38, [0, 0, 0])
my_model.set_hop(t1, 29, 10, [1, 0, 0])
my_model.set_hop(t2, 29, 27, [1, 0, 0])
my_model.set_hop(t2, 29, 34, [0, 0, 0])
my_model.set_hop(t1, 29, 36, [0, 0, 0])
my_model.set_hop(t1, 30, 39, [0, 0, 0])
my_model.set_hop(t1, 31, 37, [0, 0, 0])
my_model.set_hop(t1, 32, 0, [1, 0, 0])
my_model.set_hop(t1, 32, 44, [0, 0, 0])
my_model.set_hop(t1, 33, 2, [0, 1, 0])
my_model.set_hop(t2, 33, 31, [0, 1, 0])
my_model.set_hop(t1, 33, 46, [0, 0, 0])
my_model.set_hop(t1, 34, 45, [0, 0, 0])
my_model.set_hop(t1, 35, 47, [0, 0, 0])
my_model.set_hop(t2, 36, 10, [1, 0, 0])
my_model.set_hop(t2, 36, 44, [0, 0, 0])
my_model.set_hop(t2, 37, 42, [0, 0, 0])
my_model.set_hop(t2, 37, 47, [0, 0, 0])
my_model.set_hop(t2, 38, 43, [0, 0, 0])
my_model.set_hop(t2, 41, 46, [0, 0, 0])
my_model.set_hop(t2, 42, 47, [0, 0, 0])
my_model.set_hop(t2, 43, 6, [0, 0, 1])
my_model.set_hop(t2, 44, 0, [1, 0, 0])
my_model.set_hop(t2, 46, 2, [0, 1, 0])

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

ax.set_title('bccG8-0 schwarzite band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(48):
	ax.plot(k_dist, evals[xband])

# fig.tight_layout()
fig.savefig('bccG8-0.pdf')

print('Done.')