# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[7.0830001831, 0.0000000000, 0.0000000000], \
	[-3.5415000916, 6.1340580936, 0.0000000000], \
	[0.0000000000, 0.0000000000, 5.4200000763]]

orb = [[0.000000000, 0.000000000, 0.000000000], \
	[0.179930001, 0.431919992, 0.000000000], \
	[0.568080008, 0.748010039, 0.000000000], \
	[0.251989990, 0.820070028, 0.000000000], \
	[0.431919992, 0.179930001, 0.000000000], \
	[0.748010039, 0.568080008, 0.000000000], \
	[0.820070028, 0.251989990, 0.000000000], \
	[0.111520000, 0.590279996, 0.000000000], \
	[0.409720004, 0.521239996, 0.000000000], \
	[0.478760004, 0.888480008, 0.000000000], \
	[0.590279996, 0.111520000, 0.000000000], \
	[0.521239996, 0.409720004, 0.000000000], \
	[0.888480008, 0.478760004, 0.000000000], \
	[0.000000000, 0.209619999, 0.000000000], \
	[0.790380001, 0.790380001, 0.000000000], \
	[0.209619999, 0.000000000, 0.000000000]]

my_model = tb_model(3, 3, lat, orb)

# set model parameters
delta = 0.0
t1 = -1.
#t2 = -0.09
#t3 = -.3
t4 = 0.0
t5 = 0.0

my_model.set_onsite([delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta])

# set hopping parameters for connected orbitals
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(t1, 0, 13, [0, 0, 0])
my_model.set_hop(t1, 0, 15, [0, 0, 0])
my_model.set_hop(t1, 1, 7, [0, 0, 0])
my_model.set_hop(t1, 1, 8, [0, 0, 0])
my_model.set_hop(t1, 1, 13, [0, 0, 0])
my_model.set_hop(t1, 2, 8, [0, 0, 0])
my_model.set_hop(t1, 2, 9, [0, 0, 0])
my_model.set_hop(t1, 2, 14, [0, 0, 0])
my_model.set_hop(t1, 3, 7, [0, 0, 0])
my_model.set_hop(t1, 3, 9, [0, 0, 0])
my_model.set_hop(t1, 3, 15, [0, 1, 0])
my_model.set_hop(t1, 4, 10, [0, 0, 0])
my_model.set_hop(t1, 4, 11, [0, 0, 0])
my_model.set_hop(t1, 4, 15, [0, 0, 0])
my_model.set_hop(t1, 5, 11, [0, 0, 0])
my_model.set_hop(t1, 5, 12, [0, 0, 0])
my_model.set_hop(t1, 5, 14, [0, 0, 0])
my_model.set_hop(t1, 6, 10, [0, 0, 0])
my_model.set_hop(t1, 6, 12, [0, 0, 0])
my_model.set_hop(t1, 6, 13, [1, 0, 0])
my_model.set_hop(t1, 8, 11, [0, 0, 0])
my_model.set_hop(t1, 9, 10, [0, 1, 0])
my_model.set_hop(t1, 12, 7, [1, 0, 0])
my_model.set_hop(t1, 14, 0, [1, 1, 0])

my_model.display()

path = [ [0., 0., 0.], [1./3., 1./3., 0.], [0.5, 0., 0.], [0., 0., 0.] ]
label = (r'$\Gamma $',r'$K$', r'$M$', r'$\Gamma$')
nk = 1000

(k_vec, k_dist, k_node) = my_model.k_path(path, nk)

print('-'*20)
print('starting calculation')
print('-'*20)
print('Calculating bands . . .')

out = open('haeck_bands.txt', 'w')
nodes = open('haeck_nodes.txt', 'w')
evals = my_model.solve_all(k_vec)

fig, ax = plt.subplots(figsize=(4,5))

# ax.set_xlim([0, k_node[-1]])
ax.set_xticks(k_node)
ax.set_xticklabels(label)

for n in range(len(k_node)):
	ax.axvline(x=k_node[n], linewidth=0.3, color='black')

ax.set_title('haeck band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(16):
	for k in range(len(k_dist)):
		out.write(str(k_dist[k]) + '\t' + str(evals[xband][k]) + '\n')
	ax.plot(k_dist, evals[xband], linewidth=.5, color='black')
	out.write('\n')

# fig.tight_layout()
for n in range(len(k_node)):
	nodes.write(label[n] + '\t' + str(k_node[n]) + '\n')
fig.savefig('haeck.eps')

k_vec = my_model.k_uniform_mesh([7, 7, 7])
evals2 = my_model.solve_all(k_vec)

evals2 = evals2.flatten()

print('Plotting DoS...')

fig, ax = plt.subplots()
ax.hist(evals2, 80,range=(-4.,4.))
#ax.set_ylim(0.0,80.0)
# put title
ax.set_title('haeck density of states')
ax.set_xlabel('Band energy')
ax.set_ylabel('Number of states')
# make an PDF figure of a plot
fig.tight_layout()
fig.savefig('haeck_dos.pdf')
print('Done.')