# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import plotly
import plotly.plotly as plt
import plotly.graph_objs as go

lat = [[11.2060003281, 0.0000000000, 0.0000000000], \
	[0.0000000000, 11.2060003281, 0.0000000000], \
	[0.0000000000, 0.0000000000, 4.2670001984]]

orb = [[0.850627899, 0.500000000, 0.832919121], \
	[0.850627899, 0.500000000, 0.499929070], \
	[0.833493590, 0.608338594, 0.332990050], \
	[0.833493590, 0.608338594, 0.000000000], \
	[0.783697128, 0.706146717, 0.832919121], \
	[0.783697128, 0.706146717, 0.499929070], \
	[0.706146717, 0.783697128, 0.332990050], \
	[0.706146717, 0.783697128, 0.000000000], \
	[0.608338594, 0.833493590, 0.832919121], \
	[0.608338594, 0.833493590, 0.499929070], \
	[0.500000000, 0.850627899, 0.332990050], \
	[0.500000000, 0.850627899, 0.000000000], \
	[0.391661435, 0.833493590, 0.832919121], \
	[0.391661435, 0.833493590, 0.499929070], \
	[0.293853283, 0.783697128, 0.332990050], \
	[0.293853283, 0.783697128, 0.000000000], \
	[0.216302872, 0.706146717, 0.832919121], \
	[0.216302872, 0.706146717, 0.499929070], \
	[0.166506380, 0.608338594, 0.332990050], \
	[0.166506380, 0.608338594, 0.000000000], \
	[0.149372131, 0.500000000, 0.832919121], \
	[0.149372131, 0.500000000, 0.499929070], \
	[0.166506380, 0.391661435, 0.332990050], \
	[0.166506380, 0.391661435, 0.000000000], \
	[0.216302872, 0.293853283, 0.832919121], \
	[0.216302872, 0.293853283, 0.499929070], \
	[0.293853283, 0.216302872, 0.332990050], \
	[0.293853283, 0.216302872, 0.000000000], \
	[0.391661435, 0.166506380, 0.832919121], \
	[0.391661435, 0.166506380, 0.499929070], \
	[0.500000000, 0.149372131, 0.332990050], \
	[0.500000000, 0.149372131, 0.000000000], \
	[0.608338594, 0.166506380, 0.832919121], \
	[0.608338594, 0.166506380, 0.499929070], \
	[0.706146717, 0.216302872, 0.332990050], \
	[0.706146717, 0.216302872, 0.000000000], \
	[0.783697128, 0.293853283, 0.832919121], \
	[0.783697128, 0.293853283, 0.499929070], \
	[0.833493590, 0.391661435, 0.332990050], \
	[0.833493590, 0.391661435, 0.000000000]]

my_model = tb_model(1, 3, lat, orb, per=[2])

# set model parameters
delta = 0.0
t1 = -2.9
t2 = 0.
t3 = 0.

my_model.set_onsite([delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta])

# set hopping parameters for connected orbitals
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(t3, 0, 1, [0, 0, 1])
my_model.set_hop(t2, 0, 2, [0, 0, 1])
my_model.set_hop(t1, 0, 3, [0, 0, 1])
my_model.set_hop(t3, 0, 5, [0, 0, 1])
my_model.set_hop(t3, 0, 7, [0, 0, 1])
my_model.set_hop(t3, 0, 35, [0, 0, 1])
my_model.set_hop(t3, 0, 37, [0, 0, 1])
my_model.set_hop(t2, 0, 38, [0, 0, 1])
my_model.set_hop(t1, 0, 39, [0, 0, 1])
my_model.set_hop(t1, 0, 1, [0, 0, 0])
my_model.set_hop(t2, 0, 2, [0, 0, 0])
my_model.set_hop(t3, 0, 3, [0, 0, 0])
my_model.set_hop(t2, 0, 4, [0, 0, 0])
my_model.set_hop(t3, 0, 5, [0, 0, 0])
my_model.set_hop(t2, 0, 36, [0, 0, 0])
my_model.set_hop(t3, 0, 37, [0, 0, 0])
my_model.set_hop(t2, 0, 38, [0, 0, 0])
my_model.set_hop(t3, 0, 39, [0, 0, 0])
my_model.set_hop(t1, 1, 2, [0, 0, 0])
my_model.set_hop(t2, 1, 3, [0, 0, 0])
my_model.set_hop(t3, 1, 4, [0, 0, 0])
my_model.set_hop(t2, 1, 5, [0, 0, 0])
my_model.set_hop(t3, 1, 6, [0, 0, 0])
my_model.set_hop(t3, 1, 34, [0, 0, 0])
my_model.set_hop(t3, 1, 36, [0, 0, 0])
my_model.set_hop(t2, 1, 37, [0, 0, 0])
my_model.set_hop(t1, 1, 38, [0, 0, 0])
my_model.set_hop(t2, 1, 39, [0, 0, 0])
my_model.set_hop(t1, 2, 3, [0, 0, 0])
my_model.set_hop(t2, 2, 4, [0, 0, 0])
my_model.set_hop(t1, 2, 5, [0, 0, 0])
my_model.set_hop(t2, 2, 6, [0, 0, 0])
my_model.set_hop(t3, 2, 7, [0, 0, 0])
my_model.set_hop(t3, 2, 9, [0, 0, 0])
my_model.set_hop(t3, 2, 37, [0, 0, 0])
my_model.set_hop(t2, 2, 38, [0, 0, 0])
my_model.set_hop(t3, 2, 39, [0, 0, 0])
my_model.set_hop(t3, 3, 4, [0, 0, 0])
my_model.set_hop(t2, 3, 5, [0, 0, 0])
my_model.set_hop(t3, 3, 6, [0, 0, 0])
my_model.set_hop(t2, 3, 7, [0, 0, 0])
my_model.set_hop(t3, 3, 38, [0, 0, 0])
my_model.set_hop(t2, 3, 39, [0, 0, 0])
my_model.set_hop(t3, 4, 1, [0, 0, 1])
my_model.set_hop(t2, 4, 2, [0, 0, 1])
my_model.set_hop(t1, 4, 3, [0, 0, 1])
my_model.set_hop(t3, 4, 5, [0, 0, 1])
my_model.set_hop(t2, 4, 6, [0, 0, 1])
my_model.set_hop(t1, 4, 7, [0, 0, 1])
my_model.set_hop(t3, 4, 9, [0, 0, 1])
my_model.set_hop(t3, 4, 11, [0, 0, 1])
my_model.set_hop(t3, 4, 39, [0, 0, 1])
my_model.set_hop(t1, 4, 5, [0, 0, 0])
my_model.set_hop(t2, 4, 6, [0, 0, 0])
my_model.set_hop(t3, 4, 7, [0, 0, 0])
my_model.set_hop(t2, 4, 8, [0, 0, 0])
my_model.set_hop(t3, 4, 9, [0, 0, 0])
my_model.set_hop(t1, 5, 6, [0, 0, 0])
my_model.set_hop(t2, 5, 7, [0, 0, 0])
my_model.set_hop(t3, 5, 8, [0, 0, 0])
my_model.set_hop(t2, 5, 9, [0, 0, 0])
my_model.set_hop(t3, 5, 10, [0, 0, 0])
my_model.set_hop(t3, 5, 38, [0, 0, 0])
my_model.set_hop(t1, 6, 7, [0, 0, 0])
my_model.set_hop(t2, 6, 8, [0, 0, 0])
my_model.set_hop(t1, 6, 9, [0, 0, 0])
my_model.set_hop(t2, 6, 10, [0, 0, 0])
my_model.set_hop(t3, 6, 11, [0, 0, 0])
my_model.set_hop(t3, 6, 13, [0, 0, 0])
my_model.set_hop(t3, 7, 8, [0, 0, 0])
my_model.set_hop(t2, 7, 9, [0, 0, 0])
my_model.set_hop(t3, 7, 10, [0, 0, 0])
my_model.set_hop(t2, 7, 11, [0, 0, 0])
my_model.set_hop(t3, 8, 3, [0, 0, 1])
my_model.set_hop(t3, 8, 5, [0, 0, 1])
my_model.set_hop(t2, 8, 6, [0, 0, 1])
my_model.set_hop(t1, 8, 7, [0, 0, 1])
my_model.set_hop(t3, 8, 9, [0, 0, 1])
my_model.set_hop(t2, 8, 10, [0, 0, 1])
my_model.set_hop(t1, 8, 11, [0, 0, 1])
my_model.set_hop(t3, 8, 13, [0, 0, 1])
my_model.set_hop(t3, 8, 15, [0, 0, 1])
my_model.set_hop(t1, 8, 9, [0, 0, 0])
my_model.set_hop(t2, 8, 10, [0, 0, 0])
my_model.set_hop(t3, 8, 11, [0, 0, 0])
my_model.set_hop(t2, 8, 12, [0, 0, 0])
my_model.set_hop(t3, 8, 13, [0, 0, 0])
my_model.set_hop(t1, 9, 10, [0, 0, 0])
my_model.set_hop(t2, 9, 11, [0, 0, 0])
my_model.set_hop(t3, 9, 12, [0, 0, 0])
my_model.set_hop(t2, 9, 13, [0, 0, 0])
my_model.set_hop(t3, 9, 14, [0, 0, 0])
my_model.set_hop(t1, 10, 11, [0, 0, 0])
my_model.set_hop(t2, 10, 12, [0, 0, 0])
my_model.set_hop(t1, 10, 13, [0, 0, 0])
my_model.set_hop(t2, 10, 14, [0, 0, 0])
my_model.set_hop(t3, 10, 15, [0, 0, 0])
my_model.set_hop(t3, 10, 17, [0, 0, 0])
my_model.set_hop(t3, 11, 12, [0, 0, 0])
my_model.set_hop(t2, 11, 13, [0, 0, 0])
my_model.set_hop(t3, 11, 14, [0, 0, 0])
my_model.set_hop(t2, 11, 15, [0, 0, 0])
my_model.set_hop(t3, 12, 7, [0, 0, 1])
my_model.set_hop(t3, 12, 9, [0, 0, 1])
my_model.set_hop(t2, 12, 10, [0, 0, 1])
my_model.set_hop(t1, 12, 11, [0, 0, 1])
my_model.set_hop(t3, 12, 13, [0, 0, 1])
my_model.set_hop(t2, 12, 14, [0, 0, 1])
my_model.set_hop(t1, 12, 15, [0, 0, 1])
my_model.set_hop(t3, 12, 17, [0, 0, 1])
my_model.set_hop(t3, 12, 19, [0, 0, 1])
my_model.set_hop(t1, 12, 13, [0, 0, 0])
my_model.set_hop(t2, 12, 14, [0, 0, 0])
my_model.set_hop(t3, 12, 15, [0, 0, 0])
my_model.set_hop(t2, 12, 16, [0, 0, 0])
my_model.set_hop(t3, 12, 17, [0, 0, 0])
my_model.set_hop(t1, 13, 14, [0, 0, 0])
my_model.set_hop(t2, 13, 15, [0, 0, 0])
my_model.set_hop(t3, 13, 16, [0, 0, 0])
my_model.set_hop(t2, 13, 17, [0, 0, 0])
my_model.set_hop(t3, 13, 18, [0, 0, 0])
my_model.set_hop(t1, 14, 15, [0, 0, 0])
my_model.set_hop(t2, 14, 16, [0, 0, 0])
my_model.set_hop(t1, 14, 17, [0, 0, 0])
my_model.set_hop(t2, 14, 18, [0, 0, 0])
my_model.set_hop(t3, 14, 19, [0, 0, 0])
my_model.set_hop(t3, 14, 21, [0, 0, 0])
my_model.set_hop(t3, 15, 16, [0, 0, 0])
my_model.set_hop(t2, 15, 17, [0, 0, 0])
my_model.set_hop(t3, 15, 18, [0, 0, 0])
my_model.set_hop(t2, 15, 19, [0, 0, 0])
my_model.set_hop(t3, 16, 11, [0, 0, 1])
my_model.set_hop(t3, 16, 13, [0, 0, 1])
my_model.set_hop(t2, 16, 14, [0, 0, 1])
my_model.set_hop(t1, 16, 15, [0, 0, 1])
my_model.set_hop(t3, 16, 17, [0, 0, 1])
my_model.set_hop(t2, 16, 18, [0, 0, 1])
my_model.set_hop(t1, 16, 19, [0, 0, 1])
my_model.set_hop(t3, 16, 21, [0, 0, 1])
my_model.set_hop(t3, 16, 23, [0, 0, 1])
my_model.set_hop(t1, 16, 17, [0, 0, 0])
my_model.set_hop(t2, 16, 18, [0, 0, 0])
my_model.set_hop(t3, 16, 19, [0, 0, 0])
my_model.set_hop(t2, 16, 20, [0, 0, 0])
my_model.set_hop(t3, 16, 21, [0, 0, 0])
my_model.set_hop(t1, 17, 18, [0, 0, 0])
my_model.set_hop(t2, 17, 19, [0, 0, 0])
my_model.set_hop(t3, 17, 20, [0, 0, 0])
my_model.set_hop(t2, 17, 21, [0, 0, 0])
my_model.set_hop(t3, 17, 22, [0, 0, 0])
my_model.set_hop(t1, 18, 19, [0, 0, 0])
my_model.set_hop(t2, 18, 20, [0, 0, 0])
my_model.set_hop(t1, 18, 21, [0, 0, 0])
my_model.set_hop(t2, 18, 22, [0, 0, 0])
my_model.set_hop(t3, 18, 23, [0, 0, 0])
my_model.set_hop(t3, 18, 25, [0, 0, 0])
my_model.set_hop(t3, 19, 20, [0, 0, 0])
my_model.set_hop(t2, 19, 21, [0, 0, 0])
my_model.set_hop(t3, 19, 22, [0, 0, 0])
my_model.set_hop(t2, 19, 23, [0, 0, 0])
my_model.set_hop(t3, 20, 15, [0, 0, 1])
my_model.set_hop(t3, 20, 17, [0, 0, 1])
my_model.set_hop(t2, 20, 18, [0, 0, 1])
my_model.set_hop(t1, 20, 19, [0, 0, 1])
my_model.set_hop(t3, 20, 21, [0, 0, 1])
my_model.set_hop(t2, 20, 22, [0, 0, 1])
my_model.set_hop(t1, 20, 23, [0, 0, 1])
my_model.set_hop(t3, 20, 25, [0, 0, 1])
my_model.set_hop(t3, 20, 27, [0, 0, 1])
my_model.set_hop(t1, 20, 21, [0, 0, 0])
my_model.set_hop(t2, 20, 22, [0, 0, 0])
my_model.set_hop(t3, 20, 23, [0, 0, 0])
my_model.set_hop(t2, 20, 24, [0, 0, 0])
my_model.set_hop(t3, 20, 25, [0, 0, 0])
my_model.set_hop(t1, 21, 22, [0, 0, 0])
my_model.set_hop(t2, 21, 23, [0, 0, 0])
my_model.set_hop(t3, 21, 24, [0, 0, 0])
my_model.set_hop(t2, 21, 25, [0, 0, 0])
my_model.set_hop(t3, 21, 26, [0, 0, 0])
my_model.set_hop(t1, 22, 23, [0, 0, 0])
my_model.set_hop(t2, 22, 24, [0, 0, 0])
my_model.set_hop(t1, 22, 25, [0, 0, 0])
my_model.set_hop(t2, 22, 26, [0, 0, 0])
my_model.set_hop(t3, 22, 27, [0, 0, 0])
my_model.set_hop(t3, 22, 29, [0, 0, 0])
my_model.set_hop(t3, 23, 24, [0, 0, 0])
my_model.set_hop(t2, 23, 25, [0, 0, 0])
my_model.set_hop(t3, 23, 26, [0, 0, 0])
my_model.set_hop(t2, 23, 27, [0, 0, 0])
my_model.set_hop(t3, 24, 19, [0, 0, 1])
my_model.set_hop(t3, 24, 21, [0, 0, 1])
my_model.set_hop(t2, 24, 22, [0, 0, 1])
my_model.set_hop(t1, 24, 23, [0, 0, 1])
my_model.set_hop(t3, 24, 25, [0, 0, 1])
my_model.set_hop(t2, 24, 26, [0, 0, 1])
my_model.set_hop(t1, 24, 27, [0, 0, 1])
my_model.set_hop(t3, 24, 29, [0, 0, 1])
my_model.set_hop(t3, 24, 31, [0, 0, 1])
my_model.set_hop(t1, 24, 25, [0, 0, 0])
my_model.set_hop(t2, 24, 26, [0, 0, 0])
my_model.set_hop(t3, 24, 27, [0, 0, 0])
my_model.set_hop(t2, 24, 28, [0, 0, 0])
my_model.set_hop(t3, 24, 29, [0, 0, 0])
my_model.set_hop(t1, 25, 26, [0, 0, 0])
my_model.set_hop(t2, 25, 27, [0, 0, 0])
my_model.set_hop(t3, 25, 28, [0, 0, 0])
my_model.set_hop(t2, 25, 29, [0, 0, 0])
my_model.set_hop(t3, 25, 30, [0, 0, 0])
my_model.set_hop(t1, 26, 27, [0, 0, 0])
my_model.set_hop(t2, 26, 28, [0, 0, 0])
my_model.set_hop(t1, 26, 29, [0, 0, 0])
my_model.set_hop(t2, 26, 30, [0, 0, 0])
my_model.set_hop(t3, 26, 31, [0, 0, 0])
my_model.set_hop(t3, 26, 33, [0, 0, 0])
my_model.set_hop(t3, 27, 28, [0, 0, 0])
my_model.set_hop(t2, 27, 29, [0, 0, 0])
my_model.set_hop(t3, 27, 30, [0, 0, 0])
my_model.set_hop(t2, 27, 31, [0, 0, 0])
my_model.set_hop(t3, 28, 23, [0, 0, 1])
my_model.set_hop(t3, 28, 25, [0, 0, 1])
my_model.set_hop(t2, 28, 26, [0, 0, 1])
my_model.set_hop(t1, 28, 27, [0, 0, 1])
my_model.set_hop(t3, 28, 29, [0, 0, 1])
my_model.set_hop(t2, 28, 30, [0, 0, 1])
my_model.set_hop(t1, 28, 31, [0, 0, 1])
my_model.set_hop(t3, 28, 33, [0, 0, 1])
my_model.set_hop(t3, 28, 35, [0, 0, 1])
my_model.set_hop(t1, 28, 29, [0, 0, 0])
my_model.set_hop(t2, 28, 30, [0, 0, 0])
my_model.set_hop(t3, 28, 31, [0, 0, 0])
my_model.set_hop(t2, 28, 32, [0, 0, 0])
my_model.set_hop(t3, 28, 33, [0, 0, 0])
my_model.set_hop(t1, 29, 30, [0, 0, 0])
my_model.set_hop(t2, 29, 31, [0, 0, 0])
my_model.set_hop(t3, 29, 32, [0, 0, 0])
my_model.set_hop(t2, 29, 33, [0, 0, 0])
my_model.set_hop(t3, 29, 34, [0, 0, 0])
my_model.set_hop(t1, 30, 31, [0, 0, 0])
my_model.set_hop(t2, 30, 32, [0, 0, 0])
my_model.set_hop(t1, 30, 33, [0, 0, 0])
my_model.set_hop(t2, 30, 34, [0, 0, 0])
my_model.set_hop(t3, 30, 35, [0, 0, 0])
my_model.set_hop(t3, 30, 37, [0, 0, 0])
my_model.set_hop(t3, 31, 32, [0, 0, 0])
my_model.set_hop(t2, 31, 33, [0, 0, 0])
my_model.set_hop(t3, 31, 34, [0, 0, 0])
my_model.set_hop(t2, 31, 35, [0, 0, 0])
my_model.set_hop(t3, 32, 27, [0, 0, 1])
my_model.set_hop(t3, 32, 29, [0, 0, 1])
my_model.set_hop(t2, 32, 30, [0, 0, 1])
my_model.set_hop(t1, 32, 31, [0, 0, 1])
my_model.set_hop(t3, 32, 33, [0, 0, 1])
my_model.set_hop(t2, 32, 34, [0, 0, 1])
my_model.set_hop(t1, 32, 35, [0, 0, 1])
my_model.set_hop(t3, 32, 37, [0, 0, 1])
my_model.set_hop(t3, 32, 39, [0, 0, 1])
my_model.set_hop(t1, 32, 33, [0, 0, 0])
my_model.set_hop(t2, 32, 34, [0, 0, 0])
my_model.set_hop(t3, 32, 35, [0, 0, 0])
my_model.set_hop(t2, 32, 36, [0, 0, 0])
my_model.set_hop(t3, 32, 37, [0, 0, 0])
my_model.set_hop(t1, 33, 34, [0, 0, 0])
my_model.set_hop(t2, 33, 35, [0, 0, 0])
my_model.set_hop(t3, 33, 36, [0, 0, 0])
my_model.set_hop(t2, 33, 37, [0, 0, 0])
my_model.set_hop(t3, 33, 38, [0, 0, 0])
my_model.set_hop(t1, 34, 35, [0, 0, 0])
my_model.set_hop(t2, 34, 36, [0, 0, 0])
my_model.set_hop(t1, 34, 37, [0, 0, 0])
my_model.set_hop(t2, 34, 38, [0, 0, 0])
my_model.set_hop(t3, 34, 39, [0, 0, 0])
my_model.set_hop(t3, 35, 36, [0, 0, 0])
my_model.set_hop(t2, 35, 37, [0, 0, 0])
my_model.set_hop(t3, 35, 38, [0, 0, 0])
my_model.set_hop(t2, 35, 39, [0, 0, 0])
my_model.set_hop(t3, 36, 1, [0, 0, 1])
my_model.set_hop(t3, 36, 3, [0, 0, 1])
my_model.set_hop(t3, 36, 31, [0, 0, 1])
my_model.set_hop(t3, 36, 33, [0, 0, 1])
my_model.set_hop(t2, 36, 34, [0, 0, 1])
my_model.set_hop(t1, 36, 35, [0, 0, 1])
my_model.set_hop(t3, 36, 37, [0, 0, 1])
my_model.set_hop(t2, 36, 38, [0, 0, 1])
my_model.set_hop(t1, 36, 39, [0, 0, 1])
my_model.set_hop(t1, 36, 37, [0, 0, 0])
my_model.set_hop(t2, 36, 38, [0, 0, 0])
my_model.set_hop(t3, 36, 39, [0, 0, 0])
my_model.set_hop(t1, 37, 38, [0, 0, 0])
my_model.set_hop(t2, 37, 39, [0, 0, 0])
my_model.set_hop(t1, 38, 39, [0, 0, 0])

my_model.display()

path = 'fullc'
nk = 800

(k_vec, k_dist, k_node) = my_model.k_path(path, nk)

print('-'*20)
print('starting calculation')
print('-'*20)
print('Calculating bands . . .')

evals = my_model.solve_all(k_vec)

plotly.tools.set_credentials_file(username='boiled-greggs', api_key='vLwVsrAXGmlgJzki3Q4a')
data = []
for band in evals:
    trace = go.Scatter(
        x = k_dist,
        y = band,
        mode='lines',
        showlegend = False
    )
    data.append(trace)

labels = [r'$\pi/a$',r'$\Gamma $',r'$\pi/a$']
vlines = list() 
for i, labels in enumerate(labels):
    vlines.append(
        go.Scatter(
            x = [k_node[i], k_node[i]],
            y = [-10, 10],
            mode = 'lines',
            line = go.Line(color='#111111', width=1),
            showlegend=False
        )
    )
bandxaxis = go.XAxis(
    title = 'k path',
    range=[0, k_dist[-1]],
    showgrid=True,
    showline=True,
    ticks="", 
    showticklabels=True,
    mirror=True,
    linewidth=2,
    ticktext=[r'$\pi/a$',r'$\Gamma $',r'$\pi/a$'],
    tickvals=k_node
    )
bandyaxis = go.YAxis(
    title="$E - E_f \quad (\\text{eV})$",
    range = [-7, 7],
    showgrid=True,
    showline=True,
    zeroline=True,
    mirror="ticks",
    ticks="inside",
    linewidth=2,
    tickwidth=2,
    zerolinewidth=2
    )
layout = go.Layout(
    title = '(10,0) tube bands',
    width = 500,
    height = 700,
    xaxis = bandxaxis,
    yaxis = bandyaxis
    )
fig = go.Figure(data=data + vlines, layout = layout) 
plt.iplot(fig, filename='plotly (10,0) tube bands')




"""
fig, ax = plt.subplots(figsize=(3,4))

# ax.set_xlim([0, k_node[-1]])
ax.set_xticks(k_node)
ax.set_xticklabels(labels)

for n in range(len(k_node)):
	ax.axvline(x=k_node[n], linewidth=0.3, color='black')

ax.set_title('10-0-tube band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(40):
	ax.plot(k_dist, evals[xband], linewidth=.5, color='black')

# fig.tight_layout()
fig.savefig('10-0-tube.eps')
"""
print('Done.')
