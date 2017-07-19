# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[16.9340000153, 0.0000000000, 0.0000000000], \
	[0.0000000000, 16.9340000153, 0.0000000000], \
	[0.0000000000, 0.0000000000, 2.4619998932]]

orb = [[0.901153803, 0.500000000, 0.499903798], \
	[0.892413855, 0.583443284, 0.499903798], \
	[0.881547928, 0.623954356, 0.000000000], \
	[0.847414672, 0.700606465, 0.000000000], \
	[0.824560761, 0.735802650, 0.499903798], \
	[0.768400431, 0.798163652, 0.499903798], \
	[0.735802650, 0.824560761, 0.000000000], \
	[0.663107157, 0.866489112, 0.000000000], \
	[0.623954356, 0.881547928, 0.499903798], \
	[0.541869283, 0.898968816, 0.499903798], \
	[0.500000000, 0.901153803, 0.000000000], \
	[0.416556716, 0.892413855, 0.000000000], \
	[0.376045644, 0.881547928, 0.499903798], \
	[0.299393564, 0.847414672, 0.499903798], \
	[0.264197379, 0.824560761, 0.000000000], \
	[0.201836377, 0.768400431, 0.000000000], \
	[0.175439239, 0.735802650, 0.499903798], \
	[0.133510888, 0.663107157, 0.499903798], \
	[0.118452102, 0.623954356, 0.000000000], \
	[0.101031184, 0.541869283, 0.000000000], \
	[0.098846167, 0.500000000, 0.499903798], \
	[0.107586175, 0.416556716, 0.499903798], \
	[0.118452102, 0.376045644, 0.000000000], \
	[0.152585328, 0.299393564, 0.000000000], \
	[0.175439239, 0.264197379, 0.499903798], \
	[0.231599569, 0.201836377, 0.499903798], \
	[0.264197379, 0.175439239, 0.000000000], \
	[0.336892843, 0.133510888, 0.000000000], \
	[0.376045644, 0.118452102, 0.499903798], \
	[0.458130717, 0.101031184, 0.499903798], \
	[0.500000000, 0.098846167, 0.000000000], \
	[0.583443284, 0.107586175, 0.000000000], \
	[0.623954356, 0.118452102, 0.499903798], \
	[0.700606465, 0.152585328, 0.499903798], \
	[0.735802650, 0.175439239, 0.000000000], \
	[0.798163652, 0.231599569, 0.000000000], \
	[0.824560761, 0.264197379, 0.499903798], \
	[0.866489112, 0.336892843, 0.499903798], \
	[0.881547928, 0.376045644, 0.000000000], \
	[0.898968816, 0.458130717, 0.000000000]]

my_model = tb_model(1, 3, lat, orb, per=[2])

# set model parameters
delta = 0.0
t1 = -1.8
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
my_model.set_hop(t2, 0, 0, [0, 0, 1])
my_model.set_hop(t3, 0, 1, [0, 0, 1])
my_model.set_hop(t2, 0, 2, [0, 0, 1])
my_model.set_hop(t3, 0, 3, [0, 0, 1])
my_model.set_hop(t3, 0, 37, [0, 0, 1])
my_model.set_hop(t2, 0, 38, [0, 0, 1])
my_model.set_hop(t1, 0, 39, [0, 0, 1])
my_model.set_hop(t1, 0, 1, [0, 0, 0])
my_model.set_hop(t2, 0, 2, [0, 0, 0])
my_model.set_hop(t3, 0, 3, [0, 0, 0])
my_model.set_hop(t3, 0, 37, [0, 0, 0])
my_model.set_hop(t2, 0, 38, [0, 0, 0])
my_model.set_hop(t1, 0, 39, [0, 0, 0])
my_model.set_hop(t3, 1, 0, [0, 0, 1])
my_model.set_hop(t2, 1, 1, [0, 0, 1])
my_model.set_hop(t1, 1, 2, [0, 0, 1])
my_model.set_hop(t2, 1, 3, [0, 0, 1])
my_model.set_hop(t3, 1, 4, [0, 0, 1])
my_model.set_hop(t3, 1, 38, [0, 0, 1])
my_model.set_hop(t2, 1, 39, [0, 0, 1])
my_model.set_hop(t1, 1, 2, [0, 0, 0])
my_model.set_hop(t2, 1, 3, [0, 0, 0])
my_model.set_hop(t3, 1, 4, [0, 0, 0])
my_model.set_hop(t3, 1, 38, [0, 0, 0])
my_model.set_hop(t2, 1, 39, [0, 0, 0])
my_model.set_hop(t1, 2, 3, [0, 0, 0])
my_model.set_hop(t2, 2, 4, [0, 0, 0])
my_model.set_hop(t3, 2, 5, [0, 0, 0])
my_model.set_hop(t3, 2, 39, [0, 0, 0])
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
my_model.set_hop(t3, 17, 14, [0, 0, 1])
my_model.set_hop(t2, 17, 15, [0, 0, 1])
my_model.set_hop(t3, 17, 16, [0, 0, 1])
my_model.set_hop(t2, 17, 17, [0, 0, 1])
my_model.set_hop(t1, 17, 18, [0, 0, 1])
my_model.set_hop(t2, 17, 19, [0, 0, 1])
my_model.set_hop(t3, 17, 20, [0, 0, 1])
my_model.set_hop(t1, 17, 18, [0, 0, 0])
my_model.set_hop(t2, 17, 19, [0, 0, 0])
my_model.set_hop(t3, 17, 20, [0, 0, 0])
my_model.set_hop(t1, 18, 19, [0, 0, 0])
my_model.set_hop(t2, 18, 20, [0, 0, 0])
my_model.set_hop(t3, 18, 21, [0, 0, 0])
my_model.set_hop(t1, 19, 20, [0, 0, 0])
my_model.set_hop(t2, 19, 21, [0, 0, 0])
my_model.set_hop(t3, 19, 22, [0, 0, 0])
my_model.set_hop(t3, 20, 17, [0, 0, 1])
my_model.set_hop(t2, 20, 18, [0, 0, 1])
my_model.set_hop(t1, 20, 19, [0, 0, 1])
my_model.set_hop(t2, 20, 20, [0, 0, 1])
my_model.set_hop(t3, 20, 21, [0, 0, 1])
my_model.set_hop(t2, 20, 22, [0, 0, 1])
my_model.set_hop(t3, 20, 23, [0, 0, 1])
my_model.set_hop(t1, 20, 21, [0, 0, 0])
my_model.set_hop(t2, 20, 22, [0, 0, 0])
my_model.set_hop(t3, 20, 23, [0, 0, 0])
my_model.set_hop(t3, 21, 18, [0, 0, 1])
my_model.set_hop(t2, 21, 19, [0, 0, 1])
my_model.set_hop(t3, 21, 20, [0, 0, 1])
my_model.set_hop(t2, 21, 21, [0, 0, 1])
my_model.set_hop(t1, 21, 22, [0, 0, 1])
my_model.set_hop(t2, 21, 23, [0, 0, 1])
my_model.set_hop(t3, 21, 24, [0, 0, 1])
my_model.set_hop(t1, 21, 22, [0, 0, 0])
my_model.set_hop(t2, 21, 23, [0, 0, 0])
my_model.set_hop(t3, 21, 24, [0, 0, 0])
my_model.set_hop(t1, 22, 23, [0, 0, 0])
my_model.set_hop(t2, 22, 24, [0, 0, 0])
my_model.set_hop(t3, 22, 25, [0, 0, 0])
my_model.set_hop(t1, 23, 24, [0, 0, 0])
my_model.set_hop(t2, 23, 25, [0, 0, 0])
my_model.set_hop(t3, 23, 26, [0, 0, 0])
my_model.set_hop(t3, 24, 21, [0, 0, 1])
my_model.set_hop(t2, 24, 22, [0, 0, 1])
my_model.set_hop(t1, 24, 23, [0, 0, 1])
my_model.set_hop(t2, 24, 24, [0, 0, 1])
my_model.set_hop(t3, 24, 25, [0, 0, 1])
my_model.set_hop(t2, 24, 26, [0, 0, 1])
my_model.set_hop(t3, 24, 27, [0, 0, 1])
my_model.set_hop(t1, 24, 25, [0, 0, 0])
my_model.set_hop(t2, 24, 26, [0, 0, 0])
my_model.set_hop(t3, 24, 27, [0, 0, 0])
my_model.set_hop(t3, 25, 22, [0, 0, 1])
my_model.set_hop(t2, 25, 23, [0, 0, 1])
my_model.set_hop(t3, 25, 24, [0, 0, 1])
my_model.set_hop(t2, 25, 25, [0, 0, 1])
my_model.set_hop(t1, 25, 26, [0, 0, 1])
my_model.set_hop(t2, 25, 27, [0, 0, 1])
my_model.set_hop(t3, 25, 28, [0, 0, 1])
my_model.set_hop(t1, 25, 26, [0, 0, 0])
my_model.set_hop(t2, 25, 27, [0, 0, 0])
my_model.set_hop(t3, 25, 28, [0, 0, 0])
my_model.set_hop(t1, 26, 27, [0, 0, 0])
my_model.set_hop(t2, 26, 28, [0, 0, 0])
my_model.set_hop(t3, 26, 29, [0, 0, 0])
my_model.set_hop(t1, 27, 28, [0, 0, 0])
my_model.set_hop(t2, 27, 29, [0, 0, 0])
my_model.set_hop(t3, 27, 30, [0, 0, 0])
my_model.set_hop(t3, 28, 25, [0, 0, 1])
my_model.set_hop(t2, 28, 26, [0, 0, 1])
my_model.set_hop(t1, 28, 27, [0, 0, 1])
my_model.set_hop(t2, 28, 28, [0, 0, 1])
my_model.set_hop(t3, 28, 29, [0, 0, 1])
my_model.set_hop(t2, 28, 30, [0, 0, 1])
my_model.set_hop(t3, 28, 31, [0, 0, 1])
my_model.set_hop(t1, 28, 29, [0, 0, 0])
my_model.set_hop(t2, 28, 30, [0, 0, 0])
my_model.set_hop(t3, 28, 31, [0, 0, 0])
my_model.set_hop(t3, 29, 26, [0, 0, 1])
my_model.set_hop(t2, 29, 27, [0, 0, 1])
my_model.set_hop(t3, 29, 28, [0, 0, 1])
my_model.set_hop(t2, 29, 29, [0, 0, 1])
my_model.set_hop(t1, 29, 30, [0, 0, 1])
my_model.set_hop(t2, 29, 31, [0, 0, 1])
my_model.set_hop(t3, 29, 32, [0, 0, 1])
my_model.set_hop(t1, 29, 30, [0, 0, 0])
my_model.set_hop(t2, 29, 31, [0, 0, 0])
my_model.set_hop(t3, 29, 32, [0, 0, 0])
my_model.set_hop(t1, 30, 31, [0, 0, 0])
my_model.set_hop(t2, 30, 32, [0, 0, 0])
my_model.set_hop(t3, 30, 33, [0, 0, 0])
my_model.set_hop(t1, 31, 32, [0, 0, 0])
my_model.set_hop(t2, 31, 33, [0, 0, 0])
my_model.set_hop(t3, 31, 34, [0, 0, 0])
my_model.set_hop(t3, 32, 29, [0, 0, 1])
my_model.set_hop(t2, 32, 30, [0, 0, 1])
my_model.set_hop(t1, 32, 31, [0, 0, 1])
my_model.set_hop(t2, 32, 32, [0, 0, 1])
my_model.set_hop(t3, 32, 33, [0, 0, 1])
my_model.set_hop(t2, 32, 34, [0, 0, 1])
my_model.set_hop(t3, 32, 35, [0, 0, 1])
my_model.set_hop(t1, 32, 33, [0, 0, 0])
my_model.set_hop(t2, 32, 34, [0, 0, 0])
my_model.set_hop(t3, 32, 35, [0, 0, 0])
my_model.set_hop(t3, 33, 30, [0, 0, 1])
my_model.set_hop(t2, 33, 31, [0, 0, 1])
my_model.set_hop(t3, 33, 32, [0, 0, 1])
my_model.set_hop(t2, 33, 33, [0, 0, 1])
my_model.set_hop(t1, 33, 34, [0, 0, 1])
my_model.set_hop(t2, 33, 35, [0, 0, 1])
my_model.set_hop(t3, 33, 36, [0, 0, 1])
my_model.set_hop(t1, 33, 34, [0, 0, 0])
my_model.set_hop(t2, 33, 35, [0, 0, 0])
my_model.set_hop(t3, 33, 36, [0, 0, 0])
my_model.set_hop(t1, 34, 35, [0, 0, 0])
my_model.set_hop(t2, 34, 36, [0, 0, 0])
my_model.set_hop(t3, 34, 37, [0, 0, 0])
my_model.set_hop(t1, 35, 36, [0, 0, 0])
my_model.set_hop(t2, 35, 37, [0, 0, 0])
my_model.set_hop(t3, 35, 38, [0, 0, 0])
my_model.set_hop(t3, 36, 33, [0, 0, 1])
my_model.set_hop(t2, 36, 34, [0, 0, 1])
my_model.set_hop(t1, 36, 35, [0, 0, 1])
my_model.set_hop(t2, 36, 36, [0, 0, 1])
my_model.set_hop(t3, 36, 37, [0, 0, 1])
my_model.set_hop(t2, 36, 38, [0, 0, 1])
my_model.set_hop(t3, 36, 39, [0, 0, 1])
my_model.set_hop(t1, 36, 37, [0, 0, 0])
my_model.set_hop(t2, 36, 38, [0, 0, 0])
my_model.set_hop(t3, 36, 39, [0, 0, 0])
my_model.set_hop(t3, 37, 0, [0, 0, 1])
my_model.set_hop(t3, 37, 34, [0, 0, 1])
my_model.set_hop(t2, 37, 35, [0, 0, 1])
my_model.set_hop(t3, 37, 36, [0, 0, 1])
my_model.set_hop(t2, 37, 37, [0, 0, 1])
my_model.set_hop(t1, 37, 38, [0, 0, 1])
my_model.set_hop(t2, 37, 39, [0, 0, 1])
my_model.set_hop(t1, 37, 38, [0, 0, 0])
my_model.set_hop(t2, 37, 39, [0, 0, 0])
my_model.set_hop(t1, 38, 39, [0, 0, 0])

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

ax.set_title('10-10 band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(40):
	ax.plot(k_dist, evals[xband], linewidth=.5, color='black')

# fig.tight_layout()
fig.savefig('10-10.eps')

print('Done.')