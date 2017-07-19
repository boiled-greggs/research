# Greg Stewart 2017

from __future__ import print_function
from pythtb import *
import numpy as np
import pylab as plt

lat = [[13.0763549805, 0.0000000000, 0.0000000000], \
	[-4.3587852709, 12.3285056084, 0.0000000000], \
	[-4.3587852709, -6.1642533927, 10.6767987078]]

orb = [[0.329797983, 0.558250964, 0.679593027], \
	[0.320436746, 0.441821247, 0.670235932], \
	[0.373703063, 0.596954465, 0.596916139], \
	[0.626296759, 0.223213106, 0.223251432], \
	[0.000038147, 0.776786685, 0.403083593], \
	[0.000038147, 0.596954465, 0.223251432], \
	[0.403045267, 0.403083593, 0.626296759], \
	[0.403045267, 0.776748359, 0.999961615], \
	[0.223213106, 0.223251432, 0.626296759], \
	[0.223213106, 0.596916139, 0.999961615], \
	[0.403083593, 0.626296759, 0.403045267], \
	[0.403083593, 0.000038147, 0.776786685], \
	[0.776748359, 0.999961615, 0.403045267], \
	[0.776748359, 0.373703063, 0.776786685], \
	[0.223251432, 0.626296759, 0.223213106], \
	[0.223251432, 0.000038147, 0.596954465], \
	[0.596916139, 0.999961615, 0.223213106], \
	[0.596916139, 0.373703063, 0.596954465], \
	[0.776786685, 0.403083593, 0.000038147], \
	[0.776786685, 0.776748359, 0.373703063], \
	[0.596954465, 0.223251432, 0.000038147], \
	[0.596954465, 0.596916139, 0.373703063], \
	[0.999961615, 0.403045267, 0.776748359], \
	[0.999961615, 0.223213106, 0.596916139], \
	[0.373703063, 0.776786685, 0.776748359], \
	[0.626296759, 0.403045267, 0.403083593], \
	[0.329763860, 0.679563046, 0.558178484], \
	[0.670235932, 0.228414744, 0.349799216], \
	[0.121384487, 0.771585107, 0.441821247], \
	[0.121384487, 0.679563046, 0.349799216], \
	[0.320436746, 0.650200546, 0.878615320], \
	[0.228414744, 0.349799216, 0.670235932], \
	[0.228414744, 0.558178484, 0.878615320], \
	[0.441821247, 0.670235932, 0.320436746], \
	[0.441821247, 0.121384487, 0.771585107], \
	[0.650200546, 0.878615320, 0.320436746], \
	[0.650200546, 0.329763860, 0.771585107], \
	[0.349799216, 0.670235932, 0.228414744], \
	[0.349799216, 0.121384487, 0.679563046], \
	[0.558178484, 0.878615320, 0.228414744], \
	[0.558178484, 0.329763860, 0.679563046], \
	[0.771585107, 0.441821247, 0.121384487], \
	[0.771585107, 0.650200546, 0.329763860], \
	[0.679563046, 0.349799216, 0.121384487], \
	[0.679563046, 0.558178484, 0.329763860], \
	[0.878615320, 0.320436746, 0.650200546], \
	[0.878615320, 0.228414744, 0.558178484], \
	[0.329763860, 0.771585107, 0.650200546], \
	[0.670235932, 0.320436746, 0.441821247], \
	[0.670201838, 0.349795043, 0.228453010], \
	[0.878657758, 0.650204718, 0.320406765], \
	[0.878657758, 0.558250964, 0.228453010], \
	[0.441748828, 0.320406765, 0.670201838], \
	[0.441748828, 0.771546841, 0.121342048], \
	[0.349795043, 0.228453010, 0.670201838], \
	[0.349795043, 0.679593027, 0.121342048], \
	[0.320406765, 0.670201838, 0.441748828], \
	[0.320406765, 0.878657758, 0.650204718], \
	[0.771546841, 0.121342048, 0.441748828], \
	[0.771546841, 0.329797983, 0.650204718], \
	[0.228453010, 0.670201838, 0.349795043], \
	[0.228453010, 0.878657758, 0.558250964], \
	[0.679593027, 0.121342048, 0.349795043], \
	[0.679593027, 0.329797983, 0.558250964], \
	[0.650204718, 0.320406765, 0.878657758], \
	[0.650204718, 0.771546841, 0.329797983], \
	[0.558250964, 0.228453010, 0.878657758], \
	[0.558250964, 0.679593027, 0.329797983], \
	[0.121342048, 0.441748828, 0.771546841], \
	[0.121342048, 0.349795043, 0.679593027], \
	[0.329797983, 0.650204718, 0.771546841], \
	[0.670201838, 0.441748828, 0.320406765], \
	[0.422558755, 0.538737476, 0.538703382], \
	[0.577441037, 0.116144642, 0.116178736], \
	[0.000033975, 0.883855164, 0.461296380], \
	[0.000033975, 0.538737476, 0.116178736], \
	[0.461262316, 0.461296380, 0.577441037], \
	[0.461262316, 0.883821070, 0.999965847], \
	[0.116144642, 0.116178736, 0.577441037], \
	[0.116144642, 0.538703382, 0.999965847], \
	[0.461296380, 0.577441037, 0.461262316], \
	[0.461296380, 0.000033975, 0.883855164], \
	[0.883821070, 0.999965847, 0.461262316], \
	[0.883821070, 0.422558755, 0.883855164], \
	[0.116178736, 0.577441037, 0.116144642], \
	[0.116178736, 0.000033975, 0.538737476], \
	[0.538703382, 0.999965847, 0.116144642], \
	[0.538703382, 0.422558755, 0.538737476], \
	[0.883855164, 0.461296380, 0.000033975], \
	[0.883855164, 0.883821070, 0.422558755], \
	[0.538737476, 0.116178736, 0.000033975], \
	[0.538737476, 0.538703382, 0.422558755], \
	[0.999965847, 0.461262316, 0.883821070], \
	[0.999965847, 0.116144642, 0.538703382], \
	[0.422558755, 0.883855164, 0.883821070], \
	[0.577441037, 0.461262316, 0.461296380]]

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
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta, delta, delta, delta, delta, \
	delta])

# set hopping parameters for connected orbitals
# (amplitude, i, j, [lattice vector to cell containing j])
my_model.set_hop(2.7222, 0, 1, [0, 0, 0])
my_model.set_hop(2.2543, 0, 2, [0, 0, 0])
my_model.set_hop(0.3503, 0, 6, [0, 0, 0])
my_model.set_hop(0.3714, 0, 24, [0, 0, 0])
my_model.set_hop(0.3162, 0, 26, [0, 0, 0])
my_model.set_hop(0.3375, 0, 30, [0, 0, 0])
my_model.set_hop(0.3370, 0, 31, [0, 0, 0])
my_model.set_hop(0.2222, 0, 47, [0, 0, 0])
my_model.set_hop(3.0981, 0, 70, [0, 0, 0])
my_model.set_hop(0.3507, 0, 72, [0, 0, 0])
my_model.set_hop(0.2367, 0, 76, [0, 0, 0])
my_model.set_hop(0.3502, 1, 2, [0, 0, 0])
my_model.set_hop(2.2564, 1, 6, [0, 0, 0])
my_model.set_hop(0.3712, 1, 8, [0, 0, 0])
my_model.set_hop(3.0908, 1, 31, [0, 0, 0])
my_model.set_hop(0.3162, 1, 52, [0, 0, 0])
my_model.set_hop(0.2222, 1, 54, [0, 0, 0])
my_model.set_hop(0.3370, 1, 69, [0, 0, 0])
my_model.set_hop(0.3375, 1, 70, [0, 0, 0])
my_model.set_hop(0.2366, 1, 72, [0, 0, 0])
my_model.set_hop(0.3509, 1, 76, [0, 0, 0])
my_model.set_hop(0.2476, 2, 6, [0, 0, 0])
my_model.set_hop(0.2476, 2, 10, [0, 0, 0])
my_model.set_hop(0.2737, 2, 24, [0, 0, 0])
my_model.set_hop(2.2564, 2, 26, [0, 0, 0])
my_model.set_hop(0.3712, 2, 47, [0, 0, 0])
my_model.set_hop(0.3503, 2, 56, [0, 0, 0])
my_model.set_hop(0.3714, 2, 70, [0, 0, 0])
my_model.set_hop(3.3626, 2, 72, [0, 0, 0])
my_model.set_hop(0.4115, 2, 76, [0, 0, 0])
my_model.set_hop(0.4113, 2, 80, [0, 0, 0])
my_model.set_hop(0.2476, 3, 20, [0, 0, 0])
my_model.set_hop(0.2737, 3, 25, [0, 0, 0])
my_model.set_hop(2.2564, 3, 27, [0, 0, 0])
my_model.set_hop(0.3502, 3, 43, [0, 0, 0])
my_model.set_hop(0.3712, 3, 48, [0, 0, 0])
my_model.set_hop(2.2543, 3, 49, [0, 0, 0])
my_model.set_hop(0.3503, 3, 62, [0, 0, 0])
my_model.set_hop(0.3714, 3, 71, [0, 0, 0])
my_model.set_hop(3.3626, 3, 73, [0, 0, 0])
my_model.set_hop(0.4115, 3, 90, [0, 0, 0])
my_model.set_hop(0.2737, 4, 5, [0, 0, 0])
my_model.set_hop(2.2564, 4, 28, [0, 0, 0])
my_model.set_hop(0.3712, 4, 29, [0, 0, 0])
my_model.set_hop(0.3503, 4, 61, [0, 0, 0])
my_model.set_hop(3.3626, 4, 74, [0, 0, 0])
my_model.set_hop(0.2476, 5, 14, [0, 0, 0])
my_model.set_hop(0.3712, 5, 28, [0, 0, 0])
my_model.set_hop(2.2564, 5, 29, [0, 0, 0])
my_model.set_hop(0.3503, 5, 60, [0, 0, 0])
my_model.set_hop(3.3626, 5, 75, [0, 0, 0])
my_model.set_hop(0.4113, 5, 84, [0, 0, 0])
my_model.set_hop(0.2737, 6, 8, [0, 0, 0])
my_model.set_hop(0.2476, 6, 17, [0, 0, 0])
my_model.set_hop(0.3712, 6, 31, [0, 0, 0])
my_model.set_hop(0.3502, 6, 40, [0, 0, 0])
my_model.set_hop(2.2543, 6, 52, [0, 0, 0])
my_model.set_hop(0.3714, 6, 54, [0, 0, 0])
my_model.set_hop(0.4113, 6, 72, [0, 0, 0])
my_model.set_hop(3.3626, 6, 76, [0, 0, 0])
my_model.set_hop(0.4115, 6, 87, [0, 0, 0])
my_model.set_hop(0.2476, 7, 16, [0, 0, 1])
my_model.set_hop(0.3502, 7, 39, [0, 0, 1])
my_model.set_hop(2.2543, 7, 53, [0, 0, 1])
my_model.set_hop(0.3714, 7, 55, [0, 0, 1])
my_model.set_hop(0.4115, 7, 86, [0, 0, 1])
my_model.set_hop(0.2737, 7, 9, [0, 0, 0])
my_model.set_hop(0.2476, 7, 24, [0, 0, 0])
my_model.set_hop(2.2564, 7, 30, [0, 0, 0])
my_model.set_hop(0.3712, 7, 32, [0, 0, 0])
my_model.set_hop(0.3503, 7, 70, [0, 0, 0])
my_model.set_hop(3.3626, 7, 77, [0, 0, 0])
my_model.set_hop(0.4113, 7, 94, [0, 0, 0])
my_model.set_hop(0.2476, 8, 15, [0, 0, 0])
my_model.set_hop(2.2564, 8, 31, [0, 0, 0])
my_model.set_hop(0.3502, 8, 38, [0, 0, 0])
my_model.set_hop(0.3714, 8, 52, [0, 0, 0])
my_model.set_hop(2.2543, 8, 54, [0, 0, 0])
my_model.set_hop(0.3503, 8, 69, [0, 0, 0])
my_model.set_hop(3.3626, 8, 78, [0, 0, 0])
my_model.set_hop(0.4115, 8, 85, [0, 0, 0])
my_model.set_hop(0.2476, 9, 14, [0, 0, 1])
my_model.set_hop(0.3502, 9, 37, [0, 0, 1])
my_model.set_hop(0.3714, 9, 53, [0, 0, 1])
my_model.set_hop(2.2543, 9, 55, [0, 0, 1])
my_model.set_hop(0.4115, 9, 84, [0, 0, 1])
my_model.set_hop(0.3712, 9, 30, [0, 0, 0])
my_model.set_hop(2.2564, 9, 32, [0, 0, 0])
my_model.set_hop(0.3503, 9, 68, [0, 0, 0])
my_model.set_hop(3.3626, 9, 79, [0, 0, 0])
my_model.set_hop(0.2737, 10, 14, [0, 0, 0])
my_model.set_hop(0.2476, 10, 21, [0, 0, 0])
my_model.set_hop(0.3502, 10, 26, [0, 0, 0])
my_model.set_hop(2.2564, 10, 33, [0, 0, 0])
my_model.set_hop(0.3712, 10, 37, [0, 0, 0])
my_model.set_hop(2.2543, 10, 56, [0, 0, 0])
my_model.set_hop(0.3714, 10, 60, [0, 0, 0])
my_model.set_hop(0.3503, 10, 67, [0, 0, 0])
my_model.set_hop(0.4115, 10, 72, [0, 0, 0])
my_model.set_hop(3.3626, 10, 80, [0, 0, 0])
my_model.set_hop(0.4113, 10, 91, [0, 0, 0])
my_model.set_hop(0.2737, 11, 15, [0, 0, 0])
my_model.set_hop(2.2564, 11, 34, [0, 0, 0])
my_model.set_hop(0.3712, 11, 38, [0, 0, 0])
my_model.set_hop(0.3503, 11, 66, [0, 0, 0])
my_model.set_hop(3.3626, 11, 81, [0, 0, 0])
my_model.set_hop(0.2476, 12, 23, [0, 1, 0])
my_model.set_hop(0.3502, 12, 46, [0, 1, 0])
my_model.set_hop(2.2543, 12, 58, [0, 1, 0])
my_model.set_hop(0.3714, 12, 62, [0, 1, 0])
my_model.set_hop(0.4115, 12, 93, [0, 1, 0])
my_model.set_hop(0.2737, 12, 16, [0, 0, 0])
my_model.set_hop(0.2476, 12, 19, [0, 0, 0])
my_model.set_hop(2.2564, 12, 35, [0, 0, 0])
my_model.set_hop(0.3712, 12, 39, [0, 0, 0])
my_model.set_hop(0.3503, 12, 65, [0, 0, 0])
my_model.set_hop(3.3626, 12, 82, [0, 0, 0])
my_model.set_hop(0.4113, 12, 89, [0, 0, 0])
my_model.set_hop(0.2737, 13, 17, [0, 0, 0])
my_model.set_hop(0.2476, 13, 22, [0, 0, 0])
my_model.set_hop(2.2564, 13, 36, [0, 0, 0])
my_model.set_hop(0.3712, 13, 40, [0, 0, 0])
my_model.set_hop(0.3502, 13, 45, [0, 0, 0])
my_model.set_hop(2.2543, 13, 59, [0, 0, 0])
my_model.set_hop(0.3714, 13, 63, [0, 0, 0])
my_model.set_hop(0.3503, 13, 64, [0, 0, 0])
my_model.set_hop(3.3626, 13, 83, [0, 0, 0])
my_model.set_hop(0.4115, 13, 92, [0, 0, 0])
my_model.set_hop(0.3502, 14, 29, [0, 0, 0])
my_model.set_hop(0.3712, 14, 33, [0, 0, 0])
my_model.set_hop(2.2564, 14, 37, [0, 0, 0])
my_model.set_hop(0.3503, 14, 55, [0, 0, 0])
my_model.set_hop(0.3714, 14, 56, [0, 0, 0])
my_model.set_hop(2.2543, 14, 60, [0, 0, 0])
my_model.set_hop(0.4115, 14, 75, [0, 0, 0])
my_model.set_hop(3.3626, 14, 84, [0, 0, 0])
my_model.set_hop(0.3712, 15, 34, [0, 0, 0])
my_model.set_hop(2.2564, 15, 38, [0, 0, 0])
my_model.set_hop(0.3503, 15, 54, [0, 0, 0])
my_model.set_hop(0.4113, 15, 78, [0, 0, 0])
my_model.set_hop(3.3626, 15, 85, [0, 0, 0])
my_model.set_hop(0.2476, 16, 3, [0, 1, 0])
my_model.set_hop(0.3502, 16, 27, [0, 1, 0])
my_model.set_hop(0.3714, 16, 58, [0, 1, 0])
my_model.set_hop(2.2543, 16, 62, [0, 1, 0])
my_model.set_hop(0.4115, 16, 73, [0, 1, 0])
my_model.set_hop(0.3712, 16, 35, [0, 0, 0])
my_model.set_hop(2.2564, 16, 39, [0, 0, 0])
my_model.set_hop(0.3503, 16, 53, [0, 0, 0])
my_model.set_hop(3.3626, 16, 86, [0, 0, 0])
my_model.set_hop(0.2476, 17, 25, [0, 0, 0])
my_model.set_hop(0.3712, 17, 36, [0, 0, 0])
my_model.set_hop(2.2564, 17, 40, [0, 0, 0])
my_model.set_hop(0.3502, 17, 48, [0, 0, 0])
my_model.set_hop(0.3503, 17, 52, [0, 0, 0])
my_model.set_hop(0.3714, 17, 59, [0, 0, 0])
my_model.set_hop(2.2543, 17, 63, [0, 0, 0])
my_model.set_hop(0.4113, 17, 76, [0, 0, 0])
my_model.set_hop(3.3626, 17, 87, [0, 0, 0])
my_model.set_hop(0.4115, 17, 95, [0, 0, 0])
my_model.set_hop(0.2737, 18, 20, [0, 0, 0])
my_model.set_hop(2.2564, 18, 41, [0, 0, 0])
my_model.set_hop(0.3712, 18, 43, [0, 0, 0])
my_model.set_hop(0.3503, 18, 51, [0, 0, 0])
my_model.set_hop(3.3626, 18, 88, [0, 0, 0])
my_model.set_hop(0.2737, 19, 21, [0, 0, 0])
my_model.set_hop(0.3502, 19, 35, [0, 0, 0])
my_model.set_hop(2.2564, 19, 42, [0, 0, 0])
my_model.set_hop(0.3712, 19, 44, [0, 0, 0])
my_model.set_hop(0.3503, 19, 50, [0, 0, 0])
my_model.set_hop(2.2543, 19, 65, [0, 0, 0])
my_model.set_hop(0.3714, 19, 67, [0, 0, 0])
my_model.set_hop(0.4115, 19, 82, [0, 0, 0])
my_model.set_hop(3.3626, 19, 89, [0, 0, 0])
my_model.set_hop(0.3712, 20, 41, [0, 0, 0])
my_model.set_hop(2.2564, 20, 43, [0, 0, 0])
my_model.set_hop(0.3503, 20, 49, [0, 0, 0])
my_model.set_hop(0.4113, 20, 73, [0, 0, 0])
my_model.set_hop(3.3626, 20, 90, [0, 0, 0])
my_model.set_hop(0.2476, 21, 25, [0, 0, 0])
my_model.set_hop(0.3502, 21, 33, [0, 0, 0])
my_model.set_hop(0.3712, 21, 42, [0, 0, 0])
my_model.set_hop(2.2564, 21, 44, [0, 0, 0])
my_model.set_hop(0.3714, 21, 65, [0, 0, 0])
my_model.set_hop(2.2543, 21, 67, [0, 0, 0])
my_model.set_hop(0.3503, 21, 71, [0, 0, 0])
my_model.set_hop(0.4115, 21, 80, [0, 0, 0])
my_model.set_hop(3.3626, 21, 91, [0, 0, 0])
my_model.set_hop(0.4113, 21, 95, [0, 0, 0])
my_model.set_hop(0.2476, 22, 9, [1, 0, 0])
my_model.set_hop(0.3502, 22, 32, [1, 0, 0])
my_model.set_hop(2.2543, 22, 68, [1, 0, 0])
my_model.set_hop(0.3714, 22, 69, [1, 0, 0])
my_model.set_hop(0.4115, 22, 79, [1, 0, 0])
my_model.set_hop(0.2737, 22, 23, [0, 0, 0])
my_model.set_hop(2.2564, 22, 45, [0, 0, 0])
my_model.set_hop(0.3712, 22, 46, [0, 0, 0])
my_model.set_hop(0.3503, 22, 59, [0, 0, 0])
my_model.set_hop(0.4113, 22, 83, [0, 0, 0])
my_model.set_hop(3.3626, 22, 92, [0, 0, 0])
my_model.set_hop(0.2476, 23, 8, [1, 0, 0])
my_model.set_hop(0.3502, 23, 31, [1, 0, 0])
my_model.set_hop(0.3714, 23, 68, [1, 0, 0])
my_model.set_hop(2.2543, 23, 69, [1, 0, 0])
my_model.set_hop(0.4115, 23, 78, [1, 0, 0])
my_model.set_hop(0.3712, 23, 45, [0, 0, 0])
my_model.set_hop(2.2564, 23, 46, [0, 0, 0])
my_model.set_hop(0.3503, 23, 58, [0, 0, 0])
my_model.set_hop(3.3626, 23, 93, [0, 0, 0])
my_model.set_hop(0.3712, 24, 26, [0, 0, 0])
my_model.set_hop(0.3502, 24, 30, [0, 0, 0])
my_model.set_hop(2.2564, 24, 47, [0, 0, 0])
my_model.set_hop(0.3503, 24, 57, [0, 0, 0])
my_model.set_hop(2.2543, 24, 70, [0, 0, 0])
my_model.set_hop(0.4115, 24, 77, [0, 0, 0])
my_model.set_hop(3.3626, 24, 94, [0, 0, 0])
my_model.set_hop(0.3712, 25, 27, [0, 0, 0])
my_model.set_hop(0.3502, 25, 44, [0, 0, 0])
my_model.set_hop(2.2564, 25, 48, [0, 0, 0])
my_model.set_hop(0.3714, 25, 49, [0, 0, 0])
my_model.set_hop(0.3503, 25, 63, [0, 0, 0])
my_model.set_hop(2.2543, 25, 71, [0, 0, 0])
my_model.set_hop(0.4113, 25, 87, [0, 0, 0])
my_model.set_hop(0.4115, 25, 91, [0, 0, 0])
my_model.set_hop(3.3626, 25, 95, [0, 0, 0])
my_model.set_hop(3.0908, 26, 47, [0, 0, 0])
my_model.set_hop(2.7222, 26, 56, [0, 0, 0])
my_model.set_hop(0.3370, 26, 57, [0, 0, 0])
my_model.set_hop(0.3375, 26, 60, [0, 0, 0])
my_model.set_hop(0.2222, 26, 70, [0, 0, 0])
my_model.set_hop(0.3509, 26, 72, [0, 0, 0])
my_model.set_hop(0.2366, 26, 80, [0, 0, 0])
my_model.set_hop(3.0908, 27, 48, [0, 0, 0])
my_model.set_hop(0.3162, 27, 49, [0, 0, 0])
my_model.set_hop(0.3375, 27, 58, [0, 0, 0])
my_model.set_hop(2.7222, 27, 62, [0, 0, 0])
my_model.set_hop(0.3370, 27, 63, [0, 0, 0])
my_model.set_hop(0.2222, 27, 71, [0, 0, 0])
my_model.set_hop(0.3509, 27, 73, [0, 0, 0])
my_model.set_hop(3.0908, 28, 29, [0, 0, 0])
my_model.set_hop(0.3375, 28, 57, [0, 0, 0])
my_model.set_hop(0.3370, 28, 60, [0, 0, 0])
my_model.set_hop(2.7222, 28, 61, [0, 0, 0])
my_model.set_hop(0.3509, 28, 74, [0, 0, 0])
my_model.set_hop(0.3375, 29, 56, [0, 0, 0])
my_model.set_hop(2.7222, 29, 60, [0, 0, 0])
my_model.set_hop(0.3370, 29, 61, [0, 0, 0])
my_model.set_hop(0.3509, 29, 75, [0, 0, 0])
my_model.set_hop(0.2366, 29, 84, [0, 0, 0])
my_model.set_hop(3.0908, 30, 32, [0, 0, 0])
my_model.set_hop(0.3370, 30, 68, [0, 0, 0])
my_model.set_hop(2.7222, 30, 70, [0, 0, 0])
my_model.set_hop(0.3509, 30, 77, [0, 0, 0])
my_model.set_hop(0.2366, 30, 94, [0, 0, 0])
my_model.set_hop(0.2222, 31, 52, [0, 0, 0])
my_model.set_hop(0.3162, 31, 54, [0, 0, 0])
my_model.set_hop(0.3375, 31, 68, [0, 0, 0])
my_model.set_hop(2.7222, 31, 69, [0, 0, 0])
my_model.set_hop(0.3509, 31, 78, [0, 0, 0])
my_model.set_hop(2.7222, 32, 68, [0, 0, 0])
my_model.set_hop(0.3375, 32, 69, [0, 0, 0])
my_model.set_hop(0.3370, 32, 70, [0, 0, 0])
my_model.set_hop(0.3509, 32, 79, [0, 0, 0])
my_model.set_hop(3.0908, 33, 37, [0, 0, 0])
my_model.set_hop(0.3370, 33, 55, [0, 0, 0])
my_model.set_hop(0.3162, 33, 56, [0, 0, 0])
my_model.set_hop(0.2222, 33, 60, [0, 0, 0])
my_model.set_hop(0.3375, 33, 65, [0, 0, 0])
my_model.set_hop(2.7222, 33, 67, [0, 0, 0])
my_model.set_hop(0.3509, 33, 80, [0, 0, 0])
my_model.set_hop(0.2366, 33, 91, [0, 0, 0])
my_model.set_hop(3.0908, 34, 38, [0, 0, 0])
my_model.set_hop(0.3370, 34, 54, [0, 0, 0])
my_model.set_hop(0.3375, 34, 64, [0, 0, 0])
my_model.set_hop(2.7222, 34, 66, [0, 0, 0])
my_model.set_hop(0.3509, 34, 81, [0, 0, 0])
my_model.set_hop(3.0908, 35, 39, [0, 0, 0])
my_model.set_hop(0.3370, 35, 53, [0, 0, 0])
my_model.set_hop(2.7222, 35, 65, [0, 0, 0])
my_model.set_hop(0.3375, 35, 67, [0, 0, 0])
my_model.set_hop(0.3509, 35, 82, [0, 0, 0])
my_model.set_hop(0.2366, 35, 89, [0, 0, 0])
my_model.set_hop(3.0908, 36, 40, [0, 0, 0])
my_model.set_hop(0.3370, 36, 52, [0, 0, 0])
my_model.set_hop(0.3162, 36, 59, [0, 0, 0])
my_model.set_hop(0.2222, 36, 63, [0, 0, 0])
my_model.set_hop(2.7222, 36, 64, [0, 0, 0])
my_model.set_hop(0.3375, 36, 66, [0, 0, 0])
my_model.set_hop(0.3509, 36, 83, [0, 0, 0])
my_model.set_hop(0.3375, 37, 53, [0, 0, 0])
my_model.set_hop(2.7222, 37, 55, [0, 0, 0])
my_model.set_hop(0.2222, 37, 56, [0, 0, 0])
my_model.set_hop(0.3162, 37, 60, [0, 0, 0])
my_model.set_hop(0.3370, 37, 67, [0, 0, 0])
my_model.set_hop(0.3509, 37, 84, [0, 0, 0])
my_model.set_hop(0.3375, 38, 52, [0, 0, 0])
my_model.set_hop(2.7222, 38, 54, [0, 0, 0])
my_model.set_hop(0.3370, 38, 66, [0, 0, 0])
my_model.set_hop(0.2366, 38, 78, [0, 0, 0])
my_model.set_hop(0.3509, 38, 85, [0, 0, 0])
my_model.set_hop(2.7222, 39, 53, [0, 0, 0])
my_model.set_hop(0.3375, 39, 55, [0, 0, 0])
my_model.set_hop(0.3370, 39, 65, [0, 0, 0])
my_model.set_hop(0.3509, 39, 86, [0, 0, 0])
my_model.set_hop(2.7222, 40, 52, [0, 0, 0])
my_model.set_hop(0.3375, 40, 54, [0, 0, 0])
my_model.set_hop(0.2222, 40, 59, [0, 0, 0])
my_model.set_hop(0.3162, 40, 63, [0, 0, 0])
my_model.set_hop(0.3370, 40, 64, [0, 0, 0])
my_model.set_hop(0.2366, 40, 76, [0, 0, 0])
my_model.set_hop(0.3509, 40, 87, [0, 0, 0])
my_model.set_hop(3.0908, 41, 43, [0, 0, 0])
my_model.set_hop(0.3370, 41, 49, [0, 0, 0])
my_model.set_hop(0.3375, 41, 50, [0, 0, 0])
my_model.set_hop(2.7222, 41, 51, [0, 0, 0])
my_model.set_hop(0.3509, 41, 88, [0, 0, 0])
my_model.set_hop(3.0908, 42, 44, [0, 0, 0])
my_model.set_hop(2.7222, 42, 50, [0, 0, 0])
my_model.set_hop(0.3375, 42, 51, [0, 0, 0])
my_model.set_hop(0.3162, 42, 65, [0, 0, 0])
my_model.set_hop(0.2222, 42, 67, [0, 0, 0])
my_model.set_hop(0.3370, 42, 71, [0, 0, 0])
my_model.set_hop(0.3509, 42, 89, [0, 0, 0])
my_model.set_hop(2.7222, 43, 49, [0, 0, 0])
my_model.set_hop(0.3370, 43, 51, [0, 0, 0])
my_model.set_hop(0.3375, 43, 71, [0, 0, 0])
my_model.set_hop(0.2366, 43, 73, [0, 0, 0])
my_model.set_hop(0.3509, 43, 90, [0, 0, 0])
my_model.set_hop(0.3375, 44, 49, [0, 0, 0])
my_model.set_hop(0.3370, 44, 50, [0, 0, 0])
my_model.set_hop(0.2222, 44, 65, [0, 0, 0])
my_model.set_hop(0.3162, 44, 67, [0, 0, 0])
my_model.set_hop(2.7222, 44, 71, [0, 0, 0])
my_model.set_hop(0.3509, 44, 91, [0, 0, 0])
my_model.set_hop(0.2366, 44, 95, [0, 0, 0])
my_model.set_hop(3.0908, 45, 46, [0, 0, 0])
my_model.set_hop(0.3370, 45, 58, [0, 0, 0])
my_model.set_hop(2.7222, 45, 59, [0, 0, 0])
my_model.set_hop(0.3375, 45, 63, [0, 0, 0])
my_model.set_hop(0.2366, 45, 83, [0, 0, 0])
my_model.set_hop(0.3509, 45, 92, [0, 0, 0])
my_model.set_hop(2.7222, 46, 58, [0, 0, 0])
my_model.set_hop(0.3370, 46, 59, [0, 0, 0])
my_model.set_hop(0.3375, 46, 62, [0, 0, 0])
my_model.set_hop(0.3509, 46, 93, [0, 0, 0])
my_model.set_hop(0.3370, 47, 56, [0, 0, 0])
my_model.set_hop(2.7222, 47, 57, [0, 0, 0])
my_model.set_hop(0.3375, 47, 61, [0, 0, 0])
my_model.set_hop(0.3162, 47, 70, [0, 0, 0])
my_model.set_hop(0.3509, 47, 94, [0, 0, 0])
my_model.set_hop(0.2222, 48, 49, [0, 0, 0])
my_model.set_hop(0.3375, 48, 59, [0, 0, 0])
my_model.set_hop(0.3370, 48, 62, [0, 0, 0])
my_model.set_hop(2.7222, 48, 63, [0, 0, 0])
my_model.set_hop(0.3162, 48, 71, [0, 0, 0])
my_model.set_hop(0.2366, 48, 87, [0, 0, 0])
my_model.set_hop(0.3509, 48, 95, [0, 0, 0])
my_model.set_hop(3.0981, 49, 71, [0, 0, 0])
my_model.set_hop(0.3507, 49, 73, [0, 0, 0])
my_model.set_hop(0.2367, 49, 90, [0, 0, 0])
my_model.set_hop(3.0981, 50, 51, [0, 0, 0])
my_model.set_hop(0.2367, 50, 89, [0, 0, 0])
my_model.set_hop(0.2367, 51, 88, [0, 0, 0])
my_model.set_hop(3.0981, 52, 54, [0, 0, 0])
my_model.set_hop(0.3507, 52, 76, [0, 0, 0])
my_model.set_hop(0.2367, 52, 87, [0, 0, 0])
my_model.set_hop(3.0981, 53, 55, [0, 0, 0])
my_model.set_hop(0.2367, 53, 86, [0, 0, 0])
my_model.set_hop(0.3507, 54, 78, [0, 0, 0])
my_model.set_hop(0.2367, 54, 85, [0, 0, 0])
my_model.set_hop(0.2367, 55, 84, [0, 0, 0])
my_model.set_hop(3.0981, 56, 60, [0, 0, 0])
my_model.set_hop(0.2367, 56, 72, [0, 0, 0])
my_model.set_hop(0.3507, 56, 80, [0, 0, 0])
my_model.set_hop(3.0981, 57, 61, [0, 0, 0])
my_model.set_hop(0.2367, 57, 94, [0, 0, 0])
my_model.set_hop(3.0981, 58, 62, [0, 0, 0])
my_model.set_hop(0.2367, 58, 93, [0, 0, 0])
my_model.set_hop(3.0981, 59, 63, [0, 0, 0])
my_model.set_hop(0.3507, 59, 83, [0, 0, 0])
my_model.set_hop(0.2367, 59, 92, [0, 0, 0])
my_model.set_hop(0.2367, 60, 75, [0, 0, 0])
my_model.set_hop(0.3507, 60, 84, [0, 0, 0])
my_model.set_hop(0.2367, 61, 74, [0, 0, 0])
my_model.set_hop(0.2367, 62, 73, [0, 0, 0])
my_model.set_hop(0.3507, 63, 87, [0, 0, 0])
my_model.set_hop(0.2367, 63, 95, [0, 0, 0])
my_model.set_hop(3.0981, 64, 66, [0, 0, 0])
my_model.set_hop(0.2367, 64, 83, [0, 0, 0])
my_model.set_hop(3.0981, 65, 67, [0, 0, 0])
my_model.set_hop(0.2367, 65, 82, [0, 0, 0])
my_model.set_hop(0.3507, 65, 89, [0, 0, 0])
my_model.set_hop(0.2367, 66, 81, [0, 0, 0])
my_model.set_hop(0.2367, 67, 80, [0, 0, 0])
my_model.set_hop(0.3507, 67, 91, [0, 0, 0])
my_model.set_hop(3.0981, 68, 69, [0, 0, 0])
my_model.set_hop(0.2367, 68, 79, [0, 0, 0])
my_model.set_hop(0.2367, 69, 78, [0, 0, 0])
my_model.set_hop(0.2367, 70, 77, [0, 0, 0])
my_model.set_hop(0.3507, 70, 94, [0, 0, 0])
my_model.set_hop(0.2367, 71, 91, [0, 0, 0])
my_model.set_hop(0.3507, 71, 95, [0, 0, 0])
my_model.set_hop(2.8059, 72, 76, [0, 0, 0])
my_model.set_hop(2.8059, 72, 80, [0, 0, 0])
my_model.set_hop(0.3668, 72, 87, [0, 0, 0])
my_model.set_hop(0.3668, 72, 91, [0, 0, 0])
my_model.set_hop(0.2366, 72, 95, [0, 0, 0])
my_model.set_hop(2.8059, 73, 90, [0, 0, 0])
my_model.set_hop(2.8059, 75, 84, [0, 0, 0])
my_model.set_hop(0.3668, 76, 80, [0, 0, 0])
my_model.set_hop(2.8059, 76, 87, [0, 0, 0])
my_model.set_hop(0.2366, 76, 91, [0, 0, 0])
my_model.set_hop(0.3668, 76, 95, [0, 0, 0])
my_model.set_hop(0.4113, 77, 16, [0, 0, 1])
my_model.set_hop(0.2366, 77, 39, [0, 0, 1])
my_model.set_hop(0.3507, 77, 53, [0, 0, 1])
my_model.set_hop(2.8059, 77, 86, [0, 0, 1])
my_model.set_hop(2.8059, 77, 94, [0, 0, 0])
my_model.set_hop(2.8059, 78, 85, [0, 0, 0])
my_model.set_hop(0.4113, 79, 14, [0, 0, 1])
my_model.set_hop(0.2366, 79, 37, [0, 0, 1])
my_model.set_hop(0.3507, 79, 55, [0, 0, 1])
my_model.set_hop(0.3668, 79, 75, [0, 0, 1])
my_model.set_hop(2.8059, 79, 84, [0, 0, 1])
my_model.set_hop(0.2366, 80, 87, [0, 0, 0])
my_model.set_hop(2.8059, 80, 91, [0, 0, 0])
my_model.set_hop(0.3668, 80, 95, [0, 0, 0])
my_model.set_hop(0.4113, 82, 23, [0, 1, 0])
my_model.set_hop(0.2366, 82, 46, [0, 1, 0])
my_model.set_hop(0.3507, 82, 58, [0, 1, 0])
my_model.set_hop(2.8059, 82, 93, [0, 1, 0])
my_model.set_hop(2.8059, 82, 89, [0, 0, 0])
my_model.set_hop(2.8059, 83, 92, [0, 0, 0])
my_model.set_hop(0.4113, 86, 3, [0, 1, 0])
my_model.set_hop(0.2366, 86, 27, [0, 1, 0])
my_model.set_hop(0.3507, 86, 62, [0, 1, 0])
my_model.set_hop(2.8059, 86, 73, [0, 1, 0])
my_model.set_hop(0.3668, 86, 90, [0, 1, 0])
my_model.set_hop(0.3668, 87, 91, [0, 0, 0])
my_model.set_hop(2.8059, 87, 95, [0, 0, 0])
my_model.set_hop(2.8059, 91, 95, [0, 0, 0])
my_model.set_hop(0.4113, 92, 9, [1, 0, 0])
my_model.set_hop(0.2366, 92, 32, [1, 0, 0])
my_model.set_hop(0.3507, 92, 68, [1, 0, 0])
my_model.set_hop(2.8059, 92, 79, [1, 0, 0])
my_model.set_hop(0.4113, 93, 8, [1, 0, 0])
my_model.set_hop(0.2366, 93, 31, [1, 0, 0])
my_model.set_hop(0.3507, 93, 69, [1, 0, 0])
my_model.set_hop(2.8059, 93, 78, [1, 0, 0])
my_model.set_hop(0.3668, 93, 85, [1, 0, 0])

my_model.display()

path = [ [0., 0., 0.], [-0.5, 0.5, 0.5], [0., .5, 0.], [.25, .25, .25], [0., 0., 0.], [0., 0.5, 0.] ]
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

ax.set_title('bccP8-1 band structure')
ax.set_xlabel('Path in k-space')
ax.set_ylabel('Band energy')

for xband in range(96):
	ax.plot(k_dist, evals[xband], linewidth=.5, color='black')

# fig.tight_layout()
fig.savefig('bccP8-1.eps')

print('Done.')