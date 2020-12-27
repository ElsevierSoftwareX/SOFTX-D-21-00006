# -*- coding: utf-8 -*-
"""
generate_distribution.py

Module to generate lognormal, Gauss normal and gamma distribution.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

@authors: Serrao Prince Henry, Arun Prakash
@email: 
created: 24 January 2020
Copyright © 2020 by Serrao Prince Henry, Dr. Arun Prakash

This file is part of Optimized Micro-structure Generator.

Optimized Micro-structure Generator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Optimized Micro-structure Generator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Optimized Micro-structure Generator.  If not, see <https://www.gnu.org/licenses/>.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm
from scipy.stats import norm
from scipy.stats import gamma

## Parameters for distribution
s = 0.945
a = 2
loc = 40
scale = 4
x = np.linspace(0, 10, 1000)
x_disorientation = np.linspace(0, 60, 1000)

## Generating distribution
distribution_lognormal = lognorm.pdf(x, s, loc=1)
distribution_gaussian = norm.pdf(x_disorientation, loc = loc, scale = scale)
distribution_gamma = gamma.pdf(x, a)

## Storing lognorm distribution to an array
data_array = np.zeros([len(x), 2])
data_array[:, 0] = np.around(x[:], decimals=4)
data_array[:, 1] = np.around(distribution_lognormal, decimals=4)

## Storing gaussian distribution to array
data_array_gauss = np.zeros([len(x), 2])
data_array_gauss[:, 0] = np.around(x_disorientation[:], decimals=4)
data_array_gauss[:, 1] = np.around(distribution_gaussian, decimals=4)

## Storing gamma distribution to array
data_array_gamma = np.zeros([len(x), 2])
data_array_gamma[:, 0] = np.around(x[:], decimals=4)
data_array_gamma[:, 1] = np.around(distribution_gamma, decimals=4)

## Checking if the total probability (area under curve) is 1
print((np.sum((data_array[1:, 0] - data_array[:-1, 0]) * data_array[:-1, 1])))
assert np.isclose((np.sum((data_array[1:, 0] - data_array[:-1, 0]) * data_array[:-1, 1])), 1, atol=1e-1) 
assert np.isclose((np.sum((data_array_gauss[1:, 0] - data_array_gauss[:-1, 0]) * data_array_gauss[:-1, 1])), 1, atol=1e-1)
assert np.isclose((np.sum((data_array_gamma[1:, 0] - data_array_gamma[:-1, 0]) * data_array_gamma[:-1, 1])), 1, atol=1e-1)


## Writing the generated distribution to file
with open('user_distribution.txt', 'a+') as f:
    f.truncate(0)
    f.write('# Grain Sizes, Frequency of Grain Sizes \n')
    np.savetxt(f, data_array, delimiter=',', fmt= '%.4f')
    f.write('end\n')
    f.write('# Disorientation_angles, Frequency of disorientation angles \n')
    np.savetxt(f, data_array_gauss, delimiter=',', fmt= '%.4f')
    f.write('end\n')
    f.write('# Distance between grains, Frequency of distances \n')
    np.savetxt(f, data_array_gamma, delimiter=',', fmt= '%.4f')
    f.write('end\n')

## Plotting
fig, ax = plt.subplots()
ax.plot(x, distribution_lognormal,
       'r-', lw=5, alpha=0.6, label='lognorm pdf')
ax.plot(x_disorientation, distribution_gaussian,
       'y-', lw=5, alpha=0.6, label='gaussian pdf')
ax.plot(x, distribution_gamma,
       'b-', lw=5, alpha=0.6, label='gamma pdf')

plt.legend()
plt.show()