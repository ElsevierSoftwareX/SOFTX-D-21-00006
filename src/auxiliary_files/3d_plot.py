# -*- coding: utf-8 -*-
"""
3d_plot.py

Module to generate 3d plot. It has to be modified based on needs

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

@authors: Serrao Prince Henry, Arun Prakash
@email: 
created: 19 March 2020
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
from mpl_toolkits.mplot3d import Axes3D

size_of_simulation_box = 2.6
a = 1.3
limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
seed_array = np.zeros([2*((int(limit[0])) * (int(limit[1])) + (int(limit[1]) * int(limit[2])) + (int(limit[0]) * int(limit[2])) + (int(limit[0])*int(limit[1])*int(limit[2])*3)), 3])
dimension = 3
material = "FCC_3D"

seed_array = np.zeros([2*((int(limit[0])) * (int(limit[1])) + (int(limit[1]) * int(limit[2])) + (int(limit[0]) * int(limit[2])) + (int(limit[0])*int(limit[1])*int(limit[2])*3)), 3])

## Generating regular grid seeds array
counter = 0
for x in np.linspace(0, limit[0], limit[0]/a +1):
    for y in np.linspace(0, limit[1], limit[1]/a +1):
        for z in np.linspace(0, limit[2], limit[2]/a +1):
            
            if a*x >= limit[0] or a*y >= limit[1] or a*z >= limit[2]:
                continue
            seed_array[counter, 0] = x
            seed_array[counter, 1] = y
            seed_array[counter, 2] = z
            counter += 1

            if x + a/2 <= limit[0] and y + a/2 <= limit[1]:
                seed_array[counter, 0] = x + a/2
                seed_array[counter, 1] = y + a/2
                seed_array[counter, 2] = z
                counter += 1

            if y + a/2 <= limit[1] and z + a/2 <= limit[2]:
                seed_array[counter, 0] = x
                seed_array[counter, 1] = y + a/2
                seed_array[counter, 2] = z + a/2
                counter += 1

            if x + a/2 <= limit[0] and z + a/2 <= limit[2]:
                seed_array[counter, 0] = x + a/2
                seed_array[counter, 1] = y
                seed_array[counter, 2] = z + a/2
                counter += 1 

new_array = [tuple(row) for row in seed_array]
seed_array_unique = np.unique(new_array, axis = 0)

fig = plt.figure(figsize=(10,10))
ax = Axes3D(fig)

ax.scatter(seed_array_unique[:, 0], seed_array_unique[:, 1], seed_array_unique[:, 2])

fig, ax = plt.subplots()
ax.scatter(seed_array_unique[:, 0], seed_array_unique[:, 1])
plt.show()