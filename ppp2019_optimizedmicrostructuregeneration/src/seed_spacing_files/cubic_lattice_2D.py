# -*- coding: utf-8 -*-
"""
cubic_lattice_2D.py

Module to generate seed in Cubic lattice type of spacing in 2D.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

@authors: Serrao Prince Henry, Arun Prakash
@email: 
created: 04 January 2020
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

from ppp2019_optimizedmicrostructuregeneration.src.seed_spacing_files.import_statements import *

from ppp2019_optimizedmicrostructuregeneration.src.set_logger import set_logger as set_logger
name_str = __name__

def cubic_lattice_2D(limit, a, log_level):
    """
    Input: The function requires an array specifying the limits along X, Y, Z 
            directions respectively and the seed spacing length (lattice constant).
    
    Process: The function generates seeds in cubic lattice spacing in 2D and 
            ensures that there is no overlapping of seeds at the extreme ends of 
            simulation box.

    Returns: The function returns an array of regularly spaced seeds and also 
            ensures that the seeds are unique.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Starting to generate Cubic 2D type seed spacing')

    seed_array = np.zeros([(int(limit[0])) * (int(limit[1])), 3])
    seed_array[:, 2] = limit[2]												# to initiate first seed with same Z coordinate

    counter = 0
    for x in np.linspace(0, limit[0], limit[0]/a + 1):
        for y in np.linspace(0, limit[1], limit[1]/a + 1):
                if x >= limit[0] or y >= limit[1]:
                    continue
                seed_array[counter, 0] = x
                seed_array[counter, 1] = y
                seed_array[counter, 2] = limit[2]
                counter += 1

    new_array = [tuple(row) for row in seed_array]                          # generating tuple of each row of the seeds array
    seed_array_unique = np.unique(new_array, axis = 0)                      # removing duplicate rows to avoid overlapping seeds

    log.debug('Generated seeds array is ' + str(seed_array_unique))
    log.info('Cubic 2D type seed array was successfully generated')

    return seed_array_unique
