# -*- coding: utf-8 -*-
"""
bcc_lattice_3D.py

Module to generate seed in Body Centered Cubic (BCC) type of spacing.

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

def bcc_lattice_3D(limit, a, log_level):
    """
    Generate seeds with BCC lattice type spacing in 3D.

    Processing
    ----------
    The function generates seeds in BCC spacing and ensures that there is no 
    overlapping of seeds at the extreme ends of simulation box.

    Parameters
    ----------
    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

    a: float
        Spacing between seeds along X, Y & Z in 3D case and along X & Y
        directions in Quasi-2D case. Also the spacing_length must be a perfect
        divisor of size of simulation box along all three directions.

    log_level: string
        Logger level to be used.

    Returns
    ------- 
    The function returns an array of regularly spaced seeds and also ensures 
    that the seeds are unique.
    """

    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Starting to generate BCC type seed spacing')
    
    ## Generating an empty list to store seed coordinates
    seed_array = []

    ## Generating regular grid seeds array
    ## Iterating in all directions
    for x in np.arange(0, limit[0], a):
        for y in np.arange(0, limit[1], a):
            for z in np.arange(0, limit[2], a):
                if x >= limit[0] or y >= limit[1] or z >= limit[2]:       # Checking if the current position is outside limits
                    continue
                seed_array.append([x, y, z])

                if x < limit[0] and y < limit[1] and z < limit[2]:        # Placing the body centered seed
                    seed_array.append([x + a/2, y + a/2, z + a/2])

    seed_array = np.array(seed_array)
    ## Checking for uniqueness of all seeds
    new_array = [tuple(row) for row in seed_array]
    seed_array_unique = np.unique(new_array, axis = 0)

    log.debug('Generated seeds array is ' + str(seed_array_unique))
    log.info('BCC type seed array was successfully generated')

    return seed_array_unique