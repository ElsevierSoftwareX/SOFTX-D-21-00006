# -*- coding: utf-8 -*-
"""
hcp_lattice_3D.py

Module to generate seed in HCP lattice type of spacing in 3D.

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

def hcp_lattice_3D(limit, spacing_length, log_level):
    """
    Generate seeds with HCP lattice type spacing in 3D.

    Processing
    ----------
    The function generates seeds in HCP lattice spacing in 3D and ensures 
    that there is no overlapping of seeds at the extreme ends of simulation box.

    Parameters
    ----------
    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

    spacing_length: float
        Spacing between seeds along X direction. Spacing_length is a factor 
        multiplied to SQRT(3) in Y & Z directions in 3D case and along Y 
        directions in Quasi-2D case.

    log_level: string
        Logger level to be used.

    Returns
    ------- 
    The function returns an array of regularly spaced seeds and also ensures 
    that the seeds are unique.
    """

    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Starting to generate HCP 3D type seed spacing')

    x_limit, y_limit, z_limit = limit 
    seed_array = np.zeros([int(np.ceil(x_limit/spacing_length) * (2*np.ceil(y_limit/(np.sqrt(3)*spacing_length))) * (2*np.ceil(z_limit/(np.sqrt(3)*spacing_length)))), 3]) # sqrt(3) is due to spacing of seeds in HCP

    counter = 0
    for x in np.arange(0, x_limit, spacing_length):
        for y in np.arange(0, y_limit, np.sqrt(3)*spacing_length):
            for z in np.arange(0, z_limit, np.sqrt(3)*spacing_length):
                if (x >= x_limit) or (y >= y_limit) or (y >= z_limit):
                    continue
                seed_array[counter, :] = x, y, z
                counter += 1

                seed_array[counter, 0] = x + (spacing_length/2)
                seed_array[counter, 1] = y + (np.sqrt(3)*spacing_length/2)
                seed_array[counter, 2] = z + (np.sqrt(3)*spacing_length/2)
                counter += 1

    new_array = [tuple(row) for row in seed_array]                          # generating tuple of each row of the seeds array
    seed_array_unique = np.unique(new_array, axis = 0)                      # removing duplicate rows to avoid overlapping seeds

    log.debug('Generated seeds array is ' + str(seed_array_unique))
    log.info('HCP 3D type seed array was successfully generated')

    return seed_array_unique
