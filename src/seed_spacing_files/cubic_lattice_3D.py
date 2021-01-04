# -*- coding: utf-8 -*-
"""
cubic_lattice_3D.py

Module to generate seed in Cubic lattice type of spacing in 3D.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
created: 04 January 2020
Copyright © 2020 by Serrao Prince Henry, Dr. Arun Prakash

This file is part of OptiMic.

OptiMic is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OptiMic is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OptiMic.  If not, see <https://www.gnu.org/licenses/>.

"""

from src.seed_spacing_files.import_statements import *

from src.set_logger import set_logger 
name_str = __name__

def cubic_lattice_3D(limit, a, log_level):
    """
    Generate seeds with Cubic lattice type spacing in 3D.

    Processing
    ----------
    The function generates seeds in cubic lattice spacing in 3D and ensures 
    that there is no overlapping of seeds at the extreme ends of simulation box.

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
    log.info('Starting to generate Cubic 3D type seed spacing')

    ## Generating an array to store seed coordinates
    seed_array = np.zeros([(int(np.ceil(limit[0]/a)) + 1) * (int(np.ceil(limit[1]/a)) + 1) * (int(np.ceil(limit[2]/a)) + 1), 3])
    counter = 0
    
    ## Iterating in all directions
    for x in np.arange(0, limit[0], a):
        for y in np.arange(0, limit[1], a):
            for z in np.arange(0, limit[2], a):
                
                if x >= limit[0] or y >= limit[1] or z >= limit[2]:
                    continue
                seed_array[counter, 0] = x
                seed_array[counter, 1] = y
                seed_array[counter, 2] = z
                counter += 1
    
    ## Checking for uniqueness of all seeds
    new_array = [tuple(row) for row in seed_array]
    seed_array_unique = np.unique(new_array, axis = 0)

    log.debug('Generated seeds array is ' + str(seed_array_unique))
    log.info('Cubic 3D type seed array was successfully generated')
    
    return seed_array_unique