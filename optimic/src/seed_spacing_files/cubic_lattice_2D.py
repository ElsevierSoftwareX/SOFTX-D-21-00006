# -*- coding: utf-8 -*-
"""
cubic_lattice_2D.py

Module to generate seed in Cubic lattice type of spacing in 2D.

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

def cubic_lattice_2D(limit, a, log_level):
    """
    Generate seeds with Cubic lattice type spacing in 2D.

    Processing
    ----------
    The function generates seeds in cubic lattice spacing in 2D and ensures 
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
    log.info('Starting to generate Cubic 2D type seed spacing')

    seed_array = []

    for x in np.arange(0, limit[0], a):
        for y in np.arange(0, limit[1], a):
                if x >= limit[0] or y >= limit[1]:
                    continue
                seed_array.append([x, y, limit[2]/2.0])

    seed_array = np.array(seed_array)
    new_array = [tuple(row) for row in seed_array]                          # generating tuple of each row of the seeds array
    seed_array_unique = np.unique(new_array, axis = 0)                      # removing duplicate rows to avoid overlapping seeds

    log.debug('Generated seeds array is ' + str(seed_array_unique))
    log.info('Cubic 2D type seed array was successfully generated')

    return seed_array_unique