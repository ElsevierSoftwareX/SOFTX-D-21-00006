# -*- coding: utf-8 -*-
"""
random_generator.py

Module to generate random seed coordinates.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
created: 16 November 2019
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

from optimic.src.main_import_statements import *

from optimic.src.set_logger import set_logger 

name_str = __name__

def random_generator(number_of_seeds, dimension, limits, log_level):       
    """
    Generates random numbers within the ranges from 0 to respective length of 
    simulation box along specific axis.

    Parameters
    ----------
    number_of_seeds: integer
        Number of seeds/grains

    dimension: integer 
        Dimension of study. (2 or 3)

    limits: array
        Size of simulation box (array of length along X, Y, Z directions)
    
    log_level: string
        Logger level to be used.

    Returns
    -------
    The function returns a list containing the coordinates of each seed.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    ## Generating random numbers based on dimensions (3D or 2D) of simulation box
    seeds = []
    length_x_axis, length_y_axis, length_z_axis = limits[:]
    
    for i in range(number_of_seeds):
               
        if dimension == 2:
            rand_x = np.around(np.random.uniform(low = 0, high = length_x_axis), decimals=4)
            rand_y = np.around(np.random.uniform(low = 0, high = length_y_axis), decimals=4)
            random_number_list = [rand_x, rand_y, limits[2]/2.0]
            seeds.append(random_number_list)
        
        elif dimension == 3:
            rand_x = np.around(np.random.uniform(low = 0, high = length_x_axis), decimals=4)
            rand_y = np.around(np.random.uniform(low = 0, high = length_y_axis), decimals=4)
            rand_z = np.around(np.random.uniform(low = 0, high = length_z_axis), decimals=4)
            random_number_list = [rand_x, rand_y, rand_z]
            seeds.append(random_number_list)

    log.info('Successfully generated random seeds')
    return seeds