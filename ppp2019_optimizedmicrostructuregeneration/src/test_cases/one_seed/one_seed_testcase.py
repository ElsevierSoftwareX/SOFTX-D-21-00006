# -*- coding: utf-8 -*-
"""
one_seed_testcase.py

Module to test some features of program for one seed.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

@authors: Serrao Prince Henry, Arun Prakash
@email: 
created: 31 March 2020
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

from ppp2019_optimizedmicrostructuregeneration.src.main_import_statements import *

from ppp2019_optimizedmicrostructuregeneration.src.set_logger import set_logger as set_logger
name_str = __name__

def one_seed_testcase(tessellation, cell_number, size_of_simulation_box, length_z, log_level):

    log = set_logger(name_str, 'log_data.log', log_level)
    total_volume = 0                                                        # initializing the variable
    neighbors = tessellation['neighbors_list'][cell_number]                       # storing the neighbors of the specified cell
    for v in range(copy.deepcopy(tessellation['number_of_grains'])):
        total_volume += copy.deepcopy(tessellation['volume_list'][v])                                          # updating the total_volume variable
    
    try:
        assert np.isclose(total_volume, size_of_simulation_box*size_of_simulation_box*length_z)  ## Testing for total volume of all cells                         
        assert np.isclose(np.count_nonzero(neighbors), 0)                       ## Testing periodicity
    except AssertionError:
        log.exception('one_seed_testcase failed !!')

    log.info('one_seed_testcase passed !!')