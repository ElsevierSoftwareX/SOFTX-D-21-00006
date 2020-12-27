# -*- coding: utf-8 -*-
"""
one_seed_testcase.py

Module to test some features of program for one seed.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
created: 31 March 2020
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

from src.main_import_statements import *

from src.set_logger import set_logger as set_logger
name_str = __name__

def one_seed_testcase(tessellation, cell_number, size_of_simulation_box, length_z, log_level):
    """
    Execute all assert statements related to one_seed_testcase.

    Parameters
    ----------
    tessellation: dictionary
        Dictionary of tessellations data with following keys:
            1. number_of_grains
            2. number_of_faces_list
            3. vertices_list
            4. face_vertices_list
            5. centroid_list
            6. volume_list
            7. normals_list
            8. neighbors_list
            9. face_area_list
            10. number_of_edges_list

    cell_number: integer
        Grain number.(Since only 1 grain, it must be always 0)

    size_of_simulation_box: array of length 3
        Size of simulation box (array of length along X, Y, Z directions)

    length_z: float
        Size of simulation box along Z axis.

    log_level: string
        Logger level to be used.

    Returns
    -------
    Function returns nothing.
    """

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