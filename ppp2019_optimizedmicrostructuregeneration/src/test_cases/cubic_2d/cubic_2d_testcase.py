# -*- coding: utf-8 -*-
"""
cubic_2d_testcase.py

Module to test some features of program for CUBIC 2D type spacing.

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

def cubic_2d_testcase(tessellation, dimension, size_of_simulation_box, \
    spacing_length, length_z, grain_size_distributions, number_of_neighbor, \
    grain_boundary_area_distribution, junction_lengths, \
    junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, \
    disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level):

    log = set_logger(name_str, 'log_data.log', log_level)
    try:
        for v in range(copy.deepcopy(tessellation['number_of_grains'])):
            assert np.isclose(copy.deepcopy(tessellation['number_of_faces_list'][v]) - dimension, 4)                           # ignoring surfaces with surface normals along z axis. REMEMBER !!!!!!
        assert np.all([grain_size_distributions[:, 1] == grain_size_distributions[0, 1]])     # The grain sizes should be the same
        assert np.allclose(grain_size_distributions[:, 1], 1.1284*spacing_length, atol=1e-2)	# Please refer documentation for details regarding 1.1284
        assert np.all([np.around(grain_size_distributions[:, 1], decimals=2) == np.around(np.sqrt((size_of_simulation_box**2 * 4)/(np.pi * copy.deepcopy(tessellation['number_of_grains'])) ), decimals=2)])                               # The grain sizes should be the same
        assert np.all([np.isclose(num[1], number_of_neighbor[0][1]) for num in number_of_neighbor])   # The no. of neighbors should be the same                
        assert np.all([np.isclose(num[1], 4) for num in number_of_neighbor])
        assert np.all([np.isclose(area[3], grain_boundary_area_distribution[0][3]) for area in grain_boundary_area_distribution])
        assert np.all([np.isclose(area[3], length_z*spacing_length) for area in grain_boundary_area_distribution])
        #assert np.all([length[2] == junction_lengths[0][2] for length in junction_lengths])   # The junctions length should be the same
        assert np.all([((np.around(length[2], decimals=2) == np.around(1 * spacing_length, decimals=2)) or (length[2] == length_z)) for length in junction_lengths])
        assert np.all(np.concatenate(np.array([np.equal(angles[2::2], 90.0) for angles in junction_angles_degrees])).flatten())           # All angles should be the same
    except AssertionError:
        log.exception('cubic_2d_testcase failed !!')

    log.info('cubic_2d_testcase passed !!')