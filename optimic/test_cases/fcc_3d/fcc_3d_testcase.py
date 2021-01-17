# -*- coding: utf-8 -*-
"""
fcc_3d_testcase.py

Module to test some features of program for FCC 3D type spacing.

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

from src.set_logger import set_logger 
name_str = __name__

def fcc_3d_testcase(tessellation, dimension, size_of_simulation_box, \
    spacing_length, grain_size_distributions, number_of_neighbor, \
    grain_boundary_area_distribution, junction_lengths, \
    junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, \
    disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level):
    """
    Execute all assert statements related to fcc_3d_testcase.

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

    dimension: integer    
        Dimension of study (2 or 3)

    size_of_simulation_box: array of length 3
        Size of simulation box (array of length along X, Y, Z directions)

    spacing_length: float
        Spacing between seeds along X, Y & Z in 3D case and along X & Y
        directions in Quasi-2D case. Also the spacing_length must be a perfect
        divisor of size of simulation box along all three directions.

    grain_size_distribution: 2D array
        An array of all the grain sizes in terms of radius. The
        first column is the grain number and second column consists of grain 
        sizes.

    number_of_neighbor: 2D array
        An array comprising of data related to the number of 
        neighbors of each grains. Column names are: Grain number, number of 
        neighbors, grain indexes of neighbors.

    grain_boundary_area_distribution_ 2D array
        An array of consisting of columns Sr. No., Grain 1, 
        Grain 2 and grain boundary area.

    junction_lengths: List of lists
        List of lists consisting of column names: Sr. no., Junction type, 
        Junction lengths, grains with this junction.

    junction_angles_degrees: List of lists
        List of lists consisting of column names: Sr. No., Junction Type, 
        1st Junction angle, Grain containing the 1st junction angle, 
        2nd Junction angle, so on...

    distance_btw_grain_array: 2D array
        Symmetric array with rows and columns represented by grain numbers in 
        ascending order and each element of the array representing distance 
        between respective grain numbers.

    distance_btw_grain_1d: 1D array
        1D array of distances between grains

    disorientation_angle: 2D array
        An array consisting of columns grain 1, grain 2, disorientation 
        angle, disorientation axis.

    schmid_factors: 2D array
        An array consisting of Schmid Factor related data with column names as 
        Grain no., Schmid factor, Slip plane and direction (Slip System).

    type_of_grain_boundaries: 2D array
        An array with the column names as grain 1, grain 2, rotation angle, 
        rotation axis, type of csl (0 indicates normal grain boundary).

    log_level: string
        Logger level to be used.

    Returns
    -------
    Function returns nothing.

    """

    log = set_logger(name_str, 'log_data.log', log_level)
    try:
        for v in range(copy.deepcopy(tessellation['number_of_grains'])):
            assert np.isclose(copy.deepcopy(tessellation['number_of_faces_list'][v]), 12)
        assert np.all([grain_size_distributions[:, 1] == grain_size_distributions[0, 1]])     # The grain sizes should be the same
        assert np.all([np.around(grain_size_distributions[:, 1], decimals=2) == np.around(0.7816*spacing_length, decimals=2)])     # The grain sizes should be the same
        assert np.all([np.around(grain_size_distributions[:, 1], decimals=2) == np.around(np.cbrt((size_of_simulation_box**3 * 6)/(np.pi * copy.deepcopy(tessellation['number_of_grains'])) ), decimals=2)])                               # The grain sizes should be the same
        assert np.all([np.isclose(num[2], number_of_neighbor[0][2]) for num in number_of_neighbor])   # The no. of neighbors should be the same 
        assert np.all([np.isclose(num[2], 12) for num in number_of_neighbor])   # The no. of neighbors should be the same               
        assert np.all([np.isclose(area[3], grain_boundary_area_distribution[0][3]) for area in grain_boundary_area_distribution])   # The GB areas should be the same
        assert np.all([np.isclose(area[3], 0.1768*(spacing_length**2), atol=1e-2) for area in grain_boundary_area_distribution])   # The GB areas should be the same
        assert np.all([np.isclose(length[2], junction_lengths[0][2]) for length in junction_lengths])   # The junctions length should be the same
        assert np.all([np.isclose(length[2], 0.4330*spacing_length, atol=1e-2) for length in junction_lengths])   # The junctions length should be the same
        assert np.all(np.array([np.allclose(angles[2::2], 120.0) for angles in junction_angles_degrees]))           # All angles should be the same
    except AssertionError:
        log.exception('fcc_3d_testcase failed !!')

    log.info('fcc_3d_testcase passed !!')