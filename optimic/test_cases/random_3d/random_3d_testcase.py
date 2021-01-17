# -*- coding: utf-8 -*-
"""
random_3d_testcase.py

Module to test some features of program for RANDOM 3D type spacing.

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

from src.textural_characteristics import sharp_texture_quaternions 
from src.textural_characteristics import random_quaternions_generator 
from src.textural_characteristics import disorientation_angles 
from src.textural_characteristics import type_of_grain_boundary 
from src.textural_characteristics import schmid_factor 
from src.textural_characteristics import available_required_texture 

def random_3d_testcase(store_folder, version, now, material, tessellation, \
    dimension, size_of_simulation_box, spacing_length, seed_array_unique, \
    required_texture, number_of_seeds, grain_size_distributions, \
    number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
    junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, \
    disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level):
    """
    Execute all test conditions related to random_3d_testcase.

    Parameters
    ----------
    store_folder: string
        Name of directory where output files are to be stored.

    version: string
        Current version of code.

    now: string
        Current time and date.

    material: string
        Material name of which microstructure is being tested.

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

    seed_array_unique: array of shape (number of grains, 3)
        Unique seed coordinates

    required_texture: array of length 3
        Specific texture to be used for all grains

    number_of_seeds: integer
        Number of seeds/grains.

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

    Output
    ------
        1. Seeds_data_without_orientations.txt in current working directory.
        2. Seeds_data_with_orientations.txt in current working directory.
        3. array_verification_file.txt in appropriate directory within specified 
            storage_folder

    """
    log = set_logger(name_str, 'log_data.log', log_level)
    ##Writing Seeds information without orientations to a file
    with open('seeds_data_without_orientations.txt', 'a+') as f:
        f.truncate(0)
        f.write("# X, Y, Z, W, OX, OY, OZ \n")
        np.savetxt(f, seed_array_unique, fmt='%.4f', delimiter= ',')

    ## Writing Seeds information with orientations to a file
    orientation_quaternion = sharp_texture_quaternions(1, required_texture, log_level)                       # Generating a random quaternion
    #orientation_quaternion = np.array([1, 0, 0, 0])
    orientation_array = np.broadcast_to(orientation_quaternion, (seed_array_unique.shape[0], 4))   # Assigning same orientation to each grain
    seeds_data = np.concatenate((seed_array_unique, orientation_array), axis=1)         # Concatenating seed coorinates and orientations

    with open('seeds_data_with_orientations.txt', 'a+') as f:
        f.truncate(0)
        f.write("# X, Y, Z, W, OX, OY, OZ \n")
        np.savetxt(f, seeds_data, fmt='%.4f', delimiter= ',')
    
    ## Reading Seeds data from the files
    with open('seeds_data_without_orientations.txt', 'r') as f:
        data = np.loadtxt(f, delimiter=',', comments='#')
        
        ## Extracting orientations data if available
        if data.shape[1] == 7:
            orientation_data = data[:, 3:7]
        else: orientation_data = None
    
    with open('seeds_data_with_orientations.txt', 'r') as f:
        data_orientation = np.loadtxt(f, delimiter=',', comments='#')
        seeds_data = data_orientation[:, :3]
        
        ## Extracting orientations data if available
        if data_orientation.shape[1] == 7:
            orientation_data = data_orientation[:, 3:7]
        else: orientation_data = None
    
    total_volume = 0
    for v in range(copy.deepcopy(tessellation['number_of_grains'])):
        total_volume += copy.deepcopy(tessellation['volume_list'][v]) 
    
    try:
        assert seed_array_unique.shape[0] == number_of_seeds                ## Checking if the number of random unique seeds are same as required
        assert np.isclose(total_volume, size_of_simulation_box*size_of_simulation_box*size_of_simulation_box)
        assert np.array_equal(seed_array_unique, data)
        assert np.array_equal(seed_array_unique, seeds_data)
        assert np.allclose(orientation_data, orientation_array, atol=1e-4)
    except AssertionError:
        log.exception('random_3d_testcase failed !!')

    ## Saving to array_verification_file.txt 
    output_file_path = Path(store_folder, material, now, "Text_output", "array_verification_file.txt")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    with open(str(output_file_path), 'a+') as f:
        f.truncate(0)
        f.write("Verification matrix for seed_array read from file without orientations \n")
        np.savetxt(f, seed_array_unique - data, newline='\n', delimiter=',', fmt="%.4f")
        f.write("\n \n Verification matrix for seed_array read from file with orientations \n")
        np.savetxt(f, seed_array_unique - seeds_data, newline='\n', delimiter=',', fmt="%.4f")
        f.write("\n \n Verification matrix for orientations array read from file with orientations \n")
        np.savetxt(f, orientation_data - orientation_array, newline='\n', delimiter=',', fmt="%.4f")

    log.info('random_3d_testcase passed !!')
