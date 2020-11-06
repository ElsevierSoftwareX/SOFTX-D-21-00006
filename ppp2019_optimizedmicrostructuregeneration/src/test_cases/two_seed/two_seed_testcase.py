# -*- coding: utf-8 -*-
"""
two_seed_testcase.py

Module to test some features of program for TWO SEEDS.

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

from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import disorientation_angles as disorientation_angles
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import type_of_grain_boundary as type_of_grain_boundary

from ppp2019_optimizedmicrostructuregeneration.src.set_logger import set_logger as set_logger
name_str = __name__

def two_seed_testcase(store_folder, version, now, material, tessellation, \
    dimension, size_of_simulation_box, spacing_length, length_z, \
    required_texture, rand_quat_flag, grain_size_distributions, \
    number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
    junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, \
    disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level):
    """
    Execute all test conditions related to two_seed_testcase.

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

    size_of_simulation_box: float
        Size of simulation box (array of length along X, Y, Z directions)

    spacing_length: float
        Spacing between seeds along X, Y & Z in 3D case and along X & Y
        directions in Quasi-2D case. Also the spacing_length must be a perfect
        divisor of size of simulation box along all three directions.

    length_z: float
        Size of simulation box along Z axis.

    required_texture: array of length 3
        Specific texture to be used for all grains

    rand_quat_flag: boolean
        Flag to indicate that random orientations are to be assigned to all 
        grains.

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
    array_verification_file.txt in appropriate directory within specified 
    storage_folder
    """
    log = set_logger(name_str, 'log_data.log', log_level)
    limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
    
    total_volume = 0    
    assert np.isclose(copy.deepcopy(tessellation['volume_list'][0]) , (size_of_simulation_box*size_of_simulation_box*length_z)/2) ## Testing the volume of each cell
    assert np.isclose(copy.deepcopy(tessellation['volume_list'][1]) , (size_of_simulation_box*size_of_simulation_box*length_z)/2) ## Testing the volume of each cell
    for v in range(copy.deepcopy(tessellation['number_of_grains'])):
        total_volume += copy.deepcopy(tessellation['volume_list'][v]) 
    assert np.isclose(total_volume, size_of_simulation_box*size_of_simulation_box*length_z)                                    ## Testing the total volume of the cells
    assert np.isclose(distance_btw_grain_1d, size_of_simulation_box/2)

    #tessellation = create_tessellations(seed_array_unique, limit)
    
    # Σ = 15 grain boundary orientations
    orientation_data = np.array([[1, 0, 0, 0], [0.7745967, 0.5163978, 0.2581989, 0.2581989]])
    skewed_boundary_flag = False                               
    disorientation_angle, orientation_data = disorientation_angles(dimension, limit, skewed_boundary_flag, required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    type_of_grain_boundaries_sigma_15, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, dimension, limit, skewed_boundary_flag, log_level)
    assert (np.all(type_of_grain_boundaries_sigma_15[:, 6] == 15))
    assert (np.all(np.isclose(disorientation_angle[:, 2], 48.2, atol=1e-1)))

    ## Saving to array_verification_file.txt
    output_file_path = Path(store_folder, material, now, "Text_output", "array_verification_file.txt")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    with open(str(output_file_path), 'a+') as f:
        f.truncate(0)
        f.write("Verification matrix for disorientation matrix for Sigma 15 GB type \n")
        f.write("# Grain 1, Grain 2, Disorientation angle, Disorientation axis \n")
        np.savetxt(f, disorientation_angle, newline='\n', delimiter=',', fmt="%.4f")
        f.write("\n \n Verification matrix for Type of GB for Sigma 15 GB type \n")
        f.write("# Grain 1, Grain 2, Misorientation angle, Misorientation axis, Type of Grain Boundary \n")
        np.savetxt(f, type_of_grain_boundaries_sigma_15, newline='\n', delimiter=',', fmt="%.4f")
    
    # Σ = 5 grain boundary orientations
    orientation_data = np.array([[1, 0, 0, 0], [0.6710739, 0.6706129, 0.2235376, 0.2235376]])
    skewed_boundary_flag = False 
    disorientation_angle, orientation_data = disorientation_angles(dimension, limit, skewed_boundary_flag, required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    type_of_grain_boundaries_sigma_5, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, dimension, limit, skewed_boundary_flag, log_level)
    assert (np.all(type_of_grain_boundaries_sigma_5[:, 6] == 5))
    assert (np.all(np.isclose(disorientation_angle[:, 2], 36.9, atol=1e-1)))

    ## Saving to array_verification_file.txt
    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for disorientation matrix for Sigma 5 GB type \n")
        f.write("# Grain 1, Grain 2, Disorientation angle, Disorientation axis \n")
        np.savetxt(f, disorientation_angle, newline='\n', delimiter=',', fmt="%.4f")
        f.write("\n \n Verification matrix for Type of GB for Sigma 5 GB type \n")
        f.write("# Grain 1, Grain 2, Misorientation angle, Misorientation axis, Type of Grain Boundary \n")
        np.savetxt(f, type_of_grain_boundaries_sigma_5, newline='\n', delimiter=',', fmt="%.4f")

    # Σ = 3 grain boundary orientations
    orientation_data = np.array([[1, 0, 0, 0], [0.8166416, 0.4081033, 0.4081033, 0]])
    skewed_boundary_flag = False 
    disorientation_angle, orientation_data = disorientation_angles(dimension, limit, skewed_boundary_flag, required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    type_of_grain_boundaries_sigma_3, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, dimension, limit, skewed_boundary_flag, log_level)
    assert (np.all(type_of_grain_boundaries_sigma_3[:, 6] == 3))
    assert (np.all(np.isclose(disorientation_angle[:, 2], 60.0, atol=1e-1)))

    ## Saving to array_verification_file.txt
    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for disorientation matrix for Sigma 3 GB type \n")
        f.write("# Grain 1, Grain 2, Disorientation angle, Disorientation axis \n")
        np.savetxt(f, disorientation_angle, newline='\n', delimiter=',', fmt="%.4f")
        f.write("\n \n Verification matrix for Type of GB for Sigma 3 GB type \n")
        f.write("# Grain 1, Grain 2, Misorientation angle, Misorientation axis, Type of Grain Boundary \n")
        np.savetxt(f, type_of_grain_boundaries_sigma_3, newline='\n', delimiter=',', fmt="%.4f")

    # Σ = 7 grain boundary orientations
    orientation_data = np.array([[1, 0, 0, 0], [0.8017756, 0.5345322, 0.2672661, 0]])
    skewed_boundary_flag = False 
    disorientation_angle, orientation_data = disorientation_angles(dimension, limit, skewed_boundary_flag, required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    type_of_grain_boundaries_sigma_7, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, dimension, limit, skewed_boundary_flag, log_level)
    assert (np.all(type_of_grain_boundaries_sigma_7[:, 6] == 7))
    assert (np.all(np.isclose(disorientation_angle[:, 2], 38.2, atol=1e-1)))

    ## Saving to array_verification_file.txt
    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for disorientation matrix for Sigma 7 GB type \n")
        f.write("# Grain 1, Grain 2, Disorientation angle, Disorientation axis \n")
        np.savetxt(f, disorientation_angle, newline='\n', delimiter=',', fmt="%.4f")
        f.write("\n \n Verification matrix for Type of GB for Sigma 7 GB type \n")
        f.write("# Grain 1, Grain 2, Misorientation angle, Misorientation axis, Type of Grain Boundary \n")
        np.savetxt(f, type_of_grain_boundaries_sigma_7, newline='\n', delimiter=',', fmt="%.4f")

    log.info('two_seed_testcase passed !!')