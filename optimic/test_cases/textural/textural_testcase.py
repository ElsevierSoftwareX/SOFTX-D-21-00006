# -*- coding: utf-8 -*-
"""
textural_testcase.py

Module to test functionalities related to textural characteristics.

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

from optimic.src.main_import_statements import *

from optimic.src.seed_spacing_files.cubic_lattice_2D import cubic_lattice_2D 

from optimic.src.set_logger import set_logger 
name_str = __name__

from optimic.src.textural_characteristics import sharp_texture_quaternions 
from optimic.src.textural_characteristics import random_quaternions_generator 
from optimic.src.textural_characteristics import disorientation_angles 
from optimic.src.textural_characteristics import type_of_grain_boundary 
from optimic.src.textural_characteristics import schmid_factor 
from optimic.src.textural_characteristics import available_required_texture
from optimic.src.textural_characteristics import slip_system_generator 

def textural_testcase(store_folder, tessellation, dimension, \
    size_of_simulation_box, spacing_length, seed_array_unique, \
    required_texture, now, material, grain_size_distributions, \
    number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
    junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, \
    disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level):
    """
    Execute all test conditions related to textural_testcase.

    Parameters
    ----------
    store_folder: string
        Name of directory where output files are to be stored.

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

    seed_array_unique: array of shape (number of grains, 3)
        Unique seed coordinates

    required_texture: array of length 3
        Specific texture to be used for all grains

    now: string
        Current time and date.

    material: string
        Material name of which microstructure is being tested.

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
    Following files in appropriate directory within specified storage_folder:
        1. array_verification_file.txt
        2. textural_characteristics.txt
    """

    log = set_logger(name_str, 'log_data.log', log_level)
    #size_of_simulation_box = 10.0
    #spacing_length = 1
    limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
    #seed_array = np.zeros([(int(limit[0]) + 1) * (int(limit[1]) + 1) * (int(limit[2]) + 1), 3])
    #dimension = 2
    #material = "Textural_Test"
    #required_texture = np.array([1, 1, 1])
    rand_quat_flag = True
    skewed_boundary_flag = False

    #seed_array_unique = cubic_lattice_2D(limit, spacing_length)
    
    ## Creating tessellations
    #tessellation = create_tessellations(seed_array_unique, limit)

    ## Textural Characteristics
    ####################################################################################################################
    #Testing random quaternions function
    ####################################################################################################################
    number_of_quaternions = 1
    quaternion_array = sharp_texture_quaternions(number_of_quaternions, required_texture, log_level)
    assert np.isclose(np.linalg.norm(quaternion_array), 1)

    ####################################################################################################################
    ## Testing Disorientation_angles function
    ## Major testing is done in 'two_seed' case
    ####################################################################################################################
    ## Due to cubic symmetry all disorientation angles should be less than 62.8 degrees
    
    orientation_data = None
    disorientation_angle, orientation_data = disorientation_angles(dimension, limit, skewed_boundary_flag, required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    assert np.all([angle[2] <= 62.8 for angle in disorientation_angle])

    ## Saving to array_verification_file.txt
    output_file_path = Path(store_folder, material, now, "Text_output", "array_verification_file.txt")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    ## Testing with random orientation but same random orientation for each grain
    orientation_quaternion = sharp_texture_quaternions(1, required_texture, log_level)                       # Generating a random quaternion
    orientation_data = np.broadcast_to(orientation_quaternion, (seed_array_unique.shape[0], 4))   # Assigning same orientation to each grain
    disorientation_angle, orientation_data = disorientation_angles(dimension, limit, skewed_boundary_flag, required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    assert np.all([angle[2] <= 62.8 for angle in disorientation_angle])
    assert np.all([angle[2] == 0 for angle in disorientation_angle])

    with open(str(output_file_path), 'a+') as f:
        f.truncate(0)
        f.write("Verification matrix for disorientation matrix for same random orientation for all grains \n")
        f.write("# Grain 1, Grain 2, Disorientation angle, Disorientation axis \n")
        np.savetxt(f, disorientation_angle, newline='\n', delimiter=',', fmt="%.4f")

    ## Testing for all grains having same crystal orientation
    orientation_quaternion = np.array([1, 0, 0, 0])
    orientation_data = np.broadcast_to(orientation_quaternion, (seed_array_unique.shape[0], 4))   # Assigning same orientation to each grain
    disorientation_angle, orientation_data = disorientation_angles(dimension, limit, skewed_boundary_flag, required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    assert np.all([angle[2] <= 62.8 for angle in disorientation_angle])
    assert np.all([angle[2] == 0 for angle in disorientation_angle])

    ####################################################################################################################
    #Testing slip system generator function of number of slip systems generated
    #and orthogonality check
    ####################################################################################################################
    
    ## FCC has 12 slip systems [111]<110>
    slip_system_family = np.array([1,1,1,1,1,0])
    crystal_symmetry_type = 'CUBIC'
    slip_systems = slip_system_generator(slip_system_family, crystal_symmetry_type)
    assert len(slip_systems)==12                                                ## FCC has 12 slip systems
    all_normals = slip_systems[:, :3]
    all_dir = slip_systems[:,3:]
    all_dot_products = np.einsum('ij,ij->i', all_normals, all_dir)
    assert np.all(np.isclose(all_dot_products,0.0))                                          ## for orthogonality, dot product must be zero

    ## BCC has 12 slip systems [110]<-111>
    slip_system_family = np.array([1,1,0,-1,1,1])
    crystal_symmetry_type = 'CUBIC'
    slip_systems = slip_system_generator(slip_system_family, crystal_symmetry_type)
    assert len(slip_systems)==12                                                ## BCC has 12 slip systems
    all_normals = slip_systems[:, :3]
    all_dir = slip_systems[:,3:]
    all_dot_products = np.einsum('ij,ij->i', all_normals, all_dir)
    assert np.all(np.isclose(all_dot_products,0.0))                                          ## for orthogonality, dot product must be zero

    ## BCC has 12 slip systems [211]<-111>
    slip_system_family = np.array([2,1,1,-1,1,1])
    crystal_symmetry_type = 'CUBIC'
    slip_systems = slip_system_generator(slip_system_family, crystal_symmetry_type)
    assert len(slip_systems)==12                                                ## BCC has 12 slip systems
    all_normals = slip_systems[:, :3]
    all_dir = slip_systems[:,3:]
    all_dot_products = np.einsum('ij,ij->i', all_normals, all_dir)
    assert np.all(np.isclose(all_dot_products,0.0))                                          ## for orthogonality, dot product must be zero

    ## BCC has 24 slip systems [321]<-111>
    slip_system_family = np.array([3,2,1,-1,1,1])
    crystal_symmetry_type = 'CUBIC'
    slip_systems = slip_system_generator(slip_system_family, crystal_symmetry_type)
    assert len(slip_systems)==24                                                ## BCC has 24 slip systems
    all_normals = slip_systems[:, :3]
    all_dir = slip_systems[:,3:]
    all_dot_products = np.einsum('ij,ij->i', all_normals, all_dir)
    assert np.all(np.isclose(all_dot_products,0.0))                                          ## for orthogonality, dot product must be zero

    ## Body centered orthorhombic crystal
    slip_system_family = np.array([3,2,1,-1,1,1])
    crystal_symmetry_type = ['ORTHORHOMBIC', 'TETRAGONAL', 'HEXAGONAL']
    for sym_type in crystal_symmetry_type:
        slip_systems = slip_system_generator(slip_system_family, sym_type)
        #assert len(slip_systems)==12                                                ## FCC has 12 slip systems
        all_normals = slip_systems[:, :3]
        all_dir = slip_systems[:,3:]
        all_dot_products = np.einsum('ij,ij->i', all_normals, all_dir)
        assert np.all(np.isclose(all_dot_products,0.0))                                          ## for orthogonality, dot product must be zero


    ####################################################################################################################
    #Testing for known inputs of stress directions and known schmid factors
    ####################################################################################################################
    ## Schmid Factors can have a maximum value of 0.5
    ## Testing for various stress directions

    slip_system_family = np.array([1,1,1,1,1,0])
    crystal_symmetry_type = 'CUBIC'

    ## With random orientations
    orientation_data = None
    stress_direction = np.array([1, 0, 0])
    schmid_factors, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, limit, stress_direction, slip_system_family, crystal_symmetry_type, orientation_data, tessellation, log_level)
    assert np.all([factors[2] <= 0.5 for factors in schmid_factors])

    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for schmid factor for stress direction (1, 0, 0) \n")
        f.write("\n# Grain number, Maximum schmid factor, Respective slip plane and slip direction (Slip System) \n")
        np.savetxt(f, schmid_factors, newline='\n', delimiter=',', fmt="%.4f")

    ## With same orientations but with different stress directions
    orientation_quaternion = np.array([1, 0, 0, 0])
    orientation_data = np.broadcast_to(orientation_quaternion, (seed_array_unique.shape[0], 4))   # Assigning same orientation to each grain

    stress_direction = np.array([1, -1, 0])
    schmid_factors, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, limit, stress_direction, slip_system_family, crystal_symmetry_type, orientation_data, tessellation, log_level)
    assert np.all([factors[2] <= 0.5 for factors in schmid_factors])
    assert np.all([np.isclose(factor[2], 0.408, atol=1e-3) for factor in schmid_factors])

    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for schmid factor for stress direction (1, -1, 0) \n")
        f.write("\n# Grain number, Maximum schmid factor, Respective slip plane and slip direction (Slip System) \n")
        np.savetxt(f, schmid_factors, newline='\n', delimiter=',', fmt="%.4f")


    stress_direction = np.array([1, 0, 0])
    schmid_factors, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, limit, stress_direction, slip_system_family, crystal_symmetry_type, orientation_data, tessellation, log_level)
    assert np.all([factors[2] <= 0.5 for factors in schmid_factors])
    assert np.all([np.isclose(factor[2], 0.408, atol=1e-3) for factor in schmid_factors])

    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for schmid factor for stress direction (1, 0, 0) \n")
        f.write("\n# Grain number, Maximum schmid factor, Respective slip plane and slip direction (Slip System) \n")
        np.savetxt(f, schmid_factors, newline='\n', delimiter=',', fmt="%.4f")


    stress_direction = np.array([1, 1, 0])
    schmid_factors, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, limit, stress_direction, slip_system_family, crystal_symmetry_type, orientation_data, tessellation, log_level)
    assert np.all([factors[2] <= 0.5 for factors in schmid_factors])
    assert np.all([np.isclose(factor[2], 0.408, atol=1e-3) for factor in schmid_factors])

    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for schmid factor for stress direction (1, 1, 0) \n")
        f.write("\n# Grain number, Maximum schmid factor, Respective slip plane and slip direction (Slip System) \n")
        np.savetxt(f, schmid_factors, newline='\n', delimiter=',', fmt="%.4f")

    ####################################################################################################################
    ## Testing Type of Grain Boundaries
    ## Major testing is done in 'two_seed' case
    ####################################################################################################################

    required_texture = np.array([2, 1, 0])
    rand_quat_flag = True
    type_of_grain_boundaries, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, dimension, limit, skewed_boundary_flag, log_level)

    from optimic.src.textural_characteristics import available_required_texture
    from optimic.src.textural_characteristics import type_of_csl_data
    available_required_texture_info = copy.deepcopy(available_required_texture)
    type_of_csl_dict = copy.deepcopy(type_of_csl_data)
    for row in available_required_texture_info:
        required_texture = np.around((1/np.linalg.norm(row))*row, decimals=4)
        key = "".join(str(list(required_texture)))
        assert (key in type_of_csl_dict)

    ####################################################################################################################
    ## Saving Textural Characteristics
    output_file_path = Path(store_folder, material, now, "Text_output", "textural_characteristics.txt")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    with open(str(output_file_path), 'a+') as f:
        f.truncate(0)
        f.write(now)
        
        f.write("\n\nSize of Simulation Box: " + str(size_of_simulation_box))
        #f.write("\nNumber of Seeds: " + str(number_of_seeds))
        f.write("\nDimension: " + str(dimension))
        f.write("\nMaterial Name: " + material)
        f.write("\n \n# Disorientation Angles \n")
        f.write("\n# Grain 1, Grain 2, Disorientation angle, Disorientation axis, GB area \n")
        np.savetxt(f, disorientation_angle, delimiter=',', fmt="%.4f")
        f.write("\n \n# Type of Grain Boundaries \n")
        f.write("\n# Grain 1, Grain 2, Misorientation angle, Misorientation axis, Type of Grain Boundary, Grain Boundary area \n")
        np.savetxt(f, type_of_grain_boundaries, delimiter=',', fmt='%.14f,%.14f,%.14f,%.14f,%.14f,%.14f,%s,%.14f')
        f.write("\n# Schmid Factors")
        f.write("\n# Grain number, Grain sizes, Maximum schmid factor, Respective slip plane and slip direction (Slip System), 2nd Maximum schmid factor, Respective Slip System, 3rd Maximum schmid factor, Respective Slip System \n")
        np.savetxt(f, schmid_factors, delimiter=',', fmt="%.4f")

    log.info('textural_testcase passed !!')