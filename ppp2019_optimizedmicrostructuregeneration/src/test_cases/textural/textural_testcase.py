# -*- coding: utf-8 -*-
"""
textural_testcase.py

Module to test functionalities related to textural characteristics.

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

from ppp2019_optimizedmicrostructuregeneration.src.seed_spacing_files.cubic_lattice_2D import cubic_lattice_2D as cubic_lattice_2D

from ppp2019_optimizedmicrostructuregeneration.src.set_logger import set_logger as set_logger
name_str = __name__

from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import sharp_texture_quaternions as sharp_texture_quaternions
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import random_quaternions_generator as random_quaternions_generator
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import disorientation_angles as disorientation_angles
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import type_of_grain_boundary as type_of_grain_boundary
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import schmid_factor as schmid_factor
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import available_required_texture as available_required_texture

def textural_testcase(store_folder, tessellation, dimension, \
    size_of_simulation_box, spacing_length, seed_array_unique, \
    required_texture, now, material, grain_size_distributions, \
    number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
    junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, \
    disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level):
    
    """
    Testing the textural characteristics function with seeds at regular grid in 3D
    Input: Tessellations data and orientations data
    Output: Known disorientation angles and schmid factors
    """

    log = set_logger(name_str, 'log_data.log', log_level)
    #size_of_simulation_box = 10.0
    #spacing_length = 1
    #limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
    #seed_array = np.zeros([(int(limit[0]) + 1) * (int(limit[1]) + 1) * (int(limit[2]) + 1), 3])
    #dimension = 2
    #material = "Textural_Test"
    #required_texture = np.array([1, 1, 1])
    rand_quat_flag = True

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
    disorientation_angle, orientation_data = disorientation_angles(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    assert np.all([angle[2] <= 62.8 for angle in disorientation_angle])

    ## Saving to array_verification_file.txt
    output_file_path = Path("visualization_files", store_folder, now, material, "Text_output", "array_verification_file.txt")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    ## Testing with random orientation but same random orientation for each grain
    orientation_quaternion = sharp_texture_quaternions(1, required_texture, log_level)                       # Generating a random quaternion
    orientation_data = np.broadcast_to(orientation_quaternion, (seed_array_unique.shape[0], 4))   # Assigning same orientation to each grain
    disorientation_angle, orientation_data = disorientation_angles(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
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
    disorientation_angle, orientation_data = disorientation_angles(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    assert np.all([angle[2] <= 62.8 for angle in disorientation_angle])
    assert np.all([angle[2] == 0 for angle in disorientation_angle])
    
    ####################################################################################################################
    #Testing for known inputs of stress directions and known schmid factors
    ####################################################################################################################
    ## Schmid Factors can have a maximum value of 0.5
    ## Testing for various stress directions

    ## With random orientations
    orientation_data = None
    stress_direction = np.array([1, 0, 0])
    schmid_factors, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, stress_direction, orientation_data, tessellation, log_level)
    assert np.all([factors[1] <= 0.5 for factors in schmid_factors])

    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for schmid factor for stress direction (1, 0, 0) \n")
        f.write("\n# Grain number, Maximum schmid factor, Respective slip plane and slip direction (Slip System) \n")
        np.savetxt(f, schmid_factors, newline='\n', delimiter=',', fmt="%.4f")

    ## With same orientations but with different stress directions
    orientation_quaternion = np.array([1, 0, 0, 0])
    orientation_data = np.broadcast_to(orientation_quaternion, (seed_array_unique.shape[0], 4))   # Assigning same orientation to each grain

    stress_direction = np.array([1, -1, 0])
    schmid_factors, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, stress_direction, orientation_data, tessellation, log_level)
    assert np.all([factors[1] <= 0.5 for factors in schmid_factors])
    assert np.all([np.isclose(factor[1], 0.408, atol=1e-3) for factor in schmid_factors])

    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for schmid factor for stress direction (1, -1, 0) \n")
        f.write("\n# Grain number, Maximum schmid factor, Respective slip plane and slip direction (Slip System) \n")
        np.savetxt(f, schmid_factors, newline='\n', delimiter=',', fmt="%.4f")


    stress_direction = np.array([1, 0, 0])
    schmid_factors, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, stress_direction, orientation_data, tessellation, log_level)
    assert np.all([factors[1] <= 0.5 for factors in schmid_factors])
    assert np.all([np.isclose(factor[1], 0.408, atol=1e-3) for factor in schmid_factors])

    with open(str(output_file_path), 'a+') as f:
        f.write("\n \n Verification matrix for schmid factor for stress direction (1, 0, 0) \n")
        f.write("\n# Grain number, Maximum schmid factor, Respective slip plane and slip direction (Slip System) \n")
        np.savetxt(f, schmid_factors, newline='\n', delimiter=',', fmt="%.4f")


    stress_direction = np.array([1, 1, 0])
    schmid_factors, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, stress_direction, orientation_data, tessellation, log_level)
    assert np.all([factors[1] <= 0.5 for factors in schmid_factors])
    assert np.all([np.isclose(factor[1], 0.408, atol=1e-3) for factor in schmid_factors])

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
    type_of_grain_boundaries, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)

    from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import available_required_texture
    from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import type_of_csl_data
    available_required_texture_info = copy.deepcopy(available_required_texture)
    type_of_csl_dict = copy.deepcopy(type_of_csl_data)
    for row in available_required_texture_info:
        required_texture = np.around((1/np.linalg.norm(row))*row, decimals=4)
        key = "".join(str(list(required_texture)))
        assert (key in type_of_csl_dict)

    ####################################################################################################################
    ## Saving Textural Characteristics
    output_file_path = Path("visualization_files", store_folder, now, material, "Text_output", "textural_characteristics.txt")
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
        f.write("\n# Grain 1, Grain 2, Disorientation angle, Disorientation axis \n")
        np.savetxt(f, disorientation_angle, delimiter=',', fmt="%.4f")
        f.write("\n \n# Type of Grain Boundaries \n")
        f.write("\n# Grain 1, Grain 2, Misorientation angle, Misorientation axis, Type of Grain Boundary \n")
        np.savetxt(f, type_of_grain_boundaries, delimiter=',', fmt='%.4f %.4f %.4f %.4f %.4f %.4f %.4f')
        f.write("\n# Schmid Factors")
        f.write("\n# Grain number, Maximum schmid factor, Respective slip plane and slip direction (Slip System) \n")
        np.savetxt(f, schmid_factors, delimiter=',', fmt="%.4f")

    log.info('textural_testcase passed !!')