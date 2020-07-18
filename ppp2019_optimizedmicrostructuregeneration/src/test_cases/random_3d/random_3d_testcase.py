# -*- coding: utf-8 -*-
"""
random_3d_testcase.py

Module to test some features of program for RANDOM 3D type spacing.

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

from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import sharp_texture_quaternions as sharp_texture_quaternions
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import random_quaternions_generator as random_quaternions_generator
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import disorientation_angles as disorientation_angles
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import type_of_grain_boundary as type_of_grain_boundary
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import schmid_factor as schmid_factor
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import available_required_texture as available_required_texture

def random_3d_testcase(store_folder, version, now, material, tessellation, \
    dimension, size_of_simulation_box, spacing_length, seed_array_unique, \
    required_texture, number_of_seeds, grain_size_distributions, \
    number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
    junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, \
    disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level):

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
    output_file_path = Path("visualization_files", store_folder, now, material, "Text_output", "array_verification_file.txt")
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
