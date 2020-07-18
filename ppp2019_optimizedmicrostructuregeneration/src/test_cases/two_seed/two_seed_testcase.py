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

    log = set_logger(name_str, 'log_data.log', log_level)
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
    disorientation_angle, orientation_data = disorientation_angles(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    type_of_grain_boundaries_sigma_15, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    assert (np.all(type_of_grain_boundaries_sigma_15[:, 6] == 15))
    assert (np.all(np.isclose(disorientation_angle[:, 2], 48.2, atol=1e-1)))

    ## Saving to array_verification_file.txt
    output_file_path = Path("visualization_files", store_folder, now, material, "Text_output", "array_verification_file.txt")
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
    disorientation_angle, orientation_data = disorientation_angles(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    type_of_grain_boundaries_sigma_5, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
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
    disorientation_angle, orientation_data = disorientation_angles(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    type_of_grain_boundaries_sigma_3, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
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
    disorientation_angle, orientation_data = disorientation_angles(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    type_of_grain_boundaries_sigma_7, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
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