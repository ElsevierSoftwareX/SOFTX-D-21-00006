# -*- coding: utf-8 -*-
"""
test.py

Module to test structural as well as textural characteristics.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de 
created: 22 January 2020
Initial separated test cases created on 16 November 2019
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

from src.__version__ import __version__ as version

from src.create_tessellations import create_tessellations as create_tessellations

from src.create_obj_files import create_obj_file_all_grains as create_obj_file_all_grains
from src.create_obj_files import create_obj_file_individual_grains as create_obj_file_individual_grains
from src.create_vtk_files import create_vtk_file_all_grains as create_vtk_file_all_grains
from src.create_vtk_files import create_vtk_file_individual_grains as create_vtk_file_individual_grains

from src.random_generator import random_generator as random_generator
from src.structural_characteristics import grain_size_distribution as grain_size_distribution
from src.structural_characteristics import number_of_neighbors as number_of_neighbors
from src.structural_characteristics import grain_boundary_areas as grain_boundary_areas
from src.structural_characteristics import junction_length as junction_length
from src.structural_characteristics import junction_angle as junction_angle
from src.textural_characteristics import sharp_texture_quaternions as sharp_texture_quaternions
from src.textural_characteristics import random_quaternions_generator as random_quaternions_generator
from src.textural_characteristics import disorientation_angles as disorientation_angles
from src.textural_characteristics import type_of_grain_boundary as type_of_grain_boundary
from src.textural_characteristics import schmid_factor as schmid_factor
from src.textural_characteristics import available_required_texture as available_required_texture

from src.seed_spacing_files.bcc_lattice_3D import bcc_lattice_3D as bcc_lattice_3D
from src.seed_spacing_files.cubic_lattice_2D import cubic_lattice_2D as cubic_lattice_2D
from src.seed_spacing_files.cubic_lattice_3D import cubic_lattice_3D as cubic_lattice_3D
from src.seed_spacing_files.fcc_lattice_2D import fcc_lattice_2D as fcc_lattice_2D
from src.seed_spacing_files.fcc_lattice_3D import fcc_lattice_3D as fcc_lattice_3D
from src.seed_spacing_files.hcp_lattice_2D import hcp_lattice_2D as hcp_lattice_2D
from src.seed_spacing_files.hcp_lattice_3D import hcp_lattice_3D as hcp_lattice_3D

from test_cases.cubic_2d.cubic_2d_testcase import cubic_2d_testcase as cubic_2d_testcase
from test_cases.cubic_3d.cubic_3d_testcase import cubic_3d_testcase as cubic_3d_testcase
from test_cases.bcc_3d.bcc_3d_testcase import bcc_3d_testcase as bcc_3d_testcase
from test_cases.fcc_2d.fcc_2d_testcase import fcc_2d_testcase as fcc_2d_testcase
from test_cases.fcc_3d.fcc_3d_testcase import fcc_3d_testcase as fcc_3d_testcase
from test_cases.hcp_2d.hcp_2d_testcase import hcp_2d_testcase as hcp_2d_testcase
from test_cases.hcp_3d.hcp_3d_testcase import hcp_3d_testcase as hcp_3d_testcase
from test_cases.random_3d.random_3d_testcase import random_3d_testcase as random_3d_testcase
from test_cases.one_seed.one_seed_testcase import one_seed_testcase as one_seed_testcase
from test_cases.two_seed.two_seed_testcase import two_seed_testcase as two_seed_testcase
from test_cases.textural.textural_testcase import textural_testcase as textural_testcase

from src.execute_func import execute_func as execute_func

## Defining Pytest fixture for providing input arguments for main test function
@pytest.fixture
def name():
    return 'all'

@pytest.fixture
def f():
    return False

@pytest.fixture
def log_level():
    return 'INFO'

def test_func(name, f, log_level):
    """
    Function to iteratively execute required test case/s.

    Processing
    ---------- 
    The function gets the required data from the respective test case 
    function such as unique seed coordinates, etc and then computes and stores 
    generated data. The function also checks various assert statements based on 
    the name of the test case. 

    Parameters
    ----------
    
    name: string
        Name of the test case to be executed. Following are the available test cases:
            1. 'all'
            2. 'cubic_2d'
            3. 'cubic_3d'
            4. 'fcc_2d'
            5. 'fcc_3d'
            6. 'bcc_3d'
            7. 'hcp_2d'
            8. 'hcp_3d'
            9. 'random_3d'
            10. 'one_seed'
            11. 'two_seed'
            12. 'textural'

    f: boolean
        Flag to indicate that opaque surface is to be used instead of transparent.

    log_level: string
        Logger level to be used.

    Output
    -------
    Function stores appropriate outputs based on individual test cases.

    Returns
    -------
    Function returns nothing.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)    
    
    if log_level == 'DEBUG':
        ## starting yappi profiler
        yappi.start()
 
    
    ## Current time and date
    now = datetime.now()
    now = now.strftime("%Y%m%d_%H%M%S")
    
    store_folder = "visualization_output_tests"    
    face_flag = f                                                               ## Flag to indicate if closed surface is to be used or line
    skewed_boundary_flag = False
    number_of_bins = 10

    ## Dictionary with all test cases
    test_names = {'cubic_2d': ['cubic_2d'],
                    'cubic_3d': ['cubic_3d'],
                    'bcc_3d': ['bcc_3d'],
                    'fcc_2d': ['fcc_2d'],
                    'fcc_3d': ['fcc_3d'],
                    'hcp_2d': ['hcp_2d'],
                    'hcp_3d': ['hcp_3d'],
                    'random_3d': ['random_3D'],
                    'one_seed': ['one_seed'],
                    'two_seed': ['two_seed'],
                    'textural': ['textural'],
                    'all': ['textural', 'cubic_2d', 'cubic_3d', 'bcc_3d', 'fcc_2d', 'fcc_3d', 'hcp_2d', 'hcp_3d', 'random_3D', 'one_seed', 'two_seed']
                    }
    
    ## Identifies the tests to be executed based on input from user
    test_to_be_executed = test_names[name.lower()]
    
    ## Looping through all the tests
    for case_name in test_to_be_executed:
        log.info("Currently test case "  + case_name + " is running.")
        
        ## Common Parameters
                
        if case_name in ('cubic_2d', 'cubic_3d', 'bcc_3d', 'fcc_2d', 'fcc_3d', 'hcp_2d', 'hcp_3d'):
            size_of_simulation_box = 10.0
            spacing_lengths = [2.5, 5]                                       #1, 2, 2.5, 5]
        else:
            size_of_simulation_box = 10.0
            spacing_lengths = [1]

        required_texture = np.array([1, 1, 1])
        rand_quat_flag = True
        mesh_flag = 'TET'
        global_mesh_size = 0.5
        
        for spacing_length in spacing_lengths:
            ## Extracting required data based on name of the test case        
            if str(case_name).lower() == 'cubic_2D'.lower():
                log.info("Current spacing length = " + str(spacing_length))
                length_z = 1
                dimension = 2
                limit = np.array([size_of_simulation_box, size_of_simulation_box, length_z])
                material = "Cubic_2D_spacing_len_" + str(spacing_length)
                seed_array_unique = cubic_lattice_2D(limit, spacing_length, log_level)
                orientation_data = None
            
            elif str(case_name).lower() == 'cubic_3D'.lower():
                log.info("Current spacing length = " + str(spacing_length))
                limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
                dimension = 3
                material = "Cubic_3D_spacing_len_" + str(spacing_length)
                seed_array_unique = cubic_lattice_3D(limit, spacing_length, log_level)
                orientation_data = None

            elif str(case_name).lower() == 'bcc_3D'.lower():
                log.info("Current spacing length = " + str(spacing_length))
                limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
                dimension = 3
                material = "BCC_3D_spacing_len_" + str(spacing_length)
                seed_array_unique = bcc_lattice_3D(limit, spacing_length, log_level)
                orientation_data = None

            elif str(case_name).lower() == 'fcc_2D'.lower():
                log.info("Current spacing length = " + str(spacing_length))
                length_z = 1
                limit = np.array([size_of_simulation_box, size_of_simulation_box, length_z])
                dimension = 2
                material = "FCC_2D_spacing_len_" + str(spacing_length)
                seed_array_unique = fcc_lattice_2D(limit, spacing_length, log_level)
                orientation_data = None
            
            elif str(case_name).lower() == 'fcc_3D'.lower():
                log.info("Current spacing length = " + str(spacing_length))
                limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
                dimension = 3
                material = "FCC_3D_spacing_len_" + str(spacing_length)
                seed_array_unique = fcc_lattice_3D(limit, spacing_length, log_level)
                orientation_data = None

            elif str(case_name).lower() == 'hcp_2D'.lower():
                log.info("Current spacing length = " + str(spacing_length))
                length_z = 1
                limit = np.array([size_of_simulation_box, size_of_simulation_box, length_z])
                ## modifying limits
                if not np.isclose((limit[0] % spacing_length), 0):
                    limit[0] = limit[0] + (spacing_length - (limit[0]%spacing_length))
                    assert np.isclose(limit[0]%spacing_length, 0)
                if not np.isclose(limit[1] % (np.sqrt(3)*spacing_length), 0):
                    hcp_y_spacing = np.sqrt(3)*spacing_length
                    limit[1] = limit[1] + (hcp_y_spacing - (limit[1] % hcp_y_spacing))
                    assert np.isclose(limit[1]%hcp_y_spacing, 0)
                dimension = 2
                material = "HCP_2D_spacing_len_" + str(spacing_length)
                seed_array_unique = hcp_lattice_2D(limit, spacing_length, log_level)
                orientation_data = None            

            elif str(case_name).lower() == 'hcp_3D'.lower():
                log.info("Current spacing length = " + str(spacing_length))
                limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
                ## modifying limits
                if not np.isclose((limit[0] % spacing_length), 0):
                    limit[0] = limit[0] + (spacing_length - (limit[0]%spacing_length))
                    assert np.isclose(limit[0]%spacing_length, 0)
                if not np.isclose(limit[1] % (np.sqrt(3)*spacing_length), 0):
                    hcp_y_spacing = np.sqrt(3)*spacing_length
                    limit[1] = limit[1] + (hcp_y_spacing - (limit[1] % hcp_y_spacing))
                    assert np.isclose(limit[1]%hcp_y_spacing, 0)
                if not np.isclose((limit[2] % spacing_length), 0):
                    limit[2] = limit[2] + (spacing_length - (limit[2]%spacing_length))
                    assert np.isclose(limit[2]%spacing_length, 0)
                dimension = 3
                material = "HCP_3D_spacing_len_" + str(spacing_length)
                seed_array_unique = hcp_lattice_3D(limit, spacing_length, log_level)
                orientation_data = None 

            elif str(case_name).lower() == 'random_3D'.lower():
                number_of_seeds = 100
                limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
                dimension = 3
                material = "random_3D"
                
                ##Generating random seeds
                seed_coordinates_list = random_generator(number_of_seeds, dimension, limit, log_level)
                seed_array = np.array(seed_coordinates_list)
                new_array = [tuple(row) for row in seed_array]
                seed_array_unique = np.unique(new_array, axis = 0)

                orientation_data = None

            elif str(case_name).lower() == 'one_seed'.lower():
                length_z = 1
                dimension = 2
                material = "one_seed_2D"
                cell_number = 0
                seed_array_unique = np.array([[size_of_simulation_box/2, size_of_simulation_box/2, length_z]])                           # defining the seed array
                limit = np.array([size_of_simulation_box, size_of_simulation_box, length_z])                                      # defining the limits
                orientation_data = None

            elif str(case_name).lower() == 'two_seed'.lower():
                length_z = 1
                dimension = 2
                material = "two_seed_2D"
                seed_array_unique = np.array([[size_of_simulation_box/4, size_of_simulation_box/2, length_z], [0.75 * size_of_simulation_box, size_of_simulation_box/2, length_z]])
                limit = np.array([size_of_simulation_box, size_of_simulation_box, length_z])
                orientation_data = None
            
            elif str(case_name).lower() == 'textural'.lower():
                log.info("Current spacing length = " + str(spacing_length))
                limit = np.array([size_of_simulation_box, size_of_simulation_box, size_of_simulation_box])
                seed_array = np.zeros([(int(limit[0]) + 1) * (int(limit[1]) + 1) * (int(limit[2]) + 1), 3])
                dimension = 2
                material = "Textural_Test"
                seed_array_unique = cubic_lattice_2D(limit, spacing_length, log_level)
                orientation_data = None

            else:
                log.critical('Please enter appropriate test keyword mentioned in the documentation.')
                exit()
            
            ## Creating tessellations
            tessellation = create_tessellations(seed_array_unique, limit, log_level)

            ## Calling another function which executes the characteristics functions and store the computed data
            stress_direction = np.array([0, 1, 0])
            grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries = \
                execute_func(size_of_simulation_box, dimension, limit, material, orientation_data, required_texture, rand_quat_flag, \
                seed_array_unique, stress_direction, store_folder, face_flag, now, number_of_bins, skewed_boundary_flag, mesh_flag, global_mesh_size, log_level)

            ## Executing various assert statements based on the name of the test case
            if str(case_name).lower() == 'cubic_2D'.lower():
                cubic_2d_testcase(tessellation, dimension, size_of_simulation_box, spacing_length, length_z, grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level)        
            
            elif str(case_name).lower() == 'cubic_3D'.lower():
                cubic_3d_testcase(tessellation, dimension, size_of_simulation_box, spacing_length, grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level)         
            
            elif str(case_name).lower() == 'bcc_3D'.lower():
                bcc_3d_testcase(tessellation, dimension, size_of_simulation_box, spacing_length, grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level)
            
            elif str(case_name).lower() == 'fcc_2D'.lower():
                fcc_2d_testcase(tessellation, dimension, size_of_simulation_box, spacing_length, length_z, grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level)
            
            elif str(case_name).lower() == 'fcc_3D'.lower():
                fcc_3d_testcase(tessellation, dimension, size_of_simulation_box, spacing_length, grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level)
            
            elif str(case_name).lower() == 'hcp_2D'.lower():
                hcp_2d_testcase(tessellation, dimension, size_of_simulation_box, spacing_length, length_z, grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level)

            elif str(case_name).lower() == 'random_3D'.lower():
                random_3d_testcase(store_folder, version, now, material, tessellation, dimension, size_of_simulation_box, spacing_length, seed_array_unique, required_texture, number_of_seeds, grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level)
            
            elif str(case_name).lower() == 'one_seed'.lower():
                one_seed_testcase(tessellation, cell_number, size_of_simulation_box, length_z, log_level)
            
            elif str(case_name).lower() == 'two_seed'.lower():
                two_seed_testcase(store_folder, version, now, material, tessellation, dimension, size_of_simulation_box, spacing_length, length_z, required_texture, rand_quat_flag, grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level)
            
            elif str(case_name).lower() == 'textural'.lower():
                textural_testcase(store_folder, tessellation, dimension, size_of_simulation_box, spacing_length, seed_array_unique, required_texture, now, material, grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, \
                junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries, log_level)

    log.info('Successfully executed testcase/s')

    if log_level == 'DEBUG':    
        ## Saving profiler output
        stats = yappi.get_func_stats()

        output_file_path = Path(store_folder, material, now, "profiler_output", "callgrind.yappi_profiler_output.prof")
        output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

        stats.save(str(output_file_path), type='callgrind')
        
        log.info('Successfully saved profiler data file')
        yappi.get_thread_stats().print_all()

    ## Moving log file to appropriate directory in visualization_files
    output_file_path = Path(store_folder, material, now, "log_directory", "log_data.log")
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    os.system('mv log_data.lo* ' + str(output_file_path))

