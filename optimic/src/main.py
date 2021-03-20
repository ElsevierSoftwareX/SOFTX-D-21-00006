# -*- coding: utf-8 -*-
"""
main.py

Module to execute main program.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash, Stefan Sandfeld
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
created: 22 December 2019
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

## Below commented lines are command line input for specified output

## INPUT FROM BASH COMMAND IS TO BE PROVIDED AS for running main for no optimization: $ ./optimic main -s 18 12 15 -sl 1 -d 2 -n 122 -m steel -st 0 1 0 -ss cubic_2d -rs 1 -noopti
## INPUT FROM BASH COMMAND IS TO BE PROVIDED AS for running main for optimization: $ ./optimic main -s 180 120 15 -sl 1 -d 2 -n 122 -t user_grain_size_distribution.txt -c 0 -m steel -st 0 1 0 -ss random_3d -rs 10 -om COBYLA -mf 2000
## INPUT FROM BASH COMMAND IS TO BE PROVIDED AS for running main for optimization using USER defined COST Function: $ ./optimic main -s 180 120 15 -sl 1 -d 2 -n 122 -t user_grain_size_distribution.txt -c 3 -c 5 -m steel -st 0 1 0 -ss random_3d -rs 10 -om COBYLA -mf 2000 -ucf user_cost_function_2.py
## INPUT FROM BASH COMMAND IS TO BE PROVIDED AS for running tests: $ ./optimic test --name cubic_2d
## INPUT FROM BASH COMMAND IS TO BE PROVIDED AS for running all tests: $ ./optimic test --name all
## INPUT FROM BASH COMMAND IS TO BE PROVIDED AS for optimizing grain size, disorientation angles and distance between grains: $ ./optimic main -s 180 120 15 -sl 1 -d 2 -n 122 -t user_grain_size_distribution.txt -c 0 -c 6 -c 5 -m steel -st 0 1 0 -ss random_3d -rs 10 -om COBYLA -mf 2000
## For help: $ ./optimic --help
## For help: $ ./optimic test --help
## For help: $ ./optimic main --help

## Using faulthandler incase segmentation fault occurs

# import faulthandler
# faulthandler.enable()
from optimic.src.main_import_statements import *
from optimic.src.__version__ import __version__ as version

from optimic.src.create_tessellations import create_tessellations

"""
Import statements are written in this fashion due to philosophical reason. 
Importing functions in this way helps me to quickly identify the module 
containing a specific function. All modules in this package use this style for
absolute import of functions from another modules.
"""
from optimic.src.create_obj_files import create_obj_file_all_grains
from optimic.src.create_obj_files import create_obj_file_individual_grains
from optimic.src.create_vtk_files import create_vtk_file_all_grains 
from optimic.src.create_vtk_files import create_vtk_file_individual_grains 

from optimic.src.random_generator import random_generator 
from optimic.src.structural_characteristics import grain_size_distribution 
from optimic.src.structural_characteristics import number_of_neighbors 
from optimic.src.structural_characteristics import grain_boundary_areas 
from optimic.src.structural_characteristics import junction_length 
from optimic.src.structural_characteristics import junction_angle 
from optimic.src.structural_characteristics import distance_btw_grains 
from optimic.src.textural_characteristics import sharp_texture_quaternions 
from optimic.src.textural_characteristics import random_quaternions_generator 
from optimic.src.textural_characteristics import disorientation_angles 
from optimic.src.textural_characteristics import type_of_grain_boundary 
from optimic.src.textural_characteristics import schmid_factor 

from optimic.src.seed_spacing_files.bcc_lattice_3D import bcc_lattice_3D 
from optimic.src.seed_spacing_files.cubic_lattice_2D import cubic_lattice_2D 
from optimic.src.seed_spacing_files.cubic_lattice_3D import cubic_lattice_3D 
from optimic.src.seed_spacing_files.fcc_lattice_2D import fcc_lattice_2D 
from optimic.src.seed_spacing_files.fcc_lattice_3D import fcc_lattice_3D 
from optimic.src.seed_spacing_files.hcp_lattice_2D import hcp_lattice_2D 
from optimic.src.seed_spacing_files.hcp_lattice_3D import hcp_lattice_3D 

from optimic.src.test import test_func 
from optimic.src.execute_func import execute_func 

from optimic.src.check_libraries import check_libraries 

from optimic.src.set_logger import set_logger 

name_str = __name__

def cost_function_general(combined_user_data, combined_predicted_data, start_row_combined_data, data_dictionary, args):
    """
    Calculating cost value using Sum of Squared Error (SSE) between the 
    predicted data and the user data.

    Parameters
    ----------
    combined_user_data: array
        Distribution of user defined data for all characteristics.

    combined_predicted_data: array
        Distribution of predicted data for all characteristics.

    start_row_combined_data: list
        List comprising of row indexes of user and predicted data associated 
        with each characteristic to be optimized.

    data_dictionary: dictionary
        Dictionary consisting of complete data pertaining to characteristic to 
        be optimized. Keys used in the dictionary are integers used for 
        identifying appropriate characteristics as below:
        '0': Grain sizes
        '1': Number of neighbors
        '2': Grain boundary areas
        '3': Junction lengths
        '4': Junction angles in degrees
        '5a': Distance between grains as 1D array
        '5b': Distance between grains as matrix array
        '6': Disorientation angles
        '7': Type of grain boundaries
        '8': Schmid factors

    args: list
        Common arguments consisting of following:
            1. parameter_list
            2. dimension
            3. user_data
            4. start_row_of_parameter
            5. limit
            6. number_of_bins
            7. fig_animate
            8. line, ax_animate
            9. cost_function_names
            10. func_name_key
            11. required_texture
            12. rand_quat_flag
            13. stress_direction
            14. orientation_data
            15. skewed_boundary_flag
            16. tessellation (DICTIONARY OF ALL TESSELLATION DATA)
            17. limits_factor
            18. slip_system_family
            19. crystal_symmetry_type

    Returns
    -------
    Scalar SSE cost value.
    """
    
    ## Extracting USER and PREDICTED Y data as 1D array
    user_y = combined_user_data[:, 1]
    predicted_y = combined_predicted_data[:, 1]
    
    ############################################################################
    ## Below commented lines are kept to remind myself about why 'start_row_combined_data' was created 
    ############################################################################
    # ## Iterating to accumulate cost function value for each parameter
    # value = 0    
    # for i, j in zip(start_row_combined_data[:-1], start_row_combined_data[1:]):
    #     value += (np.sum(np.square(user_y[i:j] - predicted_y[i:j])))
    
    value = (np.sum(np.square(user_y - predicted_y)))

    return value

class optimize_class():
    """
    Using classes so that use of Global variables can be avoided
    """
    def __init__(self, store_folder, now, material, save_interval, log_level):
        self.initial_distribution = None
        self.final_user_distribution_data = None
        self.final_predicted_data = None
        self.final_start_row_combined_data = None
        self.smallest_cost_function_value = np.inf 
        self.iteration_count = [0]
        self.func_eval_count = [0]
        self.current_cost_function_value = []
        self.seeds_array_all_iterations = []
        self.log_level = log_level
        self.store_folder = store_folder
        self.now = now
        self.material = material
        self.save_interval = save_interval

    def callback_f(self, x):
        self.iteration_count.append(self.iteration_count[-1] + 1)
        print("Finished with iteration: " + str(self.iteration_count[-1]))

    def cost_function(self, seed_p, args):                                            # *args is a tuple containing parameter, dimension, user_required_distribution, limit
        """
        The function represents the objective function for the optimization 
        process. The function interpolates data for the provided user data and 
        the predicted data such that there are equal number of entries in both
        user and predicted data array. This helps in comparison of both data.
        The function then calls the cost function containing the evaluation 
        steps to determine MSE. The function then updates a set of common 
        variables and updates the live plot.

        Parameters
        ----------
        seed_p: 1D array
            Predicted seeds coordinates in 3D

        args: list
            Common arguments consisting of following:
                1. parameter_list
                2. dimension
                3. user_data
                4. start_row_of_parameter
                5. limit
                6. number_of_bins
                7. fig_animate
                8. line, ax_animate
                9. cost_function_names
                10. func_name_key
                11. required_texture
                12. rand_quat_flag
                13. stress_direction
                14. orientation_data
                15. skewed_boundary_flag
                16. tessellation_dummy_variable
                17. limits_factor
                18. slip_system_family
                19. crystal_symmetry_type

        Returns 
        -------
        The function returns scalar cost function value.
        """

        log = set_logger(name_str, 'log_data.log', self.log_level)
        ## Extracting the data from input arguments
        ## Unique seed coordinates are not isolated as this would lead to inappropriate
        ## results as the number of seeds would be less compared to requirement.
        seed_array_wo_factors = np.around(np.reshape(seed_p, (-1, 3)), decimals=14)

        # new_array = [tuple(row) for row in seed_array]
        # seed_array_unique = np.unique(new_array, axis = 0)

        #args_list = [parameter_list, dimension, user_data, start_row_of_parameter, limit, number_of_bins, fig_animate, ax_animate, cost_function_names, func_name_key, required_texture, rand_quat_flag, stress_direction, orientation_data, skewed_boundary_flag]
        parameter_list = args[0]
        dimension = args[1]
        user_data = args[2]
        start_row_of_parameter = args[3]
        limit = args[4]
        number_of_bins = args[5]
        fig_animate = args[6]
        line, ax_animate = args[7]
        cost_function_names = args[8]
        func_name_key = args[9]
        required_texture = args[10]
        rand_quat_flag = args[11]
        stress_direction = args[12]
        orientation_data = args[13]
        skewed_boundary_flag = args[14]
        limits_factor = args[16]
        slip_system_family = args[17]
        crystal_symmetry_type = args[18]

        log.debug('Successfully extracted all information from args in the cost function')

        ## Scaling seed coordinates with respective factors    
        seed_array = seed_array_wo_factors*limits_factor

        ## Adjusting seed array based on dimension
        if dimension == 2:
            seed_array[:, 2] = limit[2]

        ## Setting log_level to critical while evaluating charateristic functions
        original_log_level = copy.deepcopy(self.log_level)
        self.log_level = 'CRITICAL'

        ## Creating tessellations
        ## Using try and return to skip the iteration if any error occurs during 
        ## tessellations.
        try:
            tessellation = create_tessellations(seed_array, limit, self.log_level)
            args[15] = tessellation                                             ## Dummy variable is replaced with actual tessellation data
        except ValueError:
            log.exception('tessellation creation failed in cost function and infinity would be returned as cost function value')
            return 1e8

        ## Creating an empty array for storing the USER and PREDICTED data as an array
        combined_user_data = np.empty([0, 2])
        combined_predicted_data = np.empty([0, 2])
        ## Creating a list to indicate the start and end row numbers for each parameter
        start_row_combined_data = [0]

        ## Parsing data array into a dictionary

        data_dictionary = {}

        ## Iterating through each parameter
        for index, parameter in enumerate(parameter_list):
            ## Based on parameter finding the distribution using histogram
            if parameter is 'grain_size_distribution':
                grain_size_distributions = np.around(grain_size_distribution(dimension, tessellation, limit, self.log_level), decimals=14)
                hist, bins = np.histogram(grain_size_distributions[:, 1], bins= number_of_bins, density= True)
                data_dictionary['0'] = grain_size_distributions                     # Key as integer refers to the integer representing the characteristic feature

            elif parameter is 'number_of_neighbors':
                number_of_neighbor = number_of_neighbors(dimension, tessellation, limit, self.log_level)
                neighbors_array = np.array([v[2] for v in number_of_neighbor])
                hist, bins = np.histogram(neighbors_array, bins= number_of_bins, density= True)
                data_dictionary['1'] = number_of_neighbor                     # Key as integer refers to the integer representing the characteristic feature

            elif parameter is 'grain_boundary_areas':
                parent_function_name = inspect.stack()[1][3]
                grain_boundary_area_distribution, all_vertices_list = grain_boundary_areas(dimension, limit, tessellation, parent_function_name, skewed_boundary_flag, self.log_level)
                all_areas = [v[3] for v in grain_boundary_area_distribution]
                hist, bins = np.histogram(all_areas, bins= number_of_bins, density= True)
                data_dictionary['2'] = grain_boundary_area_distribution                     # Key as integer refers to the integer representing the characteristic feature

            elif parameter is 'junction_length':
                junction_lengths = junction_length(dimension, tessellation, self.log_level)
                all_lengths = [v[2] for v in junction_lengths]
                hist, bins = np.histogram(all_lengths, bins= number_of_bins, density= True)
                data_dictionary['3'] = junction_lengths                     # Key as integer refers to the integer representing the characteristic feature

            elif parameter is 'junction_angle':
                junction_angles_degrees = junction_angle(dimension, tessellation, self.log_level)
                all_angles = [v[2::2] for v in junction_angles_degrees]
                all_angles_flatten = np.around(np.array([w for v in all_angles for w in v]), decimals=8)
                hist, bins = np.histogram(all_angles_flatten, bins= number_of_bins, density= True)
                data_dictionary['4'] = junction_angles_degrees                     # Key as integer refers to the integer representing the characteristic feature

            elif parameter is 'distance_btw_grains':
                distance_btw_grain_array, distance_btw_grain_1d, smallest_distance_btw_neighbors_1d = distance_btw_grains(dimension, limit, tessellation, self.log_level)
                all_distances = distance_btw_grain_1d
                hist, bins = np.histogram(all_distances, bins= number_of_bins, density= True)
                data_dictionary['5a'] = distance_btw_grain_1d                     # Key as integer refers to the integer representing the characteristic feature
                data_dictionary['5b'] = distance_btw_grain_array                     # Key as integer refers to the integer representing the characteristic feature
                data_dictionary['5c'] = smallest_distance_btw_neighbors_1d      # smallest distance between neighbors

            elif parameter is 'disorientation_angles':
                disorientation_angle, orientation_data = disorientation_angles(dimension, limit, skewed_boundary_flag, required_texture, rand_quat_flag, orientation_data, tessellation, self.log_level)
                all_disorientation_angles = [v[2] for v in disorientation_angle]
                hist, bins = np.histogram(all_disorientation_angles, bins= number_of_bins, density= True)
                data_dictionary['6'] = disorientation_angle                     # Key as integer refers to the integer representing the characteristic feature

            elif parameter is 'type_of_grain_boundary':
                parent_function_name = inspect.stack()[1][3]
                type_of_grain_boundaries, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, dimension, limit, skewed_boundary_flag, self.log_level)
                all_grain_boundaries = type_of_grain_boundaries[:, 6]
                hist, bins = np.histogram(all_grain_boundaries, bins = number_of_bins, density= True)
                data_dictionary['7'] = type_of_grain_boundaries                     # Key as integer refers to the integer representing the characteristic feature

            elif parameter is 'schmid_factor':
                schmid_factors, single_gr_single_sys_schmid_factor, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, limit, stress_direction, slip_system_family, crystal_symmetry_type, orientation_data, tessellation, self.log_level)
                all_schmid_factors = schmid_factors[:, 1]
                hist, bins = np.histogram(all_schmid_factors, bins = number_of_bins, density= True)
                data_dictionary['8'] = schmid_factors                     # Key as integer refers to the integer representing the characteristic feature

            ## Extracting User Distribution data
            user_x_data = user_data[start_row_of_parameter[index]:start_row_of_parameter[index + 1], 0]
            user_y_data = user_data[start_row_of_parameter[index]:start_row_of_parameter[index + 1], 1]

            ## Extracting predicted distribution data
            predicted_x_data = bins[: -1]
            predicted_y_data = hist

            ## Finding minimum and maximum values for both User and predicted data along X axis
            min_user_x, max_user_x = np.min(user_x_data), np.max(user_x_data)    
            min_predicted_x, max_predicted_x = np.min(predicted_x_data), np.max(predicted_x_data)

            ## Finding the lowest and Highest values along X axis
            min_x = np.min([min_user_x, min_predicted_x])
            max_x = np.max([max_user_x, max_predicted_x])

            ############################################################################
            # Setting cost function to be evaluated at 1000 points
            no_evaluation_points = 1000
            ############################################################################
            
            ## Creating equally spaced evaluation points along X axis
            range_of_x = np.linspace(min_x, max_x, no_evaluation_points)

            ## Interpolating user and predicted data based on Existing values for values of Y axis along equally spaced X axis range
            user_y_data_interpolated = np.interp(range_of_x, user_x_data, user_y_data, left=0, right=0)
            predicted_y_data_interpolated = np.interp(range_of_x, predicted_x_data, predicted_y_data, left=0, right=0)

            ## Checking if size of User and predicted data is the same
            assert len(user_y_data_interpolated) == no_evaluation_points
            assert len(predicted_y_data_interpolated) == no_evaluation_points
            
            ## Combining vectors of X axis and Y axis into 2D array
            combined_user_local_array = np.hstack((np.reshape(range_of_x, (-1, 1)), np.reshape(user_y_data_interpolated, (-1, 1))))
            combined_predicted_local_array = np.hstack((np.reshape(range_of_x, (-1, 1)), np.reshape(predicted_y_data_interpolated, (-1, 1))))
            
            ## appending to main array
            combined_user_data = np.concatenate((combined_user_data, combined_user_local_array))
            combined_predicted_data = np.concatenate((combined_predicted_data, combined_predicted_local_array))
            start_row_combined_data.append(start_row_combined_data[-1] + no_evaluation_points)
        
        ## Re-setting log_level to original 
        self.log_level = original_log_level
        
        ## Finding the cost function value using the function
        ## Function name is identified based on the key of the dictionary 'cost_function_names'
        cost_function_value = cost_function_names[func_name_key](combined_user_data, combined_predicted_data, start_row_combined_data, data_dictionary, args)

        ## Storing the arrays of data for smallest cost function value achieved
        if cost_function_value < self.smallest_cost_function_value:
            self.final_user_distribution_data = combined_user_data
            self.final_predicted_data = combined_predicted_data
            self.final_start_row_combined_data = start_row_combined_data
            self.smallest_cost_function_value = cost_function_value
        
        if len(self.func_eval_count) == 1:                                      # again changed to 1 as with 2 results are not appropriate # set to 2 since the bins are of same values when equally spaced seeds are used at func_eval_count=1
            self.initial_distribution = combined_predicted_data
            print("Iteration number, Cost function value: \n")

        ## Assigning values for Dynamic Plot
        self.func_eval_count.append(self.func_eval_count[-1] + 1)
        self.current_cost_function_value.append(cost_function_value)
        
        print(self.func_eval_count[-1], cost_function_value)

        ## Adding current iteration seeds array to a list of seeds array of all iteration
        self.seeds_array_all_iterations.append(seed_array)

        log.debug('Iteration: ' +str(self.func_eval_count[-1]) + ', Cost function value: ' +str(self.current_cost_function_value[-1]))

        ## Updating plot and saving seeds data
        previous_iter = self.func_eval_count[-1]-1
        if previous_iter % self.save_interval == 0:                               # self.save_interval denotes the interval in which axes is updated
            ax_animate.set_xlim(0, self.func_eval_count[-1] + self.save_interval)                # X axis is extended with more self.save_interval units
            ax_animate.set_ylim(0, np.max(self.current_cost_function_value) + 10)                 # Y axis is updated with max value of cost value achieved where 10 is used randomly so that the upper portion of plot is not congested
            fig_animate.canvas.draw()
            fig_animate.canvas.flush_events()

            ## saving all iterations seed data current interval
            output_file_path = Path(self.store_folder, self.material, self.now, "Text_output", "seed_data_all_func_eval.txt")
            
            if (previous_iter) == 0:
                output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
                output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
                output_file_path.parent.parent.parent.mkdir(exist_ok=True)
                output_file_path.parent.parent.mkdir(exist_ok=True)
                output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

                with open(str(output_file_path), 'a+') as f:
                    f.truncate(0)

                    f.write("""# seed_data_all_func_eval.txt
# 
# This was created as part of "Personal Programming Project (PPP)" coursework in 
# "Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
# Germany.
# 
# For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>
# 
# @authors: Serrao Prince Henry, Arun Prakash
# @email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
# Copyright © 2020 by Serrao Prince Henry, Dr. Arun Prakash
# 
# This file is part of OptiMic.
# 
# OptiMic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# OptiMic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with OptiMic.  If not, see <https://www.gnu.org/licenses/>. \n""")
        
                    f.write("################################################################################\n\n")
        
                    f.write('# Date and time [YYYYMMDD_HHMMSS]: ' + self.now)

                    f.write('\n\n# Command line input: ' + " ".join(sys.argv))

                    f.write("""\n\n# This file contains data related to:
#        1. Grain seeds coordinates and grain orientation Quaternions predicted 
#           during each objective function calls\n\n""")
                    
                    f.write("# Seeds data (Seed coordinates + Orientation as Quaternions)\n")
                    f.write("# X coordinate, Y coordinate, Z coordinate, W, q1, q2, q3 \n")
            else:
                with open(str(output_file_path), 'a+') as f:
                    for i in range(previous_iter - self.save_interval, previous_iter):
                        f.write("\n# Objective function Evaluation count No.: " +str(i+1) + ", Cost Function value: " + str(self.current_cost_function_value[i]) + "\n")
                        f.write("# X coordinate, Y coordinate, Z coordinate, W, q1, q2, q3 \n")
                        seed_data = np.array(self.seeds_array_all_iterations[i])
                        temp_write_data = np.around(np.concatenate((seed_data, orientation_data), axis=1), decimals=14)
                        np.savetxt(f, temp_write_data, delimiter=',', comments='#', fmt='%.14f')

        line.set_ydata(self.current_cost_function_value)                        # Updating Y data
        line.set_xdata(self.func_eval_count[1:])                               # Updating X data
        ax_animate.draw_artist(ax_animate.patch)                                # Restoring plot region
        ax_animate.draw_artist(line)                                            # Plotting line
        fig_animate.canvas.blit()                                               # Updates only the line instead of entire plot comprising of title, axes, etc.
        fig_animate.canvas.flush_events()

        log.debug('Cost function was successfully evaluated at iteration no. ' + str(self.func_eval_count[-1]))
        return cost_function_value

def constraints_func(dimension, limits_wo_factors, log_level):
    """
    Define a list of constraints for each element of the prediction array.

    Parameter
    ---------
    dimension: integer
        Dimension of study. (2 or 3)

    limits_wo_factors: array
        Size of simulation box along X, Y & Z direction without factors.

    log_level: string
        Logger level to be used

    Returns
    -------
    List of all constraints.
    """
    log = set_logger(name_str, 'log_data.log', log_level)
    constraints = []
    '''
    Constraining the seed coordinates along all three directions such that 
    difference between the Upper limit of the dimension along particular direction 
    and the seed coordinate along that direction is greater than equal to 0. Apply 
    these constraints along all the directions respectively 
    '''
    if dimension == 3:
        lower = 0
        upper_x = limits_wo_factors[0]
        upper_y = limits_wo_factors[1]
        upper_z = limits_wo_factors[2]
        lower_bound = {'type': 'ineq', 'fun': lambda x: np.array(x[:]) - lower}
        upper_bound_x = {'type': 'ineq', 'fun': lambda x: upper_x - np.array(x[0::3])}
        upper_bound_y = {'type': 'ineq', 'fun': lambda x: upper_y - np.array(x[1::3])}
        upper_bound_z = {'type': 'ineq', 'fun': lambda x: upper_z - np.array(x[2::3])}
        constraints.append(lower_bound)
        constraints.append(upper_bound_x)
        constraints.append(upper_bound_y)
        constraints.append(upper_bound_z)
    
        '''
        When the dimension is 2, it is required that the coordinate along Z axis of 
        the seeds remain the same. Since COBYLA only supports inequality constraints,
        we have modified to inequality constraint as
        0 <= (difference of limit along Z and the seed coordinate along Z) <= 0 and 
        this is possible only when difference is zero.
        '''
    elif dimension == 2:
        lower = 0
        upper_x = limits_wo_factors[0]
        upper_y = limits_wo_factors[1]
        upper_z = limits_wo_factors[2]
        lower_bound = {'type': 'ineq', 'fun': lambda x: np.array(x[:]) - lower}
        upper_bound_x = {'type': 'ineq', 'fun': lambda x: upper_x - np.array(x[0::3])}
        upper_bound_y = {'type': 'ineq', 'fun': lambda x: upper_y - np.array(x[1::3])}
        upper_bound_z_1 = {'type': 'ineq', 'fun': lambda x: upper_z - np.array(x[2::3])}
        upper_bound_z_2 = {'type': 'ineq', 'fun': lambda x: np.array(x[2::3]) - upper_z}
        constraints.append(lower_bound)
        constraints.append(upper_bound_x)
        constraints.append(upper_bound_y)
        constraints.append(upper_bound_z_1)
        constraints.append(upper_bound_z_2)
    
    log.info('Constraints were successfully assigned')

    return constraints

# ## Array of all the available texture information from the research paper
# available_required_texture_array = copy.deepcopy(available_required_texture)

## Using Click group to execute either the tests or the main program
## This function is just to guide the execution of required function based on input argument
@click.group()
def guide():
    print("==================================")
    print("  ___        _   _ __  __ _   \n\
 / _ \ _ __ | |_(_)  \/  (_) ___  \n\
| | | | '_ \| __| | |\/| | |/ __| \n\
| |_| | |_) | |_| | |  | | | (__  \n\
 \___/| .__/ \__|_|_|  |_|_|\___| \n\
      |_|                        ")
    print("==================================")
    print("Optimized Microstructure Generator")
    print("Copyright © 2020-2021")
    print("Authors: P. H. Serrao, S. Sandfeld, A. Prakash")
    print("Version: "+ str(version)+"\n\n")
    log = set_logger(name_str, 'log_data.log', 'INFO')
    #log.info('\n')
    log.info("Welcome to Optimized Microstructure Generator !!")
    log.info("Version: " + str(version))
    log.info("Command-line input: " + " ".join(sys.argv))
    pass

############# REMEMBER PROMPT HAS AN ISSUE WITH NARGS != 1################
## The main program which would optimise the random seed coordinates based on user requirements
@guide.command()
@click.option('-s', '--size', required=True, help='Size of simulation box along X, Y & Z direction in the format n n n', type= float, nargs= 3)
@click.option('-dim', '--dimension', help='Dimension of study ie; 2D or 3D in the format n where n is an integer', type= int, nargs= 1)
@click.option('-n', '--number_seed', help='Number of seeds/grains in the format n', type= int, nargs= 1)
@click.option('-t', '--target', help='Target distribution file name as a string stored in the same directory where the package is executed from', type= str, nargs= 1)
@click.option('-c', '--characteristic', help='The characteristic that has to be optimized in the format n where n corresponds to the integer number corresponding to the characteristic', type= int, multiple = True)#nargs= 1)
@click.option('-m', '--material', help='The name of the material as a string', type=str, nargs=1)
@click.option('-sdir', '--stress_direction', help='The direction of stress for computing the Schmid Factors', default= [1, 0, 0], show_default=True, type=int, nargs=3)
@click.option('-slip', '--slip_family', help='The slip system family having first 3 components denoting family of slip plane normal and remaining 3 components denoting family of slip direction', default= [1, -1, 1, 1, 1, 0], show_default=True, type=int, nargs=6)
@click.option('-ctype', '--crystal_type', help='The crystal symmetry type to be used for Schmid factor calculations. Either Cubic, Orthorhombic, Hexagonal or Tetragonal', default='CUBIC', show_default=True, type=str, nargs=1)
@click.option('-so', '--sharp_orientation', help='Required texture common to each grain in the format n n n as provided in the documentation', type=float, nargs=3)
# @click.option('-r', help='Flag to indicate if Random orientations is to be generated', is_flag=True)
@click.option('-noopti', '--no_optimization', help='Flag to indicate if optimization is not to be performed', is_flag=True)
@click.option('-f', '--face_flag', is_flag=True, help= 'This flag is to be used to indicate if a closed surface is to be used for visualization files for all grain in one file')
@click.option('-ss', '--seed_spacing', help='Option to indicate if randomly placed seeds are required or regularly spaced seeds ie; cubic_2d, cubic_3d, etc', type=str, nargs=1, show_default=True, default= 'random_3d')
@click.option('-sl', '--spacing_length', help='Option to specify the spacing length such that it is exactly a multiple of size of simulation box along all three directions', show_default=True, default= 1.0, type=float, nargs=1)
@click.option('-om', '--optimization_method', help='Method to be used for optimization', show_default=True, default='COBYLA', type=str, nargs=1)
@click.option('-skw', '--skew_boundary', help='Flag to indicate skewed boundary requirement', is_flag=True)
@click.option('-ucf', '--user_cost_func', help='Specify the user defined cost function file name. Refer documentation for file specifications', type=str, nargs=1)
@click.option('-msh', '--mesh', help='Flag to indicate type of meshing required of the simulation box. Eg. hex (for Hexahedral), tet (for Tetrahedral), vis (for Visualization)', type=str, nargs=1)
@click.option('-gms', '--mesh_size', help='Provide global mesh size', type=float, nargs=1, show_default=True, default= 0.5)
@click.option('-mi', '--max_iter', help='Provide maximum number of iterations during optimizytion. The objective function might be evaluated multiple times during each iteration.', show_default=True, default=200, type=int, nargs=1)
@click.option('-rs', '--rand_seed', help='Enter the seed value for Numpy random function', show_default=True, default=None, type=int)
@click.option('-nb', '--number_bins', help='Specify the number of bins', show_default=True, default=10, type=int)
@click.option('-si', '--save_interval', help='Intervals of iterations in which the seeds data is to be saved and live plot is to be extended', show_default=True, default=100, type=int)
@click.option('-r', '--restart', help='Restart optimization using seed positions of specified function evaluation count number. Specify -1 to use optimized seed positions.', show_default=True, default=0, type=int)
@click.option('-deb', '--debug', help='Flag to activate Debug mode', is_flag=True)
def main(size, dimension, number_seed, target, characteristic, material, stress_direction, slip_family, crystal_type, sharp_orientation, no_optimization, face_flag, seed_spacing, spacing_length, optimization_method, skew_boundary, user_cost_func, mesh, mesh_size, max_iter, rand_seed, number_bins, save_interval, restart, debug):
    """
    \f
    Function to parse command-line inputs of Click.

    Parameter \n
    --------- \n
    size: list of 3 elements \n
    \t Size of simulation box along X, Y & Z direction \n \n

    dimension: integer \n
    \t Dimension of study. (2 or 3)\n \n

    number_seed: integer \n
    \t Number of seeds/grains.\n \n

    target: string \n
    \t Target distribution filename.\n \n

    characteristic: list \n
    \t Single or multiple characteristics to be optimized. Characteristics are
    \t identified based on integers as follows: \n
    \t\t 0: Grain sizes \n
    \t\t 1: Number of neighbors \n
    \t\t 2: Grain boundary areas \n
    \t\t 3: Junction lengths \n
    \t\t 4: Junction angles in degrees \n
    \t\t 5: Distance between grains \n
    \t\t 6: Disorientation angles \n
    \t\t 7: Type of grain boundaries \n
    \t\t 8: Schmid factors \n
    
    material: string \n
    \t Name of the material of which microstructure is being generated. \n

    stress_direction: list of 3 elements \n
    \t Direction of loading \n

    slip_family: list of shape (6,) \n
    \t First three components denote family of slip plane normal and remaining
    \t three components denote family of slip direction. \n

    crystal_type: string \n
    \t Crystal symmetry to be used.(CUBIC, ORTHORHOMBIC, HEXAGONAL, TETRAGONAL) \n

    sharp_orientation: list of 3 elements \n
    \t Specific texture to be used for all grains \n
    
    no_optimization: boolean \n
    \t Flag to indicate that optimization is NOT required. \n
    
    face_flag: boolean \n
    \t Flag to indicate that opaque surface is to be used instead of transparent. \n
    
    seed_spacing: string \n
    \t Type of seed spacing to be used from the following:\n
    \t\t 1. cubic_2d \n
    \t\t 2. cubic_3d \n
    \t\t 3. fcc_2d \n
    \t\t 4. fcc_3d \n
    \t\t 5. bcc_3d \n
    \t\t 6. random_3d \n
    \t\t 7. Any filename to specify file containing user defined seeds \n
    
    spacing_length: float \n
    \t Spacing between seeds along X, Y & Z in 3D case and along X & Y  \n
    \t directions in Quasi-2D case. Also the spacing_length must be a perfect \n
    \t divisor of size of simulation box along all three directions.\n
    
    optimization_method: string \n
    \t Optimization algorithm to be used. (COBYLA, SLSQP, POWELL)
    
    skew_boundary: boolean \n
    \t Flag to indicate that skewed grain boundaries are required. Applicable
    \t only in Quasi-2D case. \n
    
    user_cost_func: string \n
    \t Filename of user defined cost function. \n
    
    mesh: string \n
    \t Specify the type of elements to be used from the following: \n
    \t\t 1. 'HEX' for hexahedral elements \n
    \t\t 2. 'TET' for tetrahedral elements \n
    \t\t 3. 'VIS' for visualization using tetrahedral elements for grains not  \n
    \t\t\t confined to simulation box.
    
    mesh_size: float \n
    \t Global mesh size to be used for meshing.\n

    max_iter: integer \n
    \t Maximum number of iterations during optimization. The objective \n
    \t function might be evaluated multiple times during each iteration \n
    
    rand_seed: integer \n
    \t Seed to be used for Numpy random so that same output/results can be repeated. \n
    
    number_bins: integer \n
    \t Total number of bins to be used while computing distribution. \n
    
    save_interval: integer \n
    \t The intervals in iteration at which the seeds data is to be stored and \n
    \t the live plot axes is to extended. \n

    restart: integer \n
    \t Restart optimization using the seed positions of the specified function \n
    \t evaluation count number. Specify -1 to use optimized seed positions. \n

    debug: boolean \n
    \t Flag to indicate if DEBUG mode is to be activated. \n

    Returns \n
    ------- \n
    Function returns nothing. 
    """

    ## This is done so that the function 'main_run()' can be imported in some another python script
    main_run(size, dimension, number_seed, target, characteristic, material, stress_direction, slip_family, crystal_type, sharp_orientation, no_optimization, face_flag, seed_spacing, spacing_length, optimization_method, skew_boundary, user_cost_func, mesh, mesh_size, max_iter, rand_seed, number_bins, save_interval, restart, debug)

def main_run(size, dimension, number_seed, target, characteristic, material, stress_direction, slip_family, crystal_type, sharp_orientation, no_optimization, face_flag, seed_spacing, spacing_length, optimization_method, skew_boundary, user_cost_func, mesh, mesh_size, max_iter, rand_seed, number_bins, save_interval, restart, debug):
    """
    Function to execute main tasks.

    Parameter \n
    --------- \n
    size: list of 3 elements \n
    \t Size of simulation box along X, Y & Z direction \n \n

    dimension: integer \n
    \t Dimension of study. (2 or 3)\n \n

    number_seed: integer \n
    \t Number of seeds/grains.\n \n

    target: string \n
    \t Target distribution filename.\n \n

    characteristic: list \n
    \t Single or multiple characteristics to be optimized. Characteristics are
    \t identified based on integers as follows: \n
    \t\t 0: Grain sizes \n
    \t\t 1: Number of neighbors \n
    \t\t 2: Grain boundary areas \n
    \t\t 3: Junction lengths \n
    \t\t 4: Junction angles in degrees \n
    \t\t 5: Distance between grains \n
    \t\t 6: Disorientation angles \n
    \t\t 7: Type of grain boundaries \n
    \t\t 8: Schmid factors \n
    
    material: string \n
    \t Name of the material of which microstructure is being generated. \n

    stress_direction: list of 3 elements \n
    \t Direction of loading \n

    slip_family: list of shape (6,) \n
    \t First three components denote family of slip plane normal and remaining
    \t three components denote family of slip direction. \n

    crystal_type: string \n
    \t Crystal symmetry to be used.(CUBIC, ORTHORHOMBIC, HEXAGONAL, TETRAGONAL) \n

    sharp_orientation: list of 3 elements \n
    \t Specific texture to be used for all grains \n
    
    no_optimization: boolean \n
    \t Flag to indicate that optimization is NOT required. \n
    
    face_flag: boolean \n
    \t Flag to indicate that opaque surface is to be used instead of transparent. \n
    
    seed_spacing: string \n
    \t Type of seed spacing to be used from the following:\n
    \t\t 1. cubic_2d \n
    \t\t 2. cubic_3d \n
    \t\t 3. fcc_2d \n
    \t\t 4. fcc_3d \n
    \t\t 5. bcc_3d \n
    \t\t 6. random_3d \n
    \t\t 7. Any filename to specify file containing user defined seeds \n
    
    spacing_length: float \n
    \t Spacing between seeds along X, Y & Z in 3D case and along X & Y  \n
    \t directions in Quasi-2D case. Also the spacing_length must be a perfect \n
    \t divisor of size of simulation box along all three directions.\n
    
    optimization_method: string \n
    \t Optimization algorithm to be used. (COBYLA, SLSQP, POWELL)
    
    skew_boundary: boolean \n
    \t Flag to indicate that skewed grain boundaries are required. Applicable
    \t only in Quasi-2D case. \n
    
    user_cost_func: string \n
    \t Filename of user defined cost function. \n
    
    mesh: string \n
    \t Specify the type of elements to be used from the following: \n
    \t\t 1. 'HEX' for hexahedral elements \n
    \t\t 2. 'TET' for tetrahedral elements \n
    \t\t 3. 'VIS' for visualization using tetrahedral elements for grains not  \n
    \t\t\t confined to simulation box.
    
    mesh_size: float \n
    \t Global mesh size to be used for meshing.\n

    max_iter: integer \n
    \t Maximum number of iterations during optimization. The objective \n
    \t function might be evaluated multiple times during each iteration \n
    
    rand_seed: integer \n
    \t Seed to be used for Numpy random so that same output/results can be repeated. \n
    
    number_bins: integer \n
    \t Total number of bins to be used while computing distribution. \n
    
    save_interval: integer \n
    \t The intervals in iteration at which the seeds data is to be stored and \n
    \t the live plot axes is to extended. \n

    restart: integer \n
    \t Restart optimization using the seed positions of the specified function \n
    \t evaluation count number. Specify -1 to use optimized seed positions. \n

    debug: boolean \n
    \t Flag to indicate if DEBUG mode is to be activated. \n

    Returns \n
    ------- \n
    Function returns nothing. 
    """

    ## Current time and date
    now = datetime.now()
    now = now.strftime("%Y%m%d_%H%M%S")

    log_level = 'INFO'
    if debug:
        log_level = 'DEBUG'
        
        ## Checking libraries
        check_libraries(log_level)
        
        ## starting yappi profiler
        yappi.start()

    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Executing MAIN program')
    log.info('Following input parameters will be used: \n' \
            + 'Size of simulation box: ' + str(size) + '\n' \
            + 'Dimension: ' +str(dimension) + '\n' \
            + 'Number of seeds: ' +str(number_seed) + '\n' \
            + 'Target distribution file name: ' +str(target) + '\n' \
            + 'Characteristic to be optimized: ' +str(characteristic) + '\n' \
            + 'Name of the material: ' +str(material) + '\n' \
            + 'Direction vector for stress applied: ' +str(stress_direction) + '\n' \
            + 'Slip system family (Slip plane normal + Slip direction): ' +str(slip_family) + '\n' \
            + 'Crystal symmetry type: ' +str(crystal_type) + '\n' \
            + 'Required sharp texture direction: ' +str(sharp_orientation) + '\n' \
            + 'Flag for NO optimization: ' + str(no_optimization) + '\n' \
            + 'Flag for opaque surface: ' + str(face_flag) + '\n' \
            + 'Seed spacing type or file name: ' +str(seed_spacing) + '\n' \
            + 'Spacing length between seeds: ' +str(spacing_length) + '\n' \
            + 'Optimization method/algorithm: ' +str(optimization_method) + '\n' \
            + 'Flag for skewed grain boundary: ' +str(skew_boundary) + '\n' \
            + 'Name of user defined cost function file: ' +str(user_cost_func) + '\n' \
            + 'Required meshing type: ' +str(mesh) + '\n' \
            + 'Global mesh size: ' +str(mesh_size) + '\n' \
            + 'Maximum number of function evaluations: ' +str(max_iter) + '\n' \
            + 'Seed of random generator: ' +str(rand_seed) + '\n' \
            + 'Number of bins for histogram: ' +str(number_bins) + '\n' \
            + 'Save interval: ' + str(save_interval) + '\n' \
            + 'Restart position: ' + str(restart) + '\n' \
            + 'Debug mode activated: ' + str(debug) + '\n')

    np.random.seed(rand_seed)                                                                    # Initializing seed to None
    log.debug("Rand seed set to " + str(rand_seed))
    
    ## The name of the folder where the output files would be stored
    store_folder = "visualization_output"
       
    ## Defining the input parameters obtained from bash command line
    limit = tuple(np.array(size))
    dimension = dimension
    number_of_seeds = number_seed
    target_distribution_file_name = target
    characteristic_to_be_optimized = np.array(characteristic)
    material = material
    stress_direction = np.array(stress_direction)
    
    ## Checking if the slip system family provided is orthogonal or not 
    slip_system_family = np.array(slip_family)
    dot_family = np.dot(slip_system_family[0:3], slip_system_family[3:])
    assert np.isclose(dot_family, 0.0), 'Please provide a slip system such that \
family of slip plane normal is orthogonal to family of slip direction. For example: --slip_family 1 -1 1 1 1 0'

    crystal_symmetry_type = crystal_type.upper()

    required_texture = np.array(sharp_orientation)
#   rand_quat_flag = r

    ## Checking if random orientations are to be assigned
    if len(required_texture) <= 0:
        rand_quat_flag = True
    else:
        rand_quat_flag = None

    no_optimization_flag = no_optimization
    face_flag = face_flag
    seed_spacing_type = seed_spacing
    spacing_length = spacing_length
    optimization_method = optimization_method.upper()
    skewed_boundary_flag = skew_boundary
    user_cost_function_name = user_cost_func
    
    if mesh:
        mesh_flag = mesh.upper()
    else: 
        mesh_flag = mesh

    if skew_boundary and dimension == 3:
        log.critical('You have entered dimension = 3 with skew grain boundary flag. Please check your inputs and try again. Exiting now...')
        raise ValueError('You have entered dimension = 3 with skew grain boundary flag. Please check your inputs and try again. Exiting now...')
		
    global_mesh_size = mesh_size
    max_iterations = max_iter

    number_of_bins = number_bins

    ## Adjusting size of simulation box as per the requirement
    assert (np.array(limit) > 0.0).all(), 'Size of simulation box cannot be negative or zero.'
    limit = list(limit)
    if str(seed_spacing_type).lower() in ['cubic_3d', 'bcc_3d', 'fcc_3d']:
        
        assert (restart == 0), 'Restart option can only be used with user specified seed positions'
        if not np.isclose((limit[0] % spacing_length), 0):
            limit[0] = limit[0] + (spacing_length - (limit[0]%spacing_length))
            assert np.isclose(limit[0]%spacing_length, 0)
        if not np.isclose(limit[1] % spacing_length, 0):
            limit[1] = limit[1] + (spacing_length - (limit[1]%spacing_length))
            assert np.isclose(limit[1]%spacing_length, 0)
        if not np.isclose((limit[2] % spacing_length), 0):
            limit[2] = limit[2] + (spacing_length - (limit[2]%spacing_length))
            assert np.isclose(limit[2]%spacing_length, 0)

    elif str(seed_spacing_type).lower() in ['cubic_2d', 'fcc_2d']:
        
        assert (restart == 0), 'Restart option can only be used with user specified seed positions'
        if not np.isclose((limit[0] % spacing_length), 0):
            limit[0] = limit[0] + (spacing_length - (limit[0]%spacing_length))
            assert np.isclose(limit[0]%spacing_length, 0)
        if not np.isclose(limit[1] % spacing_length, 0):
            limit[1] = limit[1] + (spacing_length - (limit[1]%spacing_length))
            assert np.isclose(limit[1]%spacing_length, 0)

    elif str(seed_spacing_type).lower() in ['hcp_2d']:
        
        assert (restart == 0), 'Restart option can only be used with user specified seed positions'
        if not np.isclose((limit[0] % spacing_length), 0):
            limit[0] = limit[0] + (spacing_length - (limit[0]%spacing_length))
            assert np.isclose(limit[0]%spacing_length, 0)
        if not np.isclose(limit[1] % (np.sqrt(3)*spacing_length), 0):
            hcp_y_spacing = np.sqrt(3)*spacing_length
            limit[1] = limit[1] + (hcp_y_spacing - (limit[1] % hcp_y_spacing))
            assert np.isclose(limit[1]%hcp_y_spacing, 0)

    elif str(seed_spacing_type).lower() in ['hcp_3d']:

        assert (restart == 0), 'Restart option can only be used with user specified seed positions'
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
    limit = tuple(limit)
    log.warning("Size of simulation box is set to " + str(limit))

    ## Creating list of all characteristics to be oprimized
    parameter_dict = {'0': 'grain_size_distribution', '1': 'number_of_neighbors', '2': 'grain_boundary_areas', '3': 'junction_length', '4': 'junction_angle', '5': 'distance_btw_grains', '6': 'disorientation_angles', '7': 'type_of_grain_boundary', '8': 'schmid_factor'}
    parameter_list = [parameter_dict[str(i)] for i in characteristic_to_be_optimized]
    
    log.info('Gathering seeds and orientation information')

    ##Generating seeds
    if str(seed_spacing_type).lower() == 'cubic_2D'.lower():
        #seed_array = np.zeros([(int(limit[0])) * (int(limit[1])), 3])
        seed_array_unique = np.around(cubic_lattice_2D(limit, spacing_length, log_level), decimals=14)
        orientation_data = None
        assert dimension == 2, 'Please enter dimension as 2 (--d 2 in command line input)'
        
    elif str(seed_spacing_type).lower() == 'cubic_3D'.lower():
        #seed_array = np.zeros([(int(limit[0]) + 1) * (int(limit[1]) + 1) * (int(limit[2]) + 1), 3])
        seed_array_unique = np.around(cubic_lattice_3D(limit, spacing_length, log_level), decimals=14)
        orientation_data = None
        assert dimension == 3, 'Please enter dimension as 3 (--d 3 in command line input)'

    elif str(seed_spacing_type).lower() == 'bcc_3D'.lower():
        #seed_array = np.zeros([((int(limit[0]) + 1) * (int(limit[1]) + 1) * (int(limit[2]) + 1)) + (int(limit[0]) * int(limit[1]) * int(limit[2])), 3])
        seed_array_unique = np.around(bcc_lattice_3D(limit, spacing_length, log_level), decimals=14)
        orientation_data = None
        assert dimension == 3, 'Please enter dimension as 3 (--d 3 in command line input)'

    elif str(seed_spacing_type).lower() == 'fcc_2D'.lower():
        #seed_array = np.zeros([((int(limit[0])) * (int(limit[1]))) + (int(limit[0]) * int(limit[1])), 3])
        seed_array_unique = np.around(fcc_lattice_2D(limit, spacing_length, log_level), decimals=14)
        orientation_data = None
        assert dimension == 2, 'Please enter dimension as 2 (--d 2 in command line input)'
    
    elif str(seed_spacing_type).lower() == 'fcc_3D'.lower():
        #seed_array = np.zeros([2*((int(limit[0])) * (int(limit[1])) + (int(limit[1]) * int(limit[2])) + (int(limit[0]) * int(limit[2])) + (int(limit[0])*int(limit[1])*int(limit[2])*3)), 3])
        seed_array_unique = np.around(fcc_lattice_3D(limit, spacing_length, log_level), decimals=14)
        orientation_data = None
        assert dimension == 3, 'Please enter dimension as 3 (--d 3 in command line input)'

    elif str(seed_spacing_type).lower() == 'hcp_2D'.lower():
        seed_array_unique = np.around(hcp_lattice_2D(limit, spacing_length, log_level), decimals=14)
        orientation_data = None
        assert dimension == 2, 'Please enter dimension as 2 (--d 2 in command line input)'

    elif str(seed_spacing_type).lower() == 'hcp_3D'.lower():
        seed_array_unique = np.around(hcp_lattice_3D(limit, spacing_length, log_level), decimals=14)
        orientation_data = None
        assert dimension == 3, 'Please enter dimension as 3 (--d 3 in command line input)'
    
    elif str(seed_spacing_type).lower() == 'random_3d'.lower():
        seed_coordinates_list = random_generator(number_of_seeds, dimension, limit, log_level)
        seed_array = np.array(seed_coordinates_list)
        new_array = [tuple(row) for row in seed_array]
        seed_array_unique = np.around(np.unique(new_array, axis = 0), decimals=14)
        orientation_data = None
        ## Checking if the number of random unique seeds are same as required
        assert seed_array_unique.shape[0] == number_of_seeds
        assert (restart == 0), 'Restart option can only be used with user specified seed positions'


    else:
        with open(seed_spacing_type, 'r') as f:
            all_data = np.loadtxt(f, delimiter=',', comments='#')
            assert ((len(all_data)%number_of_seeds) == 0), 'There is mis-match in the number of seed positions provided. Probably total required seeds number is changed.'
            reshape_data = np.reshape(all_data, (-1, number_of_seeds, all_data.shape[1]))
            
            data = reshape_data[restart]
            
            seed_array_unique = np.around(data[:, :3], decimals=14)
            
            ## Extracting orientations data if available
            if data.shape[1] == 7:
                orientation_data = np.around(data[:, 3:7], decimals=14)
            else: orientation_data = None

    ## Generating/assigning orientations for grains
    if orientation_data is None:
        ## Checking if random orientations are required or sharp texture is required
        if rand_quat_flag:
            orientation_data = np.around(random_quaternions_generator(len(seed_array_unique), log_level), decimals=14)
        else:
            orientation_data = np.around(sharp_texture_quaternions(len(seed_array_unique), required_texture, log_level), decimals=14)    # Assigning random quaternions to each grain

    ############ Optimization Flag ###############
    if no_optimization_flag:
        log.info('Continuing without performing optimization')
        execute_func(limit, dimension, limit, material, orientation_data, required_texture, rand_quat_flag, seed_array_unique, stress_direction, slip_system_family, crystal_symmetry_type, store_folder, face_flag, now, number_of_bins, skewed_boundary_flag, mesh_flag, global_mesh_size, log_level)
    else:
        log.info('Starting optimization process')
        ## Reading required user distribution
        start_row_of_parameter = [0]
        row_counter = 0        #  Using row_counter since enumerate will take into account the comments and the 'end' line as well
        user_data = []
        
        ## Extracting target distribution
        log.debug('Reading target distribution')
        with open(target_distribution_file_name, 'r') as f:
            target_data = f.readlines()
            for line in target_data:
                if '#' in line:
                    continue
                elif 'end' in line:
                    start_row_of_parameter.append(row_counter)
                else:
                    user_data.append(list(map(float, line.rstrip('\n').split(','))))
                    row_counter += 1
        
        ## Converting to array
        user_data = np.array(user_data) 

        ## Normalizing target distribution
        for index, parameter in enumerate(parameter_list):
            ## Extracting User Distribution data
            user_x_data = user_data[start_row_of_parameter[index]:start_row_of_parameter[index + 1], 0]                 # extracting x data for a particular characteristic feature
            user_y_data = user_data[start_row_of_parameter[index]:start_row_of_parameter[index + 1], 1]                 # extracting y data for a particular characteristic feature

            area_under_curve = np.sum((user_x_data[1:]-user_x_data[:-1])*user_y_data[:-1])                              # calculating area under curve
            normalized_y_data = user_y_data/area_under_curve                                                            # normalizing Y data

            user_data[start_row_of_parameter[index]:start_row_of_parameter[index + 1], 1] = normalized_y_data           # replacing og data with normalized data

            assert np.isclose(np.sum((user_x_data[1:]-user_x_data[:-1])*normalized_y_data[:-1]), 1.0)                   # Ensuring area under curve is 1.0

        log.info('Target distribution read successfully')

        ## Importing user defined cost function if name is specified else general cost function is used
        ########################################################################
        # Importing user defined cost function formula if specified. The import
        # is done based on the file name specified by user. The user has to store
        # the function formula within a python file (.py extension). Please refer
        # documentation for more details.
        ########################################################################
        #if user_cost_function_name is None:
        
        func_name_key = 'general_cost_func'
        ## Defining the dictionary consisting of the function call names to selected based on required cost function
        cost_function_names = {'general_cost_func': cost_function_general}
        
        if not (user_cost_function_name is None):
            
            log.debug('Importing user defined cost function')
            ## The user function is read from the current working directory
            current_working_directory = os.getcwd()                             # Path of the current working directory is identified

            ## Getting module spec (specifications such as name, loader, origin) based on path
            spec = importlib.util.spec_from_file_location(user_cost_function_name, current_working_directory + '/' + user_cost_function_name)
            ## Creating a new module based on spec
            user_module = importlib.util.module_from_spec(spec)
            ## Executing the module
            spec.loader.exec_module(user_module)
            ## Defining the dictionary key
            func_name_key = 'user_cost_func'
        
            ## Adding new key to the dictionary consisting of the function call names to selected based on required cost function
            cost_function_names['user_cost_func'] = user_module.function_formula

            log.info('Successfully imported user defined cost function')
        
        log.debug('Defining constraints')

        ## Extracting order of magnitude from limits
        limits_factors = 10**(np.floor(np.log10(limit)))
        limits_wo_factors = limit/limits_factors
        seed_array_unique_wo_factors = seed_array_unique/limits_factors

        ## Defining constraints for the seeds coordinates for optimizer        
        constraints_list = constraints_func(dimension, limits_wo_factors, log_level)

        seed_array_unique_flatten = seed_array_unique_wo_factors.flatten()
        
        ## For Dynamic Plot of evolution of cost function
        font_size_value = 40
        label_size = 25
        #plt.ion()                                                               # turn interactive mode on                               
        fig_animate, ax_animate = plt.subplots(figsize=(10,10))
        line, = ax_animate.plot([], c='C2', lw=2)
        ax_animate.set_xlabel('Number of Iterations', fontsize=0.5*font_size_value)
        ax_animate.set_ylabel('Cost Function Value', fontsize=0.5*font_size_value)
        ax_animate.title.set_text('Evolution of Cost Function Value [' + optimization_method + ' Algorithm]')
        ax_animate.tick_params(labelsize=0.5*label_size)
        fig_animate.canvas.draw()
        plt.show(block=False)   
        
        log.info('Initialized dynamic plot successfully')

        ## Defining args for cost function
        tessellation_dummy_variable = []                                        ## dummy variable that will be replaced in args with tessellations data during optimization
        args_list = [parameter_list, dimension, user_data, start_row_of_parameter, limit, number_of_bins, fig_animate, [line, ax_animate], cost_function_names, func_name_key, required_texture, rand_quat_flag, stress_direction, orientation_data, skewed_boundary_flag, tessellation_dummy_variable, limits_factors, slip_system_family, crystal_symmetry_type]

        log.info('Calling optimizer')

        ## Optimization Algorithms options dictionary
        ## np.random.uniform is important in rhobeg argument to ensure that initial translations of seeds do not lead to overlapping seeds
        algo_options_dict = {'COBYLA': {'rhobeg': np.random.uniform(low=0.4, high=0.6), 'maxiter': max_iterations, 'disp': True},
                            'SLSQP': {'maxiter': max_iterations, 'ftol': 1e-6, 'disp': True, 'eps': 0.1},  ##### eps of 1e-2
                            'POWELL':{'maxiter': max_iterations},
                            'NELDER-MEAD': {'maxiter': max_iterations},
                            'BFGS': {'maxiter': max_iterations, 'eps': 0.1},
                            'L-BFGS-B': {'maxiter': max_iterations, 'eps': 0.1},
                            'TRUST-CONSTR': {'maxiter': max_iterations},
                            'CG': {'maxiter': max_iterations, 'eps': 0.1},
                            'TNC': {'maxiter': max_iterations, 'eps': 0.1}}

        ## Calling Optimizer
        optimize_class_instance = optimize_class(store_folder, now, material, save_interval, log_level)
        #print(optimize_class_instance.cost_function(seed_array_unique_flatten, args_list)); exit()                    # command: time ./optimic main -s 10. 10. 10. -sl 1 -dim 3 -n 100 -t user_grain_size_distribution.txt -c 0 -m steel -sdir 0 1 0 -ss seed_data.txt -rs 1 -om cobyla -mi 10 -nb 10
        optimization_result = minimize(optimize_class_instance.cost_function, seed_array_unique_flatten, args=(args_list), callback=optimize_class_instance.callback_f, method=optimization_method, constraints=constraints_list, options=algo_options_dict[optimization_method])
        
        log.info('Finished with optimization')

        plt.close(fig_animate)                                                  # Animation plot is closed
        
        initial_distribution            = optimize_class_instance.initial_distribution
        final_user_distribution_data    = optimize_class_instance.final_user_distribution_data
        final_predicted_data            = optimize_class_instance.final_predicted_data
        final_start_row_combined_data   = optimize_class_instance.final_start_row_combined_data
        smallest_cost_function_value    = optimize_class_instance.smallest_cost_function_value
        func_eval_count                = optimize_class_instance.func_eval_count
        current_cost_function_value     = optimize_class_instance.current_cost_function_value
        seeds_array_all_iterations      = optimize_class_instance.seeds_array_all_iterations                                                             # turn interactive mode off
        
        ##Saving Cost function evolution Plot
        output_file_path = Path(store_folder, material, now, "Plots", "evolution_of_cost_function.png")
        output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

        ## Plotting the graph of evolution of cost function again to save the plot
        fig_animate, ax_animate = plt.subplots(figsize=(10,10))
        line, = ax_animate.plot(func_eval_count[1:], current_cost_function_value, c='C2', lw=2)
        ax_animate.set_xlabel('Number of Iterations', fontsize=0.5*font_size_value)
        ax_animate.set_ylabel('Cost Function Value', fontsize=0.5*font_size_value)
        ax_animate.title.set_text('Evolution of Cost Function Value [' + optimization_method + ' Algorithm]')
        ax_animate.tick_params(labelsize=0.5*label_size)
        fig_animate.savefig(str(output_file_path))

        print(optimization_result.message)
        optimized_seed_array_flatten = optimization_result.x
        optimized_seed_array = np.reshape(optimized_seed_array_flatten, (-1, 3)) * limits_factors
        
        ## Adjusting seed array based on dimension
        if dimension == 2:
            optimized_seed_array[:, 2] = limit[2]/2.0

        execute_func(limit, dimension, limit, material, orientation_data, required_texture, rand_quat_flag, optimized_seed_array, stress_direction, slip_system_family, crystal_symmetry_type, store_folder, face_flag, now, number_of_bins, skewed_boundary_flag, mesh_flag, global_mesh_size, log_level)

        ## Plotting Expected and predicted data of each parameter on same plot
        
        fig, ax = plt.subplots(ncols = len(parameter_list), figsize=(30,15), squeeze=False)
        for plot_number, (i, j) in enumerate(zip(final_start_row_combined_data[:-1], final_start_row_combined_data[1:])):
            ax[0, plot_number].plot(initial_distribution[i:j, 0], initial_distribution[i:j, 1], label='Initial distribution')
            ax[0, plot_number].plot(final_user_distribution_data[i:j, 0], final_user_distribution_data[i:j, 1], label='User required distribution')
            ax[0, plot_number].plot(final_predicted_data[i:j, 0], final_predicted_data[i:j, 1], label='Optimized distribution')
            ax[0, plot_number].set_xlabel(parameter_list[plot_number], fontsize=font_size_value)
            ax[0, plot_number].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
            ax[0, plot_number].legend(prop={'size': 0.5*font_size_value})
            ax[0, plot_number].tick_params(labelsize=label_size)

        ##Saving Plot
        output_file_path = Path(store_folder, material, now, "Plots", "optimization_result.png")
        output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one
        
        
        plt.subplots_adjust(wspace=0.5)
        plt.suptitle("User required and Optimized distribution", fontsize=60)
        plt.tight_layout
        fig.savefig(str(output_file_path))

        ##Saving seeds data of all iterations

        log.debug('Saving seeds data at all iterations into a text file')

        ## Saving seeds data
        output_file_path = Path(store_folder, material, now, "Text_output", "seed_data_all_func_eval.txt")
        output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

        with open(str(output_file_path), 'a+') as f:
            ## Writing seed data for iterations that were pending during optimization
            if func_eval_count[-1] % save_interval != 0:
                for i in range(func_eval_count[-1] - (func_eval_count[-1] % save_interval), func_eval_count[-1]):
                    f.write("\n# Objective function Evaluation count No.: " +str(i+1) + ", Cost Function value: " + str(current_cost_function_value[i]) + "\n")
                    f.write("# X coordinate, Y coordinate, Z coordinate, W, q1, q2, q3 \n")
                    seed_data = np.array(seeds_array_all_iterations[i])
                    temp_write_data = np.around(np.concatenate((seed_data, orientation_data), axis=1), decimals=14)
                    np.savetxt(f, temp_write_data, delimiter=',', comments='#', fmt='%.14f')
            
            ## Writing optimization data
            f.write("\n\n\n## Optimization Information \n")
            f.write("# Optimized cost function value = " + str(optimization_result.fun) + "\n")
            f.write("# message = " + str(optimization_result.message) + "\n")
            #f.write("# nit = " + str(optimization_result.nit) + "\n")
            f.write("# nfev = " + str(optimization_result.nfev)+ "\n")
            f.write("# Optimized Seeds data (Seed coordinates + Orientation as Quaternions)\n")
            f.write("# X coordinate, Y coordinate, Z coordinate, W, q1, q2, q3 \n")
            temp_write_data = np.around(np.concatenate((optimized_seed_array, orientation_data), axis=1), decimals=14)
            np.savetxt(f, temp_write_data, delimiter=',', comments='#', fmt='%.14f')



        log.info('Successfully saved seeds data at all iterations into a text file')

        log.info('Successfully completed optimization process')

    if debug:
        ## Saving profiler output
        stats = yappi.get_func_stats()

        output_file_path = Path(store_folder, material, now, "Text_output", "callgrind.yappi_profiler_output.prof")
        output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one


        stats.save(str(output_file_path), type='callgrind')

        log.info('Successfully saved profiler data file')
        yappi.get_thread_stats().print_all()


    ## Moving log file to appropriate directory in visualization_files
    output_file_path = Path(store_folder, material, now, "Text_output")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    os.system('mv log_data.lo* ' + str(output_file_path) + '/')

@guide.command()
@click.option('-nt', '--name', help = 'Enter the test case name that you want to execute. For eg; Cubic_2d, Cubic_3d, fcc_2d, fcc_3d, bcc_3d, random_3d, textural or ALL', prompt='Enter test case name')
@click.option('-f', '--face_flag', is_flag=True, help= 'This flag is to be used to indicate if a closed surface is to be used for visualization files for all grain in one file')
@click.option('-rs', '--rand_seed', help='Enter the seed value for Numpy random function', type=int, nargs=1, show_default=True, default=None)
@click.option('-deb', '--debug', help='Flag to activate Debug mode', is_flag=True)
def test(name, face_flag, rand_seed, debug):
    """
    \f
    Function to parse command-line input arguments of Click.

    Parameters \n
    ---------- \n
    
    name: string \n
    \t Name of the test case to be executed. Following are the available test case: \n
    \t\t 1. 'all' \n
    \t\t 2. 'cubic_2d' \n
    \t\t 3. 'cubic_3d' \n
    \t\t 4. 'fcc_2d' \n
    \t\t 5. 'fcc_3d' \n
    \t\t 6. 'bcc_3d' \n
    \t\t 7. 'random_3d' \n
    \t\t 8. 'one_seed' \n
    \t\t 9. 'two_seed' \n
    \t\t 10. 'textural' \n

    face_flag: boolean \n
    \t Flag to indicate that opaque surface is to be used instead of transparent. \n

    rand_seed: integer \n
    \t Seed to be used for Numpy random so that same output/results can be repeated. \n

    debug: boolean \n
    \t Flag to indicate if DEBUG mode is to be activated. \n

    Returns \n
    ------- \n
    Function returns nothing.
    """
    
    ## This is done so that the function 'test_run()' can be imported in some another python script
    test_run(name, face_flag, rand_seed, debug)

def test_run(name, face_flag, rand_seed, debug):
    """
    Function to execute test cases.

    Parameters \n
    ---------- \n
    
    name: string \n
    \t Name of the test case to be executed. Following are the available test cases: \n
    \t\t 1. 'all' \n
    \t\t 2. 'cubic_2d' \n
    \t\t 3. 'cubic_3d' \n
    \t\t 4. 'fcc_2d' \n
    \t\t 5. 'fcc_3d' \n
    \t\t 6. 'bcc_3d' \n
    \t\t 7. 'random_3d' \n
    \t\t 8. 'one_seed' \n
    \t\t 9. 'two_seed' \n
    \t\t 10. 'textural' \n

    face_flag: boolean \n
    \t Flag to indicate that opaque surface is to be used instead of transparent. \n

    rand_seed: integer \n
    \t Seed to be used for Numpy random so that same output/results can be repeated. \n

    debug: boolean \n
    \t Flag to indicate if DEBUG mode is to be activated. \n

    Returns \n
    ------- \n
    Function returns nothing.
    """
 
    log_level = 'INFO'
    if debug:
        log_level = 'DEBUG'
        check_libraries(log_level)

    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Executing TEST case/s')
    log.info('Following input parameters will be used: \n' \
            + 'Name of specific test case or all test cases: ' + str(name) + '\n' \
            + 'Flag for opaque surfaces: ' +str(face_flag) + '\n' \
            + 'Seed of random generator: ' +str(rand_seed))

    np.random.seed(rand_seed)                                                                    # Initializing seed to None
    test_func(name, face_flag, log_level)



if __name__ == '__main__':
    guide()
