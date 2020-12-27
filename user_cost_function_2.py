# -*- coding: utf-8 -*-
"""
user_cost_function_2.py

Module to demonstrate an example of user-defined cost function.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
created: 16 May 2020
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

import numpy as np
import copy

## User-defined cost function designed to be used for default modules ##
## Designed only for 2D ##

def mathematical_function(C, M, distArr, CritDist):
    """
    Compute cost based on mathematical formula dependent on exponential and 
    linear functions

    Parameters
    ----------
    C: float
        Model constant term (used for scaling data).

    M: float
        Model constant term (used for scaling data).

    distArr: array
        Array of distances or lengths.

    CritDist: float
        Critical threshold value.

    Returns
    -------
    Cost as sum of both exponential and linear terms.

    """
    distCost=np.exp(-C * (distArr-CritDist)) + M * (distArr-CritDist)
    #print('Critdist: ', CritDist, distCost); exit()
    return np.sum(distCost)

def function_formula(combined_user_data, combined_predicted_data, start_row_combined_data, data_dictionary, args):
    """
    Computing cost function by penalizing seeds & vertices that are 
    very close and attracts them when they are very far away.

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
            8. ax_animate
            9. cost_function_names
            10. func_name_key
            11. required_texture
            12. rand_quat_flag
            13. stress_direction
            14. orientation_data
            15. skewed_boundary_flag
            16. tessellation (DICTIONARY OF ALL TESSELLATION DATA)

    Returns
    -------
    Sum of cost value obtained by considering both distance between grains
    and edge length.

    """
    C = 0.6
    M = 0.6
    
    #args_list = [parameter_list, dimension, user_data, start_row_of_parameter, limit, number_of_bins, fig_animate, ax_animate, cost_function_names, func_name_key, required_texture, rand_quat_flag, stress_direction, orientation_data, skewed_boundary_flag, tessellation]
    limit = args[4]
    tessellation = copy.deepcopy(args[15])

    CritDist_dist = 0.1*min(args[4][:args[1]])
    all_distances_array = np.array(data_dictionary['5a'])

    cost_value_dist = mathematical_function(C, M, all_distances_array, CritDist_dist)

    CritDist_edge = 0.05*min(args[4][:args[1]])
    all_ridge_lengths = np.array([v[3]/limit[2] for v in data_dictionary['2']]) #######REMEMBER#######
    cost_value_edge = mathematical_function(C, M, all_ridge_lengths, CritDist_edge)

    return (cost_value_dist + cost_value_edge)