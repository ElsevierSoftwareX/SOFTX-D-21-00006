import numpy as np
import copy

## User-defined cost function designed to be used for default modules ##
## Designed only for 2D ##

def mathematical_function(C, M, distArr, CritDist):
    """
    Input: Constants C & M, Array of lengths or distances, Critical length or 
    distance.

    Processing: Computes cost value based on exponential and linear functions. 
    Basically penalizes seeds or vertices too close to each other and rewards 
    seeds or vertices too far away from each other.

    Returns: Sum of both exponential and linear terms
    """
    distCost=np.exp(-C * (distArr-CritDist)) + M * (distArr-CritDist)
    #print('Critdist: ', CritDist, distCost); exit()
    return np.sum(distCost)

def function_formula(combined_user_data, combined_predicted_data, start_row_combined_data, data_dictionary, args):
    """
    Input: As specified in the documentation of Optimized Microstructure Generator
    under section user-defined cost function.

    Processing: Computes cost function by penalizing seeds & vertices that are 
    very close and attracts them when they are very far away.

    Returns: Sum of cost value obtained by considering both distance between grains
    and junction length.
    """
    C = 0.6
    M = 0.6
    
    #args_list = [parameter_list, dimension, user_data, start_row_of_parameter, limit, number_of_bins, fig_animate, ax_animate, cost_function_names, func_name_key, required_texture, rand_quat_flag, stress_direction, orientation_data, skewed_boundary_flag, tessellation]
    limit = args[4]
    tessellation = copy.deepcopy(args[15])

    CritDist_dist = 0.1*min(args[4][:args[1]])
    all_distances_array = np.array(data_dictionary['5b'])
    neighbors_list_data = tessellation['neighbors_list']

    ## Computing distance between neighbors
    all_distances = []                                                          # empty list where distance between neighbors would be stored
    for index_grain, grain in enumerate(neighbors_list_data):
        for index_neighbor, neighbor_number in enumerate(grain):
            if (index_grain != neighbor_number):
                all_distances.append(all_distances_array[index_grain, neighbor_number])

    cost_value_dist = mathematical_function(C, M, all_distances, CritDist_dist)

    CritDist_edge = 0.05*min(args[4][:args[1]])
    all_junction_lengths = np.array([v[2] for v in data_dictionary['3'] if v[2] != limit[2]]) #######REMEMBER#######
    cost_value_edge = mathematical_function(C, M, all_junction_lengths, CritDist_edge)

    return (cost_value_dist + cost_value_edge)