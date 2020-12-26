import numpy as np
import copy

## User-defined cost function designed to be used for non-default structural_characteristics module and Scipy-Voronoi library ##

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
    
    CritDist_dist = 0.1*min(args[4][:args[1]])
    CritDist_edge = 0.05*min(args[4][:args[1]])

    ## Extracting data from the data_dictionary
    all_distances = np.array(data_dictionary['5b'])
    all_junction_lengths = np.array(data_dictionary['3'][:, 2])

    cost_value_dist = mathematical_function(C, M, all_distances, CritDist_dist)
    cost_value_edge = mathematical_function(C, M, all_junction_lengths, CritDist_edge)

    return 0.5*(cost_value_dist + cost_value_edge)
