import numpy as np
import copy

## User-defined cost function designed to be used for default modules ##

def function_formula(combined_user_data, combined_predicted_data, start_row_combined_data, data_dictionary, args):
    """
    Input: As specified in the documentation of Optimized Microstructure Generator
    under section user-defined cost function.

    Processing: Computes cost function by penalizing prediction with volume fraction 
    higher than 50% and varies linearly for volume fraction between 50% to 100%.

    Returns: Sum of cost value based on linear function.
    """
    threshold = 50.0
    C = 10
    M = 0.1
    
    ## extracting data
    #args_list = [parameter_list, dimension, user_data, start_row_of_parameter, limit, number_of_bins, fig_animate, ax_animate, cost_function_names, func_name_key, required_texture, rand_quat_flag, stress_direction, orientation_data, skewed_boundary_flag, tessellation]
    all_CSL_array = np.array(data_dictionary['7'])
    CSL_type_array = all_CSL_array[:, 6]

    ## Computing volume fraction of Special Grain Boundaries
    number_sgb = np.count_nonzero(CSL_type_array)
    volume_fraction = (number_sgb/len(CSL_type_array))*100

    ## Computing cost function based on exponential function
    if volume_fraction > threshold:
        return C*np.exp(M*(100.0 - volume_fraction))
    return 10*C*(np.exp(M * threshold) + np.exp(M*volume_fraction))