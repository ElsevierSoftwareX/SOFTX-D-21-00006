import numpy as np
import copy

## User-defined cost function designed to be used for default modules ##

def function_formula(combined_user_data, combined_predicted_data, start_row_combined_data, data_dictionary, args):
    """
    Input: As specified in the documentation of Optimized Microstructure Generator
    under section user-defined cost function.

    Processing: Computes cost function by penalizing prediction with area fraction 
    lower than 50% and varies linearly for area fraction between 50% to 100%.

    Returns: Sum of cost value based on linear function.
    """
    threshold = 50.0
    C = 10
    M = 0.1
    
    ## extracting data
    #args_list = [parameter_list, dimension, user_data, start_row_of_parameter, limit, number_of_bins, fig_animate, ax_animate, cost_function_names, func_name_key, required_texture, rand_quat_flag, stress_direction, orientation_data, skewed_boundary_flag, tessellation]
    all_CSL_array = np.array(data_dictionary['7'])
    CSL_type_array = all_CSL_array[:, 6]
    CSL_gb_area_array = all_CSL_array[:, 7]
    total_gb_area = np.sum(CSL_gb_area_array)

    ## Computing volume fraction of Special Grain Boundaries
    indices_sgb = np.nonzero(CSL_type_array)
    sgb_gb_area = np.sum(CSL_gb_area_array[indices_sgb])
    area_fraction = (sgb_gb_area/total_gb_area)*100

    ## Computing cost function based on exponential function
    if area_fraction > threshold:
        return C*np.exp(M*(100.0 - area_fraction))
    return 10*C*(np.exp(M * threshold) + np.exp(M*(100-area_fraction)))