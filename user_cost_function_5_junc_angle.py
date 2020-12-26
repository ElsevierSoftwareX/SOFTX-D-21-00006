import numpy as np

def function_formula(combined_user_data, combined_predicted_data, start_row_combined_data, data_dictionary, args):
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
    Scalar SSE cost value.
    """
    user_y = combined_user_data[:, 1]
    predicted_y = combined_predicted_data[:, 1]
    value = (np.sum(np.square((user_y - predicted_y)*100)))            # 100 is multiplied so that the data is scaled for appropriate optimization
    return value