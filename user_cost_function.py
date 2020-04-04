import numpy as np

def function_formula(combined_user_data, combined_predicted_data, start_row_combined_data):
    user_y = combined_user_data[:, 1]
    predicted_y = combined_predicted_data[:, 1]
    value = 0
    for i, j in zip(start_row_combined_data[:-1], start_row_combined_data[1:]):
        value += (np.sum(np.square(user_y[i:j] - predicted_y[i:j])))
    return value