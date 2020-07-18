# -*- coding: utf-8 -*-
"""
stereographic_projections.py

Module to generate pole figures.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

@authors: Serrao Prince Henry, Arun Prakash
@email: 
created: 27 February 2020
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

## INPUT FROM BASH COMMAND IS TO BE PROVIDED AS: $ python stereographic_projections.py --p 1 1 1  --n 50
## For help: $ python stereographic_projections.py --help

import numpy as np
import math
import matplotlib.pyplot as plt
import click

from ..textural_characteristics import sharp_texture_quaternions as sharp_texture_quaternions
from ..textural_characteristics import random_quaternions_generator as random_quaternions_generator

import quaternion

## INPUT FROM BASH COMMAND IS TO BE PROVIDED AS: $ python stereographic_projections.py --p 1 1 1  --n 50

def euler (thet_1_deg, phhi_deg, thet_2_deg):
    """
    Input Parameters: All 3 Bunge angles in degrees in the order of rotation about z axis, Rotated x axis and rotated z axis respectively.
    Processing: The matrices of rotation of respective axis are multiplied together in the order of z , x, z. 
    Returns: Transformation matrix from Crystal to System coordinates  
    """
    
    thet_1 = math.radians(thet_1_deg)
    phhi = math.radians(phhi_deg)
    thet_2 = math.radians(thet_2_deg)
    
    z_1 = np.array([[math.cos(thet_1), math.sin(thet_1), 0 ], 
                [- math.sin(thet_1), math.cos(thet_1), 0], 
                [0, 0, 1]])
    
    x_2 = np.array([[1, 0, 0], 
                    [0, math.cos(phhi), math.sin(phhi)], 
                    [0, - math.sin(phhi), math.cos(phhi)]])
    
    z_2 = np.array([[math.cos(thet_2), math.sin(thet_2), 0 ], 
                    [- math.sin(thet_2), math.cos(thet_2), 0], 
                    [0, 0, 1]])
    
    rotation_matrix = z_2.dot(x_2).dot(z_1)
    return rotation_matrix

def bunge_angles (transformation_matrix):
    """
    Input parameters: Transformation Matrix
    Processing: Required elements are extracted from the transformation matrix and the corresponding bunge angles are calculated.
    Returns: All 3 Bunge angles in degrees in the order of rotation about z axis, Rotated x axis and rotated z axis respectively.
    """
    r_13 = transformation_matrix[0, 2]
    r_23 = transformation_matrix[1, 2]
    r_31 = transformation_matrix[2, 0]
    r_32 = transformation_matrix[2, 1]
    r_33 = transformation_matrix[2, 2]
    
    phhi = math.acos(r_33)
    phhi_degrees = math.degrees(phhi)
    thet_1 = math.degrees(math.atan((r_31 /(math.sin(phhi)))/(-r_32/math.sin(phhi))))
    thet_2 = math.degrees(math.atan((r_13 /(math.sin(phhi)))/(r_23/math.sin(phhi))))
    return thet_1, phhi_degrees, thet_2

def test__bunge_angles():
    
    transformation_matrix = euler(atheta_1_deg, aphi_1_deg, atheta_2_deg)
    
    bunge_angle_1, bunge_angle_2, bunge_angle_3 = bunge_angles (transformation_matrix)
    assert np.isclose(bunge_angle_1, atheta_1_deg)
    assert np.isclose(bunge_angle_2, aphi_1_deg)
    assert np.isclose(bunge_angle_3, atheta_2_deg)

def symmetric_poles(symmetry_operator, normalized_pole):
    """
    Input parameters: array of symmetric operators and the array of normalized pole
    Processing: The normalized pole is multipled with each symmetric operator and return as an array.
    Returns: All the symmetric normalized poles as an array.
    """    
        
    sym_poles = []
    for i in range(symmetry_operator.shape[0]):
        sym_variable = symmetry_operator[i].dot(normalized_pole.T)
        sym_poles.append(sym_variable)
    return sym_poles

def poles_in_sys_axis(transformation_matrix, sym_normalized_poles):
    """
    Input Parameters: array of transformation matrix and of all symmetric normalized poles
    Processing: all the normalized poles are multiplied with the transformation matrix.
    Returns: an array of all the transformed symmetric normalized poles.
    """
    rotate_sym_poles = []
    for i in range(len(sym_normalized_poles)):
        rotated_variable = (transformation_matrix.dot(sym_normalized_poles[i]))
        rotate_sym_poles.append(rotated_variable)
    return(rotate_sym_poles)

def stereographic_proj (rotated_sym_poles):
    """
    Input Parameters: an array of all the transformed symmetric poles.
    Processing: the coordinates of the poles on the projected planes are calculated using vector calculation. 
    Returns: an array of the coordinates of all the points of the stereographic projections.
    """
    x = []
    y = []
    for i in range(len(rotated_sym_poles)):
        value = rotated_sym_poles[i]
        for j in range(len(value)):
            p_x = value[j][0]
            p_y = value[j][1]
            p_z = value[j][2]
            if p_z < 0: 
                continue
            var_x = (p_x/(p_z+1))
            var_y = (p_y/(p_z+1))
            x.append(var_x)
            y.append(var_y)
    return x, y

def random_angles_generator(i):
    """
    Input Parameters: -
    Processing: random bunge angles are generated within the respective range using numpy.
    Returns: Randomly generated bunge angles.
    """
    angle_1 = np.random.uniform(low = 0, high = 360)
    angle_2 = np.random.uniform(low = 0, high = 180)
    angle_3 = np.random.uniform(low = 0, high = 360)
    print("Bunge Angles for Orientation ", i+1, " are", angle_1, angle_2, angle_3)
    return angle_1, angle_2, angle_3

@click.command()
@click.option('--p', help = 'pole to be plotted in the format n n n', type = int, nargs = 3)
@click.option('--n', help ='number of orientations', type = int)
def parameters(p, n):
    
    ## Taking input from the command like and generating random euler angles based on number of orientations
    pole = np.array([p])
    number_of_orientations = n
    print("Pole is: ", pole[0])
    print("Number of orientations is/are: ", n)
    
    ## Normalizing Pole
    normalized_pole = (1/np.linalg.norm(pole))*pole
    
    ## Defining an array of the symmetric operators.
    symmetry_operator = np.array([[[1, 0, 0],
                                    [0, 1, 0], 
                                    [0, 0, 1]], 
                                    
                                    [[0, 0, 1], 
                                    [1, 0, 0], 
                                    [0, 1, 0]], 
                                    
                                    [[0, 1, 0], 
                                    [0, 0, 1], 
                                    [1, 0, 0]],
                                    
                                    [[0, -1, 0], 
                                    [0, 0, 1], 
                                    [-1, 0, 0]], 
                                    
                                    [[0, -1, 0], 
                                    [0, 0, -1], 
                                    [1, 0, 0]], 
                                    
                                    [[0, 1, 0], 
                                    [0, 0, -1], 
                                    [-1, 0, 0]], 
                                    
                                    [[0, 0, -1], 
                                    [1, 0, 0], 
                                    [0, -1, 0]], 
                                    
                                    [[0, 0, -1], 
                                    [-1, 0, 0], 
                                    [0, 1, 0]],  
                                    
                                    [[0, 0, 1], 
                                    [-1, 0, 0], 
                                    [0, -1, 0]], 
                                    
                                    [[-1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, -1]],
                                    
                                    [[-1, 0, 0], 
                                    [0, -1, 0],
                                    [0, 0, 1]], 
                                    
                                    [[1,0, 0], 
                                    [0, -1, 0], 
                                    [0, 0, -1]], 
                                    
                                    [[0, 0, -1], 
                                    [0, -1, 0], 
                                    [-1, 0, 0]], 
                                    
                                    [[0, 0, 1], 
                                    [0, -1, 0], 
                                    [1, 0, 0]], 
                                    
                                    [[0, 0, 1], 
                                    [0, 1, 0], 
                                    [-1, 0, 0]], 
                                    
                                    [[0, 0, -1], 
                                    [0, 1, 0], 
                                    [1, 0, 0]], 
                                    
                                    [[-1, 0, 0], 
                                    [0, 0, -1], 
                                    [0, -1, 0]], 
                                    
                                    [[1, 0, 0], 
                                    [0, 0, -1], 
                                    [0, 1, 0]], 
                                    
                                    [[1, 0, 0], 
                                    [0, 0, 1], 
                                    [0, -1, 0]], 
                                    
                                    [[-1, 0, 0], 
                                    [0, 0, 1], 
                                    [0, 1, 0]], 
                                    
                                    [[0, -1, 0], 
                                    [-1, 0, 0], 
                                    [0, 0, -1]], 
                                    
                                    [[0, 1, 0], 
                                    [-1, 0, 0], 
                                    [0, 0, 1]],
                                    
                                    [[0, 1, 0], 
                                    [1, 0, 0], 
                                    [0, 0, -1]], 
                                    
                                    [[0, -1, 0], 
                                    [1, 0, 0], 
                                    [0, 0, 1]]])
    
    ## Calculate symmetric pole array
    sym_normalized_poles = symmetric_poles(symmetry_operator, normalized_pole)
    
    ## Rotating all symmetric Poles
    rotated_sym_poles = []
    #quaternions_array = random_quaternions_generator(n)
    quaternions_array = sharp_texture_quaternions(n, np.array([1, 1, 1]))
    for i in range(number_of_orientations):
        quaternion_grain = quaternion.from_float_array(quaternions_array[i, :])
        transformation_matrix  = quaternion.as_rotation_matrix(quaternion_grain)
        rotated_sym_poles_iter  = poles_in_sys_axis(transformation_matrix , sym_normalized_poles)
        rotated_sym_poles.append(rotated_sym_poles_iter)
    
    ## Calculating stereographic projections
    x_plot , y_plot  = stereographic_proj (rotated_sym_poles)
    
    ## Plotting
    fig, ax = plt.subplots(figsize = (12, 12))
    r = 1
    t = np.arange(0, 2*math.pi, 0.01)                                             # For plotting a Circle      
    x_c = r*np.sin(t)                                                             # For plotting a Circle 
    y_c = r*np.cos(t)                                                             # For plotting a Circle
    ax.plot(x_c, y_c, linewidth = 5, color ='green', label = 'Plane (001)')
    ax.scatter(x_plot , y_plot  , label = 'Stereographic Projections')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    plt.legend()
    plt.show()

if __name__ == '__main__':
        parameters()
