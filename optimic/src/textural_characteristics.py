# -*- coding: utf-8 -*-
"""
textural_characteristics.py

Module to compute following textural characteristics:
1. Disorientation Angles
2. Type of Grain boundary
3. Schmid factors

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
created: 07  December 2019
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

from optimic.src.main_import_statements import *
from scipy.spatial.transform import Rotation as R

from optimic.src.structural_characteristics import *
from optimic.src.set_logger import set_logger 
name_str = __name__
   
## Array of all the available texture information from the research paper
available_required_texture =  np.array([[1, 0, 0],
                                        [1, 1, 0],
                                        [1, 1, 1],
                                        [2, 1, 0],
                                        [2, 1, 1],
                                        [2, 2, 1],
                                        [3, 1, 0],
                                        [3, 1, 1],
                                        [3, 2, 0],
                                        [3, 2, 1],
                                        [3, 2, 2],
                                        [4, 1, 0],
                                        [4, 1, 1],
                                        [3, 3, 1],
                                        [4, 2, 1],
                                        [3, 3, 2],
                                        [4, 3, 0],
                                        [4, 3, 1],
                                        [5, 1, 0],
                                        [5, 1, 1],
                                        [4, 3, 2],
                                        [5, 2, 0],
                                        [5, 2, 1],
                                        [4, 4, 1],
                                        [5, 2, 2],
                                        [4, 3, 3],
                                        [5, 3, 0],
                                        [6, 1, 0],
                                        [5, 3, 2],
                                        [6, 1, 1],
                                        [4, 4, 3],
                                        [5, 4, 0],
                                        [6, 2, 1],
                                        [5, 3, 1],
                                        [5, 3, 3],
                                        [5, 5, 1],
                                        [5, 4, 1],
                                        [5, 4, 2],
                                        [6, 3, 1],
                                        [6, 3, 2],
                                        [5, 4, 3],
                                        [7, 1, 0],
                                        [7, 1, 1],
                                        [5, 5, 3],
                                        [7, 3, 1],
                                        [6, 4, 1],
                                        [7, 2, 0],
                                        [5, 5, 2],
                                        [7, 2, 1],
                                        [5, 4, 4],
                                        [7, 2, 2],
                                        [7, 3, 0],
                                        [7, 3, 3],
                                        [7, 5, 1],
                                        [7, 5, 3],
                                        [9, 1, 1],
                                        [9, 3, 1],
                                        [7, 5, 5],
                                        [7, 7, 1],
                                        [7, 7, 3],
                                        [6, 4, 5],
                                        [9, 5, 1],
                                        [6, 5, 0],
                                        [9, 5, 3],
                                        [7, 7, 5],
                                        [7, 3, 2],
                                        [11, 1, 1],
                                        [6, 5, 1]])

## Array of all the available texture information from the book by V. Randle in the normalized form 
norm_vectors_available_data = np.transpose([np.linalg.norm(available_required_texture, axis=1)])
broadcast_norm_vectors_available_data = np.tile(norm_vectors_available_data, (1, 3))
norm_available_required_texture = available_required_texture / broadcast_norm_vectors_available_data

## Defining the Cubic symmetry operators
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

## Defining the orthorhombic symmetry operators
orthorhombic_symmetry_operators = np.array([[[1, 0, 0],
                                             [0, 1, 0],
                                             [0, 0, 1]],

                                            [[-1, 0, 0],
                                             [0, 1, 0],
                                             [0, 0, -1]],

                                            [[-1, 0, 0],
                                             [0, -1, 0],
                                             [0, 0, 1]],

                                            [[1, 0, 0],
                                             [0, -1, 0],
                                             [0, 0, -1]]])

## Defining the hexagonal symmetry operators
ha = np.sqrt(3)/2.0
hexagonal_symmetry_operators = np.array([[[1, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1]],
                                         
                                         [[-.5, ha, 0],
                                          [-ha, -.5, 0],
                                          [0, 0, 1]],

                                         [[-.5, -ha, 0],
                                          [ha, -.5, 0],
                                          [0, 0, 1]],

                                         [[.5, ha, 0],
                                          [-ha, .5, 0],
                                          [0, 0, 1]],

                                         [[-1, 0, 0],
                                          [0, -1, 0],
                                          [0, 0, 1]],

                                         [[.5, -ha, 0],
                                          [ha, .5, 0],
                                          [0, 0, 1]],

                                         [[-.5, -ha, 0],
                                          [-ha, .5, 0],
                                          [0, 0, -1]],

                                         [[1, 0, 0],
                                          [0, -1, 0],
                                          [0, 0, -1]],

                                         [[-.5, ha, 0],
                                          [ha, .5, 0],
                                          [0, 0, -1]],

                                         [[.5, ha, 0],
                                          [ha, -.5, 0],
                                          [0, 0, -1]],

                                         [[-1, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, -1]],

                                         [[.5, -ha, 0],
                                          [-ha, -.5, 0],
                                          [0, 0, -1]]])

## Defining tetragonal symmetry operators
tetragonal_symmetry_operators = np.array([[[1, 0, 0],
                                           [0, 1, 0],
                                           [0, 0, 1]],

                                          [[-1, 0, 0],
                                           [0, 1, 0],
                                           [0, 0, -1]],

                                          [[1, 0, 0],
                                           [0, -1, 0],
                                           [0, 0, -1]],

                                          [[-1, 0, 0],
                                           [0, -1, 0],
                                           [0, 0, 1]],

                                          [[0, 1, 0],
                                           [-1, 0, 0],
                                           [0, 0, 1]],

                                          [[0, -1, 0],
                                           [1, 0, 0],
                                           [0, 0, 1]],

                                          [[0, 1, 0],
                                           [1, 0, 0],
                                           [0, 0, -1]],

                                          [[0, -1, 0],
                                           [-1, 0, 0],
                                           [0, 0, -1]]])


## Transforming rotation matrices of each symmetry operator to quaternions
symmetric_operators_quaternions = np.array([quaternion.from_rotation_matrix(r) for r in symmetry_operator])

type_of_csl_data = {'[1.0, 0.0, 0.0]': np.array([[5, 13, 17, 25, 29],
                                            [36.9, 22.6, 28.1, 16.3, 43.6]]),
        '[0.7071, 0.7071, 0.0]': np.array([[3, 9, 11, 17, 19, 27],
                                        [70.5, 38.9, 50.5, 86.6, 26.5, 31.6]]),
        '[0.5774, 0.5774, 0.5774]': np.array([[3, 7, 13, 19, 21, 31],
                                        [60, 38.2, 27.8, 46.6, 21.8, 17.9]]),
        '[0.8944, 0.4472, 0.0]': np.array([[3, 5, 7, 9, 15, 21, 23, 27, 29],
                                        [131.8, 180, 73.4, 96.4, 48.2, 58.4, 163.0, 35.4, 112.3]]),
        '[0.8165, 0.4082, 0.4082]': np.array([[3, 5, 7, 11, 15, 21, 25, 29, 31],
                                                [180, 101.5, 135.6, 63.0, 78.5, 44.4, 156.9, 149.6, 52.2]]),
        '[0.6667, 0.6667, 0.3333]': np.array([[5, 9, 9, 13, 17, 25, 29],
                                            [143.1, 90, 180, 112.2, 61.9, 73.7, 46.4]]),
        '[0.9487, 0.3162, 0.0]': np.array([[5, 7, 11, 13, 19, 23],
                                    [180, 115.4, 144.9, 76.7, 93.0, 55.6]]),
        '[0.9045, 0.3015, 0.3015]': np.array([[3, 5, 9, 11, 15, 15, 23, 25, 27, 31],
                                        [146.4, 95.7, 67.1, 180, 50.7, 117.8, 40.5, 168.3, 79.3, 126.6]]),
        '[0.8321, 0.5547, 0.0]': np.array([[7, 11, 13, 17, 19, 29, 31],
                                [149.0, 100.5, 180, 122.0, 71.6, 84.1, 54.5]]),
        '[0.8018, 0.5345, 0.2673]': np.array([[7, 9, 15, 15, 23, 25],
                                    [180.0, 123.8, 86.2, 150.1, 102.6, 63.9]]),
        '[0.7276, 0.4851, 0.4851]': np.array([[9, 13, 17, 21, 21],
                                    [152.7, 107.9, 180, 128.3, 79.0]]),
        '[0.9701, 0.2425, 0.0]': np.array([[9, 13, 17, 21, 21],
                                    [152.7, 107.9, 180, 79.0, 128.3]]),
        '[0.9428, 0.2357, 0.2357]': np.array([[9, 11, 17, 19, 27, 27],
                                [180, 129.5, 93.4, 153.5, 109.5, 70.5]]),
        '[0.6882, 0.6882, 0.2294]': np.array([[5, 7, 11, 17, 19, 23, 25],
                                    [154.2, 110.9, 82.2, 63.8, 180, 130.7, 51.7]]),
        '[0.8729, 0.4364, 0.2182]': np.array([[11, 15, 21, 23, 25],
                                [155.4, 113.6, 180, 85.0, 132.8]]),
        '[0.6396, 0.6396, 0.4264]': np.array([[11, 13, 19, 23, 29, 31],
                                    [180, 133.8, 99.1, 155.9, 76.0, 114.8]]),
        '[0.8, 0.6, 0.0]': np.array([[13, 17, 25, 25, 29],
                            [157.4, 118.1, 180, 90, 136.4]]),
        '[0.7845, 0.5883, 0.1961]': np.array([[13, 15, 21, 27, 31],
                            [180, 137.2, 103.8, 157.8, 80.7]]),
        '[0.9806, 0.1961, 0.0]': np.array([[13, 15, 21, 27, 31],
                            [180, 137.2, 103.8, 157.8, 80.7]]),
        '[0.9623, 0.1925, 0.1925]': np.array([[7, 9, 13, 19, 27, 27, 31],
                                    [158.2, 120.0, 92.2, 73.2, 60, 180, 137.9]]),
        '[0.7428, 0.5571, 0.3714]': np.array([[15, 19, 27, 29],
                                [159.0, 121.8, 94.3, 180]]),
        '[0.9285, 0.3714, 0.0]': np.array([[15, 19, 27, 29],
                                [159.0, 121.8, 94.3, 180]]),
        '[0.9129, 0.3651, 0.1826]': np.array([[15, 17, 23, 31],
                                [180, 139.9, 107.7, 159.3]]),
        '[0.6963, 0.6963, 0.1741]': np.array([[17, 21, 29],
                            [160.3, 124.9, 97.9]]),
        '[0.8704, 0.3482, 0.3482]': np.array([[17, 21, 29],
                                    [160.3, 124.9, 97.9]]),
        '[0.686, 0.5145, 0.5145]': np.array([[17, 19, 25],
                                        [180, 142.1, 111.1]]),
        '[0.8575, 0.5145, 0.0]': np.array([[17, 19, 25],
                                        [180, 142.1, 111.1]]),
        '[0.9864, 0.1644, 0.0]': np.array([[19, 23, 31],
                                [161.3, 127.5, 101.2]]),
        '[0.8111, 0.4867, 0.3244]': np.array([[19, 21, 27],
                                [180, 144.1, 114.0]]),
        '[0.9733, 0.1622, 0.1622]': np.array([[19, 21, 27],
                                [180, 144.1, 114.0]]),
        '[0.6247, 0.6247, 0.4685]': np.array([[21, 25],
                                [162.3, 129.8]]),
        '[0.7809, 0.6247, 0.0]': np.array([[21, 25],
                                [162.3, 129.8]]),
        '[0.937, 0.3123, 0.1562]': np.array([[21, 25],
                                [162.3, 129.8]]),
        '[0.8452, 0.5071, 0.169]': np.array([[9, 11, 15, 21, 29],
                                [160.8, 126.2, 99.6, 80.4, 66.6]]),
        '[0.7625, 0.4575, 0.4575]': np.array([[11, 13, 17, 23, 31],
                                [162.7, 130.8, 105.3, 86.3, 72.2]]),
        '[0.7001, 0.7001, 0.14]': np.array([[13, 15, 19, 25],
                                [164.1, 134.4, 110.0, 91.2]]),
        '[0.7715, 0.6172, 0.1543]': np.array([[21, 23, 29],
                                [180, 145.7, 116.6]]),
        '[0.7454, 0.5963, 0.2981]': np.array([[23, 27],
                                [163, 131.8]]),
        '[0.8847, 0.4423, 0.1474]': np.array([[23, 25, 31],
                                    [180, 147.1, 118.9]]),
        '[0.8571, 0.4286, 0.2857]': np.array([[25, 29],
                                    [163.7, 133.6]]),
        '[0.7071, 0.5657, 0.4243]': np.array([[25, 27],
                                    [180, 148.4]]),
        '[0.9899, 0.1414, 0.0]': np.array([[25, 27],
                                    [180, 148.4]]),
        '[0.9802, 0.14, 0.14]': np.array([[13, 15, 19, 25],
                                    [164.1, 134.4, 110.0, 91.2]]),
        '[0.6509, 0.6509, 0.3906]': np.array([[15, 17, 21, 27],
                                    [165.2, 137.3, 113.9, 95.3]]),
        '[0.9113, 0.3906, 0.1302]': np.array([[15, 17, 21, 27],
                                    [165.2, 137.3, 113.9, 95.3]]),
        '[0.8242, 0.5494, 0.1374]': np.array([[27, 31],
                                [164.4, 135.2]]),
        '[0.9615, 0.2747, 0.0]': np.array([[27, 31],
                                [164.4, 135.2]]),
        '[0.6804, 0.6804, 0.2722]': np.array([[27, 29], 
                                    [180, 149.6]]),
        '[0.9526, 0.2722, 0.1361]': np.array([[27, 29],
                                [180, 149.6]]),
        '[0.6623, 0.5298, 0.5298]': np.array([[29],
                                        [164.9]]),
        '[0.9272, 0.2649, 0.2649]': np.array([[29],
                                        [164.9]]),
        '[0.9191, 0.3939, 0.0]': np.array([[29, 31],
                                        [180, 150.6]]),
        '[0.8552, 0.3665, 0.3665]': np.array([[17, 19, 23, 29], 
                                    [166.1, 139.7, 117.2, 98.9]]),
        '[0.8083, 0.5774, 0.1155]': np.array([[19, 21, 25, 31],
                                            [166.8, 141.8, 120, 102.1]]),
        '[0.7683, 0.5488, 0.3293]': np.array([[21, 23, 27],
                                        [167.5, 143.6, 122.5]]),
        '[0.9879, 0.1098, 0.1098]': np.array([[21, 23, 27],
                                    [167.5, 143.6, 122.5]]),
        '[0.9435, 0.3145, 0.1048]': np.array([[23, 25, 29], 
                                        [168.0, 145.1, 124.7]]),
        '[0.7035, 0.5025, 0.5025]': np.array([[25, 27, 31], 
                                            [168.5, 146.4, 126.6]]),
        '[0.7035, 0.7035, 0.1005]': np.array([[25, 27, 31],
                                    [168.5, 146.4, 126.6]]),
        '[0.6767, 0.6767, 0.29]': np.array([[27, 29],
                                        [169.0, 147.7]]),
        '[0.6838, 0.4558, 0.5698]': np.array([[31],
                                    [165.4]]),
        '[0.7682, 0.6402, 0.0]': np.array([[31],
                                    [165.4]]),
        '[0.8701, 0.4834, 0.0967]': np.array([[27, 29],
                                        [169.0, 147.7]]),
        '[0.8393, 0.4663, 0.2798]': np.array([[29, 31],
                                        [169.4, 148.7]]),
        '[0.6312, 0.6312, 0.4508]': np.array([[31],
                                    [169.7]]),
        '[0.889, 0.381, 0.254]': np.array([[31], 
                                            [180]]),
        '[0.9918, 0.0902, 0.0902]': np.array([[31],
                                    [169.7]]),
        '[0.762, 0.635, 0.127]': np.array([[31], 
                                            [180]]),
        }

"""
https://stackoverflow.com/questions/33658620/generating-two-orthogonal-vectors-that-are-orthogonal-to-a-particular-direction
"""
def sharp_texture_quaternions(number_of_grains, required_texture, log_level):
    """
    Compute sharp textured quaternions in Scalar-first format.

    Processing
    ----------
    The function generates 2 orthogonal vectors based on the required texture 
    at random and then generates a rotation matrix by the combination of these 
    vectors. These rotation matrices are then transformed to unit quaternions.

    Parameters
    ----------
    number_of_grains: integer
        Number of grains.

    required_texture: 1D array of length 3
        Direction axis along which texture is to be generated.

    log_level: string
        Logger level to be used.

    Returns
    -------
    The function returns an array of sharp textured quaternions in scalar first
    format for the required number of grains.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Generating sharp texture quaternions')

    ## Normalizing the required texture direction
    vector_n = required_texture
    vector_n = vector_n * (1/np.linalg.norm(vector_n))
    
    all_orientation_of_grains = np.zeros([number_of_grains, 4])                 # Initializing an array to store all orientations
    
    ############################################################################
    # Source: https://www.ucl.ac.uk/~ucahmdl/LessonPlans/Lesson10.pdf
    # Based on Gram-Schmidt Process
    ############################################################################

    for i in range(number_of_grains):       
        vector_b = np.random.randn(3)                                           # Generating a vector at random
        vector_b -= (np.dot(vector_b, vector_n) * vector_n) / np.linalg.norm(vector_n)**2   # Making orthogonal to required texture
        vector_b /= np.linalg.norm(vector_b)                                    # Normalizing the vector

        vector_t = np.cross(vector_n, vector_b)                                 # Determining the third vector

        transpose_rotation_matrix = np.vstack((vector_b, vector_t, vector_n))             # Defininng the rotation matrix
        rotation_matrix = np.transpose(transpose_rotation_matrix)
        
        ## Testing for orthogonality of each row pairs
        assert np.isclose(np.dot(vector_b, vector_t), 0)
        assert np.isclose(np.dot(vector_t, vector_n), 0)
        assert np.isclose(np.dot(vector_b, vector_n), 0)

        ## Testing for orthogonality of each column pairs
        assert np.isclose(np.dot(rotation_matrix[:, 0], rotation_matrix[:, 1]), 0)
        assert np.isclose(np.dot(rotation_matrix[:, 1], rotation_matrix[:, 2]), 0)
        assert np.isclose(np.dot(rotation_matrix[:, 0], rotation_matrix[:, 2]), 0)

        ## Testing for orthogonality of the rotation matrix
        orthogonal_test = np.around((rotation_matrix @ np.transpose(rotation_matrix)), decimals=2)
        assert np.all(np.array_equal(orthogonal_test, np.eye(3)))
        
        ## Testing for the determinant of rotation matrix
        assert np.isclose(np.linalg.det(rotation_matrix), 1)
        
        ## Testing for reverse transformation from Rotation matrix - quaternion - rotation matrix
        roundtrip = quaternion.as_rotation_matrix(quaternion.from_rotation_matrix(rotation_matrix, nonorthogonal=False))
        assert np.all(np.allclose(roundtrip, rotation_matrix))
        
        ## Transforming rotation matrix to quaternion
        orientation_grain = quaternion.from_rotation_matrix(rotation_matrix)
        
        ## Assigning grain with the quaternion
        all_orientation_of_grains[i, :] = quaternion.as_float_array(orientation_grain)

    log.info('Completed generating sharp texture quaternions')    
    return all_orientation_of_grains

def random_quaternions_generator(number_of_grains, log_level):
    """
    Generate random Quaternions in Scalar-first format.

    Processing
    ----------
    The function generates random quaternions using Scipy library. Scipy uses 
    scalar last format, hence it is necessary to convert to scalar first format.
    
    Parameters
    ----------
    number_of_grains: integer
        Number of grains.

    required_texture: 1D array of length 3
        Direction axis along which texture is to be generated.

    log_level: string
        Logger level to be used.

    Returns
    -------
    The function returns an array of random quaternions in scalar first
    format for the required number of grains.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Generating random orientation quaternions')
    ## REMEMBER SCIPY USES SCALAR LAST FORMAT AND NUMPY-QUATERNION USES SCALAR FIRST FORMAT
    ## Random Quaternions in Scalar-last format
    random_quats_array = R.random(number_of_grains, random_state = None).as_quat()
    
    ## Converting to Scalar-first Format
    random_quats_reperesentation = np.zeros_like(random_quats_array)
    random_quats_reperesentation[:, 0] = random_quats_array[:, 3]
    random_quats_reperesentation[:, 1:4] = random_quats_array[:, 0:3]
    
    log.info('Completed generating random orientation quaternions')
    return random_quats_reperesentation

def disorientation_angles(dimension, limit, skewed_boundary_flag, required_texture, rand_quat_flag, orientation_data, tessellation_og, log_level):
    """
    Compute disorientation angles between grains.

    Processing
    ----------
    The misorientation angle are computed between two grains along with its 
    symmetric orientations. The minimum misorientation angle is considered to 
    be the disorientation angle and its corresponding axis as its 
    disorientation axis.

    Parameters
    ----------
    dimension: integer    
        Dimension of study (2 or 3)

    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

    skewed_boundary_flag: boolean 
        Flag to specify if skewed grain boundaries are required. Only functional
        in quasi-2D case.
    
    required_texture: 1D array of length 3
        Direction axis along which texture is to be generated.

    rand_quat_flag: boolean
        Flag to indicate that random orientations are to be assigned to all 
        grains.

    orientation_data: array of shape (number of grains, 4)
        Orientation data of each grain as rows. Orientations to be specified as
        Quaternions with Scalar-first format.

    tessellation_og: dictionary
        Dictionary of tessellations data with following keys:
            1. number_of_grains
            2. number_of_faces_list
            3. vertices_list
            4. face_vertices_list
            5. centroid_list
            6. volume_list
            7. normals_list
            8. neighbors_list
            9. face_area_list
            10. number_of_edges_list

    log_level: string
        Logger level to be used.

    Returns
    -------
        1. An array consisting of columns grain 1, grain 2, disorientation 
            angle, disorientation axis, GB area.
        2. An array of generated or inherited orientations data. This is returned to
            ensure that when random or sharp texture orientations are generated, they 
            are generated only once and same are used in other functions such as 
            disorientation_angles, type_of_grain_boundary, etc.  
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Started computing disorientation angles')
    tessellation = copy.deepcopy(tessellation_og)

    parent_function_name = inspect.stack()[1][3]
    grain_boundary_area, all_vertices_list = grain_boundary_areas(dimension, limit, tessellation, parent_function_name, skewed_boundary_flag, log_level)

    ## Again creating a deepcopy to avoid call by reference issues
    tessellation = copy.deepcopy(tessellation_og)

    number_of_grains = copy.deepcopy(tessellation['number_of_grains'])                                        # extracting the total number of grains    
    quaternions_of_grains = np.zeros([number_of_grains, 5])                     # Array with column names: Gr. no., quaternions
    quaternions_of_grains[:, 0] = range(number_of_grains)                       # assigning gr. no. to first column
    
    if orientation_data is None:
        if rand_quat_flag:
            quaternions_of_grains[:, 1:5] = random_quaternions_generator(number_of_grains, log_level)
        else:
            quaternions_of_grains[:, 1:5] = sharp_texture_quaternions(number_of_grains, required_texture, log_level)    # Assigning random quaternions to each grain
    else: 
        quaternions_of_grains[:, 1:5] = orientation_data

    neighbors_each_grain = copy.deepcopy(tessellation['neighbors_list']) #[v.neighbors() for v in tessellation]                # extracting the neighbors of each grain

    disorientation_data = []                                                    ## column names: grain 1, grain 2, disorientation angle, disorientation axis, GB area
    for grain in range(number_of_grains):
        
        for neighbor in neighbors_each_grain[grain]:                            # iterating through each neighbor of respective grain
            ## Checking if CSL type data already exists for current grain-neighbor pair
            if len(disorientation_data) != 0:
                match_index = np.where((np.array(disorientation_data)[:, 1] == grain) & (np.array(disorientation_data)[:, 0] == neighbor))[0]
                if len(match_index) != 0:
                    #print(grain, neighbor, np.array(csl_data)[match_index])
                    continue
            
            ####### Second condition to accomodate for one seed test ###########
            if (grain == neighbor) and (np.count_nonzero(neighbors_each_grain) != 0):  # Due to periodicity in 2D case, self is also considered as neighbor 
                continue
            disorientation_data_single_neighbor = []                            # Disorientation data like gr. no., neighbor gr., disorientation angle and axis will be stored
            disorientation_data_single_neighbor.append(grain)
            disorientation_data_single_neighbor.append(neighbor)
            
            first_grain_quaternion = quaternion.as_quat_array(quaternions_of_grains[grain, 1:5])        # extracting quaternion of the respective grain
            second_grain_quaternion = quaternion.as_quat_array(quaternions_of_grains[neighbor, 1:5])    # extracting qauternion of neighbour grain
            inverse_second_grain_quaternion = second_grain_quaternion.inverse()                         # doing inverse of the quternion of neighbour grain
            
            misorientation = first_grain_quaternion * inverse_second_grain_quaternion

            ####################################################################
            # if using the below misorientation quaternion then the disorientation
            # is achieved as 48.1897 degrees with 0.8944, 0, 0.4472 as disorientation
            # angle and axis respectively. This is in agreement with the result shown 
            # in the Book 'Measurement of Grain Boundary Geometry' by V Randle
            # Page no. 44. The same has been done by converting axis angle to
            # quaternion for Σ = 5 and Σ = 3.
            # ##################################################################   
            
            # misorientation = quaternion.as_quat_array([0.7745967, 0.5163978, 0.2581989, 0.2581989])         # Σ = 15         
            # misorientation = quaternion.as_quat_array([0, 0.9486773, 0.3162458, 0])                         # Σ = 5                     
            # misorientation = quaternion.as_quat_array([0, 0.8164966, 0.4082483, 0.4082483])                # Σ = 3

            ## initializing lists where all the misorientation angles and axes would be stored for finding minimum angle
            all_misorientation_angles =[]
            misorientation_axis = []

            ## Iterating over each symmetric operator
            for operator in symmetric_operators_quaternions:
                
                ## Product of symmetry operator and misorientation
                operator_quaternion = operator
                product_operator_misorientation = operator_quaternion * misorientation

                product_operator_misorientation_array = quaternion.as_float_array(product_operator_misorientation)
                trace_product_operator_misorientation = 4 * (product_operator_misorientation_array[0])**2 - 1
                
                ## Determining Misorientation angle               
                rotation_angle = np.degrees(np.arccos((trace_product_operator_misorientation-1)/2))
                
                ## Determining Misorientation axis
                denominator = (product_operator_misorientation_array[1]**2 + product_operator_misorientation_array[2]**2 + product_operator_misorientation_array[3]**2)**0.5 
                rotation_axis = [product_operator_misorientation_array[1]/denominator, product_operator_misorientation_array[2]/denominator, product_operator_misorientation_array[3]/denominator]
              
                ## Appending the misorientation angle and axis to storage lists
                all_misorientation_angles.append(rotation_angle)
                misorientation_axis.append(rotation_axis)
           
            ## Determing the Minimum misorientation angle and its corresponding axis (Disorientation angle and axis)
            index_minimum_misorientation_angle = np.argmin(all_misorientation_angles)
            disorientation_angle = np.around(all_misorientation_angles[index_minimum_misorientation_angle], decimals=4)
            disorientation_axis = misorientation_axis[index_minimum_misorientation_angle]
            required_texture_test = required_texture * (1/np.linalg.norm(required_texture))
            
            if rand_quat_flag is False:
                assert np.all(np.allclose(np.abs(required_texture_test), np.abs(disorientation_axis)))
            
            ## Finding GB area for this current pair of grain and neighbor            
            ## index of grain-neighbor pair in grain_boundary_area
            ## It is to be noted that the sequence of GB area is exactly the same as 
            ## that of Type of GB. Hence it can also be horizontally stacked instead
            ## of finding matching pair. Matching is done to be on safer side.
            index_pair_gb_area = np.where((grain_boundary_area[:, 1] == grain) & (grain_boundary_area[:, 2] == neighbor))[0]
            current_pair_gb_area = grain_boundary_area[index_pair_gb_area, 3]   # 4th column is the area

            ## Storing the corresponding disorientation angle and axis to main array
            disorientation_data_single_neighbor.append(disorientation_angle)
            disorientation_data_single_neighbor = disorientation_data_single_neighbor + disorientation_axis
            disorientation_data_single_neighbor.append(current_pair_gb_area[0])
            disorientation_data.append(disorientation_data_single_neighbor)

    log.info('Completed computing disorientation angles')
    return np.array(disorientation_data), np.array(quaternions_of_grains[:, 1:5])

def type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation_og, dimension, limit, skewed_boundary_flag, log_level):
    """
    Compute type of grain boundaries between each pair of grains.

    Processing
    ----------
    The function computes misorientations between grains and its symmetries. The 
    function then comparess the misorientation angle and axis and based on the 
    Brandons Criterion decides if the Grain boundary is a specific Special Grain
    boundary or not and then stores it into the main list.

    Parameters
    ----------
    required_texture: 1D array of length 3
        Direction axis along which texture is to be generated.

    rand_quat_flag: boolean
        Flag to indicate that random orientations are to be assigned to all 
        grains.

    orientation_data: array of shape (number of grains, 4)
        Orientation data of each grain as rows. Orientations to be specified as
        Quaternions with Scalar-first format.

    tessellation_og: dictionary
        Dictionary of tessellations data with following keys:
            1. number_of_grains
            2. number_of_faces_list
            3. vertices_list
            4. face_vertices_list
            5. centroid_list
            6. volume_list
            7. normals_list
            8. neighbors_list
            9. face_area_list
            10. number_of_edges_list

    dimension: integer    
        Dimension of study (2 or 3)

    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

    skewed_boundary_flag: boolean 
        Flag to specify if skewed grain boundaries are required. Only functional
        in quasi-2D case.

    log_level: string
        Logger level to be used.

    Returns
    -------
        1. An array with the column names as grain 1, grain 2, rotation angle, 
            rotation axis, type of csl (0 indicates normal grain boundary), GB
            area.
        2. An array of generated or inherited orientations data. This is returned to
            ensure that when random or sharp texture orientations are generated, they 
            are generated only once and same are used in other functions such as 
            disorientation_angles, type_of_grain_boundary, etc. 
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Started computing type of grain boundaries')
    tessellation = copy.deepcopy(tessellation_og)

    parent_function_name = inspect.stack()[1][3]
    grain_boundary_area, all_vertices_list = grain_boundary_areas(dimension, limit, tessellation, parent_function_name, skewed_boundary_flag, log_level)

    ## Again creating a deepcopy to avoid call by reference issues
    tessellation = copy.deepcopy(tessellation_og)

    number_of_grains = copy.deepcopy(tessellation['number_of_grains'])                                          # extracting the total number of grains    
    quaternions_of_grains = np.zeros([number_of_grains, 5])                     # Array with column names: Gr. no., quaternions
    quaternions_of_grains[:, 0] = range(number_of_grains)                       # assigning gr. no. to first column 

    ## Checking if orientation data is provided and if not then generating random
    if orientation_data is None:
        ## Checking if random orientations are required or sharp texture is required
        if rand_quat_flag:
            quaternions_of_grains[:, 1:5] = random_quaternions_generator(number_of_grains, log_level)
        else:
            quaternions_of_grains[:, 1:5] = sharp_texture_quaternions(number_of_grains, required_texture, log_level)    # Assigning random quaternions to each grain
    else: 
        quaternions_of_grains[:, 1:5] = orientation_data
    
    neighbors_each_grain = copy.deepcopy(tessellation['neighbors_list']) #[v.neighbors() for v in tessellation]                # extracting the neighbors of each grain

    csl_data = []                                                    ## column names: grain 1, grain 2, rotation angle, rotation axis, type of csl, GB area
    for grain in range(number_of_grains):
        
        for neighbor in neighbors_each_grain[grain]:                            # iterating through each neighbor of respective grain
            ## Checking if CSL type data already exists for current grain-neighbor pair
            if len(csl_data) != 0:
                match_index = np.where((np.array(csl_data)[:, 1] == grain) & (np.array(csl_data)[:, 0] == neighbor))[0]
                if len(match_index) != 0:
                    #print(grain, neighbor, np.array(csl_data)[match_index])
                    continue

            ####### Second condition to accomodate for one seed test ###########
            if (grain == neighbor) and (np.count_nonzero(neighbors_each_grain) != 0):  # Due to periodicity in 2D case, self is also considered as neighbor 
                continue
            misorientation_data_single_neighbor = []                            # misorientation data pertaining to single combination of grain and its neighbor would be stored
            misorientation_data_single_neighbor.append(grain)
            misorientation_data_single_neighbor.append(neighbor)
            
            first_grain_quaternion = quaternion.as_quat_array(quaternions_of_grains[grain, 1:5])        # extracting quaternion of the respective grain
            second_grain_quaternion = quaternion.as_quat_array(quaternions_of_grains[neighbor, 1:5])    # extracting qauternion of neighbour grain
            inverse_second_grain_quaternion = second_grain_quaternion.inverse()                         # doing inverse of the quternion of neighbour grain
            
            ## Finding Misorientation between grain and its neighbor
            misorientation = first_grain_quaternion * inverse_second_grain_quaternion   

            ####################################################################
            # The quaternion representing the misorientation of Σ = 15 is 
            # (3, 2, 1, 1) in the real part first format. Using this quaternion
            # as the misorientation between the grains, the resulting type of 
            # grain boundary  should be Σ = 15. This has been observed by using 
            # the misorientation quaternion commented in the below statement and
            # the type of GB were all Σ = 15. This is in agreement with the result shown 
            # in the Book 'Measurement of Grain Boundary Geometry' by V Randle
            # Page no. 44. The same has been done by converting axis angle to
            # quaternion for Σ = 5 and Σ = 3.
            # ##################################################################

            # misorientation = quaternion.as_quat_array([0.7745967, 0.5163978, 0.2581989, 0.2581989])         # Σ = 15         
            # misorientation = quaternion.as_quat_array([0, 0.9486773, 0.3162458, 0])                         # Σ = 5                     
            # misorientation = quaternion.as_quat_array([0, 0.8164966, 0.4082483, 0.4082483])                # Σ = 3
            
            ## Broadcasting misorientation to the shape (-1, 1) and symmetric operators from (24,) to (24, 1)
            broadcast_misorientation = np.tile(misorientation, (symmetric_operators_quaternions.shape[0], 1))
            broadcast_symmetric_operator_quaternions = symmetric_operators_quaternions.reshape((-1, 1))
            
            ## Multiplying the misorientation with individual symmetric operators
            product_operator_misorientation = np.array([x * y for x, y in zip(broadcast_symmetric_operator_quaternions, broadcast_misorientation)])
            
            ## Converting the quaternion to an array
            product_operator_misorientation_array = quaternion.as_float_array(product_operator_misorientation)
            
            ## While converting from quat to array it was a 3d array, hence converting to 2d array
            product_operator_misorientation_array = np.array([x[0] for x in product_operator_misorientation_array])
            
            ## Finding trace of the quaternion using formula 4*q_0**2 - 1
            trace_product_operator_misorientation = 4 * (product_operator_misorientation_array[:, 0])**2 - 1
            
            ## Determining Misorientation angle               
            rotation_angle = np.degrees(np.arccos((trace_product_operator_misorientation-1)/2))

            ## Determining Misorientation axis
            denominator = ((product_operator_misorientation_array[:, 1]**2 + product_operator_misorientation_array[:, 2]**2 + product_operator_misorientation_array[:, 3]**2)**0.5).reshape((-1, 1))
            rotation_axis_column_1 = (product_operator_misorientation_array[:, 1]).reshape(-1, 1)/denominator[:]
            rotation_axis_column_2 = (product_operator_misorientation_array[:, 2]).reshape(-1, 1)/denominator[:]
            rotation_axis_column_3 = (product_operator_misorientation_array[:, 3]).reshape(-1, 1)/denominator[:]
            rotation_axis = np.hstack((rotation_axis_column_1, rotation_axis_column_2, rotation_axis_column_3))

            ## Initializing the final Rotation axis and angle            
            smallest_rotation_axis = rotation_axis[np.argmin(rotation_angle), :]
            smallest_rotation_angle = rotation_angle[np.argmin(rotation_angle)]
            
            ## Initializing the Smallest dot product angle discovered and the type of CSL
            smallest_dot_product_angle_checker = 360
            type_of_csl = 0
            
            ## Iterating through each of the symmetric rotation axis
            for row_number, row in enumerate(rotation_axis):
                
                ## Broadcasting it to the number of rows of the available rotation axis array
                broadcast_each_rotation_axis = np.tile(row, (len(norm_available_required_texture), 1))
                
                ## Finding dot product of each of the available rotation axis with obtained rotation axis and converting to angles
                ################################################################
                # Here rounding off is very important else dot product returns nan
                ################################################################
                dot_product_rows = np.around(np.array([np.dot(x, y) for x, y in zip(broadcast_each_rotation_axis, norm_available_required_texture)]), decimals=2)   ## Dot product of row with each available texture info row
                
                ################################################################
                # Here it is necessary to find absolute of the dot product and 
                # then convert it to degrees since the range of cosine is from
                # -1 to 1 ie; from 180° to 0° respectively and our condition is
                # based on dot product to be 0° which is not trapped if the 
                # dot product is -1. It is necessary to assume that the rotation
                # axis is the same even if the rotation axes are at 180° to each
                # other.
                ################################################################
                angle_dot_product = np.degrees(np.arccos(np.abs(dot_product_rows)))             # Converting to angle ##############################################################################################remove abs#######################################
                
                ## Finding the minimum of the dot product angles
                min_angle_dot_product = np.min(angle_dot_product)
                
                ## Checking if the min angle is less than 2 degrees and also 
                # if the type of csl is already discovered and if the min angle 
                # is less than the smallest dot product angle discovered till now
                if min_angle_dot_product < 2 and ((min_angle_dot_product < smallest_dot_product_angle_checker) or (type_of_csl == 0)):           ############################ Arbitrarily choosing 2 degree to be the threshold
                    index_min_angle = np.argmin(angle_dot_product)                      # Finding index of lowest angle
                    corresponding_rotation_angle = rotation_angle[row_number]
                    matching_rotation_axis = np.around(norm_available_required_texture[index_min_angle], decimals= 4)

                    ## Extracting data from the dictionary based on key
                    key_direction = "".join(str(list(matching_rotation_axis)))
                    type_of_csl_array = type_of_csl_data[key_direction]
                    type_of_csl_angles = type_of_csl_array[1, :]

                    # Finding the index of closest angle
                    index_of_closest_angle = (np.abs(type_of_csl_angles - corresponding_rotation_angle)).argmin()
                    
                    ## Using Brandons Criterion to decide if the GB is of special type
                    brandon_criteria_angle_min = type_of_csl_array[1, index_of_closest_angle] - (15 * (1/(type_of_csl_array[0, index_of_closest_angle])**0.5))      # 15 is constant and formula is [constant * (sigma_value)**(-0.5)]
                    brandon_criteria_angle_max = type_of_csl_array[1, index_of_closest_angle] + (15 * (1/(type_of_csl_array[0, index_of_closest_angle])**0.5))
                    
                    ## Extracting the type of CSL if within bounds
                    if ((corresponding_rotation_angle > brandon_criteria_angle_min) & (corresponding_rotation_angle < brandon_criteria_angle_max)):
                        type_of_csl = type_of_csl_array[0, index_of_closest_angle]
                    
                    ## Defining final rotation axis and angle
                    smallest_rotation_axis = row
                    smallest_rotation_angle = corresponding_rotation_angle

            ## Finding GB area for this current pair of grain and neighbor            
            ## index of grain-neighbor pair in grain_boundary_area
            ## It is to be noted that the sequence of GB area is exactly the same as 
            ## that of Type of GB. Hence it can also be horizontally stacked instead
            ## of finding matching pair. Matching is done to be on safer side.
            index_pair_gb_area = np.where((grain_boundary_area[:, 1] == grain) & (grain_boundary_area[:, 2] == neighbor))[0]
            current_pair_gb_area = grain_boundary_area[index_pair_gb_area, 3]   # 4th column is the area

            ## Updating the data to main array
            misorientation_data_single_neighbor.append(smallest_rotation_angle)
            misorientation_data_single_neighbor += list(smallest_rotation_axis)
            misorientation_data_single_neighbor.append(type_of_csl)
            misorientation_data_single_neighbor.append(current_pair_gb_area[0]) # slicing or else appended as 1 element array object
            csl_data.append(misorientation_data_single_neighbor)
    
    log.info('Completed computing type of grain boundaries')
    return np.array(csl_data), np.array(quaternions_of_grains[:, 1:5])
                

#######################################################################################
# normalize each slip plane and direction or else it gives a huge difference in results
#######################################################################################
def slip_system_generator(slip_system_family, crystal_symmetry_type):
    """
    Generates all slip systems

    Processing
    ----------
    The function checks all permutations of symmetries of slip systems for its
    orthogonality. It also checks if the same parallel slip system exists in the
    list.

    Parameters
    ----------
    slip_system_family: 1D array of shape (6,)
        First three components denote family of slip plane normal and remaining
        three components denote family of slip direction.

    crystal_symmetry_type: string
        Crystal symmetry to be used.(CUBIC, ORTHORHOMBIC, HEXAGONAL, TETRAGONAL)

    Returns
    -------
        An array of all slip systems. First three columns denote slip plane 
        normal and remaining three columns denote slip direction.
    """
    slip_normal_family = slip_system_family[:3]
    slip_direction_family = slip_system_family[3:]
    
    symmetry_operator_types = {'CUBIC': symmetry_operator, 
                               'ORTHORHOMBIC': orthorhombic_symmetry_operators,
                               'HEXAGONAL': hexagonal_symmetry_operators,
                               'TETRAGONAL': tetragonal_symmetry_operators}
    crystal_symmetry_operator = symmetry_operator_types[crystal_symmetry_type.upper()]

    all_slip_system = [[0, 0, 0, 0, 0, 0]]                                      ## adding a row so that first time checking of rows is made possible. If this array is empty then checking is not possible.
    ## Creating permutations of symmetries of slip plane normal and slip direction
    for index, operator in enumerate(crystal_symmetry_operator):
        rot_norm = operator@slip_normal_family
        for index_dir, operator_dir in enumerate(crystal_symmetry_operator):
            rot_dir = operator_dir@slip_direction_family
            dot = np.dot(rot_norm, rot_dir)
            
            ## Slip system must be orthogonal in this case
            if np.isclose(dot, 0.0):
                check_list = []
                
                ## The permutations of negatives of symmetric normal and directions are checked if they already exist in the slip system list
                for i, j in zip([1, 1, -1, -1], [1, -1, 1, -1]):
                    
                    ## Creating a row having slip plane normal and direction
                    match = np.hstack([i * rot_norm, j * rot_dir])
                    
                    ## Checking if row exists
                    check = any(np.equal(np.array(all_slip_system), match).all(1))
                    
                    ## variable 'Check' boolean result is added to a list
                    check_list.append(check)
                    
                    ## If a match is found then the loop is broken immediately
                    if check == True:
                        break
                
                ## Adding the row to all slip system list only if no match is found for any permutation
                if np.all(~(np.array(check_list, dtype=bool))):
                    all_slip_system.append(np.hstack([rot_norm, rot_dir]))

    ## deleting first dummy zeros row
    all_slip_system = np.delete(np.array(all_slip_system), 0, 0)
    return all_slip_system

def schmid_factor(required_texture, rand_quat_flag, dimension, limit, stress_direction, slip_system_family, crystal_symmetry_type, orientation_data, tessellation_og, log_level):
    """
    Compute Schmid factors.

    Processing
    ----------
    The function computes Schmid factor for each grain for each slip system and 
    then stores the maximum of schmid factor for each grain as the slip system 
    having maximum schmid factor starts to deform first.
    
    Parameters
    ----------
    required_texture: 1D array of length 3
        Direction axis along which texture is to be generated.

    rand_quat_flag: boolean
        Flag to indicate that random orientations are to be assigned to all 
        grains.

    dimension: integer    
        Dimension of study (2 or 3)

    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

    stress_direction: array of length 3    
        Direction of loading

    slip_system_family: 1D array of shape (6,)
        First three components denote family of slip plane normal and remaining
        three components denote family of slip direction.

    crystal_symmetry_type: string
        Crystal symmetry to be used.(CUBIC, ORTHORHOMBIC, HEXAGONAL, TETRAGONAL)

    orientation_data: array of shape (number of grains, 4)
        Orientation data of each grain as rows. Orientations to be specified as
        Quaternions with Scalar-first format.

    tessellation_og: dictionary
        Dictionary of tessellations data with following keys:
            1. number_of_grains
            2. number_of_faces_list
            3. vertices_list
            4. face_vertices_list
            5. centroid_list
            6. volume_list
            7. normals_list
            8. neighbors_list
            9. face_area_list
            10. number_of_edges_list

    log_level: string
        Logger level to be used.

    Returns
    -------
        1. An array consisting of Schmid Factor related data with column names as 
            Grain no., Schmid factor, Slip plane and direction (Slip System), 
            2nd max Schmid factor, Slip system and 3rd max Schmid factor, Slip 
            system.
        2. An array of generated or inherited orientations data. This is returned to
            ensure that when random or sharp texture orientations are generated, they 
            are generated only once and same are used in other functions such as 
            disorientation_angles, type_of_grain_boundary, etc.
    """

    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Started computing Schmid factors')
    tessellation = copy.deepcopy(tessellation_og)

    ## Generating slip systems
    slip_systems = slip_system_generator(slip_system_family, crystal_symmetry_type)

    ## Computing grain size distribution
    grain_sizes_data = grain_size_distribution(dimension, tessellation, limit, log_level)
    grain_sizes = grain_sizes_data[:, 1]

    ## copying tessellation to avoid errors due to call by reference
    tessellation = copy.deepcopy(tessellation_og)

    ## Checking the function input arguments
    if dimension == 2 and stress_direction[2] != 0:
        print("Stress cannot be along Z-direction for 2-D case")
        exit()
    
    stress_direction_index = [np.where(stress_direction == 0)]

    ## Testing if direction [0, 0, 0] is provided. Second condition to check if proper dimension of vector is provided
    if (len(stress_direction_index[0][0]) > 2) or (len(stress_direction) > 3):
        print("Stress Vector should have definite direction")
        exit()

    #stress_direction_value = stress_direction_index[0][0][0]
    
    stress_direction = (1/np.linalg.norm(stress_direction)) * stress_direction

    number_of_grains = copy.deepcopy(tessellation['number_of_grains'])                                         # extracting the total number of grains    
    quaternions_of_grains = np.zeros([number_of_grains, 5])                     # Array with column names: Gr. no., quaternions
    quaternions_of_grains[:, 0] = range(number_of_grains)                       # assigning gr. no. to first column
    
    if orientation_data is None:
        if rand_quat_flag:
            quaternions_of_grains[:, 1:5] = random_quaternions_generator(number_of_grains, log_level)
        else:
            quaternions_of_grains[:, 1:5] = sharp_texture_quaternions(number_of_grains, required_texture, log_level)    # Assigning random quaternions to each grain
    else: 
        quaternions_of_grains[:, 1:5] = orientation_data
    
    ## Storing inverse of quaternions
    inverse_quaternions_of_grains = np.zeros([number_of_grains, 5])                     
    inverse_quaternions_of_grains[:, 0] = range(number_of_grains)
    for i in range(number_of_grains):
        quaternion_to_be_inversed = quaternion.as_quat_array(quaternions_of_grains[i, 1:5])
        inverse_quaternions_of_grains[i, 1:5] = quaternion.as_float_array(quaternion_to_be_inversed.inverse())
    
    ## Initializing a list where all data related to Schmid Factor would be stored
    schmid_factors = []                                                         # Column names are: grain, grain size, 1st max schmid factor, respective slip plane and slip direction, 2nd max schmid factor, respective slip plane and slip direction, 3rd max schmid factor, respective slip plane and slip direction
    single_gr_single_sys_schmid_factor = []                                     # Column names are: grain, grain size, Schmid factor, slip system

    for grain in range(number_of_grains):
        
        ## Extracting quaternion of grain and its inverse
        quaternion_of_grain = quaternion.as_quat_array(quaternions_of_grains[grain, 1:5])
        inverse_quaternion_of_grain = quaternion.as_quat_array(inverse_quaternions_of_grains[grain, 1:5])

        all_schmid_factors_grain = []

        ## Calculating Schmid Factor for each slip system
        for slip_system in slip_systems:
            
            ## Normalizing Slip Plane
            slip_plane =  [0] + list(slip_system[0:3])
            slip_plane = (1/np.linalg.norm(slip_plane))*np.array(slip_plane)

            ## Normalizing Slip Direction
            slip_direction =  [0] + list(slip_system[3:6])
            slip_direction = (1/np.linalg.norm(slip_direction))*np.array(slip_direction)

            ## Orientaing Slip System to Crystal Coordinates (Passive Rotation)
            rotated_slip_plane = (inverse_quaternion_of_grain * quaternion.as_quat_array(slip_plane)) * quaternion_of_grain
            rotated_slip_direction = (inverse_quaternion_of_grain * quaternion.as_quat_array(slip_direction)) * quaternion_of_grain
            
            ## Converting the rotated Slip Sytem to an array
            rotated_slip_plane_array = quaternion.as_float_array(rotated_slip_plane)
            rotated_slip_direction_array = quaternion.as_float_array(rotated_slip_direction)

            ## Computing Schmid Factor based on the axis of applied uniaxial stress
            rotated_slip_plane_vector = rotated_slip_plane_array[1:]
            rotated_slip_direction_vector = rotated_slip_direction_array[1:]
            
            cosine_lambda = np.dot(stress_direction, rotated_slip_plane_vector)/(np.linalg.norm(stress_direction) * np.linalg.norm(rotated_slip_plane_vector))
            cosine_theta = np.dot(stress_direction, rotated_slip_direction_vector)/(np.linalg.norm(stress_direction) * np.linalg.norm(rotated_slip_direction_vector))
            
            schmid_factor_grain = cosine_theta * cosine_lambda

            all_schmid_factors_grain.append(schmid_factor_grain)                # used for main schmid factor array
            single_gr_single_sys_schmid_factor.append([grain] + [grain_sizes[grain]] \
                + [schmid_factor_grain] + list(slip_system))                    # for writing to a separate file other than textural_characteristics.txt
        

        all_schmid_factors_grain = np.abs(np.array(all_schmid_factors_grain))   # Finding absolute values of each schmid factor
        # maximum_schmid_factor_index = np.argmax(all_schmid_factors_grain)       # Finding index of Maximum schid factor
        # second_max_schmid_factor_index = 
        ## based on answer "https://stackoverflow.com/questions/6910641/how-do-i-get-indices-of-n-maximum-values-in-a-numpy-array"
        ## by Fred Foo

        no_sorted_indices_max_schmid_factor = np.argpartition(all_schmid_factors_grain, -3)[-3:]       # Indices of max 3 schmid factors
        sorted_indices_max_schmid_factor = no_sorted_indices_max_schmid_factor[np.argsort(all_schmid_factors_grain[no_sorted_indices_max_schmid_factor])]   # Indices are in ascending order of schmid factor values
        indices_schmid_descending = sorted_indices_max_schmid_factor[::-1]      # Indices are in descending order of schmid factor values
        
        
        ## Appending all the Schmid Factor related data to main list
        schmid_factor_data_individual_grain = [grain] + [grain_sizes[grain]] \
            + [all_schmid_factors_grain[indices_schmid_descending[0]]] + list(slip_systems[indices_schmid_descending[0]]) \
            + [all_schmid_factors_grain[indices_schmid_descending[1]]] + list(slip_systems[indices_schmid_descending[1]]) \
            + [all_schmid_factors_grain[indices_schmid_descending[2]]] + list(slip_systems[indices_schmid_descending[2]])
        schmid_factors.append(schmid_factor_data_individual_grain)
    
    log.info('Completed computing Schmid factors')

    return np.array(schmid_factors), np.array(single_gr_single_sys_schmid_factor), np.array(quaternions_of_grains[:, 1:5])






