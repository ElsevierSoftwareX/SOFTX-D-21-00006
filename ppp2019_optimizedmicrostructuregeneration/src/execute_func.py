# -*- coding: utf-8 -*-
"""
execute_func.py

Module to perform tasks that are common to tests as well as main program.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

@authors: Serrao Prince Henry, Arun Prakash
@email: 
created: 22 January 2020
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

from ppp2019_optimizedmicrostructuregeneration.src.main_import_statements import *

from ppp2019_optimizedmicrostructuregeneration.src.set_logger import set_logger as set_logger
name_str = __name__

from ppp2019_optimizedmicrostructuregeneration.src.__version__ import __version__ as version

from ppp2019_optimizedmicrostructuregeneration.src.create_tessellations import create_tessellations as create_tessellations

from ppp2019_optimizedmicrostructuregeneration.src.create_obj_files import create_obj_file_all_grains as create_obj_file_all_grains
from ppp2019_optimizedmicrostructuregeneration.src.create_obj_files import create_obj_file_individual_grains as create_obj_file_individual_grains
from ppp2019_optimizedmicrostructuregeneration.src.create_vtk_files import create_vtk_file_all_grains as create_vtk_file_all_grains
from ppp2019_optimizedmicrostructuregeneration.src.create_vtk_files import create_vtk_file_individual_grains as create_vtk_file_individual_grains

from ppp2019_optimizedmicrostructuregeneration.src.random_generator import random_generator as random_generator
from ppp2019_optimizedmicrostructuregeneration.src.structural_characteristics import grain_size_distribution as grain_size_distribution
from ppp2019_optimizedmicrostructuregeneration.src.structural_characteristics import number_of_neighbors as number_of_neighbors
from ppp2019_optimizedmicrostructuregeneration.src.structural_characteristics import grain_boundary_areas as grain_boundary_areas
from ppp2019_optimizedmicrostructuregeneration.src.structural_characteristics import junction_length as junction_length
from ppp2019_optimizedmicrostructuregeneration.src.structural_characteristics import junction_angle as junction_angle
from ppp2019_optimizedmicrostructuregeneration.src.structural_characteristics import distance_btw_grains as distance_btw_grains
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import sharp_texture_quaternions as sharp_texture_quaternions
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import random_quaternions_generator as random_quaternions_generator
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import disorientation_angles as disorientation_angles
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import type_of_grain_boundary as type_of_grain_boundary
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import schmid_factor as schmid_factor
from ppp2019_optimizedmicrostructuregeneration.src.textural_characteristics import available_required_texture as available_required_texture

from ppp2019_optimizedmicrostructuregeneration.src.mesh import mesh_hex as mesh_hex
from ppp2019_optimizedmicrostructuregeneration.src.mesh import mesh_tetra as mesh_tetra
from ppp2019_optimizedmicrostructuregeneration.src.mesh import mesh_visualization as mesh_visualization

def execute_func(size_of_simulation_box, dimension, limit, material, orientation_data, required_texture, rand_quat_flag, seed_array_unique, stress_direction, store_folder, face_flag, now, number_of_bins, skewed_boundary_flag, mesh_flag, global_mesh_size, log_level):
    """
    The execute_func is a function which computes data based on the input parameters.
    The function then stores the computed data in the form of text files and plots.
    The function can be used for tests as well as for the main program

    Parameters
    ----------

    size_of_simulation_box: array of length 3
        Size of simulation box (array of length along X, Y, Z directions)

    dimension: integer    
        Dimension of study (2 or 3)

    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

    material: string
        Material name

    orientation_data: array of shape (number of grains, 4)
        Orientation data of each grain as rows. Orientations to be specified as
        Quaternions with Scalar-first format.

    required_texture: array of length 3
        Specific texture to be used for all grains

    rand_quat_flag: boolean
        Flag to indicate if random quaternions are to be assigned

    seed_array_unique: array of shape (number of grains, 3)
        Unique seed coordinates
        
    stress_direction: array of length 3    
        Direction of loading

    store_folder: string
        Name of directory where output files are to be stored.
    
    face_flag: boolean
        Flag to indicate that opaque surface is to be used instead of transparent.
    
    now: string
        Current time and date.
    
    number_of_bins: integer
        Total number of bins to be used while computing distribution

    skewed_boundary_flag: boolean 
        Flag to specify if skewed grain boundaries are required. Only functional
        in quasi-2D case.

    mesh_flag: string (Hex, quat or vis)
        String to indicate type of mesh elements required.
        
    global_mesh_size: float   
        Glomal mesh size to be used

    log_level: string
        Logger level to be used.

    Returns
    ------- 
    The function returns following:
        1. grain size distribution
        2. Information of number of neighbors 
        3. Information of grain boundary area
        4. Information of junction lengths
        5. Information of junction angles
        6. Information of distance between grains
        7. Information of distance between grains as 1D array
        8. Information of disorientation angles
        9. Information of Schmid factors
        10. Information of type of grain boundaries.

    Output
    ------
    The function creates following files:
        1. Single VTK file containing information of all grains.
        2. Individual VTK files containing information of individual grains.
        3. Single OBJ file containing information of all grains.
        4. Individual OBJ files containing information of individual grains.
        5. INP file consisting of mesh information
        6. Text file of structural characteristics
        7. Text file of textural characteristics
        8. Distribution plots 
    """
    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Executing common execute_func module')

    ## Identify the name of function that has called this function
    parent_function_name = inspect.stack()[1][3]

    ## Creating tessellations
    tessellation = create_tessellations(seed_array_unique, limit, log_level)

    log.info('Starting to compute structural characteristics')

    ## Structural Characteristics
    grain_size_distributions = np.around(grain_size_distribution(dimension, tessellation, limit, log_level), decimals=10)
    number_of_neighbor = number_of_neighbors(dimension, tessellation, log_level)
    grain_boundary_area_distribution, all_vertices_list = grain_boundary_areas(dimension, limit, tessellation, parent_function_name, skewed_boundary_flag, log_level)
    junction_lengths = junction_length(tessellation, log_level)
    junction_angles_degrees = junction_angle(tessellation, log_level)
    distance_btw_grain_array, distance_btw_grain_1d = distance_btw_grains(dimension, limit, tessellation, log_level)
    
    log.info('Successfully computed all structural characteristics')

    log.info('Starting to compute textural characteristics')

    ## Textural Characteristics
    disorientation_angle, orientation_data = disorientation_angles(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    schmid_factors, orientation_data = schmid_factor(required_texture, rand_quat_flag, dimension, stress_direction, orientation_data, tessellation, log_level)
    type_of_grain_boundaries, orientation_data = type_of_grain_boundary(required_texture, rand_quat_flag, orientation_data, tessellation, log_level)
    
    log.info('Successfully computed all textural characteristics')

    log.info('Starting to create visualization files in VTK and OBJ formats')

    ## Writing visualization files (VTK and OBJ files)
    create_vtk_file_all_grains(material, tessellation, store_folder, face_flag, now, skewed_boundary_flag, all_vertices_list, log_level)
    create_vtk_file_individual_grains(material, tessellation, store_folder, now, skewed_boundary_flag, all_vertices_list, log_level)
    create_obj_file_all_grains(material, tessellation, store_folder, face_flag, now, skewed_boundary_flag, all_vertices_list, log_level)
    create_obj_file_individual_grains(material, tessellation, store_folder, now, skewed_boundary_flag, all_vertices_list, log_level)

    log.info('Successfully created visualization files in VTK and OBJ formats')

    ## Meshing
    ## If mesh_flag is None then meshing will not be performed
    if mesh_flag:
        
        log.info('Meshing configuration')

        ## Depending on the required mesh type, function call name would be identified from dict
        mesh_type_dictionary ={'HEX': mesh_hex, 'TET': mesh_tetra, 'VIS': mesh_visualization}
        
        if (mesh_flag == 'HEX') or (mesh_flag == 'TET'):
            mesh_type_dictionary[mesh_flag](limit, global_mesh_size)

            '''
            The mesh function stores the nodes and elements data into the inp file
            hence in order to do further processing we need to extract nodes and 
            volume elements data from inp file.
            '''
            nodes = []
            elements = []
            node_flag = False
            element_flag = False
            
            node_flag = False                                                       # Flags to indicate lines of inp file belonging to nodes data
            element_flag = False                                                    # Flags to indicate lines of inp file belonging to volume element data
            ## read nodes from inp file
            with open('mesh_data.inp', 'r') as f:
                data = f.readlines()
                for line in data:
                    ## As soon as '*NODE' is identified in the lines the node flag is activated
                    if '*NODE' in line:                                     
                        node_flag = True
                        continue
                    
                    ## As soon as '**' is identified in the lines the node flag is deactivated
                    if '**' in line:
                        node_flag = False
                    
                    ## C3D8 and C3D4 are the volume element types for hexahedral and tetrahedral mesh respectively
                    ## As soon as '*ELEMENT, type=C3D8' is identified in the lines the elements flag is activated
                    if ('*ELEMENT, type=C3D8' in line) or ('*ELEMENT, type=C3D4' in line):
                        element_flag = True
                        continue

                    ## Storing nodes and elements data to the list
                    if node_flag:
                        nodes.append(list(map(float, list((line.rstrip('\n')).split(',')))))
                    
                    if element_flag:
                        elements.append(list(map(int, list((line.rstrip('\n')).split(',')))))
            
            ## Converting to arrays
            nodes = np.array(nodes)
            elements = np.array(elements)

            '''
            Since we are interested only in volume elements, we will be deleting the line
            and surface elements in future. Therefore we need to change the element 
            numbers to make it start from proper sequence. We are subtracting 9 since
            8 lines are reserved for 8 corners of the cube and elements numbering
            should start from 9.
            '''
            
            elements[:,0] -= (elements[0,0] - 9)

            ## Determine centroids of elements
            elements_centroid = []
            for row in elements:
                # First column is the element number
                node_numbers = row[1:]-1                                            # Subtracting by 1 since in inp file node numbering starts from 1

                element_nodal_array = nodes[node_numbers, :]                        # Array of all nodal coordinates of the particular element
                
                ## Determining centroid
                centroid_coordinate = np.mean(element_nodal_array[:, 1:4], axis= 0)
                elements_centroid.append(centroid_coordinate)
            
            ## Converting to array
            elements_centroid = np.array(elements_centroid)
            
            ## Extracting centroids of each grains
            grain_centroids = np.array(copy.deepcopy(tessellation['centroid_list'])) #[v.centroid() for v in tessellation])
            
            ## Adding Grain Numbers to grain centroids array
            grain_numbers = np.arange(grain_centroids.shape[0])
            grain_numbers = grain_numbers[:, np.newaxis]
            grain_centroids = np.hstack((grain_numbers, grain_centroids))

            ## Creating a Permutation matrix for considering Periodicity
            x_column = np.repeat(np.repeat(np.array([0, limit[0], -limit[0]]), 3), 3)
            y_column = np.tile(np.repeat(np.array([0, limit[1], -limit[1]]), 3), 3)
            z_column = np.tile(np.array([0, limit[2], -limit[2]]), 9)
            periodicity_matrix = np.vstack((x_column, y_column, z_column)).T
            
            ## Checking for grains within simulation box
            ## Here we are trying to isolate grain numbers that have vertices outside the simulation box
            protruding_grains = []
            for grain_index in range(copy.deepcopy(tessellation['number_of_grains'])):
                vertices_simulation_box = np.array([[0, 0, 0],
                                                    [limit[0], 0, 0],
                                                    [limit[0], limit[1], 0],
                                                    [0, limit[1], 0],
                                                    [0, 0, limit[2]],
                                                    [limit[0], 0, limit[2]],
                                                    [limit[0], limit[1], limit[2]],
                                                    [0, limit[1], limit[2]]])
                
                # Create a convex hull of the simulation box vertices
                hull = ConvexHull(vertices_simulation_box, incremental=True)
                old_hull_vertices = hull.vertices                               # extract vertices of the generated hull

                ## Add the vertices of the grain
                grain_vertices = np.array(copy.deepcopy(tessellation['vertices_list'][grain_index]))
                hull.add_points(grain_vertices, restart=False)
                new_hull_vertices = hull.vertices
                
                ## Check if the vertices of the hull is different that the old vertices of hull
                # If yes, then it means that the grain has vertices outside simulation box
                if np.array_equal(old_hull_vertices, new_hull_vertices) is False:
                    protruding_grains.append(grain_index)
            
            ## Translating the centroids of enlisted grains in all possible 
            # permutations and appending it separately to the centroids array
            for protruding_grain in protruding_grains:
                centroid = grain_centroids[protruding_grain, 1:]
                broadcast_centroid = np.broadcast_to(centroid, (periodicity_matrix.shape))        
                translated_periodic_centroids = broadcast_centroid + periodicity_matrix
                
                current_grain_numbers = np.repeat(protruding_grain, translated_periodic_centroids.shape[0])
                current_grain_numbers = current_grain_numbers[:, np.newaxis]
                translated_periodic_centroids = np.hstack((current_grain_numbers, translated_periodic_centroids))
                
                grain_centroids = np.vstack((grain_centroids, translated_periodic_centroids))
                
            ## Removing any duplicate centroids
            new_centroids_array = [tuple(row) for row in grain_centroids]
            grain_centroids = np.unique(new_centroids_array, axis = 0)
                    
            ## Identifying the grains of each element
            grain_number_of_each_element = np.zeros([elements.shape[0], 1]).flatten()
            for element_number in range(elements.shape[0]):
                broadcast_element_centroid = np.tile(elements_centroid[element_number, :], (grain_centroids.shape[0], 1))   ## To match the number of rows
                subtract_array = broadcast_element_centroid - grain_centroids[:, 1:]
                #distance_between_centroids = np.linalg.norm(subtract_array, axis=1)
                distance_between_centroids = np.einsum('ij,ij->i', subtract_array, subtract_array)
                grain_number_of_each_element[element_number] = grain_centroids[np.argmin(distance_between_centroids), 0]    ## updating the grain number with the shortest distance
                
            ## Writing to new INP file
            output_file_path = Path(store_folder, material, now, "Text_output", "mesh_" + mesh_flag + ".inp")
            output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
            output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
            output_file_path.parent.parent.parent.mkdir(exist_ok=True)
            output_file_path.parent.parent.mkdir(exist_ok=True)
            output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one
            with open(str(output_file_path), 'a+') as f:
                f.truncate(0)
                
                ## Writing header
                f.write('*heading \n')
                f.write('   abaqus_input_file.inp \n')
                
                ## Writing Nodes
                f.write('*NODE \n')
                for row_number in range(nodes.shape[0]):
                    f.write(str(int(nodes[row_number, 0])) + ', ' + str(list(nodes[row_number, 1:])).lstrip('[').rstrip(']') + '\n')
                
                ## Writing elements
                f.write('******* E L E M E N T S ************* \n')
                
                ## Writing element type based on type of meshing
                if (mesh_flag == 'HEX'):
                    f.write('*ELEMENT, type=C3D8, ELSET=Volume34 \n')
                else:
                    f.write('*ELEMENT, type=C3D4, ELSET=Volume34 \n')
                
                
                for row_number in range(elements.shape[0]):
                    f.write(str(list(elements[row_number, :])).lstrip('[').rstrip(']') + '\n')
                
                ## Creating Element sets
                for i in grain_numbers:
                    f.write('\n')
                    # f.write('*ELEMENT, type=C3D8, ELSET=Grain' + str(i) + '\n')
                    f.write('*ELSET, ELSET=Grain' + str(i) + '\n')
                    
                    ## Identifying elements belonging to particular grain
                    elements_index_for_grain = np.where(grain_number_of_each_element == i)
                    elements_number_list = elements[elements_index_for_grain[0], 0]
                    
                    ## An elset line can have only 16 entities 
                    elsets_line = []
                    for count in range(len(elements_number_list)):
                        elsets_line.append(elements_number_list[count])
                        if (count + 1)%16 == 0 and count != 0:
                            f.write(str(list(elsets_line)).lstrip('[').rstrip(']') + ', ' + '\n')
                            elsets_line = []
                            continue
                        elif count >= (len(elements_number_list) - 1):
                            f.write(str(list(elsets_line)).lstrip('[').rstrip(']') + '\n')
                            f.write('*Solid Section, ELSET=Grain' + str(i) + ', material=material_Grain' + str(i) + '\n')   ## Creating solid sections
            
            os.remove('mesh_data.inp')

            log.info('Successfully completed with meshing')

        elif (mesh_flag == 'VIS'):
            mesh_type_dictionary[mesh_flag](tessellation, global_mesh_size)

            ## Saving Mesh file as 'MSH' and 'INP'
            for extension in ['inp', 'msh']:
                output_file_path = Path(store_folder, material, now, "Text_output", "mesh_" + mesh_flag + "." + extension)
                output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
                output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
                output_file_path.parent.parent.parent.mkdir(exist_ok=True)
                output_file_path.parent.parent.mkdir(exist_ok=True)
                output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one 

                with open('mesh_data_file.' + extension, 'r') as old_file, open(str(output_file_path), 'a+') as new_file:
                    new_file.truncate(0)
                    data = old_file.read()
                    new_file.write(data)

                os.remove('mesh_data_file.' + extension)
            
            log.info('Successfully completed with meshing')

        else:
            log.error("Please enter appropriate mesh type. Refer documentation for more details. \n")


    ## Saving Structural and Txtural characteristics to a file
    
    log.debug('Saving structural characteristics into a text file')

    ## Saving Structural CHaracteristics
    output_file_path = Path(store_folder, material, now, "Text_output", "structural_characteristics.txt")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    with open(str(output_file_path), 'a+') as f:
        f.truncate(0)
        f.write(now)
        
        f.write("\n\nSize of Simulation Box: " + str(size_of_simulation_box))
        #f.write("\nNumber of Seeds: " + str(number_of_seeds))
        f.write("\nDimension: " + str(dimension))
        f.write("\nMaterial Name: " + material)
        f.write("\n \n# Grain sizes \n")
        f.write('# Grain Number, Grain Size \n')
        np.savetxt(f, grain_size_distributions, delimiter=',', fmt="%.4f")
        f.write("\n# Number of Neighbors \n")
        f.write("# Grain Number, Number of neighbors, List of indices of all neighboring grains \n")
        for line in number_of_neighbor:
            np.savetxt(f, np.array(line), newline=' ', delimiter=',', fmt="%.4f")
            f.write("\n")
        f.write("\n# Grain Boundary Areas \n")
        f.write("# Sr. No., Grain Number 1, Grain Number 2, Area \n")
        np.savetxt(f, grain_boundary_area_distribution, delimiter=',', fmt="%.4f")
        f.write("\n# Junction Lengths \n")
        f.write("# Sr. no., Junction type, Junction lengths, Grains with this junction \n")
        for line in junction_lengths:
            np.savetxt(f, np.array(line), newline=' ', delimiter=',', fmt="%.4f")
            f.write("\n")
        f.write("\n# Junction angles degrees \n")
        f.write('# Sr. No., Junction Type, 1st Junction angle, Grain containing the 1st junction angle, 2nd Junction angle, so on... \n')
        for line in junction_angles_degrees:
            np.savetxt(f, np.array(line), newline=' ', delimiter=',', fmt="%.4f")
            f.write("\n")
        f.write('\n# Distance between grains \n')
        f.write('# Sr. No., Distance between grains in ascending orders \n')
        np.savetxt(f, distance_btw_grain_array, delimiter=',', fmt="%.4f")

    log.info('Successfully saved all structural characteristics into a text file')

    log.debug('Saving textural characteristics into a text file')

    ## Saving Textural CHaracteristics
    output_file_path = Path(store_folder, material, now, "Text_output", "textural_characteristics.txt")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    with open(str(output_file_path), 'a+') as f:
        f.truncate(0)
        f.write(now)
        
        f.write("\n\nSize of Simulation Box: " + str(size_of_simulation_box))
        #f.write("\nNumber of Seeds: " + str(number_of_seeds))
        f.write("\nDimension: " + str(dimension))
        f.write("\nMaterial Name: " + material)
        f.write("\n \n# Disorientation Angles \n")
        f.write("\n# Grain 1, Grain 2, Disorientation angle, Disorientation axis \n")
        np.savetxt(f, disorientation_angle, delimiter=',', fmt="%.4f")
        f.write("\n \n# Type of Grain Boundaries \n")
        f.write("\n# Grain 1, Grain 2, Misorientation angle, Misorientation axis, Type of Grain Boundary \n")
        np.savetxt(f, type_of_grain_boundaries, delimiter=',', fmt='%.4f %.4f %.4f %.4f %.4f %.4f %s')
        f.write("\n# Schmid Factors")
        f.write("\n# Grain number, Maximum schmid factor, Respective slip plane and slip direction (Slip System) \n")
        np.savetxt(f, schmid_factors, delimiter=',', fmt="%.4f")

    log.info('Successfully saved all textural characteristics into a text file')

    log.debug('Saving seeds data into a text file')

    ## Saving seeds data
    output_file_path = Path(store_folder, material, now, "Text_output", "seed_data.txt")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    with open(str(output_file_path), 'a+') as f:
        f.truncate(0)
        f.write("# " + now + " \n")
        f.write("# Seeds data (Seed coordinates + Orientation as Quaternions)\n")
        f.write("# X coordinate, Y coordinate, Z coordinate, W, q1, q2, q3 \n")
        seed_data = np.around(np.concatenate((seed_array_unique, orientation_data), axis=1), decimals=14)
        np.savetxt(f, seed_data, delimiter=',', comments='#', fmt='%.14f')

    log.info('Successfully saved seeds data into a text file')

    log.debug('Starting to save all plots')

    ## Plotting Structural Characteristics
    fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize = (25,25))
    font_size_value = 40
    label_size = 25
    
    ## Plotting Grain SIze Distribution
    hist, bins = np.histogram(grain_size_distributions[:, 1], bins= number_of_bins, density= True)
    #bins = 0.5 * (bins[1:] + bins[:-1])
    ax[0][0].plot(bins[:-1], hist)
    ax[0][0].scatter(bins[:-1], hist)
    ax[0][0].set_xlabel("Grain Sizes", fontsize=font_size_value)
    ax[0][0].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
    ax[0][0].tick_params(labelsize=label_size)
    
    ## Plotting Number of Neighbors
    neighbors_array = np.array([v[1] for v in number_of_neighbor])
    hist, bins = np.histogram(neighbors_array, bins= number_of_bins, density= True)
    #bins = 0.5 * (bins[1:] + bins[:-1])
    ax[0][1].plot(bins[:-1], hist)
    ax[0][1].scatter(bins[:-1], hist)
    ax[0][1].set_xlabel("Number of Neighbors", fontsize=font_size_value)
    ax[0][1].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
    ax[0][1].tick_params(labelsize=label_size)
    
    ## Plotting Grain Boundary Areas
    all_areas = [v[3] for v in grain_boundary_area_distribution]
    hist, bins = np.histogram(all_areas, bins= number_of_bins, density= True)
    #bins = 0.5 * (bins[1:] + bins[:-1])
    ax[0][2].plot(bins[:-1], hist)
    ax[0][2].scatter(bins[:-1], hist)
    ax[0][2].set_xlabel("Grain Boundary Areas", fontsize=font_size_value)
    ax[0][2].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
    ax[0][2].tick_params(labelsize=label_size)
    
    ## Plotting Junction Lengths
    all_lengths = [v[2] for v in junction_lengths]
    hist, bins = np.histogram(all_lengths, bins= number_of_bins, density= True)
    #bins = 0.5 * (bins[1:] + bins[:-1])
    ax[1][0].plot(bins[:-1], hist)
    ax[1][0].scatter(bins[:-1], hist)
    ax[1][0].set_xlabel("Junction Lengths", fontsize=font_size_value)
    ax[1][0].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
    ax[1][0].tick_params(labelsize=label_size)
    
    ## Plotting Junction angles
    all_angles = [v[2::2] for v in junction_angles_degrees]
    all_angles_flatten = [w for v in all_angles for w in v]
    hist, bins = np.histogram(all_angles_flatten, bins= number_of_bins, density= True)
    #bins = 0.5 * (bins[1:] + bins[:-1])
    ax[1][1].plot(bins[:-1], hist)
    ax[1][1].scatter(bins[:-1], hist)
    ax[1][1].set_xlabel("Junction angles", fontsize=font_size_value)
    ax[1][1].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
    ax[1][1].tick_params(labelsize=label_size)

    ## Plotting Distance between grains
    all_distances = distance_btw_grain_1d
    hist, bins = np.histogram(all_distances, bins= number_of_bins, density= True)
    ax[1][2].plot(bins[:-1], hist)
    ax[1][2].scatter(bins[:-1], hist)
    ax[1][2].set_xlabel("Distance between grains", fontsize=font_size_value)
    ax[1][2].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
    ax[1][2].tick_params(labelsize=label_size)

    plt.suptitle("Structural Characteristics", fontsize=60)
    
    ##Saving Plot
    output_file_path = Path(store_folder, material, now, "Plots", "structural_characteristics.png")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one
    
    
    plt.subplots_adjust(wspace=0.5)
    plt.tight_layout
    fig.savefig(str(output_file_path))
    
    ## Plotting Textural Characteristics
    fig, ax = plt.subplots(ncols = 3, figsize = (30,15))

    ## Plotting Disorientation angles
    all_disorientation_angles = [v[2] for v in disorientation_angle]
    hist, bins = np.histogram(all_disorientation_angles, bins= number_of_bins, density= True)
    #bins = 0.5 * (bins[1:] + bins[:-1])
    ax[0].plot(bins[:-1], hist)
    ax[0].scatter(bins[:-1], hist)
    ax[0].set_xlabel("Disorientation angles", fontsize=font_size_value)
    ax[0].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
    ax[0].tick_params(labelsize=label_size)
    
    ## Plotiing Schmid Factors
    all_schmid_factors = schmid_factors[:, 1]
    hist, bins = np.histogram(all_schmid_factors, bins = number_of_bins, density= True)
    #bins = 0.5 * (bins[1:] + bins[:-1])
    ax[1].plot(bins[:-1], hist)
    ax[1].scatter(bins[:-1], hist)
    ax[1].set_xlabel("Schmid Factors", fontsize=font_size_value)
    ax[1].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
    ax[1].tick_params(labelsize=label_size)

    ## Plotting Type of Grain Boundaries
    all_grain_boundaries = type_of_grain_boundaries[:, 6]
    hist, bins = np.histogram(all_grain_boundaries, bins = number_of_bins, density= True)
    #bins = 0.5 * (bins[1:] + bins[:-1])
    ax[2].plot(bins[:-1], hist)
    ax[2].scatter(bins[:-1], hist)
    ax[2].set_xlabel("Type of Grain Boundary (CSL type)", fontsize=font_size_value)
    ax[2].set_ylabel("Frequency of occurrences", fontsize=font_size_value)
    ax[2].tick_params(labelsize=label_size)

    plt.suptitle("Textural Characteristics", fontsize=60)
    
    ##Saving Plot
    output_file_path = Path(store_folder, material, now, "Plots", "textural_characteristics.png")
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one
    
    plt.tight_layout
    plt.subplots_adjust(wspace=0.5)
    fig.savefig(str(output_file_path))

    log.info('Successfully saved all required plots')

    return grain_size_distributions, number_of_neighbor, grain_boundary_area_distribution, junction_lengths, junction_angles_degrees, distance_btw_grain_array, distance_btw_grain_1d, disorientation_angle, schmid_factors, type_of_grain_boundaries