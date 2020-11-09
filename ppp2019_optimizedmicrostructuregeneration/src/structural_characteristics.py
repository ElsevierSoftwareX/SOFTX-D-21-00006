# -*- coding: utf-8 -*-
"""
structural_characteristics.py

Module to compute the following structural characteristics:
1. Grain Size
2. Number of Neighbors
3. Grain Boundary Area
4. Junction length
5. Junction Angles
6. Distance between grains

@author: Serrao Prince Henry, Dr. Arun Prakash

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

created: 23 November 2019
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

#from memory_profiler import profile
#unit normal vector of plane defined by points a, b, and c
## THIS FUNCTION IS PART OF FACE AREA CALCULATION
## SOURCE: https://stackoverflow.com/questions/12642256/python-find-area-of-polygon-from-xyz-coordinates @ 'Jamie Bull'
## aLGORITHM SOURCE: http://geomalgorithms.com/a01-_area.html#3D%20Polygons by Dan Sunday
def unit_normal(a, b, c):
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#area of polygon poly
## THIS FUNCTION IS PART OF FACE AREA CALCULATION
## SOURCE: https://stackoverflow.com/questions/12642256/python-find-area-of-polygon-from-xyz-coordinates @ 'Jamie Bull'
## aLGORITHM SOURCE: http://geomalgorithms.com/a01-_area.html#3D%20Polygons by Dan Sunday
def poly_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)

def grain_size_distribution(dimension, tessellation_og, limit, log_level):
    """
    Compute grain sizes.

    Processing
    ----------
    In case of 2D, the volume is the area of the cell as the length along z axis
    is 1. Assuming the grains to be a circle in 2D and sphere in 3D, the radius 
    of the grain is determined using the area and volume formula respectively.

    Parameter
    ---------
    dimension: integer    
        Dimension of study (2 or 3).

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

    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

    log_level: string
        Logger level to be used.

    Returns
    -------
    The function returns an array of all the grain sizes in terms of radius. The
    first column is the grain number, second column consists of grain sizes and
    third column consists of grain volume/area.
    """

    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Started computing grain sizes')

    tessellation = copy.deepcopy(tessellation_og)
    grain_sizes = np.zeros([tessellation['number_of_grains'], 3])               # Columns are: gr. number, gr. size, gr. volume/area
    grain_numbers = np.array([num for num in range(tessellation['number_of_grains'])])
    grain_volumes = np.array(tessellation['volume_list'])                # Extracting the volume of each grain from the cells data
    
    grain_sizes[:, 0] = grain_numbers[:]
    grain_sizes[:, 1] = grain_volumes[:]
    grain_sizes[:, 2] = grain_volumes[:]

    ## Calculating grain sizes  based on dimension 
    if dimension == 2:
        grain_sizes[:, 1] = grain_sizes[:, 1]/limit[2]                          # Grain size is considered to be based on XY plane, hence removing effect of length along Z axis
        grain_sizes[:, 1] = 2*(grain_sizes[:, 1]/np.pi)**0.5                    # Multiplying by 2 in order to get diameter
        grain_sizes[:, 2] = grain_sizes[:, 2]/limit[2]                          # grain area in case of quasi-2D

    elif dimension == 3:
        grain_sizes[:, 1] = 2*(grain_sizes[:, 1] * 3 / (np.pi * 4))**(1/3)      # Multiplying by 2 in order to get diameter

    log.info('Completed computing grain sizes')    
    return np.around(grain_sizes, decimals=12)

def number_of_neighbors(dimension, tessellation_og, limit, log_level):
    """
    Compute number of neighbors for all grains.

    Processing
    ----------
    In case of 2D, the number of neighbors along z axis are ignored.

    Parameter
    ---------
    dimension: integer    
        Dimension of study (2 or 3).

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

    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

    log_level: string
        Logger level to be used.

    Returns
    -------
    The function returns an array comprising of data related to the number of 
    neighbors of each grains. Column names are: Grain number, number of 
    neighbors, grain indexes of neighbors.
    """

    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Started computing number of neighbors')
    tessellation = copy.deepcopy(tessellation_og)

    ## Calculating the number of neighbors based on dimension
    grain_neighbors = []                                                        # column names are: Grain number, grain volume/area, number of neighbors, grain indexes of neighbors
    if dimension == 2:
        for grain_index in range(tessellation['number_of_grains']):
            individual_gr_neighbors = []                                        # Initializing an empty list where individual grain data would be stored
            individual_gr_neighbors.append(grain_index)
            individual_gr_neighbors.append(tessellation['volume_list'][grain_index]/limit[2])   # Grain area in case of quasi-2D
            individual_gr_neighbors.append(len(tessellation['neighbors_list'][grain_index])- dimension)   # Ignoring self grain number from list of neighbors
            individual_gr_neighbors = individual_gr_neighbors + [val for val in tessellation['neighbors_list'][grain_index] if val != grain_index] # concatenating lists of grain indices of neighboring grains
            grain_neighbors.append(individual_gr_neighbors)                     # Appending individual grain data to main list
        
    elif dimension == 3:
        for grain_index in range(tessellation['number_of_grains']):
            individual_gr_neighbors = []
            individual_gr_neighbors.append(grain_index)
            individual_gr_neighbors.append(tessellation['volume_list'][grain_index])   # Grain volume
            individual_gr_neighbors.append(len(tessellation['neighbors_list'][grain_index]))
            individual_gr_neighbors = individual_gr_neighbors + [val for val in tessellation['neighbors_list'][grain_index]]
            grain_neighbors.append(individual_gr_neighbors)

    log.info('Completed computing number of neighbors')
    return grain_neighbors

def grain_boundary_areas(dimension, limit, tessellation_og, parent_function_name, skewed_boundary_flag, log_level):
    """
    Compute grain boundary areas.

    Processing
    ----------
    In case of 2D, it is required to ignore the faces perpendicular to the 
    z axis. In order to achieve this the absolute value along the z axis of the 
    normal of the face is checked. If the normal has unit value along z axis 
    then that face is ignored. Also we consider the 2D case to be a quasi 2D 
    and displace the vertices on the negative z plane by a small shift based on 
    randoms. It is important to note that even in case of Quasi-2D grain 
    boundary has an area. 

    Parameters
    ----------
    dimension: integer    
        Dimension of study (2 or 3)

    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

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

    parent_function_name: string
        Name of the function that is calling this function.
        
    skewed_boundary_flag: boolean 
        Flag to specify if skewed grain boundaries are required. Only functional
        in quasi-2D case.
        
    log_level: string
        Logger level to be used.
    
    Returns
    -------
    The function returns an array of consisting of columns Sr. No., Grain 1, 
    Grain 2 and grain boundary area.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Started computing grain boundary areas')

    tessellation = copy.deepcopy(tessellation_og)

    if dimension == 2:                                                          # Checking for the dimension
        if ('test' in parent_function_name):
            value = False
        else: 
            value = skewed_boundary_flag
        
        if value  == False:

            row_counter = 0
            grain_boundary_area = []                                                     # Initializing array
            for grain_number in range(tessellation['number_of_grains']):
                normal_of_face = tessellation['normals_list'][grain_number]                                        # extracting the normals of faces of cell

                for face_number, normals in enumerate(normal_of_face):
                    single_pair_data = []
                    if np.isclose(np.linalg.norm(normals[2]), 1.0):                 # checking the absolute value of normal along z axis
                        continue                                                    # ignoring face if absolute value along z axis is 1
                    else:
                        single_pair_data.append(row_counter)                        # Appending the sr. no.
                        single_pair_data.append(grain_number)                       # Appending first grain
                        single_pair_data.append((tessellation['neighbors_list'][grain_number])[face_number])   # Appending second grain
                        single_pair_data.append((tessellation['face_area_list'][grain_number])[face_number])  # Appending the GB area
                        grain_boundary_area.append(single_pair_data)                # Appending the row to main array
                        row_counter += 1

            grain_boundary_area = np.around(np.array(grain_boundary_area), decimals=14)
            all_vertices_list = None
            
        else:
            face_area_list = []                                                     # Initializing array
            all_gb_face_vertices = []
            all_neighbors_list = []
            
            for v in range(tessellation['number_of_grains']):
                normal_of_face = tessellation['normals_list'][v]                                        # extracting the normals of faces of cell
                gb_face_vertices = []
                grain_face_area_list =[]                                            # face area of each grain
                current_grain_neighbors = []

                ## Extracting the face vertices of faces having normals not parallel to z axis
                for face_number, normals in enumerate(normal_of_face):
                    if np.isclose(np.linalg.norm(normals[2]), 1.0):                 # checking the absolute value of normal along z axis
                        continue                                                    # ignoring face if absolute value along z axis is 1
                    else:
                        grain_face_area_list.append((tessellation['face_area_list'][v])[face_number])
                        gb_face_vertices.append((tessellation['face_vertices_list'][v])[face_number])
                        current_grain_neighbors.append((tessellation['neighbors_list'][v])[face_number])
                face_area_list.append(grain_face_area_list)
                all_gb_face_vertices.append(gb_face_vertices)
                all_neighbors_list.append(current_grain_neighbors)

            ## Find the quasi vertex shift parameneter considering it to be 10% of the minimum grain boundary length of each grain
            all_vertices_list = tessellation['vertices_list']
            min_gb_length_list = np.array([np.min(areas/limit[2]) for areas in face_area_list])
            
            #average_gb_length = np.average(face_area_list)
            quasi_variable_d = 0.2 * min_gb_length_list

            changed_vertices_grain_vertex_index = []                                 # List to store tuple of grain and vertex number of all changed vertices 
            for grain_index, grain in enumerate(all_vertices_list):
                for vertex_index, vertex in enumerate(grain):
                    vertex = np.array(vertex)
                    all_vertices_list[grain_index][vertex_index] = list(vertex)     # Changing type from tuple to list
                    
                    ## Fing the vertices having negative z coordinate and checking if the vertex is already shifted
                    ############################################################
                    # Since it is periodic the grains are extended outside of simulation box
                    # to represent periodicity. The simulation box is still in
                    # positive values.
                    ############################################################
                    if (vertex[2] <=0 and (~((grain_index, vertex_index) in changed_vertices_grain_vertex_index))):#((vertex[0] > 0 and vertex[0]<limit[0]) and (vertex[1] > 0 and vertex[1] < limit[1]))):   
                        all_vertices_list[grain_index][vertex_index][0] += quasi_variable_d[grain_index]
                        all_vertices_list[grain_index][vertex_index][1] += quasi_variable_d[grain_index]
                        
                        ## Changing all matching vertices by same amount
                        for grain_index_match, grain_match in enumerate(all_vertices_list):
                            for vertex_index_match, vertex_match in enumerate(grain_match):
                                
                                if np.allclose(vertex_match, vertex):
                                    changed_vertices_grain_vertex_index.append((grain_index_match, vertex_index_match))
                                    all_vertices_list[grain_index_match][vertex_index_match] = list(vertex)
                                    all_vertices_list[grain_index_match][vertex_index_match][0] += quasi_variable_d[grain_index]
                                    all_vertices_list[grain_index_match][vertex_index_match][1] += quasi_variable_d[grain_index]

            ## Calculating the total number of grain boundaries
            total_gb = 0
            for grain in all_gb_face_vertices:
                for faces in grain:
                    total_gb += 1

            ## Initializing an array where all the GB related information will be stored
            grain_boundary_area = np.zeros([total_gb, 4])                           # column names: Sr. no, grain 1, grain 2, area
            grain_boundary_area[:, 0] = range(total_gb)
            
            ## Computing the GB areas and updating the storage array
            row_counter = 0
            for grain_index, grain in enumerate(all_gb_face_vertices):

                neighbors = all_neighbors_list[grain_index]#[val for val in tessellation['neighbors_list'][grain_index] if val != grain_index]

                for face_index, face in enumerate(grain):
                    face_vertex_coordinates = []
                    for vertex in face:
                        face_vertex_coordinates.append(all_vertices_list[grain_index][vertex])
                    face_area = poly_area(face_vertex_coordinates)                  # Computing the area based on the 3D coordinates of vertices
                    grain_boundary_area[row_counter, 1] = grain_index
                    grain_boundary_area[row_counter, 2] = neighbors[face_index]
                    grain_boundary_area[row_counter, 3] = face_area
                    row_counter += 1

            grain_boundary_area = np.around(np.array(grain_boundary_area), decimals=14)

    elif dimension == 3:
        all_vertices_list = None
        face_area_list = tessellation['face_area_list']                 # Storing all the face areas as grain boundary areas for 3D case
        
        ## Calculating the total number of grain boundaries
        total_gb = 0
        for grain in face_area_list:
            total_gb += len(grain)
        
        ## Initializing an array where all the GB related information will be stored
        grain_boundary_area = np.zeros([total_gb, 4])                           # column names: Sr. no, grain 1, grain 2, area
        grain_boundary_area[:, 0] = range(total_gb)

        ## Updating the data into storage array
        row_counter = 0
        for grain_index, grain in enumerate(face_area_list):
            neighbors = tessellation['neighbors_list'][grain_index]

            for face_index, area in enumerate(grain):
                grain_boundary_area[row_counter, 1] = grain_index
                grain_boundary_area[row_counter, 2] = neighbors[face_index]
                grain_boundary_area[row_counter, 3] = area
                row_counter += 1

        grain_boundary_area  = np.around(np.array(grain_boundary_area), decimals=14)
    
    log.info('Completed computing grain boundary areas')
    return grain_boundary_area, all_vertices_list

def junction_length(dimension, tessellation_og, log_level):
    """
    Compute junction lengths.

    Processing
    ----------
    The lengths of edges of grains are considered as junction lengths. It is 
    important to note that in case of Quasi-2D, junction lengths along Z axis 
    are also considered. Junction lengths are stored as list of lists as the
    lists are of irregular sizes.

    Parameters
    ----------
    dimension: integer    
        Dimension of study (2 or 3)

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
    List of lists consisting of column names: Sr. no., Junction type, Junction 
    lengths, grains with this junction.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Started computing junction lengths')

    tessellation = copy.deepcopy(tessellation_og)
    ## Extracting all the required informations from cells data
    """
    I dont know why but I had to deep copy all the lists of dictionary since when
    this function was called again, all face vertices list somehow retained the
    appended face vertex which existed even if I stripped last element of 'face'
    list.
    """
    all_vertices_list = tessellation['vertices_list']
    number_of_edges = tessellation['number_of_edges_list'] #[v.number_of_edges() for v in tessellation]
    total_number_of_edges = np.sum(number_of_edges)
    all_face_vertices_list = tessellation['face_vertices_list'] #[v.face_vertices() for v in tessellation]
    all_face_normals = tessellation['normals_list'] #[v.normals() for v in tessellation]
    all_neighbors_of_each_grain = tessellation['neighbors_list'] #[v.neighbors() for v in tessellation]
    
    ## Initializing arrays and lists
    junction_lengths = []                                                       # Column names: Sr. no., Junction type, Junction lengths, grains with this junction
    edge_vertices_array = np.zeros([int(total_number_of_edges), 7])             # First six columns for edge vertex coordinates, 7th column for grain number
    grains_with_junction = []
    edges_of_junction = []

    two_D_flag = False
    if dimension == 2:
        two_D_flag = True
        z_plane_normal = np.array([0, 0, 1])

    ## Adding Unique vertices of each edge
    """
    Each row of the 2D array consists of 3 columns of coordinates of first vertex 
    of edge and remaining 3 columns as coordinates of the other vertex of edge.
    """
    row_counter = 0
    for grain_index, grain in enumerate(all_face_vertices_list):
        grain_row_counter = row_counter
        for face_index, face in enumerate(grain):
            face.append(face[0])
            ## iterating through the adjacent face vertices
            for i, j in zip(face[:], face[1:]):
                row = all_vertices_list[grain_index][i] + all_vertices_list[grain_index][j]
                swapped_row = all_vertices_list[grain_index][j] + all_vertices_list[grain_index][i]
                
                ## Checking if the row of edge data exists in the array by swapping the vertices of edge
                no_of_edges_for_current_grain = int(grain_row_counter + number_of_edges[grain_index])                     # Checking is done only within the edges of particular grain 

                if two_D_flag:
                    check_row = np.array(row)
                    edge_vector = np.abs(check_row[:3] - check_row[3:6])
                    if np.isclose(np.dot(edge_vector, z_plane_normal), 0):
                        continue

                if any((np.equal(row, edge_vertices_array[grain_row_counter:no_of_edges_for_current_grain, :6])).all(1)):
                    continue

                if any((np.equal(swapped_row, edge_vertices_array[grain_row_counter:no_of_edges_for_current_grain, :6])).all(1)):
                    continue
                edge_vertices_array[row_counter, :6] = np.array(row[:])
                edge_vertices_array[row_counter, 6] = grain_index
                row_counter += 1

    edge_vertices_array = np.around(edge_vertices_array[~(edge_vertices_array==0).all(1)], decimals=10)

    ## Counting the number of repetitions of particular edge
    """
    If the edge is shared by 3 grains then it is considered as triple junction and length is calculated by 
    using the vertices of the edge.
    """

    i = 0
    while i < len(edge_vertices_array):
        
        row = edge_vertices_array[i]                                            # assigning row of an array to a variable                                  
        
        ## Extracting the coordinates of each vertex of the edge
        edge_1 = row[:3]                                                    
        edge_2 = row[3:6]
        
        ## Concatenating the edges in both possible ways
        row_array = np.concatenate((edge_1, edge_2), axis=None)
        swapped_row_array = np.concatenate((edge_2, edge_1), axis=None)

        ## Finding index numbers of the matching edges
        comparing_row_array = np.all((edge_vertices_array[:, :6] == row_array), axis = 1)
        comparing_swapped_row_array = np.all((edge_vertices_array[:, :6] == swapped_row_array), axis = 1)
        logical_or_operation = np.logical_or(comparing_row_array, comparing_swapped_row_array)
        repeating_edge_indices = np.where(logical_or_operation)
        
        ## Appending the length of the required type of junction to the list 'junction_lengths'
        single_junction_data = []
        single_junction_data.append(i)
        single_junction_data.append(len(repeating_edge_indices[0]))
        single_junction_data.append(np.linalg.norm(np.array(edge_2)-np.array(edge_1)))
        grain_numbers = [edge_vertices_array[row_index, 6] for row_index in repeating_edge_indices[0]]
        single_junction_data = single_junction_data + grain_numbers
        junction_lengths.append(single_junction_data)
        
        #grain_numbers = [edge_vertices_array[row_index, 6] for row_index in repeating_edge_indices[0]]
        grains_with_junction.append(grain_numbers)
        edges_of_junction.append(row_array)


        ## Deleting all the rows which have been already accounted for
        edge_vertices_array = np.delete(edge_vertices_array, repeating_edge_indices[0][1:], axis=0)
        
        i += 1                                                                  

    ## Returning additional parameters when this function is called by 'junction_angle' function
    if inspect.stack()[1][3] == 'junction_angle':                               # Identifies the function name calling this
        log.info('Completed computing junction lengths')
        return grains_with_junction, edges_of_junction
    
    log.info('Completed computing junction lengths')
    return junction_lengths

def junction_angle(dimension, tessellation_og, log_level):
    """
    Compute junction angles

    Parameters
    ----------
        Parameters
    ----------
    dimension: integer    
        Dimension of study (2 or 3)

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
    List of lists consisting of column names: Sr. No., Junction Type, 
    1st Junction angle, Grain containing the 1st junction angle, 
    2nd Junction angle, so on...
    """

    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Started computing junction angles')

    tessellation = copy.deepcopy(tessellation_og)
    ## Calling 'junction_length' function for identifying the required edges with the grains having these edges
    grains_with_junction, edges_of_junction = junction_length(dimension, tessellation, log_level)
    
    ## Creating again a copy in order to avoid error due to call by reference
    tessellation = copy.deepcopy(tessellation_og)
    ## Initializing a list where all the junction angles related information would be stored
    junction_angles = []                                                        ## Column names: Sr. No., Junction Type, 1st Junction angle, Grain containing the 1st junction angle, 2nd Junction angle, so on...           

    ## iterating over the Data to find angles between the face normals
    for row_index, row in enumerate(np.array(edges_of_junction)):
        
        ## splitting the row data into edges with coordinates
        edge_1 = row[:3]
        edge_2 = row[3:6]
        
        ## Initializing a local list where the data pertaining to a junction would be stored
        single_junction_related_data = []                                       ## List where all info related to a particular junction will be stored
        single_junction_related_data.append(row_index)
        single_junction_related_data.append(len(grains_with_junction[row_index]))

        for grain_index, grain in enumerate(np.array(grains_with_junction[row_index])):
            
            ## Extracting data regarding specific grains
            vertices_individual_grain = np.around(np.array(tessellation['vertices_list'][int(grain)]), decimals=10) #tessellation[int(grain)].vertices()), decimals=4)
            faces_individual_grain = tessellation['face_vertices_list'][int(grain)] #tessellation[int(grain)].face_vertices()
            
            ## Initializing an empty list where the normals of faces forming the required junction edges will be stored
            normals_of_junction_faces = []
            
            ## Identifying the index numbers of the edge vertices in the array of individual grain vertices
            comparing_edge_1 = np.all(vertices_individual_grain == edge_1, axis = 1)
            comparing_edge_2 = np.all(vertices_individual_grain == edge_2, axis=1)
            logical_or_operation = np.logical_or(comparing_edge_1, comparing_edge_2)
            junction_vertex_indexes = np.where(logical_or_operation)
            
            ## Iterating over all the faces of the grain to identify faces having the required junction vertices
            for face_index, face in enumerate(faces_individual_grain):
                face.append(face[0])
                for i, j in zip(face[:], face[1:]):
                    
                    ## When the face is identified containing the required edge vertices, normals of these faces are appended
                    if (i == junction_vertex_indexes[0][0] and j == junction_vertex_indexes[0][1]) or (i == junction_vertex_indexes[0][1] and j == junction_vertex_indexes[0][0]):
                        normal_of_face = (tessellation['normals_list'][int(grain)])[face_index] #(tessellation[int(grain)].normals())[face_index]
                        normals_of_junction_faces.append(normal_of_face)
            
            ## Rounding off the normal vectors to 14 decimals
            normals_of_junction_faces = np.around(np.array(normals_of_junction_faces), decimals=14)
            if len(normals_of_junction_faces) == 0: 
                continue
            
            ## Calculating the cosine theta between the normal vectors
            cosine_theta = np.dot(normals_of_junction_faces[0],normals_of_junction_faces[1])/np.linalg.norm(normals_of_junction_faces[0])/np.linalg.norm(normals_of_junction_faces[1])
            
            ## Storing angles in the storage array
            angle_inscribed_by_grain_at_junction = 180 - np.degrees(np.arccos(cosine_theta))    # Subtracting to get angles between the faces
            single_junction_related_data.append(angle_inscribed_by_grain_at_junction)
            single_junction_related_data.append(grain)
        junction_angles.append(single_junction_related_data)                    ## Appending the junction information to the main list    

    log.info('Completed computing junction angles')
    return junction_angles

def distance_btw_grains(dimension, limit, tessellation_og, log_level):
    """
    Compute distance of each grain from all grains.

    Processing
    ----------
    Centroid of all grains are replicated along all directions based on 
    dimension of study. The shortest distance between certain pair of grains out
    of all the symmetric combinations is considered as the final distance for 
    that pair of grains.
    
    Parameters
    ----------
    dimension: integer    
        Dimension of study (2 or 3)

    limit: array
        Size of simulation box (array of length along X, Y, Z directions)

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
    1. Symmetric array with rows and columns represented by grain numbers in 
    ascending order and each element of the array representing distance between 
    respective grain numbers.
    2. 1D array of distances between grains
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.debug('Started computing distance between grains')

    tessellation = copy.deepcopy(tessellation_og)
    ## Initialize array to store distance information
    distance_array = np.zeros([tessellation['number_of_grains'], tessellation['number_of_grains']])

    ## Extracting centroids of each grains
    centroid_grains = np.array(tessellation['centroid_list'])

    ## Creating periodic centroids list
    periodic_centroid_list = [centroid_grains]

    if dimension == 2:
        for ix in [-1, 0, 1]:
            for iy in [-1, 0, 1]:
                if (not((ix ==0) and (iy ==0))):
                    
                    current_centroid_grain = np.zeros_like(centroid_grains)
                    current_centroid_grain[:, 0] = centroid_grains[:, 0] + (ix*limit[0])
                    current_centroid_grain[:, 1] = centroid_grains[:, 1] + (iy*limit[1])
                    current_centroid_grain[:, 2] = centroid_grains[:, 2]

                    periodic_centroid_list.append(current_centroid_grain)

    elif dimension == 3:
        for ix in [-1, 0, 1]:
            for iy in [-1, 0, 1]:
                for iz in [-1, 0, 1]:
                    if (not((ix ==0) and (iy ==0) and (iz ==0))):
                        current_centroid_grain = np.zeros_like(centroid_grains)
                        current_centroid_grain[:, 0] = centroid_grains[:, 0] + (ix*limit[0])
                        current_centroid_grain[:, 1] = centroid_grains[:, 1] + (iy*limit[1])
                        current_centroid_grain[:, 2] = centroid_grains[:, 2] + (iz*limit[2])

                        periodic_centroid_list.append(current_centroid_grain)

    ## Computing distance between centroids
    for row_index, row in enumerate(centroid_grains):

        distance_array_all_periodic = np.zeros((tessellation['number_of_grains'], len(periodic_centroid_list)))
        
        for replicated_index, replicated_array in enumerate(periodic_centroid_list):
            
            ## Broadcasting row
            broadcast_row = np.tile(row, (tessellation['number_of_grains'], 1))
            ## Subtracting arrays row wise
            subtract_arrays = broadcast_row - replicated_array
            ## Computing Euclidean distance
            distance = np.linalg.norm(subtract_arrays, axis=1)
            distance_array_all_periodic[:, replicated_index] = distance

        distance_array[row_index, :] = np.min(distance_array_all_periodic, axis=1)
    
    ## Computing distance between neighbors
    neighbors_list_data = tessellation['neighbors_list']
    distance_btw_neighbors_1d = []                                                          # empty list where distance between neighbors would be stored
    for index_grain, grain in enumerate(neighbors_list_data):
        for index_neighbor, neighbor_number in enumerate(grain):
            if (index_grain != neighbor_number):
                distance_btw_neighbors_1d.append(distance_array[index_grain, neighbor_number])

    log.info('Completed computing distance between grains')
    return distance_array, distance_btw_neighbors_1d
    