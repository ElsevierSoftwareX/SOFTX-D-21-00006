# -*- coding: utf-8 -*-
"""
mesh.py

Module to perform meshing.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

@authors: Serrao Prince Henry, Arun Prakash
@email: 
created: 28 January 2020
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

import gmsh
import numpy as np

from ppp2019_optimizedmicrostructuregeneration.src.set_logger import set_logger as set_logger

def mesh_hex(limit, global_mesh_size):
    '''
    Using gmsh-sdk module to generate hexahedral mesh for the simulation box.

    Input: The function requires array of limits along X, Y, Z directions and 
            the global mesh size.

    Process: The function generates hexahedral mesh in the simulation box such 
            that it would be assigned to the respective grains in the later 
            stages.

    Output: The function writes an INP file with the mesh information for 
            further processing.
    '''

    ## Extracting size of simulation box
    x_limit = limit[0]
    y_limit = limit[1]
    z_limit = limit[2]
    
    ## Setting global mesh size to 0.5
    lcar = global_mesh_size

    ## Defining the number of nodes along all directions to generate hexahedral mesh
    nodes_x = int(x_limit/lcar) + 1
    nodes_y = int(y_limit/lcar) + 1
    nodes_z = int(z_limit/lcar) + 1
    
    ## Initializing GMSH
    gmsh.initialize()

    ## Print message on terminal
    gmsh.option.setNumber("General.Terminal", 1)

    ## Setting model name
    gmsh.model.add("hexahedral_mesh")

    ## Defining the 8 corners of the simulation box
    # First three args are the coordinates, then the mesh size around that point 
    # and then the Point Tag number.
    gmsh.model.geo.addPoint(0, 0, 0, lcar, 1)
    gmsh.model.geo.addPoint(x_limit, 0, 0, lcar, 2)
    gmsh.model.geo.addPoint(x_limit, y_limit, 0, lcar, 3)
    gmsh.model.geo.addPoint(0, y_limit, 0, lcar, 4)
    gmsh.model.geo.addPoint(0, 0, z_limit, lcar, 5)
    gmsh.model.geo.addPoint(x_limit, 0, z_limit, lcar, 6)
    gmsh.model.geo.addPoint(x_limit, y_limit, z_limit, lcar, 7)
    gmsh.model.geo.addPoint(0, y_limit, z_limit, lcar, 8)

    ## Defining the lines and the start and end points of the lines
    # First two args are the Tags of points connected by the line and then the 
    # Line tag number 
    gmsh.model.geo.addLine(1, 2, 9)
    gmsh.model.geo.addLine(2, 3, 10)
    gmsh.model.geo.addLine(3, 4, 11)
    gmsh.model.geo.addLine(4, 1, 12)
    gmsh.model.geo.addLine(5, 6, 13)
    gmsh.model.geo.addLine(6, 7, 14)
    gmsh.model.geo.addLine(7, 8, 15)
    gmsh.model.geo.addLine(8, 5, 16)
    gmsh.model.geo.addLine(1, 5, 17)
    gmsh.model.geo.addLine(2, 6, 18)
    gmsh.model.geo.addLine(3, 7, 19)
    gmsh.model.geo.addLine(4, 8, 20)

    ## Defining each surface of the cube
    # It is important to note that the direction of lines while specifying them
    # in the loop. First argument is the List of lines in the loop and then the
    # CurveLoop tag number
    gmsh.model.geo.addCurveLoop([9, 10, 11, 12], 21)
    gmsh.model.geo.addPlaneSurface([21], 22)                                    # First arg is the Curve loop and then the surface Tag number

    gmsh.model.geo.addCurveLoop([13, 14, 15, 16], 23)
    gmsh.model.geo.addPlaneSurface([23], 24)

    gmsh.model.geo.addCurveLoop([9, 18, -13, -17], 25)
    gmsh.model.geo.addPlaneSurface([25], 26)

    gmsh.model.geo.addCurveLoop([10, 19, -14, -18], 27)
    gmsh.model.geo.addPlaneSurface([27], 28)

    gmsh.model.geo.addCurveLoop([11, 20, -15, -19], 29)
    gmsh.model.geo.addPlaneSurface([29], 30)

    gmsh.model.geo.addCurveLoop([12, 17, -16, -20], 31)
    gmsh.model.geo.addPlaneSurface([31], 32)

    ## Defining the loop of the surfaces to enclose a volume
    # First argument is the list of the surfaces and then the SurfaceLoop Tag number.
    gmsh.model.geo.addSurfaceLoop([26, 28, 30, 32, 22, 24], 33)
    gmsh.model.geo.addVolume([33], 34)                                          # Defining the volume. First arg is the Surface Loop tag number and then the Volume Tag number
    
    ## Defining the Transfinite Lines with specific number of nodes 
    # to generate hexahedral mesh
    # First arg is the line tag number, then the number of nodes and then the 
    # mesh type. It is important to note that the lines along the same direction 
    # must have same number of nodes.
    gmsh.model.geo.mesh.setTransfiniteCurve(9, nodes_x, meshType = "Progression")
    gmsh.model.geo.mesh.setTransfiniteCurve(11, nodes_x, meshType = "Progression")
    gmsh.model.geo.mesh.setTransfiniteCurve(13, nodes_x, meshType = "Progression")
    gmsh.model.geo.mesh.setTransfiniteCurve(15, nodes_x, meshType = "Progression")

    gmsh.model.geo.mesh.setTransfiniteCurve(10, nodes_y, meshType = "Progression")
    gmsh.model.geo.mesh.setTransfiniteCurve(12, nodes_y, meshType = "Progression")
    gmsh.model.geo.mesh.setTransfiniteCurve(16, nodes_y, meshType = "Progression")
    gmsh.model.geo.mesh.setTransfiniteCurve(14, nodes_y, meshType = "Progression")

    gmsh.model.geo.mesh.setTransfiniteCurve(17, nodes_z, meshType = "Progression")
    gmsh.model.geo.mesh.setTransfiniteCurve(18, nodes_z, meshType = "Progression")
    gmsh.model.geo.mesh.setTransfiniteCurve(19, nodes_z, meshType = "Progression")
    gmsh.model.geo.mesh.setTransfiniteCurve(20, nodes_z, meshType = "Progression")

    ## Defining the Transfinite Surface
    # First arg is the surface TAg number and then the corners point Tags of the points 
    # belonging to the surface
    gmsh.model.geo.mesh.setTransfiniteSurface(22, cornerTags= [1, 2, 3, 4])
    gmsh.model.geo.mesh.setTransfiniteSurface(24, cornerTags= [5, 6, 7, 8])
    gmsh.model.geo.mesh.setTransfiniteSurface(26, cornerTags= [1, 2, 6, 5])
    gmsh.model.geo.mesh.setTransfiniteSurface(28, cornerTags= [2, 3, 7, 6])
    gmsh.model.geo.mesh.setTransfiniteSurface(30, cornerTags= [3, 4, 8, 7])
    gmsh.model.geo.mesh.setTransfiniteSurface(32, cornerTags= [4, 1, 5, 8])

    ## Defining Transfinite volume with first arg as the Volume Tag and then the 
    # corners of the Volume
    gmsh.model.geo.mesh.setTransfiniteVolume(34, cornerTags= [1, 2, 3, 4, 5, 6, 7, 8])

    ## Recombine is the most important step to generate hexahedral mesh
    # Fist arg is the dimension, then is the Surface Tag and then the angles by 
    # which the triangles are to be recombined to quadrangles
    gmsh.model.geo.mesh.setRecombine(2, 22, angle = 45)
    gmsh.model.geo.mesh.setRecombine(2, 24, angle = 45)
    gmsh.model.geo.mesh.setRecombine(2, 26, angle = 45)
    gmsh.model.geo.mesh.setRecombine(2, 28, angle = 45)
    gmsh.model.geo.mesh.setRecombine(2, 30, angle = 45)
    gmsh.model.geo.mesh.setRecombine(2, 32, angle = 45)
    
    ## Synchronizing the model before generating mesh
    gmsh.model.geo.synchronize()
    
    ## Generating 3D mesh (argument 3 for 3D)
    gmsh.model.mesh.generate(3)

    ## Storing the model to MSH and INP files
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write("mesh_data.inp")

    gmsh.finalize()

def mesh_tetra(limit, global_mesh_size):
    '''
    Using gmsh-sdk module to generate tetrahedral mesh for the simulation box.

    Input: The function requires array of limits along X, Y, Z directions and 
            the global mesh size.

    Process: The function generates tetrahedral mesh in the simulation box such 
            that it would be assigned to the respective grains in the later 
            stages.

    Output: The function writes an INP file with the mesh information for 
            further processing.
    '''

    ## Extracting size of simulation box
    x_limit = limit[0]
    y_limit = limit[1]
    z_limit = limit[2]
    
    ## Setting global mesh size to 0.5
    lcar = global_mesh_size

    ## Defining the number of nodes along all directions to generate hexahedral mesh
    nodes_x = int(x_limit/lcar) + 1
    nodes_y = int(y_limit/lcar) + 1
    nodes_z = int(z_limit/lcar) + 1
    
    ## Initializing GMSH
    gmsh.initialize()

    ## Print message on terminal
    gmsh.option.setNumber("General.Terminal", 1)

    ## Setting model name
    gmsh.model.add("hexahedral_mesh")

    ## Defining the 8 corners of the simulation box
    # First three args are the coordinates, then the mesh size around that point 
    # and then the Point Tag number.
    gmsh.model.geo.addPoint(0, 0, 0, lcar, 1)
    gmsh.model.geo.addPoint(x_limit, 0, 0, lcar, 2)
    gmsh.model.geo.addPoint(x_limit, y_limit, 0, lcar, 3)
    gmsh.model.geo.addPoint(0, y_limit, 0, lcar, 4)
    gmsh.model.geo.addPoint(0, 0, z_limit, lcar, 5)
    gmsh.model.geo.addPoint(x_limit, 0, z_limit, lcar, 6)
    gmsh.model.geo.addPoint(x_limit, y_limit, z_limit, lcar, 7)
    gmsh.model.geo.addPoint(0, y_limit, z_limit, lcar, 8)

    ## Defining the lines and the start and end points of the lines
    # First two args are the Tags of points connected by the line and then the 
    # Line tag number 
    gmsh.model.geo.addLine(1, 2, 9)
    gmsh.model.geo.addLine(2, 3, 10)
    gmsh.model.geo.addLine(3, 4, 11)
    gmsh.model.geo.addLine(4, 1, 12)
    gmsh.model.geo.addLine(5, 6, 13)
    gmsh.model.geo.addLine(6, 7, 14)
    gmsh.model.geo.addLine(7, 8, 15)
    gmsh.model.geo.addLine(8, 5, 16)
    gmsh.model.geo.addLine(1, 5, 17)
    gmsh.model.geo.addLine(2, 6, 18)
    gmsh.model.geo.addLine(3, 7, 19)
    gmsh.model.geo.addLine(4, 8, 20)

    ## Defining each surface of the cube
    # It is important to note that the direction of lines while specifying them
    # in the loop. First argument is the List of lines in the loop and then the
    # CurveLoop tag number
    gmsh.model.geo.addCurveLoop([9, 10, 11, 12], 21)
    gmsh.model.geo.addPlaneSurface([21], 22)                                    # First arg is the Curve loop and then the surface Tag number

    gmsh.model.geo.addCurveLoop([13, 14, 15, 16], 23)
    gmsh.model.geo.addPlaneSurface([23], 24)

    gmsh.model.geo.addCurveLoop([9, 18, -13, -17], 25)
    gmsh.model.geo.addPlaneSurface([25], 26)

    gmsh.model.geo.addCurveLoop([10, 19, -14, -18], 27)
    gmsh.model.geo.addPlaneSurface([27], 28)

    gmsh.model.geo.addCurveLoop([11, 20, -15, -19], 29)
    gmsh.model.geo.addPlaneSurface([29], 30)

    gmsh.model.geo.addCurveLoop([12, 17, -16, -20], 31)
    gmsh.model.geo.addPlaneSurface([31], 32)

    ## Defining the loop of the surfaces to enclose a volume
    # First argument is the list of the surfaces and then the SurfaceLoop Tag number.
    gmsh.model.geo.addSurfaceLoop([26, 28, 30, 32, 22, 24], 33)
    gmsh.model.geo.addVolume([33], 34)                                          # Defining the volume. First arg is the Surface Loop tag number and then the Volume Tag number
    
    ## Synchronizing the model before generating mesh
    gmsh.model.geo.synchronize()
    
    ## Generating 3D mesh (argument 3 for 3D)
    gmsh.model.mesh.generate(3)

    ## Storing the model to MSH and INP files
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write("mesh_data.inp")

    gmsh.finalize()

def mesh_visualization(tessellation, global_mesh_size):
    '''
    Using gmsh-sdk module to generate tetrahedral mesh for the grains 
    configuration without limiting them to the simulation box.

    Input: The function requires tessellations data and the global mesh size.

    Process: The function generates tetrahedral mesh in all the grains individually
            and ensures that the surfaces are not repeated in order to ensure
            conformity between the grains.

    Output: The function writes an INP file with the mesh information for 
            configuration without limiting it to the simulation box and 
            considering periodic boundary conditions (PBC).
    '''

    ## Setting global mesh size to 0.5
    lcar = global_mesh_size

    ## Initializing GMSH
    gmsh.initialize()

    ## Print message on terminal
    gmsh.option.setNumber("General.Terminal", 1)

    ## Setting model name
    gmsh.model.add("tetra_mesh")

    ## Initialize line counter
    line_counter = 0

    ## Initializing an array where Tag of first point for that particular cell will be stored
    vertex_line_data = np.zeros(len(tessellation))
    print(vertex_line_data)
    ## Iterating to add all vertices to geometry
    for index, v in enumerate(tessellation):
        ## Extracting the array of vertices of the particular cell
        vertex_array = np.array(v.vertices())
        ## Updating the Tag of first point of the cell which will be added
        vertex_line_data[index] = line_counter
        ## Adding points to GMSH geometry and incrementing line counter
        for vertex in vertex_array:
            gmsh.model.geo.addPoint(vertex[0], vertex[1], vertex[2], lcar, line_counter)
            line_counter += 1

    ## Creating an empty list where tuples of end points of a created line would be stored
    list_created_lines = []
    ## Creating an empty list where Tag number corresponding to the tuple of line from above list, would be stored 
    line_line_data = []
    ############################################################################
    ## Initializing a list where the centroid data with Tag numbers of surface 
    ## would be stored. The first 3 elements correspond to centroid coords and 
    ## next is the PlaneSurface Tag number. Initializing with a list of zeros so 
    ## that the centroidal distance for first surface could be calculated.
    ############################################################################
    centroid_all_face_data = [[0, 0, 0, 0]]
    for index, v in enumerate(tessellation):
        ## Extracting vertices of cell and face vertices of a cell
        face_vertices_array = v.face_vertices()
        vertex_array = np.array(v.vertices())

        ## Initializing a list to store the TAG numbers of all the surfaces enclosing a volume
        plane_surface_line_counters = []
        ## Iterating through each surface(face) of the cell
        for row_num, row in enumerate(face_vertices_array):
            ## In order to ensure that the edges of the surface form a closed loop
            row.append(row[0])
            ## Creating an empty list to store the TAG numbers of all the lines 
            # belonging to the face
            lines_of_face = []

            ## Iterating through the face vertices pair wise
            for start_point, last_point in zip(row[:-1], row[1:]):
                ############################################################################
                ## Creating the tuple of the end points of a line as 'first_tuple',
                ## Creating the tuple of swapped end points of a line as 'second_tuple'.
                ## Here it is important to add the starting index of the points 
                ## of the particular cell in order to ensure that the points of 
                ## earlier cells are not considered while creating lines.
                ############################################################################
                first_tuple = tuple([start_point + int(vertex_line_data[index]), last_point+ int(vertex_line_data[index])])
                second_tuple = tuple([last_point+ int(vertex_line_data[index]), start_point + int(vertex_line_data[index])])
                
                ## Checking if the line already exists
                if (first_tuple in list_created_lines):
                    ############################################################################
                    ## Appending the TAG number of the line to the lines of the face.
                    ## First finding the index where the tuple matches in the list_created_lines
                    ## and then using that index extracting the TAG number of the
                    ## line from line_line_data.
                    ############################################################################
                    lines_of_face.append(line_line_data[list_created_lines.index(first_tuple)])
                    continue
                elif (second_tuple in list_created_lines):
                    ############################################################################
                    ## While appending to the list it is important to ensure that
                    ## the surface would be closed by the sequence of lines.
                    ## Suppose line (1, 3) already exists and 'second tuple' is (1, 3), means
                    ## that start is 3 and end is 1 that is reverse of existing line
                    ## hence multiply by (-1).
                    ############################################################################
                    lines_of_face.append((-1)*line_line_data[list_created_lines.index(second_tuple)])
                    continue
                else:
                    ## If none of the above conditions matches then create a line
                    gmsh.model.geo.addLine(start_point + int(vertex_line_data[index]), last_point + int(vertex_line_data[index]), line_counter)
                    ## append the TAG number to the lines of face list
                    lines_of_face.append(line_counter)
                    line_counter += 1
                    ## Appending the tuple of new line and the TAG number to 
                    ## list created lines list and line_line_data list respectively
                    list_created_lines.append(first_tuple)
                    line_line_data.append(line_counter-1)
            
            ## Important to ensure synchronization before extracting boundary details
            gmsh.model.geo.synchronize()
            ## Extracting the TAGs of points of the lines.
            ## Tuple consists of dimension of entity and TAG of point
            point_tags_of_face = np.array([gmsh.model.getBoundary((1, np.abs(i)), recursive=True) for i in lines_of_face])
            ## Extracting only the TAG numbers from the 3d array
            flatten_point_tags_of_face = (point_tags_of_face[:, :, 1]).flatten()
            ## Creating an array of unique TAG of points and deducing the row of 
            # the vertex array by subtracting the starting TAG of point of each cell
            unique_point_tags_of_face = np.array(list(set(flatten_point_tags_of_face))) - int(vertex_line_data[index])
            ## Extracting the coordinates the vertices of the face
            face_vertices_coordinates_array = vertex_array[unique_point_tags_of_face, :]
            ## Finding the centroid of the surface
            centroid_face = np.mean(face_vertices_coordinates_array, axis=0)
            ## Broadcasting the centroid of face 
            broadcast_centroid_face = np.tile(centroid_face, (len(np.array(centroid_all_face_data)), 1))
            ## Subtracting the current face centroid with all the other face centroids
            subtract_arrays = broadcast_centroid_face - np.array(centroid_all_face_data)[:, :-1]
            ## Finding the Euclideandistance between the centroids
            distance_between_centroids = np.linalg.norm(subtract_arrays, axis=1)
            ## FInding the index where the distance is the lowest
            lowest_distance_index = np.argmin(distance_between_centroids)

            ## For the surfaces to be overlapping it is important that the distance 
            # between centroids is zero and append plane surface TAG number
            if np.isclose(distance_between_centroids[lowest_distance_index], 0):
                plane_surface_line_counters.append(np.array(centroid_all_face_data)[lowest_distance_index, 3])
                continue
            ## If it is a uniques surface then Create a surface
            gmsh.model.geo.addCurveLoop(lines_of_face, line_counter)
            line_counter += 1
            gmsh.model.geo.addPlaneSurface([line_counter - 1], line_counter)
            ## appending the centroid and the Plane Surface TAG number of the unique face
            centroid_all_face_data.append(list(centroid_face) + [line_counter])
            plane_surface_line_counters.append(line_counter)
            line_counter += 1
        
        ## Creating the surface loop of all the surfaces enclosing a volume
        gmsh.model.geo.addSurfaceLoop(plane_surface_line_counters, line_counter)
        line_counter += 1
        ## Creating a volume from the surface loop
        gmsh.model.geo.addVolume([line_counter - 1], line_counter)
        line_counter += 1
        ## Creating a SET based on the volume
        gmsh.model.addPhysicalGroup(3, [line_counter-1], line_counter)
        line_counter += 1
        gmsh.model.setPhysicalName(3, line_counter - 1, "Grain_" + str(index))
    
    ## Remove all the duplicate entities(nodes) and synchornize before generating mesh
    gmsh.model.geo.removeAllDuplicates()
    gmsh.model.geo.synchronize()    
    
    ## Generating 3D mesh (argument 3 for 3D)
    #gmsh.model.mesh.setOrder(2)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen", force=False)

    ## Storing the model to MSH and INP files
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write("mesh_data_file.msh")
    gmsh.write("mesh_data_file.inp")

    gmsh.finalize()