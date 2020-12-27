# -*- coding: utf-8 -*-
"""
create_obj_files.py

Module to create OBJ files for visualization.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

@authors: Serrao Prince Henry, Arun Prakash
@email: 
created: 16 November 2019
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

def create_obj_file_all_grains(material_name, tessellation, store_folder, face_flag, now, skewed_boundary_flag, all_vertices_list, log_level):
    """
    The function creates a directory "obj_file_all_grains" if it doesnt 
    exist in the 'visualization_files' directorywithin appropriate 
    sub-directories. An OBJ file consisting data related to all the grains is 
    saved in this directory.

    Parameters
    ----------
    material_name: string
        Name of the material of which microstructure is being generated.
    
    tessellation: dictionary
        Dictionary consisting of data related to tessellations generated.
    
    store_folder: string
        Sub-directory name (output/output_test) in which output data would be
        stored. 
    
    face_flag: boolean
        Flag to specify if opaque surfaces are to be used instead of transparent.

    now: string
        Current time and date.
        
    skewed_boundary_flag: boolean
        Flag to specify if skewed grain boundaries are required. Only functional
        in quasi-2D case.
        
    all_vertices_list: list
        List of vertices of skewed faces in appropriate sequence.
        
    log_level: string
        Logger level to be used.

    Output
    ------
    The functions saves OBJ file consisting of all the vertices and the line 
    connectivity data for all the vertices in the respective folder.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Creating OBJ file of entire configuration')
    ## Checking if opaque surfaces are required or transparent
    if face_flag:
        type_of_surface = "f "                  # 'f' indicates opaque face
    else:
        type_of_surface = "l "                  # 'l' indicates connection of vertices by polyline

    ## Creating Folder and assigning the output file path
    output_file_path = Path(store_folder, material_name, now, "obj_file_all_grains", "{:s}.obj".format(material_name + "_morphology"))
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one
    

    ## Using generator statement to extract data from all cells of the tessellation
    if skewed_boundary_flag:
        grain_vertices = all_vertices_list
    else:
        grain_vertices = copy.deepcopy(tessellation['vertices_list']) #[v.vertices() for v in tessellation]                   # all the vertices of each grains
    grain_faces = copy.deepcopy(tessellation['face_vertices_list']) #[v.face_vertices() for v in tessellation]
    
    ## Writing the required commands of the OBJ format to the file
    with open(str(output_file_path), 'a+') as f:
        
        f.truncate(0)                                                               # Deleting all contents of the file
        f.write("# All grains morphology \n")
        
        ## Writing all vertices to the file
        vertices_indices = [0]
        for vertex in grain_vertices:
            vertices_indices.append(len(vertex) + vertices_indices[-1])                        
            for vertex_coords_list in vertex:
                vertex_coords = np.array(vertex_coords_list)
                f.write ("v " + str(vertex_coords).strip('[').strip(']') + "\n")
        
        ## Writing line connectivity data for all the faces of each cell
        for grain_number, faces in enumerate(grain_faces):
            for face_indices in faces:
                face_indice = np.array(face_indices)
                f.write(type_of_surface + str(face_indice + vertices_indices[grain_number] + 1).strip('[').strip(']') + "\n")
    
    log.info('Successfully created OBJ file of entire configuration')

def create_obj_file_individual_grains(material_name, tessellation, store_folder, now, skewed_boundary_flag, all_vertices_list, log_level):
    """
    The function creates a directory "obj_file_individual_grains" if it doesnt 
    exist in the 'visualization_files' directory within appropriate 
    sub-directories. Individual OBJ files consisting data related to respective 
    grain is saved in this folder. By default only opaque surfaces are created 
    for individual grains.

    Parameters
    ----------
    material_name: string
        Name of the material of which microstructure is being generated.
    
    tessellation: dictionary
        Dictionary consisting of data related to tessellations generated.
    
    store_folder: string
        Sub-directory name (output/output_test) in which output data would be
        stored. 
    
    now: string
        Current time and date.
        
    skewed_boundary_flag: boolean
        Flag to specify if skewed grain boundaries are required. Only functional
        in quasi-2D case.
        
    all_vertices_list: list
        List of vertices of skewed faces in appropriate sequence.
        
    log_level: string
        Logger level to be used.

    Output
    ------
    The functions saves individual OBJ files consisting of all the vertices and 
    the line connectivity data for all the vertices in the respective folder.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Creating OBJ files of individual grains')

    ## If skewed boundary is required then modified vertices list is used
    if skewed_boundary_flag:
        grain_vertices = all_vertices_list
    else:
        grain_vertices = copy.deepcopy(tessellation['vertices_list']) #[v.vertices() for v in tessellation]                   # all the vertices of each grains
    grain_faces = copy.deepcopy(tessellation['face_vertices_list']) #[v.face_vertices() for v in tessellation]
    
    for i in range(copy.deepcopy(tessellation['number_of_grains'])):
        
        log.debug('Creating OBJ file of grain no. ' + str(i))
        
        ## Creating Folder and assigning the output file path
        output_file_path = Path(store_folder, material_name, now, "obj_file_individual_grains", "{:s}.obj".format("grain_" + str(i+1)))
        output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.mkdir(exist_ok=True)

        with open(str(output_file_path), 'a+') as f:
            f.truncate(0)
            f.write("# Individual grain morphology \n")

            for vertex_list in grain_vertices[i]:
                vertex = np.array(vertex_list)
                f.write ("v " + str(vertex).strip('[').strip(']') + "\n")
            f.write("\n \n")

            for faces in grain_faces[i]:
                face_indice = np.array(faces)
                f.write("f " + str(face_indice + 1).strip('[').strip(']') + "\n")
            f.write("\n \n")

    log.info('Successfully created OBJ files of all grains individually')