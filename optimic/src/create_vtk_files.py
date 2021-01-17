# -*- coding: utf-8 -*-
"""
create_vtk_files.py

Module to create VTK files for visualization.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
created: 16 November 2019
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

from src.main_import_statements import *

from src.set_logger import set_logger 
name_str = __name__

def create_vtk_file_all_grains(material_name, tessellation, store_folder, face_flag, now, skewed_boundary_flag, all_vertices_list, log_level):
    """
    The function creates a directory "vtk_file_all_grains" if it doesnt 
    exist in the 'visualization_files' directorywithin appropriate 
    sub-directories. An VTK file consisting data related to all the grains is 
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
    The functions saves VTK file consisting of all the vertices and the line 
    connectivity data for all the vertices in the respective folder.
    """
    
    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Creating VTK file of entire configuration')

    ## Checking if opaque surfaces are required or transparent
    if face_flag:
        type_of_surface = 7                             # 7 indicates opaque surface
    else:
        type_of_surface = 4                             # 4 indicates connection of all the vertices by polyline

    ## Creating Folder and assigning the output file path
    output_file_path = Path(store_folder, material_name, now, "vtk_file_all_grains", "{:s}.vtk".format(material_name + "_morphology"))
    output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.parent.mkdir(exist_ok=True)
    output_file_path.parent.mkdir(exist_ok=True)                                # checks if folder exists, if not then creates one

    ## Using generator statement to extract data from all cells of the tessellation
    if skewed_boundary_flag:
        grain_vertices = all_vertices_list                            # Using skewed vertices list
    else:
        grain_vertices = copy.deepcopy(tessellation['vertices_list']) #[v.vertices() for v in tessellation]         # all the vertices of each grains

    number_of_grain_vertices = np.sum([len(i) for i in grain_vertices])
    grain_faces = copy.deepcopy(tessellation['face_vertices_list']) #[v.face_vertices() for v in tessellation]
    number_of_grain_faces = np.sum([len(i) for i in grain_faces])
    
    ## Writing the required commands to the file
    with open(str(output_file_path), 'a+') as f:
        f.truncate(0)
        
        ## Writing the Header lines
        f.write("# vtk DataFile Version 3.0 \n")
        f.write("Morphology of all grains \n")
        f.write("ASCII \n")
        f.write("DATASET UNSTRUCTURED_GRID \n")
        f.write("POINTS " + str(number_of_grain_vertices) + " float \n")
        
        ## Writing all vertices to the file
        vertices_indices = [0]
        for vertex in grain_vertices:
            vertices_indices.append(len(vertex) + vertices_indices[-1])                        
            for vertex_coords_list in vertex:
                vertex_coords = np.array(vertex_coords_list)
                f.write (str(vertex_coords).strip('[').strip(']').lstrip(' ').rstrip(' ') + "\n")
        f.write("\n \n")
        
        ## Writing the Cells (Face vertices) data to file
        sum = 0
        for faces in grain_faces:
            for face_indices in faces:
                sum += len(face_indices)
        
        f.write("CELLS " + str(number_of_grain_faces) + " " + str((number_of_grain_faces) + sum) + "\n")
        
        for grain_number, faces in enumerate(grain_faces):
            for face_indices in faces:
                face_indice = np.array(face_indices)
                f.write(str(len(face_indices)) + " " + str(face_indice + vertices_indices[grain_number]).strip('[').strip(']') + "\n")
        f.write("\n \n")
        
        ## Writing Cell Type data to file
        f.write("CELL_TYPES " + str(number_of_grain_faces) + "\n")
        for i in range(number_of_grain_faces):
            f.write(str(type_of_surface) + "\n")

    log.info('Successfully created VTK file of entire configuration')


def create_vtk_file_individual_grains(material_name, tessellation, store_folder, now, skewed_boundary_flag, all_vertices_list, log_level):
    """
    The function creates a directory "vtk_file_individual_grains" if it doesnt 
    exist in the 'visualization_files' directory within appropriate 
    sub-directories. Individual VTK files consisting data related to respective 
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
    The functions saves individual VTK files consisting of all the vertices and 
    the line connectivity data for all the vertices in the respective folder.
    """

    log = set_logger(name_str, 'log_data.log', log_level)
    log.info('Creating VTK files of individual grains')

    if skewed_boundary_flag:
        grain_vertices = all_vertices_list                            # Using skewed vertices list
    else:
        grain_vertices = copy.deepcopy(tessellation['vertices_list']) #[v.vertices() for v in tessellation]         # all the vertices of each grains
    grain_faces = copy.deepcopy(tessellation['face_vertices_list']) #[v.face_vertices() for v in tessellation]
    
    for i in range(copy.deepcopy(tessellation['number_of_grains'])):
        
        log.debug('Creating VTK file of grain no. ' + str(i))
        
        ## Creating Folder and assigning the output file path
        output_file_path = Path(store_folder, material_name, now, "vtk_file_individual_grains", "{:s}.vtk".format("grain_" + str(i+1)))
        output_file_path.parent.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.parent.mkdir(exist_ok=True)
        output_file_path.parent.mkdir(exist_ok=True)
        
        number_of_grain_vertices = len(grain_vertices[i])
        number_of_grain_faces = len(grain_faces[i])
    
        with open(str(output_file_path), 'a+') as f:
            f.truncate(0)
            
            ## Writing the Header lines
            f.write("# vtk DataFile Version 3.0 \n")
            f.write("Morphology of individual grains \n")
            f.write("ASCII \n")
            f.write("DATASET UNSTRUCTURED_GRID \n")
            f.write("POINTS " + str(number_of_grain_vertices) + " float \n")
            
            ## Writing all vertices to the file
            for vertex_list in grain_vertices[i]:
                vertex = np.array(vertex_list)
                f.write (str(vertex).strip('[').strip(']') + "\n")
            f.write("\n \n")
            
            ## Writing the Cells (Face vertices) data to file
            sum = 0
            for faces in grain_faces[i]:
                sum += len(faces)
            
            f.write("CELLS " + str(number_of_grain_faces) + " " + str((number_of_grain_faces) + sum) + "\n")
            
            for faces in grain_faces[i]:
                face_indice = np.array(faces)
                f.write(str(len(faces)) + " " + str(face_indice).strip('[').strip(']') + "\n")
            f.write("\n \n")
            
            ## Writing Cell Type data to file
            f.write("CELL_TYPES " + str(number_of_grain_faces) + "\n")
            for i in range(number_of_grain_faces):
                f.write(str(7) + "\n")

    log.info('Successfully created VTK files of all grains individually')