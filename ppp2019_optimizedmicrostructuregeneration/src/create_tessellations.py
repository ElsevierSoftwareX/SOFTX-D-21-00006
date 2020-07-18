# -*- coding: utf-8 -*-
"""
create_tessellations.py

Module to perform voronoi tessellations using 'Tess' package

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

def create_tessellations(seed_array, limit, log_level):                            #x are seed coordinates, y are limits
    """
    Input: The array of all the seed coordinates and the array of the limits of the simulation box.
    Process: The function calls  the 'tess' package with the arrays of seed arrays and the limits as parameters.
            The periodicity is set to be True.
    Output: The function returns the cells as a list after performing the tessellations.
    """
    log = set_logger(__name__, 'log_data.log', log_level)
    log.debug('Tessellations are to be created with limits as ' + str(limit) + ' and seed array as:\n' + str(seed_array))

    local_tessellation = Container(seed_array, limits=limit, periodic = True)
    
    log.info('Tessellations were created successfully.')

    number_of_grains = len(local_tessellation)
    number_of_faces_list = [v.number_of_faces() for v in local_tessellation]
    vertices_list = [v.vertices() for v in local_tessellation]
    face_vertices_list = [v.face_vertices() for v in local_tessellation]
    centroid_list = [v.centroid() for v in local_tessellation]
    volume_list = [v.volume() for v in local_tessellation]
    normals_list = [v.normals() for v in local_tessellation]
    neighbors_list = [v.neighbors() for v in local_tessellation]
    face_area_list = [v.face_areas() for v in local_tessellation]
    number_of_edges_list = [v.number_of_edges() for v in local_tessellation]

    tessellation = {}
    tessellation['number_of_grains'] = number_of_grains
    tessellation['number_of_faces_list'] = number_of_faces_list
    tessellation['vertices_list'] = vertices_list
    tessellation['face_vertices_list'] = face_vertices_list
    tessellation['centroid_list'] = centroid_list
    tessellation['volume_list'] = volume_list
    tessellation['normals_list'] = normals_list
    tessellation['neighbors_list'] = neighbors_list
    tessellation['face_area_list'] = face_area_list
    tessellation['number_of_edges_list'] = number_of_edges_list
    
    log.info('All required properties of created tessellations are stored in the dictionary successfully.')

    return tessellation                                                   # returns list of cells as class instances