# -*- coding: utf-8 -*-
"""
check_libraries.py

Module to ensure installation of all required python libraries.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
created: 19 May 2020
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

import os
from src.set_logger import set_logger as set_logger
name_str = __name__

def check_import(name, install_command):
    """
    Checks if a certain library is installed and installs them if it isn't
    installed.

    Parameter
    ---------
    name: string
        Name of the library
    install_command: string 
        Command-line input required to install a library

    Returns
    -------
    Function returns nothing.  
    """
    try:
        __import__(name)
    except ImportError:
        log.exception('Failed importing ' + name + ' library.')
        input_text = input("Do you want to install " + name + " ? [Yes/No] \n").lower()
        if 'y' in input_text:
            os.system(install_command)
        else:
            log.critical("Since required library " + name +" is not installed, program will exit now. \n")
            exit()

def check_libraries(log_level):
    """
    Iterates through the list of all required libraries and installs them if 
    they aren't installed.

    Parameter
    ---------
    log_level: string
        Logger level to be used.    
    """
    log = set_logger(name_str, 'log_data.log', log_level)

    library_list = ['pip', 'numpy', 'scipy', 'matplotlib', 'numba', 'pytest', 'quaternion', 'tess', 'click', 'gmsh', 'yappi']
    install_command = {'os': ' ', 
                        'pip': 'curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py; python get-pip.py', 
                        'numpy': 'pip install numpy',
                        'scipy': 'pip install scipy',
                        'matplotlib': 'pip install matplotlib',
                        'numba': 'pip install numba',
                        'pytest': 'pip install pytest==4.4.1',
                        'quaternion': 'pip install numpy-quaternion==2019.12.11.22.25.52',
                        'tess': 'pip install tess',
                        'click': 'pip install click==7.0',
                        'gmsh': 'pip install gmsh-sdk==4.5.0-1; echo -e "Please refer http://gmsh.info/ for installation of GMSH software. GMSH software is required for working of meshing functionality. \nFor ubuntu based systems, GMSH can be installed using following command: \n $ sudo apt-get update; sudo apt-get install gmsh"',
                        'yappi': 'pip install yappi==1.2.3'}
    for name in library_list:
        log.info("Checking for " + name + " library.")
        check_import(name, install_command[name])
    
    log.info("All required libraries are installed.")
   
if __name__ == '__main__':
    check_libraries('DEBUG')