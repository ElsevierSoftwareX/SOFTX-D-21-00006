# -*- coding: utf-8 -*-
"""
set_logger.py

Module to setup and configure Python logger.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

For reporting bugs/issues: <https://gitlab.com/arun.prakash.mimm/optimic>

@authors: Serrao Prince Henry, Arun Prakash
@email: prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de
created: 21 May 2020
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

import logging
from logging.handlers import RotatingFileHandler

def set_logger(name, file_name, level_name):
    """
    Instantiate logger module.

    Parameter
    ---------
    name: string
        Name of the current module.

    file_name: string
        Filename into which log data is to be stored.

    level_name: string
        Logger level to be used.

    Returns
    -------
    Logger object.

    """
    log = logging.getLogger(name)
    
    level_value = logging.getLevelName(level_name)
    log.setLevel(level_value)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    log_file = RotatingFileHandler(file_name, maxBytes=5*1024*1024, backupCount=3, encoding=None)
    log_file.setFormatter(formatter)

    print_handler = logging.StreamHandler()
    print_handler.setFormatter(formatter)

    if not len(log.handlers):
        log.addHandler(log_file)
        log.addHandler(print_handler)

    return log