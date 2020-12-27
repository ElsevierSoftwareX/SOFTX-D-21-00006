# -*- coding: utf-8 -*-
"""
main_import_statements.py

Module to import all required libraries for other modules.

This was created as part of "Personal Programming Project (PPP)" coursework in 
"Computation Materials Science (CMS)" M. Sc program at TU Bergakademie Freiberg,
Germany.

@authors: Serrao Prince Henry, Arun Prakash
@email: 
created: 20 December 2019
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

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from tess import Container
import inspect
from src.__version__ import __version__ as version
from datetime import datetime
import quaternion
import numba
import click
import copy
import scipy
from scipy.stats import lognorm
from scipy.optimize import minimize
import matplotlib.animation as animation
from time import sleep
import os
import importlib.util
from scipy.spatial import ConvexHull
import pytest
import yappi
import shutil
import subprocess
import logging

