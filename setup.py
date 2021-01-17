#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from setuptools import setup, find_packages
import ast
import re

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__)))

_version_re = re.compile(r'__version__\s+=\s+(.*)')
with open(os.path.join(str(root_dir), 'optimic', 'src', '__version__.py'), 'r') as f:
    version = str(ast.literal_eval(_version_re.search(f.read()).group(1)))

with open(os.path.join(str(root_dir), "README.md"), "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='optimic',
    version=version,
    packages=find_packages(),
    #package_dir={'':'optimic'},
    #package_data = {'': ['images/*', 'test_cases/*', 'util/*']},
    url='https://gitlab.com/arun.prakash.mimm/optimic.git',
    license='GPLv3',
    author='P. H. Serrao, A. Prakash',
    author_email='prince.serrao.code@gmail.com, arun.prakash@imfd.tu-freiberg.de',
    description='A tool to generate optimized polycrystalline microstructures for materials simulations.',
    long_description=long_description,
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'numpy==1.17.4',
        'scipy==1.4.1',
        'matplotlib==3.0.3',
        'llvmlite==0.32.1',
        'numba==0.46.0',
        'pytest==4.4.1',
        'numpy-quaternion==2019.12.11.22.25.52',
        'tess',
        'click==7.0',
        'gmsh-sdk==4.5.0-1',
        'yappi'
    ],
    entry_points={
        'console_scripts': [
            'optimic = optimic.src.main:guide',
        ]
    })