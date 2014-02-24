#!/usr/bin/env python
import os
from distutils.core import setup
from setuptools import find_packages

VERSION_FILE_NAME = os.path.join(
    os.path.sep.join(
        os.path.abspath(__file__).split(os.path.sep)[:-1]),
    'VERSION')

with open(VERSION_FILE_NAME, 'r') as version_file:
    VERSION = version_file.readline().strip()

PACKAGES = find_packages()
print PACKAGES

setup(name='Monte Python',
      version=VERSION,
      description='Cosmological Monte Carlo parameter extraction code',
      author='Benjamin Audren',
      author_email='benjamin.audren@epfl.ch',
      url='http://www.montepython.net/',
      packages= PACKAGES,
      )
