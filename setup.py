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

# Recover all needed data files from the likelihoods packages
LIKELIHOODS = [
        elem.replace('.', '/') for elem in PACKAGES
        if elem.find('montepython.likelihoods.') != -1]
DATA_FILES = [('', ['VERSION']),
              ('', ['default.conf.template'])]
for likelihood in LIKELIHOODS:
    data = [
        os.path.join(likelihood, elem) for elem in os.listdir(likelihood) if
        elem.find('.py') == -1]
    DATA_FILES.append((likelihood, data))

setup(name='montepython',
      version=VERSION,
      description='Cosmological Monte Carlo parameter extraction code',
      author='Benjamin Audren',
      author_email='benjamin.audren@epfl.ch',
      url='http://www.montepython.net/',
      packages=PACKAGES,
      install_requires=['argparse'],
      data_files=DATA_FILES,
      )
