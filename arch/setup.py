#from distutils.core import setup
from setuptools import setup, find_packages

print('\nInstalling filament module...\n')

import filament

setup(name='filament',
      description='Description.',
      version='0.1.0',
      author='Team PIE',
      packages=find_packages(),
      requires=['numpy', 'netCDF4']
      )