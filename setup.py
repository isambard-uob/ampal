"""Setup script for the AMPAL framework."""

from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize


def readme():
    """Loads the readme file for AMPAL."""
    with open('README.md', 'r') as inf:
        return inf.read()


setup(name='AMPAL',
      version='1.0.0',
      description='A simple framework for representing biomolecular structure.',
      long_description=readme(),
      url='https://github.com/isambard-uob/ampal',
      author='Woolfson Group, University of Bristol',
      license='MIT',
      packages=['ampal'],
      # This code automatically builds the Cython extensions.
      ext_modules=cythonize(
          [Extension(
              "ampal.geometry",
              ["ampal/geometry.pyx"]),
           ]
      ),
      install_requires=[
          'Cython'
      ],
      zip_safe=False,
      )
