"""Setup script for the AMPAL framework."""

from setuptools import setup

def readme():
    """Loads the readme file for AMPAL."""
    with open('README.md', 'r') as inf:
        return inf.read()

setup(name='AMPAL',
      version='0.1',
      description='A simple framework for representing biomolecular structure.',
      long_description=readme(),
      url='https://github.com/isambard-uob/ampal',
      author='Woolfson Group, University of Bristol',
      license='MIT',
      packages=['ampal'],
      zip_safe=False,
     )
