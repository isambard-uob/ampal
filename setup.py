"""Setup script for the AMPAL framework."""

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize


def readme():
    """Loads the readme file for AMPAL."""
    with open('README.md', 'r') as inf:
        return inf.read()


setup(name='AMPAL',
      version='1.1.0',
      description='A simple framework for representing biomolecular structure.',
      long_description=readme(),
      long_description_content_type='text/markdown; charset=UTF-8; variant=GFM',
      url='https://github.com/isambard-uob/ampal',
      author='Woolfson Group, University of Bristol',
      author_email='chris.wood@bristol.ac.uk',
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ],
      packages=find_packages('src'),
      package_dir={'': 'src'},
      include_package_data=True,
      # This code automatically builds the Cython extensions.
      ext_modules=cythonize(
          [Extension(
              "ampal.geometry",
              ["src/ampal/geometry.pyx"]),
           ]
      ),
      install_requires=[
          'Cython',
          'networkx',
          'numpy',
          'requests'
      ],
      zip_safe=False,
      )
