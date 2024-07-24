"""
polyhedral_analysis
"""

from setuptools import setup, find_packages
from polyhedral_analysis import __version__ as VERSION

readme = 'README.md'
long_description = open( readme ).read()

setup(
    name='polyhedral_analysis',
    version=VERSION,
    description='A library for analysis of coordination polyhedra from molecular dynamics trajectories',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Benjamin J. Morgan',
    author_email='bjm42@bath.ac.uk',
    url='https://github.com/bjmorgan/polyhedral-analysis', 
    download_url="https://github.com/bjmorgan/polyhedral-analysis/archive/%s.tar.gz" % (VERSION),
    packages=find_packages( exclude=['docs', 'tests*'] ),
    license='MIT',
    python_requires='>=3.7',
    install_requires=['numpy',
                      'pymatgen>=2024.7.18',
                      'scipy', 
                      'coverage>=6.3',
                      'vg',
                      'monty',
                      'tqdm']
    )
