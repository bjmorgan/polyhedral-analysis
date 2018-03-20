import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

readme = 'README.md'
try:
    import pypandoc
    long_description = pypandoc.convert( readme, 'rst')
except ImportError:
    long_description = open( readme ).read()

from polyhedral_analysis import __version__ as VERSION

config = {
    'name': 'polyhedral_analysis',
    'description': 'TODO',
    'long_description': long_description,
    'author': 'Benjamin J. Morgan',
    'author_email':'bjm42@bath.ac.uk',
    'url': 'https://github.com/bjmorgan/polyhedral-analysis', 
    'download_url': "https://github.com/bjmorgan/polyhedral-analysis/archive/%s.tar.gz" % (VERSION),
    'version': VERSION,
    'packages': ['polyhedral_analysis'],
    'license': 'MIT',
    'install_requires': [ 'numpy',
                          'pymatgen',
                          'scipy', 
                          'coverage==4.3.4',
                          'codeclimate-test-reporter' ]
}

setup(**config)
