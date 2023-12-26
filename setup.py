import sys
import os, subprocess, shutil
import platform

from distutils.command.clean import clean as _clean

from setuptools import Extension, setup, Command
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.install import install
from setuptools import find_packages

requires = []

def find_files(dirname, relpath=None):
    def find_paths(dirname):
        items = []
        for fname in os.listdir(dirname):
            path = os.path.join(dirname, fname)
            if os.path.isdir(path):
                items += find_paths(path)
            elif not path.endswith(".py") and not path.endswith(".pyc"):
                items.append(path)
        return items
    items = find_paths(dirname)
    if relpath is None:
        relpath = dirname
    return [os.path.relpath(path, relpath) for path in items]

def get_version():
    """Get the version info from the __version_ file"""

    with open(os.path.join('__version__')) as version_file:
        for line in version_file:
            if line.startswith('__version__'):
                return eval(line.split('=')[-1])

def get_requirements():
    """Get the requirements from the requirements file"""

    with open(os.path.join('requirements.txt')) as requirements_file:
        requirements = requirements_file.readlines()
    return [r.replace('\n', '') for r in requirements]

class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        

VERSION = get_version()
setup_requires = []
install_requires = setup_requires + get_requirements()


setup(
    name = 'scidata',
    version = VERSION,
    description = 'Library to analyze data from simulation of core collapse supernovae done with Aenus-ALCAR.',
    long_description = open('README.md').read(),
    long_description_content_type='text/markdown',
    author = 'Marco Cusinato',
    author_email = 'marco.cusinato@uv.es',
    url = 'https://github.com/MarcoCusinato/scidata',
    download_url = f'https://github.com/MarcoCusinato/scidata',
    keywords = [
        'Aenus-ALCAR',
        'astrophysics',
        'supernovae'
    ],
    setup_requires = setup_requires,
    install_requires = install_requires,
    packages = find_packages(),
    package_data = {
        'scidata.cell': find_files('scidata/cell'),
        'scidata.file_manipulation_functions': find_files('scidata/file_manipulation_functions'),
        'scidata.grid': find_files('scidata/grid'),
        'scidata.GWs_strain': find_files('scidata/GWs_strain'),
        'scidata.legendre': find_files('scidata/legendre'),
        'scidata.math_functions': find_files('scidata/math_functions'),
        'scidata.par': find_files('scidata/par'),
        'scidata.paths': find_files('scidata/paths'),
        'scidata.platform': find_files('scidata/platform'),
        'scidata.quantities': find_files('scidata/quantities'),
        'scidata.spherical_harmonics': find_files('scidata/spherical_harmonics'),
        'scidata.Tools': find_files('scidata/Tools'),
        'scidata.units': find_files('scidata/units')
    },
    python_requires='>=3.6',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Astrophysics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
