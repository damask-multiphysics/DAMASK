import setuptools
from pathlib import Path
import re

# https://www.python.org/dev/peps/pep-0440
with open(Path(__file__).parent/'damask/VERSION') as f:
    version = re.sub(r'(-([^-]*)).*$',r'.\2',re.sub(r'^v(\d+\.\d+(\.\d+)?)',r'\1',f.readline().strip()))

setuptools.setup(
    name='damask',
    version=version,
    author='The DAMASK team',
    author_email='damask@mpie.de',
    description='DAMASK library',
    long_description='Python library for managing DAMASK simulations',
    url='https://damask.mpie.de',
    packages=setuptools.find_packages(),
    include_package_data=True,
    python_requires = '>=3.6',
    install_requires = [
        'pandas>=0.24',                                                                             # requires numpy
        'scipy>=1.2',
        'h5py>=2.9',                                                                                # requires numpy
        'vtk>=8.1',
        'matplotlib>=3.0',                                                                          # requires numpy, pillow
        'pyaml>=3.12'
    ],
    classifiers = [
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
    ],
)
