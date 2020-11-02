import setuptools
from pathlib import Path
import re

with open(Path(__file__).parent/'damask/VERSION') as f:
  version = re.sub(r'(-([^-]*)).*$',r'.\2',re.sub(r'^v(\d+\.\d+(\.\d+)?)',r'\1',f.readline().strip()))

setuptools.setup(
    name="damask",
    version=version,
    author="The DAMASK team",
    author_email="damask@mpie.de",
    description="DAMASK library",
    long_description="Python library for pre and post processing of DAMASK simulations",
    url="https://damask.mpie.de",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires = [
        "pandas",                                                                                   # requires numpy
        "scipy",
        "h5py",                                                                                     # requires numpy
        "vtk",
        "matplotlib",                                                                               # requires numpy, pillow
        "pyaml"
    ],
    classifiers = [
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
)
