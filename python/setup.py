import setuptools
import os

with open(os.path.join(os.path.dirname(__file__),'damask/VERSION')) as f:
  version = f.readline()

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
        "numpy",
        "scipy",
        "pandas",
        "h5py",
        "vtk",
    ],
    classifiers = [
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
)
