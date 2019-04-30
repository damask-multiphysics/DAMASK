import setuptools
import os

with open(os.path.join(os.path.dirname(__file__),'damask/VERSION')) as f:
  version = f.readline()[1:-1]

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
        "scipy",
        "h5py",
        "vtk"
    ],
    license = 'GPL3',
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL3",
        "Operating System :: OS Independent",
    ],
)
