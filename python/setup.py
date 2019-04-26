import setuptools
import os

with open(os.path.join(os.path.dirname(__file__),'damask/VERSION')) as f:
  version = f.readline()[1:-1]

setuptools.setup(
    name="damask",
    version=version,
    author="The DAMASK team",
    author_email="damask@mpie.de",
    description="Python library for DAMASK",
    long_description='test',
    #long_description_content_type="text/markdown",
    url="https://damask.mpie.de",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GLPv3",
        "Operating System :: OS Independent",
    ],
)
