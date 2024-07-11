---
title: 'A Python Library for Pre- and Post-Processing of DAMASK Simulations'
tags:
  - Crystal Plasticity
  - Python
authors:
- name: Daniel Otto De Mentock
  affiliation: 1
- name: Sharan Roongta
  orcid: 0000-0001-5754-068X
  affiliation: 1
- name: Franz Roters
  orcid: 0000-0002-9098-9566
  affiliation: 1
- name: Philip Eisenlohr
  orcid: 0000-0002-8220-5995
  affiliation: 2
- name: Martin Diehl
  orcid: 0000-0002-3738-7363
  affiliation: "3, 4"
  corresponding: true

affiliations:
- name: Max-Planck-Institute for Sustainable Materials, 40237 Düsseldorf, Germany
  index: 1
- name: Chemical Engineering and Materials Science, Michigan State University, East Lansing, MI 48824, USA
  index: 2
- name: Department of Materials Engineering, KU Leuven, Kasteelpark Arenberg 44, 3001 Leuven, Belgium
  index: 3
- name: Department of Computer Science, KU Leuven, Celestijnenlaan 200a, 3001 Leuven, Belgium
  index: 4
bibliography: paper.bib

---

# Summary

DAMASK, the Düsseldorf Advanced Material Simulation Kit, is a modular multi-physics crystal plasticity simulation package distributed as free and open-source software under the GNU Affero General Public License (AGPL) [@roters2019damask].
To facilitate easy pre- and post- processing of DAMASK simulations and simplify the creation of custom workflows for Integrated Computational Materials Engineering (ICME), an accompanying Python library has been developed.
This library, which is also AGPL-licensed, is introduced here.

# Statement of need

Multi-physics-enriched crystal plasticity simulations are used to establish processing--structure--property relationships of crystalline materials at engineering length and time scales.
Setting up a simulation requires parameterization of the constitutive models, description of microstructure and geometry, and definition of boundary and initial conditions.
The interpretation of the resulting data requires tools for statistical analysis, plotting, and 3D visualization.
Moreover, the design and study of complex materials using various computational techniques requires interoperability between different software packages in ICME worflows [@ShahEtAl2022].
These needs are best addressed by a modular set of routines for pre- and post-processing that integrate seamlessly into an existing ecosystem.

# Features

The materialpoint models implemented in DAMASK can be used in conjunction with different solvers: `DAMASK_grid`, `DAMASK_mesh`, `DAMASK_Marc`.
The grid solver (`DAMASK_grid`) operates on periodically repeated hexahedral domains discretized by a structured grid, whereas the two other solvers are based on the finite element method and can be used for unstructured meshes and allow for more complex geometries.
Hence, the definition of the geometry as well as boundary conditions and initial conditions depends on the selected solver and requires different pre-processing tools.
In contrast, the configuration of the materialpoint model and the DAMASK-specific HDF5 file format that is used to store the simulation results [@DiehlEtAl2017b] are solver-agnostic.

The DAMASK Python library is a package called `damask` and contains functionality for pre-processing tasks, such as the generation and manipulation of microstructures, load cases, material definition, or numerical parameters, as well as functionality for post-processing that enables analysis and visualization of DAMASK results.
A particular focus is laid on finite-strain continuum mechanics and crystallography.
The routines for conversion between the different kinds of orientation representations, such as Euler angles, rotation matrices, unit quaternions, or axis-angle pairs, are based a consistent set of conventions [@rowenhorst2015consistent].
The provided routines and datas tructures interoperate seamlessly with other libraries from the Python ecosystem, such as `NumPy` [@harris2020array], `pandas` [@mckinney-proc-scipy-2010], `Matplotlib` [@hunter2007matplotlib], `SciPy` [@2020SciPy-NMeth], `VTK`/`PyVista` [@ahrens200536] [@Sullivan_PyVista_3D_plotting_2019], `PyYAML`, `h5py` [@collette_python_hdf5_2014] [@FolkEtAl2011], and `orix` [@johnstone2020density] to facilitate the definition of custom ICME workflows.

Apart from its DAMASK-specific processing capabilities, many routines of the library can be used in other material science and continuum mechanics applications.

## Pre-Processing

The main goal of the pre-processing tools is to provide users with the possibility to incorporate data from various sources into their simulation setup.
For example, the users can integrate domain-specific software such as `Neper` [@quey2011large], `DREAM.3D` [@groeber2014dream], and `Gmsh` [@geuzaine2009gmsh] to define microstructures and employ the provided Python routines to convert them into DAMASK-compatible input files.
Similarly, the materialpoint configuration and load case definition are internally represented as particular Python classes that simplify the creation, modification, and export to YAML files.
In addition to the three mandatory input files, i.e., geometry definition, material configuration, and the load case description, the creation and modification of an optional YAML file to fine-tune numerical parameters is also supported.

## Post-Processing

The post-processing tools are centered around the HDF5 output file resulting from a DAMASK simulation.
This flexible file format is called DADF5, short for "DAMASK HDF5" ensures FAIR data storage [@WilkinsonEtAl2016].
The `Result` class provides custom views on the hierarchical data layout of the DADF5 file and enables the computation of derived quantities that are stored alongside with automatically created metadata in the output file.
The data can be further processed within Python or exported to various file formats for analysis using third-party tools such as `DREAM.3D`, `ParaView`, or `MTEX` [@BachmannEtAl2010].


# Alternatives

Some of the functionality provided by the `damask` package overlaps with other packages, such as `orix` (crystallography and rotations), `PyVista` (3D visualization), and `SciPy` (rotations).
Because the `orix` and `PyVista` packages are currently in an early development status and not available via native package managers on popular Linux distributions, the corresponding `damask` functionality has been implemented independently of both.
The rotation functionality of  the `damask` package does not rely on `SciPy` in order to reproduce one-to-one the conventions used in the main simulation code of DAMASK written in Fortran, despite using `SciPy` routines for other purposes.


# Availability

The `damask` Python package is developed within the DAMASK main repository, but is also available as a separate package via multiple channels.
For documentation and installation options we refer to the [DAMASK website](https://damask.mpie.de).


# Acknowledgements
This research was financially supported by Internal Funds KU Leuven.


# References
