# -*- coding: UTF-8 no BOM -*-

# $Id$

import damask.geometry

class Geometry():
  '''
     General class for geometry parsing.
     Sub-classed by the individual solvers.
  '''
  
  def __init__(self,solver=''):
    solverClass = {
                      'spectral': damask.geometry.Spectral,
                      'marc':     damask.geometry.Marc,
                    }
    if solver.lower() in solverClass.keys():
      self.__class__=solverClass[solver.lower()]
      self.__init__()

