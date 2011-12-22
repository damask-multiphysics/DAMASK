# $Id$

from .geometry import Geometry

class Spectral(Geometry):

  def __init__(self):
    self.solver='Spectral'
