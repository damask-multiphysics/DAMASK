import damask.solver

class Solver():
  '''
     General class for solver specific functionality.
     Sub-classed by the individual solvers.
  '''
  
  def __init__(self,solver=''):
    solverClass = {
                      'spectral': damask.solver.Spectral,
                      'marc':     damask.solver.Marc,
                    }
    if solver.lower() in solverClass.keys():
      self.__class__=solverClass[solver.lower()]
      self.__init__()

