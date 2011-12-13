import damask.result

class Result():
  '''
     General class for result parsing.
     Sub-classed by the individual solvers.
  '''
  
  def __init__(self,solver=''):
    solverClass = {
                      'spectral': damask.result.Spectral,
                      'marc':     damask.result.Marc,
                    }
    if solver.lower() in solverClass.keys():
      self.__class__=solverClass[solver.lower()]
      self.__init__()

