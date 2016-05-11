#!/usr/bin/env python2
# -*- coding: UTF-8 no BOM -*-

import os,sys
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Generate histogram of N bins in given data range.

""", version = scriptID)

parser.add_option('-d','--data',
                  dest = 'data',
                  type = 'string', metavar = 'string',
                  help = 'column heading for data')
parser.add_option('-w','--weights',
                  dest = 'weights',
                  type = 'string', metavar = 'string',
                  help = 'column heading for weights')
parser.add_option('--range',
                  dest = 'range',
                  type = 'float', nargs = 2, metavar = 'float float',
                  help = 'data range of histogram [min - max]')
parser.add_option('-N',
                  dest = 'N',
                  type = 'int', metavar = 'int',
                  help = 'number of bins')
parser.add_option('--density',
                  dest = 'density',
                  action = 'store_true',
                  help = 'report probability density')
parser.add_option('--logarithmic',
                  dest = 'log',
                  action = 'store_true',
                  help = 'logarithmically spaced bins')
parser.set_defaults(data    = None,
                    weights = None,
                    range   = None,
                    N       = None,
                    density = False,
                    log     = False,
                   )

(options,filenames) = parser.parse_args()

if not options.data: parser.error('no data specified.')
if not options.N:    parser.error('no bin number specified.')

if options.log:
  def forward(x):
    return np.log(x)
  def reverse(x):
    return np.exp(x)
else:
  def forward(x):
    return x
  def reverse(x):
    return x


# --- loop over input files ------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name     = name,
                                    buffered = False,
                                    readonly = True)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors  = []
  remarks = []
  
  if table.label_dimension(options.data) != 1:  errors.append('data {} are not scalar.'.format(options.data))
  if options.weights and \
     table.label_dimension(options.data) != 1:  errors.append('weights {} are not scalar.'.format(options.weights))

  if remarks != []: damask.util.croak(remarks)
  if errors  != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# --------------- read data ----------------------------------------------------------------

  table.data_readArray([options.data,options.weights])

# --------------- auto range ---------------------------------------------------------------

  if options.range is None:
    rangeMin,rangeMax = min(table.data[:,0]),max(table.data[:,0])
  else:
    rangeMin,rangeMax = min(options.range),max(options.range)

# --------------- bin data ----------------------------------------------------------------

  count,edges = np.histogram(table.data[:,0],
                             bins = reverse(forward(rangeMin) + np.arange(options.N+1) *
                                           (forward(rangeMax)-forward(rangeMin))/options.N),
                             range = (rangeMin,rangeMax),
                             weights = None if options.weights is None else table.data[:,1],
                             density = options.density,
                            )
  bincenter = reverse(forward(rangeMin) + (0.5+np.arange(options.N)) *
                     (forward(rangeMax)-forward(rangeMin))/options.N)                             # determine center of bins
  
# ------------------------------------------ assemble header ---------------------------------------
  
  table.info_clear()
  table.info_append([scriptID + '\t' + ' '.join(sys.argv[1:]),
                     scriptID + ':\t' +
                                'data range {} -- {}'.format(rangeMin,rangeMax) +
                                (' (log)' if options.log else ''),
                     ])
  table.labels_clear()
  table.labels_append(['bincenter','count'])
  table.head_write()

# ------------------------------------------ output result -----------------------------------------  

  table.data = np.squeeze(np.dstack((bincenter,count)))
  table.data_writeArray()
  
# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
