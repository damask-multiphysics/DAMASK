#!/usr/bin/env python2.7
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
Add data of selected column(s) from (first) row of linked ASCIItable that shares the linking column value.

""", version = scriptID)

parser.add_option('--link',
                  dest = 'link', nargs = 2,
                  type = 'string', metavar = 'string string',
                  help = 'column labels containing linked values')
parser.add_option('-l','--label',
                  dest = 'label',
                  action = 'extend', metavar = '<string LIST>',
                  help = 'column label(s) to add from linked ASCIItable')
parser.add_option('-a','--asciitable',
                  dest = 'asciitable',
                  type = 'string', metavar = 'string',
                  help = 'linked ASCIItable')

parser.set_defaults()

(options,filenames) = parser.parse_args()

if options.label is None:
  parser.error('no data columns specified.')
if options.link is None:
  parser.error('no linking columns given.')

# -------------------------------------- process linked ASCIItable --------------------------------

if options.asciitable is not None and os.path.isfile(options.asciitable):

  linkedTable = damask.ASCIItable(name = options.asciitable,
                                  buffered = False,
                                  readonly = True) 
  linkedTable.head_read()                                                                           # read ASCII header info of linked table
  linkDim = linkedTable.label_dimension(options.link[1])                                            # dimension of linking column

  missing_labels = linkedTable.data_readArray([options.link[1]]+options.label)                      # try reading linked ASCII table
  linkedTable.close()                                                                               # close linked ASCII table

  if len(missing_labels) > 0:
    damask.util.croak('column{} {} not found...'.format('s' if len(missing_labels) > 1 else '',', '.join(missing_labels)))

  index = linkedTable.data[:,:linkDim]
  data  = linkedTable.data[:,linkDim:]
else:
  parser.error('no linked ASCIItable given.')

# --- loop over input files -----------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,"{} {} <== {} {}".format(name,damask.util.deemph('@ '+options.link[0]),
                                                         options.asciitable,damask.util.deemph('@ '+options.link[1])))

# ------------------------------------------ read header ------------------------------------------

  table.head_read()

# ------------------------------------------ sanity checks ----------------------------------------

  errors = []

  myLink    = table.label_index    (options.link[0])  
  myLinkDim = table.label_dimension(options.link[0])  
  if myLink <  0:          errors.append('linking column {} not found.'.format(options.link[0]))
  if myLinkDim != linkDim: errors.append('dimension mismatch for column {}.'.format(options.link[0]))

  if errors != []:
    damask.util.croak(errors)
    table.close(dismiss = True)
    continue

# ------------------------------------------ assemble header --------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_append(linkedTable.labels(raw = True)[linkDim:])                                     # extend with new labels (except for linked column)
  
  table.head_write()

# ------------------------------------------ process data ------------------------------------------

  outputAlive = True
  while outputAlive and table.data_read():                                                          # read next data line of ASCII table
    try:
      table.data_append(data[np.argwhere(np.all((map(float,table.data[myLink:myLink+myLinkDim]) - index)==0,axis=1))[0]])  # add data of first matching line
    except IndexError:
      table.data_append(np.nan*np.ones_like(data[0]))                                               # or add NaNs
    outputAlive = table.data_write()                                                                # output processed line

# ------------------------------------------ output finalization -----------------------------------  

  table.close()                                                                                     # close ASCII tables
