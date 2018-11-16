#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,re,sys,fnmatch
import math                                                                                         # noqa
import numpy as np
from optparse import OptionParser
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def sortingList(labels,whitelistitems):

  indices = []
  names   = []

  for label in labels:
    if re.match('^\d+_',label):
      indices.append(int(label.split('_',1)[0]))
      names.append(label.split('_',1)[1])
    else:
      indices.append(0)
      names.append(label)
      
  return [indices,names,whitelistitems]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Filter rows according to condition and columns by either white or black listing.

Examples:
Every odd row if x coordinate is positive -- " #ip.x# >= 0.0 and #_row_#%2 == 1 ).
All rows where label 'foo' equals 'bar' -- " #s#foo# == 'bar' "

""", version = scriptID)

parser.add_option('-w','--white',
                  dest   = 'whitelist',
                  action = 'extend', metavar = '<string LIST>',
                  help   = 'whitelist of column labels (a,b,c,...)')
parser.add_option('-b','--black',
                  dest   = 'blacklist',
                  action = 'extend', metavar='<string LIST>',
                  help   = 'blacklist of column labels (a,b,c,...)')
parser.add_option('-c','--condition',
                  dest   = 'condition', metavar='string',
                  help   = 'condition to filter rows')

parser.set_defaults(condition = None,
                   )

(options,filenames) = parser.parse_args()

# --- loop over input files -------------------------------------------------------------------------

if filenames == []: filenames = [None]

for name in filenames:
  try:    table = damask.ASCIItable(name = name,
                                    buffered = False)
  except: continue
  damask.util.report(scriptName,name)

# ------------------------------------------ assemble info ---------------------------------------  

  table.head_read()

# ------------------------------------------ process data ---------------------------------------  

  specials = { \
               '_row_': 0,
             }
  labels = []
  positions = []

  for position,label in enumerate(table.labels(raw = True)):
    if    (options.whitelist is None or     any([   position in table.label_indexrange(needle) \
                                                 or fnmatch.fnmatch(label,needle) for needle in options.whitelist])) \
      and (options.blacklist is None or not any([   position in table.label_indexrange(needle) \
                                                 or fnmatch.fnmatch(label,needle) for needle in options.blacklist])):  # a label to keep?
      labels.append(label)                                                                          # remember name...
      positions.append(position)                                                                    # ...and position

  if len(labels) > 0 and options.whitelist is not None and options.blacklist is None:               # check whether reordering is possible
    whitelistitem = np.zeros(len(labels),dtype=int)
    for i,label in enumerate(labels):                                                               # check each selected label
      match = [   positions[i] in table.label_indexrange(needle) \
               or fnmatch.fnmatch(label,needle) for needle in options.whitelist]                    # which whitelist items do match it
      whitelistitem[i] = match.index(True) if np.sum(match) == 1 else -1                            # unique match to a whitelist item --> store which

    order =      range(len(labels)) if np.any(whitelistitem < 0) \
            else np.lexsort(sortingList(labels,whitelistitem))                                      # reorder if unique, i.e. no "-1" in whitelistitem
  else:
    order = range(len(labels))                                                                      # maintain original order of labels
  
# --------------------------------------- evaluate condition ---------------------------------------
  if options.condition is not None:
    condition = options.condition                                                                   # copy per file, since might be altered inline
    breaker = False
  
    for position,(all,marker,column) in enumerate(set(re.findall(r'#(([s]#)?(.+?))#',condition))):              # find three groups
      idx = table.label_index(column)
      dim = table.label_dimension(column)
      if idx < 0 and column not in specials:
        damask.util.croak('column "{}" not found.'.format(column))
        breaker = True
      else:
        if column in specials:
          replacement = 'specials["{}"]'.format(column)
        elif dim == 1:                                                                                # scalar input
          replacement = '{}(table.data[{}])'.format({  '':'float',
                                                        's#':'str'}[marker],idx)                      # take float or string value of data column
        elif dim > 1:                                                                                 # multidimensional input (vector, tensor, etc.)
          replacement = 'np.array(table.data[{}:{}],dtype=float)'.format(idx,idx+dim)                 # use (flat) array representation
       
        condition = condition.replace('#'+all+'#',replacement)
    
    if breaker: continue                                                                              # found mistake in condition evaluation --> next file
  
# ------------------------------------------ assemble header ---------------------------------------

  table.info_append(scriptID + '\t' + ' '.join(sys.argv[1:]))
  table.labels_clear()
  table.labels_append(np.array(labels)[order])                                                        # update with new label set
  table.head_write()

# ------------------------------------------ process and output data ------------------------------------------

  positions = np.array(positions)[order]
  
  atOnce = options.condition is None
  if atOnce:                                                                                          # read full array and filter columns
    try:
      table.data_readArray(positions+1)                                                               # read desired columns (indexed 1,...)
      table.data_writeArray()                                                                         # directly write out
    except:
      table.data_rewind()
      atOnce = False                                                                                  # data contains items that prevent array chunking

  if not atOnce:                                                                                      # read data line by line
    outputAlive = True
    while outputAlive and table.data_read():                                                          # read next data line of ASCII table
      specials['_row_'] += 1                                                                          # count row
      if options.condition is None or eval(condition):                                                # valid row ?
        table.data = [table.data[position] for position in positions]                                 # retain filtered columns
        outputAlive = table.data_write()                                                              # output processed line

# ------------------------------------------ finalize output -----------------------------------------

  table.close()                                                                                     # close input ASCII table (works for stdin)
