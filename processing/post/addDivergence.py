#!/usr/bin/python

import os,re,sys,math,string,numpy,DAMASK
from optparse import OptionParser, Option

# -----------------------------
class extendableOption(Option):
# -----------------------------
# used for definition of new option parser action 'extend', which enables to take multiple option arguments
# taken from online tutorial http://docs.python.org/library/optparse.html
  
  ACTIONS = Option.ACTIONS + ("extend",)
  STORE_ACTIONS = Option.STORE_ACTIONS + ("extend",)
  TYPED_ACTIONS = Option.TYPED_ACTIONS + ("extend",)
  ALWAYS_TYPED_ACTIONS = Option.ALWAYS_TYPED_ACTIONS + ("extend",)

  def take_action(self, action, dest, opt, value, values, parser):
    if action == "extend":
      lvalue = value.split(",")
      values.ensure_value(dest, []).extend(lvalue)
    else:
      Option.take_action(self, action, dest, opt, value, values, parser)

def location(idx,res):
  return ( idx  % res[0], \
         (idx // res[0]) % res[1], \
         (idx // res[0] // res[1]) % res[2] )

def index(location,res):
  return ( location[0] % res[0]                   + \
         (location[1] % res[1]) * res[0]          + \
         (location[2] % res[2]) * res[0] * res[1]   )

def prefixMultiply(what,len):
  return {True: ['%i_%s'%(i+1,what) for i in range(len)],
         False:[what]}[len>1]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

FDcoefficients = [ \
                   [1.0/2.0,      0.0,      0.0,       0.0],\
                   [2.0/3.0,-1.0/12.0,      0.0,       0.0],\
                   [3.0/4.0,-3.0/20.0,1.0/ 60.0,       0.0],\
                   [4.0/5.0,-1.0/ 5.0,4.0/105.0,-1.0/280.0],\
                 ]
parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing divergence of requested column(s).
Operates on periodic ordered three-dimensional data sets.
Deals with both vector- and tensor-valued fields.

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('--fdm',              dest='accuracy', action='extend', type='string', \
                                        help='degree of central difference accuracy')
parser.add_option('--fft',              dest='fft', action='store_true', \
                                        help='calculate divergence in Fourier space [%default]')
parser.add_option('-m','--memory',      dest='memory', action='store_true', \
                                        help='memory efficient calculation (not possible for FFT based divergency [%default]')
parser.add_option('-v','--vector',      dest='vector', action='extend', type='string', \
                                        help='heading of columns containing vector field values')
parser.add_option('-t','--tensor',      dest='tensor', action='extend', type='string', \
                                        help='heading of columns containing tensor field values')
parser.add_option('-d','--dimension',   dest='dim', type='float', nargs=3, \
                                        help='physical dimension of data set in x (fast) y z (slow) [%default]')
parser.add_option('-r','--resolution',  dest='res', type='int', nargs=3, \
                                        help='resolution of data set in x (fast) y z (slow)')
parser.add_option('-s','--skip',        dest='skip', type='int', nargs=3, \
                                        help='items skipped due to periodicity in x (fast) y z (slow)')

parser.set_defaults(accuracy = [])
parser.set_defaults(memory = False)
parser.set_defaults(fft = False)
parser.set_defaults(vector = [])
parser.set_defaults(tensor = [])
parser.set_defaults(dim = [])
parser.set_defaults(skip = [0,0,0])
accuracyChoices = [2,4,6,8]

(options,filenames) = parser.parse_args()

if len(options.vector) + len(options.tensor) == 0:
  parser.error('no data column specified...')
if len(options.dim) < 3:
  parser.error('improper dimension specification...')
if not options.res or len(options.res) < 3:
  parser.error('improper resolution specification...')
for choice in options.accuracy:
  if int(choice) not in accuracyChoices:
    parser.error('accuracy must be chosen from %s...'%(', '.join(accuracyChoices)))
if options.fft: options.accuracy.append('fft')

resSkip = [lambda a,b: a+b for a,b in zip(options.res,options.skip)]
datainfo = {                                                               # list of requested labels per datatype
             'vector':     {'len':3,
                            'label':[]},
             'tensor':     {'len':9,
                            'label':[]},
           }

if options.vector != None:    datainfo['vector']['label'] += options.vector
if options.tensor != None:    datainfo['tensor']['label'] += options.tensor

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'handle':sys.stdin})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'handle':open(name)})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  print file['name']

  content = file['handle'].readlines()
  file['handle'].close()
  
  #  get labels by either read the first row, or - if keyword header is present - the last line of the header

  headerlines = 1
  m = re.search('(\d+)\s*head', content[0].lower())
  if m:
    headerlines = int(m.group(1))
  passOn  = content[1:headerlines]
  headers = content[headerlines].split()
  data    = content[headerlines+1:]
    
  regexp = re.compile('1_\d+_')
  for i,l in enumerate(headers):
    if regexp.match(l):
      headers[i] = l[2:]

  active = {}
  column = {}
  values = {}
  div_field ={}
  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = {True :'1_%s',
             False:'%s'   }[info['len']>1]%label
      if key not in headers:
        print 'column %s not found...'%key
      else:
        if datatype not in active: active[datatype] = []
        if datatype not in column: column[datatype] = {}
        if datatype not in values: values[datatype] = {}
        if datatype not in div_field: div_field[datatype] = {}
        active[datatype].append(label)
        column[datatype][label] = headers.index(key)
        values[datatype][label] = numpy.array([0.0 for i in xrange(datainfo[datatype]['len']*\
                                          options.res[0]*options.res[1]*options.res[2])]).\
                                          reshape((options.res[0],options.res[1],options.res[2],\
                                                     3,datainfo[datatype]['len']//3))

        for what in options.accuracy:                                                               # loop over all requested degrees of accuracy (plus potentially fft)
          if not options.memory or what != 'fft':                                                   # FFT divergence excluded in memory saving mode
            head += prefixMultiply('div%s(%s)'%(what,label),datainfo[datatype]['len']//3)
        
# ------------------------------------------ assemble header ---------------------------------------  

  output = '%i\theader'%(headerlines+1) + '\n' + \
           ''.join(passOn)                 + \
           string.replace('$Id$','\n','\\n')+ '\t' + \
           ' '.join(sys.argv[1:]) + '\n' + \
           '\t'.join(headers + head) + '\n'                              # build extended header

# ------------------------------------------ read value field ---------------------------------------  

  idx = 0
  for line in data:
    items = line.split()[:len(headers)]
    if len(items) < len(headers):                             # skip too short lines (probably comments or invalid)
      continue
    locSkip = location(idx,resSkip)
    if (    locSkip[0] < options.res[0]
        and locSkip[1] < options.res[1]
        and locSkip[2] < options.res[2] ):                    # only take values that are not periodic images
      for datatype,labels in active.items():
        for label in labels:
          values[datatype][label][locSkip[0]][locSkip[1]][locSkip[2]]\
              = numpy.reshape(items[column[datatype][label]:
                                    column[datatype][label]+datainfo[datatype]['len']],(3,datainfo[datatype]['len']//3))
    idx += 1
# ------------------------------------------ read file --------------------------------------- 
  if options.memory:
    idx = 0
    for line in data:
      items = line.split()[:len(headers)]
      if len(items) < len(headers):
        continue
    
      output += '\t'.join(items)
      (x,y,z) = location(idx,options.res)
    
      for datatype,labels in active.items():
        for label in labels:
          for accuracy in options.accuracy:
            if accuracy == 'fft': continue
            for k in range(datainfo[datatype]['len']/3):                                  # formulas from Mikhail Itskov: Tensor Algebra and Tensor Analysis for Engineers, Springer 2009, p 52
              theDiv = 0.0
              for a in range(int(accuracy)//2):

                theDiv += FDcoefficients[int(accuracy)//2-1][a] * \
                          ( \
                           (values[datatype][label][location(index([x+1+a,y,z],options.res),options.res)[0]] \
                                                   [location(index([x+1+a,y,z],options.res),options.res)[1]] \
                                                   [location(index([x+1+a,y,z],options.res),options.res)[2]][k][0] - \
                            values[datatype][label][location(index([x-1-a,y,z],options.res),options.res)[0]] \
                                                   [location(index([x-1-a,y,z],options.res),options.res)[1]] \
                                                   [location(index([x-1-a,y,z],options.res),options.res)[2]][k][0]) * options.res[0] / options.dim[0] + \
                           (values[datatype][label][location(index([x,y+1+a,z],options.res),options.res)[0]] \
                                                   [location(index([x,y+1+a,z],options.res),options.res)[1]] \
                                                   [location(index([x,y+1+a,z],options.res),options.res)[2]][k][1] - \
                            values[datatype][label][location(index([x,y-1-a,z],options.res),options.res)[0]] \
                                                   [location(index([x,y-1-a,z],options.res),options.res)[1]] \
                                                   [location(index([x,y-1-a,z],options.res),options.res)[2]][k][1]) * options.res[1] / options.dim[1] + \
                           (values[datatype][label][location(index([x,y,z+1+a],options.res),options.res)[0]] \
                                                   [location(index([x,y,z+1+a],options.res),options.res)[1]] \
                                                   [location(index([x,y,z+1+a],options.res),options.res)[2]][k][2]- \
                            values[datatype][label][location(index([x,y,z-1-a],options.res),options.res)[0]] \
                                                   [location(index([x,y,z-1-a],options.res),options.res)[1]] \
                                                   [location(index([x,y,z-1-a],options.res),options.res)[2]][k][2]) * options.res[2] / options.dim[2] \
                          )
              output += '\t%f'%theDiv
      
      output += '\n'
      idx += 1
  
  else:
    for datatype,labels in active.items():
      for label in labels:
        if label not in div_field[datatype]: div_field[datatype][label] = {}
        for accuracy in options.accuracy:
          div_field[datatype][label][accuracy] = numpy.array([0.0 for i in range((datainfo[datatype]['len'])//3*\
                                                                                  options.res[0]*options.res[1]*options.res[2])]).\
                                                                                  reshape((options.res[0],options.res[1],options.res[2],\
                                                                                  datainfo[datatype]['len']//3))
          if accuracy == 'fft':
            div_field[datatype][label][accuracy] = DAMASK.math.divergence_fft(options.res,options.dim,datainfo[datatype]['len']//3,values[datatype][label])
          else:
            div_field[datatype][label][accuracy] = DAMASK.math.divergence_fdm(options.res,options.dim,datainfo[datatype]['len']//3,eval(accuracy)//2-1,values[datatype][label])

    idx = 0
    for line in data:
      items = line.split()[:len(headers)]
      if len(items) < len(headers):
        continue
    
      output += '\t'.join(items)
      
    
      for datatype,labels in active.items():
        for label in labels:
          for accuracy in options.accuracy:
            for i in range(datainfo[datatype]['len']/3):    
              output += '\t%f'%div_field[datatype][label][accuracy][location(idx,options.res)[0]][location(idx,options.res)[1]][location(idx,options.res)[2]][i]
      output += '\n'
      idx += 1
 
  
# ------------------------------------------ output result ---------------------------------------  

  if file['name'] == 'STDIN':
    print output
  else:
    file['handle'] = open(file['name']+'_tmp','w')
    try:
      file['handle'].write(output)
      file['handle'].close()
      os.rename(file['name']+'_tmp',file['name'])
    except:
      print 'error during writing',file['name']+'_tmp'
