#!/usr/bin/env python

import os,re,sys,math,numpy,string
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



def Mises(what,tensor):

	dev = tensor - numpy.trace(tensor)/3.0*numpy.eye(3)
	symdev = 0.5*(dev+dev.T)
	return math.sqrt(numpy.sum(symdev*symdev.T)*
 	       {
	        'stress': 3.0/2.0,
	        'strain': 2.0/3.0,
	       }[what])
	

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add vonMises equivalent values for symmetric part of requested strains and/or stresses.

""" + string.replace('$Id: addMises 966 2011-08-18 08:00:19Z MPIE\p.eisenlohr $','\n','\\n')
)


parser.add_option('-m','--memory',      dest='memory', action='store_true', \
                                        help='load complete file into memory [%default]')
parser.add_option('-e','--strain',      dest='strain', action='extend', type='string', \
                                        help='heading(s) of columns containing strain tensors')
parser.add_option('-s','--stress',      dest='stress', action='extend', type='string', \
                                        help='heading(s) of columns containing stress tensors')

parser.set_defaults(memory = False)
parser.set_defaults(strain = [])
parser.set_defaults(stress = [])

(options,filenames) = parser.parse_args()

if len(options.strain) + len(options.stress) == 0:
  parser.error('no data column specified...')

datainfo = {                                                               # list of requested labels per datatype
             'strain':     {'len':9,
                            'label':[]},
             'stress':     {'len':9,
                            'label':[]},
           }


if options.strain != None:    datainfo['strain']['label'] += options.strain
if options.stress != None:    datainfo['stress']['label'] += options.stress


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name, 'input':open(name), 'output':open(name+'_tmp','w')})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  print file['name']

  #  get labels by either read the first row, or - if keyword header is present - the last line of the header

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    passOn  = [file['input'].readline() for i in range(1,headerlines)]
    headers = file['input'].readline().split()
  else:
    headerlines = 1
    passOn  = []
    headers = firstline.split()

  if options.memory:
    data = file['input'].readlines()
  else:
    data = []

  for i,l in enumerate(headers):
    if l.startswith('1_'):
      if re.match('\d+_',l[2:]) or i == len(headers)-1 or not headers[i+1].endswith(l[2:]):
        headers[i] = l[2:]

  active = {}
  column = {}
  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = {True :'1_%s',
             False:'%s'   }[info['len']>1]%label
      if key not in headers:
        sys.stderr.write('column %s not found...\n'%key)
      else:
        if datatype not in active: active[datatype] = []
        if datatype not in column: column[datatype] = {}
        active[datatype].append(label)
        column[datatype][label] = headers.index(key)
        head.append('Mises(%s)'%label)

# ------------------------------------------ assemble header ---------------------------------------  

  output = '%i\theader'%(headerlines+1) + '\n' + \
           ''.join(passOn) + \
           string.replace('$Id: addMises 966 2011-08-18 08:00:19Z MPIE\p.eisenlohr $','\n','\\n')+ '\t' + \
           ' '.join(sys.argv[1:]) + '\n' + \
           '\t'.join(headers + head) + '\n'                              # build extended header

  if not options.memory:
    file['output'].write(output)
    output = ''

# ------------------------------------------ read file ---------------------------------------  

  for line in {True  : data,
               False : file['input']}[options.memory]:
    items = line.split()[:len(headers)]
    if len(items) < len(headers):
      continue
  
    output += '\t'.join(items)

    for datatype,labels in active.items():
      for label in labels:
        theMises = Mises(datatype,
                         numpy.array(map(float,items[column[datatype][label]:
                                                     column[datatype][label]+datainfo[datatype]['len']]),'d').reshape(3,3))
        output += '\t%f'%theMises

    output += '\n'

    if not options.memory:
      file['output'].write(output)
      output = ''

  file['input'].close()
  
# ------------------------------------------ output result ---------------------------------------  

  if options.memory:
    file['output'].write(output)

  if file['name'] != 'STDIN':
    file['output'].close
    os.rename(file['name']+'_tmp',file['name'])
