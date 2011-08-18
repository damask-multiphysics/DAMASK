#!/usr/bin/python

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


def prefixMultiply(what,len):

  return {True: ['%i_%s'%(i+1,what) for i in range(len)],
          False:[what]}[len>1]




# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing Cauchy stress based on given column(s) of
deformation gradient and first Piola--Kirchhoff stress.

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('-m','--memory',      dest='memory', action='store_true', \
                                        help='load complete file into memory [%default]')
parser.add_option('-f','--defgrad',     dest='defgrad', type='string', \
                                        help='heading of columns containing deformation gradient [%default]')
parser.add_option('-p','--stress',      dest='stress', type='string', \
                                        help='heading of columns containing first Piola--Kirchhoff stress [%default]')

parser.set_defaults(memory = False)
parser.set_defaults(defgrad = 'f')
parser.set_defaults(stress = 'p')

(options,filenames) = parser.parse_args()

if options.defgrad == None or options.stress == None:
  parser.error('missing data column...')

datainfo = {                                                               # list of requested labels per datatype
             'defgrad':    {'len':9,
                            'label':[]},
             'stress':     {'len':9,
                            'label':[]},
           }


datainfo['defgrad']['label'].append(options.defgrad)
datainfo['stress']['label'].append(options.stress)


# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'handle':sys.stdin})
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
  head += prefixMultiply('Cauchy',datainfo[datatype]['len'])

# ------------------------------------------ assemble header ---------------------------------------  

  output = '%i\theader'%(headerlines+1) + '\n' + \
           ''.join(passOn) + \
           string.replace('$Id$','\n','\\n')+ '\t' + \
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

    F = numpy.array(map(float,items[column['defgrad'][active['defgrad'][0]]:
                                    column['defgrad'][active['defgrad'][0]]+datainfo['defgrad']['len']]),'d').reshape(3,3)
    P = numpy.array(map(float,items[column['stress'][active['stress'][0]]:
                                    column['stress'][active['stress'][0]]+datainfo['stress']['len']]),'d').reshape(3,3)
    output += '\t'+'\t'.join(map(str,1.0/numpy.linalg.det(F)*numpy.dot(P,F.T).reshape(9)))  # [Cauchy] = (1/det(F)) * [P].[F_transpose]

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
