#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,numpy
from optparse import OptionParser, OptionGroup, Option, SUPPRESS_HELP

# -----------------------------
class extendedOption(Option):
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



def outStdout(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    print cmd
  else:
    print cmd
  return

def outFile(cmd,locals):
  if cmd[0:3] == '(!)':
    exec(cmd[3:])
  elif cmd[0:3] == '(?)':
    cmd = eval(cmd[3:])
    locals['filepointer'].write(cmd+'\n')
  else:
    locals['filepointer'].write(cmd+'\n')
  return


def output(cmds,locals,dest):
  for cmd in cmds:
    if isinstance(cmd,list):
      output(cmd,locals,dest)
    else:
      {\
      'File': outFile,\
      'Stdout': outStdout,\
      }[dest](str(cmd),locals)
  return


# +++++++++++++++++++++++++++++++++++++++++++++++++++
def vtk_writeASCII_mesh(dim,res,origin,data):
# +++++++++++++++++++++++++++++++++++++++++++++++++++
  """ function writes data array defined on a rectilinear grid """
  N  = res[0]*res[1]*res[2]
  
  cmds = [\
          '# vtk DataFile Version 3.1',
          string.replace('powered by $Id$','\n','\\n'),
          'ASCII',
          'DATASET RECTILINEAR_GRID',
          'DIMENSIONS %i %i %i'%(res[0]+1,res[1]+1,res[2]+1),
          'X_COORDINATES %i float'%(res[0]+1),
          ' '.join(map(str,[i*dim[0]/res[0]+origin[0] for i in range(res[0]+1)])),
          'Y_COORDINATES %i float'%(res[1]+1),
          ' '.join(map(str,[i*dim[1]/res[1]+origin[1] for i in range(res[1]+1)])),
          'Z_COORDINATES %i float'%(res[2]+1),
          ' '.join(map(str,[i*dim[2]/res[2]+origin[2] for i in range(res[2]+1)])),
          'CELL_DATA %i'%N,
         ]
  
  for datatype in data:
    for item in data[datatype]:
      cmds += [\
               '%s %s float'%(datatype.upper()+{True:'',False:'S'}[datatype.lower().endswith('s')],item),
               'LOOKUP_TABLE default',
               [[['\t'.join(map(str,data[datatype][item][:,j,k]))] for j in range(res[1])] for k in range(res[2])]
              ]

  return cmds


# ----------------------- MAIN -------------------------------

identifiers = {
        'resolution': ['a','b','c'],
        'dimension':  ['x','y','z'],
        'origin':     ['x','y','z'],
          }
mappings = {
        'resolution': lambda x: int(x),
        'dimension':  lambda x: float(x),
        'origin':     lambda x: float(x),
          }

parser = OptionParser(option_class=extendedOption, usage='%prog [geomfile[s]]', description = """
Produce VTK rectilinear mesh of structure data from geom description

""" + string.replace('$Id$','\n','\\n')
)

(options, filenames) = parser.parse_args()

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN',
                'input':sys.stdin,
                'output':sys.stdout,
                'croak':sys.stderr,
               })
else:
  for name in filenames:
    if os.path.exists(name):
      files.append({'name':name,
                    'input':open(name),
                    'output':open(name+'_tmp','w'),
                    'croak':sys.stdout,
                    })

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': file['croak'].write(file['name']+'\n')

  #  get labels by either read the first row, or - if keyword header is present - the last line of the header

  firstline = file['input'].readline()
  m = re.search('(\d+)\s*head', firstline.lower())
  if m:
    headerlines = int(m.group(1))
    headers  = [file['input'].readline() for i in range(headerlines)]
  else:
    headerlines = 1
    headers = firstline

  content = file['input'].readlines()
  file['input'].close()

  info = {'resolution': [0,0,0],
          'dimension':  [0.0,0.0,0.0],
          'origin':     [0.0,0.0,0.0],
         }
  for header in headers:
    headitems = map(str.lower,header.split())
    if headitems[0] in identifiers.keys():
      for i in xrange(len(identifiers[headitems[0]])):
        info[headitems[0]][i] = mappings[headitems[0]](headitems[headitems.index(identifiers[headitems[0]][i])+1])

  if numpy.all(info['resolution'] == 0):
    file['croak'].write('no resolution info found.\n')
    continue

  if numpy.all(info['dimension'] == 0.0):
    file['croak'].write('no dimension info found.\n')
    continue

  file['croak'].write('resolution:     %s\n'%(' x '.join(map(str,info['resolution']))) + \
                      'dimension:      %s\n'%(' x '.join(map(str,info['dimension']))) + \
                      'origin:         %s\n'%(' : '.join(map(str,info['origin'])))
                     )
  
  data = {'scalar':{'structure':numpy.zeros(info['resolution'],'i')}}
  i = 0
  for line in content:  
    for item in map(int,line.split()):
      data['scalar']['structure'][i%info['resolution'][0],(i/info['resolution'][0])%info['resolution'][1],i/info['resolution'][0]/info['resolution'][1]] = item
      i += 1


  out = {}
  out['mesh'] = vtk_writeASCII_mesh(info['dimension'],info['resolution'],info['origin'],data)
  
  for what in out.keys():
    if file['name'] == 'STDIN':
      output(out[what],{},'Stdout')
    else:
      (head,tail) = os.path.split(file['name'])
      vtk = open(os.path.join(head,what+'_'+os.path.splitext(tail)[0]+'.vtk'), 'w')
      output(out[what],{'filepointer':vtk},'File')
      vtk.close()
