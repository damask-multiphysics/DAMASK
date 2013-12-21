#!/usr/bin/env python

import os,sys,string,re,platform
from optparse import OptionParser

validShell = {\
             'bash':    'DAMASK_env.sh' ,
             'csh':     'DAMASK_env.csh',
             'windows': 'DAMASK_env.bat',
              }

env_order = ['DAMASK_ROOT','DAMASK_BIN','PATH','PYTHONPATH','LD_LIBRARY_PATH']
environment = { 'DAMASK_ROOT':
                              [{'delete':'',
                                'substitute':'[DamaskRoot]',
                                'append': False,
                               },
                              ],
                'DAMASK_BIN':
                              [{'delete':'',
                                'substitute':'[os.path.join(DamaskRoot,"bin")]',
                                'append': False,
                               },
                              ],
                'PATH':
                              [{'delete':'"${DAMASK_BIN}"',
                                'substitute':'["${DAMASK_BIN}"]',
                                'append': True,
                               },
                              ],
                'PYTHONPATH':
                              [{'delete':'"${DAMASK_ROOT}/lib"',
                                'substitute':'["${DAMASK_ROOT}/lib"]',
                                'append': True,
                               },
                              ], 
                'LD_LIBRARY_PATH':
                              [{'activate':  'pathInfo["acml"]',
                                'delete':    '"acml"',                                             # what keywords trigger item deletion from existing path 
                                'substitute':'[os.path.join(pathInfo["acml"],"ifort64_mp/lib"),\
                                               os.path.join(pathInfo["acml"],"ifort64/lib"),\
                                               os.path.join(pathInfo["acml"],"gfortran64_mp/lib"),\
                                               os.path.join(pathInfo["acml"],"gfortran64/lib")]',  # what to substitute for deleted path items
                                'append': True,                                                    # whether new entries append to existing ${env}
                               },
                               {'activate':  'pathInfo["lapack"]',
                                'delete':    '[os.path.join(pathInfo["lapack"],"lib"),\
                                               os.path.join(pathInfo["lapack"],"lib64")]',         # deleted current (same) entry
                                'substitute':'[os.path.join(pathInfo["lapack"],"lib"),\
                                               os.path.join(pathInfo["lapack"],"lib64")]',         # what to substitute for deleted path
                                'append': True,                                                    # whether new entries append to existing ${env}
                               },
                               {'activate':  'pathInfo["fftw"]',
                                'delete':    '[os.path.join(pathInfo["fftw"],"lib")]',             # deleted current (same) entry
                                'substitute':'[os.path.join(pathInfo["fftw"],"lib")]',             # what to substitute for deleted path items
                                'append': True,                                                    # whether new entries append to existing ${env}
                               },
                              ],
               }

parser = OptionParser(usage="%prog [options]", description = """
Sets up your shell resource to be compatible with DAMASK. 
""" + string.replace('$Id: setup_shellrc.py 2126 2013-02-01 13:22:40Z MPIE\m.diehl $','\n','\\n')
)

parser.add_option("-s","--shell", type="string", 
                  dest = "shell", \
                  help = "type of shell, e.g. "+', '.join(validShell.keys())+" [%default]")

parser.set_defaults(shell = 'windows' if platform.system() == 'Windos' else 'bash')
(options, args) = parser.parse_args()

pathInfo = {}

DamaskRoot = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),'../'))
try:                                                                                              # check for user-defined pathinfo
  file = open(os.path.join(DamaskRoot,'lib/pathinfo'))
  content = map(lambda string: string.strip(),file.readlines())
  file.close()
  for line in content:
    if not (line.startswith('#') or line == ''):
      items = line.split() + ['','']
      pathInfo[items[0].lower()] = {False: os.path.normpath(os.path.join(DamaskRoot,'lib/',items[1])),
                                    True: items[1]}[items[1] == '']
except:
  pass

theShell = options.shell.lower()
if theShell not in validShell:
  print('invalid shell %s, assuming bash', theShell)
  theShell='bash'

theSourceFile = validShell[theShell]
thePath = os.path.join(DamaskRoot,theSourceFile)
if os.path.exists(thePath):
  print('parsing existing %s'%thePath)

  currentSourceFile = open(thePath)
  content = map(string.strip,currentSourceFile.readlines())
  currentSourceFile.close()

  match = {}
  for var in env_order: match[var] = False
  
  output = []
      
  for line in content:
    for var in env_order:
      m = re.search(\
                    {\
                     'bash':    r'^(.*? %s=)([^;]*)(.*)$'%var,
                     'csh':     r'^(\s*setenv\s+%s\s+)([^;]*)(.*)$'%var,
                     'windows': r'^(\s*SET\s+%s\s+)([^;]*)(.*)$'%var                              #case insesitiv search + second % sign
                    }[theShell],line)
      if m:
        n = re.search(r'(%s)'%var, line)
        o = re.search(r'(#)', line)
        if o:
          if o.start() < n.start():
            print 'skipped potential comment line, please check!'
            continue
        match[var] = True
        items = m.group(2).split(':')
        for piece in environment[var]:
          try:
            if 'activate' not in piece or eval(piece['activate']) != '':
              if piece['append']:
                if 'delete' in piece:
                  killer = eval(piece['delete'])
                  if type(killer) == str:
                    items = [path for path in items if killer not in path]
                  if type(killer) == list:
                    items = [path for path in items if path not in killer]
                items += eval(piece['substitute'])
              else:
                items = eval(piece['substitute'])
          except:
            pass
        line = m.group(1)+':'.join(items)+m.group(3)

    output.append(line)

for var in env_order:
  if not match[var]:
    items = ['${%s}'%var]
    for piece in environment[var]:
      try:
        if 'activate' not in piece or eval(piece['activate']) != '':
          if piece['append']:
            items += eval(piece['substitute'])
          else:
            items =  eval(piece['substitute'])
      except:
        pass
    output.append({\
       'bash':    'export %s=%s'%(var,':'.join(items)),
       'csh':     'setenv %s %s'%(var,':'.join(items)),
       'windows': 'SET %s %s'%(var,':'.join(items)),
      }[theShell])

print('writing to %s'%thePath)
newSourceFile = open(thePath,'w')
newSourceFile.write('\n'.join(output)+'\n')
newSourceFile.close()


