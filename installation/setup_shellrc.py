#!/usr/bin/env python

import os,sys,string,re
from optparse import OptionParser

validShells = {\
               'bash':['.bashrc','.bash_profile','.bash_login','.profile'],
               'csh': ['.cshrc'],
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
                                'delete':    '"acml"',                                            # what keywords trigger item deletion from existing path 
                                'substitute':'[os.path.join(pathInfo["acml"],"ifort64_mp/lib"),\
                                               os.path.join(pathInfo["acml"],"ifort64/lib")]',    # what to substitute for deleted path items
                                'append': True,                                                   # whether new entries append to existing ${env}
                               },
                               {'activate':  'pathInfo["lapack"]',
                                'delete':    '[os.path.join(pathInfo["lapack"],"lib"),\
                                               os.path.join(pathInfo["lapack"],"lib64")]',        # deleted current (same) entry
                                'substitute':'[os.path.join(pathInfo["lapack"],"lib"),\
                                               os.path.join(pathInfo["lapack"],"lib64")]',        # what to substitute for deleted path
                                'append': True,                                                   # whether new entries append to existing ${env}
                               },
                               {'activate':  'pathInfo["fftw"]',
                                'delete':    '[os.path.join(pathInfo["fftw"],"lib")]',            # deleted current (same) entry
                                'substitute':'[os.path.join(pathInfo["fftw"],"lib")]',            # what to substitute for deleted path items
                                'append': True,                                                   # whether new entries append to existing ${env}
                               },
                              ],
               }

parser = OptionParser(usage="%prog [options]", description = """
Sets up your shell resource to be compatible with DAMASK. 
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option("-s","--shell", type="string", 
                  dest = "shell", \
                  help = "type of shell, e.g. "+', '.join(validShells.keys())+" [%default]")

parser.set_defaults(shell = 'bash')

(options, args) = parser.parse_args()

pathInfo = {}

DamaskRoot = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),'../'))
try:                                        # check for user-defined pathinfo
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
theHome = os.getenv('USERPROFILE') or os.getenv('HOME')

if theShell == 'bash':
  for theRC in validShells[theShell]:
    thePath = os.path.join(theHome,theRC)
    if os.path.exists(thePath):
      rc = open(os.path.join(theHome,theRC))
      content = map(string.strip,rc.readlines())
      rc.close()

      match = {}
      for var in env_order: match[var] = False
  
      output = []
      
      for line in content:
        for var in env_order:
          m = re.search(r'^(.*? %s=)([^;]*)(.*)$'%var,line)
          if m:
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
          output.append('export %s=%s'%(var,':'.join(items)))

      rc = open(os.path.join(theHome,theRC),'w')
      rc.write('\n'.join(output)+'\n')
      rc.close()

elif theShell == 'csh':
  print 'csh not supported yet...'
