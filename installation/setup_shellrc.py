#!/usr/bin/env python

import os,sys,string,re
from optparse import OptionParser

pathInfo = {\
            'acml': 'lib/acml4.4.0',
            'fftw': 'lib/fftw',
            'msc':  '/msc',
           }

validShells = {\
               'bash':['.bashrc','.bash_profile'],
               'csh': ['.cshrc'],
              }

environment = {\
                'LD_LIBRARY_PATH': {'delete':'"acml"',                                                 # what keywords trigger item deletion from existing path 
                                    'substitute':'[os.path.join(pathInfo["acml"],"ifort64_mp/lib"),\
                                                   os.path.join(pathInfo["acml"],"ifort64/lib")]',     # what to substitute for deleted path items
                                    'append': True,                                                    # whether new entries append to existing ${env}
                                   },
                'PYTHONPATH':      {'delete':'os.path.join(DamaskRoot,"lib")',
                                    'substitute':'[os.path.join(DamaskRoot,"lib")]',
                                    'append': True,
                                   },
                'DAMASK_ROOT':     {'delete':'DamaskRoot',
                                    'substitute':'[DamaskRoot]',
                                    'append': False,
                                   },
                'PATH':            {'delete':'os.path.join(DamaskRoot,"bin")',
                                    'substitute':'[os.path.join(DamaskRoot,"bin")]',
                                    'append': True,
                                   },
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

DamaskRoot = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),'../'))

try:                                        # check for user-defined pathinfo
  file = open(os.path.join(DamaskRoot,'lib/pathinfo'))
  content = file.readlines()
  file.close()
  for line in content:
    pathinfo[line.split()[0].lower()] = os.path.normpath(line.split()[1])
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

      output = []
      for envVar in environment.keys(): environment[envVar]['matched'] = False
  
      for line in content:
        for envVar,data in environment.items():
          m = re.search(r'^(.*? %s=)([^;]*)(.*)$'%envVar,line)
          if m:
            substitute = eval(data['substitute'])+[path for path in m.group(2).split(':') if eval(data['delete']) not in path]
            line = m.group(1)+':'.join(substitute)+m.group(3)
            environment[envVar]['matched'] = True

        output.append(line)

      for envVar,data in environment.items():
        if not data['matched']:
          output.append('export %s=%s'%(envVar,':'.join(eval(data['substitute'])+{True:['${%s}'%envVar],False:[]}[data['append']])))

      rc = open(os.path.join(theHome,theRC),'w')
      rc.write('\n'.join(output)+'\n')
      rc.close()

elif theShell == 'csh':
  print 'csh not supported yet...'
