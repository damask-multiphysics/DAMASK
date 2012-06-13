#!/usr/bin/env python
import os,string,re,damask
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

    
parser = OptionParser(option_class=extendedOption, usage='%prog [options] datafile[s]', description = """
Writes version specific files for different MARC releases, adjustes the make file for the spectral solver and optionally compiles the spectral solver

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-c', '--compile', dest='spectralCompile', action='store_true', \
                  help='compiles sthe spectral solver [%default]')
parser.add_option('-m', '--make', dest='makeOptions', type='string', \
                  help='options parsed to make file when compiling spectral code [%default]')
parser.set_defaults(spectralCompile = False)
parser.set_defaults(makeOptions = 'F90=ifort')

(options, args) = parser.parse_args()

architectures = { 
                 'marc': { 
                          'parent': 'DAMASK_marc.f90', 
                          'versions' : ['%%MARCVERSION%%','2007r1','2008r1','2010','2011','2012'],
                         }, 
                }

bin_link = { \
        '' : [
                'DAMASK_spectral.exe',
             ],
            }
            
damaskEnv = damask.Environment('../../')          # script location relative to root
baseDir = damaskEnv.relPath('code/')

for arch in architectures:
	me = architectures[arch]
	try:
		parentFile = open(baseDir+os.sep+me['parent'])
		parentContent = parentFile.readlines()
		parentFile.close()
	except IOError:
		print 'unable to open',me['parent']
		continue

	
	for version in me['versions'][1:]:
		childFile = open(baseDir+os.sep+version.join(os.path.splitext(me['parent'])),'w')
		for line in parentContent:
			childFile.write(line.replace(me['versions'][0],version))
		childFile.close()

# changing dirs in makefile
makefile = open(os.path.join(baseDir,'Makefile'))
content = makefile.readlines()
makefile.close()
makefile = open(os.path.join(baseDir,'Makefile'),'w')
for line in content:
  m = re.match(r'(FFTW|IKML|ACML|LAPACK)ROOT\s*:?=',line)
  if m:
    if m.group(1).lower() in damaskEnv.pathInfo:
      substitution = damaskEnv.pathInfo[m.group(1).lower()]
    else:
      substitution = ''
    line = '%sROOT := %s\n'%(m.group(1),substitution)
  makefile.write(line)
makefile.close()

# compiling spectral code
if (options.spectralCompile):
  os.system('make --directory %s clean'%(baseDir))
  os.system('make --directory %s %s'%(baseDir,options.makeOptions))

for dir in bin_link:
  for file in bin_link[dir]:
    src = os.path.abspath(os.path.join(baseDir,dir,file))
    if (file == ''):
      sym_link = os.path.abspath(os.path.join(damaskEnv.binDir(),dir))
    else:
      sym_link = os.path.abspath(os.path.join(damaskEnv.binDir(),os.path.splitext(file)[0]))
    print sym_link,'-->',src
    if os.path.lexists(sym_link):
      os.remove(sym_link)    
    os.symlink(src,sym_link)
