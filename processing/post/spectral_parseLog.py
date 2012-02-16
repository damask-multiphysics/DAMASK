#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,string,re,damask
from optparse import OptionParser, Option

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



# ----------------------- MAIN -------------------------------


parsing = [ \
           ['__CLOCK__',                   r'^Loadcase \d+ Increment \d+/\d+ @ Iteration \d+/\d+'],
           ['increment',                   r'^Loadcase \d+ Increment (\d+)/\d+ @ Iteration \d+/\d+'],
           ['iteration',                   r'^Loadcase \d+ Increment \d+/\d+ @ Iteration (\d+)/\d+'],
           ['error_divergence',            r'^error divergence = [ +-eE0123456789.]+ \( ([+-eE0123456789.]+)'],
           ['error_stress',                r'^error stress = [ +-eE0123456789.]+ \( ([+-eE0123456789.]+)'],
           ['crit_divergence',             r'^error divergence = ([ +-eE0123456789.]+) \( [+-eE0123456789.]+'],
           ['crit_stress',                 r'^error stress = ([ +-eE0123456789.]+) \( [+-eE0123456789.]+'],
           ['max(sym(deltaF))',            r'^max symmetrix correction of deformation: ([ +-eE0123456789.]+)'],
           ['max(skew(deltaF))',           r'^max skew      correction of deformation: ([ +-eE0123456789.]+)'],
           ['max(sym/skew(avg(deltaF)))',  r'^max sym/skew of avg correction:         ([ +-eE0123456789.]+)'],
           ['det(Fbar)',                   r'^determinant of new deformation: ([ +-eE0123456789.]+)'],
           ['max(det(F))',                 r'^max determinant of deformation:([ +-eE0123456789.]+)'],
           ['min(det(F))',                 r'^min determinant of deformation:([ +-eE0123456789.]+)'],
           ['div_FT_max',                  r'^error divergence  FT  max = ([ +-eE0123456789.]+)'],
           ['div_Real_RMS',                r'^error divergence Real RMS = ([ +-eE0123456789.]+)'],
           ['div_Real_max',                r'^error divergence Real max = ([ +-eE0123456789.]+)'],
           ['real(error_FT)',              r'^max FT relative error ([ +-eE0123456789.]+) [ +-eE0123456789.]+'],
           ['img(error_FT)',               r'^max FT relative error [ +-eE0123456789.]+ ([ +-eE0123456789.]+)'],
           ['real(error_iFT)',             r'^max iFT relative error ([ +-eE0123456789.]+)'],
           ['error_postProc',              r'^max deviat. from postProc =([ +-eE0123456789.]+)'],
           ['loadcase',                    r'^Loadcase (\d+) Increment \d+/\d+ @ Iteration \d+/\d+']
          ]

parser = OptionParser(option_class=extendedOption, usage='%prog [options] spectralOut[s]', description = """
Generate ASCIItable with data per iteration.

""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-p','--prefix', dest='prefix', type='string',\
                                   help='prefix of output ASCIItable [%default]')
parser.set_defaults(prefix = 'logParsing_')

(options, filenames) = parser.parse_args()


# ------------------------------------------ pre-compile regular expressions ---------------------------------------  

for i,[what,regexp] in enumerate(parsing):
  parsing[i][1] = re.compile(regexp)

# ------------------------------------------ setup file handles ---------------------------------------  

files = []
if filenames == []:
  files.append({'name':'STDIN', 'input':sys.stdin, 'output':sys.stdout})
else:
  for name in filenames:
    if os.path.exists(name):
      (path,file) = os.path.split(name)
      files.append({'name':name, 'input':open(name),
                                'output':open(os.path.join(path,options.prefix+os.path.splitext(file)[0]+'.txt'), 'w')})

# ------------------------------------------ loop over input files ---------------------------------------  

for file in files:
  if file['name'] != 'STDIN': print file['name']
  
  table = damask.ASCIItable(fileOut=file['output'],buffered=False)            # make unbuffered ASCII_table
  table.info_append(string.replace('$Id$','\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))
  table.labels_append([x[0] for x in parsing[1:]])                            # exclude __CLOCK__
  
# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ parse log ---------------------------------------  

  for line in file['input']:
    for what,regexp in parsing:
      m = regexp.search(line)
      if m:
        if what != '__CLOCK__': table.data_set(m.group(1),what)             # unless clock ticking detected...
        else:
          table.data_write()
          table.data_clear()

  if table.data != []: table.data_write()
  
# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close                                                    # close output ASCII table