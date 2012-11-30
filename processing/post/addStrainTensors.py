#!/usr/bin/env python

import os,re,sys,math,numpy,string,damask
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



def operator(stretch,strain,eigenvalues):
  return { \
    'V#ln':    numpy.log(eigenvalues)                                 ,
    'U#ln':    numpy.log(eigenvalues)                                 ,
    'V#Biot':  ( numpy.ones(3,'d') - 1.0/eigenvalues )                ,
    'U#Biot':  ( eigenvalues - numpy.ones(3,'d') )                    ,
    'V#Green': ( numpy.ones(3,'d') - 1.0/eigenvalues*eigenvalues) *0.5,
    'U#Green': ( eigenvalues*eigenvalues - numpy.ones(3,'d'))     *0.5,
         }[stretch+'#'+strain]



# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=extendableOption, usage='%prog options [file[s]]', description = """
Add column(s) containing given strains based on given stretches of requested deformation gradient column(s).

""" + string.replace('$Id$','\n','\\n')
)


parser.add_option('-u','--right',       action='store_true', dest='right', \
                                        help='material strains based on right Cauchy--Green deformation, i.e., C and U')
parser.add_option('-v','--left',        action='store_true', dest='left', \
                                        help='spatial strains based on left Cauchy--Green deformation, i.e., B and V')
parser.add_option('-l','-0','--logarithmic', action='store_true', dest='logarithmic', \
                                        help='calculate logarithmic strain tensor')
parser.add_option('-b','-1','--biot',   action='store_true', dest='biot', \
                                        help='calculate biot strain tensor')
parser.add_option('-g','-2','--green',  action='store_true', dest='green', \
                                        help='calculate green strain tensor')
parser.add_option('-f','--deformation', dest='defgrad', action='extend', type='string', \
                                        help='heading(s) of columns containing deformation tensor values %default')

parser.set_defaults(right       = False)
parser.set_defaults(left        = False)
parser.set_defaults(logarithmic = False)
parser.set_defaults(biot        = False)
parser.set_defaults(green       = False)
parser.set_defaults(defgrad     = ['f'])

(options,filenames) = parser.parse_args()

stretches = []
stretch = {}
strains = []

if options.right: stretches.append('U')
if options.left:  stretches.append('V')
if options.logarithmic: strains.append('ln')
if options.biot: strains.append('Biot')
if options.green: strains.append('Green')

datainfo = {                                                               # list of requested labels per datatype
             'defgrad':     {'len':9,
                             'label':[]},
           }


if options.defgrad != None:   datainfo['defgrad']['label'] += options.defgrad


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
  if file['name'] != 'STDIN': print file['name']

  table = damask.ASCIItable(file['input'],file['output'],False)             # make unbuffered ASCII_table
  table.head_read()                                                         # read ASCII header info
  table.info_append(string.replace('$Id$','\n','\\n') + \
                    '\t' + ' '.join(sys.argv[1:]))

  active = {}
  column = {}
  head = []

  for datatype,info in datainfo.items():
    for label in info['label']:
      key = {True :'1_%s',
             False:'%s'   }[info['len']>1]%label
      if key not in table.labels:
        sys.stderr.write('column %s not found...\n'%key)
      else:
        if datatype not in active: active[datatype] = []
        if datatype not in column: column[datatype] = {}
        active[datatype].append(label)
        column[datatype][label] = table.labels.index(key)
        for theStretch in stretches:
          for theStrain in strains:
            table.labels_append(['%i_%s(%s)'%(i+1,theStrain,theStretch) 
                                for i in xrange(datainfo['defgrad']['len'])])         # extend ASCII header with new labels

# ------------------------------------------ assemble header ---------------------------------------  

  table.head_write()

# ------------------------------------------ process data ---------------------------------------  

  while table.data_read():                                                  # read next data line of ASCII table
  
    for datatype,labels in active.items():                                  # loop over vector,tensor
      for label in labels:                                                  # loop over all requested norms
        F = numpy.array(map(float,table.data[column['defgrad'][active['defgrad'][0]]:
                                             column['defgrad'][active['defgrad'][0]]+datainfo['defgrad']['len']]),'d').reshape(3,3)
        (U,S,Vh) = numpy.linalg.svd(F)
        R = numpy.dot(U,Vh)
        stretch['U'] = numpy.dot(numpy.linalg.inv(R),F)
        stretch['V'] = numpy.dot(F,numpy.linalg.inv(R))
        for theStretch in stretches:
          for i in range(9):
            if abs(stretch[theStretch][i%3,i//3]) < 1e-15:    # kill nasty noisy data
              stretch[theStretch][i%3,i//3] = 0.0
          (D,V) = numpy.linalg.eig(stretch[theStretch])       # eigen decomposition (of symmetric matrix)
          for i,eigval in enumerate(D):
            if eigval < 0.0:                                  # flip negative eigenvalues
              D[i] = -D[i]
              V[:,i] = -V[:,i]
          if numpy.dot(V[:,i],V[:,(i+1)%3]) != 0.0:           # check each vector for orthogonality
              V[:,(i+1)%3] = numpy.cross(V[:,(i+2)%3],V[:,i]) # correct next vector
              V[:,(i+1)%3] /= numpy.sqrt(numpy.dot(V[:,(i+1)%3],V[:,(i+1)%3].conj()))  # and renormalize (hyperphobic?)
          for theStrain in strains:
            d = operator(theStretch,theStrain,D)              # operate on eigenvalues of U or V
            eps = (numpy.dot(V,numpy.dot(numpy.diag(d),V.T)).real).reshape(9)  # build tensor back from eigenvalue/vector basis

            table.data_append(list(eps))

    table.data_write()                                                      # output processed line

# ------------------------------------------ output result ---------------------------------------  

  table.output_flush()                                                      # just in case of buffered ASCII table

  file['input'].close()                                                     # close input ASCII table
  if file['name'] != 'STDIN':
    file['output'].close                                                    # close output ASCII table
    os.rename(file['name']+'_tmp',file['name'])                             # overwrite old one with tmp new
