#!/usr/bin/env python3

import os
import sys
from io import StringIO
from optparse import OptionParser

import damask


scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def parameters(stretch,strain):
  """Albrecht Bertram: Elasticity and Plasticity of Large Deformations An Introduction (3rd Edition, 2012), p. 102."""
  return {
    'V#ln':    ('V',0.0),
    'U#ln':    ('U',0.0),
    'V#Biot':  ('V',-.5),
    'U#Biot':  ('U',+.5),
    'V#Green': ('V',-1.),
    'U#Green': ('U',+1.),
         }[stretch+'#'+strain]


# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [ASCIItable(s)]', description = """
Add column(s) containing given strains based on given stretches of requested deformation gradient column(s).

""", version = scriptID)

parser.add_option('-u','--right',
                  dest = 'right',
                  action = 'store_true',
                  help = 'material strains based on right Cauchy--Green deformation, i.e., C and U')
parser.add_option('-v','--left',
                  dest = 'left',
                  action = 'store_true',
                  help = 'spatial strains based on left Cauchy--Green deformation, i.e., B and V')
parser.add_option('-0','--logarithmic',
                  dest = 'logarithmic',
                  action = 'store_true',
                  help = 'calculate logarithmic strain tensor')
parser.add_option('-1','--biot',
                  dest = 'biot',
                  action = 'store_true',
                  help = 'calculate biot strain tensor')
parser.add_option('-2','--green',
                  dest = 'green',
                  action = 'store_true',
                  help = 'calculate green strain tensor')
parser.add_option('-f','--defgrad',
                  dest = 'defgrad',
                  action = 'extend',
                  metavar = '<string LIST>',
                  help = 'heading(s) of columns containing deformation tensor values [%default]')

parser.set_defaults(
                    defgrad     = ['f'],
                   )

(options,filenames) = parser.parse_args()
if filenames == []: filenames = [None]

if len(options.defgrad) > 1:
    options.defgrad = options.defgrad[1:]

stretches = []
strains = []

if options.right: stretches.append('U')
if options.left:  stretches.append('V')
if options.logarithmic: strains.append('ln')
if options.biot:        strains.append('Biot')
if options.green:       strains.append('Green')

if options.defgrad is None:
    parser.error('no data column specified.')

for name in filenames:
    damask.util.report(scriptName,name)
    
    table = damask.Table.from_ASCII(StringIO(''.join(sys.stdin.read())) if name is None else name)

    for defgrad in options.defgrad:
        F = table.get(defgrad).reshape((-1,3,3))
        for theStretch in stretches:
            for theStrain in strains:
                (t,m) = parameters(theStretch,theStrain)
                label = '{}({}){}'.format(theStrain,theStretch,defgrad if defgrad != 'f' else '')
                table.add(label,
                          damask.mechanics.strain_tensor(F,t,m).reshape((-1,9)),
                          scriptID+' '+' '.join(sys.argv[1:]))

    table.to_ASCII(sys.stdout if name is None else name)
