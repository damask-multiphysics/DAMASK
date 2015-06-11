#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os,sys,math,string
import numpy as np
from optparse import OptionParser
from PIL import Image
import damask

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]

#--------------------------------------------------------------------------------------------------
#                                MAIN
#--------------------------------------------------------------------------------------------------
parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """

Generate geometry description and from stacked 8bit tiff images. No orientation info is used.
Microstructure index is based on gray scale value (1...256).

""", version = scriptID)

parser.add_option('--homogenization',      dest='homogenization', type='int', metavar = 'int',
                  help='homogenization index for <microstructure> configuration [%default]')
              
parser.set_defaults(homogenization = 1)

(options,filenames) = parser.parse_args()


#--- setup file handles ---------------------------------------------------------------------------
files = []
for name in filenames:
  if os.path.exists(name):
    files.append({'name':name,
                  'input':open(name),
                  'output':open(name+'_tmp','w'),
                  'croak':sys.stdout,
                  })

#--- loop over input files ------------------------------------------------------------------------
for file in files:
  file['croak'].write('\033[1m' + scriptName + '\033[0m: ' + (file['name'] if file['name'] != 'STDIN' else '') + '\n')

  info = {
          'grid':   np.ones (3,'i'),
          'size':   np.zeros(3,'d'),
          'origin': np.zeros(3,'d'),
          'microstructures': 0,
          'homogenization':  options.homogenization,
         }

  ImageFile=Image.open(file['input'])
  i=0
  a=[]
  while True:
    try:
      ImageFile.seek(i)
      a.append(np.array(ImageFile.transpose(Image.FLIP_TOP_BOTTOM)))
      i+=1
    except EOFError:
      break
# http://docs.scipy.org/doc/scipy/reference/ndimage.html
# http://scipy-lectures.github.io/advanced/image_processing/
  imarray=np.asarray(a)
  info['grid'] = np.array(imarray.shape)[::-1]
  info['size'] = np.array(imarray.shape,'d')[::-1]
  info['microstructures']=len(np.unique(imarray))
#--- report ---------------------------------------------------------------------------------------
  file['croak'].write('grid     a b c:  %s\n'%(' x '.join(map(str,info['grid']))) +
                      'size     x y z:  %s\n'%(' x '.join(map(str,info['size']))) +
                      'origin   x y z:  %s\n'%(' : '.join(map(str,info['origin']))) +
                      'homogenization:  %i\n'%info['homogenization'] +
                      'microstructures: %i\n\n'%info['microstructures'])

  if np.any(info['grid'] < 1):
    file['croak'].write('invalid grid a b c.\n')
    continue
  if np.any(info['size'] <= 0.0):
    file['croak'].write('invalid size x y z.\n')
    continue


#--- write data -----------------------------------------------------------------------------------
  header = [' '.join([scriptID] + sys.argv[1:]),
            "grid\ta %i\tb %i\tc %i"%(info['grid'][0],info['grid'][1],info['grid'][2]),
            "size\tx %f\ty %f\tz %f"%(info['size'][0],info['size'][1],info['size'][2]),
            "origin\tx %f\ty %f\tz %f"%(info['origin'][0],info['origin'][1],info['origin'][2]),
            "microstructures\t%i"%info['microstructures'],
            "homogenization\t%i"%info['homogenization'],
            ]
  file['output'].write('\n'.join(['%i\theader'%(len(header))] + header) + '\n')
  np.savetxt(file['output'],imarray.reshape([info['grid'][1]*info['grid'][2],info['grid'][0]]),
                            fmt='%0'+str(1+int(math.log10(np.amax(imarray))))+'d'
  
#--- output finalization -------------------------------------------------------------------------- 
  if file['name'] != 'STDIN':
    file['output'].close()
    os.rename(file['name']+'_tmp', os.path.splitext(file['name'])[0] +'.geom')
