#!/usr/bin/python
# -*- coding: UTF-8 no BOM -*-

import threading,os,string
import numpy as np
from optparse import OptionParser
from shutil import copy2
from re import split
import damask

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID   = ' '.join([scriptName,damask.version])

def list_split(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

#---------------------------------------------------------------------------------------------------
class myThread (threading.Thread):
  """Runner"""

  def __init__(self, threadID):
    threading.Thread.__init__(self)
    self.threadID = threadID
  def run(self):
    s.acquire()
    conv=isFinished()
    s.release()
    while conv:
      doSim(4.,self.name)
      s.acquire()
      conv=isFinished()
      s.release()

def doSim(delay,thread):
  global dirCurrent
  s.acquire()
  delta_angle = offsetPhi()
  if len(str(delta_angle)) > 5:
    file_angle = str(delta_angle)[:5]
  else:
    file_angle = str(delta_angle)
  dire = dirCurrent+'/'+file_angle

  if not os.path.isdir(dire):
    os.mkdir(dire,0755)
  for file in [options.geometry+'.geom',options.load+'.load','numerics.config']:
    copy2(dirCurrent+'/'+file, dire)
  newMaterialConfig(dirCurrent,delta_angle)

  os.chdir(dire)
  if not os.path.isfile('%s_%s.spectralOut'%(options.geometry,options.load)):
    print('starting uniaxial tension in direction of angle %s from %s'%(file_angle,thread))
    s.release()
    damask.util.execute('DAMASK_spectral -g %s -l %s'%(options.geometry,options.load))
  else: s.release()

  s.acquire()
  if not os.path.isfile('./%s/%s_%s.txt'%('Rvalues',options.geometry,options.load)):
    print('starting post processing for angle %s from %s'%(file_angle,thread))
    s.release()
    damask.util.execute('postResults --cr f,p -d %s %s_%s.spectralOut'%('Rvalues',options.geometry,options.load))
    damask.util.execute('addCauchy ./%s/%s_%s.txt'%('Rvalues',options.geometry,options.load))
    damask.util.execute('addStrainTensors -l -v ./%s/%s_%s.txt'%('Rvalues',options.geometry,options.load))
    print('post processing for angle %s from %s is finished'%(file_angle,thread))

  else:
    s.release()
  os.chdir(dirCurrent)

def isFinished():
  global N_simulations

  if N_simulations < options.number:
    return True
  else:
    return False

def offsetPhi():
  global N_simulations #, N_tensile
  N_simulations+=1
  return np.linspace(0,90,options.number)[N_simulations-1]

def newMaterialConfig(dire,angle):
  filename   = '/material.config'
  if len(str(angle)) > 5:
    file_angle = str(angle)[:5]
  else:
    file_angle = str(angle)
  f = open(dire+'/'+file_angle+filename,'w')
  data = open(dire+filename, 'r').readlines()

  for line in data:
    if '(gauss)' in line:
      linesplit = split(r'[;,\s,\t]\s*',line)
      index = linesplit.index('phi1')
      phi = float(linesplit[index+1]) - angle
      if phi > 360. : 
        phi = phi-360.0
      elif phi < 0.0: 
        phi = phi+360.0
      index2 = line.index(linesplit[index+1])
      line2 = line[:index2]+str(phi)+line[index2+len(str(float(linesplit[index+1]))):]
    else:
      line2 = line
    f.write(line2)
  f.close()

# --------------------------------------------------------------------
#                                MAIN
# --------------------------------------------------------------------

parser = OptionParser(option_class=damask.extendableOption, usage='%prog options [file[s]]', description = """
Calculate the deformation anisotropic coefficients (r values) and 
strength anisotropic coefficients (normalized yield stress)

""", version=string.replace(scriptID,'\n','\\n')
)

parser.add_option('-l','--load' ,
                  dest='load', type='string', 
                  help='name of the load file [%default]', metavar='string')
parser.add_option('-g','--geometry',
                  dest='geometry', type='string',
                  help='name of the geometry file [%default]', metavar='string')
parser.add_option('-s', '--strain',
                  dest='strain', type='string', action='callback', callback=list_split,
                  help='threshold strains, using comma to seperate multiple strains [%default]', metavar='string')
parser.add_option('-t','--threads',
                  dest='threads', type='int',
                  help='number of parallel executions [%default]',  metavar='int')
parser.add_option('-n','--number',
                  dest='number', type='int',
                  help='Number of uni-axial tensile tests [%default]',  metavar='int')

parser.set_defaults(geometry = '20grains16x16x16')
parser.set_defaults(load     = 'tensionX')
parser.set_defaults(threads  = 1)
parser.set_defaults(number   = 7)
parser.set_defaults(strain   = ('0.0025','0.025'))

options = parser.parse_args()[0]

if not os.path.isfile(options.geometry+'.geom'):
  parser.error('geometry file %s.geom not found'%options.geometry)
if not os.path.isfile('material.config'):
  parser.error('material.config file not found')

thresStrain = [float(i) for i in options.strain]

if not os.path.isfile(options.load+'.load'):
  maxstrain = max(thresStrain)
  f11 = str(1.2*maxstrain + 1.0)
  if maxstrain < 0.05:
    maxincs = '60'; maxtime = '30'
  else:
    maxincs = str(int(maxstrain*3000)); maxtime = str(maxstrain*600)

  f=open('tensionX.load','w')
  uniaxial_load = 'f '      +f11+ ' 0 0 0 * 0 0 0 *  ' + \
                  'stress ' +'* * *  * 0 *  * * 0  '   + \
                  'time '+ maxtime + ' incs ' + maxincs + ' freq 3'
  f.write(uniaxial_load)
  f.close()

N_simulations=0
s=threading.Semaphore(1)
threads=[]
dirCurrent = os.getcwd()
for i in range(options.threads):
  threads.append(myThread(i))
  threads[i].start()

for i in range(options.threads):
  threads[i].join()

anisotropy = {}
for i,delta_angle in enumerate(np.linspace(0,90,options.number)):
  if len(str(delta_angle)) > 5:
    file_angle = str(delta_angle)[:5]
  else:
    file_angle = str(delta_angle)
  refFile = dirCurrent+'/'+file_angle+'/Rvalues/%s_%s.txt'%(options.geometry,options.load)
  if not os.path.isfile(refFile):
    print('post processing file is not found for the angle %s'%file_angle)
  else:
    File = open(refFile)
    table = damask.ASCIItable(File)
    table.head_read()
    if not set(['1_ln(V)','5_ln(V)','9_ln(V)','1_Cauchy']).issubset(set(table.labels)):
      print 'data missing in direction %s'%str(delta_angle)
    table.data_readArray(['%i_Cauchy'%(i+1) for i in xrange(9)]+['%i_ln(V)'%(i+1) for i in xrange(9)])
    line, lines = 0, np.shape(table.data)[0]
    aniso_thres = {}
    for threshold in thresStrain:
      while line < lines:
        if abs(table.data[line,9])>= threshold: # 1_ln(V), e_xx
          upper,lower = abs(table.data[line,9]),abs(table.data[line-1,9]) # values for linear interpolation
          interplate = lambda x_d, x_u : x_d * (upper-threshold)/(upper-lower) + x_u * (threshold-lower)/(upper-lower)
          e_yy = interplate (table.data[line-1,13], table.data[line,13])
          e_zz = interplate (table.data[line-1,17], table.data[line,17])
          s_xx = interplate (table.data[line-1,0 ], table.data[line,0 ])
          aniso_thres[str(threshold)] = [e_yy/e_zz, s_xx]
          break
        else:
          line+=1
    anisotropy[file_angle] = aniso_thres

f = open('./anisotropy.txt','w')
f.write('4    header \n')
f.write('#the first row mean the threshold strain e_xx \n#the first column means loading directions \n')
f.write('#none means the data is unavailable\n# \n')
title = ['*R values (Lankford coefficients) \n', '# \n*Normalized yield stress \n']
for i in xrange(2):
  f.write(title[i])
  f.write(' '*10+len(thresStrain)*'%-12.4f'%tuple(thresStrain)+'\n')
  for j,delta_angle in enumerate(np.linspace(0,90,options.number)):
    if len(str(delta_angle)) > 5:
      file_angle = str(delta_angle)[:5]
    else:
      file_angle = str(delta_angle)

    if file_angle in anisotropy.keys():
      aniso_dic_ang = anisotropy[file_angle]
      aniso_list_ang_strain = []; writeformat = ''
      for threshold in thresStrain:
        if str(threshold) in aniso_dic_ang.keys():
          if i == 1: # extract the normalized stress
            aniso = aniso_dic_ang[str(threshold)][i]/anisotropy['0.0'][str(threshold)][i]
          else:      # extract r value
            aniso = aniso_dic_ang[str(threshold)][i]
          aniso_list_ang_strain.append(aniso)
          writeformat = writeformat+'%-12.6f'
        else:
          aniso_list_ang_strain.append('none')
          writeformat = writeformat+'%-12s'
      f.write('%-10s'%file_angle + writeformat%(tuple(aniso_list_ang_strain))+'\n')
f.close()
