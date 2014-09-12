#!/usr/bin/env python
# -*- coding: UTF-8 no BOM -*-

import os, sys, math, string, re, time
import damask
from optparse import OptionParser, OptionGroup, Option



# -----------------------------
class MyOption(Option):
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

      

# -----------------------------
def ParseOutputFormat(filename,homogID,crystID,phaseID):
#
# parse .output* files in order to get a list of outputs 
# -----------------------------
  
  myID = {  
          'Homogenization': homogID,
          'Crystallite':    crystID,
          'Constitutive':   phaseID,
         }
  format = {}
  
  for what in ['Homogenization','Crystallite','Constitutive']:
    content = []
    format[what] = {'outputs':{},'specials':{'brothers':[]}}
    for prefix in ['']+map(str,range(1,17)):
      if os.path.exists(prefix+filename+'.output'+what):
        try:
          file = open(prefix+filename+'.output'+what)
          content = file.readlines()
          file.close()
          break
        except:
          pass
    
    if content == []: continue                                              # nothing found...
    
    tag = ''
    tagID = 0
    for line in content:
      if re.match("\s*$",line) or re.match("#",line):                       # skip blank lines and comments
        continue
      m = re.match("\[(.+)\]",line)                                         # look for block indicator
      if m:                                                                 # next section
        tag = m.group(1)
        tagID += 1
        format[what]['specials']['brothers'].append(tag)
        if tag == myID[what] or (myID[what].isdigit() and tagID == int(myID[what])):
          format[what]['specials']['_id'] = tagID
          format[what]['outputs'] = []
          tag = myID[what]
      else:                                           # data from section
        if tag == myID[what]:
          (output,length) = line.split()
          output.lower()
          if length.isdigit():
            length = int(length)
          if re.match("\((.+)\)",output):                           # special data, e.g. (Ngrains)
            format[what]['specials'][output] = length
          elif length > 0:
            format[what]['outputs'].append([output,length])
    
    if not '_id' in format[what]['specials']:
      print "\nsection '%s' not found in <%s>"%(myID[what], what)
      print '\n'.join(map(lambda x:'  [%s]'%x, format[what]['specials']['brothers']))
  
  return format


# -----------------------------
def ParsePostfile(p,filename, outputFormat, legacyFormat):
#
# parse postfile in order to get position and labels of outputs
# needs "outputFormat" for mapping of output names to postfile output indices
# -----------------------------

  startVar = {True: 'GrainCount',
              False:'HomogenizationCount'}

  # --- build statistics

  stat = { \
  'IndexOfLabel': {}, \
  'Title': p.title(), \
  'Extrapolation': p.extrapolate, \
  'NumberOfIncrements': p.increments() - 1, \
  'NumberOfNodes': p.nodes(), \
  'NumberOfNodalScalars': p.node_scalars(), \
  'LabelOfNodalScalar': [None]*p.node_scalars() , \
  'NumberOfElements': p.elements(), \
  'NumberOfElementalScalars': p.element_scalars(), \
  'LabelOfElementalScalar': [None]*p.element_scalars() , \
  'NumberOfElementalTensors': p.element_tensors(), \
  'LabelOfElementalTensor': [None]*p.element_tensors(), \
  }

  # --- find labels 

  for labelIndex in range(stat['NumberOfNodalScalars']):
    label =  p.node_scalar_label(labelIndex)
    stat['IndexOfLabel'][label] = labelIndex
    stat['LabelOfNodalScalar'][labelIndex] = label

  for labelIndex in range(stat['NumberOfElementalScalars']):
    label =  p.element_scalar_label(labelIndex)
    stat['IndexOfLabel'][label] = labelIndex
    stat['LabelOfElementalScalar'][labelIndex] = label

  for labelIndex in range(stat['NumberOfElementalTensors']):
    label =  p.element_tensor_label(labelIndex)
    stat['IndexOfLabel'][label] = labelIndex
    stat['LabelOfElementalTensor'][labelIndex] = label
  
  if 'User Defined Variable 1' in stat['IndexOfLabel']:       # output format without dedicated names?
    stat['IndexOfLabel'][startVar[legacyFormat]] = stat['IndexOfLabel']['User Defined Variable 1']  # adjust first named entry
  
  if startVar[legacyFormat] in stat['IndexOfLabel']:          # does the result file contain relevant user defined output at all?
    startIndex = stat['IndexOfLabel'][startVar[legacyFormat]]
    stat['LabelOfElementalScalar'][startIndex] = startVar[legacyFormat]
    
    # We now have to find a mapping for each output label as defined in the .output* files to the output position in the post file
    # Since we know where the user defined outputs start ("startIndex"), we can simply assign increasing indices to the labels
    # given in the .output* file  

    offset = 1
    if legacyFormat:
      stat['LabelOfElementalScalar'][startIndex + offset] = startVar[not legacyFormat]    # add HomogenizationCount as second
      offset += 1
    
    for (name,N) in outputFormat['Homogenization']['outputs']:
      for i in range(N):
        label = {False:   '%s'%(    name),
                  True:'%i_%s'%(i+1,name)}[N > 1]
        stat['IndexOfLabel'][label] = startIndex + offset
        stat['LabelOfElementalScalar'][startIndex + offset] = label
        offset += 1
    
    if not legacyFormat:
      stat['IndexOfLabel'][startVar[not legacyFormat]] = startIndex + offset
      stat['LabelOfElementalScalar'][startIndex + offset] = startVar[not legacyFormat]        # add GrainCount
      offset += 1

    if '(ngrains)' in outputFormat['Homogenization']['specials']:
      for grain in range(outputFormat['Homogenization']['specials']['(ngrains)']):

        stat['IndexOfLabel']['%i_CrystalliteCount'%(grain+1)] = startIndex + offset              # report crystallite count
        stat['LabelOfElementalScalar'][startIndex + offset] = '%i_CrystalliteCount'%(grain+1)    # add GrainCount
        offset += 1

        for (name,N) in outputFormat['Crystallite']['outputs']:                           # add crystallite outputs
          for i in range(N):
            label = {False:   '%i_%s'%(grain+1,    name),
                      True:'%i_%i_%s'%(grain+1,i+1,name)}[N > 1]
            stat['IndexOfLabel'][label] = startIndex + offset
            stat['LabelOfElementalScalar'][startIndex + offset] = label
            offset += 1

        stat['IndexOfLabel']['%i_ConstitutiveCount'%(grain+1)] = startIndex + offset      # report constitutive count
        stat['LabelOfElementalScalar'][startIndex + offset] = '%i_ConstitutiveCount'%(grain+1)    # add GrainCount
        offset += 1

        for (name,N) in outputFormat['Constitutive']['outputs']:                               # add constitutive outputs
          for i in range(N):
            label = {False:   '%i_%s'%(grain+1,    name),
                      True:'%i_%i_%s'%(grain+1,i+1,name)}[N > 1]
            stat['IndexOfLabel'][label] = startIndex + offset
            try:
              stat['LabelOfElementalScalar'][startIndex + offset] = label
            except IndexError:
              print 'trying to assign %s at position %i+%i'%(label,startIndex,offset)
              sys.exit(1)
            offset += 1
  
  return stat


# -----------------------------
def GetIncrementLocations(p,Nincrements,options):
#
# get mapping between positions in postfile and increment number
# -----------------------------
  
  incAtPosition = {}
  positionOfInc = {}
  
  for position in range(Nincrements):
    p.moveto(position+1)
    incAtPosition[position] = p.increment            # remember "real" increment at this position
    positionOfInc[p.increment] = position            # remember position of "real" increment
  
  if not options.range:
    options.getIncrements = False
    locations = range(Nincrements)                   # process all positions
  else:
    options.range = list(options.range)              # convert to list
    if options.getIncrements:
      locations = [positionOfInc[x] for x in range(options.range[0],options.range[1]+1,options.range[2])
                                     if x in positionOfInc]
    else:
      locations = range( max(0,options.range[0]),
                         min(Nincrements,options.range[1]+1),
                         options.range[2] )
  
  increments = [incAtPosition[x] for x in locations] # build list of increments to process
  
  return [increments,locations]


# -----------------------------
def SummarizePostfile(stat,where=sys.stdout):
# -----------------------------

  where.write('\n\n')
  where.write('title:\t%s'%stat['Title'] + '\n\n')
  where.write('extraplation:\t%s'%stat['Extrapolation'] + '\n\n')
  where.write('increments:\t%i'%(stat['NumberOfIncrements']) + '\n\n')
  where.write('nodes:\t%i'%stat['NumberOfNodes'] + '\n\n')
  where.write('elements:\t%i'%stat['NumberOfElements'] + '\n\n')
  where.write('nodal scalars:\t%i'%stat['NumberOfNodalScalars'] + '\n\n  ' + '\n  '.join(stat['LabelOfNodalScalar']) + '\n\n')
  where.write('elemental scalars:\t%i'%stat['NumberOfElementalScalars'] + '\n\n  ' + '\n  '.join(stat['LabelOfElementalScalar']) + '\n\n')
  where.write('elemental tensors:\t%i'%stat['NumberOfElementalTensors'] + '\n\n  ' + '\n  '.join(stat['LabelOfElementalTensor']) + '\n\n')
  
  return True


# -----------------------------
def SummarizeOutputfile(format,where=sys.stdout):
# -----------------------------

  where.write('\nUser Defined Outputs')
  for what in format.keys():
    where.write('\n\n  %s:'%what)
    for output in format[what]['outputs']:
      where.write('\n    %s'%output)
  
  return True


# -----------------------------
def writeHeader(myfile,stat,geomtype):
# -----------------------------
    
  myfile.write('2\theader\n')
  myfile.write(string.replace('$Id$','\n','\\n')+
           '\t' + ' '.join(sys.argv[1:]) + '\n')
  if geomtype == 'nodebased':
    myfile.write('node')
    for i in range(stat['NumberOfNodalScalars']):
      myfile.write('\t%s'%''.join(stat['LabelOfNodalScalar'][i].split()))
    
  elif geomtype == 'ipbased':
    myfile.write('elem\tip')
    for i in range(stat['NumberOfElementalScalars']):
      myfile.write('\t%s'%''.join(stat['LabelOfElementalScalar'][i].split()))
  
  myfile.write('\n')
   
  return True



# -----------------------------
# MAIN FUNCTION STARTS HERE
# -----------------------------

# --- input parsing

parser = OptionParser(option_class=MyOption, usage='%prog [options] resultfile', description = """
Extract data from a .t16 (MSC.Marc) results file. 
""" + string.replace('$Id$','\n','\\n')
)

parser.add_option('-i','--info', action='store_true', dest='info', \
                  help='list contents of resultfile [%default]')
parser.add_option('-l','--legacy', action='store_true', dest='legacy', \
                  help='legacy user result block (starts with GrainCount) [%default]')
parser.add_option('-d','--dir', dest='dir', \
                  help='name of subdirectory to hold output [%default]')
parser.add_option('-r','--range', dest='range', type='int', nargs=3, \
                  help='range of positions (or increments) to output (start, end, step) [all]')
parser.add_option('--increments', action='store_true', dest='getIncrements', \
                  help='switch to increment range [%default]')
parser.add_option('-t','--type', dest='type', type='choice', choices=['ipbased','nodebased'], \
                  help='processed geometry type [ipbased and nodebased]')

group_material = OptionGroup(parser,'Material identifier')

group_material.add_option('--homogenization', dest='homog', \
                          help='homogenization identifier (as string or integer [%default])', metavar='<ID>')
group_material.add_option('--crystallite', dest='cryst', \
                          help='crystallite identifier (as string or integer [%default])', metavar='<ID>')
group_material.add_option('--phase', dest='phase', \
                          help='phase identifier (as string or integer [%default])', metavar='<ID>')

parser.add_option_group(group_material)

parser.set_defaults(info = False)
parser.set_defaults(legacy = False)
parser.set_defaults(dir = 'vtk')
parser.set_defaults(getIncrements= False)
parser.set_defaults(homog = '1')
parser.set_defaults(cryst = '1')
parser.set_defaults(phase = '1')

(options, files) = parser.parse_args()


# --- sanity checks

if files == []:
  parser.print_help()
  parser.error('no file specified...')

filename = os.path.splitext(files[0])[0]
if not os.path.exists(filename+'.t16'):
  parser.print_help()
  parser.error('invalid file "%s" specified...'%filename+'.t16')

sys.path.append(damask.solver.Marc().libraryPath('../../'))
try:
  from py_post import *
except:
  print('error: no valid Mentat release found')
  sys.exit(-1)

if not options.type :
  options.type = ['nodebased', 'ipbased']
else: 
  options.type = [options.type]


# --- initialize mesh data 

if damask.core.mesh.mesh_init_postprocessing(filename+'.mesh'):
  print('error: init not successful')
  sys.exit(-1)


# --- check if ip data available for all elements; if not, then .t19 file is required

p = post_open(filename+'.t16')
asciiFile = False
p.moveto(1)
for e in range(p.elements()):
  if not damask.core.mesh.mesh_get_nodeAtIP(str(p.element(e).type),1):
    if os.path.exists(filename+'.t19'):
      p.close()
      p = post_open(filename+'.t19')
      asciiFile = True
      break


# --- parse *.output and *.t16 file

outputFormat = ParseOutputFormat(filename,options.homog,options.cryst,options.phase)
p.moveto(1)
p.extrapolation('translate')
stat = ParsePostfile(p,filename,outputFormat,options.legacy)


# --- output info

if options.info:
  print '\n\nMentat release %s'%damask.solver.Marc().version('../../')
  SummarizePostfile(stat)
  SummarizeOutputfile(outputFormat)
  sys.exit(0)


# --- create output dir

dirname = os.path.abspath(os.path.join(os.path.dirname(filename),options.dir))
if not os.path.isdir(dirname):
  os.mkdir(dirname,0755)


# --- get positions

[increments,locations] = GetIncrementLocations(p,stat['NumberOfIncrements'],options)


# --- loop over positions

time_start = time.time()
for incCount,position in enumerate(locations):     # walk through locations
  p.moveto(position+1)                             # wind to correct position
  time_delta = (float(len(locations)) / float(incCount+1) - 1.0) * (time.time() - time_start)
  sys.stdout.write("\r(%02i:%02i:%02i) processing increment %i of %i..."%(time_delta//3600,time_delta%3600//60,time_delta%60,incCount+1,len(locations)))
  sys.stdout.flush()
  
  # --- write header 
  
  outFilename = {}
  for geomtype in options.type:
    outFilename[geomtype] = eval('"'+eval("'%%s_%%s_inc%%0%ii.txt'%(math.log10(max(increments+[1]))+1)")+'"%(dirname + os.sep + os.path.split(filename)[1],geomtype,increments[incCount])')
    with open(outFilename[geomtype],'w') as myfile:
      writeHeader(myfile,stat,geomtype)
      
      # --- write node based data
      
      if geomtype == 'nodebased':
        for n in range(stat['NumberOfNodes']):
          myfile.write(str(n))
          for l in range(stat['NumberOfNodalScalars']):
            myfile.write('\t'+str(p.node_scalar(n,l)))
          myfile.write('\n')
    
      # --- write ip based data
    
      elif geomtype == 'ipbased':
        for e in range(stat['NumberOfElements']):
          if asciiFile:
            print 'ascii postfile not yet supported'
            sys.exit(-1)
          else:
            ipData = [[]]
            for l in range(stat['NumberOfElementalScalars']):
              data = p.element_scalar(e,l)
              for i in range(len(data)):                                               # at least as many nodes as ips
                node = damask.core.mesh.mesh_get_nodeAtIP(str(p.element(e).type),i+1)  # fortran indexing starts at 1
                if not node: break                                                     # no more ips
                while i >= len(ipData): ipData.append([])
                ipData[i].extend([data[node-1].value])                                 # python indexing starts at 0
            for i in range(len(ipData)):
              myfile.write('\t'.join(map(str,[e,i]+ipData[i]))+'\n')
            
p.close()
sys.stdout.write("\n")

# ---------------------------       DONE     --------------------------------
