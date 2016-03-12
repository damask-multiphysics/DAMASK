# -*- coding: UTF-8 no BOM -*-

# This tool converts a msc.marc result file into the vtk format that 
# can be viewed by Paraview software (Kitware), or MayaVi (needs xml-vtk, or ...
# 
# About the vtk format: http://www.vtk.org/VTK/project/about.html
# Some example vtk files: http://people.sc.fsu.edu/~jburkardt/data/vtk/vtk.html
# www.paraview.org

import os,sys,re
import numpy as np

import py_post  # MSC closed source module to access marc result files

class MARC_POST():
  import re
  def __init__(self):
      self.projdir='./'

  def opent16(self,incr=None):
      self.fpath=os.path.join(self.projdir,self.postname)
      print('Trying to open ',self.fpath,' ...')
      self.p=py_post.post_open(self.fpath)
      if self.p is None:
        print('Could not open %s.'%self.postname); #return 'err'#; sys.exit(1)
        raise Exception('Could not open t16')
      print('Postfile %s%s is open ...'%(self.projdir,self.postname))
      self.maxincr=self.p.increments()       
      print('and has %i increments'%self.maxincr)
      if incr is None:
        self.p.moveto(self.maxincr-1)   
      else:
        self.p.moveto(incr+1)
      print('moved to increment ', self.p.increment)       
      self.p.extrapolation('translate') # linear, translate, average. query with p.extrapolate
      print('extrapolation method is ', self.p.extrapolate) 
      print('collecting model information')
      self.t16info(printFlag=0)
      print('t16 is open')
      self.p

  def t16info(self, printFlag=1
              ):
      if not self.p:
        self.p=self.opent16()#post_open(self.postname)
      print(self.p)
      oldincr=self.p.position
      if oldincr==0:
        self.p.moveto(1)
        
      self.nnodes=self.p.nodes() #self.p.node(self.nnodes) crashes; only to nnodes-1 possible
      self.nodes=range(0,self.nnodes)
      self.nscals=self.p.node_scalars(); #print 'nscals', nscals
      self.nscal_list=['']*self.nscals
      self.nel=self.p.elements()
      self.elscals=self.p.element_scalars()
      self.elscal_list=['']*self.elscals
      self.eltens=self.p.element_tensors()
      self.elten_list=['']*self.eltens
      for i in range(0,self.nscals):
          self.nscal_list[i]=self.p.node_scalar_label(i)
      for i in range (0,self.elscals):
          self.elscal_list[i]=self.p.element_scalar_label(i)
          if printFlag==1: print(i, self.elscal_list[i])
      for i in range (0,self.eltens):
          self.elten_list[i]=self.p.element_tensor_label(i)
          if printFlag==1: print(i, self.elten_list[i])
      for i in range(0,self.p.element_tensors()):
        if printFlag==1: print('Element Tensor: ', i, self.p.element_tensor_label(i))
      if printFlag==1: print('')
      for i in range(0,self.p.element_scalars()):
        if printFlag==1: print('Element Scalar: ', i, self.p.element_scalar_label(i))
      if oldincr==0:
        self.p.moveto(0)     

  def closet16(self):
      print('Trying to close FEM result file ...')
      try:
        if self.p:
          self.p.close()
          print('FEM result file closed.')
          self.p=None
        else:
          print('post object not open?')
      except:
        print('ERROR. Could not close FEM result file.')

  def getLabelNr(self, label=None, type='Scalar'): 
      if type[0]=='S' or type[0]=='s':   # element scalar
         labelNr=self.elscal_list.index(label)
      elif type[0]=='N':                 # node scalar
         labelNr=self.nscal_list.index(label)
      elif type[0]=='T' or type[0]=='t': # tensor
         labelNr=self.elten_list.index(label)
      print('Found label %s at index %i'%(label,labelNr))
      return labelNr
        
        
  def writeNodes2VTK(self, fobj):
      self.points=[]
      self.VTKcnt=200 # number of values per line in vtk file
      fobj.write('POINTS %i'%self.p.nodes()+' float\n')
      self.nodes_dict={} # store the node IDs in case of holes in the numbering
      for iNd in self.nodes:
        nd=self.p.node(iNd)
        disp=self.p.node_displacement(iNd)
        nd_xyz=[nd.x+disp[0], nd.y+disp[1], nd.z+disp[2]]
        self.points.append(nd_xyz) # for pyvtk
        fobj.write('%f %f %f \n'%
        (nd.x+disp[0], nd.y+disp[1], nd.z+disp[2]))
        self.nodes_dict[nd.id-1]=iNd
      fobj.write('\n')  
      print('Nodes written to VTK: %i'%self.p.nodes())
        
  def writeElements2VTK(self, fobj):
      fobj.write('\nCELLS %i %i'%(self.p.elements(),self.p.elements()*9)+'\n')
      self.cells=[] #for pyvtk
      for iEl in range(0,self.nel):
        el=self.p.element(iEl)
        cell_nodes=[] # for pyvtk
        ndlist=el.items
        for k in [0, 1, 2, 3, 4, 5, 6, 7]: # FOR CELL TYPE VTK_HEXAHEDRON
         node=ndlist[k]-1
         cell_nodes.append(self.nodes_dict[node])
        self.cells.append(cell_nodes) # for pyvtk
      for e in self.cells:
         fobj.write('8 ')
         for n in e:
           fobj.write('%6i '%n)
         fobj.write('\n')
      fobj.write('\nCELL_TYPES %i'%self.p.elements()+'\n')
      cnt=0
      for iEl in range(0,self.nel):
        cnt=cnt+1
        #fobj.write('11\n') #VTK_VOXEL
        fobj.write('12 ') #VTK_HEXAHEDRON
        if cnt>self.VTKcnt:
          fobj.write('\n');cnt=0
      fobj.write('\n')    
      print('Elements written to VTK: %i'%self.p.elements())
        
  def writeElScalars2NodesVTK(self,fobj):
      fobj.write('\nPOINT_DATA %i\n'%self.p.nodes())
      nScal=12
      nComponents=1+nScal
      fobj.write('SCALARS scalars float %i\n'%nComponents)
      fobj.write('LOOKUP_TABLE default\n')  
      idxScal=self.nscal_list.index('Displacement Z')    
      for iNd in self.nodes:
        fobj.write('%f '%(self.p.node_scalar(iNd,idxScal)))
        for iEl in range(0,self.nel):
          el=self.p.element(iEl)
          ndlist=el.items
          if (iNd+1) in ndlist:
            idx=ndlist.index(iNd+1)
            for iV in range(0,nScal):
              elData=self.p.element_scalar(iEl,35+iV)
              fobj.write('%f '%(elData[idx].value))
            break
        fobj.write('\n')  
      fobj.write('\n')  
        
  def writeNodeScalars2VTK(self,fobj):
      fobj.write('\nPOINT_DATA %i\n'%self.p.nodes())
      self.pointDataScalars=[]
      for idxNdScal in range(-3,self.nscals): #now include node x,y,z
        if idxNdScal>=0:
          datalabel=self.nscal_list[idxNdScal]
          datalabel=re.sub("\s",'_',datalabel)
        else:
          if idxNdScal==-3: datalabel='node.x'  
          if idxNdScal==-2: datalabel='node.y'
          if idxNdScal==-1: datalabel='node.z'
        fobj.write('SCALARS %s float %i\n'%(datalabel,1))#nComponents))
        fobj.write('LOOKUP_TABLE default\n')  
        cnt=0
        for iNd in range(0,self.nnodes):
          cnt=cnt+1
          if idxNdScal>=0:
            ndData=self.p.node_scalar(iNd,idxNdScal)
          else:
             nd=self.p.node(iNd)
             if idxNdScal==-3: ndData=nd.x
             if idxNdScal==-2: ndData=nd.y
             if idxNdScal==-1: ndData=nd.z
          fobj.write('%E '%(ndData))
          if cnt>self.VTKcnt:
            fobj.write('\n')
            cnt=0
        fobj.write('\n')              
      fobj.write('\n')        

  def writeElementData2VTK(self,fobj):
      self.sig_vMises=[]
      self.sig33=[]
      idx_sig_vMises=self.getLabelNr('Equivalent Von Mises Stress')
      idx_sig33=self.getLabelNr('Comp 33 of Cauchy Stress')
      fobj.write('\nCELL_DATA %i\n'%self.p.elements())
      for idxElScal in range(0,self.elscals):
        datalabel=self.elscal_list[idxElScal]
        datalabel=re.sub("\s",'_',datalabel)
        fobj.write('\n\nSCALARS %s float %i\n'%(datalabel,1))
        fobj.write('LOOKUP_TABLE default\n')  
        cnt=0
        for iEl in range(0,self.nel):
          cnt=cnt+1
          elData=self.p.element_scalar(iEl,idxElScal)
          avgScal=0.0
          if datalabel in ['phi1', 'PHI','phi2']: # Euler angles should not be averaged
            avgScal=avgScal+elData[0].value
          else:
            for IP in range(0,8):
              avgScal=avgScal+elData[IP].value
            avgScal=avgScal/8.
          fobj.write('%E '%(avgScal))
          if idxElScal==idx_sig_vMises:
            self.sig_vMises.append(avgScal)
          elif idxElScal==idx_sig33:  
            self.sig33.append(avgScal)
          if cnt>self.VTKcnt:
            fobj.write('\n')
            cnt=0
      fobj.write('\n')      

  def get_avg_el_scal(self,idxElScal):
      result=[]
      datalabel=self.elscal_list[idxElScal]
      print('Collecting %s from all elements'%datalabel)
      for iEl in range(0,self.nel):
        elData=self.p.element_scalar(iEl,idxElScal)
        avgScal=0.0
        for IP in range(0,8):
          avgScal=avgScal+elData[IP].value
        avgScal=avgScal/8. 
        result.append(avgScal)
      return result
      
  def writeUniaxiality2VTK(self,fobj):     
      datalabel='uniaxiality_sig_vMises_durch_sig33'
      fobj.write('SCALARS %s float %i\n'%(datalabel,1))
      fobj.write('LOOKUP_TABLE default\n')  
      cnt=0
      for iEl in range(0,self.nel):
        cnt=cnt+1                  
        if abs(self.sig_vMises[iEl])<1e-5:
          datum=0.
        else:  
          datum=self.sig33[iEl]/self.sig_vMises[iEl]
        fobj.write('%E '%(datum))
        if cnt>self.VTKcnt:
          fobj.write('\n')
          cnt=0
      fobj.write('\n')     

  def stress_per_element(self):
      self.stress=[]
      for iEl in range(0,self.nel):
        sig=self.avg_elten(2,elID=iEl)
        self.stress.append(sig[0])
        
  def mean_stress_per_element(self):
      self.mean_stress=[]
      for iEl in range(0,self.nel):
        sig=self.stress[iEl]
        self.mean_stress.append(self.meanStress(sig))

  def triaxiality_per_element(self):
    # classical triaxiality 
    # 1/3 : uniax tension
      self.triaxiality=[]
      for iEl in range(0,self.nel):
        t=self.mean_stress[iEl]/self.sig_vMises[iEl]
        self.triaxiality.append(t)
      
  def moreElData2VTK(self,fobj,data=[],label='datalabel'):
      fobj.write('SCALARS %s float %i\n'%(label,1))
      fobj.write('LOOKUP_TABLE default\n')  
      cnt=0
      for iEl in range(0,self.nel):
        cnt=cnt+1                  
        fobj.write('%E '%(data[iEl]))
        if cnt>self.VTKcnt:
          fobj.write('\n')
          cnt=0
      fobj.write('\n')             
      
  def calc_lode_parameter(self):
     self.lode=[]
     try:
       self.stress
     except:
       self.stress_per_element()
     for iEl in range(0,self.nel):
       sig=self.stress[iEl]
       lode=self.stress2lode(sig)
       self.lode.append(lode)
       
  def stress2lode(self,stress):
     [pStress,pAxes]=self.princStress(stress)
     s1=pStress[0]
     s2=pStress[1]
     s3=pStress[2]
     lode=(2*s2-s1-s3) / ( s1 - s3 )
     return lode
     
  def princStress(self, stress):
    """
    Function to compute 3D principal stresses and sort them.
    
    from: http://geodynamics.org/svn/cig/short/3D/PyLith/trunk/playpen/postproc/vtkcff.py
    """
    stressMat=np.array(stress)
    (princStress, princAxes) = np.linalg.eigh(stressMat)
    idx = princStress.argsort()
    princStressOrdered = princStress[idx]
    princAxesOrdered = princAxes[:,idx]
    return princStressOrdered, princAxesOrdered     

  def avg_elten(self,
     idxElTen, mat=0, elID=None):
     tensum=np.zeros((3,3));
     T=np.zeros((3,3));
     pts=0;
     avg=np.zeros((3,3));

     if elID is None:
       averaged_elements=range(0,self.nel)
     else:
       averaged_elements=[elID]
     for i in averaged_elements:
        if mat==0 or int(self.p.element_scalar(i,4)[0].value)==mat:
          T=self.p.element_tensor(i,idxElTen)
          for k in range (0,8):
             tensum[0][0] = tensum[0][0] + T[k].t11
             tensum[0][1] = tensum[0][1] + T[k].t12
             tensum[0][2] = tensum[0][2] + T[k].t13
             tensum[1][1] = tensum[1][1] + T[k].t22
             tensum[1][2] = tensum[1][2] + T[k].t23
             tensum[2][2] = tensum[2][2] + T[k].t33
             pts=pts+1
     avg=tensum/pts
     avg=self.fillComponents(avg)
     del [T]
     return (avg,tensum,pts)

  def fillComponents(self,
                     halftensor
                     ):
     halftensor[1][0]=halftensor[0][1]
     halftensor[2][0]=halftensor[0][2]
     halftensor[2][1]=halftensor[1][2]
     return halftensor
     
  def vMises(self,tensor33):
     t=tensor33
     s=(t[0,0]-t[1,1])**2+(t[1,1]-t[2,2])**2+(t[0,0]-t[2,2])**2+\
       6*(t[0,1]**2+t[1,2]**2+t[2,0]**2)
     vM=np.sqrt(s/2.)
     return vM
     
  def meanStress(self,tensor33):
     t=tensor33
     s=t[0,0]+t[1,1]+t[2,2]
     ms=s/3.
     return ms
     
  def invariants(self,tensor33):
     t=tensor33
     I1=t[0,0]+t[1,1]+t[2,2]
     I2=t[0,0]*t[1,1]+t[1,1]*t[2,2]+t[0,0]*t[2,2]-\
        t[0,1]**2-t[1,2]**2-t[0,2]**2 
     I3=t[0,0]*t[1,1]*t[2,2]+\
        2*t[0,1]*t[1,2]*t[2,0]-\
        t[2,2]*t[0,1]**2-t[0,0]*t[1,2]**2-t[1,1]*t[0,2]**2
     return [ I1, I2, I3 ]

      
class VTK_WRITER():
  """
  The resulting vtk-file can be imported in Paraview 3.12
 
  Then use Filters: Cell Data to Point Data + Contour 
  to plot semi-transparent iso-surfaces.
  """

  import re
  def __init__(self):
    self.p=MARC_POST() # self.p
    
  def openFile(self, filename='test.vtp'):
      self.f=open(filename,'w+')
      self.fname=filename
      
  def writeFirstLines(self,
      vtkFile=None,
      version='2.0', 
      comment='Test',
      dformat='ASCII', # BINARY | [ASCII]
      dtype='UNSTRUCTURED_GRID' # UNSTRUCTURED GRID
      ):
    if vtkFile is None:
      vtkFile=self.f    
    # First Line contains Data format version
    self.versionVTK=version
    vtkFile.write('# vtk DataFile Version %s\n'%self.versionVTK)
    # Comment goes to 2nd line and has maximum 256 chars
    vtkFile.write(comment+'\n')
    vtkFile.write(dformat+'\n')
    vtkFile.write('DATASET '+dtype+'\n')
    

  def marc2vtkBatch(self):
    for iori in range(1,63):
      self.p.postname='indent_fric0.3_R2.70_cA146.0_h0.320_ori%03i_OST_h19d.t16'%(iori)
      if os.path.exists(self.p.postname):
        self.marc2vtk(mode='fast', batchMode=1)
    
  def marc2vtk(self, label=None, mode='fast', 
               batchMode=0, 
               incRange=None,
               incStepMult=1.):
    if batchMode==0:
      try:
        self.p
      except:
        self.p=MARC_POST()
    ### ---- CHANGE dir/job/model to process here    
    os.chdir('M:/nicu');  
    jobname='ori001'
    self.p.postname='indent_fric0.3_R0.25_cA90.0_h0.010_4320els_'+jobname+'.t16'
    ### ----

    self.p.opent16()  
    self.p.t16info()

    incMax=self.p.p.increments(); 
    if incRange is None:
      incStep=5
      incRange=range(0,incMax+1,incStep)
    self.vtkPath=os.getcwd()+'/vtk_%s/'%self.p.postname
    if not os.path.exists(self.vtkPath):
      os.mkdir(self.vtkPath)
    for inc in incRange:
        print('Increment: %i'%inc)
        self.p.p.moveto(inc)
        t=self.p.p.time
        sys.stdout.write('inc:%i, time:%.3f\n'%(self.p.p.increment,t))
        self.incStr='inc%03i'%(inc*incStepMult)
        self.openFile(filename=self.vtkPath+self.p.postname[0:-4]+'_'+
        self.incStr+'.vtk')
        self.writeFirstLines(comment=self.p.postname,
          dtype='UNSTRUCTURED_GRID')
        self.p.writeNodes2VTK(fobj=self.f)
        self.p.writeElements2VTK(fobj=self.f)
        self.p.writeNodeScalars2VTK(fobj=self.f)
        self.p.writeElementData2VTK(fobj=self.f)
        # insert generation and writing of derived post values 
        # *here*
        
        self.f.close()
        print('Increment (self.p.p.increment): %i'%self.p.p.increment)
        print('Data written.')
    print(self.p.postname+' ready.')     
            
    
  def scaleBar(self, length=1.0, posXYZ=[0., 0., 0.]):
      self.fsb=open('micronbar_l%.1f.vtp'%length,'w+')
      self.writeFirstLines(self.fsb, comment='micronbar')
      pts=np.array([])
      width=length*1.
      height=length*1.
      wVec=np.array([0., width, 0.])
      lVec=np.array([length,0.,0.])
      hVec=np.array([0.,0.,height])
      posXYZ=posXYZ-0.5*wVec-0.5*lVec#-0.5*hVec # CENTERING Y/N
      posXYZ=np.array(posXYZ)
      pts=[posXYZ, posXYZ+lVec,
      posXYZ+wVec,
      posXYZ+wVec+lVec]
      pts.extend([pts[0]+hVec,pts[1]+hVec,pts[2]+hVec,pts[3]+hVec])
      print(len(pts), pts)
      self.fsb.write('POINTS %i float\n'%len(pts))
      for npts in range(0,len(pts)):
        self.fsb.write('%f %f %f\n'%(pts[npts][0], pts[npts][1], pts[npts][2]))
      if 1: #Triad
        nCells=3
        ptsPerCell=2 # Lines (Type=3)
        cellSize=(ptsPerCell+1)*nCells    
        self.fsb.write('CELLS %i %i\n'%(nCells,cellSize))
        self.fsb.write('2 0 1\n')     #X-Line
        self.fsb.write('2 0 2\n')     #Y-Line
        self.fsb.write('2 0 4\n')     #Z-Line
        self.fsb.write('CELL_TYPES %i\n'%(nCells))
        self.fsb.write('3\n3\n3\n')#Line
      else: # Cube, change posXYZ      
        nCells=1
        ptsPerCell=2 # Lines (Type=3)
        cellSize=(ptsPerCell+1)*nCells    
        self.fsb.write('CELLS %i %i\n'%(nCells,cellSize))
        self.fsb.write('2 0 1\n')     #Line
        self.fsb.write('CELL_TYPES %i\n'%(nCells))
        self.fsb.write('3\n')#Line
      
      self.fsb.write('\n')
      self.fsb.close()
      print(self.fsb)

  def example_unstructured(self):
    self.openFile(filename='example_unstructured_grid.vtk')
    self.f.write("""
# vtk DataFile Version 2.0
example_unstruct_grid
ASCII
    
POINTS 12 float
0 0 0
1 0 0
1 1 0
0 1 0
0 0 1
1 0 1
1 1 1
0 1 1
0 0 1.9
1 0 1.9
1 1 1.9
0 1 1.9


CELLS 2 18
8   0 1 2 3 4 5 6 7
8   4 5 6 
7 8 9 10 11

CELL_TYPES 2
12
12

POINT_DATA 12
SCALARS nodex float 1
LOOKUP_TABLE default
2.34E+12
2.00
0.00
1.62
5.03
1.02
1.50
0.00
3 5 6 23423423423423423423.23423423""")
    self.f.close()



  def writeNodes2VTK(self, fobj):
      self.VTKcnt=200 # how many numbers per line in vtk file
      fobj.write('POINTS %i'%self.p.nodes()+' float\n')
      for iNd in self.nodes:
        nd=self.p.node(iNd)
        disp=self.p.node_displacement(iNd)
        fobj.write('%f %f %f \n'%
        (nd.x+disp[0], nd.y+disp[1], nd.z+disp[2]))
      fobj.write('\n')  
      print('Nodes written to VTK: %i'%self.p.nodes())
      
  def writeElements2VTK(self, fobj):
      fobj.write('\nCELLS %i %i'%(self.p.elements(),self.p.elements()*9)+'\n')
      for iEl in range(0,self.nel):
        el=self.p.element(iEl)
        fobj.write('8 ')
        ndlist=el.items
        for k in [0, 1, 2, 3, 4, 5, 6, 7]: # FOR CELL TYPE VTK_HEXAHEDRON
         fobj.write('%6i '%(ndlist[k]-1))
        fobj.write('\n')
      fobj.write('\nCELL_TYPES %i'%self.p.elements()+'\n')
      cnt=0
      for iEl in range(0,self.nel):
        cnt=cnt+1
        fobj.write('12 ') #VTK_HEXAHEDRON
        if cnt>self.VTKcnt:
          fobj.write('\n');cnt=0
      fobj.write('\n')    
      print('Elements written to VTK: %i'%self.p.elements())
      
  def writeElScalars2NodesVTK(self,fobj):
      fobj.write('\nPOINT_DATA %i\n'%self.p.nodes())
      nScal=12
      nComponents=1+nScal
      fobj.write('SCALARS scalars float %i\n'%nComponents)
      fobj.write('LOOKUP_TABLE default\n')  
      idxScal=self.nscal_list.index('Displacement Z')    
      for iNd in self.nodes:
        fobj.write('%f '%(self.p.node_scalar(iNd,idxScal)))
        for iEl in range(0,self.nel):
          el=self.p.element(iEl)
          ndlist=el.items
          if (iNd+1) in ndlist:
            idx=ndlist.index(iNd+1)
            for iV in range(0,nScal):
              elData=self.p.element_scalar(iEl,35+iV)
              fobj.write('%f '%(elData[idx].value))
            break
        fobj.write('\n')  
      fobj.write('\n')  

  def writeNodeScalars2VTK(self,fobj):
      fobj.write('\nPOINT_DATA %i\n'%self.p.nodes())
      for idxNdScal in range(-3,self.nscals): # include node x,y,z
        if idxNdScal>=0:
          datalabel=self.nscal_list[idxNdScal]
          datalabel=re.sub("\s",'_',datalabel)
        else:
          if idxNdScal==-3: datalabel='node.x'  
          if idxNdScal==-2: datalabel='node.y'
          if idxNdScal==-1: datalabel='node.z'
        fobj.write('SCALARS %s float %i\n'%(datalabel,1))#nComponents))
        fobj.write('LOOKUP_TABLE default\n')  
        cnt=0
        for iNd in range(0,self.nnodes):
          cnt=cnt+1
          if idxNdScal>=0:
            ndData=self.p.node_scalar(iNd,idxNdScal)
          else:
             nd=self.p.node(iNd)
             if idxNdScal==-3: ndData=nd.x
             if idxNdScal==-2: ndData=nd.y
             if idxNdScal==-1: ndData=nd.z
          fobj.write('%E '%(ndData))
          if cnt>self.VTKcnt:
            fobj.write('\n')
            cnt=0
        fobj.write('\n')              
      fobj.write('\n')
        
  def writeElementData2VTK(self,fobj):
      fobj.write('\nCELL_DATA %i\n'%self.p.elements())
      for idxElScal in range(0,self.elscals):
        datalabel=self.elscal_list[idxElScal]
        datalabel=re.sub("\s",'_',datalabel)
        fobj.write('\n\nSCALARS %s float %i\n'%(datalabel,1))#nComponents))
        fobj.write('LOOKUP_TABLE default\n')  
        cnt=0
        for iEl in range(0,self.nel):
          cnt=cnt+1
          elData=self.p.element_scalar(iEl,idxElScal)
          avgScal=0.0
          if datalabel in ['phi1', 'PHI','phi2']:
            avgScal=avgScal+elData[0].value
          else:
            for IP in range(0,8):
              avgScal=avgScal+elData[IP].value
            avgScal=avgScal/8.
          fobj.write('%E '%(avgScal))
          if cnt>self.VTKcnt:
            fobj.write('\n')
            cnt=0
      fobj.write('\n')
      
      
  def example1(self):
    self.openFile()
    self.writeFirstLines()
    self.f.write("""DATASET POLYDATA
POINTS 8 float
0.0 0.0 0.0
1.0 0.0 0.0
1.0 1.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
1.0 0.0 1.0
1.0 1.0 1.0
0.0 1.0 1.0
POLYGONS 6 30
4 0 1 2 3
4 4 5 6 7
4 0 1 5 4
4 2 3 7 6
4 0 4 7 3
4 1 2 6 5

CELL_DATA 6
SCALARS cell_scalars int 1
LOOKUP_TABLE default
0
1
2
3
4
5
NORMALS cell_normals float
0 0 -1
0 0 1
0 -1 0
0 1 0
-1 0 0
1 0 0
FIELD FieldData 2
cellIds 1 6 int
0 1 2 3 4 5
faceAttributes 2 6 float
0.0 1.0 1.0 2.0 2.0 3.0 3.0 4.0 4.0 5.0 5.0 6.0

POINT_DATA 8
SCALARS sample_scalars float 1
LOOKUP_TABLE my_table
0.0
1.0
2.0
3.0
4.0
5.0
6.0
7.0
LOOKUP_TABLE my_table 8
0.0 0.0 0.0 1.0
1.0 0.0 0.0 1.0
0.0 1.0 0.0 1.0
1.0 1.0 0.0 1.0
0.0 0.0 1.0 1.0
1.0 0.0 1.0 1.0
0.0 1.0 1.0 1.0
1.0 1.0 1.0 1.0""")
    self.f.close()

    
import pyvtk
class marc_to_vtk():
    """
    Anybody wants to implement it with pyvtk?
    
    The advantage would be that pyvtk can also wirte the 
    <xml>-VTK format and binary. 
    These can be plotted with mayavi.
    """   

    def __init__(self):
      self.p=[]#MARC_POST() # self.p

    def run(self):      
      vtk = pyvtk.VtkData(\
          pyvtk.UnstructuredGrid(self.p.points,
                         hexahedron=self.p.cells),
          'm2v output')               
      vtk.tofile('m2v_file')

