#!/usr/bin/env python

# Makes postprocessing routines acessible from everywhere.
import os
from damask import Environment

damaskEnv = Environment()
baseDir = damaskEnv.relPath('processing/')
codeDir = damaskEnv.relPath('code/')

#define ToDo list
bin_link = { \
        'pre' : [
                'marc_addUserOutput.py',
                'mentat_pbcOnBoxMesh.py',
                'mentat_spectralBox.py',
                'OIMang_hex2cub.py',
                'patchFromReconstructedBoundaries.py',
                'seeds_fromRandom.py',
                'seeds_fromGeom.py',
                'seeds_check.py',
                'geom_fromAng.py',
                'geom_fromVPSC.py',
                'geom_fromEuclideanDistance.py',
                'geom_fromMinimalSurface.py',
                'geom_fromVoronoiTessellation.py',
                'geom_fromOsteonGeometry.py',
                'geom_canvas.py',
                'geom_check.py',
                'geom_rescale.py',
                'geom_pack.py',
                'geom_unpack.py',
                'geom_translate.py',
                'geom_vicinityOffset.py',
                'geom_grainGrowth.py',
                'geom_poke.py',
                ],
        'post' : [
                '3Dvisualize.py',
                'permuteData.py',
                'addCalculation.py',
                'addCauchy.py',
                'addCompatibilityMismatch.py',
                'addCurl.py',
                'addDeformedConfiguration.py',
                'addDeterminant.py',
                'addDeviator.py',
                'addDivergence.py',
                'addEhkl.py',
                'addEuclideanDistance.py',
                'addGrainID.py',
                'addMises.py',
                'addNorm.py',
                'addOrientations.py',
                'addIPFcolor.py',
                'addPK2.py',
                'addSchmidfactors.py',
                'addSpectralDecomposition.py',
                'addStrainTensors.py',
                'averageDown.py',
                'binXY.py',
                'blowUp.py',
                'stddevDown.py',
                'deleteColumn.py',
                'deleteInfo.py',
                'filterTable.py',
                'sortTable.py',
                'marc_deformedGeometry.py',
                'marc_extractData.py',
                'mentat_colorMap.py',
                'nodesFromCentroids.py',
                'perceptualUniformColorMap.py',
                'postResults.py',
                'showTable.py',
                'tagLabel.py',
                'vtk2ang.py',
                'vtk_addData.py',
                'vtk_pointcloud.py',
                'vtk_addPointcloudData.py',
                'vtk_voxelcloud.py',
                'vtk_addVoxelcloudData.py',
                ],
            }

root=os.access('/usr/local/bin', os.W_OK)
if root:
  binDir = '/usr/local/bin'
else:
  binDir = os.path.join(os.getenv('HOME'),'bin')
  if not os.path.isdir(binDir):
    os.mkdir(binDir)
            
for dir in bin_link:
  for file in bin_link[dir]:
    src = os.path.abspath(os.path.join(baseDir,dir,file))
    if (file == ''):
      sym_link = os.path.abspath(os.path.join(binDir,dir))
    else:
      sym_link = os.path.abspath(os.path.join(binDir,os.path.splitext(file)[0]))
    print sym_link,'-->',src
    if os.path.lexists(sym_link):
      os.remove(sym_link)    
    os.symlink(src,sym_link)            



