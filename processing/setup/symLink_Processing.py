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
                'randomSeeding.py',
                'geom_fromAng.py',
                'geom_fromVPSC.py',
                'geom_fromMinimalSurface.py',
                'geom_fromVoronoiTessellation.py',
                'geom_Osteon.py',
                'geom_canvas.py',
                'geom_check.py',
                'geom_rescale.py',
                'geom_pack.py',
                'geom_unpack.py',
                'geom_translate.py',
                'geom_vicinityOffset.py',
                'geom_euclideanDistance.py'
                ],
        'post' : [
                '3Dvisualize.py',
                'addCalculation.py',
                'addCauchy.py',
                'addCompatibilityMismatch.py',
                'addCurl.py',
                'addDataToGeometry.py',
                'addDeformedConfiguration.py',
                'addDeterminant.py',
                'addDeviator.py',
                'addDivergence.py',
                'addEhkl.py',
                'addEuclideanDistance.py',
                'addMises.py',
                'addNorm.py',
                'addPK2.py',
                'addSpectralDecomposition.py',
                'addStrainTensors.py',
                'averageDown.py',
                'binXY.py',
                'blowUp.py',
                'stddevDown.py',
                'deleteColumn.py',
                'deleteInfo.py',
                'filterTable.py',
                'marc_deformedGeometry.py',
                'marc_extractData.py',
                'mentat_colorMap.py',
                'nodesFromCentroids.py',
                'perceptualUniformColorMap.py',
                'postResults.py',
                'showTable.py',
                'table2ang',
                'tagLabel.py',
                ],
            }
            
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



