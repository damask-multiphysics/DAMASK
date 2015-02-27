#!/usr/bin/env python

##
# This script will read in all the seeds and partition the space
# using scipy.spatial.Delaunay triangulation.
# The unknown location will be then interpolated through Barycentric
# interpolation method, which relies on the triangulation.
# A rim will be automatically added to the patch, which will help
# improve the compatibility with the spectral solver as well as
# maintain meaningful microstructure(reduce artifacts).


import sys, os
import numpy as np
import argparse
from scipy.spatial import Delaunay

scriptID   = string.replace('$Id$','\n','\\n')
scriptName = os.path.splitext(scriptID.split()[1])[0]


OFFSET = 0.1  #resize the seeded volume to give space for rim/pan
PHANTOM_ID = -1  #grain ID for phantom seeds

def d_print(info, data, separator=False):
    '''quickly print debug information'''
    if(separator): print "*"*80
    print info
    print data


#prepare command line interface
parser = argparse.ArgumentParser(prog="geoFromBarycentic",
                                 description='''Generate geom file through \
                                   Barycentric interpolating seeds file.''',
                                 epilog="requires numpy, scipy, and \
                                   ezxtal(https://github.com/KedoKudo/ezxtal.git).")
parser.add_argument("seeds",
                    help="seeds file in DAMASK format:\
                          http://damask.mpie.de/Documentation/AsciiTableFormat",
                    default="test.seeds")
parser.add_argument("-v", "--version",
                    action="version",
                    version="%(prog)s 0.1")
parser.add_argument("-g", "--grid",
                    nargs=3,
                    help="grid size(mesh resolution, recommend using 2^x)",
                    default=[32,32,32],
                    type=int)
parser.add_argument("-s", "--size",
                    help="physical size of the target volume.",
                    nargs=3,
                    default=[1.0,1.0,1.0],
                    type=float)
parser.add_argument("-o", "--origin",
                    help="lower left corner of the patch.",
                    nargs=3,
                    default=[0.0,0.0,0.0],
                    type=float)
parser.add_argument('-m', '--homogenization',
                    help='homogenization index to be used',
                    default=1,
                    type=int)
parser.add_argument('-c', '--crystallite',
                    help='crystallite index to be used',
                    default=1,
                    type=int)
parser.add_argument('-p', '--phase',
                    help='phase index to be used',
                    default=1,
                    type=int)
parser.add_argument('-F', '--Favg',
                    help='reshape the periodicity, not useful for RIM method',
                    nargs=9,
                    default=[1.0,0.0,0.0,
                             0.0,1.0,0.0,
                             0.0,0.0,1.0],
                    type=float)
parser.add_argument("-G", "--geomFile",
                    help='the name of the output geom file',
                    default='seeds.geom',
                    type=str)
parser.add_argument("-C", "--configFile",
                    help='output dummy material.config file',
                    action='store_true',
                    default=False)
parser.add_argument("-d", "--debug",
                    help="start debugging script",
                    action='store_true',
                    default=False)
parser.add_argument("-S", "--seedsFile",
                    help="write out resized seeds file",
                    action='store_true',
                    default=False)
parser.add_argument("-r", '--addRim',
                    help="add rim and provide control of face lifting point",
                    action='store_true',
                    default=False)
args = parser.parse_args()  # get all the arguments right after

#quick help to user
print "*"*80
parser.print_help()
print '''Sample usage:
./geoFromBarycentic.py 20grains.seeds -g 128 128 128 -S -r; geom_check seeds.geom; seeds_check new_seed.seeds.
'''
print "*"*80
if (args.debug):
    d_print("args are:", parser.parse_args(),separator=True)

#/\/\/\/\/#
# m a i n #
#\/\/\/\/\#
print "only work for 3D case now, 2D support coming soon..."
print "reading seeds file: {}".format(args.seeds)

with open(args.seeds, 'r') as f:
    rawtext = f.readlines()
    n_header = int(rawtext.pop(0).split()[0])
    #record all the seeds position
    if (args.addRim):
        grid_shift = np.array(args.size) * np.array([OFFSET,OFFSET,OFFSET*2])
        s_coords = np.array([[np.array(float(item))*(1 - OFFSET*2)
                                for item in line.split()[:3]] + grid_shift
                                for line in rawtext[n_header:]])
    else:
        #no need for shifting with periodicity
        s_coords = np.array([[np.array(float(item))
                                for item in line.split()[:3]]
                                for line in rawtext[n_header:]])

    #record ID of the seeds: int/EulerAngles
    if 'microstructure' in rawtext[n_header-1]:
        s_id = [int(line.split()[-1]) for line in rawtext[n_header:]]
    else:
        print "WARNING:"
        print "THIS SCRIPT DOES NOT UNDERSTAND HOW TO GROUP CRYSTALLITES."
        print "ALL CRYSTAL ORIENTATIONS ARE CONSIDERED TO BE UNIQUE."
        print "FOR MORE ACCURATE CONTROL OF SEEDS GROUPING, USE MICROSTRUCTURE ID."
        s_id = range(len(s_coords))
        #s_eulers here is just a quick book keeping
        s_eulers = np.array([[float(item) for item in line.split()[3:]]
                                          for line in rawtext[n_header:]])

if(args.debug):
    print d_print("resize point cloud to make space for rim/pan:",
                  s_coords)

if(args.addRim):
    #add binding box to create rim/pan for the volume where the ID of the seeds is
    #unknown
    print "Shrining the seeds to {}x in each direction".format(1 - OFFSET*2)
    x,y,z = args.size[0],args.size[1],args.size[2]
    print "Use circumscribed sphere to place phantom seeds."
    r = np.sqrt(x**2+y**2+z**2)/2.0
    BINDBOX = [[0,0,0],[x,0,0],[0,y,0],[x,y,0],
               [0,0,z],[x,0,z],[0,y,z],[x,y,z],
               [x/2.0+r,y/2,     z/2],     [x/2.0-r, y/2,     z/2],
               [x/2,    y/2.0+r, z/2],     [x/2,     y/2.0-r, z/2],
               [x/2,    y/2,     z/2.0-r]]  #8 corners + 5 face centers (no top)
    print "Adding phantom seeds for RIM generation:"
    for point in BINDBOX:
        print point
        s_coords = np.vstack([s_coords,point])
        s_id.append(PHANTOM_ID)
else:
    #The idea here is that we read in each seed point, than duplicate in 3D (make a few copies),
    #move on to the next seed point, repeat the same procedure. As for the ID list, we can just use the
    #same one. The trick here is use the floor division to find the correct id since we pretty much duplicate
    #the same point several times.
    Favg = np.array(args.Favg).reshape((3,3))
    x,y,z = args.size[0],args.size[1],args.size[2]
    tmp = []
    for seed in s_coords:
        tmp += [np.dot(Favg, np.array(seed) + np.array([dx,dy,dz]))
            for dz in [-z, 0, z]
            for dy in [-y, 0, y]
            for dx in [-x, 0, x]]
    s_coords = tmp
    for item in tmp:
        print item

if (args.seedsFile):
    with open("new_seed.seeds", "w") as f:
        outstr = "4\theader\n"
        outstr += "grid\ta {}\tb {}\tc {}\n".format(args.grid[0],
                                                    args.grid[1],
                                                    args.grid[2])
        outstr += "microstructures {}\n".format(len(set(s_id)))
        outstr += "randomSeed 0\n"
        outstr += "x\ty\tz\tmicrostructure"
        if (args.addRim):
            for i in range(len(s_id)):
                outstr += "{}\t{}\t{}\t{}\n".format(s_coords[i][0],
                                                    s_coords[i][1],
                                                    s_coords[i][2],
                                                    s_id[i])
        else:
            for i in range(len(s_coords)):
                outstr += "{}\t{}\t{}\t{}\n".format(s_coords[i][0],
                                                    s_coords[i][1],
                                                    s_coords[i][2],
                                                    s_id[i//3**3])
        f.write(outstr)

#triangulate space with given point-clouds
tri = Delaunay(s_coords)

if(args.debug):
    d_print("simplices:", tri.simplices, separator=True)
    d_print("vertices:", s_coords[tri.simplices])

#populate grid points (only 3D for now)
mesh_pts = [[(i+0.5)*args.size[0]/args.grid[0],
             (j+0.5)*args.size[1]/args.grid[1],
             (k+0.5)*args.size[2]/args.grid[2]]
                for k in range(args.grid[2])
                for j in range(args.grid[1])
                for i in range(args.grid[0])]

mesh_ids = [PHANTOM_ID*2]*len(mesh_pts)  #initialize grid

#search ID for each grid point
s_id = np.array(s_id)  #allow multi-indexing
mesh_idx = tri.find_simplex(mesh_pts)

for i, pt in enumerate(mesh_pts):
    if mesh_idx[i] < 0:
        continue  #didn't find any envelop tetrahedron --> something wrong!
    #calculate Barycentric coordinates
    bary_c = tri.transform[mesh_idx[i],:3,:3].dot(pt-tri.transform[mesh_idx[i],3,:])
    bary_c = np.append(bary_c, 1 - bary_c.sum())

    if (args.addRim):
        tmp_ids = s_id[tri.simplices[mesh_idx[i]]]  #rim method
    else:
        tmp_ids = np.array(s_id[tri.simplices[mesh_idx[i]]//(3**3)])  #kill periodicity through floor division
        #print tmp_ids
        #print tri.simplices[mesh_idx[i]]//(3**3)

    max_weight = -1960
    for this_id in tmp_ids:
        msk = [item==this_id for item in tmp_ids]  #find vertex with the same id
        tmp_weight = sum([bary_c[j] for j in range(len(bary_c)) if msk[j]])
        if tmp_weight > max_weight:
            max_weight = tmp_weight
            mesh_ids[i] = this_id
    if (args.debug):
        d_print("bary_c:",bary_c,separator=True)
        d_print("vertex ID:", tmp_ids)
        d_print("final ID:", mesh_ids[i])

mesh_ids = np.reshape(mesh_ids, (-1, args.grid[0]))

#write to file
with open(args.geomFile, "w") as f:
    outstr = "5\theader\n"
    outstr += "grid\ta {}\tb {}\tc {}\n".format(args.grid[0],
                                                args.grid[1],
                                                args.grid[2])
    outstr += "size\tx {}\ty {}\tz {}\n".format(args.size[0],
                                                args.size[1],
                                                args.size[2])
    outstr += "origin\tx {}\ty {}\tz {}\n".format(args.origin[0],
                                                  args.origin[1],
                                                  args.origin[2])
    outstr += "homogenization\t{}\nmicrostructure\t{}\n".format(args.homogenization,
                                                                len(set(s_id)))
    for row in mesh_ids:
        row = [str(item) for item in list(row)]
        outstr += "\t".join(row) + "\n"
    f.write(outstr)