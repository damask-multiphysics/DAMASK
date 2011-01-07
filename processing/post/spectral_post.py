#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

# This script is used for the post processing of the results achieved by the spectral method.
# As it reads in the data coming from "materialpoint_results", it can be adopted to the data
# computed using the FEM solvers. Until now, its capable to handle elements with one IP in a regular order
# written by M. Diehl, m.diehl@mpie.de

import array
import struct
import numpy
import reconstruct

#this funtion finds a string value in the given file for a given keyword, also returns the position
def searchNameForKeyword(searchstring, file):
    file.seek(0,0)
    begin = file.read(2048).find(searchstring) + len(searchstring) # function will return -1 if string is not found
    file.seek(begin -len(searchstring) -4)                         # position of header (Fortran specific, 4 bit long on msuws4
    header = file.read(4)                                          # store header. header will be repeated, thus it allows us to determine the length of the searchstring
    file.seek(begin,0)                                             # go back to starting postion
    length = file.read(2048).find(header)                          # find header(footer) a second time
    file.seek(begin,0) 
    return file.read(length), results.tell()                       # return searchstring and maximum position in file

#this funtion finds an integer value in the given file for a given keyword, works similar as "searchNameForKeyword"
def searchValueForKeyword(searchstring, file):
    results.seek(0,0)
    begin = file.read(2048).find(searchstring) + len(searchstring) 
    file.seek(begin,0)                                                                     
    value = struct.unpack('i',results.read(4))[0]
    return value, results.tell()

#finds the value for the three integers (a,b,c) following the given keyword
def searchIntegerArrayForKeyword(searchstring, file):
	values = array.array('i',[0,0,0])
	file.seek(0,0)
	begin = file.read(2048).find(searchstring) + len(searchstring) #end position of string "searchstring"
	pos = file.read(60).find('a')
	file.seek(begin+pos+2,0)                                       #why 2, not 1??
	values[0]=struct.unpack('i',file.read(4))[0]
	maxpos=file.tell()
	file.seek(begin,0)
	pos = file.read(60).find('b')
	file.seek(begin+pos+1,0)
	values[1]=struct.unpack('i',file.read(4))[0]
	maxpos = max(maxpos,file.tell())
	file.seek(begin,0)
	pos = file.read(60).find('c')
	file.seek(begin+pos+1,0)
	values[2]=struct.unpack('i',file.read(4))[0]
	maxpos=max(maxpos,file.tell())
	return values, maxpos

#finds the value for the three doubles (x,y,z) following the given keyword
def searchDoubleArrayForKeyword(searchstring, file):
	values = array.array('d',[0,0,0])
	file.seek(0,0)
	begin = file.read(2048).find(searchstring) + len(searchstring) 
	pos = file.read(60).find('x')
	file.seek(begin+pos+2,0) 
	values[0]=struct.unpack('d',file.read(8))[0]
	maxpos=file.tell()
	file.seek(begin,0)
	pos = file.read(60).find('y')
	file.seek(begin+pos+1,0)
	values[1] = struct.unpack('d',file.read(8))[0]
	maxpos = max(maxpos,file.tell())
	file.seek(begin,0)
	pos = file.read(60).find('z')
	file.seek(begin+pos+1,0)
	values[2] = struct.unpack('d',file.read(8))[0]
	maxpos = max(maxpos,file.tell())
	return values, maxpos
print '*********************************************************************************'	
print 'Post Processing for Material subroutine for BVP solution using spectral method'
print '*********************************************************************************\n'

#reading in the header of the results file
filename ='32x32x32x100.out'
print 'Results Filename:', filename 
results = open(filename, 'rb')

loadcase, position = searchNameForKeyword('Loadcase', results)
workingdir, temp = searchNameForKeyword('Workingdir', results)
position = max(temp, position)
jobname, temp = searchNameForKeyword('JobName', results)
position = max(temp, position)
totalincs, temp =  searchValueForKeyword('totalincs', results)
position = max(temp, position)
materialpoint_sizeResults, temp =  searchValueForKeyword('materialpoint_sizeResults', results)
position = max(temp, position)
resolution, temp = searchIntegerArrayForKeyword('resolution', results)
position = max(temp, position)
geomdimension, temp = searchDoubleArrayForKeyword('geomdimension', results)

print 'Load case:', loadcase
print 'Workingdir:', workingdir
print 'Job Name:', jobname
print 'Total No. of Increments:', totalincs 
print 'Materialpoint_sizeResults:', materialpoint_sizeResults
print 'Resolution:', resolution
print 'Geomdimension', geomdimension 
print 'Position in File:', position

# Ended reading of header
# Now starting to read information concerning output of materialpoint_results(:,1,:)

filename = jobname[0:len(jobname)-5]+'.outputCrystallite'
outputCrystallite = open(filename, 'r')


# some funtions that might be useful
def InCharCount(location, character):
    subj = file(location, "r")
    body = subj.read()
    subj.close()
    return body.count(character)


def InCharCount(location, character):
    subj = file(location, "r")

    nbr_of_char = 0
    for line in subj:
        nbr_of_char = nbr_of_char + line.count(character)

    return nbr_of_char


