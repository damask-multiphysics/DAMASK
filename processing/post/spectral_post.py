import array
import struct
print('post processing for mpie_spectral')

filename ='results.out'
results = open(filename, 'rb')
print('filename:', filename)

#find load case name
searchstring = b'Loadcase '

begin = results.read(2048).find(searchstring) + 9 # function will return -1 if string isnt found
results.seek(begin -13)  #position of header
header = results.read(4) #store header
results.seek(begin,0)
length = results.read(2048).find(header)
results.seek(begin,0)

loadcase = results.read(length)
print('loadcase:', loadcase)


#find workingdir name
results.seek(0,0)
searchstring = b'Workingdir '

begin = results.read(2048).find(searchstring) + 11 
results.seek(begin -15)  #position of header
header = results.read(4) #store header
results.seek(begin,0)
length = results.read(2048).find(header)
results.seek(begin,0)

workingdir = results.read(length)
print('workingdir:', workingdir)

#find job name
results.seek(0,0)
searchstring = b'JobName '

begin = results.read(2048).find(searchstring) + 8 
results.seek(begin -12)  #position of header
header = results.read(4) #store header
results.seek(begin,0)
length = results.read(2048).find(header)
results.seek(begin,0)

jobname = results.read(length)
print('jobname:', jobname)

#find size of materialpoint results
results.seek(0,0)
searchstring = b'materialpoint_sizeResults'

begin = results.read(2048).find(searchstring) + 24 
results.seek(begin,0)
materialpoint_sizeResults = array.array('i',[0])
materialpoint_sizeResults[0] = struct.unpack('i',results.read(4))[0]
print('materialpoint_sizeResults:', materialpoint_sizeResults)

#read in resolution (FPs) of geometry
resolution = array.array('i',[0,0,0])
results.seek(0,0)
searchstring = b'resolution'
begin = results.read(2048).find(searchstring) + 10 #end position of string "resolution"
results.seek(begin,0)
pos = results.read(60).find(b'a')
results.seek(begin+pos+1,0)
resolution[0]=struct.unpack('i',results.read(4))[0]
results.seek(begin,0)
pos = results.read(60).find(b'b')
results.seek(begin+pos+1,0)
resolution[1]=struct.unpack('i',results.read(4))[0]
results.seek(begin,0)
pos = results.read(60).find(b'c')
results.seek(begin+pos+1,0)
resolution[2]=struct.unpack('i',results.read(4))[0]

print('resolution:',resolution)

#read in dimension of geometry
geomdimension = array.array('d',[0,0,0])
results.seek(0,0)
searchstring = b'geomdimension'
begin = results.read(2048).find(searchstring) + 13
results.seek(begin,0)
pos = results.read(60).find(b'x')
results.seek(begin+pos+1,0)
geomdimension[0]=struct.unpack('d',results.read(8))[0]
results.seek(begin,0)
pos = results.read(60).find(b'y')
results.seek(begin+pos+1,0)
geomdimension[1]=struct.unpack('d',results.read(8))[0]
results.seek(begin,0)
pos = results.read(60).find(b'z')
results.seek(begin+pos+1,0)
geomdimension[2]=struct.unpack('d',results.read(8))[0]

print('geomdimension:',geomdimension)

