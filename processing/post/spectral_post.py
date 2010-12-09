import array
import struct
#import numpy
print('post processing for mpie_spectral')

filename ='results.out'
results = open(filename, 'rb')
print('filename:', filename)
position=0

#this funtion finds a string value in the given file for a given keyword
def searchNameForKeyword(searchstring, file, maxpos):
    results.seek(0,0)
    begin = file.read(2048).find(searchstring) + len(searchstring) # function will return -1 if string isnt found
    file.seek(begin -len(searchstring) -4)  #position of header
    header = file.read(4) #store header
    file.seek(begin,0)
    length = file.read(2048).find(header)
    file.seek(begin,0)
    name = file.read(length)
    return name, max(maxpos,results.tell())

#this funtion finds an integer value in the given file for a given keyword
def searchValueForKeyword(searchstring, file, maxpos):
    results.seek(0,0)
    begin = file.read(2048).find(searchstring) + len(searchstring) # function will return -1 if string isnt found
    file.seek(begin,0)
    value = array.array('i',[0])
    value = struct.unpack('i',results.read(4))[0]
    return value, max(maxpos,results.tell())

#finds the value for the three integers followed the given keyword
def searchArrayForKeyword(searchstring, file, maxpos):
    values = array.array('i',[0,0,0])
    results.seek(0,0)
    begin = results.read(2048).find(searchstring) + len(searchstring) #end position of string "resolution"
    pos = results.read(60).find(b'a')
    results.seek(begin+pos+2,0)  #why two, not 1??
    values[0]=struct.unpack('i',results.read(4))[0]
    maxpos=max(maxpos,results.tell())
    results.seek(begin,0)
    pos = results.read(60).find(b'b')
    results.seek(begin+pos+1,0)
    values[1]=struct.unpack('i',results.read(4))[0]
    maxpos=max(maxpos,results.tell())
    results.seek(begin,0)
    pos = results.read(60).find(b'c')
    results.seek(begin+pos+1,0)
    values[2]=struct.unpack('i',results.read(4))[0]
    maxpos=max(maxpos,results.tell())
    return values, maxpos

#finds the value for the three doubles followed the given keyword
def searchFarrayForKeyword(searchstring, file, maxpos):
    values = array.array('d',[0,0,0])
    results.seek(0,0)
    begin = results.read(2048).find(searchstring) + len(searchstring) #end position of string "resolution"
    pos = results.read(60).find(b'x')
    results.seek(begin+pos+2,0)  #why two, not 1??
    values[0]=struct.unpack('d',results.read(8))[0]
    maxpos=max(maxpos,results.tell())
    results.seek(begin,0)
    pos = results.read(60).find(b'y')
    results.seek(begin+pos+1,0)
    values[1]=struct.unpack('d',results.read(8))[0]
    maxpos=max(maxpos,results.tell())
    results.seek(begin,0)
    pos = results.read(60).find(b'z')
    results.seek(begin+pos+1,0)
    values[2]=struct.unpack('d',results.read(8))[0]
    maxpos=max(maxpos,results.tell())
    return values, maxpos

loadcase, position = searchNameForKeyword(b'Loadcase ', results, position)
workingdir, position = searchNameForKeyword(b'Workingdir ', results, position)
jobname, position = searchNameForKeyword(b'JobName ', results, position)
totalincs, position =  searchValueForKeyword(b'totalincs ', results, position)
materialpoint_sizeResults, position =  searchValueForKeyword(b'materialpoint_sizeResults ', results, position)
resolution, position = searchArrayForKeyword(b'resolution ', results, position)
geomdimension, position = searchFarrayForKeyword(b'geomdimension ', results, position)
position=position+4 +13*8
results.seek(position,0)
tnsr = array.array('d',[0,0,0,0,0,0,0,0,0])
tnsr[0]=struct.unpack('d',results.read(8))[0]
tnsr[1]=struct.unpack('d',results.read(8))[0]
tnsr[2]=struct.unpack('d',results.read(8))[0]
tnsr[3]=struct.unpack('d',results.read(8))[0]
tnsr[4]=struct.unpack('d',results.read(8))[0]
print(tnsr)


# ended reading of header. now starting to read information concerning output
#filename = jobname[0:len(jobname)-5]+b'.outputCrystallite'
#outputCrystallite = open(filename, 'r')

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


