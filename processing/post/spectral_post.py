import array
import struct
print('post processing for mpie_spectral')

filename ='results.out'
results = open(filename, 'rb')
print('filename:', filename)
header = results.read(4) #read header
begin = results.tell()

end = results.read(1024).find(header) #find header (second time)

results.seek(begin+9) #position: header + string "loadcase "
#loadcase = binascii.b2a_base64(results.read(end-9)) #conversion doesnt work
loadcase = results.read(end -9)
print('loadcase:', loadcase)


results.seek(begin+end+4)

header = results.read(4)
begin = results.tell()
end = results.read(1024).find(header) #find header (second time)
results.seek(begin+11) #position: header + string "workingdir "
workingdir = results.read(end -11)
print('workingdir:', workingdir)

results.seek(begin+end+4)

header = results.read(4)
begin = results.tell()
end = results.read(1024).find(header) #find header (second time)
results.seek(begin+8) #position: header + string "workingdir "
jobname = results.read(end -8)
print('jobname:', jobname)

results.seek(results.tell()+4)
header = results.read(4)


resolution = array.array('i',[0,0,0])

header1 = b'a'
begin = results.seek(results.tell()+11)
begin = begin + results.read(1024).find(header1)
results.seek(begin+1)
resolution[0]=struct.unpack('i',results.read(4))[0]

print(resolution)

