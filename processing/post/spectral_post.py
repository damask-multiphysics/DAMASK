import binascii
print('post processing for mpie_spectral')

results = open('results.out', 'rb')

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
#
results.seek(begin+end+4)

header = results.read(4)
begin = results.tell()
end = results.read(1024).find(header) #find header (second time)
results.seek(begin+8) #position: header + string "workingdir "
jobname = results.read(end -8)
print('jobname:', jobname)

header = b'a'
results.seek(begin+end+19)
print(results.seek(begin+end+19))
begin = results.read(1024).find(header)
print(begin)
