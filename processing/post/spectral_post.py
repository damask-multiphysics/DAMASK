print('post processing for mpie_spectral')

results = open('results.out', 'rb')

header = results.read(4)
end = results.read(1024).find(header)
results.seek(4)
loadcase = results.read(end)
print('loadcase:', loadcase)
#
begin = end
results.seek(end+4)
header = results.read(4)
print(header)
end = results.read(1024).find(header)
#results.seek(4)
#workingdir = results.read(end)
#print('workingdir:', workingdir)
