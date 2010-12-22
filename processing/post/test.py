import numpy
import reconstruct
# dummy values for resolution and geomdimension (formerly called meshdimension)
resolution=numpy.zeros((3),'i')
resolution[0]=2
resolution[1]=2
resolution[2]=2
geomdimension=numpy.zeros((3),'d')
geomdimension[0]=1.0
geomdimension[1]=1.0
geomdimension[2]=1.0
defgrad=numpy.zeros((resolution[0],resolution[1],resolution[2],3,3),'d')
for i in range (0, resolution[0]):
	for j in range (0, resolution[1]):
		for k in range (0, resolution[2]):
			defgrad[i][j][k][0][0]=1.0
			defgrad[i][j][k][1][1]=1.0
			defgrad[i][j][k][2][2]=1.0
print defgrad
current_configuration=reconstruct.simple(defgrad,resolution[0],resolution[1],resolution[2],geomdimension) # don't know how to pass arrays for the resolution
print current_configuration