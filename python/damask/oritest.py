import damask
import numpy as np


rot0= damask.Rotation.from_random()
rot1= damask.Rotation.from_random()
rot2= damask.Rotation.from_random()

ori0=damask.Orientation(rot0,'fcc')
ori1=damask.Orientation(rot1,'fcc')
ori2=damask.Orientation(rot2,'fcc')




quat=np.array([rot0.as_quaternion(),rot1.as_quaternion(),rot2.as_quaternion()])
rot=damask.Rotation.from_quaternion(quat)

ori=damask.Orientation(rot,'fcc')

ori.equivalent()



# doesn't work this way, don't know why
#ori.equivalent()[:,0][0] == ori0.equivalentOrientations()[0]

for s in range(24):
    print(ori.equivalent()[s,0].rotation.as_Eulers() == ori0.equivalentOrientations()[s].rotation.as_Eulers())
    print(ori.equivalent()[s,1].rotation.as_Eulers() == ori1.equivalentOrientations()[s].rotation.as_Eulers())
    print(ori.equivalent()[s,2].rotation.as_Eulers() == ori2.equivalentOrientations()[s].rotation.as_Eulers())


