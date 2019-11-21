import numpy as np
from damask import Rotation

class TestRotation:
   
    n = 1000

    def test_Eulers(self):
        for r in range(self.n):
            rot = Rotation.fromRandom()
            assert  np.allclose(rot.asQuaternion(), 
                                Rotation.fromEulers(rot.asEulers()).asQuaternion())


    def test_AxisAngle(self):
        for r in range(self.n):
            rot = Rotation.fromRandom()
            assert  np.allclose(rot.asEulers(),
                                Rotation.fromAxisAngle(rot.asAxisAngle()).asEulers())


    def test_Matrix(self):
        for r in range(self.n):
            rot = Rotation.fromRandom()
            assert  np.allclose(rot.asAxisAngle(),
                                Rotation.fromMatrix(rot.asMatrix()).asAxisAngle())


    def test_Rodriques(self):
        for r in range(self.n):
            rot = Rotation.fromRandom()
            assert  np.allclose(rot.asMatrix(),
                                Rotation.fromRodrigues(rot.asRodrigues()).asMatrix())


    def test_Homochoric(self):
        for r in range(self.n):
            rot = Rotation.fromRandom()
            assert  np.allclose(rot.asRodrigues(),
                                Rotation.fromHomochoric(rot.asHomochoric()).asRodrigues())


    def test_Cubochoric(self):
        for r in range(self.n):
            rot = Rotation.fromRandom()
            assert  np.allclose(rot.asHomochoric(),
                                Rotation.fromCubochoric(rot.asCubochoric()).asHomochoric())


    def test_Quaternion(self):
        for r in range(self.n):
            rot = Rotation.fromRandom()
            assert  np.allclose(rot.asCubochoric(),
                                Rotation.fromQuaternion(rot.asQuaternion()).asCubochoric())
