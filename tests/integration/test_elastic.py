# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest
import numpy as np

import damask

@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'elastic'


def add_derived_quantities(fname):
    result = damask.Result(fname).view(increments=20)
    result.add_stress_second_Piola_Kirchhoff()
    result.add_strain(F='F',t='U',m=1)

    result.add_calculation('#F#-1','H')
    result.add_calculation('#S#[:,2,2]/#epsilon_U^1(F)#[:,2,2]','E','youngs modulus in z direction','Pa')

def Youngs_modulus(fname):
    result = damask.Result(fname).view(increments=20)
    return np.average(result.place('E'))

def strain_energy_density_ratio(fname):
    result = damask.Result(fname).view(increments=20)
    P      = np.average(result.place('P'),axis=0)
    H      = np.average(result.place('H'),axis=0)
    eps_GL = np.average(result.place('epsilon_U^1(F)'),axis=0)
    S      = np.average(result.place('S'),axis=0)
    return np.tensordot(S,eps_GL)/(H[2,2]*P[2,2])


@pytest.mark.parametrize('E,eulers',[( 673932.0e5, [  0.0000,  0.0000,  0.0000]),
                                        ( 706810.0e5, [270.0000, 9.00000, 90.0000]),
                                        ( 740384.0e5, [315.0000, 12.6253, 45.0000]),
                                        ( 810305.0e5, [270.0000, 18.0000, 90.0000]),
                                        ( 845783.0e5, [295.9873, 19.8733, 64.0127]),
                                        ( 956312.0e5, [315.0000, 24.6791, 45.0000]),
                                        ( 989373.0e5, [270.0000, 27.0000, 90.0000]),
                                        (1027028.0e5, [287.2677, 28.0832, 72.7323]),
                                        (1143138.0e5, [302.5253, 31.1449, 57.4747]),
                                        (1335172.0e5, [315.0000, 35.7757, 45.0000]),
                                        (1204766.0e5, [270.0000, 36.0000, 90.0000]),
                                        (1242176.0e5, [282.2979, 36.6347, 77.7021]),
                                        (1356802.0e5, [294.0948, 38.5158, 65.9052]),
                                        (1545069.0e5, [305.0420, 41.5859, 54.9580]),
                                        (1751587.0e5, [315.0000, 45.7767, 45.0000]),
                                        (1314037.0e5, [270.0000, 45.0000, 90.0000]),
                                        (1345101.0e5, [279.0000, 45.3549, 81.0000]),
                                        (1440960.0e5, [288.0000, 46.4370, 72.0000]),
                                        (1603736.0e5, [297.0000, 48.2988, 63.0000]),
                                        (1806050.0e5, [306.0000, 51.0266, 54.0000]),
                                        (1922800.0e5, [315.0000, 54.7356, 45.0000])])
@pytest.mark.parametrize('lattice',['cI','cF'])
def test_cubic(res_path,tmp_path,copy_files,assert_allclose,E,eulers,lattice):
    grid = 'singleCrystal'
    load = 'tensionZ'
    material = 'material'
    job = f'{grid}_{load}'

    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{grid}.msh'])

    config = damask.ConfigMaterial.load(res_path/'material.yaml')
    config['phase']['cubic']['lattice']=lattice
    config['material'][0]['constituents'][0]['phase']='cubic'
    config['material'][0]['constituents'][0]['O']=damask.Rotation.from_Euler_angles(eulers,degrees=True)
    config.save(tmp_path/f'{material}.yaml')

    damask.util.run(f'damask_mesh -l {load}.yaml -g {grid}.msh -m {material}.yaml -j {job}',
                    wd=tmp_path)
    add_derived_quantities(tmp_path/f'{job}.hdf5')

    assert 0.995 < strain_energy_density_ratio(tmp_path/f'{job}.hdf5') <= 1.0
    assert_allclose(E,Youngs_modulus(tmp_path/f'{job}.hdf5'),atol=0,rtol=5e-3)


@pytest.mark.parametrize('E,eulers',[(1432708.0e5, [  0.0000,  0.0000,  0.0000]),
                                        (1400580.0e5, [270.0000, 15.0000, 90.0000]),
                                        (1373331.0e5, [300.0000, 20.7536, 60.0000]),
                                        (1317927.0e5, [270.0000, 30.0000, 90.0000]),
                                        (1301514.0e5, [286.5974, 32.4767, 73.4026]),
                                        (1255337.0e5, [300.0000, 39.2315, 60.0000]),
                                        (1215974.0e5, [270.0000, 45.0000, 90.0000]),
                                        (1209335.0e5, [280.0000, 45.9930, 80.0000]),
                                        (1188944.0e5, [290.0000, 49.1066, 70.0000]),
                                        (1154240.0e5, [300.0000, 54.7356, 60.0000]),
                                        (1043717.0e5, [270.0000, 90.0000, 90.0000]),
                                        (1043717.0e5, [280.0000, 90.0000, 80.0000]),
                                        (1063321.0e5, [280.0000, 75.4892, 80.0000]),
                                        (1043717.0e5, [290.0000, 90.0000, 70.0000]),
                                        (1059631.0e5, [290.0000, 76.9357, 70.0000]),
                                        (1108077.0e5, [290.0000, 63.4349, 70.0000]),
                                        (1043717.0e5, [300.0000, 90.0000, 60.0000]),
                                        (1054472.0e5, [300.0000, 79.2714, 60.0000]),
                                        (1089112.0e5, [300.0000, 67.7923, 60.0000]),
                                        (1064653.0e5, [270.0000, 75.0000, 90.0000]),
                                        (1125068.0e5, [270.0000, 60.0000, 90.0000]),
                                        (1120691.0e5, [280.0000, 60.8526, 80.0000])])
def test_hex(res_path,tmp_path,copy_files,assert_allclose,E,eulers):
    grid = 'singleCrystal'
    load = 'tensionZ'
    material = 'material'
    job = f'{grid}_{load}'

    copy_files(res_path,tmp_path,[f'{load}.yaml',f'{grid}.msh'])

    config = damask.ConfigMaterial.load(res_path/'material.yaml')
    config['material'][0]['constituents'][0]['phase']='hex'
    config['material'][0]['constituents'][0]['O']=damask.Rotation.from_Euler_angles(eulers,degrees=True)
    config.save(tmp_path/f'{material}.yaml')

    damask.util.run(f'damask_mesh -l {load}.yaml -g {grid}.msh -m {material}.yaml -j {job}',wd=tmp_path)
    add_derived_quantities(tmp_path/f'{job}.hdf5')

    assert 0.995 < strain_energy_density_ratio(tmp_path/f'{job}.hdf5') <= 1.0
    assert_allclose(E,Youngs_modulus(tmp_path/f'{job}.hdf5'),atol=0,rtol=5e-3)
