import pytest
import numpy as np
from pathlib import Path

import damask
from damask import Table
from damask import Crystal
from damask import util


@pytest.fixture
def res_path(res_path_base):
    """Directory containing testing resources."""
    return res_path_base/'Crystal'

@pytest.mark.parametrize('lattice,family',[('aP','cubic'),('xI','cubic')])
def test_invalid_init(lattice,family):
    with pytest.raises(KeyError):
        Crystal(family=family,lattice=lattice)

def test_eq(np_rng):
    family = np_rng.choice(list(damask._crystal.lattice_symmetries.values()))
    assert Crystal(family=family) == Crystal(family=family)

def test_double_to_lattice():
    c = Crystal(lattice='cF')
    with pytest.raises(KeyError):
        c.to_lattice(direction=np.ones(3),plane=np.ones(3))

def test_double_to_frame():
    c = Crystal(lattice='cF')
    with pytest.raises(KeyError):
        c.to_frame(uvw=np.ones(3),hkl=np.ones(3))


@pytest.mark.parametrize('lattice,a,b,c,alpha,beta,gamma',
                        [
                         ('aP',0.5,2.0,3.0,0.8,0.5,1.2),
                         ('mP',1.0,2.0,3.0,np.pi/2,0.5,np.pi/2),
                         ('oI',0.5,1.5,3.0,np.pi/2,np.pi/2,np.pi/2),
                         ('tP',0.5,0.5,3.0,np.pi/2,np.pi/2,np.pi/2),
                         ('hP',1.0,None,1.6,np.pi/2,np.pi/2,2*np.pi/3),
                         ('cF',1.0,1.0,None,np.pi/2,np.pi/2,np.pi/2),
                        ])
def test_bases_contraction(lattice,a,b,c,alpha,beta,gamma):
    c = Crystal(lattice=lattice,
                a=a,b=b,c=c,
                alpha=alpha,beta=beta,gamma=gamma)
    assert np.allclose(np.eye(3),np.einsum('ik,jk',c.basis_real,c.basis_reciprocal))

def test_basis_invalid():
    with pytest.raises(KeyError):
        Crystal(family='cubic').basis_real

def test_basis_real(np_rng):
    for gamma in np_rng.random(2**8)*np.pi:
        basis = np.tril(np_rng.random((3,3))+1e-6)
        basis[1,:2] = basis[1,1]*np.array([np.cos(gamma),np.sin(gamma)])
        basis[2,:2] = basis[2,:2]*2-1
        lengths = np.linalg.norm(basis,axis=-1)
        cosines = np.roll(np.einsum('ij,ij->i',basis,np.roll(basis,1,axis=0))/lengths/np.roll(lengths,1),1)
        o = Crystal(lattice='aP',
                    **dict(zip(['a','b','c'],lengths)),
                    **dict(zip(['alpha','beta','gamma'],np.arccos(cosines))),
                    )
        assert np.allclose(o.to_frame(uvw=np.eye(3)),basis,rtol=1e-4), 'Lattice basis disagrees with initialization'

@pytest.mark.parametrize('keyFrame,keyLattice',[('uvw','direction'),('hkl','plane'),])
@pytest.mark.parametrize('vector',np.array([
                                            [1.,1.,1.],
                                            [-2.,3.,0.5],
                                            [0.,0.,1.],
                                            [1.,1.,1.],
                                            [2.,2.,2.],
                                            [0.,1.,1.],
                                            ]))
@pytest.mark.parametrize('lattice,a,b,c,alpha,beta,gamma',
                        [
                         ('aP',0.5,2.0,3.0,0.8,0.5,1.2),
                         ('mP',1.0,2.0,3.0,np.pi/2,0.5,np.pi/2),
                         ('oI',0.5,1.5,3.0,np.pi/2,np.pi/2,np.pi/2),
                         ('tP',0.5,0.5,3.0,np.pi/2,np.pi/2,np.pi/2),
                         ('hP',1.0,1.0,1.6,np.pi/2,np.pi/2,2*np.pi/3),
                         ('cF',1.0,1.0,1.0,np.pi/2,np.pi/2,np.pi/2),
                        ])
def test_to_frame_to_lattice(lattice,a,b,c,alpha,beta,gamma,vector,keyFrame,keyLattice):
    c = Crystal(lattice=lattice,
                a=a,b=b,c=c,
                alpha=alpha,beta=beta,gamma=gamma)
    assert np.allclose(vector,
                        c.to_frame(**{keyFrame:c.to_lattice(**{keyLattice:vector})}))

@pytest.mark.parametrize('lattice,a,b,c,alpha,beta,gamma,points',
                        [
                         ('aP',0.5,2.0,3.0,0.8,0.5,1.2,[[0,0,0]]),
                         ('mS',1.0,2.0,3.0,np.pi/2,0.5,np.pi/2,[[0,0,0],[0.5,0.5,0.0]]),
                         ('oI',0.5,1.5,3.0,np.pi/2,np.pi/2,np.pi/2,[[0,0,0],[0.5,0.5,0.5]]),
                         ('hP',1.0,1.0,1.6,np.pi/2,np.pi/2,2*np.pi/3,[[0,0,0],[2./3.,1./3.,0.5]]),
                         ('cF',1.0,1.0,1.0,np.pi/2,np.pi/2,np.pi/2,[[0,0,0],[0.0,0.5,0.5],[0.5,0.0,0.5],[0.5,0.5,0.0]]),
                        ])
def test_lattice_points(lattice,a,b,c,alpha,beta,gamma,points):
    c = Crystal(lattice=lattice,
                a=a,b=b,c=c,
                alpha=alpha,beta=beta,gamma=gamma)
    assert np.allclose(points,c.lattice_points)

@pytest.mark.parametrize('crystal,length',
                        [(Crystal(lattice='cF'),[12,6]),
                         (Crystal(lattice='cI'),[12,12,24]),
                         (Crystal(lattice='hP'),[3,3,6,12,6]),
                         (Crystal(lattice='tI',c=1.2),[2,2,2,4,2,4,2,2,4,8,4,8,8])
                        ])
def test_N_slip(crystal,length):
    assert [len(s) for s in crystal.kinematics('slip')['direction']] == length
    assert [len(s) for s in crystal.kinematics('slip')['plane']] == length

@pytest.mark.parametrize('crystal,length',
                        [(Crystal(lattice='cF'),[12]),
                         (Crystal(lattice='cI'),[12]),
                         (Crystal(lattice='hP'),[6,6,6,6]),
                        ])
def test_N_twin(crystal,length):
    assert [len(s) for s in crystal.kinematics('twin')['direction']] == length
    assert [len(s) for s in crystal.kinematics('twin')['plane']] == length

@pytest.mark.parametrize('lattice',['hP','cI','cF','tI'])
def test_Schmid(update,res_path,lattice):
    C = Crystal(lattice=lattice,c=(1.2 if lattice == 'tI' else None))                               # noqa
    for mode in ['slip']+([] if lattice == 'tI' else ['twin']):
        reference = res_path/f'{lattice}_{mode}.txt'
        P = C.Schmid(N_slip='*') if mode == 'slip' else C.Schmid(N_twin='*')
        if update:
            Table({'Schmid':(3,3,)},P.reshape(-1,9)).save(reference)
        assert np.allclose(P,Table.load(reference).get('Schmid'))

def test_Schmid_invalid():
    with pytest.raises(KeyError):
        Crystal(lattice='fcc').Schmid()

# https://doi.org/10.1016/0079-6425(94)00007-7, Fig. 22
@pytest.mark.parametrize('c_a,mode',
                        [(np.sqrt(2)*0.99,['c','c','c','c']),
                         (np.sqrt(2)*1.01,['c','c','c','t']),
                         (1.5*0.99,['c','c','c','t']),
                         (1.5*1.01,['c','c','t','t']),
                         (np.sqrt(3)*0.99,['c','c','t','t']),
                         (np.sqrt(3)*1.01,['t','c','t','t'])])
def test_characteristic_twin_direction(c_a,mode):
    C = Crystal(lattice='hP',c=c_a)
    assert (np.where(C.characteristic_shear_twin([1,1,1,1]).flatten()>0,'c','t')==mode).all()

@pytest.mark.parametrize('crystal', [Crystal(lattice='cF'),
                                     Crystal(lattice='cI'),
                                     Crystal(lattice='hP'),
                                     Crystal(lattice='tI',c=1.2)])
def test_related(crystal):
    for r in crystal.orientation_relationships:
        crystal.relation_operations(r)

@pytest.mark.parametrize('crystal', [Crystal(lattice='cF'),
                                     Crystal(lattice='cI'),
                                     Crystal(lattice='hP')])
def test_related_invalid_target(np_rng,crystal):
    relationship = np_rng.choice(crystal.orientation_relationships)
    with pytest.raises(ValueError):
        crystal.relation_operations(relationship,crystal)

@pytest.mark.parametrize('crystal', [Crystal(lattice='cF'),
                                     Crystal(lattice='cI'),
                                     Crystal(lattice='hP',c=np.sqrt(2.)*.99),
                                     Crystal(lattice='tI',c=1.2)])
@pytest.mark.parametrize('mode',['slip','twin'])
@pytest.mark.need_damask_root
def test_system_match(crystal,mode,damask_root):
    if crystal.lattice == 'tI' and mode == 'twin': return

    raw = []
    name = f'{crystal.lattice.upper()}_SYSTEM{mode.upper()}'
    with open(damask_root/'src'/'crystal.f90') as f:
        in_matrix = False
        for line in [l for l in f if l.split('!')[0].strip()]:
            if f'shape({name})' in line: break
            if in_matrix:
                entries = line.split('&')[0].strip().split(',')
                raw.append(list(filter(None, entries)))
            in_matrix |= line.strip().startswith(f'{name}') and 'reshape' in line

    d_fortran,p_fortran = np.hsplit(np.array(raw).astype(int),2)
    if crystal.lattice == 'hP':
        d_fortran = util.Bravais_to_Miller(uvtw=d_fortran)
        p_fortran = util.Bravais_to_Miller(hkil=p_fortran)

    python = crystal.kinematics(mode)

    assert np.array_equal(d_fortran,np.vstack(python['direction']))
    assert np.array_equal(p_fortran,np.vstack(python['plane']))
