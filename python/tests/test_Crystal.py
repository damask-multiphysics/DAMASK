import pytest
import numpy as np
from pathlib import Path

import damask
from damask import Crystal
from damask import util

class TestCrystal:

    @pytest.mark.parametrize('lattice,family',[('aP','cubic'),('xI','cubic')])
    def test_invalid_init(self,lattice,family):
        with pytest.raises(KeyError):
            Crystal(family=family,lattice=lattice)

    def test_eq(self):
        family = np.random.choice(list(damask._crystal.lattice_symmetries.values()))
        assert Crystal(family=family) == Crystal(family=family)

    def test_double_to_lattice(self):
        c = Crystal(lattice='cF')
        with pytest.raises(KeyError):
            c.to_lattice(direction=np.ones(3),plane=np.ones(3))

    def test_double_to_frame(self):
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
    def test_bases_contraction(self,lattice,a,b,c,alpha,beta,gamma):
        c = Crystal(lattice=lattice,
                    a=a,b=b,c=c,
                    alpha=alpha,beta=beta,gamma=gamma)
        assert np.allclose(np.eye(3),np.einsum('ik,jk',c.basis_real,c.basis_reciprocal))

    def test_basis_invalid(self):
        with pytest.raises(KeyError):
            Crystal(family='cubic').basis_real

    def test_basis_real(self):
        for gamma in np.random.random(2**8)*np.pi:
            basis = np.tril(np.random.random((3,3))+1e-6)
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
    def test_to_frame_to_lattice(self,lattice,a,b,c,alpha,beta,gamma,vector,keyFrame,keyLattice):
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
    def test_lattice_points(self,lattice,a,b,c,alpha,beta,gamma,points):
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
    def test_N_slip(self,crystal,length):
        assert [len(s) for s in crystal.kinematics('slip')['direction']] == length
        assert [len(s) for s in crystal.kinematics('slip')['plane']] == length

    @pytest.mark.parametrize('crystal,length',
                            [(Crystal(lattice='cF'),[12]),
                             (Crystal(lattice='cI'),[12]),
                             (Crystal(lattice='hP'),[6,6,6,6]),
                            ])
    def test_N_twin(self,crystal,length):
        assert [len(s) for s in crystal.kinematics('twin')['direction']] == length
        assert [len(s) for s in crystal.kinematics('twin')['plane']] == length

    @pytest.mark.parametrize('c_a,mode',
                            [(np.sqrt(2)*0.99,['c','c','c','c']),
                             (np.sqrt(2)*1.01,['c','c','c','t']),
                             (1.5*0.99,['c','c','c','t']),
                             (1.5*1.01,['c','c','t','t']),
                             (np.sqrt(3)*0.99,['c','c','t','t']),
                             (np.sqrt(3)*1.01,['t','c','t','t'])])
    def test_characteristic_twin_direction(self,c_a,mode):
        C = Crystal(lattice='hP',c=c_a)
        assert (np.where(C.characteristic_shear_twin()>0,'c','t')==mode).all()

    @pytest.mark.parametrize('crystal', [Crystal(lattice='cF'),
                                         Crystal(lattice='cI'),
                                         Crystal(lattice='hP'),
                                         Crystal(lattice='tI',c=1.2)])
    def test_related(self,crystal):
        for r in crystal.orientation_relationships:
            crystal.relation_operations(r)

    @pytest.mark.parametrize('crystal', [Crystal(lattice='cF'),
                                         Crystal(lattice='cI'),
                                         Crystal(lattice='hP')])
    def test_related_invalid_target(self,crystal):
        relationship = np.random.choice(crystal.orientation_relationships)
        with pytest.raises(ValueError):
            crystal.relation_operations(relationship,crystal)

    @pytest.mark.parametrize('crystal', [Crystal(lattice='cF'),
                                         Crystal(lattice='cI'),
                                         Crystal(lattice='hP',c=np.sqrt(2.)*.99),
                                         Crystal(lattice='tI',c=1.2)])
    @pytest.mark.parametrize('mode',['slip','twin'])
    @pytest.mark.need_damaskroot
    def test_system_match(self,crystal,mode,damaskroot):
        if crystal.lattice == 'tI' and mode == 'twin': return

        raw = []
        name = f'{crystal.lattice.upper()}_SYSTEM{mode.upper()}'
        with open(Path(damaskroot).expanduser()/'src'/'crystal.f90') as f:
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
