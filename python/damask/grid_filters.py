from scipy import spatial
import numpy as np

def __ks(size,field,first_order=False):
    """Get wave numbers operator."""
    grid = np.array(np.shape(field)[:3])

    k_sk = np.where(np.arange(grid[0])>grid[0]//2,np.arange(grid[0])-grid[0],np.arange(grid[0]))/size[0]
    if grid[0]%2 == 0 and first_order: k_sk[grid[0]//2] = 0                                         # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

    k_sj = np.where(np.arange(grid[1])>grid[1]//2,np.arange(grid[1])-grid[1],np.arange(grid[1]))/size[1]
    if grid[1]%2 == 0 and first_order: k_sj[grid[1]//2] = 0                                         # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

    k_si = np.arange(grid[2]//2+1)/size[2]

    kk, kj, ki = np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij')
    return np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3)


def curl(size,field):
    """Calculate curl of a vector or tensor field in Fourier space."""
    n = np.prod(field.shape[3:])
    k_s = __ks(size,field,True)

    e = np.zeros((3, 3, 3))
    e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = +1.0                                                     # Levi-Civita symbol 
    e[0, 2, 1] = e[2, 1, 0] = e[1, 0, 2] = -1.0

    field_fourier = np.fft.rfftn(field,axes=(0,1,2))
    curl = (np.einsum('slm,ijkl,ijkm ->ijks', e,k_s,field_fourier)*2.0j*np.pi if n == 3 else        # vector, 3   -> 3
            np.einsum('slm,ijkl,ijknm->ijksn',e,k_s,field_fourier)*2.0j*np.pi)                      # tensor, 3x3 -> 3x3

    return np.fft.irfftn(curl,axes=(0,1,2),s=field.shape[:3])


def divergence(size,field):
    """Calculate divergence of a vector or tensor field in Fourier space."""
    n = np.prod(field.shape[3:])
    k_s = __ks(size,field,True)

    field_fourier = np.fft.rfftn(field,axes=(0,1,2))
    divergence = (np.einsum('ijkl,ijkl ->ijk', k_s,field_fourier)*2.0j*np.pi if n == 3 else         # vector, 3   -> 1
                  np.einsum('ijkm,ijklm->ijkl',k_s,field_fourier)*2.0j*np.pi)                       # tensor, 3x3 -> 3

    return np.fft.irfftn(divergence,axes=(0,1,2),s=field.shape[:3])


def gradient(size,field):
    """Calculate gradient of a vector or scalar field in Fourier space."""
    n = np.prod(field.shape[3:])
    k_s = __ks(size,field,True)

    field_fourier = np.fft.rfftn(field,axes=(0,1,2))
    gradient = (np.einsum('ijkl,ijkm->ijkm', field_fourier,k_s)*2.0j*np.pi if n == 1 else           # scalar, 1 -> 3
                np.einsum('ijkl,ijkm->ijklm',field_fourier,k_s)*2.0j*np.pi)                         # vector, 3 -> 3x3

    return np.fft.irfftn(gradient,axes=(0,1,2),s=field.shape[:3])


def cell_coord0(grid,size):
    """Cell center positions (undeformed)."""
    delta = size/grid*0.5
    x, y, z = np.meshgrid(np.linspace(delta[2],size[2]-delta[2],grid[2]),
                          np.linspace(delta[1],size[1]-delta[1],grid[1]),
                          np.linspace(delta[0],size[0]-delta[0],grid[0]),
                          indexing = 'ij')

    return np.concatenate((z[:,:,:,None],y[:,:,:,None],x[:,:,:,None]),axis = 3) 

def cell_displacement_fluct(size,F):
    """Cell center displacement field from fluctuation part of the deformation gradient field."""
    integrator = 0.5j*size/np.pi

    k_s = __ks(size,F,False)
    k_s_squared = np.einsum('...l,...l',k_s,k_s)
    k_s_squared[0,0,0] = 1.0

    displacement = -np.einsum('ijkml,ijkl,l->ijkm',
                              np.fft.rfftn(F,axes=(0,1,2)),
                              k_s,
                              integrator,
                             ) / k_s_squared[...,np.newaxis]

    return np.fft.irfftn(displacement,axes=(0,1,2),s=F.shape[:3])

def cell_displacement_avg(size,F):
    """Cell center displacement field from average part of the deformation gradient field."""
    F_avg = np.average(F,axis=(0,1,2))
    return np.einsum('ml,ijkl->ijkm',F_avg-np.eye(3),cell_coord0(F.shape[:3],size))


def node_coord0(grid,size):
    """Nodal positions (undeformed)."""
    x, y, z = np.meshgrid(np.linspace(0,size[2],1+grid[2]),
                          np.linspace(0,size[1],1+grid[1]),
                          np.linspace(0,size[0],1+grid[0]),
                          indexing = 'ij')
    
    return np.concatenate((z[:,:,:,None],y[:,:,:,None],x[:,:,:,None]),axis = 3) 

def node_displacement_fluct(size,F):
    """Nodal displacement field from fluctuation part of the deformation gradient field."""
    return cell_2_node(cell_displacement_fluct(size,F))

def node_displacement_avg(size,F):
    """Nodal displacement field from average part of the deformation gradient field."""
    F_avg = np.average(F,axis=(0,1,2))
    return np.einsum('ml,ijkl->ijkm',F_avg-np.eye(3),node_coord0(F.shape[:3],size))


def cell_2_node(cell_data):
    """Interpolate cell data to nodal data."""
    n = (  cell_data + np.roll(cell_data,1,(0,1,2))
         + np.roll(cell_data,1,(0,))  + np.roll(cell_data,1,(1,))  + np.roll(cell_data,1,(2,))
         + np.roll(cell_data,1,(0,1)) + np.roll(cell_data,1,(1,2)) + np.roll(cell_data,1,(2,0)))*0.125
  
    return np.pad(n,((0,1),(0,1),(0,1))+((0,0),)*len(cell_data.shape[3:]),mode='wrap')

def node_2_cell(node_data):
    """Interpolate nodal data to cell data."""
    c = (  node_data + np.roll(node_data,1,(0,1,2))
         + np.roll(node_data,1,(0,))  + np.roll(node_data,1,(1,))  + np.roll(node_data,1,(2,))
         + np.roll(node_data,1,(0,1)) + np.roll(node_data,1,(1,2)) + np.roll(node_data,1,(2,0)))*0.125
  
    return c[:-1,:-1,:-1]

def regrid(size,F,new_grid):
    """tbd."""
    c = cell_coord0(F.shape[:3][::-1],size) \
      + cell_displacement_avg(size,F) \
      + cell_displacement_fluct(size,F)

    outer = np.dot(np.average(F,axis=(0,1,2)),size)
    for d in range(3):
        c[np.where(c[:,:,:,d]<0)]        += outer[d]
        c[np.where(c[:,:,:,d]>outer[d])] -= outer[d]

    tree = spatial.cKDTree(c.reshape((-1,3)),boxsize=outer)
    return tree.query(cell_coord0(new_grid,outer))[1]
