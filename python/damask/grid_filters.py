import numpy as np

def __ks(size,field):
    """Get differential operator."""
    grid = np.array(np.shape(field)[0:3])

    k_sk = np.where(np.arange(grid[0])>grid[0]//2,np.arange(grid[0])-grid[0],np.arange(grid[0]))/size[0]
    if grid[0]%2 == 0: k_sk[grid[0]//2] = 0                                                         # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

    k_sj = np.where(np.arange(grid[1])>grid[1]//2,np.arange(grid[1])-grid[1],np.arange(grid[1]))/size[1]
    if grid[1]%2 == 0: k_sj[grid[1]//2] = 0                                                         # Nyquist freq=0 for even grid (Johnson, MIT, 2011)

    k_si = np.arange(grid[2]//2+1)/size[2]

    kk, kj, ki = np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij')
    return np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3)


def curl(size,field):
    """Calculate curl of a vector or tensor field in Fourier space."""
    n = np.prod(field.shape[3:])
    k_s = __ks(size,field)

    e = np.zeros((3, 3, 3))
    e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = +1.0                                                     # Levi-Civita symbol 
    e[0, 2, 1] = e[2, 1, 0] = e[1, 0, 2] = -1.0

    field_fourier = np.fft.rfftn(field,axes=(0,1,2))
    curl = (np.einsum('slm,ijkl,ijkm ->ijks', e,k_s,field_fourier)*2.0j*np.pi if n == 3 else        # vector, 3   -> 3
            np.einsum('slm,ijkl,ijknm->ijksn',e,k_s,field_fourier)*2.0j*np.pi)                      # tensor, 3x3 -> 3x3

    return np.fft.irfftn(curl,axes=(0,1,2),s=field.shape[0:3])


def divergence(size,field):
    """Calculate divergence of a vector or tensor field in Fourier space."""
    n = np.prod(field.shape[3:])
    k_s = __ks(size,field)

    field_fourier = np.fft.rfftn(field,axes=(0,1,2))
    divergence = (np.einsum('ijkl,ijkl ->ijk', k_s,field_fourier)*2.0j*np.pi if n == 3 else         # vector, 3   -> 1
                  np.einsum('ijkm,ijklm->ijkl',k_s,field_fourier)*2.0j*np.pi)                       # tensor, 3x3 -> 3

    return np.fft.irfftn(divergence,axes=(0,1,2),s=field.shape[0:3])


def gradient(size,field):
    """Calculate gradient of a vector or scalar field in Fourier space."""
    n = np.prod(field.shape[3:])
    k_s = __ks(size,field)

    field_fourier = np.fft.rfftn(field,axes=(0,1,2))
    gradient = (np.einsum('ijkl,ijkm->ijkm', field_fourier,k_s)*2.0j*np.pi if n == 1 else           # scalar, 1 -> 3
                np.einsum('ijkl,ijkm->ijklm',field_fourier,k_s)*2.0j*np.pi)                         # vector, 3 -> 3x3

    return np.fft.irfftn(gradient,axes=(0,1,2),s=field.shape[0:3])


#--------------------------------------------------------------------------------------------------
def displacementFluctFFT(F,size):
  """Calculate displacement field from deformation gradient field."""
  integrator = 0.5j * size / np.pi

  kk, kj, ki = np.meshgrid(np.where(np.arange(grid[2])>grid[2]//2,np.arange(grid[2])-grid[2],np.arange(grid[2])),
                           np.where(np.arange(grid[1])>grid[1]//2,np.arange(grid[1])-grid[1],np.arange(grid[1])),
                                    np.arange(grid[0]//2+1),
                           indexing = 'ij')
  k_s = np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3) 
  k_sSquared = np.einsum('...l,...l',k_s,k_s)
  k_sSquared[0,0,0] = 1.0                                                                           # ignore global average frequency

#--------------------------------------------------------------------------------------------------
# integration in Fourier space

  displacement_fourier = -np.einsum('ijkml,ijkl,l->ijkm',
                                    np.fft.rfftn(F,axes=(0,1,2)),
                                    k_s,
                                    integrator,
                                   ) / k_sSquared[...,np.newaxis]

#--------------------------------------------------------------------------------------------------
# backtransformation to real space

  return np.fft.irfftn(displacement_fourier,grid[::-1],axes=(0,1,2))
