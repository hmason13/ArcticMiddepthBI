import numpy as np
import matplotlib.pyplot as plt
import pyqg
from pyqg import diagnostic_tools as tools


    ### Grid ###
nx = 64
nz = 50
L = 30e3
H = np.load(‘/scratch/projects/shaferlab/hassan/dz.npy’)
H = np.array(H.tolist()) # removes the dtype associated with the above dz.npy


    ### Planetary parameters ###
f0 = 1.4e-4
beta = 0


    ### Background shear and density profiles ###
U = np.load(‘/scratch/projects/shaferlab/hassan/velocities.npy’)
V = 0 * U
rho = np.load(‘/scratch/projects/shaferlab/hassan/rho.npy’)


    ### Time for pyqg’s sake ###
dt = 3600
tmax = dt
m = pyqg.LayeredModel(nx = nx, nz = nz, L = L, H = H, U = U, V = V,
                      rho = rho, f = f0, beta = beta, dt = dt, tmax = tmax)


# Do the stability analysis
evals, evecs = m.stability_analysis()


# Etc. etc.
evals = np.fft.fftshift(evals.imag,axes=(0,))
k, l = m.k * m.radii[1], np.fft.fftshift(m.l, axes=(0,)) * m.radii[1]
argmax = evals[int(m.ny / 2),:].argmax()
evec = np.fft.fftshift(evecs,axes=(1))[:, int(m.ny / 2), argmax]
kmax = k[int(m.ny / 2), argmax]
mag, phase = np.abs(evec), np.arctan2(evec.imag, evec.real)


# Plots
plt.figure(figsize=(8,5))
plt.contourf(k, l, evals, levels = 100)
plt.colorbar()
plt.xlabel(r’$k \, L_d$‘); plt.ylabel(r’$l \, L_d$‘)
plt.title(f’fastest growing mode’)


# Vertical structure associated with most unstable mode
zc = -np.cumsum(H)
plt.plot(mag, zc)
plt.grid()
plt.xlabel(f’eigenvector amplitude $|\psi|$’)
plt.ylabel(‘depth (m)’)
plt.title(‘vertical structure associated with most unstable mode’)