# hassan m, nyu
# Jul 2023
#
# This script also generates other fields through thermal wind balance.
#
# Important note for these files! When converting to bigEndian format,
# ensure that the dimensions are ordered z,y,x! (or y,x etc.)
#
# data.cal? .diagnostics? .exf!! .kpp? .pkg! .seaice!
 
import sys
import numpy as np
import scipy.interpolate as terp
import netCDF4 as nc
import matplotlib.pyplot as plt

# file i/o stuff
locLoad = '/home/hcm7920/data/'
locOther= '/home/hcm7920/MITgcm/experiments/e110-0-modified/other/'
month = 12            # jan=1, dec=12
fname = f"OCEAN_DENS_STRAT_PRESS_mon_mean_2017-{month:02}_ECCO_V4r4_latlon_0p50deg.nc"

# some constants
nx        = 250  # number of cells
ny        = 300  # ""
nz        = 50   # ""
dx        = 2000  # cell size in meters
dy        = 2000  # ""
rhoConst  = 1029. # ref density
s0        = 33    # ref salinity
g         = 9.81
tAlpha    = 0
sBeta     = 7.4e-4
f0        = 1.4e-4
halfDeg   = 222000    # 0.5deg latitude in meters 
pdelta    = 4.4e-10   # pressure parameter for linear equation of state

# indices for profile location
dVI = -8              # deepValueIndex
myLat = 324           # gives 71.25N, +1 index gives +0.5 deg
myLon = 59            # gives 150.25W +1 index gives -0.5 deg

# load ecco data
data = nc.Dataset(locLoad+fname)
lat = data['latitude']
lon = data['longitude']
depth = data['Z'][:dVI]
myRho1 = data['RHOAnoma'][0,:dVI,myLat,myLon]
myRho2 = data['RHOAnoma'][0,:dVI,myLat+1,myLon]
vertiDG1 = data['DRHODR'][0,:dVI,myLat,myLon]
vertiDG2 = data['DRHODR'][0,:dVI,myLat+1,myLon]

# some necessary depth fields
depthDeltas = (depth[:-1]-depth[1:])
depthBetween = (depth[1:]+depth[:-1])/2

# calculate desired profiles
horizDGData     = (myRho2-myRho1) / halfDeg # @ 150.25W 71.5N
vertiDGData     = (vertiDG1+vertiDG2) / 2   # same location, hopefully prudent??
nSquare     = vertiDGData * g / rhoConst
vertiShear  = horizDGData * g / rhoConst / f0

horizDGData.fill_value=0

# make spline for grid conversion
horizDGSpline = terp.splrep(-depth,horizDGData)
vertiDGSpline = terp.splrep(-depth,vertiDGData)
rho1Spline    = terp.splrep(-depth,myRho1)

# save as bin vertical cell thickness
# 10 in surface layer, increasing to 300 in abyss,
# same spacing as ECCO 0.5deg with extra 10m cells near surface
padSize = 50 - depthDeltas.shape[0] + 1
dz = np.pad(depthDeltas,(padSize,0),constant_values=10.)[:-1]
dz.astype('>f4').tofile('eccoZCoordSpacingPadded.bin')
np.save("eccoZCoordSpacingPadded", dz)

# determine new vertical grid
zCenters      = np.zeros((nz,))                   # depth of cell centers
zFaces        = np.zeros((nz,))                   # depth at cell faces
for i in range(50):
  zCenters[i] = dz[i]/2 + dz[:i].sum()
  zFaces[i] = dz[:i+1].sum()
zDeltas       = (zCenters[1:]-zCenters[:-1])      # distance between centers
zHalfBetween  = (zCenters[1:]+zCenters[:-1]) / 2  # depth half between centers

# determine horizontal grid
xdist = np.zeros((nx,)) # get horizontal position
for i in range(nx): xdist[i] = dx*i + (dx/2) # from m to km
xdist /= 1000

# use spline to interpolate fields onto new vertical grid
horizDGHalfs  = terp.splev(zHalfBetween,horizDGSpline)
horizDG       = terp.splev(zCenters,horizDGSpline)
vertiDGHalfs  = terp.splev(zHalfBetween,vertiDGSpline)
vertiDG       = terp.splev(zCenters,vertiDGSpline)
rhoProfile    = terp.splev(zCenters,rho1Spline)
nSquare       = vertiDG * g / rhoConst
vertiShear    = horizDGHalfs * g / rhoConst / f0

# save some plots of the data for viewing later
plt.figure(figsize=(10,10))
plt.plot(horizDG,-zCenters)
plt.title(f"Horizontal density gradient v depth, month {month}", fontsize=18)
plt.ylabel("Depth (m)", fontsize=14)
plt.xlabel("(kg m-4)", fontsize=14)
plt.savefig(locOther+f"horizDGProfile-month{month}")

plt.figure(figsize=(10,10))
plt.plot(nSquare,-zCenters)
plt.title(f"N Square v depth, month {month}", fontsize=18)
plt.ylabel("Depth (m)", fontsize=14)
plt.xlabel("(s -2) ", fontsize=14)
plt.savefig(locOther+f"nSquareProfile-month{month}")

plt.figure(figsize=(10,10))
plt.plot(vertiShear,-zHalfBetween)
plt.title(f"Mean shear v depth, month {month}", fontsize=18)
plt.ylabel("Depth (m)", fontsize=14)
plt.xlabel("(s -1) ", fontsize=14)
plt.savefig(locOther+f"meanShearProfile-month{month}")


# generate bathymetry field
bathy = -1*np.ones((ny,nx))*dz.sum()
bathy[:,-1] = 0                     # wall
bathy.astype('>f4').tofile('bathy-2km.bin')

# save u: 0 m/s
u = np.zeros((nz,ny,nx))
u.astype('>f4').tofile('u-2km.bin')

# generate v from thermal wind balance: 0m/s at domain floor
velocities = np.zeros(50)
for i in range(49,0,-1):
  velocities[i-1] = velocities[i] + vertiShear[i-1]*zDeltas[i-1]
v = np.ones((nz,ny,nx))*velocities[:,np.newaxis,np.newaxis]
v.astype('>f4').tofile('v-2km.bin')

# generate initial density field: RHOAnoma
rho = np.ones((nz,ny,nx))*rhoProfile[:,np.newaxis,np.newaxis]
for ix in range(1,nx):
  rho[:,:,ix] = rho[:,:,ix-1] + (horizDG*2000)[:,np.newaxis]
rho.astype('>f4').tofile('rhoAnoma-2km.bin')
rhoTotal = rho+rhoConst

# generate initial SSH
eta = (-1*(rhoTotal)*dz[:,np.newaxis,np.newaxis]).sum(0) / rhoTotal[0,:,:]
eta -= eta.mean((0,1))      # make sure average eta is 0
eta.astype('>f4').tofile('eta-2km.bin')

# generate initial temperature field
# tempStepping = FALSE ==> temp identically 1
temp = np.ones((nz,ny,nx))
temp.astype('>f4').tofile('temp-2km.bin')

# generate initial hydrostatic pressure field
pressure = np.zeros((nz,ny,nx))
for i in range(49):
  pressure[i+1,:,:] = g * (rhoConst) * zCenters[i]

# generate initial salinity field from rho field
# note: linear eq of state 
# note: must consider the pressure contribution
# note: add a perturbation in this field
salt = (rho+rhoConst) / rhoConst
salt = (salt-1-(pdelta*pressure)) / sBeta
salt = salt + s0
salt.astype('>f4').tofile('saltUnperturbed-2km.bin')
salt += np.fromfunction(lambda z,y,x: 0.001*np.sin(y/10)*np.exp((-(125-x)**2)/1000),
			 (nz,ny,nx))
salt.astype('>f4').tofile('salt-2km.bin')

# generate restoring mask
mask = np.ones((nz,ny,nx))
mask.astype('>f4').tofile('restore.bin')

