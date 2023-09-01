# hassan m, nyu
# jan 2023
#
# This script generates the initial Cice file. This file varies with parameters
# and so takes in an argument through the shell.
#
# Important note for these files! When converting to bigEndian format,
# -- ensure that the dimensions are ordered z,y,x! (or y,x etc.)

import sys
import numpy as np

nx = 250  # number of cells
ny = 300  # ""
nz = 50   # ""
dx = 2000  # cell size in meters
dy = 2000  # ""

# The argument should be a percentage from 0 to 100 representing the
# -- concentration of sea ice in an individual cell.
# conc = float(sys.argv[1])
# conc = np.round(conc, 0)
conc = 0

# Cice area
area = 0.01 * conc * np.ones((ny,nx))
area.astype('>f4').tofile(f'iceAreaLowRes.bin')

# Cice heff (effective thickness)
area = 2 * np.ones((ny,nx))
area.astype('>f4').tofile(f'iceHeffLowRes.bin')
