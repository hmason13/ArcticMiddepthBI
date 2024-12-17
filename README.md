## Experiment description
This experiment represents a flat-bottom channel initialized with a balanced & baroclinically unstable
density field. This initialization is generated from ECCOv4 data in the Beaufort gyre region.
Simple linear stability analysis on the raw ECCO data implies a Phillips like BI regime, one that
originates in the fluid interior. This setup is designed to explore baroclinic instability under
mechanically active sea ice. The sea ice thermodynamics here are turned off; the ice neither melts nor
grows.

The experiment is forced and dissipated. Dissipation is set by the QGLeith scheme, scaling with grid
enstrophy. The forcing is a zonal restoring everwhere determined by the difference between channel mean
fields and initial fields. This forcing scheme is not basic MITgcm functionality and was constructed
for this experiment.

## Numerical details
This experiment includes a few modified Fortran files. These files change the behavior of the model.
Specifically, the functionality of the RBCS package is changed to restore based on
the along channel mean fields instead of the actual model fields. This change means the RBCS package
**does not** function as described in the standard model documentation. The model calculates zonal
averages at every timestep. The zonal averages are then used in place of local values when
calculating the restoring terms. Diagnostics compatible with the RBCS package also work
with this modified version. The data.rbcs file still defines the target fields, restoring timescale, and
toggles the forcing for each field.

## Other details
The modified forcing is mpi friendly for my particular build choices.\
Specifically:\
#define GLOBAL_SUM_ORDER_TILES\
#undef GLOBAL_SUM_SEND_RECV\
#define ALLOW_USE_MPI\
nSx = 1\
nSy = 1\
This setup is multi process but single threaded!\

The link to the ECCO files used in initialization is below. All twelve files are necessary here
because an average is taken over the entire year.\
https://cmr.earthdata.nasa.gov/virtual-directory/collections/C1990404798-POCLOUD/temporal/2017\

Files to analyze the model output and produce the submitted figures are also in the respository. Some
notebooks present were not required to perform the final analyses.

MITgcm checkpoint67y
