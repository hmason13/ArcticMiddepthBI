# arcticMiddepthBI
Welcome! Please note this experiment is still under construction as of 1 Sep 2023.

## Experiment description
This experiment represents a flat-bottom channel initialized with a balanced & baroclinically unstable
density field. This initialization is generated from ECCOv4 data in the Beaufort gyre region.
Simple linear stability analysis on the raw ECCO data implies a Phillips like BI regime, one that
originates in the fluid interior. This setup captures some aspects of an eddy field generated through
this mechanism.

The experiment is forced and dissipated. Dissipation is standard (Leith), while the forcing is
a restoring everwhere determined by the difference between channel mean (2d) fields and initial fields.

## Numerical details
This experiment includes a few modified Fortran files. These files change the behavior of the model.
Specifically, the functionality of the RBCS package is changed to restore based on
meridional mean fields instead of the actual model fields. Note that this implies the RBCS package
**does not** function as described in the standard model documentation. A more complete
description of the behavior will be provided at some point.
