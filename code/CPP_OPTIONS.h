C $Header: /u/gcmpack/MITgcm/verification/exp4/code/CPP_OPTIONS.h,v 1.16 2008/01/23 16:42:47 jmc Exp $
C $Name:  $

#ifndef CPP_OPTIONS_H
#define CPP_OPTIONS_H

C CPP flags controlling particular source code features

C o Shortwave heating as extra term in external_forcing.F
C Note: this should be a run-time option
#undef SHORTWAVE_HEATING

C adding this option so that lwdown works?
C hassan
#undef ALLOW_DOWNWARD_RADIATION

#define DISABLE_SEAICE_GROWTH

C o Include/exclude phi_hyd calculation code
#define INCLUDE_PHIHYD_CALCULATION_CODE

C o Include/exclude call to S/R CONVECT
#define INCLUDE_CONVECT_CALL

C o Include/exclude call to S/R CALC_DIFFUSIVITY
#define INCLUDE_CALC_DIFFUSIVITY_CALL

C o Allow latitudinally varying BryanLewis79 vertical diffusivity
#undef ALLOW_BL79_LAT_VARY

C o Allow full 3D specification of vertical diffusivity
#define ALLOW_3D_DIFFKR

C o Include/exclude Implicit vertical advection code
#define INCLUDE_IMPLVERTADV_CODE

C o Include/exclude AdamsBashforth-3rd-Order code
#undef ALLOW_ADAMSBASHFORTH_3

C o Include/exclude nonHydrostatic code
#define ALLOW_NONHYDROSTATIC

C o Include pressure loading code
#define ATMOSPHERIC_LOADING

C o exclude/allow external forcing-fields load 
C   this allows to read & do simple linear time interpolation of oceanic
C   forcing fields, if no specific pkg (e.g., EXF) is used to compute them.
#undef EXCLUDE_FFIELDS_LOAD

C o Use "Exact Convervation" of fluid in Free-Surface formulation
C   so that d/dt(eta) is exactly equal to - Div.Transport
#define EXACT_CONSERV

C o Allow the use of Non-Linear Free-Surface formulation
C   this implies that surface thickness (hFactors) vary with time
#undef NONLIN_FRSURF

C o ALLOW isotropic scaling of harmonic and bi-harmonic terms when
C   using an locally isotropic spherical grid with (dlambda) x (dphi*cos(phi))
C *only for use on a lat-lon grid*
C   Setting this flag here affects both momentum and tracer equation unless
C   it is set/unset again in other header fields (e.g., GAD_OPTIONS.h).
C   The definition of the flag is commented to avoid interference with
C   such other header files.
C   The preferred method is specifying a value for viscAhGrid or viscA4Grid
C   in data which is then automatically scaled by the grid size;
C   the old method of specifying viscAh/viscA4 and this flag is provided
C   for completeness only (and for use with the adjoint).
C#define ISOTROPIC_COS_SCALING

C o This flag selects the form of COSINE(lat) scaling of bi-harmonic term.
C *only for use on a lat-lon grid*
C   Has no effect if ISOTROPIC_COS_SCALING is undefined.
C   Has no effect on vector invariant momentum equations.
C   Setting this flag here affects both momentum and tracer equation unless
C   it is set/unset again in other header fields (e.g., GAD_OPTIONS.h).
C   The definition of the flag is commented to avoid interference with
C   such other header files.
C#define COSINEMETH_III

C o Use "OLD" UV discretisation near boundaries (*not* recommended)
C   Note - only works with  #undef NO_SLIP_LATERAL  in calc_mom_rhs.F
C          because the old code did not have no-slip BCs
#undef  OLD_ADV_BCS

C o Use LONG.bin, LATG.bin, etc., initialization for ini_curviliear_grid.F
C   Default is to use "new" grid files (OLD_GRID_IO undef) but OLD_GRID_IO
C   is still useful with, e.g., single-domain curvilinear configurations.
#undef OLD_GRID_IO

C o Execution environment support options
#include "CPP_EEOPTIONS.h"

C o Include/exclude code specific to the ECCO/SEALION version.
C   AUTODIFF or EXF package.
C   Currently controled by a single header file
C   For this to work, PACKAGES_CONFIG.h needs to be included!
cph#if (defined (ALLOW_AUTODIFF) || \
cph     defined (ALLOW_ECCO) || \
cph     defined (ALLOW_EXF))
cph# include "ECCO_CPPOPTIONS.h"
cph#endif

#endif /* CPP_OPTIONS_H */

