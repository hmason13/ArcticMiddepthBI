#ifndef SEAICE_OPTIONS_H
#define SEAICE_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"


C By default, the sea-ice package uses its own integrated bulk
C formulae to compute fluxes over open ocean. When this flag is set,
C these variables are computed in a separate external package,
C for example, pkg/exf, and then modified for sea-ice effects by
C pkg/seaice.
#define SEAICE_EXTERNAL_FLUXES

#define SEAICE_CGRID

#define SEAICE_ALLOW_DYNAMICS
#define SEAICE_ALLOW_EVP
#define SEAICE_ALLOW_JFNK
#define SEAICE_ALLOW_KRYLOV
#define SEAICE_ZETA_SMOOTHREG

#endif /* SEAICE_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
