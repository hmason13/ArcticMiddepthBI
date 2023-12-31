C $Header: /u/gcmpack/MITgcm/model/src/external_forcing_surf.F,v 1.63 2014/05/21 21:48:32 heimbach Exp $
C $Name:  $
C hassan: adapted from jc restoring code

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif
#ifdef ALLOW_SALT_PLUME
# include "SALT_PLUME_OPTIONS.h"
#endif
#undef CHECK_OVERLAP_FORCING

CBOP
C     !ROUTINE: EXTERNAL_FORCING_SURF
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_SURF(
     I             iMin, iMax, jMin, jMax,
     I             myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE EXTERNAL_FORCING_SURF
C     | o Determines forcing terms based on external fields
C     |   relaxation terms etc.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "FFIELDS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "SURFACE.h"
#ifdef ALLOW_AUTODIFF
# include "tamc.h"
# include "tamc_keys.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     iMin,iMax, jMin,jMax :: Range of points for calculation
C     myTime :: Current time in simulation
C     myIter :: Current iteration number in simulation
C     myThid :: Thread no. that called this routine.
      INTEGER iMin, iMax
      INTEGER jMin, jMax
      _RL myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     === Local variables ===
C     bi,bj  :: tile indices
C     i,j    :: loop indices
C     ks     :: index of surface interface layer
C     hassan: helper variables to compute average
      INTEGER bi,bj
      INTEGER i,j,k
      INTEGER ks
      INTEGER valCol ! hassan: index for this thread's column of tiles
      _RL recip_Cp
      _RL sum_T(lshare8), sum_S(lshare8) 
      _RL sum_U(lshare8), sum_V(lshare8) ! hassan: helper variables
      _RL meridArea
      _RL T_tile(nSx,nSy), S_tile(nSx,nSy) ! hassan: more helper arrays
      _RL U_tile(nSx,nSy), V_tile(nSx,nSy) ! hassan
#ifdef ALLOW_BALANCE_FLUXES
      _RS tmpVar(1)
#endif
#ifdef CHECK_OVERLAP_FORCING
      _RS fixVal
#endif
CEOP

      IF ( usingPCoords ) THEN
       ks        = Nr
      ELSE
       ks        = 1
      ENDIF
      recip_Cp = 1. _d 0 / HeatCapacity_Cp

C hassan: calculate mean_T, mean_S, mean_U, mean_V (meridional average of theta,
C s, u, and v) GLOBAL_SUM_TILE_MERID returns nPx*nSx values, and each thread
C only needs to know one of the values. This MPI function must be called once
C for each horizontal grid cell in a tile. Inefficient, but maybe not too bad?

      DO k=1,Nr
        sum_T = 0.
        sum_S = 0.
        sum_U = 0.
        sum_V = 0.
        DO bj=myByLo(myThid), myByHi(myThid)
          DO bi=myBxLo(myThid), myBxHi(myThid)
            DO i=iMin,iMax
              T_tile(bi,bj) = 0.
              S_tile(bi,bj) = 0.
              U_tile(bi,bj) = 0.
              V_tile(bi,bj) = 0.
              meridArea = 0.
              DO j=jMin,jMax
C               hassan: modified to mean over only meridional slices
                T_tile(bi,bj) = T_tile(bi,bj)
     &            + rA(i,j,bi,bj)*theta(i,j,k,bi,bj)
                S_tile(bi,bj) = S_tile(bi,bj)
     &            + rA(i,j,bi,bj)*salt(i,j,k,bi,bj)
C               JC: The following is probably not right for uneven grids.
                U_tile(bi,bj) = U_tile(bi,bj)
     &            + rA(i,j,bi,bj)*uVel(i,j,k,bi,bj)
                V_tile(bi,bj) = V_tile(bi,bj)
     &            + rA(i,j,bi,bj)*vVel(i,j,k,bi,bj)
                meridArea = meridArea + rA(i,j,bi,bj)
              ENDDO
              T_tile(bi,bj) = T_tile(bi,bj) / meridArea
              S_tile(bi,bj) = S_tile(bi,bj) / meridArea
              U_tile(bi,bj) = U_tile(bi,bj) / meridArea
              V_tile(bi,bj) = V_tile(bi,bj) / meridArea
              CALL GLOBAL_SUM_TILE_MERID(T_tile, sum_T, myThid)
              CALL GLOBAL_SUM_TILE_MERID(S_tile, sum_S, myThid)
              CALL GLOBAL_SUM_TILE_MERID(U_tile, sum_U, myThid)
              CALL GLOBAL_SUM_TILE_MERID(V_tile, sum_V, myThid)
              valCol = (myXGlobalLo-1)/sNx + bi
              mean_T(i,k,bi,bj) = sum_T(valCol) / (nPy*nSy) ! divided by # of tiles
              mean_S(i,k,bi,bj) = sum_S(valCol) / (nPy*nSy) ! divided by # of tiles
              mean_U(i,k,bi,bj) = sum_U(valCol) / (nPy*nSy) ! divided by # of tiles
              mean_V(i,k,bi,bj) = sum_V(valCol) / (nPy*nSy) ! divided by # of tiles
              _BARRIER
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Apply adjustment (balancing forcing) and exchanges
C     to oceanic surface forcing

#ifdef ALLOW_BALANCE_FLUXES
C     balance fluxes
      tmpVar(1) = oneRS
      IF ( balanceEmPmR .AND. (.NOT.useSeaice .OR. useThSIce) ) THEN
        CALL REMOVE_MEAN_RS( 1, EmPmR, maskInC, maskInC, rA, tmpVar,
     &                       'EmPmR', myTime, myThid )
      ENDIF
      IF ( balanceQnet  .AND. (.NOT.useSeaice .OR. useThSIce) ) THEN
        CALL REMOVE_MEAN_RS( 1, Qnet,  maskInC, maskInC, rA, tmpVar,
     &                       'Qnet ', myTime, myThid )
      ENDIF
#endif /* ALLOW_BALANCE_FLUXES */

C-    Apply exchanges (if needed)

#ifdef CHECK_OVERLAP_FORCING
C     Put large value in overlap of forcing array to check if exch is needed
c     IF ( .NOT. useKPP ) THEN
       fixVal = 1.
       CALL RESET_HALO_RS ( EmPmR, fixVal, 1, myThid )
       fixVal = 400.
       CALL RESET_HALO_RS ( Qnet, fixVal, 1, myThid )
       fixVal = -200.
       CALL RESET_HALO_RS ( Qsw, fixVal, 1, myThid )
       fixVal = 40.
       CALL RESET_HALO_RS ( saltFlux, fixVal, 1, myThid )
c     ENDIF
#endif /* CHECK_OVERLAP_FORCING */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef EXACT_CONSERV
C NB: synchronous time step: PmEpR lag 1 time step behind EmPmR
C     to stay consitent with volume change (=d/dt etaH).
# ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE PmEpR = comlev1, key = ikey_dynamics,  kind = isbyte
CADJ STORE EmPmR = comlev1, key = ikey_dynamics,  kind = isbyte
# endif
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        IF ( staggerTimeStep ) THEN
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           PmEpR(i,j,bi,bj) = -EmPmR(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO
#endif /* EXACT_CONSERV */

C--   set surfaceForcingT,S to zero.
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           surfaceForcingT(i,j,bi,bj) = 0. _d 0
           surfaceForcingS(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Start with surface restoring term :
      IF ( doThetaClimRelax .OR. doSaltClimRelax ) THEN
       CALL FORCING_SURF_RELAX(
     I              iMin, iMax, jMin, jMax,
     I              myTime, myIter, myThid )
      ENDIF

#ifdef ALLOW_PTRACERS
C--   passive tracer surface forcing:
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE surfaceForcingS = comlev1, key = ikey_dynamics,
CADJ &    kind = isbyte
#endif
      IF ( usePTRACERS ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         CALL PTRACERS_FORCING_SURF(
     I        surfaceForcingS(1-OLx,1-OLy,bi,bj),
     I        bi, bj, iMin, iMax, jMin, jMax,
     I        myTime, myIter, myThid )
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_PTRACERS */

C- Notes: setting of PmEpR and pTracers surface forcing could have been
C         moved below, inside a unique bi,bj block. However this results
C         in tricky dependencies for TAF (and recomputations).
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

#ifdef ALLOW_AUTODIFF_TAMC
          act1 = bi - myBxLo(myThid)
          max1 = myBxHi(myThid) - myBxLo(myThid) + 1
          act2 = bj - myByLo(myThid)
          max2 = myByHi(myThid) - myByLo(myThid) + 1
          act3 = myThid - 1
          max3 = nTx*nTy
          act4 = ikey_dynamics - 1
          ikey = (act1 + 1) + act2*max1
     &                      + act3*max1*max2
     &                      + act4*max1*max2*max3
#endif /* ALLOW_AUTODIFF_TAMC */

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE EmPmR(:,:,bi,bj) = comlev1_bibj, key=ikey,  kind = isbyte
#ifdef EXACT_CONSERV
CADJ STORE PmEpR(:,:,bi,bj) = comlev1_bibj, key=ikey,  kind = isbyte
#endif
#endif /* ALLOW_AUTODIFF_TAMC */

C--   Surface Fluxes :
        DO j = jMin, jMax
         DO i = iMin, iMax

C     Zonal wind stress fu:
          surfaceForcingU(i,j,bi,bj) = fu(i,j,bi,bj)*mass2rUnit
C     Meridional wind stress fv:
          surfaceForcingV(i,j,bi,bj) = fv(i,j,bi,bj)*mass2rUnit
C     Net heat flux Qnet:
          surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &       - ( Qnet(i,j,bi,bj)
#ifdef SHORTWAVE_HEATING
     &          -Qsw(i,j,bi,bj)
#endif
     &         ) *recip_Cp*mass2rUnit
C     Net Salt Flux :
          surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &      -saltFlux(i,j,bi,bj)*mass2rUnit

         ENDDO
        ENDDO

#ifdef ALLOW_SALT_PLUME
C saltPlume is the amount of salt rejected by ice while freezing;
C it is here subtracted from surfaceForcingS and will be redistributed
C to multiple vertical levels later on as per Duffy et al. (GRL 1999)
C-- for the case of SALT_PLUME_VOLUME, need to call this S/R right
C-- before kpp in do_oceanic_phys.F due to recent moved of
C-- external_forcing_surf.F outside bi,bj loop.
#ifndef SALT_PLUME_VOLUME
        IF ( useSALT_PLUME ) THEN
         CALL SALT_PLUME_FORCING_SURF(
     I        bi, bj, iMin, iMax, jMin, jMax,
     I        myTime, myIter, myThid )
        ENDIF
#endif /* SALT_PLUME_VOLUME */
#endif /* ALLOW_SALT_PLUME */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Fresh-water flux

C-    Apply mask on Fresh-Water flux (if useRealFreshWaterFlux)
C     <== removed: maskInC is applied directly in S/R SOLVE_FOR_PRESSURE

#ifdef EXACT_CONSERV
      IF ( (nonlinFreeSurf.GT.0 .OR. usingPCoords)
     &     .AND. useRealFreshWaterFlux ) THEN

C--   NonLin_FrSurf and RealFreshWaterFlux : PmEpR effectively changes
C     the water column height ; temp., salt, (tracer) flux associated
C     with this input/output of water is added here to the surface tendency.

       IF (temp_EvPrRn.NE.UNSET_RL) THEN
        DO j = jMin, jMax
         DO i = iMin, iMax
          surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &      + PmEpR(i,j,bi,bj)
     &          *( temp_EvPrRn - theta(i,j,ks,bi,bj) )
     &          *mass2rUnit
         ENDDO
        ENDDO
       ENDIF

       IF (salt_EvPrRn.NE.UNSET_RL) THEN
        DO j = jMin, jMax
         DO i = iMin, iMax
          surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &      + PmEpR(i,j,bi,bj)
     &          *( salt_EvPrRn - salt(i,j,ks,bi,bj) )
     &          *mass2rUnit
         ENDDO
        ENDDO
       ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      ELSE
#else /* EXACT_CONSERV */
      IF (.TRUE.) THEN
#endif /* EXACT_CONSERV */

C--   EmPmR does not really affect the water column height (for tracer budget)
C     and is converted to a salt tendency.

       IF (convertFW2Salt .EQ. -1.) THEN
C-    use local surface tracer field to calculate forcing term:

        IF (temp_EvPrRn.NE.UNSET_RL) THEN
C     account for Rain/Evap heat content (temp_EvPrRn) using local SST
         DO j = jMin, jMax
          DO i = iMin, iMax
           surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &       + EmPmR(i,j,bi,bj)
     &           *( theta(i,j,ks,bi,bj) - temp_EvPrRn )
     &           *mass2rUnit
          ENDDO
         ENDDO
        ENDIF
        IF (salt_EvPrRn.NE.UNSET_RL) THEN
C     converts EmPmR to salinity tendency using surface local salinity
         DO j = jMin, jMax
          DO i = iMin, iMax
           surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &       + EmPmR(i,j,bi,bj)
     &           *( salt(i,j,ks,bi,bj) - salt_EvPrRn )
     &           *mass2rUnit
          ENDDO
         ENDDO
        ENDIF

       ELSE
C-    use uniform tracer value to calculate forcing term:

        IF (temp_EvPrRn.NE.UNSET_RL) THEN
C     account for Rain/Evap heat content (temp_EvPrRn) assuming uniform SST (=tRef)
         DO j = jMin, jMax
          DO i = iMin, iMax
           surfaceForcingT(i,j,bi,bj) = surfaceForcingT(i,j,bi,bj)
     &       + EmPmR(i,j,bi,bj)
     &           *( tRef(ks) - temp_EvPrRn )
     &           *mass2rUnit
          ENDDO
         ENDDO
        ENDIF
        IF (salt_EvPrRn.NE.UNSET_RL) THEN
C     converts EmPmR to virtual salt flux using uniform salinity (default=35)
         DO j = jMin, jMax
          DO i = iMin, iMax
           surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &       + EmPmR(i,j,bi,bj)
     &           *( convertFW2Salt - salt_EvPrRn )
     &           *mass2rUnit
          ENDDO
         ENDDO
        ENDIF

C-    end local-surface-tracer / uniform-value distinction
       ENDIF

      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ATMOSPHERIC_LOADING
C-- Atmospheric surface Pressure loading : added to phi0surf when using Z-coord;
C   Not yet implemented for Ocean in P: would need to be applied to the other end
C   of the column, as a vertical velocity (omega); (meaningless for Atmos in P).
C- Note:
C   Using P-coord., a hack (now directly applied from S/R INI_FORCING)
C   is sometime used to read phi0surf from a file (pLoadFile) instead
C   of computing it from bathymetry & density ref. profile.

        IF ( usingZCoords ) THEN
C   The true atmospheric P-loading is not yet implemented for P-coord
C   (requires time varying dP(Nr) like dP(k-bottom) with NonLin FS).
         IF ( useRealFreshWaterFlux ) THEN
          DO j = jMin, jMax
           DO i = iMin, iMax
            phi0surf(i,j,bi,bj) = ( pLoad(i,j,bi,bj)
     &                          +sIceLoad(i,j,bi,bj)*gravity
     &                            )*recip_rhoConst
           ENDDO
          ENDDO
         ELSE
          DO j = jMin, jMax
           DO i = iMin, iMax
            phi0surf(i,j,bi,bj) = pLoad(i,j,bi,bj)*recip_rhoConst
           ENDDO
          ENDDO
         ENDIF
c       ELSEIF ( usingPCoords ) THEN
C-- This is a hack used to read phi0surf from a file (pLoadFile)
C   instead of computing it from bathymetry & density ref. profile.
C   ==> now done only once, in S/R INI_FORCING
c         DO j = jMin, jMax
c          DO i = iMin, iMax
c           phi0surf(i,j,bi,bj) = pLoad(i,j,bi,bj)
c          ENDDO
c         ENDDO
        ENDIF
#endif /* ATMOSPHERIC_LOADING */

#ifdef ALLOW_SHELFICE
        IF ( useSHELFICE) THEN
          CALL SHELFICE_FORCING_SURF(
     I                  bi, bj, iMin, iMax, jMin, jMax,
     I                  myTime, myIter, myThid )
        ENDIF
#endif /* ALLOW_SHELFICE */

C--   end bi,bj loops.
       ENDDO
      ENDDO

      RETURN
      END
