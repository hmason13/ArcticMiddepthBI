#include "CPP_EEOPTIONS.h"

C--   File global_sum_tile_merid.F: Routines that perform global sum
C                             on a tile array
C      Contents
C      o GLOBAL_SUM_TILE_MERID
C      Author
C      Adapted from original global_sum_tile.F by hassan, NYU 
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: GLOBAL_SUM_TILE_MERID

C     !INTERFACE:
      SUBROUTINE GLOBAL_SUM_TILE_MERID(
     I                       phiTile,
     O                       sumPhi,
     I                       myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GLOBAL\_SUM\_TILE\_MERID
C     | o Handle sum for _RL data.
C     *==========================================================*
C     | Apply sum on an array of one value per tile
C     |  and operate over all tiles & all the processes.
C     *==========================================================*
C     | hassan: A modification to the original GLOBAL_SUM_TILE_RL
C     |  that returns a vector of values. This vector is the sum
C     |  over a column (meridional) of tiles. Each entry
C     |  represents the sum for each column of tiles. The total
C     |  number of entries is given by nSx*nPx.
C     |  
C     |  nSx*Npx must be less than 31! (lshare8 - 1)
C     | 
C     |  My CPP_EEOPTIONS.h file shows:
C     |    #define GLOBAL_SUM_ORDER_TILES
C     |    #undef  GLOBAL_SUM_SEND_RECV
C     |    #define ALLOW_USE_MPI
C     |  So I have modified only the relevant block of code.
C     |  Ensure your settings are the same, or write your own!
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "GLOBAL_SUM.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     phiTile :: Input array with one value per tile
C     sumPhi  :: Result of sum.
C     myThid  :: My thread id.
      _RL     phiTile(nSx,nSy)
      _RL     sumPhi(lshare8)
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     bi,bj   :: Loop counters
C     mpiRC   :: MPI return code
C- type declaration of: sumMyPr, sumAllP, localBuf and shareBufGSR8 :
C         all 4 needs to have the same length as MPI_DOUBLE_PRECISION
      INTEGER bi,bj
#ifdef ALLOW_USE_MPI
#ifdef GLOBAL_SUM_SEND_RECV
      INTEGER biG, bjG, np, pId
      INTEGER lbuff, idest, itag, ready_to_receive
      INTEGER istatus(MPI_STATUS_SIZE), ierr
      Real*8  localBuf (nSx,nSy)
      Real*8  globalBuf(nSx*nPx,nSy*nPy)
#elif defined (GLOBAL_SUM_ORDER_TILES)
      INTEGER biG, bjG, lbuff
      Real*8  localBuf (nSx*nPx,nSy*nPy)
      Real*8  globalBuf(nSx*nPx,nSy*nPy)
#endif
      INTEGER mpiRC
#endif /* ALLOW_USE_MPI */
      Real*8  sumMyPr
      Real*8  sumAllP
C     hassan: the new variable that will be returned
      Real*8  sumAllVec(nPx*nSx)
      INTEGER indexx
CEOP

C     this barrier is not necessary:
c     CALL BAR2( myThid )

C--   write local sum into shared-buffer array
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
         shareBufGSR8(bi,bj) = phiTile(bi,bj)
       ENDDO
      ENDDO

C--   Master thread cannot start until everyone is ready:
      CALL BAR2( myThid )
      _BEGIN_MASTER( myThid )

#if (defined (GLOBAL_SUM_SEND_RECV) && defined (ALLOW_USE_MPI) )
      IF ( usingMPI ) THEN

        lbuff = nSx*nSy
        idest = 0
        itag  = 0
        ready_to_receive = 0

        IF ( mpiMyId.NE.0 ) THEN

C--   All proceses except 0 wait to be polled then send local array
#ifndef DISABLE_MPI_READY_TO_RECEIVE
          CALL MPI_RECV (ready_to_receive, 1, MPI_INTEGER,
     &         idest, itag, MPI_COMM_MODEL, istatus, ierr)
#endif
          CALL MPI_SEND (shareBufGSR8, lbuff, MPI_DOUBLE_PRECISION,
     &         idest, itag, MPI_COMM_MODEL, ierr)

C--   All proceses except 0 receive result from process 0
          CALL MPI_RECV (sumAllP, 1, MPI_DOUBLE_PRECISION,
     &         idest, itag, MPI_COMM_MODEL, istatus, ierr)

        ELSE
C-      case mpiMyId = 0

C--   Process 0 fills-in its local data
         np = 1
         DO bj=1,nSy
          DO bi=1,nSx
            biG = (mpi_myXGlobalLo(np)-1)/sNx+bi
            bjG = (mpi_myYGlobalLo(np)-1)/sNy+bj
            globalBuf(biG,bjG) = shareBufGSR8(bi,bj)
          ENDDO
         ENDDO

C--   Process 0 polls and receives data from each process in turn
         DO np = 2, nPx*nPy
          pId = np - 1
#ifndef DISABLE_MPI_READY_TO_RECEIVE
          CALL MPI_SEND (ready_to_receive, 1, MPI_INTEGER,
     &           pId, itag, MPI_COMM_MODEL, ierr)
#endif
          CALL MPI_RECV (localBuf, lbuff, MPI_DOUBLE_PRECISION,
     &           pId, itag, MPI_COMM_MODEL, istatus, ierr)

C--   Process 0 gathers the local arrays into a global array.
          DO bj=1,nSy
           DO bi=1,nSx
             biG = (mpi_myXGlobalLo(np)-1)/sNx+bi
             bjG = (mpi_myYGlobalLo(np)-1)/sNy+bj
             globalBuf(biG,bjG) = localBuf(bi,bj)
           ENDDO
          ENDDO
C-       end loop on np
         ENDDO

C--   Sum over all tiles:
         sumAllP = 0.
         DO bjG = 1,nSy*nPy
          DO biG = 1,nSx*nPx
           sumAllP = sumAllP + globalBuf(biG,bjG)
          ENDDO
         ENDDO

C--   Process 0 sends result to all other processes
         lbuff = 1
         DO np = 2, nPx*nPy
          pId = np - 1
          CALL MPI_SEND (sumAllP, 1, MPI_DOUBLE_PRECISION,
     &                   pId, itag, MPI_COMM_MODEL, ierr)
         ENDDO

C       End if/else mpiMyId = 0
        ENDIF

      ELSE
#elif (defined (GLOBAL_SUM_ORDER_TILES) && defined (ALLOW_USE_MPI) )
      IF ( usingMPI ) THEN

C--   Initialise local buffer
        DO bjG=1,nSy*nPy
         DO biG=1,nSx*nPx
           localBuf(biG,bjG) = 0.
         ENDDO
        ENDDO

C--   Put my own data in local buffer
        DO bj=1,nSy
         DO bi=1,nSx
           biG = (myXGlobalLo-1)/sNx+bi
           bjG = (myYGlobalLo-1)/sNy+bj
           localBuf(biG,bjG) = shareBufGSR8(bi,bj)
         ENDDO
        ENDDO

C--   Collect data from all procs
        lbuff = nSx*nPx*nSy*nPy
        CALL MPI_Allreduce( localBuf, globalBuf, lbuff,
     &           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_MODEL, mpiRC )

C-    hassan: modified to sum over tiles in the same column
C-            luckily, the biG,bjG are based on the tile spatial location
C-            within the simulation!
C--   Sum over each column:
         DO biG = 1,nSx*nPx
          sumAllVec(biG) = 0.
         ENDDO

         DO biG = 1,nSx*nPx
          DO bjG = 1,nSy*nPy
           sumAllVec(biG) = sumAllVec(biG) + globalBuf(biG,bjG)
          ENDDO
         ENDDO

      ELSE
#else /* not ((GLOBAL_SUM_SEND_RECV | GLOBAL_SUM_ORDER_TILES) & ALLOW_USE_MPI) */
      IF ( .TRUE. ) THEN
#endif /* not ((GLOBAL_SUM_SEND_RECV | GLOBAL_SUM_ORDER_TILES) & ALLOW_USE_MPI) */

C--   Sum over all tiles (of the same process) first
        sumMyPr = 0.
        DO bj = 1,nSy
         DO bi = 1,nSx
          sumMyPr = sumMyPr + shareBufGSR8(bi,bj)
         ENDDO
        ENDDO

C     in case MPI is not used:
        sumAllP = sumMyPr

#ifdef ALLOW_USE_MPI
        IF ( usingMPI ) THEN
         CALL MPI_Allreduce(sumMyPr,sumAllP,1,MPI_DOUBLE_PRECISION,
     &                      MPI_SUM,MPI_COMM_MODEL,mpiRC)
        ENDIF
#endif /* ALLOW_USE_MPI */

      ENDIF
      
C-    hassan: This code block must also be changed.
C-            The vector sumAllVec holds max 32 values. It is passed to the
C-            preexisting array phiGSR8 which has size (lshare8,MAX_NO_THREADS).
C             lshare8 is 32, and the first slot is used by other subroutines,
C             so we limit ourselves to the following 31 slots. It may be
C             possible to use a larger shared array like phivGSR8, but then
C             the localBuf and globalBuf arrays would have to grow as well.
C             If performance is terrible (too large sNx?) then I suggest
C             using the larger shared array to reduce the # of MPI calls. Of
C             course, this isn't easy because if it were I'd have done it.
C-             
C--   Write solution to shared buffer (all threads can see it)
c     shareBufGSR8(1,1) = sumAllP
      DO indexx = 1,nSx*nPx
        phiGSR8(indexx+1,0) = sumAllVec(indexx)
      ENDDO

      _END_MASTER( myThid )
C--   Everyone wait for Master thread to be ready
      CALL BAR2( myThid )

C--   set result for every threads
c     sumPhi = shareBufGSR8(1,1)
      DO indexx = 1,nSx*nPx
        sumPhi(indexx) = phiGSR8(indexx+1,0)
      ENDDO
C--   A barrier was needed here to prevent thread 1 to modify shareBufGSR8(1,1)
C     (as it would in the following call to this S/R) before all threads get
C     their global-sum result out.
C     No longer needed since a dedicated shared var. is used to share the output
c     CALL BAR2( myThid )

      RETURN
      END
