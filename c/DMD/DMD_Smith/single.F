C<PRE>
C***  Discontinuous MD for hard or square-well chains
C***  Algorithm implements the single-event scheduling of Allen.
C***  A combination neighbor list/link list is used.
C***  Events are scheduled with a simple binary tree.
C***  Particles only experience periodic boundaries
C***  when interactions are considered. (Infinite 3D checkerboard)
C***  Densities are assumed to be liquid-like so the
C***  minimum image convention holds.
C***  Periodic Boundary condition assumes a box length of 1
C***  Box length = 1 and sigma is scaled accordingly.
C***  Box is assumed to be centered at the origin (-0.5,+0.5)
C***  Most suited for systems with more the 500 particles.
C***  Square well interaction requires an extra input file.
C***  For canonical square, a scaling technique must be added.

C*** The following defines will be expanded by your C preprocessor
C*** this usually happens automatically by naming the file with 
C*** a .F extension.

C*** C preprocessor defines
C*** Number of particles
#define NOP 512
C*** To use or not to use Belleman's type bond
#define BELMAN
C*** Maximum number of neighbors
#define MAXNAB 60
C*** Average subcell occupany for link lists
#define Cell_Occ 8.0
C*** Number of subcells for inadequate f77
#define SMARTF77
#ifndef SMARTF77
#define Number_of_cells 4	
#endif
C*** To add or not a square-well interaction, need the sw.cnf file
C#define SQWELL
C***  Variables
C***  NOP			Number of particles
C***  VX(NOP), VY(NOP), VZ(NOP)	Particle velocities
C***  RX(NOP), RY(NOP), RZ(NOP)	Particle positions
C***  COLTIM                   	Time to next collision
C***  COLTYP			Collision type
C***  PARTNR		        Collision partner
C***  SIGMA			Particle/Segment diamter
C***  ETA                       Volume Fraction 
C***  TEMP                      kT/m
C***  DEL		        Bond factor in Rapaport model 
C***  BOXL		        Box length  
  
      PROGRAM DMD_SINGLE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER PARTNR, COLTYP
      COMMON /TIMES / COLTIM(NOP), PARTNR(NOP), COLTYP(NOP)
      COMMON /BONDS / ISBND(0:NOP)

#include "header.f"

C**   Neighbor arrays and parameters
      INTEGER     MAXARR
      PARAMETER ( MAXARR = MAXNAB*NOP )
#ifdef SMARTF77
C** Some FORTRAN compilers may have trouble with this
      PARAMETER( M = NINT(NOP/Cell_Occ)**(0.3333333333) )
#else
      PARAMETER( M = Number_of_cells )
#endif
      INTEGER*2  UPNAB(MAXARR), DNNAB(MAXARR)
      INTEGER*2  NNABUP(NOP), NNABDN(NOP)
      INTEGER  INDEX(NOP)
      COMMON /NABORS / UPNAB, NNABUP, DNNAB, NNABDN, INDEX

C** Dimension array for storing which particles to update
      INTEGER*2 NN(MAXNAB)

      DOUBLE PRECISION  LSTAR, KE
      LOGICAL   CANON, LOOP, OVRLAP
      INTEGER   ACCESS, IERRNO
      REAL      TARRAY(2), DELTA, EXTIME
      CHARACTER CNFILE*15, FDATE*24, RUNID*1, CF*16
      CHARACTER CHKFIL*19, FRSTLET*1,CFTRAJ*16, NXRUN*1

#ifdef CRAY
      EXTERNAL SECOND,DATE,CLOCK,ACCESS
      CHARACTER*10 DDAY,STIM
      CHARACTER ENVNAM*4, HSTNAM*20
#else 
      INTEGER HOSTNM, IHOST
      CHARACTER HSTNAM*20
      REAL DTIME
      EXTERNAL FDATE,DTIME,ACCESS,HOSTNM,FLUSH
#endif

      OPEN(UNIT=5,FILE='md.cnf',STATUS='UNKNOWN',FORM='FORMATTED')
C***  Read Packing Fraction
      READ(5,'(/F10.5)') ETA
C***  Read Chain Length
      READ(5,'(/I3)') CHNLEN
C***  Read Total Number of Collisions Desired
      READ(5,'(/F10.5)') RNCOLL
C***  Read Total Number of Batches
      READ(5,'(/I2)') NBATCH
C***  Read Initial Configuration Filename
      READ(5,'(/A15)') CNFILE
      READ(5,'(/F10.5)') DELFAC
      READ(5,'(/I10)') NABLIM
      READ(5,'(/F10.5)') SETEMP
      READ(5,'(/I10)') MODUPD
      CLOSE(UNIT=5)

      NCOLL = DINT(RNCOLL*1000000)
      NCOLL = NCOLL/NBATCH

#ifdef SQWELL
      OPEN(UNIT=95,FILE='sw.cnf',STATUS='UNKNOWN',FORM='FORMATTED')
      READ(95,'(/ F10.5)') 	GAMMA
      READ(95,'(/ F10.5)') 	V0
      READ(95,'(/ I10)') 	MODTMP

      IF(MODTMP .EQ. 0) THEN
C** Never check
         MODTMP = NCOLL+1
      ELSE
         IF(MOD(MODTMP,MODUPD) .NE. 0) THEN
	  WRITE(6,*) " The frequency of energy calculation must be"
	  WRITE(6,*) " a multiple of the position update (MODUPD)"
	  STOP
	 ENDIF
      END IF
#endif

      ICFRST=ICHAR('A')-1
      LOOP=.TRUE.
C** Find the first unused identifier
       DO WHILE(LOOP)
	 ICFRST=ICFRST+1
         FRSTLET=CHAR(ICFRST)
         CHKFIL=CNFILE//FRSTLET//'.md'
         IERRNO=ACCESS(CHKFIL,'r')
	 IF(IERRNO .NE. 0)LOOP=.FALSE.
      END DO

C**** Do several batches at once ***************************

      DO 10000 NBAT=ICFRST,ICFRST+NBATCH-1
         RUNID=CHAR(NBAT)
         NXRUN=CHAR(NBAT+1)
         CF=CNFILE//RUNID

C***  Initialize counters
         COLL=0
         NNEIGH=0
         NINTER=0
         NINTRA=0
         NBOND=0
         NNAB=0
         NWELL=0

C***  Zero the Virial Accumulator and time
         PSTOT=0.D0
         T=0.D0
         TNEXT = 0.D0

C***  Read in initial configuration
         CFTRAJ=CNFILE//CHAR(NBAT -1)
         CALL READCF(CF)

C***  Open output files
C***  Write output as input for next run in multi-batch mode
         IF(NBATCH .EQ. 1) THEN
            OPEN(UNIT=4,FILE=CF//'.out',STATUS='UNKNOWN',
     &         	 FORM='FORMATTED')
         ELSE
            OPEN(UNIT=4,FILE=CNFILE//NXRUN//'.in',STATUS='UNKNOWN',
     &           FORM='FORMATTED')
         END IF

C**  Log file and debug info
         OPEN(UNIT=6,FILE=CF//'.md',STATUS='UNKNOWN',FORM='FORMATTED')
C         OPEN(UNIT=99,FILE=CF//'.dbg',STATUS='UNKNOWN',FORM='FORMATTED')

C***  Set chain number for each bead
         CALL CHNNUM(CHNLEN)

C***  Initial Parameters
         NCHN = NOP/CHNLEN
         RCHLEN = DBLE(CHNLEN)
         PI = 4.D0*DATAN(1.D0)
         DENSTY=6.D0/PI*ETA

         SIGMA = (DENSTY/DBLE(NOP))**(1.D0/3.D0)
         SIGSQ=SIGMA*SIGMA
         BONDL = SIGMA*(1.D0 + DELFAC)
         BNDLSQ=BONDL*BONDL
         SIGSQIJ=(SIGMA*(1.D0 - DELFAC))**2.D0

         BOXL = 1.D0
         BOXLM1 = 1.D0/BOXL
         BOXHLF = BOXL/2.D0
         VOL = BOXL*BOXL*BOXL

#ifdef SQWELL
         TSTAR=SETEMP
         SIG2=GAMMA*SIGMA
         SIG2SQ=SIG2*SIG2
         DELPE=4.D0*V0*SIG2SQ
#endif

C*** Initialize maps for link lists
         CELLI = DBLE(M)
         CALL MAPS

C***  Rescale velocities to insure zero total momentum and set temp
         CALL MOMENT(SETEMP)

C***  Unfold chains from the effects of periodic boundaries
         CALL DELPBC(CHNLEN,BOXL,BOXLM1)

C***  Check for overlap and calculate energy
         CALL CHECK(KE,OVRLAP)
         IF (OVRLAP)STOP 'Particle overlap in initial Configuration'
         EN=KE/DBLE(NOP)
         TEMP = 2.D0*EN/3.D0

#ifdef SQWELL
         OPEN(UNIT=9,FILE=CF//'.ener',STATUS='UNKNOWN',FORM='FORMATTED')
	 SCLTMP=DSQRT(V0*TSTAR*3.D0*REAL(NOP)/(2.D0*KE))
	 DO I=1,NOP
	  VX(I)=SCLTMP*VX(I)
	  VY(I)=SCLTMP*VY(I)
	  VZ(I)=SCLTMP*VZ(I)
         END DO
         TEMP = 2.D0*EN/3.D0
         TSTACT=TEMP/V0
	 TOTE = KE + PE
         WRITE(9,'("Collision   Kin. Ener.  Pot. Ener. Tot. Ener. T*")')
         WRITE(9,'(I10,4F10.3)')COLL,KE,PE,TOTE,TSTACT
#endif

         ENKT = EN/TEMP
         PSFACT=3.D0*TEMP*DBLE(NCHN)*SIGMA

#ifdef CRAY
         ENVNAM = 'HOST'
         CALL GETENV(ENVNAM, HSTNAM)
         CALL DATE(DDAY)
         CALL CLOCK(STIM)
         WRITE(6,'("  This simulation began on    ",A24)')DDAY
         WRITE(6,'("                        at    ",A24)')STIM
#else 
         WRITE(6,'("  This simulation began on    ",A24)') FDATE()
         IHOST = HOSTNM (HSTNAM)
#endif
         WRITE(6,'("  Host name is   ",A20/)')HSTNAM
         WRITE(6,'("  Reduced Density is          ",F10.7)') DENSTY
         WRITE(6,'("  Sigma is                    ",E13.6)') SIGMA
         WRITE(6,'("  Box Length is               ",E13.6)') BOXL
#ifdef BELMAN
         WRITE(6,'("  Compiled for BELLEMANS BOND !!!!")')
#endif
         WRITE(6,'("  Bond Length Factor is       ",F10.7)') DELFAC
         WRITE(6,'("  Packing Fraction is         ",F10.7)') ETA
         WRITE(6,'("  Run length (MCollisions) is ",I10)') NCOLL/1000000
         WRITE(6,'("  Configuration Filename is   ",A15)') CNFILE
         WRITE(6,'("  Temperature kT/m is         ",F20.16)')TEMP
         WRITE(6,'("  Initial E/NkT is            ",F10.7)') ENKT
#ifdef SQWELL
         WRITE(6,'(" Compiled for SQUARE WELL SIMULATION ")')
         WRITE(6,'("  Reduced Temp(kT/eps) is     ",F10.7)') TSTAR
         WRITE(6,'("  Well Depth epsilon/m is    ",F10.7)') V0
         WRITE(6,'("  Well Width gamma is        ",F10.7)') GAMMA
         WRITE(6,'("  Initial Total Energy is    ",F15.8)') TOTE
#endif

         WRITE(6,'("  Number of neighbors         ",I5)') NABLIM
         WRITE(6,'("  Number of subcells          ",I5)') M*M*M
         WRITE(6,'("  Position updates every      ",I5/)') MODUPD
	 IF(MOD(NCOLL,MODUPD) .NE. 0) THEN
         WRITE(6,*) " Number of collisions must a multiple of"
         WRITE(6,*) " position update interval"
	 Stop
	 END IF
C***  Set up neighbor lists
#ifdef SQWELL
         SIGCUT = SIG2
         RLIST  = 2.D0*SIGCUT
#else 
         SIGCUT = SIGMA
         RLIST  = 3.D0*SIGCUT
#endif
         RLSTSQ = RLIST*RLIST
         DSPLIM = ((RLIST - SIGCUT)/2.5D0)**2.D0

C***  Set up indices for down neighbor lists
         DO I = 1,NOP
            INDEX(I) = (I-1)*MAXNAB
         END DO

C***  Set up initial collision lists and Partners **
         DO 10 I=1,NOP
            COLTIM(I) = 1D10
            PARTNR(I) = NOP
            COLTYP(I) = NOP
 10      CONTINUE

C***  Set an integer to decide if the bead has a bond time
         DO I = 1,NOP
            ISBND(I) = 1
         END DO

         DO I = CHNLEN,NOP,CHNLEN
            ISBND(I) = 0
         END DO

C***  Set up neighbor lists and initial time lists
         CALL NABCHK
         CELL = BOXL/CELLI
         IF(CELL .LT. RLIST) THEN
C***  Subcell length CELL must be  greater than neighbor list radius RLIST
            WRITE(6,'("Subcell size too small",F10.4)')CELL/SIGCUT
            STOP
         END IF

         WRITE(6,'(/"  START OF DYNAMICS")')
         CALL FLUSH(6)
#ifndef CRAY
         DELTA=DTIME(TARRAY)
#endif

C=================Main Loop=========================================

         DO 1000 COLL = 1,NCOLL

            CALL NEXTEV(I)
            TNEXT  = COLTIM(I)
            J      = PARTNR(I) 
            ICOTYP = COLTYP(I)

            IF(ICOTYP .EQ. 1) THEN
               IF(ICHN(I) .NE. ICHN(J))THEN
                  NINTER=NINTER+1
               ELSEIF((J-I) .EQ. 1)THEN
                  NNEIGH=NNEIGH + 1
               ELSE
                  NINTRA=NINTRA + 1
               END IF
            ELSEIF(ICOTYP .EQ. 0)THEN
C**   A bond stretch
               NBOND = NBOND + 1
            END IF
#ifdef SQWELL
            IF(ICOTYP .GE. 2)THEN
               NWELL = NWELL + 1
            END IF
#endif

C**   A safety net, just to be sure 
            IF(TNEXT .LE. 0.D0) THEN
               WRITE(6,*)'***NEGATIVE OR ZERO COLLISION TIME***'
               WRITE(6,*)'Coll#, I,  J,  Coltype, Tnext'
               WRITE(6,*)COLL, I, J, ICOTYP, TNEXT
               STOP
            END IF

C***  Compute collision dynamics and collisional virial
            CALL BUMP(I,J)	
            PSTOT = PSTOT - W

C***  Check to see if the neighbor lists need updating
            IF(MOD(COLL,MODUPD) .EQ. 0)THEN
C***  Do a big update to put all the particles in true positions
C***  Move particles forward by TNEXT and reduce collision times in list
               T = T + TNEXT
               DO 490 K=1,NOP
                  COLTIM(K) = COLTIM(K) - TNEXT
                  RX(K) = RX(K) + VX(K) * TNEXT
                  RY(K) = RY(K) + VY(K) * TNEXT
                  RZ(K) = RZ(K) + VZ(K) * TNEXT
 490           CONTINUE
               TNEXT = 0.D0

C***  Check neighbor list expiration and update time lists
C***  if neighbors are renewed
               CALL NABCHK
#ifdef SQWELL
C   ** Calcualte the temperature every MODTMP collisions
	IF(MOD(COLL,MODTMP) .EQ. 0) THEN
	 CALL CHECK(KE,OVRLAP)
	 TSTACT=2.D0*KE/3.D0/REAL(NOP)/V0
         TOTE = KE + PE
	 WRITE(9,'(I10,4F10.3)')COLL,KE,PE,TOTE,TSTACT
	END IF
#endif
            END IF

C***  Find new core collision time for those particles affected 
C***  by the previous collision, since the time list is constructed
C***  so PARTNR(I) is greater than I we only need to check the down neighbors
C***  of I and J. I will automatically be included since it is
C***  a down neighbor of J. We add J to the list as well for
C***  convenience.

            NUP = 1
            NN(NUP) = J

C***  Check the down neighbors of I
            KBEG = INDEX(I) + 1
            KEND = KBEG + NNABDN(I) - 1

C***  Cray Directive
CDIR$ IVDEP
            DO 449 KK=KBEG, KEND
               K = DNNAB(KK)
               IF(PARTNR(K) .EQ. I)THEN
                  NUP = NUP + 1
                  NN(NUP) = K
               END IF
 449        CONTINUE

            KBEG = INDEX(J) + 1
            KEND = KBEG + NNABDN(J) - 1

C***  Check the down neighbors of J
CDIR$ IVDEP
            DO 450 KK=KBEG, KEND
               K = DNNAB(KK)
               IF(PARTNR(K) .EQ. J)THEN
                  NUP = NUP + 1
                  NN(NUP) = K
               END IF
 450        CONTINUE


C***  Now update only those who need it
            IMINUS1 = 0
            DO 451 L = 1, NUP
               K = NN(L)
               IF(K .EQ. I-1)IMINUS1 = 1
               CALL UPLIST(K)
C***  Update bond time since only one time list is used
               CALL BTIME(K)
C***  Now schedule the new event
               CALL ADDTRE(K)
 451        CONTINUE

C***  Check the down neighbors for new collisions and schedule them
            CALL DNLIST(I)
            CALL DNLIST(J)

C***  Update bond times for every particle that was
C***  affected by the velocity change in I and J. Must
C***  be done after DNLIST(I,J) since that routine 
C***  may change the next collision time for the same particles.
C***  Since this may (emphasis) have been done in
C***  loop 451 we use some logic to avoid doing it again.
            IF(IMINUS1 .EQ. 0)THEN
               CALL BTIME(I-1)
               IF(I.NE.1)CALL ADDTRE(I-1)
            END IF

C***  When I and J are bonded then J-1 equals I which has already been updated 
            IF(J-1 .NE. I) THEN
               CALL BTIME(J-1)
               CALL ADDTRE(J-1)
            END IF

 1000    CONTINUE

C=======================Main Loop Ends=================================

#ifdef CRAY
         DELTA=SECOND()
#else 
         DELTA=DTIME(TARRAY)
#endif

 1300    WRITE(6,'(/"  END OF DYNAMICS")')

C***  Check for Overlap
         CALL CHECK(KE,OVRLAP)
         IF (OVRLAP) THEN
            WRITE(6,'(" Particle overlap in final configuration")')
         END IF

C***  Compute average bond distance and scale for SIGMA=1
         IF(CHNLEN .GT. 1) THEN
            CALL BNDLNGTH(LSTAR,CHNLEN)
            LSTAR=LSTAR/SIGMA
         END IF

C***  Write final configuration
         CALL WRITCF (4, BOXL)

C***  Write out interesting info.
         PVNKT = PSTOT/DBLE(NCHN)/3.D0/T/TEMP
         KE    = KE/DBLE(NOP)
         TEMP  = 2.D0*KE/3.D0
         ENKT  = KE/TEMP
         TN    = T*DSQRT(TEMP)/SIGMA
         RATE  = DBLE(NCOLL)/TN
         TBC   = DBLE(NOP)/RATE/2.D0
         TRUTBC= 1.D0/RATE
         NTOTAL= NBOND+NINTRA+NINTER+NNEIGH

#ifdef SQWELL
         TSTAR = TEMP/V0
         NTOTAL= NTOTAL+NWELL
         REWIND 9
         READ(9,*)
         II =0
         SUMT = 0.D0
 1400    READ(9,*,END=1401)ICOLL,EN1,EN2,EN3,TACT
         II = II+1
         SUMT = SUMT + TACT
         GOTO 1400
 1401    TAVER = SUMT/DBLE(II)
         WRITE(6,'("  Average Reduced Temp is          ",F20.16)') TAVER
         WRITE(6,'("  Number of Well Collisions        ",I10)') NWELL
         WRITE(6,'("  Final Total Energy is            ",F15.8)') TOTE
#endif

         ZTOT=PVNKT + DBLE(CHNLEN)

C***  Check total momentum again
         SUMVX=0.D0
         SUMVY=0.D0
         SUMVZ=0.D0

         DO 2000 I=1,NOP
            SUMVX=SUMVX + VX(I)
            SUMVY=SUMVY + VY(I)
            SUMVZ=SUMVZ + VZ(I)
 2000    CONTINUE
         SUMV = (SUMVX + SUMVY + SUMVZ)/3.D0

         WRITE(6,'("  Total Momentum is                ",E10.4)') SUMV
         WRITE(6,'("  Temperature(kT/m) is             ",F20.16)')TEMP
         WRITE(6,'("  Final E/NkT is                   ",F15.8/)')ENKT

         WRITE(6,'("  Final Time(Reduced) is           ",E20.13)') TN
         WRITE(6,'("  Collision Rate is 	       ",F15.3)') RATE
         WRITE(6,'("  Reduced Mean Collision Time is   ",E20.13)') TBC
         WRITE(6,'("  True Mean Collision Time is      ",E20.13/)') TRUTBC

         WRITE(6,'("  Number of Inter Collisions       ",I10)') NINTER
         WRITE(6,'("  Number of Intra Collisions       ",I10)') NINTRA
         WRITE(6,'("  Number of Neighbor Collisions    ",I10)') NNEIGH
         WRITE(6,'("  Number of Bond Stretches         ",I10)') NBOND
         WRITE(6,'("  Sum of the above collisions is   ",I10/)') NTOTAL

         WRITE(6,'("  Compressibility from summing is  ",F15.8/)') ZTOT

         EXTIME=DELTA/3600
         WRITE(6,'("  Execution Time (CPU hours)       ",F15.4)') EXTIME
         COLRAT=FLOAT(COLL-1)/EXTIME/1.0E6
         WRITE(6,'("  Millions of Coll/CPU hour        ",F15.3/)') COLRAT

10000 CONTINUE
      
      STOP
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C**   Subroutine to calculate collision dynamics for particles
C**   in contact or at the bond extension. Also computes collision virial.

      SUBROUTINE BUMP(I, J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER PARTNR, COLTYP
      COMMON /TIMES / COLTIM(NOP), PARTNR(NOP), COLTYP(NOP)
#include "header.f"
      PARAMETER(TSTEP=1D-10)

      VXIJ=VX(I) - VX(J)
      VYIJ=VY(I) - VY(J)
      VZIJ=VZ(I) - VZ(J)

      RXIJ = RX(I) - RX(J)
      RYIJ = RY(I) - RY(J)
      RZIJ = RZ(I) - RZ(J)

      RXIJ = RXIJ + VXIJ*TNEXT
      RYIJ = RYIJ + VYIJ*TNEXT
      RZIJ = RZIJ + VZIJ*TNEXT

      RXIJ = RXIJ - DNINT(RXIJ)
      RYIJ = RYIJ - DNINT(RYIJ)
      RZIJ = RZIJ - DNINT(RZIJ)

      RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
      BIJ = RXIJ*VXIJ + RYIJ*VYIJ + RZIJ*VZIJ
C**   Assume a core or bond collision
      FACTOR = BIJ/RIJSQ
#ifdef SQWELL
      IF(COLTYP(I) .GT. 1)THEN 
         BIJSQ = BIJ*BIJ
         IF(COLTYP(I) .EQ. 2)THEN 
C     ** Capture
            FACTOR=0.5D0*(DSQRT(DELPE + BIJSQ) + BIJ)/SIG2SQ
         ELSEIF(COLTYP(I) .EQ. 3)THEN 
            IF(BIJSQ .GT. DELPE)THEN
C     ** Dissociation
               FACTOR=0.5D0*(BIJ - DSQRT(BIJSQ - DELPE))/SIG2SQ
            ELSE
C     ** Bounce
               FACTOR=BIJ/SIG2SQ
            END IF
         END IF
      END IF
C**   Since this is a well event, add a little extra to get off the well
#endif

      DELVX = FACTOR * RXIJ
      DELVY = FACTOR * RYIJ
      DELVZ = FACTOR * RZIJ

      VX(I) = VX(I) - DELVX
      VX(J) = VX(J) + DELVX
      VY(I) = VY(I) - DELVY
      VY(J) = VY(J) + DELVY
      VZ(I) = VZ(I) - DELVZ
      VZ(J) = VZ(J) + DELVZ

#ifdef SQWELL
      IF(COLTYP(I) .GT. 1)THEN 
C** For a well event, move I ahead slightly to get off well
         RX(I) = RX(I) + VX(I)*TSTEP
         RY(I) = RY(I) + VY(I)*TSTEP
         RZ(I) = RZ(I) + VZ(I)*TSTEP
         RX(J) = RX(J) + VX(J)*TSTEP
         RY(J) = RY(J) + VY(J)*TSTEP
         RZ(J) = RZ(J) + VZ(J)*TSTEP
      END IF
#endif
C**   Take care of false positions by rewinding
      RX(I)=RX(I)+DELVX*TNEXT
      RY(I)=RY(I)+DELVY*TNEXT
      RZ(I)=RZ(I)+DELVZ*TNEXT
      RX(J)=RX(J)-DELVX*TNEXT
      RY(J)=RY(J)-DELVY*TNEXT
      RZ(J)=RZ(J)-DELVZ*TNEXT

      W = DELVX*RXIJ + DELVY*RYIJ + DELVZ*RZIJ
      
      RETURN
      END


C***  Subroutine to calculate bond stretch times for I and J

      SUBROUTINE BTIME(K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER PARTNR, COLTYP
      COMMON /TIMES / COLTIM(NOP), PARTNR(NOP), COLTYP(NOP)
      COMMON /BONDS / ISBND(0:NOP)
#include "header.f"

      IF(ISBND(K) .EQ. 1) THEN
	 L=K+1
C***  No need to apply PBC since chains are unfolded
	 VXIJ = VX(K) - VX(L)
	 VYIJ = VY(K) - VY(L)
	 VZIJ = VZ(K) - VZ(L)

	 RXIJ = RX(K) - RX(L) + VXIJ*TNEXT
	 RYIJ = RY(K) - RY(L) + VYIJ*TNEXT
	 RZIJ = RZ(K) - RZ(L) + VZIJ*TNEXT

	 BIJ = RXIJ*VXIJ + RYIJ*VYIJ + RZIJ*VZIJ
	 RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
	 VIJSQ = VXIJ*VXIJ + VYIJ*VYIJ + VZIJ*VZIJ
	 C = RIJSQ - BNDLSQ
	 DISCR = BIJ*BIJ - VIJSQ*C

C***  Using the formula from Numerical Recipes 5.5
	 IF(BIJ .GT. 0.D0)THEN
            TIJ = -C/(DSQRT(DISCR) + BIJ)
	 ELSE
C***  This doesn't happen very often, but it does happen
            TIJ = (DSQRT(DISCR) - BIJ)/VIJSQ
	 END IF

	 TIJ = TIJ + TNEXT

	 IF(TIJ .LT. COLTIM(K)) THEN
            COLTIM(K) = TIJ
            PARTNR(K) = L
            COLTYP(K) = 0
	 END IF

      END IF

      RETURN
      END

C**   Subroutine to look for collisions of atoms J greater than I

      SUBROUTINE UPLIST(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER PARTNR, COLTYP
      COMMON /TIMES / COLTIM(NOP), PARTNR(NOP), COLTYP(NOP)

      INTEGER     MAXARR
      PARAMETER ( MAXARR = MAXNAB*NOP )
      INTEGER*2  UPNAB(MAXARR), DNNAB(MAXARR)
      INTEGER*2  NNABUP(NOP), NNABDN(NOP)
      INTEGER  INDEX(NOP)
      COMMON /NABORS / UPNAB, NNABUP, DNNAB, NNABDN, INDEX

#include "header.f"

      PARAMETER(TIMBIG = 1D10)

      IF(I .EQ. NOP) RETURN
      KBEG = INDEX(I) + 1
      KEND = KBEG + NNABUP(I) - 1

C**   Reset COLTIM for the collision
      COLTIM(I) = TIMBIG

      DO 100 K = KBEG,KEND
         J=UPNAB(K)
	 VXIJ = VX(I) - VX(J)
	 VYIJ = VY(I) - VY(J)
	 VZIJ = VZ(I) - VZ(J)
	 RXIJ = RX(I) - RX(J)
	 RYIJ = RY(I) - RY(J)
	 RZIJ = RZ(I) - RZ(J)

	 RXIJ = RXIJ + VXIJ*TNEXT
	 RYIJ = RYIJ + VYIJ*TNEXT
	 RZIJ = RZIJ + VZIJ*TNEXT

	 RXIJ = RXIJ - DNINT(RXIJ)
	 RYIJ = RYIJ - DNINT(RYIJ)
	 RZIJ = RZIJ - DNINT(RZIJ)

	 BIJ   = RXIJ*VXIJ + RYIJ*VYIJ + RZIJ*VZIJ
	 RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
	 VIJSQ = VXIJ*VXIJ + VYIJ*VYIJ + VZIJ*VZIJ
	 C     = RIJSQ - SIGSQ
#ifdef BELMAN
	 IF((J-I) .EQ. 1 .AND. ICHN(I) .EQ. ICHN(J) )C = RIJSQ - SIGSQIJ
#endif
	 DISCR = BIJ*BIJ - VIJSQ*C

#ifdef SQWELL
	 CIJ2  = RIJSQ - SIG2SQ
	 DISCR2 = BIJ*BIJ - VIJSQ*CIJ2
C** Calculate the square well time first 
	 IF(CIJ2 .LT. 0.D0)THEN
C** Rij is within the attractive well -- Bounce or dissociation
            TWELL = (DSQRT(DISCR2) - BIJ)/VIJSQ + TNEXT
            IF(TWELL .LT. COLTIM(I)) THEN
               COLTIM(I) = TWELL
               PARTNR(I) = J
               COLTYP(I) = 3
            END IF
            IF(BIJ .LT. 0.D0 .AND. DISCR .GT. 0.D0)THEN
C**   Using the formula from Numerical Recipes 5.5
               TIJ = C/(DSQRT(DISCR) - BIJ) + TNEXT
               IF(TIJ .LT. COLTIM(I)) THEN
                  COLTIM(I) = TIJ
                  PARTNR(I) = J
                  COLTYP(I) = 1
               END IF
            END IF
	 ELSE
C** Particles outside the well, Capture
            IF(BIJ .LT. 0.D0 .AND. DISCR2 .GT. 0.D0)THEN
C** Centers are approaching
               TWELL = (-BIJ - DSQRT(DISCR2))/VIJSQ + TNEXT
               IF(TWELL .LT. COLTIM(I)) THEN
                  COLTIM(I) = TWELL
                  PARTNR(I) = J
                  COLTYP(I) = 2
               END IF
            END IF
	 END IF
#else 
	 IF(BIJ .LT. 0.D0 .AND. DISCR .GT. 0.D0)THEN
C**   Using the formula from Numerical Recipes 5.5
            TIJ = C/(DSQRT(DISCR) - BIJ) + TNEXT
            IF(TIJ .LT. COLTIM(I)) THEN
               COLTIM(I) = TIJ
               PARTNR(I) = J
               COLTYP(I) = 1
            END IF
         END IF
#endif

 100  CONTINUE

      RETURN
      END

C**   Subroutine to look for collisions of atoms J less than I

      SUBROUTINE DNLIST(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER PARTNR, COLTYP
      COMMON /TIMES / COLTIM(NOP), PARTNR(NOP), COLTYP(NOP)

      INTEGER     MAXARR
      PARAMETER ( MAXARR = MAXNAB*NOP )
      INTEGER*2  UPNAB(MAXARR), DNNAB(MAXARR)
      INTEGER*2  NNABUP(NOP), NNABDN(NOP)
      INTEGER  INDEX(NOP)
      COMMON /NABORS / UPNAB, NNABUP, DNNAB, NNABDN, INDEX

#include "header.f"

      IF(I .EQ. 1) RETURN
      KBEG = INDEX(I) + 1
      KEND = KBEG + NNABDN(I) - 1

CDIR$ IVDEP
      DO 100 K=KBEG, KEND
	 J = DNNAB(K)
	 VXIJ = VX(J) - VX(I)
	 VYIJ = VY(J) - VY(I)
	 VZIJ = VZ(J) - VZ(I)
	 RXIJ = RX(J) - RX(I)
	 RYIJ = RY(J) - RY(I)
	 RZIJ = RZ(J) - RZ(I)

	 RXIJ = RXIJ + VXIJ*TNEXT
	 RYIJ = RYIJ + VYIJ*TNEXT
	 RZIJ = RZIJ + VZIJ*TNEXT

	 RXIJ = RXIJ - DNINT(RXIJ)
	 RYIJ = RYIJ - DNINT(RYIJ)
	 RZIJ = RZIJ - DNINT(RZIJ)

	 BIJ   = RXIJ*VXIJ + RYIJ*VYIJ + RZIJ*VZIJ
	 RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
	 VIJSQ = VXIJ*VXIJ + VYIJ*VYIJ + VZIJ*VZIJ
	 C     = RIJSQ - SIGSQ
#ifdef BELMAN
	 IF((I-J) .EQ. 1 .AND. ICHN(I) .EQ. ICHN(J))C = RIJSQ - SIGSQIJ
#endif
	 DISCR = BIJ*BIJ - VIJSQ*C
#ifdef SQWELL
	 CIJ2  = RIJSQ - SIG2SQ
	 DISCR2 = BIJ*BIJ - VIJSQ*CIJ2

	 IF(CIJ2 .LT. 0.D0)THEN
C** Rij is within the attractive well -- Bounce or dissociation
            TWELL = (DSQRT(DISCR2) - BIJ)/VIJSQ + TNEXT
            IF(TWELL .LT. COLTIM(J)) THEN
               COLTIM(J) = TWELL
               PARTNR(J) = I
               COLTYP(J) = 3
               CALL ADDTRE(J)
            END IF
            IF(BIJ .LT. 0.D0 .AND. DISCR .GT. 0.D0)THEN
C**   Using the formula from Numerical Recipes 5.5
               TIJ = C/(DSQRT(DISCR) - BIJ) + TNEXT
               IF(TIJ .LT. COLTIM(J)) THEN
                  COLTIM(J) = TIJ
                  PARTNR(J) = I
                  COLTYP(J) = 1
                  CALL ADDTRE(J)
               END IF
            END IF
	 ELSE
C** Particles outside the well, Capture
            IF(BIJ .LT. 0.D0 .AND. DISCR2 .GT. 0.D0)THEN
               TWELL = (-BIJ - DSQRT(DISCR2))/VIJSQ  + TNEXT
               IF(TWELL .LT. COLTIM(J)) THEN
                  COLTIM(J) = TWELL
                  PARTNR(J) = I
                  COLTYP(J) = 2
                  CALL ADDTRE(J)
               END IF
            END IF
	 END IF
#else 
	 IF(BIJ .LT. 0.D0 .AND. DISCR .GT. 0.D0)THEN
C**   Using the formula from Numerical Recipes 5.5
	    TIJ = C/(DSQRT(DISCR) - BIJ) + TNEXT
            IF(TIJ .LT. COLTIM(J)) THEN
               COLTIM(J) = TIJ
               PARTNR(J) = I
               COLTYP(J) = 1
               CALL ADDTRE(J)
            END IF
         END IF
#endif
 100  CONTINUE
      RETURN
      END

C**   Subroutine to read configuration file in the following format
C**   x,y,z
C**   ....
C**   v_x, v_y, v_z
C**   ...
      SUBROUTINE READCF (CF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PHSSPC / RX(NOP),RY(NOP),RZ(NOP),
     &     VX(NOP),VY(NOP),VZ(NOP)
      CHARACTER CF*16
      OPEN(UNIT=90,FILE=CF//'.in',STATUS='OLD',FORM='FORMATTED')	

      DO 100 I=1,NOP
         READ(90,*) RX(I), RY(I), RZ(I)
 100  CONTINUE
      DO 200 I=1,NOP
         READ(90,*) VX(I), VY(I), VZ(I)
 200  CONTINUE

      CLOSE(UNIT=90)

      RETURN
      END

C**   Subroutine to write configuration file.
C**   Periodic boundaries are applied before writing.

      SUBROUTINE WRITCF (NFILE, BOXL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PHSSPC / RX(NOP),RY(NOP),RZ(NOP),
     &     VX(NOP),VY(NOP),VZ(NOP)
      BOXLM1 = 1.D0/BOXL

      DO 100 I=1,NOP
	 RXI = RX(I) - BOXL*DNINT(RX(I)*BOXLM1)
	 RYI = RY(I) - BOXL*DNINT(RY(I)*BOXLM1)
	 RZI = RZ(I) - BOXL*DNINT(RZ(I)*BOXLM1)
         RXI=RXI/BOXL
         RYI=RYI/BOXL
         RZI=RZI/BOXL
         WRITE(NFILE,*)RXI, RYI, RZI
 100  CONTINUE
      DO 200 I=1,NOP
         WRITE(NFILE,*)VX(I), VY(I), VZ(I)
 200  CONTINUE

      RETURN
      END

C**   Subroutine to check and adjust total momentum

      SUBROUTINE MOMENT(SETEMP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PHSSPC / RX(NOP),RY(NOP),RZ(NOP),
     &     VX(NOP),VY(NOP),VZ(NOP)

      SUMVX=0.D0
      SUMVY=0.D0
      SUMVZ=0.D0
      SUMV=0.D0
      E=0.D0

      DO 100 I=1,NOP
	 SUMVX=SUMVX + VX(I)
	 SUMVY=SUMVY + VY(I)
	 SUMVZ=SUMVZ + VZ(I)
 100  CONTINUE

      SCALEX=SUMVX/DBLE(NOP)
      SCALEY=SUMVY/DBLE(NOP)
      SCALEZ=SUMVZ/DBLE(NOP)

      DO 200 I=1,NOP
         VX(I)=VX(I) - SCALEX
         VY(I)=VY(I) - SCALEY
         VZ(I)=VZ(I) - SCALEZ
 200  CONTINUE

      DO 300 I=1,NOP
         E = E + VX(I)*VX(I) + VY(I)*VY(I) + VZ(I)*VZ(I)
 300  CONTINUE

      TEMP = E/3.D0/DBLE(NOP)
      SCALET=DSQRT(SETEMP/TEMP)

      DO 400 I=1,NOP
         VX(I)=VX(I)*SCALET
         VY(I)=VY(I)*SCALET
         VZ(I)=VZ(I)*SCALET
 400  CONTINUE
      
      RETURN
      END

C**   Subroutine to unfold chains and eliminate periodic boundary effect

      SUBROUTINE DELPBC(CHNLEN,BOXL,BOXLM1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PHSSPC / RX(NOP),RY(NOP),RZ(NOP),
     &     VX(NOP),VY(NOP),VZ(NOP)
      INTEGER CHNLEN
      
      DO 300 J=1,NOP,CHNLEN
         DO 200 I=J+1,J+CHNLEN-1
            RX(I) = RX(I) - BOXL*DNINT((RX(I)-RX(I-1))*BOXLM1)
            RY(I) = RY(I) - BOXL*DNINT((RY(I)-RY(I-1))*BOXLM1)
            RZ(I) = RZ(I) - BOXL*DNINT((RZ(I)-RZ(I-1))*BOXLM1)
 200     CONTINUE
 300  CONTINUE
      
      RETURN
      END

C**   Subroutine to assign chain numbers to each atom

      SUBROUTINE CHNNUM(CHNLEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CHNUM / ICHN(NOP)
      INTEGER CHNLEN

      DO I = 1,NOP
         ICHN(I) = (I-1)/CHNLEN
      END DO
      
      RETURN
      END


C**   Subroutine to calculate average bond length

      SUBROUTINE BNDLNGTH (LSTAR,CHNLEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PHSSPC / RX(NOP),RY(NOP),RZ(NOP),
     &     VX(NOP),VY(NOP),VZ(NOP)
      DOUBLE PRECISION LSTAR
      INTEGER NBNDS,CHNLEN

      BNDDIS=0.D0
      NBNDS=NOP/CHNLEN*(CHNLEN-1)

      DO 110 K=1,NOP,CHNLEN
	 DO 100 I=K,K+CHNLEN-2,1
            J=I+1
            RXIJ = RX(I) - RX(J)
            RYIJ = RY(I) - RY(J)
            RZIJ = RZ(I) - RZ(J)
            RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
            BNDDIS=BNDDIS+DSQRT(RIJSQ)
 100     CONTINUE
 110  CONTINUE

      LSTAR=BNDDIS/DBLE(NBNDS)

      RETURN
      END

C**   Subroutine to test for overlaps and calculate kinetic energy

      SUBROUTINE CHECK(KE,OVRLAP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "header.f"
      DOUBLE PRECISION KE
      LOGICAL OVRLAP

      PARAMETER (TOL = 1.0D-6)
      OVRLAP=.FALSE.
      KE=0.D0
      PE=0.D0

      DO 100 I=1,NOP-1
         DO 99 J =I+1,NOP
            RXIJ=RX(I)-RX(J)
            RYIJ=RY(I)-RY(J)
            RZIJ=RZ(I)-RZ(J)
            RXIJ = RXIJ - DNINT(RXIJ)
            RYIJ = RYIJ - DNINT(RYIJ)
            RZIJ = RZIJ - DNINT(RZIJ)
            RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
            RIJMIN = SIGSQ
#ifdef SQWELL
            IF(RIJSQ .LT. SIG2SQ) PE = PE - V0
#endif

C**   For Bellmans type bond, check that I,J are on the same 
C**   chain, and are neighbors.
            IF((J-I) .EQ. 1 .AND. ICHN(I) .EQ. ICHN(J))RIJMIN = SIGSQIJ
            IF(RIJSQ .LT. RIJMIN) THEN
               RIJ=DSQRT(RIJSQ/RIJMIN)
               IF( (1.D0-RIJ) .GT. TOL)THEN
                  OVRLAP = .TRUE.
                  WRITE(6,'(" Overlap I,J,RIJ = ",2I7,F15.10)')I,J,RIJ
               END IF
            END IF
 99      CONTINUE
 100  CONTINUE

C**   Calculate the Kinetic Energy
      DO 200 I=1,NOP
         KE = KE + VX(I)*VX(I) + VY(I)*VY(I) + VZ(I)*VZ(I)
 200  CONTINUE
      KE=0.5D0*KE

      RETURN
      END

C************************************************************************
C**   Decide whether the neighbor list needs to be renewed

C**   Variables
c**   rx0,ry0,rz0	        coordinates at last update
c**   rlist			Radius of neighbor list
c**   dispmx			largest displacement allowed

      SUBROUTINE NABCHK()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "header.f"

      DIMENSION RX0(NOP), RY0(NOP), RZ0(NOP), RSQ(NOP)
      LOGICAL   EXCNAB

      SAVE RX0,RY0,RZ0,EXCNAB

C**   Calculate maximum displacement since last update

      DO 30 I = 1, NOP
         DX =  RX(I) - RX0(I)
         DY =  RY(I) - RY0(I)
         DZ =  RZ(I) - RZ0(I)
         RSQ(I) = DX*DX + DY*DY + DZ*DZ
 30   CONTINUE

      IMAX = ISMAX(NOP, RSQ, 1)
      DISPMX = RSQ(IMAX)

      IF(DISPMX .GT. DSPLIM .OR. COLL .EQ. 0) THEN
C**   Exceeded displacement so rehash neighbor lists

C**   Save the old positions
         DO 100 I = 1, NOP
            RX0(I) = RX(I)
            RY0(I) = RY(I)
            RZ0(I) = RZ(I)
 100     CONTINUE

C**   Use a link list to save time
         CALL CREATE_LINKS
 99      CALL LNKLST(EXCNAB)

         IF(EXCNAB)THEN
C**   Too many neighbors so reduce cutoff
            RLIST=RLIST*0.95D0
            RLSTSQ = RLIST*RLIST
#ifdef SQWELL
            SIGCUT = SIG2
#else 
            SIGCUT = SIGMA
#endif
C**   Use a conservative value for the max displacement
C**   Do this to avoid the likelihood of missing an event
C**   Since we don't recheck the displacement each time step
            DSPLIM = ((RLIST - SIGCUT)/2.5D0)**2.D0
      WRITE(6,'(" Neighbor cutoff reduced to ",F10.4)')RLIST/SIGCUT
            EXCNAB = .FALSE.
            GO TO 99
         END IF

C**   Check the new neighbors for possible events
         CALL DETREE
         DO 400 L=1,NOP
            CALL UPLIST(L)
            CALL BTIME(L)
            CALL SETREE(L)
 400     CONTINUE

      END IF

      RETURN
      END
      
C**   FUNCTION to find maximum array position
      
      FUNCTION ISMAX(N,X,INC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION X(NOP)

      DUMMY1 = 0.D0
      ISMAX  = -1

      DO 100 K=1,N,INC
         IF(X(K) .GT. DUMMY1)THEN
            DUMMY1=X(K)
            ISMAX=K
         END IF
 100  CONTINUE

      RETURN
      END

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C** Construction of cell linked-lists

C** References:
C** Quentrec and Brot, J. Comput. Phys. 13, 430, 1975.
C** Hockney and Eastwood, Comp. Sim. using Particles,McGraw Hill, 1981.

C**                
C** Subroutine MAPS
C**    Sets up map of cell structure
C** Subroutine CREATE_LINKS
C**    Sets up head of chain array and linked list
C** Subroutine LNKLST
C**    Calculates neighbors using a linked list

C** Subroutine maps is called once at the start of a simulation
C** to establish cell neighbour identities.  At each renewal,
C** subroutine create_links is called to set up the linked list
C** which is immediately used by subroutine LNKLST.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C** Routine to set up a map of neighboring cells

C** Variables:
C** INTEGER M                  Number of cells in each direction
C** INTEGER MAPSIZ             Size of cell-cell map 
C** INTEGER MAP(MAPSIZ)        List of neighboring cells

C** This subroutine sets up a list of the thirteen neighboring
C** cells for each small cell in the central box. The
C** effects of the periodic boundary conditions are included.
C** The subroutine is called once at the beginning of the
C** simulation and the map is used in the LNKLST subroutine.

        SUBROUTINE MAPS
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER     IX, IY, IZ, IMAP, ICELL
	INTEGER     MAXARR 
	PARAMETER ( MAXARR = MAXNAB*NOP )
        INTEGER     M, NCELL, MAPSIZ
#ifdef SMARTF77
C** Some FORTRAN compilers may have trouble with this
      PARAMETER( M = NINT(NOP/Cell_Occ)**(0.3333333333) )
#else
      PARAMETER( M = Number_of_cells )
#endif
        PARAMETER ( NCELL = M * M * M )
        PARAMETER ( MAPSIZ = 13 * NCELL )
        COMMON /LNKCOM / LINKS(NOP+NCELL), MAP(MAPSIZ)

C** Statement function to give cell index
        ICELL ( IX, IY, IZ) = 1 + MOD ( IX - 1 + M, M )
     :                          + MOD ( IY - 1 + M, M ) * M
     :                          + MOD ( IZ - 1 + M, M ) * M * M

C**  Find half the nearest neighbours of each cell
        DO 50 IZ = 1, M
           DO 40 IY = 1, M
              DO 30 IX = 1, M

                 IMAP = ( ICELL ( IX, IY, IZ ) - 1 ) * 13

                 MAP( IMAP + 1  ) = ICELL( IX + 1, IY    , IZ     )
                 MAP( IMAP + 2  ) = ICELL( IX + 1, IY + 1, IZ     )
                 MAP( IMAP + 3  ) = ICELL( IX    , IY + 1, IZ     )
                 MAP( IMAP + 4  ) = ICELL( IX - 1, IY + 1, IZ     )
                 MAP( IMAP + 5  ) = ICELL( IX + 1, IY    , IZ - 1 )
                 MAP( IMAP + 6  ) = ICELL( IX + 1, IY + 1, IZ - 1 )
                 MAP( IMAP + 7  ) = ICELL( IX    , IY + 1, IZ - 1 )
                 MAP( IMAP + 8  ) = ICELL( IX - 1, IY + 1, IZ - 1 )
                 MAP( IMAP + 9  ) = ICELL( IX + 1, IY    , IZ + 1 )
                 MAP( IMAP + 10 ) = ICELL( IX + 1, IY + 1, IZ + 1 )
                 MAP( IMAP + 11 ) = ICELL( IX    , IY + 1, IZ + 1 )
                 MAP( IMAP + 12 ) = ICELL( IX - 1, IY + 1, IZ + 1 )
                 MAP( IMAP + 13 ) = ICELL( IX    , IY    , IZ + 1 )

30            CONTINUE
40         CONTINUE
50      CONTINUE

        RETURN
        END

C**   Subroutine to set up linked list and the head pointers

C**   INTEGER M                  Number of cells in each direction
C**   INTEGER NCELL              Total number of cells (M**3)
C**   INTEGER LINKS(NOP+NCELL)   Linked list of atoms

C**   Each atom is sorted into one of the M**3 small cells.
C**   The first atom in each cell is placed in the head position.
C**   The head position consists of the upper NCELL positions in the
C**   links array. Subsequent atoms are placed in the linked list array.
C**   All coordinates are assumed to be between -0.5 and +0.5.

      SUBROUTINE CREATE_LINKS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER     MAXARR 
      PARAMETER ( MAXARR = MAXNAB*NOP )
      INTEGER     M, NCELL, MAPSIZ
#ifdef SMARTF77
C** Some FORTRAN compilers may have trouble with this
      PARAMETER( M = NINT(NOP/Cell_Occ)**(0.3333333333) )
#else
      PARAMETER( M = Number_of_cells )
#endif
      PARAMETER ( NCELL = M * M * M )
      PARAMETER ( MAPSIZ = 13 * NCELL )
      COMMON /LNKCOM/ LINKS(NOP+NCELL), MAP(MAPSIZ)

#include "header.f"

C**   Zero head of link list array
      DO 10 I = NOP + 1, NOP + NCELL
         LINKS(I) = 0
 10   CONTINUE

C**   Sort all atoms
      DO 20 I = 1, NOP
         RXI = RX(I) - DNINT(RX(I)) + BOXHLF
         RYI = RY(I) - DNINT(RY(I)) + BOXHLF
         RZI = RZ(I) - DNINT(RZ(I)) + BOXHLF
         ICELL = 1 + INT(RXI*CELLI) + (INT(RYI*CELLI) 
     &        + INT(RZI*CELLI)*M)*M
         LINKS(I)    = LINKS(ICELL+NOP)
         LINKS(ICELL+NOP) = I
 20   CONTINUE

      RETURN
      END

C**   Subroutine to find nearest neighbors using a link list

C**   Variables:
C**   INTEGER M                  Number of cells in each direction
C**   INTEGER NCELL              Number of small cells (M**3)
C**   INTEGER MAPSIZ             Size of cell-cell map
C**   INTEGER LINKS(NOP+NCELL)   The linked list
C**   INTEGER MAP(MAPSIZ)        List of neighboring cells

      SUBROUTINE  LNKLST(EXCNAB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "header.f"
      INTEGER     MAXARR 
      PARAMETER ( MAXARR = MAXNAB*NOP )
      INTEGER     M, NCELL, MAPSIZ
#ifdef SMARTF77
C** Some FORTRAN compilers may have trouble with this
      PARAMETER( M = NINT(NOP/Cell_Occ)**(0.3333333333) )
#else
      PARAMETER( M = Number_of_cells )
#endif
      PARAMETER ( NCELL = M * M * M )
      PARAMETER ( MAPSIZ = 13 * NCELL )
      COMMON /LNKCOM / LINKS(NOP+NCELL), MAP(MAPSIZ)

      INTEGER*2  UPNAB(MAXARR), DNNAB(MAXARR) 
      INTEGER*2  NNABUP(NOP), NNABDN(NOP)
      INTEGER  INDEX(NOP)
      COMMON /NABORS / UPNAB, NNABUP, DNNAB, NNABDN, INDEX

      DIMENSION   JNEAR(MAXNAB*Cell_Occ)
      INTEGER     ICELL, JCELL0, JCELL, I, J, NABOR
      LOGICAL     EXCNAB

C**   Zero the neighbor index arrays
      DO I = 1,NOP
         NNABDN(I)=0
         NNABUP(I)=0
      END DO
      
C**   Loop Over All Cells **
      DO 5000 ICELL = 1, NCELL
         I = LINKS(ICELL + NOP)
         
C**   Loop Over All Molecules In The Cell, Loop Over I **
         DO WHILE( I .GT. 0 )
            
C**   Loop over all molecules below I in the current cell
C**   Use DO WHILE to determine number of iterations and
C**   to gather the LINKS array
            J = LINKS(I)
            NITER = 0
            DO WHILE( J .GT. 0) 
               NITER = NITER + 1
               JNEAR(NITER) = J
               J = LINKS(J)
            END DO
            
C**   Loop over neighboring cells
            JCELL0 = 13 * (ICELL - 1)
            
            DO 4000 NABOR = 1, 13
               JCELL = MAP ( JCELL0 + NABOR )
               
C**   Loop over all molecules in neighboring cells 
               J = LINKS(JCELL + NOP)
               DO WHILE( J .NE. 0 )
                  NITER = NITER + 1
                  JNEAR(NITER) = J
                  J = LINKS(J)
               END DO
 4000       CONTINUE
            
C**   This loop vectorizes
C     DIR$ IVDEP
            DO NN = 1,NITER
               J = JNEAR(NN)
               RXIJ  = RX(I) - RX(J)
               RYIJ  = RY(I) - RY(J)
               RZIJ  = RZ(I) - RZ(J)
               RXIJ = RXIJ - DNINT(RXIJ)
               RYIJ = RYIJ - DNINT(RYIJ)
               RZIJ = RZIJ - DNINT(RZIJ)
               RSQ   = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
C**   J never equals I in this loop nor does J repeat itself
               IF (RSQ .LT. RLSTSQ)THEN
                  J = JNEAR(NN)
C**   J could be less than or greater than I
                  IF(J .LT. I) THEN
                     UPNAB(INDEX(J) + NNABUP(J) + 1) = I
                     DNNAB(INDEX(I) + NNABDN(I) + 1) = J
                     NNABDN(I)   = NNABDN(I) + 1
                     NNABUP(J)   = NNABUP(J) + 1
                  ELSE
                     UPNAB(INDEX(I) + NNABUP(I) + 1) = J
                     DNNAB(INDEX(J) + NNABDN(J) + 1) = I
                     NNABUP(I)   = NNABUP(I) + 1
                     NNABDN(J)   = NNABDN(J) + 1
                     IF(NNABDN(J) .GE. NABLIM) EXCNAB = .TRUE.
                  ENDIF
               ENDIF
            END DO

            IF(EXCNAB)RETURN 
            I = LINKS(I)
         END DO
         
 5000 CONTINUE
      
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C*** Binary tree algorithm
C*** See D. C. Rapaport, J. Comp. Phys. 34, 184, (1980)

C**   Subroutine to remove entire tree
      SUBROUTINE DETREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /EVLIST / NODEP(0:NOP),NODEL(0:NOP),NODER(0:NOP)

      DO I=0,NOP
         NODER(I) = 0
         NODEL(I) = 0
         NODEP(I) = 0
      END DO

      RETURN
      END


C**   Subroutine to search down left most branch and
C**   find the next event
      SUBROUTINE NEXTEV(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /EVLIST / NODEP(0:NOP),NODEL(0:NOP),NODER(0:NOP)

      I=NODER(0)
      DO WHILE(NODEL(I) .NE. 0)
         I = NODEL(I)
      END DO

      RETURN
      END

C**   Subroutine to delete old node and put new event in correc location
      SUBROUTINE ADDTRE(NUNODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /EVLIST / NODEP(0:NOP),NODEL(0:NOP),NODER(0:NOP)
      INTEGER PARTNR, COLTYP
      COMMON /TIMES / COLTIM(NOP),PARTNR(NOP),COLTYP(NOP)

C**   Delete a node from the tree and reconnect the pointers around
C**   where it was removed.
      ID = NUNODE
      IF(NODER(ID) .EQ. 0)THEN
         IS = NODEL(ID)
      ELSE
         IF(NODEL(ID) .EQ. 0)THEN
	    IS = NODER(ID)
         ELSE
	    IF(NODEL(NODER(ID)) .EQ. 0)THEN
               IS = NODER(ID)
	    ELSE
               IS = NODEL(NODER(ID))

               DO WHILE(NODEL(IS) .NE. 0)
                  IS = NODEL(IS)
               END DO

               NODEP(NODER(IS)) = NODEP(IS)
               NODEL(NODEP(IS)) = NODER(IS)
               NODEP(NODER(ID)) = IS
               NODER(IS) = NODER(ID)
	    END IF
            NODEP(NODEL(ID)) = IS
            NODEL(IS) = NODEL(ID)
	 END IF
      END IF

      IP = NODEP(ID)
      NODEP(IS) = IP

      IF(NODEL(IP) .EQ. ID)THEN
         NODEL(IP) = IS
      ELSE
         NODER(IP) = IS
      END IF

C**   Put the NUNODE into the binary tree event list.
C**   Note that TNEW is required but is not stored via this sub.
C**   TNEW is the time that the event will happen.
C**   Start by checking if NODER(0)=0. Then the tree must be empty 
C**   and we can just put this event in first.
      INPO = 0
      TNEW = COLTIM(NUNODE)
C**   Otherwise, we must search thru all the branches to find the proper
C**   branch-end for event insertion. The idea is that the branch to
C**   the right should always have a larger time.

      INPO = NODER(0)
      NFOUND = 0 
      DO WHILE(NFOUND .EQ. 0)
C**   if TNEW less than TEVENT(INPO) then go left.
         IF(TNEW .LE. COLTIM(INPO))THEN
            IF(NODEL(INPO) .NE. 0)THEN
               INPO = NODEL(INPO)
            ELSE
               NODEL(INPO) = NUNODE
               NFOUND = 1
            END IF
C**   if TNEW greater than TEVENT(INPO) then go right.
         ELSE
            IF(NODER(INPO) .NE. 0)THEN
               INPO=NODER(INPO)
            ELSE
               NODER(INPO) = NUNODE
               NFOUND = 1
            END IF
         END IF
      END DO

C**   That's the end of the search for the insertion position.
C**   INPO now equals the position preceding where the event was added.
C**   NODEL and NODER are set to zero since this new node is at the end of a
C**   branch.
      NODEL(NUNODE) = 0
      NODER(NUNODE) = 0
C**   NODEP is set to INPO because that is where this node should go if any
C**   future deletions in the tree cause INPO to be moved.
      NODEP(NUNODE) = INPO

      RETURN
      END

C**   Put the new time into the tree, this does not delete the node first
C**   since it is called only for a new tree

      SUBROUTINE SETREE(L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /EVLIST / NODEP(0:NOP),NODEL(0:NOP),NODER(0:NOP)
      INTEGER PARTNR, COLTYP
      COMMON /TIMES / COLTIM(NOP),PARTNR(NOP),COLTYP(NOP)

      TNEW = COLTIM(L)

      INPO = 0
      IF(NODER(0) .NE. 0)THEN
         INPO = NODER(0)
         NFOUND = 0
         DO WHILE(NFOUND .EQ. 0)
            IF(TNEW .LE. COLTIM(INPO))THEN
               IF(NODEL(INPO) .NE. 0)THEN
                  INPO = NODEL(INPO)
               ELSE
                  NODEL(INPO) = L
                  NFOUND = 1
               END IF
            ELSE
               IF(NODER(INPO) .NE. 0)THEN
                  INPO=NODER(INPO)
               ELSE
                  NODER(INPO) = L
                  NFOUND = 1
               END IF
            END IF
         END DO
      ELSE
         NODER(0) = L
      END IF

      NODEL(L) = 0
      NODER(L) = 0
      NODEP(L) = INPO

      RETURN
      END
C</PRE>

