<PRE>
<CODE>
C***  Discontinuous MD for hard or square-well chains
C***  Algorithm implements the multiple-event scheduling of Alder.
C***  A combination subcell/link list method is used.
C***  Events are scheduled with a doubly-linked binary tree.
C***  Particles only experience periodic boundaries
C***  when interactions are considered. (Infinite 3D checkerboard)
C***  Densities are assumed to be liquid-like so the
C***  minimum image convention holds.
C***  Periodic boundary condition assumes a box length of 1
C***  Box length = 1 and sigma is scaled accordingly.
C***  Box is assumed to be centered at the origin (-0.5,+0.5)

C***  The following defines will be expanded by your C preprocessor
C***  this usually happens automatically by naming the file with 
C***  a .F extension.

C***  C preprocessor defines
C***  Number of particles
#define NOP 6144
C***  To use or not to use Belleman's type bond
C#define BELMAN

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

      PROGRAM MULTI_DMD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

#include "header.f"

      INTEGER*2 ISBND 
      COMMON /BONDS / ISBND(0:NOP)

      INTEGER     M, NCELL, MAPSIZ
C**   Some Fortran compilers can't handle this
      PARAMETER( M = NINT(NOP**(1.D0/3.D0)) )
      PARAMETER ( NCELL = M * M * M )
      PARAMETER ( MAPSIZ = 26 * NCELL )
C**   This can be INTEGER*2 as long as NOP and NCELL are less than 32,768
      INTEGER*2   LINKS(NOP+NCELL), MAP(MAPSIZ), INCELL(NOP)
      COMMON / LNKCOM / LINKS, MAP, INCELL

C**   Common block for binary tree and its associated link list
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     +     NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +     JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +     COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP

      PARAMETER (TIMBIG=1D10)
      DOUBLE PRECISION   LSTAR, KE, PE
      LOGICAL   CANON, LOOP, OVRLAP
      INTEGER   ACCESS, IERRNO
      REAL      TARRAY(2), DELTA, EXTIME, COLRAT
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

C**   Read Packing Fraction
      READ(5,'(/F10.5)') ETA
C**   Read Chain Length
      READ(5,'(/I3)') CHNLEN
C**   Read Total Number of Collisions Desired
      READ(5,'(/F10.5)') RNCOLL
C**   Read Total Number of Batches
      READ(5,'(/I2)') NBATCH
C**   Read Initial Configuration Filename
C***  Any 15 characters you like, an alpha character 
C***    and apprpriate extensions will be added.
      READ(5,'(/A15)') CNFILE
      READ(5,'(/F10.5)') DELFAC
      READ(5,'(/I10)')IDUM
      READ(5,'(/F10.5)') SETEMP
      READ(5,'(/I10)') MODUPD

      NCOLL = DINT(RNCOLL*1000000)
      NCOLL=NCOLL/NBATCH
      CLOSE(UNIT=5)

#ifdef SQWELL
      OPEN(UNIT=95,FILE='sw.cnf',STATUS='UNKNOWN',FORM='FORMATTED')
      READ(95,'(/ L1)') CANON
      READ(95,'(/ F10.5)') GAMMA
      READ(95,'(/ F10.5)') V0
      READ(95,'(/ I10)') MODTMP
      READ(95,'(/ I10)') MODTNEW

      IF(MODTMP .EQ. 0) THEN
         MODTMP = NCOLL+1
      ELSE
         MODTMP = NCOLL/MODTMP
      END IF

      IF(CANON) THEN
	 MODTNEW = NCOLL/MODTNEW
      ELSE
         MODTNEW = NCOLL+1
      END IF
#endif

      ICFRST=ICHAR('A')-1
      LOOP=.TRUE.
C**   Find the first unused identifier
      DO WHILE(LOOP)
         ICFRST=ICFRST+1
         FRSTLET=CHAR(ICFRST)
         CHKFIL=CNFILE//FRSTLET//'.md'
         IERRNO=ACCESS(CHKFIL,'r')
         IF(IERRNO .NE. 0)LOOP=.FALSE.
      END DO


C***  Do several batches at once

      DO 10000 NBAT=ICFRST,ICFRST+NBATCH-1
         RUNID=CHAR(NBAT)
         NXRUN=CHAR(NBAT+1)
         CF=CNFILE//RUNID
         COLL=0
         NNEIGH=0
         NINTER=0
         NINTRA=0
         NBOND=0
         NCROSS=0
         NWELL=0

C***  Zero the Virial Accumulator and time
         PSTOT=0.D0
         T=0.D0
         TNEXT = 0.D0


C**   Read in initial configuration
         CFTRAJ=CNFILE//CHAR(NBAT -1)
         CALL READCF(CF)

C***  Put output in next run for multiple batches
         IF(NBATCH .EQ. 1) THEN
            OPEN(UNIT=4,FILE=CF//'.out',STATUS='UNKNOWN',FORM='FORMATTED')
         ELSE
            OPEN(UNIT=4,FILE=CNFILE//NXRUN//'.in',STATUS='UNKNOWN',
     &           FORM='FORMATTED')
         END IF

         OPEN(UNIT=6,FILE=CF//'.md',STATUS='UNKNOWN',FORM='FORMATTED')
C         OPEN(UNIT=99,FILE=CF//'.dbg',STATUS='UNKNOWN',FORM='FORMATTED')


C**   Set chain number for each bead
         CALL CHNNUM(CHNLEN)

C**   Initial Parameters
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
         VOL = BOXL*BOXL*BOXL

         CELLI = DBLE(M)
         BOXHLF = BOXL/2.D0

C**   Rescale velocities to insure zero total momentum and set temp
         CALL MOMENT(SETEMP)

#ifdef SQWELL
         TSTAR=SETEMP
         SIG2=GAMMA*SIGMA
         SIG2SQ=SIG2*SIG2
         DELPE=4.D0*V0*SIG2SQ
#endif

C**   Check for Overlap and Calculate Energy
         CALL DELPBC(CHNLEN,BOXL,BOXLM1)

C***  Set up link lists
         CALL MAPS
         CALL CREATE_LINKS

         CALL CHECK(SIGSQ,SIGSQIJ,KE,OVRLAP,PE)
         EN=KE/DBLE(NOP)
         TEMP = 2.D0*EN/3.D0

#ifdef SQWELL
         OPEN(UNIT=9,FILE=CF//'.ener',STATUS='UNKNOWN',FORM='FORMATTED')         
         CALL ADJSTT(KE)
         TEMP = 2.D0*EN/3.D0
         TSTACT=TEMP/V0
         WRITE(9,'("Collision   Kin. Ener.  Pot. Ener. Tot. Ener. Tstar")')
         WRITE(9,'(I10,4F10.3)')COLL,KE,PE,TOTE,TSTACT
#endif

         ENKT = EN/TEMP
         PSFACT=3.D0*TEMP*DBLE(NCHN)*SIGMA

         IF (OVRLAP) THEN
            Write(6,*)' Particle overlap in initial Configuration'
            STOP
         END IF

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
         WRITE(6,'("  Host machine was            ",A20)')HSTNAM
         WRITE(6,'("  Reduced Density is          ",F10.7)') DENSTY
         WRITE(6,'("  Sigma is                    ",E13.6)') SIGMA
         WRITE(6,'("  Box Length is               ",E13.6)') BOXL
         WRITE(6,'("  Temperature kT/m is         ",F20.16)')TEMP
         WRITE(6,'("  Initial E/NkT is            ",F10.7/)') ENKT

#ifdef BELMAN
         LAVER = 1.D0
         WRITE(6,'("  Compiled for BELLEMANS BOND !!!!")')
#endif
         WRITE(6,'("  Bond Length Factor is       ",F10.7)') DELFAC
         WRITE(6,'("  Packing Fraction is         ",F10.7)') ETA
         WRITE(6,'("  Run length (MCollisions) is ",I10)') NCOLL/1000000

         WRITE(6,'("  Configuration Filename is   ",A15)') CNFILE

#ifdef SQWELL
         WRITE(6,'(" Compiled for SQUARE WELL SIMULATION ")')
         WRITE(6,'("  Reduced Temp(kT/eps) is     ",F10.7)') TSTAR
         WRITE(6,'("  Well Depth epsilon/m is    ",F10.7)') V0
         WRITE(6,'("  Well Width gamma is        ",F10.7)') GAMMA
         WRITE(6,'("  Initial Total Energy is    ",F15.8/)') TOTE
#endif
         WRITE(6,'("  Big Updates every           ",I8/)') MODUPD


#ifdef SQWELL
         SIGCUT = SIG2
#else 
         SIGCUT = SIGMA
#endif
         CELL_LENGTH = BOXL/CELLI 
         IF(CELL_LENGTH .LT. SIGCUT) THEN
            WRITE(6,'("Subcell size too small",F10.4)')CELL_LENGTH/SIGCUT
            STOP
         END IF
         
C**   Set an integer to decide if the bead has a bond time
         DO I = 1,NOP
            ISBND(I) = 1
         END DO

         DO I = CHNLEN,NOP,CHNLEN
            ISBND(I) = 0
         END DO

         CALL INLIST
         
C**   Schedule cell crossing times 
         DO I = 1,NOP
            CALL CELCRO(TNEXT,I)
         END DO

C**   Schedule core collision times, two at a time
         DO I = 1,NOP,2
            J = I + 1
            CALL NEWTIM(I,J)
            CALL BTIME(I,J)
         END DO

         WRITE(6,'("  Using Rapaports Binary tree method ")')

         WRITE(6,'(/"  START OF DYNAMICS")')
         CALL FLUSH(6)

#ifdef CRAY
         DELTA=SECOND()
#else 
         DELTA=DTIME(TARRAY)
#endif

C=================Main Loop=========================================

         DO 1000 COLL = 1,NCOLL

            CALL NEXT(TIJ,ICOTYP,I,J)
            TNEXT=TIJ

            IF(ICOTYP .EQ. 0)THEN
C**   A bond stretch
               NBOND = NBOND + 1
            ELSEIF(ICOTYP .EQ. 1) THEN
               IF(ICHN(I) .NE. ICHN(J))THEN
                  NINTER=NINTER+1
               ELSEIF((J-I) .EQ. 1)THEN
                  NNEIGH=NNEIGH + 1
               ELSE
                  NINTRA=NINTRA + 1
               END IF
            ELSEIF(ICOTYP .GT. 90)THEN
C**   A cell crossing
               NCROSS = NCROSS + 1
            END IF
#ifdef SQWELL
            IF(ICOTYP .GE. 2)THEN
               NWELL = NWELL + 1
            END IF
#endif

C**   A safety net, just to be sure 
            IF(TNEXT .LE. 0.D0) THEN
               WRITE(6,*)'***NEGATIVE OR ZERO COLLISION TIME***'
               WRITE(6,*)'Coll#, Time,  I,  J,  Coltype'
               WRITE(6,*)COLL,TNEXT, I, J, ICOTYP
               WRITE(6,*)'Total time = ',T
               STOP
            END IF


C**   Find new core collision times for those particles affected 
C**   by the previous event
            IF(ICOTYP .LE. 1)THEN
C**   Compute collision dynamics and collisional virial
               CALL BUMP(I,J)	
               PSTOT = PSTOT - W
C**   Update cell crossing time since velocity has changed
C**   Must do this first since the cell-crossing time is the entry point
C**   into the binary tree 
               CALL CELCRO(TNEXT,I)
               CALL CELCRO(TNEXT,J)
C**   Update possible core collisions
               CALL NEWTIM(I,J)
C**   Update bond collisions
               CALL BTIME(I,J)
            ELSEIF (ICOTYP .GT. 90) THEN
C**   Cell Crossing, update possible new events
               CALL UPCROS(I,J,ICOTYP)
            END IF

            IF(MOD(COLL,MODUPD) .EQ. 0)THEN
C**   Update positions
               T = T + TNEXT
               CALL BIGUP()
            END IF

 1000    CONTINUE

C=================Main Loop Ends =====================================

#ifdef CRAY
         DELTA=SECOND()
#else 
         DELTA=DTIME(TARRAY)
#endif

 1300    WRITE(6,'(/"  END OF DYNAMICS")')
         WRITE(6,'(/"  Final Colliding Pair               ",2I7)') I, J

	IF(ICOTYP .GT. 90)THEN
C**   We need to move the particle off the subcell division
C**   in order for the configuration to be reusable.
	RX(I) = RX(I)+VX(I)*1.0D-12
	RY(I) = RY(I)+VY(I)*1.0D-12
	RZ(I) = RZ(I)+VZ(I)*1.0D-12
	END IF

C**   Check for Overlap
         CALL CHECK(SIGSQ,SIGSQIJ,KE,OVRLAP,PE)
         IF (OVRLAP) THEN
            WRITE(6,'(" Particle overlap in Final Configuration")')
         END IF

C**   Compute average bond distance and scale for SIGMA=1
         IF(CHNLEN .GT. 1) THEN
            CALL BNDLNGTH(LSTAR,CHNLEN)
            LSTAR=LSTAR/SIGMA
         END IF

C**   Write Final Configuration
         CALL WRITCF (4, BOXL)

C**   Write out interesting info.
         PVNKT = PSTOT/DBLE(NCHN)/3.D0/T/TEMP
         KE = KE/DBLE(NOP)
         TEMP = 2.D0*KE/3.D0
         ENKT = KE/TEMP
         TN = T*DSQRT(TEMP)/SIGMA
         RATE = DBLE(NCOLL-NCROSS)/TN
         TBC = DBLE(NOP)/RATE/2.D0
         TRUTBC = 1.D0/RATE
         NTOTAL = NBOND+NINTRA+NINTER+NNEIGH+NCROSS

#ifdef SQWELL
         TSTAR=TEMP/V0
         NTOTAL=NTOTAL+NWELL
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

C**   Check total momentum again
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
         WRITE(6,'("  Final Time(Reduced) is           ",E20.13)') TN
         WRITE(6,'("  Collision Rate is 	              ",F15.3)') RATE
         WRITE(6,'("  Mean Collision Time is           ",E20.13)') TBC
         WRITE(6,'("  True Mean Collision Time is      ",E20.13/)') TRUTBC
         WRITE(6,'("  Number of Inter Collisions       ",I10)') NINTER
         WRITE(6,'("  Number of Intra Collisions       ",I10)') NINTRA
         WRITE(6,'("  Number of Neighbor Collisions    ",I10)') NNEIGH
         WRITE(6,'("  Number of Cell Crossings         ",I10)') NCROSS
         WRITE(6,'("  Number of Bond Stretches         ",I10)') NBOND
         WRITE(6,'("  Sum of the above collisions is   ",I10)') NTOTAL
         WRITE(6,'("  Final E/NkT is                   ",F15.8/)') ENKT
         WRITE(6,'("  Compressibility from summing is  ",F15.8)') ZTOT
         WRITE(6,'("  Average Bond Distance            ",F15.8)') LSTAR

         EXTIME=DELTA/3600
         COLRAT=FLOAT(COLL-1-NCROSS)/EXTIME/1E6
         WRITE(6,'("  Execution Time (CPU hours)       ",F15.4)') EXTIME
         WRITE(6,'("  Millions of Coll/CPU hour        ",F15.3)') COLRAT

10000 CONTINUE
      
      STOP
      END


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SUBROUTINES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

      SUBROUTINE CHECK(SIGSQ,SIGSQIJ,KE,OVRLAP,PE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PHSSPC / RX(NOP),RY(NOP),RZ(NOP),
     &     VX(NOP),VY(NOP),VZ(NOP)
      COMMON /CHNUM / ICHN(NOP)
      DOUBLE PRECISION KE,PE,SIGSQ,SIGSQIJ
      LOGICAL OVRLAP

      PARAMETER (TOL = 1.0D-6)
      OVRLAP=.FALSE.
      KE=0.D0

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

      
C**   Subroutine to look for bond stretch times for I and J
      SUBROUTINE BTIME(I,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*2 ISBND
      COMMON /BONDS / ISBND(0:NOP)
#include "header.f"

      DO 10 K = I-1,I
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
C**   Using the formula from Numerical Recipes 5.5
            IF(BIJ .GT. 0.D0)THEN
               TIJ = -C/(DSQRT(DISCR) + BIJ)
            ELSE
C**   This doesn't happen very often, but it does happen
               TIJ = (DSQRT(DISCR) - BIJ)/VIJSQ
            END IF
            TIJ = TIJ + TNEXT
            CALL SCHED(TIJ,0,K,L)
         END IF
 10   CONTINUE

      DO 20 K = J-1,J
         IF(ISBND(K) .EQ. 1 .AND. K .NE. I) THEN
            L=K+1
C**   No need to apply PBC since chains are unfolded
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
C**   Using the formula from Numerical Recipes 5.5
            IF(BIJ .GT. 0.D0)THEN
               TIJ = -C/(DSQRT(DISCR) + BIJ)
            ELSE
C**   This doesn't happen very often, but it does happen
               TIJ = (DSQRT(DISCR) - BIJ)/VIJSQ
            END IF
            TIJ = TIJ + TNEXT
            CALL SCHED(TIJ,0,K,L)
         END IF
 20   CONTINUE

      RETURN
      END

C**   Subroutine to calculate collision dynamics for particles
C**   in contact or at the bond extension. Also computes collision virial.

      SUBROUTINE BUMP(I, J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "header.f"

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
C**   FACTOR is just Bij/Sigma^2 for interchain collision
      FACTOR = (RXIJ*VXIJ + RYIJ*VYIJ + RZIJ*VZIJ)/RIJSQ
      DELVX = FACTOR * RXIJ
      DELVY = FACTOR * RYIJ
      DELVZ = FACTOR * RZIJ

      VX(I) = VX(I) - DELVX
      VX(J) = VX(J) + DELVX
      VY(I) = VY(I) - DELVY
      VY(J) = VY(J) + DELVY
      VZ(I) = VZ(I) - DELVZ
      VZ(J) = VZ(J) + DELVZ

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


C** Subroutine to place particle in their true positions

      SUBROUTINE BIGUP()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "header.f"
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP
C** Do a big update to put all the particles in the right place
C** Move particles forward by time TNEXT
C** and Reduce Collision Times in list

      DO 100 I=1,NOP
         TEVENT(I) = TEVENT(I) - TNEXT
         RX(I) = RX(I) + VX(I) * TNEXT
         RY(I) = RY(I) + VY(I) * TNEXT
         RZ(I) = RZ(I) + VZ(I) * TNEXT

C**   Circle left on the i-site nodes finding events of this site#
C**   divide by uptime/2 because we will do this again when we pick up
C**   the site this is colliding with.
         ID=I
         DO WHILE(ILNODE(ID).NE.I)
            ID=ILNODE(ID)
            TEVENT(ID)=TEVENT(ID)-TNEXT/2.D0
         END DO
         
C**   Circle left on the j-site nodes finding events of this site#
         ID=I
         DO WHILE(JLNODE(ID).NE.I)
            ID=JLNODE(ID)
            TEVENT(ID)=TEVENT(ID)-TNEXT/2.D0
         END DO

 100  CONTINUE

C** Rezero time
      TNEXT = 0.D0

      RETURN
      END

C**   Subroutine to handle a cell crossong event.

      SUBROUTINE CELCRO(TIJNEX,I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
#include "header.f"

      INTEGER     M, NCELL, MAPSIZ
      PARAMETER( M = NINT(NOP**(1.D0/3.D0)) )
      PARAMETER ( NCELL = M * M * M )
      PARAMETER ( MAPSIZ = 26 * NCELL )
      INTEGER*2   LINKS(NOP+NCELL), MAP(MAPSIZ), INCELL(NOP)
      COMMON / LNKCOM / LINKS, MAP, INCELL

      PARAMETER( RONDUP =  1.D-12 )
C**   Compute the cell-crossing times for I. The time till
C**   I crosses one of the faces of the its cell.
C**   Schedule the cell-crossing time.  Note:  Use true positions,
C**   during the calculation of the next cell crossing, but add uptime to that
C**   cell-crossing time to get it back to the basis of the last big update.

C**   Get true positions
      RXI=RX(I)+TIJNEX*VX(I)
      RYI=RY(I)+TIJNEX*VY(I)
      RZI=RZ(I)+TIJNEX*VZ(I)

C**   Apply PBC, shift coordinates, and scale for subcell size
      RXI = (RXI - DNINT(RXI) + BOXHLF)*CELLI
      RYI = (RYI - DNINT(RYI) + BOXHLF)*CELLI
      RZI = (RZI - DNINT(RZI) + BOXHLF)*CELLI

C**   Calculate the crossing time
      RXFACI=(RXI-INT(RXI)-0.5D0*( 1+SIGN(1.D0,VX(I)) ))/CELLI
      TCROSX=-RXFACI/VX(I)
      RYFACI=(RYI-INT(RYI)-0.5D0*( 1+SIGN(1.D0,VY(I)) ))/CELLI
      TCROSY=-RYFACI/VY(I)
      RZFACI=(RZI-INT(RZI)-0.5D0*( 1+SIGN(1.D0,VZ(I)) ))/CELLI
      TCROSZ=-RZFACI/VZ(I)

C**   Now decide which face is first
      TCROS = MIN(TCROSX,TCROSY,TCROSZ)
      IF (TCROS .EQ. TCROSX )THEN
         ITYPE=97
      ELSEIF (TCROS .EQ. TCROSY) THEN
         ITYPE=98
      ELSEIF (TCROS .EQ. TCROSZ) THEN
         ITYPE=99
      END IF

C**   Get back to basis of last big update
      TCROS=TCROS+TIJNEX
      TMPTIM = TCROS + RONDUP

C**   Find which cell the partcle will arrive in
      RXI=RX(I)+TMPTIM*VX(I)
      RYI=RY(I)+TMPTIM*VY(I)
      RZI=RZ(I)+TMPTIM*VZ(I)

C**   Apply PBC, shift coordinates, and scale for subcell size
      RXI = (RXI - DNINT(RXI) + BOXHLF)*CELLI
      RYI = (RYI - DNINT(RYI) + BOXHLF)*CELLI
      RZI = (RZI - DNINT(RZI) + BOXHLF)*CELLI

C**   Store the new cell number as the partner for I
      J = 1 + INT(RXI) + (INT(RYI) + INT(RZI)*M)*M

      CALL SCHED(TCROS,ITYPE,I,J)

      RETURN
      END

C**   Update the position of I and its cell-location as it
C**   undergoes a cell-crossing event.
      SUBROUTINE UPCROS(I,J,ICOTYP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "header.f"
      INTEGER     M, NCELL, MAPSIZ
      PARAMETER( M = NINT(NOP**(1.D0/3.D0)) )
      PARAMETER ( NCELL = M * M * M )
      PARAMETER ( MAPSIZ = 26 * NCELL )
      INTEGER*2   LINKS(NOP+NCELL), MAP(MAPSIZ), INCELL(NOP)
      COMMON / LNKCOM / LINKS, MAP, INCELL
      INTEGER KNEAR(100)

      PARAMETER( RONDUP = 1.D-12 )	
      PARAMETER( TIMBIG = 1D10 )

      TMPTIM=TNEXT+RONDUP

C**   TNEXT is subject to roundoff error that causes problems
C**   because the true positions will be right on the edge of a
C**   cell.  Therefore, use a temporary value which is slightly
C**   larger when deciding on which cells are neighbors etc.
C**   CELCRO will also follow this rule, but COLT will use actual
C**   value of UPTIME so that dynamics are still correct.
C**   Note that TNEXT here must be the time just after the cell-cross.
      
      ICEL= INCELL(I)
C**   Delete i-site from the linked-list of the old cell.
      INUM=LINKS(ICEL+NOP)

      IF(INUM.EQ.I)THEN
C**   Then this is the first site in the cell, so make the 
C**   head point to the 2nd site.
         LINKS(ICEL+NOP)=LINKS(INUM)
      ELSE
C**   Otherwise, we must search thru the list until we find our site.
C**   if LINKS(INUM)=I then we found our site.
         DO WHILE(LINKS(INUM).NE.I)
            INUM=LINKS(INUM)
         END DO 
         LINKS(INUM) = LINKS(I)
      END IF

C**   The link list of the old cell has now been updated.
C**   Now update the link list for the new cell.

      NEWCEL = J
      INCELL(I)=NEWCEL
      
C**   Insert I at the head of the new cell's list
      LINKS(I)=LINKS(NEWCEL+NOP)
      LINKS(NEWCEL+NOP)=I
      
C**   Since the cell-crossing event deletes this node from the tree, we
C**   must call celcro to compute and schedule the next cell-crossing.

      CALL CELCRO(TMPTIM,I)
      
C**   Now that we have the new cell location and the cell-lists have been
C**   updated, we must compute and schedule any collisions with atoms in
C**   the cells that have just become new neighbors.
C**   The face we actually just crossed is given by the collision type
C**   variable(97,98,99).  The new neighbors are along the opposite face 
C**   of the new cell. Use true positions for finding the neighbors. 
      
C**   Look in neighboring cells for possible collisions with I.
      NITER = 0

      IF (ICOTYP .EQ. 97) THEN
C**   Crossed the x-face
         IF(VX(I) .GT. 0.0) THEN
            ILO = (NEWCEL-1)*26 + 3
         ELSE
            ILO = (NEWCEL-1)*26 + 6
         ENDIF
         ISTRD = 9
         DO 60 JX = ILO,ILO+2
            IHI = JX+18
            DO 50 KK=JX,IHI,ISTRD
               ICEL=MAP(KK)
               K=LINKS(ICEL+NOP)
C**   If the head of the link list is zero,
C**   then there is nothing in this cell so go on.
               DO WHILE(K .NE. 0)
                  NITER = NITER + 1
                  KNEAR(NITER) = K
                  K=LINKS(K)
               END DO
 50         CONTINUE
 60      CONTINUE
      END IF
      
      IF (ICOTYP .EQ. 98) THEN
C**   Crossed the y-face
         IF(VY(I) .GT. 0.D0) THEN
            ILO = (NEWCEL-1)*26 + 1
         ELSE
            ILO = (NEWCEL-1)*26 + 2
         ENDIF
         IHI = ILO + 24
         ISTRD = 3
         DO 80 KK = ILO,IHI,ISTRD
            ICEL=MAP(KK)
            K=LINKS(ICEL+NOP)
            DO WHILE (K .NE. 0)
               NITER = NITER + 1
               KNEAR(NITER) = K
               K=LINKS(K)
            END DO
 80      CONTINUE
      END IF
      
      IF (ICOTYP .EQ. 99) THEN
C**   Crossed the z-face
         IF(VZ(I) .GT. 0.0) THEN
            ILO = (NEWCEL-1)*26 + 9
         ELSE
            ILO = (NEWCEL-1)*26 + 18
         ENDIF
         IHI = ILO + 8
         ISTRD = 1
         DO 101 KK = ILO,IHI,ISTRD
            ICEL=MAP(KK)
            K=LINKS(ICEL+NOP)
            DO WHILE (K .NE. 0)
               NITER = NITER + 1
               KNEAR(NITER) = K
               K=LINKS(K)
            END DO
 101     CONTINUE
      END IF
      
      CALL COLT(I,NITER,KNEAR)
      
      RETURN
      END

      SUBROUTINE NEWTIM(I,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "header.f"
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP
        INTEGER     M, NCELL, MAPSIZ
        PARAMETER( M = NINT(NOP**(1.D0/3.D0)) )
        PARAMETER ( NCELL = M * M * M )
        PARAMETER ( MAPSIZ = 26 * NCELL )
C ** This can be INTEGER*2 as long as NOP and NCELL are less than 32,768
        INTEGER*2   LINKS(NOP+NCELL), MAP(MAPSIZ), INCELL(NOP)
        COMMON / LNKCOM / LINKS, MAP, INCELL
	INTEGER KNEAR(100)

      PARAMETER(TIMBIG = 1D10)

C**   Calculate new event times for I and J.
C**   This routine is called after every core collision between I and J.
C**   Note that collision times are also calc'd in UPCROS, but only for
C**   the 9 cells that just became new neighbors.
C**   This routine is called during initial setup and from UPCOL.
C**   Note that it may be possible to schedule the i-j collision separately
C**   from the j-i collision during startup.  
C**   The event scheduler will handle this since it will
C**   wipe out all events related to the i-site before any more can be
C**   calc'd.  Therefore we only need to make sure we don't calculate the i-j
C**   collision over when we consider J.  The initial startup
C**   will schedule some events in duplicate but the duplicates will be
C**   wiped out by the event scheduler when one comes up.
      
C**   Look in neighboring cells for possible collisions with a site.

C**   Go to the head of the cell that I is in
      K = LINKS(INCELL(I) + NOP)
      
C**   Loop over all molecules in the current cell. Use do while to
C**   determine number of iterations and gather the LINKS array
      NITER = 0

      DO WHILE( K .GT. 0 )
C**   Don't do I with itself
         IF (K .NE. I) THEN
            NITER = NITER + 1
            KNEAR(NITER) = K
         END IF
         K = LINKS(K)
      END DO
      
C** Loop Over Neighboring Cells **
      JCELL0 = 26 * (INCELL(I) - 1)
      
      DO 100 NABOR = 1, 26
         JCELL = MAP ( JCELL0 + NABOR )
         
C** Loop Over All Molecules In Neighboring Cells **
         K = LINKS(JCELL + NOP)
         DO WHILE( K .NE. 0 )
            NITER = NITER + 1
            KNEAR(NITER) = K
            K = LINKS(K)
         END DO         
 100  CONTINUE

C**   Pass the array of possible collision partners to COLT
      CALL COLT(I,NITER,KNEAR)

C**   Do the same for J but don't reschedule I
C**   Go to the head of the cell that J is in
      K = LINKS(INCELL(J) + NOP)

      NITER = 0
      DO WHILE( K .GT. 0 )
C**   Don't do J with itself, nor with I
         IF (K .NE. J .AND. K .NE. I) THEN
            NITER = NITER + 1
            KNEAR(NITER) = K
         END IF
         K = LINKS(K)	      
      END DO
      
      JCELL0 = 26 * (INCELL(J) - 1)
      
      DO 101 NABOR = 1, 26
         JCELL = MAP ( JCELL0 + NABOR )
         K = LINKS(JCELL + NOP)
         DO WHILE( K .NE. 0 )
            NITER = NITER + 1
            KNEAR(NITER) = K
            K = LINKS(K)
         END DO         
 101  CONTINUE

      CALL COLT(J,NITER,KNEAR)

      RETURN
      END
      
C**   Subroutine to calculate  collision time for atom I with J
      SUBROUTINE COLT(I,NITER,KNEAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "header.f"

      INTEGER KNEAR(100), NITER

      DO 100 NN = 1,NITER
         J = KNEAR(NN)
C**   Set TIJ to a large value in case there is no collision
         VXIJ = VX(I) - VX(J)
         VYIJ = VY(I) - VY(J)
         VZIJ = VZ(I) - VZ(J)
         RXIJ = RX(I) - RX(J)
         RYIJ = RY(I) - RY(J)
         RZIJ = RZ(I) - RZ(J)

C** Take care of false positions
         RXIJ = RXIJ + VXIJ*TNEXT
         RYIJ = RYIJ + VYIJ*TNEXT
         RZIJ = RZIJ + VZIJ*TNEXT
C** Take care of PBC
	 RXIJ = RXIJ - DNINT(RXIJ)
	 RYIJ = RYIJ - DNINT(RYIJ)
	 RZIJ = RZIJ - DNINT(RZIJ)

         BIJ   = RXIJ*VXIJ + RYIJ*VYIJ + RZIJ*VZIJ
         RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
         VIJSQ = VXIJ*VXIJ + VYIJ*VYIJ + VZIJ*VZIJ
         C = RIJSQ - SIGSQ
#ifdef BELMAN
C**   J can be less than or greater than I
         IF(ABS(J-I) .EQ. 1 .AND. ICHN(I) .EQ. ICHN(J) )C = RIJSQ - SIGSQIJ
#endif
         DISCR = BIJ*BIJ - VIJSQ*C
         IF(BIJ .LT. 0.D0 .AND. DISCR .GT. 0.D0)THEN
C**   Using the formula from Numerical Recipes 5.5
            TIJ = C/(DSQRT(DISCR) - BIJ) + TNEXT
            CALL SCHED(TIJ,1,I,J)
         END IF
 100  CONTINUE

      RETURN
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C***  Binary tree algorithm for multiple event scheduling
C***  See D. C. Rapaport, J. Comp. Phys. 34, 184, (1980)

      SUBROUTINE SCHED(TIJ,ICOTYP,I,J)
C**   rapaport's algorithm for keeping the event list up to date.
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP

C**   If this is a cell-crossing, then we must insert this first.
C**   Other events for this site may come later, but we do not need
C**   to worry about connecting pointers until those events are
C**   scheduled.  so scheduling a cell-crossing is fairly simple.
      IF(ICOTYP .GT. 90)THEN
         NUNODE=I
         ISITE(NUNODE)=I
         JSITE(NUNODE)=J
         TEVENT(NUNODE)=TIJ
         COLTYP(NUNODE)=ICOTYP
         CALL ADDTRE(TIJ,NUNODE)
         RETURN
      END IF
C**   If this event is a collision, then go to work.
      
      CALL GETNOD(NUNODE)
      TEVENT(NUNODE)= TIJ
      COLTYP(NUNODE)= ICOTYP
      CALL ADDTRE(TIJ,NUNODE)
      
      ISITE(NUNODE)=I
      JSITE(NUNODE)=J
      
C**   The left node is set to the site# of the atom, meaning that we
C**   always point back to the cell-crossing event.  The other connections
C**   are arranged so that this event is added to the end of the list.
      
      ILNODE(NUNODE)=I
      IRNODE(NUNODE)=IRNODE(I)
      ILNODE(IRNODE(I))=NUNODE
      IRNODE(I)=NUNODE
      
      JLNODE(NUNODE)=J
      JRNODE(NUNODE)=JRNODE(J)
      JLNODE(JRNODE(J))=NUNODE
      JRNODE(J)=NUNODE
      
      RETURN
      END

C**   Initialize the event-scheduling list.                            
C**   This program is called when the SCHED routine is first called.

      SUBROUTINE INLIST()
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C**   NPOOL points to the first
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP
C**   node that will be needed to schedule a collision. The nodes of sites
C**   less than nsites are reserved for home's of cell-crossings.
      NPOOL=NOP+1
      NODER(0)=0
C**   NPOOL points to the next available node.
C**   Initially, we set the available nodes in ascending order.
C**   the partner pointers of each site's home node points to the
C**   home node itself initially because the list is empty (rapaport, p191)
      DO 3 I=1,NLIST-1
         NODER(I)=I+1
         ILNODE(I)=I
         IRNODE(I)=I
         JLNODE(I)=I
         JRNODE(I)=I
 3    CONTINUE
      NODER(NLIST)=0
C**   Setting NODER(NLIST)=0 will make NPOOL=0 if we run out of nodes.
      RETURN
      END


C**   Return a node to the pool
      SUBROUTINE RENODE(NODOLD)
      NODER(NODOLD)=NPOOL
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP
      NPOOL=NODOLD
      TEVENT(NODOLD)=0
      COLTYP(NODOLD)=0
      RETURN
      END

C**   Get next event.  NODER(0) points to the "remainder of the tree".
C**   ie.  it points to the last event added. (p190).
C**   if NODER(0) is the leftmost node, then it must be the next event.
C**   Otherwise, we must search down and to the left until we get to the
C**   leftmost, lowest node of any branch.

      SUBROUTINE NEXT(TIJ,ICOTYP,I,J)

      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP
      NEXNOD=NODER(0)
      DO WHILE (NODEL(NEXNOD).NE.0)
         NEXNOD=NODEL(NEXNOD)
      END DO
      
      TIJ=TEVENT(NEXNOD)
      ICOTYP= COLTYP(NEXNOD)
      I=ISITE(NEXNOD)
      J=JSITE(NEXNOD)
      IF(ICOTYP .LT. 90)THEN
         CALL CLEAR(I)
         CALL CLEAR(J)
      ELSE
         CALL DETREE(NEXNOD)
      END IF
      RETURN
      END

C**   Clear all collision event nodes involving site I by traversing both
C**   circular lists belonging to I
C**   This routine cleans up the collision partner list and calls
C**   DETREE to clean up the time list.
      
C**   First we circle left on the i-nodes, relinking the j-node
      
      SUBROUTINE CLEAR(I)
      ID=I
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP
      DO WHILE (ILNODE(ID).NE.I)
         ID=ILNODE(ID)
         JRNODE(JLNODE(ID))=JRNODE(ID)
         JLNODE(JRNODE(ID))=JLNODE(ID)
         CALL DETREE(ID)
         CALL RENODE(ID)
      END DO
      
      ILNODE(I)=I
      IRNODE(I)=I
      
C**   Then we circle left on the j-nodes, relinking the i-nodes.
      ID=I
      DO WHILE  (JLNODE(ID).NE.I)
         ID=JLNODE(ID)
         IRNODE(ILNODE(ID))=IRNODE(ID)
         ILNODE(IRNODE(ID))=ILNODE(ID)
         
         CALL DETREE(ID)
         CALL RENODE(ID)
      END DO
      
      JLNODE(I)=I
      JRNODE(I)=I

C**   Don't forget to delete cell-crossing time of the home node.      
      CALL DETREE(I)

      RETURN
      END

C**   Delete a node from the time tree and reconnect the pointers around
C**   where it was removed.

      SUBROUTINE DETREE(ID)
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP
      IF(NODER(ID).EQ.0)THEN
         IS=NODEL(ID)
      ELSE
         IF(NODEL(ID).EQ.0)THEN
            IS=NODER(ID)
         ELSE
            IF(NODEL(NODER(ID)).EQ.0)THEN
               IS=NODER(ID)
            ELSE
               IS=NODEL(NODER(ID))
               
               DO WHILE (NODEL(IS).NE.0)
                  IS=NODEL(IS)
               END DO
               
               NODEP(NODER(IS))=NODEP(IS)
               NODEL(NODEP(IS))=NODER(IS)
               NODEP(NODER(ID))=IS
               NODER(IS)=NODER(ID)
            END IF
            NODEP(NODEL(ID))=IS
            NODEL(IS)=NODEL(ID)
         END IF
      END IF
      IP=NODEP(ID)
      NODEP(IS)=IP
      IF(NODEL(IP).EQ.ID)THEN
         NODEL(IP)=IS
      ELSE
         NODER(IP)=IS
      END IF
      RETURN
      END

      
C**   Pulls a fresh node from the lowest right corner of the binary tree
      SUBROUTINE GETNOD(NUNODE)
      IF(NPOOL.EQ.0)WRITE(6,*)'ERROR POOL IS EMPTY'
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP
      IF(NPOOL.EQ.0)WRITE(0,*)'ERROR POOL IS EMPTY'
      NUNODE=NPOOL
      NPOOL=NODER(NPOOL)
      RETURN
      END


      SUBROUTINE ADDTRE(TNEW,NUNODE)
      DOUBLE PRECISION TNEW	
      DOUBLE PRECISION TIJ,TEVENT
      PARAMETER(NLIST=7*NOP)
      INTEGER*2 COLTYP,JSITE,ISITE
      INTEGER*4 NODEL,NODEP,NODER,ILNODE,IRNODE,JLNODE,JRNODE
      DIMENSION TEVENT(0:NLIST),NODEL(0:NLIST),NODEP(0:NLIST),
     + NODER(0:NLIST),ILNODE(0:NLIST),ISITE(0:NLIST),IRNODE(0:NLIST),
     +                JLNODE(0:NLIST),JSITE(0:NLIST),JRNODE(0:NLIST),
     +                COLTYP(0:NLIST)
      COMMON /TREE/ NODEP,NODEL,NODER,ILNODE,IRNODE,JLNODE,JRNODE,NPOOL
      COMMON /EVENTS/ TEVENT,ISITE,JSITE,COLTYP
C**   Puts the new node into the binary tree event list
C**   Note that TNEW is required but is not stored via this sub.
C**   TNEW is the time of the event going to happen.
C**   Start by checking if NODER(0)=0.  Then the tree must be empty and we
C**   can just put this event in first.
      INPO=0
      IF(NODER(0).EQ.0)THEN
         NODER(0)=NUNODE
      ELSE
C**   Otherwise, we must search thru all the branches to find the proper
C**   branch-end to insert this event.  The idea is that the branch to
C**   the right should have a larger time for all branches.  Note that the
C**   ENDIF statement is way down in the program.  Do 10 goes thru the
C**   entire list in principle but the go to 11 stops us if NFOUND=1.
         INPO=NODER(0)
         NFOUND=0
         
         DO WHILE(NFOUND.EQ.0)
            
C**   if TNEW less than TEVENT(INPO) then look to the left for a branch-end.
            IF(TNEW.LE.TEVENT(INPO))THEN
               IF(NODEL(INPO).NE.0)THEN
                  INPO=NODEL(INPO)
               ELSE
                  NODEL(INPO)=NUNODE
                  NFOUND=1
               END IF
               
C**   if TNEW gretaer than TEVENT(INPO) then look to the right for a branch-end.
            ELSE
               IF(NODER(INPO).NE.0)THEN
                  INPO=NODER(INPO)
               ELSE
                  NODER(INPO)=NUNODE
                  NFOUND=1
               END IF
            END IF
         END DO
      END IF
      

C**   INPO now equals the position preceding where the event was added.
C**   NODEL and NODER are set to zero since this new node is at the end of a
C**   branch.
      NODEL(NUNODE)=0
      NODER(NUNODE)=0
C**   NODEP is set to INPO because that is where this node should go if any
C**   future deletions in the tree cause INPO to be moved.
      NODEP(NUNODE)=INPO

      RETURN
      END



C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C**   Construction of cell linked-lists

C**   References:
C**   Quentrec and Brot, J. Comput. Phys. 13, 430, 1975.
C**   Hockney and Eastwood, Comp. Sim. using Particles,McGraw Hill, 1981.

C**   
C**   Subroutine MAPS
C**   Sets up map of cell structure (26 for each cell)
C**   Subroutine CREATE_LINKS
C**   Sets up head of chain array and linked list

C**   Subroutine maps is called once at the start of a simulation
C**   to establish cell neighbour identities.  
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C**   Routine to set up a map of neighboring cells

C**   Variables:
C**   INTEGER M                  Number of cells in each direction
C**   INTEGER MAPSIZ             Size of cell-cell map 
C**   INTEGER MAP(MAPSIZ)        List of neighboring cells

C**   This subroutine sets up a list of the 26 neighboring
C**   cells for each small cell in the central box. The
C**   effects of the periodic boundary conditions are included.
C**   The subroutine is called once at the beginning of the
C**   simulation and the map is used in the cell crossing routines.


      SUBROUTINE MAPS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER     M, NCELL, MAPSIZ, NP
      PARAMETER( NP = NOP)

      PARAMETER( M = NINT(NOP**(1.D0/3.D0)) )
      PARAMETER ( NCELL = M * M * M )
      PARAMETER ( MAPSIZ = 26 * NCELL )
      INTEGER*2   LINKS(NOP+NCELL), MAP(MAPSIZ), INCELL(NOP)
      COMMON / LNKCOM / LINKS, MAP, INCELL
      INTEGER KNEAR(100)

      INTEGER     IX, IY, IZ, IMAP, ICELL

C**   STATEMENT FUNCTION TO GIVE CELL INDEX **

      ICELL ( IX, IY, IZ) = 1 + MOD ( IX - 1 + M, M )
     :     + MOD ( IY - 1 + M, M ) * M
     :     + MOD ( IZ - 1 + M, M ) * M * M

C**   FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **

      DO 50 IZ = 1, M

         DO 40 IY = 1, M

            DO 30 IX = 1, M

               IMAP = ( ICELL ( IX, IY, IZ ) - 1 ) * 26
C**   Same Plane
               MAP( IMAP + 1  ) = ICELL( IX    , IY + 1, IZ     )
               MAP( IMAP + 2  ) = ICELL( IX    , IY - 1, IZ     )
               MAP( IMAP + 3  ) = ICELL( IX + 1, IY    , IZ     )
               MAP( IMAP + 4  ) = ICELL( IX + 1, IY + 1, IZ     )
               MAP( IMAP + 5  ) = ICELL( IX + 1, IY - 1, IZ     )
               MAP( IMAP + 6  ) = ICELL( IX - 1, IY    , IZ     )
               MAP( IMAP + 7  ) = ICELL( IX - 1, IY + 1, IZ     )
               MAP( IMAP + 8  ) = ICELL( IX - 1, IY - 1, IZ     )
C**   Plane in front
               MAP( IMAP + 9  ) = ICELL( IX    , IY    , IZ + 1 )
               MAP( IMAP + 10 ) = ICELL( IX    , IY + 1, IZ + 1 )
               MAP( IMAP + 11 ) = ICELL( IX    , IY - 1, IZ + 1 )
               MAP( IMAP + 12 ) = ICELL( IX + 1, IY    , IZ + 1 )
               MAP( IMAP + 13 ) = ICELL( IX + 1, IY + 1, IZ + 1 )
               MAP( IMAP + 14 ) = ICELL( IX + 1, IY - 1, IZ + 1 )
               MAP( IMAP + 15 ) = ICELL( IX - 1, IY    , IZ + 1 )
               MAP( IMAP + 16 ) = ICELL( IX - 1, IY + 1, IZ + 1 )
               MAP( IMAP + 17 ) = ICELL( IX - 1, IY - 1, IZ + 1 )
C**   Plane in rear
               MAP( IMAP + 18 ) = ICELL( IX    , IY    , IZ - 1 )
               MAP( IMAP + 19 ) = ICELL( IX    , IY + 1, IZ - 1 )
               MAP( IMAP + 20 ) = ICELL( IX    , IY - 1, IZ - 1 )
               MAP( IMAP + 21 ) = ICELL( IX + 1, IY    , IZ - 1 )
               MAP( IMAP + 22 ) = ICELL( IX + 1, IY + 1, IZ - 1 )
               MAP( IMAP + 23 ) = ICELL( IX + 1, IY - 1, IZ - 1 )
               MAP( IMAP + 24 ) = ICELL( IX - 1, IY    , IZ - 1 )
               MAP( IMAP + 25 ) = ICELL( IX - 1, IY + 1, IZ - 1 )
               MAP( IMAP + 26 ) = ICELL( IX - 1, IY - 1, IZ - 1 )
 30         CONTINUE

 40      CONTINUE

 50   CONTINUE

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
#include "header.f"
      INTEGER     M, NCELL, MAPSIZ, NP
      PARAMETER( NP = NOP)
      PARAMETER( M = NINT(NOP**(1.D0/3.D0)) )
      PARAMETER ( NCELL = M * M * M )
      PARAMETER ( MAPSIZ = 26 * NCELL )
      INTEGER*2   LINKS(NOP+NCELL), MAP(MAPSIZ), INCELL(NOP)
      COMMON / LNKCOM / LINKS, MAP, INCELL

C**   Zero Head of Chain Array **
      DO 10 I = NOP + 1, NOP + NCELL
         LINKS(I) = 0
 10   CONTINUE

C**   Sort all Atoms **
      DO 20 I = 1, NOP
         RXI = RX(I) - DNINT(RX(I)) + BOXHLF
         RYI = RY(I) - DNINT(RY(I)) + BOXHLF
         RZI = RZ(I) - DNINT(RZ(I)) + BOXHLF
         INCELL(I) = 1 + INT(RXI*CELLI) + (INT(RYI*CELLI) 
     &        + INT(RZI*CELLI)*M)*M
         LINKS(I)    = LINKS(INCELL(I)+NOP)
         LINKS(INCELL(I)+NOP) = I
 20   CONTINUE

      RETURN
      END

</CODE>
</PRE>
