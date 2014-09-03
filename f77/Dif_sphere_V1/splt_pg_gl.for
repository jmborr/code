      SUBROUTINE SPLT(DXMIN,DXMAX,DYMIN,DYMAX,X,Y,NSPNT)
      COMMON/CALDAT/A,YCALC,WT
      COMMON/TOSPLT/RESID,R1,R2,CURSOR
      COMMON/PLTINF/SPEED,RTEMP,STEMP
      COMMON/FILES/RFILE,LENRF,SFILE,LENSF,NUMSPC

      LOGICAL CURSOR,pglog
      REAL*8 DXMIN,DXMAX,DYMIN,DYMAX,RTEMP,STEMP
      REAL*8 RESID(5000),A(100),YCALC(5000),WT(5000)
      REAL*8 X(NSPNT),Y(NSPNT)
      REAL TEMPX(5000),TEMPY(5000),tempye(5000)
      INTEGER R1,R2,SPEED,NUMPNT,TITLE(20),NUMSPC
      CHARACTER*20 RFILE,SFILE,LINE*80, Device*10, answer*1

C CONVERT TO SINGLE PRECISION
      XMIN=DXMIN
      XMAX=DXMAX
      YMIN=DYMIN
      YMAX=DYMAX

      NUMPNT=R2-R1+1
c the only thing that should need moifying in this routine is to remove the 
c pgcalls to set up the device and always open the windows graphics


c	call pgbegin(0, '/gw', 1, 1)
c	call pgbeg(0, '/gw', 1, 1)
c	pglog=.FALSE.
c	call pgask(pglog)
      call pgpage
      call pgvport(0.1, 0.95, 0.25,0.85)

      LENLIN=14+LENSF
      LINE(1:LENLIN)='SAMPLE FILE = '//SFILE(1:LENSF)
      call pgmtext('T', 4.0, 0.5, 0.5, Line(1:lenlin))
      LENLIN=18+LENRF
      LINE(1:LENLIN)='RESOLUTION FILE = '//RFILE(1:LENRF)
      call pgmtext('T', 3.0, 0.5, 0.5, Line(1:lenlin))
      LENLIN=20
      WRITE (LINE(1:LENLIN),40),NUMSPC
40    FORMAT('SPECTRUM NUMBER = ',I2)
      call pgmtext('T', 2.0, 0.5, 0.5, Line(1:lenlin))

      call pgwindow(xmin, xmax, ymin, ymax)
      call pglabel(' ', 'Intensity', Line(1:lenlin))
      call pgbox('BCST', 0.0, 0, 'BCNST', 0.0, 0)

      ii=1
      do 60 i=1,nspnt
        if (x(i).gt.xmin) then
          ii=i
          goto 80
        endif
60    continue
80    continue

      j=0
      do 100 i=ii,nspnt
        if (x(i).gt.xmax) goto 120
        j=j+1
        tempx(j)=x(i)
        tempy(j)=y(i)
        tempye(j)=wt(i)
100   continue
120   continue
      numpnt=j

      call pgpoint(numpnt, tempx, tempy, 20)
      do 140 i=1,numpnt
      ylow = tempy(i) - tempye(i)
      yhigh = tempy(i) + tempye(i)
      call pgerry(1, tempx(i), ylow, yhigh, 1.0)
140   continue

      RMAX=0
      DO 160 I=1,NUMPNT
      IPOS=R1+I-1
      if (resid(ipos).gt.rmax) rmax=resid(ipos)
      TEMPX(I)=X(IPOS)
      TEMPY(I)=YCALC(IPOS)
160   CONTINUE
      call pgsci(2)
      CALL pgline(numpnt, TEMPX,TEMPY)
      call pgsci(1)

      IF (RMAX.GT.1E-5) THEN
          call pgvport(0.1, 0.95, 0.1,0.25)
          RMAX=INT(RMAX)+1.0
          call pgwindow(xmin, xmax, -rmax, rmax)
          call pglabel('Energy Transfer, meV', 'Residuals', ' ')
          call pgbox('BCNST', 0.0, 0, 'BCNST', 0.0, 0)
          call pgsls(2)
          call pgmove(xmin, 0.0)
          call pgdraw(xmax, 0.0)
          call pgsls(1)
          IPOS=1
          DO 180 I=R1,R2
              TEMPY(IPOS)=RESID(I)
              IPOS=IPOS+1
180       CONTINUE

          call pgpoint(numpnt, tempx, tempy, 20)
          CALL pgline(numpnt, TEMPX, TEMPY)
      ENDIF

c      CALL pgend

      read (5,200) answer
200   format(a)

      END



      SUBROUTINE HPLT(DXMIN,DXMAX,DYMIN,DYMAX,X,Y,NSPNT)
c don't use this anymore
c      COMMON/CALDAT/A,YCALC,WT
c      COMMON/TOSPLT/RESID,R1,R2,CURSOR
c      COMMON/PLTINF/SPEED,RTEMP,STEMP
c      COMMON/FILES/RFILE,LENRF,SFILE,LENSF,NUMSPC

c      LOGICAL CURSOR
c      REAL*8 DXMIN,DXMAX,DYMIN,DYMAX,RTEMP,STEMP
c      REAL*8 RESID(5000),A(80),YCALC(5000),WT(5000)
c      REAL*8 X(NSPNT),Y(NSPNT)
c      REAL TEMPX(5000),TEMPY(5000),tempye(5000)
c      INTEGER R1,R2,SPEED,NUMPNT,TITLE(20),NUMSPC, lennam
c      CHARACTER*20 RFILE,SFILE,LINE*80, pltfil*20

C CONVERT TO SINGLE PRECISION
c      XMIN=DXMIN
c      XMAX=DXMAX
c      YMIN=DYMIN
c      YMAX=DYMAX

c First extract the sample filename to create the plot file
c This is required to avoid extensions to the file name such as .dat
c      do 10 i=1,lensf
c          if (sfile(i:i).eq.'.') then
c              lennam = i - 1
c              goto 15
c          endif
c10    continue
c      lennam = lensf
c15    continue
c      pltfil(1:lennam + 7) = sfile(1:lennam)//'.plt/ps'
c      lennam = lennam + 7

c      write (6,23) pltfil(1:lennam)
c23    format(1x,'Writing postscript file to ',A)

c      NUMPNT=R2-R1+1
c      CALL pgbegin(0, pltfil(1:lennam), 1, 1)
c      call pgpage
c      call pgvport(0.15, 0.85, 0.35,0.85)

c      LENLIN=14+LENSF
c      LINE(1:LENLIN)='SAMPLE FILE = '//SFILE(1:LENSF)
c      call pgmtext('T', 4.0, 0.5, 0.5, Line(1:lenlin))
c      LENLIN=18+LENRF
c      LINE(1:LENLIN)='RESOLUTION FILE = '//RFILE(1:LENRF)
c      call pgmtext('T', 3.0, 0.5, 0.5, Line(1:lenlin))
c      LENLIN=20
c      WRITE (LINE(1:LENLIN),20),NUMSPC
c20    FORMAT('SPECTRUM NUMBER = ',I2)
c      call pgmtext('T', 2.0, 0.5, 0.5, Line(1:lenlin))

c      call pgwindow(xmin, xmax, ymin, ymax)
c      call pglabel(' ', 'Intensity', Line(1:lenlin))
c      call pgbox('BCST', 0.0, 0, 'BCNST', 0.0, 0)

c      ii=1
c      do 40 i=1,nspnt
c        if (x(i).gt.xmin) then
c          ii=i
c          goto 60
c        endif
c40    continue
c60    continue
c
c      j=0
c      do 80 i=ii,nspnt
c        if (x(i).gt.xmax) goto 100
c        j=j+1
c        tempx(j)=x(i)
c        tempy(j)=y(i)
c        tempye(j)=wt(i)
c80    continue
c100   continue
c      numpnt=j

c      call pgpoint(numpnt, tempx, tempy, 20)
c      do 120 i=1,numpnt
c      ylow = tempy(i) - tempye(i)
c      yhigh = tempy(i) + tempye(i)
c      call pgerry(1, tempx(i), ylow, yhigh, 1.0)
c120   continue

c      RMAX=0
c      DO 140 I=1,NUMPNT
c      IPOS=R1+I-1
c      if (resid(ipos).gt.rmax) rmax=resid(ipos)
c      TEMPX(I)=X(IPOS)
c      TEMPY(I)=YCALC(IPOS)
c140   CONTINUE
c      call pgsci(2)
c      CALL pgline(numpnt, TEMPX,TEMPY)
c      call pgsci(1)

c      IF (RMAX.GT.1E-5) THEN
c          call pgvport(0.15, 0.85, 0.15,0.35)
c          RMAX=INT(RMAX)+1.0
c          call pgwindow(xmin, xmax, -rmax, rmax)
c          call pglabel('Energy Transfer, meV', 'Residuals', ' ')
c          call pgbox('BCNST', 0.0, 0, 'BCST', 0.0, 0)
c          call pgsls(2)
c          call pgmove(xmin, 0.0)
c          call pgdraw(xmax, 0.0)
c          call pgsls(1)
c          IPOS=1
c          DO 160 I=R1,R2
c              TEMPY(IPOS)=RESID(I)
c              IPOS=IPOS+1
c160       CONTINUE

c          call pgpoint(numpnt, tempx, tempy, 20)
c          CALL pgline(numpnt, TEMPX, TEMPY)
c      ENDIF

c      CALL pgend

      END




