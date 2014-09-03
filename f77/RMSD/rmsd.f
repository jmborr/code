*************************************************************************
*
*     pgf77 -Mextend -O -s -fast -Wl,-static -o rmsd rmsd.f
*
*     This program is to calculate RMSD of file1, file2
*     finf out common residues and then calculate RMSD
*     pdb3=file1(rotate)+file2
*************************************************************************

      program compares
      PARAMETER(nmax=3000)

      dimension xx(nmax),yy(nmax),zz(nmax)
      character*100 s,du
      character*200 fnam1,fnam2,fnam3
      character*3 aa1(nmax),aa2(nmax)

      dimension mm1(nmax),mm2(nmax) !aligned residues in origional order
      dimension x1(nmax),y1(nmax),z1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

*****instructions ----------------->
      call getarg(1,fnam1)
      if(fnam1.eq.' '.or.fnam1.eq.'?'.or.fnam1.eq.'-h')then
         write(*,*)
         write(*,*)'Brief instruction for running rmsd program:'
         write(*,*)
         write(*,*)'1. Run rmsd to compare ''model1'' and ''model2'':'
         write(*,*)'  rmsd model1 model2'
         write(*,*)
         write(*,*)'2. Run rmsd with superposition in ''rmsd.sup'':'
         write(*,*)'  rmsd model1 model2 1'
         write(*,*)
         goto 9999
      endif
      call getarg(2,fnam2)
      call getarg(3,fnam3)
      
cccccccccread data from first CA file:
      open(unit=10,file=fnam1,status='old')
      i=0
      do while (.true.)
         read(10,9001,end=1010) s
         if(i.gt.0.and.s(1:3).eq.'TER')goto 1010
         if(s(1:3).eq.'ATO')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.
     &           s(13:16).eq.'  CA')then
               i=i+1
               read(s,9000)du,aa1(i),du,mm1(i),du,
     $              x1(i),y1(i),z1(i)
            endif
         endif
      enddo
 1010 continue
 9000 format(A17,A3,A2,i4,A4,3F8.3)
 9001 format(A100)
      close(10)
      nseq1=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

cccccccccread data from second CA file:
      open(unit=10,file=fnam2,status='old')
      i=0
      do while (.true.)
         read(10,9001,end=1011) s
         if(i.gt.0.and.s(1:3).eq.'TER')goto 1011
         if(s(1:3).eq.'ATO')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.
     &           s(13:16).eq.'  CA')then
               i=i+1
               read(s,9000)du,aa2(i),du,mm2(i),du,
     $              x2(i),y2(i),z2(i)
            endif
         endif
      enddo
 1011 continue
      close(10)
      nseq2=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ccc   find the common residues for RMSD calculations -------->
      k=0
      do i=1,nseq1
         do j=1,nseq2
            if(mm1(i).eq.mm2(j))then
               k=k+1
               r_1(1,k)=x1(i)
               r_1(2,k)=y1(i)
               r_1(3,k)=z1(i)
               r_2(1,k)=x2(j)
               r_2(2,k)=y2(j)
               r_2(3,k)=z2(j)
c               write(*,*)k, r_1(1,k),r_2(1,k)
            endif
         enddo
      enddo
      nseq=k
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      call u3b(w,r_1,r_2,nseq,1,rms,u,t,ier) !u rotate r_1 to r_2
      armsd=dsqrt(rms/nseq)
***   rotate r_1:
      do j=1,nseq1
         xx(j)=t(1)+u(1,1)*x1(j)+u(1,2)*y1(j)+u(1,3)*z1(j)
         yy(j)=t(2)+u(2,1)*x1(j)+u(2,2)*y1(j)+u(2,3)*z1(j)
         zz(j)=t(3)+u(3,1)*x1(j)+u(3,2)*y1(j)+u(3,3)*z1(j)
      enddo
***   output rotated chain1----->
      if(fnam3.ne.'1')goto 999
      OPEN(unit=7,file='rmsd.sup',status='unknown') !file1(rotate)+file2
      do i=1,nseq1
c         write(7,1237)i,aa1(i),mm1(i),xx(i),yy(i),zz(i)
         write(7,1237)mm1(i),aa1(i),mm1(i),xx(i),yy(i),zz(i)
      enddo
      write(7,1238)
      do i=2,nseq1
         write(7,1239)mm1(i-1),mm1(i) !connection for atom
      enddo
***   output chain2------------->
      do i=1,nseq2
         write(7,1237)1000+mm2(i),aa2(i),mm2(i),x2(i),y2(i),z2(i)
      enddo
      write(7,1238)
      do i=2,nseq2
         write(7,1239)1000+mm2(i-1),1000+mm2(i)
      enddo
 999  continue
*^^^^^^^^^^^^^^^^^^ rotation finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 1237 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
 1238 format('TER')
 1239 format('CONECT',I5,I5)

      write(*,*)
      write(*,101)fnam1,nseq1
 101  format('Chain 1:',A10,'  Size=',I4)
      write(*,102)fnam2,nseq2
 102  format('Chain 2:',A10,'  Size=',I4)
      write(*,104)nseq
 104  format('Common length:',I4)
      write(*,*)
      write(*,103)armsd
 103  format('RMSD=',f7.3)
      write(*,*)
      write(*,*)'Rotation matrix'
      write(*,*) u(1,1),u(2,1),u(3,1)
      write(*,*) u(1,2),u(2,2),u(3,2)
      write(*,*) u(1,3),u(2,2),u(3,3)

 9999 END

cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m           (given)
c  x    - x(i,m) are coordinates of atom c m in set x       (given)
c  y    - y(i,m) are coordinates of atom c m in set y       (given)
c  n    - n is number of atom pairs                         (given)
c  mode  - 0:calculate rms only                             (given)
c          1:calculate rms,u,t                              (takes longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      double precision w(1), x(3, 1), y(3, 1), u(3, 3), t(3), rms, sigma
      double precision r(3, 3), xc(3), yc(3), wc, a(3, 3), b(3, 3), e0, 
     &e(3), e1, e2, e3, d, spur, det, cof, h, g, cth, sth, sqrth, p, tol
     &, rr(6), rr1, rr2, rr3, rr4, rr5, rr6, ss(6), ss1, ss2, ss3, ss4, 
     &ss5, ss6, zero, one, two, three, sqrt3
      equivalence (rr(1), rr1), (rr(2), rr2), (rr(3), rr3), (rr(4), rr4)
     &, (rr(5), rr5), (rr(6), rr6), (ss(1), ss1), (ss(2), ss2), (ss(3), 
     &ss3), (ss(4), ss4), (ss(5), ss5), (ss(6), ss6), (e(1), e1), (e(2)
     &, e2), (e(3), e3)
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data one / 1.0d+00 /
      data two / 2.0d+00 /
      data three / 3.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
c 156 "rms.for"
      wc = zero
      rms = 0.0
      e0 = zero
      do 1 i = 1, 3
      xc(i) = zero
      yc(i) = zero
      t(i) = 0.0
      do 1 j = 1, 3
      d = zero
      if (i .eq. j) d = one
      u(i,j) = d
      a(i,j) = d
    1 r(i,j) = zero
      ier = -1
c**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
c 170 "rms.for"
      if (n .lt. 1) return 
c 172 "rms.for"
      ier = -2
      do 2 m = 1, n
      if (w(m) .lt. 0.0) return 
      wc = wc + w(m)
      do 2 i = 1, 3
      xc(i) = xc(i) + (w(m) * x(i,m))
    2 yc(i) = yc(i) + (w(m) * y(i,m))
      if (wc .le. zero) return 
      do 3 i = 1, 3
      xc(i) = xc(i) / wc
c**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
c 182 "rms.for"
    3 yc(i) = yc(i) / wc
c 184 "rms.for"
      do 4 m = 1, n
      do 4 i = 1, 3
      e0 = e0 + (w(m) * (((x(i,m) - xc(i)) ** 2) + ((y(i,m) - yc(i)) ** 
     &2)))
c 187 "rms.for"
      d = w(m) * (y(i,m) - yc(i))
      do 4 j = 1, 3
c**** CALCULATE DETERMINANT OF R(I,J)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
c 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2
     &) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3) * ((r(2,1)
     & * r(3,2)) - (r(2,2) * r(3,1))))
c**** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
c 194 "rms.for"
      sigma = det
c 196 "rms.for"
      m = 0
      do 5 j = 1, 3
      do 5 i = 1, j
      m = m + 1
c***************** EIGENVALUES *****************************************
c**** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
c 200 "rms.for"
    5 rr(m) = ((r(1,i) * r(1,j)) + (r(2,i) * r(2,j))) + (r(3,i) * r(3,j)
     &)
c 203 "rms.for"
      spur = ((rr1 + rr3) + rr6) / three
      cof = ((((((rr3 * rr6) - (rr5 * rr5)) + (rr1 * rr6)) - (rr4 * rr4)
     &) + (rr1 * rr3)) - (rr2 * rr2)) / three
c 205 "rms.for"
      det = det * det
      do 6 i = 1, 3
    6 e(i) = spur
c**** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR
c 208 "rms.for"
      if (spur .le. zero) goto 40
c 210 "rms.for"
      d = spur * spur
      h = d - cof
c**** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
c 212 "rms.for"
      g = (((spur * cof) - det) / two) - (spur * h)
c 214 "rms.for"
      if (h .le. zero) goto 8
      sqrth = dsqrt(h)
      d = ((h * h) * h) - (g * g)
      if (d .lt. zero) d = zero
      d = datan2(dsqrt(d),- g) / three
      cth = sqrth * dcos(d)
      sth = (sqrth * sqrt3) * dsin(d)
      e1 = (spur + cth) + cth
      e2 = (spur - cth) + sth
      e3 = (spur - cth) - sth
c.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
c 224 "rms.for"
      if (mode) 10, 50, 10
c**************** EIGENVECTORS *****************************************
c 226 "rms.for"
    8 if (mode) 30, 50, 30
c 228 "rms.for"
   10 do 15 l = 1, 3, 2
      d = e(l)
      ss1 = ((d - rr3) * (d - rr6)) - (rr5 * rr5)
      ss2 = ((d - rr6) * rr2) + (rr4 * rr5)
      ss3 = ((d - rr1) * (d - rr6)) - (rr4 * rr4)
      ss4 = ((d - rr3) * rr4) + (rr2 * rr5)
      ss5 = ((d - rr1) * rr5) + (rr2 * rr4)
      ss6 = ((d - rr1) * (d - rr3)) - (rr2 * rr2)
      j = 1
      if (dabs(ss1) .ge. dabs(ss3)) goto 12
      j = 2
      if (dabs(ss3) .ge. dabs(ss6)) goto 13
   11 j = 3
      goto 13
   12 if (dabs(ss1) .lt. dabs(ss6)) goto 11
   13 d = zero
      j = 3 * (j - 1)
      do 14 i = 1, 3
      k = ip(i + j)
      a(i,l) = ss(k)
   14 d = d + (ss(k) * ss(k))
      if (d .gt. zero) d = one / dsqrt(d)
      do 15 i = 1, 3
   15 a(i,l) = a(i,l) * d
      d = ((a(1,1) * a(1,3)) + (a(2,1) * a(2,3))) + (a(3,1) * a(3,3))
      m1 = 3
      m = 1
      if ((e1 - e2) .gt. (e2 - e3)) goto 16
      m1 = 1
      m = 3
   16 p = zero
      do 17 i = 1, 3
      a(i,m1) = a(i,m1) - (d * a(i,m))
   17 p = p + (a(i,m1) ** 2)
      if (p .le. tol) goto 19
      p = one / dsqrt(p)
      do 18 i = 1, 3
   18 a(i,m1) = a(i,m1) * p
      goto 21
   19 p = one
      do 20 i = 1, 3
      if (p .lt. dabs(a(i,m))) goto 20
      p = dabs(a(i,m))
      j = i
   20 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((a(k,m) ** 2) + (a(l,m) ** 2))
      if (p .le. tol) goto 40
      a(j,m1) = zero
      a(k,m1) = - (a(l,m) / p)
      a(l,m1) = a(k,m) / p
   21 a(1,2) = (a(2,3) * a(3,1)) - (a(2,1) * a(3,3))
      a(2,2) = (a(3,3) * a(1,1)) - (a(3,1) * a(1,3))
c****************** ROTATION MATRIX ************************************
c 282 "rms.for"
      a(3,2) = (a(1,3) * a(2,1)) - (a(1,1) * a(2,3))
c 284 "rms.for"
   30 do 32 l = 1, 2
      d = zero
      do 31 i = 1, 3
      b(i,l) = ((r(i,1) * a(1,l)) + (r(i,2) * a(2,l))) + (r(i,3) * a(3,l
     &))
c 288 "rms.for"
   31 d = d + (b(i,l) ** 2)
      if (d .gt. zero) d = one / dsqrt(d)
      do 32 i = 1, 3
   32 b(i,l) = b(i,l) * d
      d = ((b(1,1) * b(1,2)) + (b(2,1) * b(2,2))) + (b(3,1) * b(3,2))
      p = zero
      do 33 i = 1, 3
      b(i,2) = b(i,2) - (d * b(i,1))
   33 p = p + (b(i,2) ** 2)
      if (p .le. tol) goto 35
      p = one / dsqrt(p)
      do 34 i = 1, 3
   34 b(i,2) = b(i,2) * p
      goto 37
   35 p = one
      do 36 i = 1, 3
      if (p .lt. dabs(b(i,1))) goto 36
      p = dabs(b(i,1))
      j = i
   36 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((b(k,1) ** 2) + (b(l,1) ** 2))
      if (p .le. tol) goto 40
      b(j,2) = zero
      b(k,2) = - (b(l,1) / p)
      b(l,2) = b(k,1) / p
   37 b(1,3) = (b(2,1) * b(3,2)) - (b(2,2) * b(3,1))
      b(2,3) = (b(3,1) * b(1,2)) - (b(3,2) * b(1,1))
      b(3,3) = (b(1,1) * b(2,2)) - (b(1,2) * b(2,1))
      do 39 i = 1, 3
      do 39 j = 1, 3
c****************** TRANSLATION VECTOR *********************************
c 320 "rms.for"
   39 u(i,j) = ((b(i,1) * a(j,1)) + (b(i,2) * a(j,2))) + (b(i,3) * a(j,3
     &))
   40 do 41 i = 1, 3
c********************** RMS ERROR **************************************
c 323 "rms.for"
   41 t(i) = ((yc(i) - (u(i,1) * xc(1))) - (u(i,2) * xc(2))) - (u(i,3)
     & * xc(3))
   50 do 51 i = 1, 3
      if (e(i) .lt. zero) e(i) = zero
   51 e(i) = dsqrt(e(i))
      ier = 0
      if (e2 .le. (e1 * 1.0d-05)) ier = -1
      d = e3
      if (sigma .ge. 0.0) goto 52
      d = - d
      if ((e2 - e3) .le. (e1 * 1.0d-05)) ier = -1
   52 d = (d + e2) + e1
      rms = (e0 - d) - d
      if (rms .lt. 0.0) rms = 0.0
      return 
c.....END U3B...........................................................
c----------------------------------------------------------
c                       THE END
c----------------------------------------------------------
c 338 "rms.for"
      end
