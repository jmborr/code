c     ###################################
      real*8 function uniform(d1,d2)
      real*8 d1,d2,rn
          uniform=(d2-d1)*rn()+d1
      return
      end      




c#################################################
c Uniformlu Distributed Random Number Generator 
c Output:
c     Double Precision RN from 0 to 1:
c Input (changed by Han 2/13/92):
c     Iseed --- random seed for initialization
c     the range should be 0<= Iseed <=31328
c
      FUNCTION RAN(ISEED)
      REAL*8 RAN,RAN1
      INTEGER*4 ISEED
      SAVE INIT
      DATA INIT /1/
      IF (INIT.EQ.1) THEN
        INIT=0
        CALL RMARIN(ISEED,23456)
      END IF

  10  CALL RANMAR(RAN1)
      IF (RAN1.LT.1D-16) GOTO 10
      RAN=RAN1
	RETURN
      END

      FUNCTION RN()
      REAL*8 RN,RAN1
      INTEGER*4 ISEED
      COMMON/SEED/ISEED
      SAVE INIT
      DATA INIT /1/
      IF (INIT.EQ.1) THEN
        INIT=0
        CALL RMARIN(ISEED,23456)
      END IF

  10  CALL RANMAR(RAN1)
      IF (RAN1.LT.1D-16) GOTO 10
      RN=RAN1
	RETURN
      END
C	##############################
      SUBROUTINE RANMAR(RVEC)
*     -----------------
* Universal random number generator proposed by Marsaglia and Zaman
* in report FSU-SCRI-87-50
* In this version RVEC is a double precision variable.
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
      COMMON/ RASET2 / IRANMR,JRANMR
C      SAVE /RASET1/,/RASET2/
      UNI = RANU(IRANMR) - RANU(JRANMR)
      IF(UNI .LT. 0D0) UNI = UNI + 1D0
      RANU(IRANMR) = UNI
      IRANMR = IRANMR - 1
      JRANMR = JRANMR - 1
      IF(IRANMR .EQ. 0) IRANMR = 97
      IF(JRANMR .EQ. 0) JRANMR = 97
      RANC = RANC - RANCD
      IF(RANC .LT. 0D0) RANC = RANC + RANCM
      UNI = UNI - RANC
      IF(UNI .LT. 0D0) UNI = UNI + 1D0
      RVEC = UNI
      END
 
      SUBROUTINE RMARIN(IJ,KL)
*     -----------------
* Initializing routine for RANMAR, must be called before generating
* any pseudorandom numbers with RANMAR. The input values should be in
* the ranges 0<=ij<=31328 ; 0<=kl<=30081
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
      COMMON/ RASET2 / IRANMR,JRANMR
      SAVE /RASET1/,/RASET2/
* This shows correspondence between the simplified input seeds IJ, KL
* and the original Marsaglia-Zaman seeds I,J,K,L.
* To get the standard values in the Marsaglia-Zaman paper (i=12,j=34
* k=56,l=78) put ij=1802, kl=9373
      I = MOD( IJ/177 , 177 ) + 2
      J = MOD( IJ     , 177 ) + 2
      K = MOD( KL/169 , 178 ) + 1
      L = MOD( KL     , 169 )
      DO 300 II = 1 , 97
        S =  0D0
        T = .5D0
        DO 200 JJ = 1 , 24
          M = MOD( MOD(I*J,179)*K , 179 )
          I = J
          J = K
          K = M
          L = MOD( 53*L+1 , 169 )
          IF(MOD(L*M,64) .GE. 32) S = S + T
          T = .5D0*T
  200   CONTINUE
        RANU(II) = S
  300 CONTINUE
      RANC  =   362436D0 / 16777216D0
      RANCD =  7654321D0 / 16777216D0
      RANCM = 16777213D0 / 16777216D0
      IRANMR = 97
      JRANMR = 33
      END
c       ###################
          subroutine calrg(rg)
          include "frg_com.f"
          real*8 rg,cx,cy,cz,xnatom

          cx=0.
          cy=0.
          cz=0.
          do i=1,natom
          cx=cx+xp(i)
          cy=cy+yp(i)
          cz=cz+zp(i)
          enddo
          xnatom=1./natom
          cx=cx*xnatom
          cy=cy*xnatom
          cz=cz*xnatom
          rg=0.0
          do i=1,natom
          rg=rg+(xp(i)-cx)**2+(yp(i)-cy)**2+(zp(i)-cz)**2
          enddo
          rg=sqrt(rg*xnatom)
          return
          end

c######################################################
      subroutine dihedral (x,y,z,dih)
      parameter (Nn=4)
      real*8  x(Nn),y(Nn), z(Nn),dih,pi,cov
      real mconst
      DATA PI/3.141592653589793d0/,cov/57.29578d0/
c
c bond angle distribution
c
c
c dihedral angle
c
         rijx = x(1) - x(2)
         rijy = y(1) - y(2)
         rijz = z(1) - z(2)
         rjkx = x(2) - x(3)
         rjky = y(2) - y(3)
         rjkz = z(2) - z(3)
         rklx = x(3) - x(4)
         rkly = y(3) - y(4)
         rklz = z(3) - z(4)
c
         ax = rijy*rjkz - rijz*rjky
         ay = rijz*rjkx - rijx*rjkz
         az = rijx*rjky - rijy*rjkx
         bx = rjky*rklz - rjkz*rkly
         by = rjkz*rklx - rjkx*rklz
         bz = rjkx*rkly - rjky*rklx
C
C     set to MCONST if smaller than MCONST
         mconst=1.0d-10
         RG=SQRT(MAX(MCONST,RJKX*RJKX+RJKY*RJKY+RJKZ*RJKZ))
c         RGR=1.d0/RG
         RA2R=1.d0/MAX(MCONST,AX*AX+AY*AY+AZ*AZ)
         RB2R=1.d0/MAX(MCONST,BX*BX+BY*BY+BZ*BZ)
         RABR=SQRT(RA2R*RB2R)
CP=COS(PHI)
         CP=RABR*(AX*BX+AY*BY+AZ*BZ)
C SP=SIN(PHI)
C which can be simplify to sin(phi)=|G|H.A/(|A|.|B|)
Cab...B950603 06/29/95
Cab        SP=RG*RABR*(AX*HX+AY*HY+AZ*HZ)
         SP=-RG*RABR*(AX*RKLX+AY*RKLY+AZ*RKLZ)
Cab...
C
         if(cp.gt.1.d0) then
            cp = 1.d0
            theta  = 0.d0
         else if(cp.lt.-1.d0) then

            cp = -1.d0
            theta  = 180.d0
         else             
            theta    = acos(cp)*cov
         endif

         if(sp.lt.0.0) theta = -theta          
         dih = theta
         return
      end
c#######################################
         real*8 function gaussbkn()
         real*8 sg,the,sgtb(100),sg2,pi,fac,fac1,ytab(100),
     &          rn,r
         data init/0/,sgtb/100*0.d0/,sg/5.0/
         DATA PI/3.141592653589793/

          if(init.eq.0) then
               init=1
               sg2=sg*sg
               fac1=sqrt(2.0)*sg
               do i=1,99
               x=i*0.01
               y=2.*x-1.0
               sgtb(i)=fac1*aerf(y)
               enddo              

          endif

          i=int(rn()*99.0)+1
          gaussbkn=sgtb(i)                          
              
c           do i=1,99
c             write(*,*) sgtb(i),i
c           enddo                              
c           stop
         
         return
         end 


         real*8 function gaussside()
         real*8 sg,the,sgtb(100),sg2,pi,fac,fac1,ytab(100),
     &          rn,r
         data init/0/,sgtb/100*0.d0/,sg/10.0/
         DATA PI/3.141592653589793/

          if(init.eq.0) then
               init=1
               sg2=sg*sg
               fac1=sqrt(2.0)*sg
               do i=1,99
               x=i*0.01
               y=2.*x-1.0
               sgtb(i)=fac1*aerf(y)
               enddo              

          endif

          i=int(rn()*99.0)+1
          gaussside=sgtb(i)                          
              
c           do i=1,99
c             write(*,*) sgtb(i),i
c           enddo                              
c             stop
         
         return
         end 

c     #############################3
         function aerf(y)
         real*8 f1,f2,ef,f3
         data ef/1.d-5/

              aerf=0.0
              a1=-50.0
              a2=50.0
              a3=0.0
              i=0
            f1=y-erf(a1)
            f2=y-erf(a2)
 10         f3=y-erf(a3)
            i=i+1

            if(abs(f3).lt.ef.or.i.gt.50000) then
               aerf=a3
c                write(*,*) 'f3',f3,i
               return
            else if((f3*f1).lt.0.) then
               a2=a3
               f2=f3
               a3=0.5*(a1+a3)
            else if((f3*f2).lt.0.) then
               a1=a3
               f1=f3
               a3=0.5*(a3+a2)

            endif 

            goto 10


         return
         end

c     ###########################################
         subroutine rotation(nX,cs,ss,xinp,xout)
         real*8 nX(*),dtphi,xinp(*),xout(*),mx(3),rn1,rm1,
     &   xpara(3),xvert(3),rnx,rmx,rx,rxvert,csphi,ssphi,tmp,
     &   cs,ss
         integer i,j,k

          
         tmp=(nX(1)*xinp(1)+nX(2)*xinp(2)+nX(3)*xinp(3))

         do i=1,3
         xpara(i)=nX(i)*tmp
         xvert(i)=xinp(i)-xpara(i)
         enddo

         mX(1)=nX(2)*xvert(3)-nX(3)*xvert(2)
         mX(2)=nX(3)*xvert(1)-nX(1)*xvert(3)
         mX(3)=nX(1)*xvert(2)-nX(2)*xvert(1)


         do i=1,3
         xout(i)=xpara(i)+cs*xvert(i)+mX(i)*ss
         enddo


         return
         end

c     #######################################3
        subroutine getang(l1,l2,l3,xp,yp,zp,thet)
        real*8 xp(*),yp(*),zp(*),xj(3),xk(3),thet,r1,r2,cov
        parameter(cov=180.0/3.1415926d0)

          xj(1)=xp(l1)-xp(l2)
          xj(2)=yp(l1)-yp(l2)
          xj(3)=zp(l1)-zp(l2)

          xk(1)=xp(l3)-xp(l2)
          xk(2)=yp(l3)-yp(l2)
          xk(3)=zp(l3)-zp(l2)

           thet=0.0
           r1=0.0
           r2=0.0

           do i=1,3
           thet=thet+xj(i)*xk(i)
           r1=r1+xj(i)*xj(i)
           r2=r2+xk(i)*xk(i)
           enddo

          thet=acos(thet/sqrt(r1*r2))*cov

         return
        end


c     #######################################################
        subroutine getdih(l1,l2,l3,l4,xp,yp,zp,phi)
       real*8 xp(*),yp(*),zp(*),xh(4),yh(4)
     & ,zh(4),phi

           xh(1)=xp(l1)
           yh(1)=yp(l1)
           zh(1)=zp(l1)

           xh(2)=xp(l2)
           yh(2)=yp(l2)
           zh(2)=zp(l2)

           xh(3)=xp(l3)
           yh(3)=yp(l3)
           zh(3)=zp(l3)

           xh(4)=xp(l4)
           yh(4)=yp(l4)
           zh(4)=zp(l4)
           call dihedral(xh,yh,zh,phi)

       return
       end
c     ###################################
          
c	###########################################
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This software package is copyrighted, and may be used subject to the following C
C conditions:                                                                    C 
C                                                                                C 
C * This software is provided free of charge to academic users, subject to the   C
C   condition that no part of it be sold or used otherwise for commercial        C
C   purposes, including, but not limited to its incorporation into commercial    C
C   software packages, without written consent from the authors. For permission  C
C   contact Prof. H. A. Scheraga, Cornell University.                            C
C                                                                                C 
C * This software package is provided on an "as is" basis. We in no way warrant  C
C   either this software or results it may produce.                              C
C                                                                                C 
C * Reports or publications using this software package must contain an          C
C   acknowledgement to the authors and the NIH Resource in the form commonly     C 
C   used in academic research.                                                   C 
C                                                                                C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fitsq(rms,x,y,nn,t,b,non_conv)
      implicit real*8 (a-h,o-z)
c      include 'COMMON.IOUNITS'
c  x and y are the vectors of coordinates (dimensioned (3,n)) of the two
c  structures to be superimposed.  nn is 3*n, where n is the number of  
c  points.   t and b are respectively the translation vector and the    
c  rotation matrix that transforms the second set of coordinates to the 
c  frame of the first set.                                              
c  eta =  machine-specific variable                                     
                                                                        
      dimension x(3*nn),y(3*nn),t(3)                                          
      dimension b(3,3),q(3,3),r(3,3),v(3),xav(3),yav(3),e(3),c(3,3)     
      logical non_conv
      eta = z00100000                                                   
c     small=25.0*rmdcon(3)                                              
c     small=25.0*eta                                                    
c     small=25.0*10.e-10                                                
c the following is a very lenient value for 'small'                     
      small = 0.0001D0                                                  
      non_conv=.false.
      fn=nn                                                             
      do 10 i=1,3                                                       
      xav(i)=0.0D0                                                      
      yav(i)=0.0D0                                                      
      do 10 j=1,3                                                       
   10 b(j,i)=0.0D0                                                      
      nc=0                                                              
c                                                                       
      do 30 n=1,nn                                                      
      do 20 i=1,3                                                       
c      write(iout,*)'x = ',x(nc+i),'  y = ',y(nc+i)                           
      xav(i)=xav(i)+x(nc+i)/fn                                          
   20 yav(i)=yav(i)+y(nc+i)/fn                                          
   30 nc=nc+3                                                           
c                                                                       
      do i=1,3
        t(i)=yav(i)-xav(i)
      enddo

      rms=0.0d0
      do n=1,nn
        do i=1,3
          rms=rms+(y(3*(n-1)+i)-x(3*(n-1)+i)-t(i))**2
        enddo
      enddo
      rms=dabs(rms/fn)

c     write(iout,*)'xav = ',(xav(j),j=1,3)                                    
c     write(iout,*)'yav = ',(yav(j),j=1,3)                                    
c     write(iout,*)'t   = ',(t(j),j=1,3)
c     write(iout,*)'rms=',rms
      if (rms.lt.small) then
        b(1,1)=1.0
        b(2,2)=1.0
        b(3,3)=1.0
        return
      endif                                                                  
                                                                        
      nc=0                                                              
      rms=0.0D0                                                         
      do 50 n=1,nn                                                      
      do 40 i=1,3                                                       
      rms=rms+((x(nc+i)-xav(i))**2+(y(nc+i)-yav(i))**2)/fn              
      do 40 j=1,3                                                       
      b(j,i)=b(j,i)+(x(nc+i)-xav(i))*(y(nc+j)-yav(j))/fn                
   40 c(j,i)=b(j,i)                                                     
   50 nc=nc+3                                                           
      call sivade(b,q,r,d,non_conv)
      sn3=dsign(1.0d0,d)                                                   
      do 120 i=1,3                                                      
      do 120 j=1,3                                                      
  120 b(j,i)=-q(j,1)*r(i,1)-q(j,2)*r(i,2)-sn3*q(j,3)*r(i,3)             
      call mvvad(b,xav,yav,t)                                           
      do 130 i=1,3                                                      
      do 130 j=1,3                                                      
      rms=rms+2.0*c(j,i)*b(j,i)                                         
  130 b(j,i)=-b(j,i)                                                    
      if (dabs(rms).gt.small) go to 140                                  
*     write (6,301)                                                     
      return                                                            
  140 if (rms.gt.0.0d0) go to 150                                         
c     write (iout,303) rms                                                 
      rms=0.0d0
*     stop                                                              
c 150 write (iout,302) dsqrt(rms)                                           
  150 continue
      return                                                            
  301 format (5x,'rms deviation negligible')                            
  302 format (5x,'rms deviation ',f14.6)                                
  303 format (//,5x,'negative ms deviation - ',f14.6)                   
      end                                                               
      subroutine sivade(x,q,r,dt,non_conv)
      implicit real*8(a-h,o-z)
c  computes q,e and r such that q(t)xr = diag(e)                        
      dimension x(3,3),q(3,3),r(3,3),e(3)                               
      dimension h(3,3),p(3,3),u(3,3),d(3)                               
      logical non_conv
      eta = z00100000                                                   
      nit = 0
      small=25.0*10.e-10                                                
c     small=25.0*eta                                                    
c     small=2.0*rmdcon(3)                                               
      xnrm=0.0d0                                                          
      do 20 i=1,3                                                       
      do 10 j=1,3                                                       
      xnrm=xnrm+x(j,i)*x(j,i)                                           
      u(j,i)=0.0d0                                                        
      r(j,i)=0.0d0                                                        
   10 h(j,i)=0.0d0                                                        
      u(i,i)=1.0                                                        
   20 r(i,i)=1.0                                                        
      xnrm=dsqrt(xnrm)                                                   
      do 110 n=1,2                                                      
      xmax=0.0d0                                                          
      do 30 j=n,3                                                       
   30 if (dabs(x(j,n)).gt.xmax) xmax=dabs(x(j,n))                         
      a=0.0d0                                                             
      do 40 j=n,3                                                       
      h(j,n)=x(j,n)/xmax                                                
   40 a=a+h(j,n)*h(j,n)                                                 
      a=dsqrt(a)                                                         
      den=a*(a+dabs(h(n,n)))                                             
      d(n)=1.0/den                                                      
      h(n,n)=h(n,n)+dsign(a,h(n,n))                                      
      do 70 i=n,3                                                       
      s=0.0d0                                                             
      do 50 j=n,3                                                       
   50 s=s+h(j,n)*x(j,i)                                                 
      s=d(n)*s                                                          
      do 60 j=n,3                                                       
   60 x(j,i)=x(j,i)-s*h(j,n)                                            
   70 continue                                                          
      if (n.gt.1) go to 110                                             
      xmax=dmax1(dabs(x(1,2)),dabs(x(1,3)))                               
      h(2,3)=x(1,2)/xmax                                                
      h(3,3)=x(1,3)/xmax                                                
      a=dsqrt(h(2,3)*h(2,3)+h(3,3)*h(3,3))                               
      den=a*(a+dabs(h(2,3)))                                             
      d(3)=1.0/den                                                      
      h(2,3)=h(2,3)+sign(a,h(2,3))                                      
      do 100 i=1,3                                                      
      s=0.0d0                                                             
      do 80 j=2,3                                                       
   80 s=s+h(j,3)*x(i,j)                                                 
      s=d(3)*s                                                          
      do 90 j=2,3                                                       
   90 x(i,j)=x(i,j)-s*h(j,3)                                            
  100 continue                                                          
  110 continue                                                          
      do 130 i=1,3                                                      
      do 120 j=1,3                                                      
  120 p(j,i)=-d(1)*h(j,1)*h(i,1)                                        
  130 p(i,i)=1.0+p(i,i)                                                 
      do 140 i=2,3                                                      
      do 140 j=2,3                                                      
      u(j,i)=u(j,i)-d(2)*h(j,2)*h(i,2)                                  
  140 r(j,i)=r(j,i)-d(3)*h(j,3)*h(i,3)                                  
      call mmmul(p,u,q)                                                 
  150 np=1                                                              
      nq=1                                                              
      nit=nit+1
      if (nit.gt.10000) then
c        print '(a)','!!!! Over 10000 iterations in SIVADE!!!!!'
        non_conv=.true.
        return
      endif
      if (dabs(x(2,3)).gt.small*(dabs(x(2,2))+abs(x(3,3)))) go to 160     
      x(2,3)=0.0d0                                                        
      nq=nq+1                                                           
  160 if (dabs(x(1,2)).gt.small*(dabs(x(1,1))+dabs(x(2,2)))) go to 180     
      x(1,2)=0.0d0                                                        
      if (x(2,3).ne.0.0d0) go to 170                                      
      nq=nq+1                                                           
      go to 180                                                         
  170 np=np+1                                                           
  180 if (nq.eq.3) go to 310                                            
      npq=4-np-nq                                                       
      if (np.gt.npq) go to 230                                          
      n0=0                                                              
      do 220 n=np,npq                                                   
      nn=n+np-1                                                         
      if (dabs(x(nn,nn)).gt.small*xnrm) go to 220                        
      x(nn,nn)=0.0d0                                                      
      if (x(nn,nn+1).eq.0.0d0) go to 220                                  
      n0=n0+1                                                           
      go to (190,210,220),nn                                            
  190 do 200 j=2,3                                                      
  200 call givns(x,q,1,j)                                               
      go to 220                                                         
  210 call givns(x,q,2,3)                                               
  220 continue                                                          
      if (n0.ne.0) go to 150                                            
  230 nn=3-nq                                                           
      a=x(nn,nn)*x(nn,nn)                                               
      if (nn.gt.1) a=a+x(nn-1,nn)*x(nn-1,nn)                            
      b=x(nn+1,nn+1)*x(nn+1,nn+1)+x(nn,nn+1)*x(nn,nn+1)                 
      c=x(nn,nn)*x(nn,nn+1)                                             
      dd=0.5*(a-b)                                                      
      xn2=c*c                                                           
      rt=b-xn2/(dd+sign(dsqrt(dd*dd+xn2),dd))                            
      y=x(np,np)*x(np,np)-rt                                            
      z=x(np,np)*x(np,np+1)                                             
      do 300 n=np,nn                                                    
      if (dabs(y).lt.dabs(z)) go to 240                                   
      t=z/y                                                             
      c=1.0/dsqrt(1.0d0+t*t)                                               
      s=c*t                                                             
      go to 250                                                         
  240 t=y/z                                                             
      s=1.0/dsqrt(1.0d0+t*t)                                               
      c=s*t                                                             
  250 do 260 j=1,3                                                      
      v=x(j,n)                                                          
      w=x(j,n+1)                                                        
      x(j,n)=c*v+s*w                                                    
      x(j,n+1)=-s*v+c*w                                                 
      a=r(j,n)                                                          
      b=r(j,n+1)                                                        
      r(j,n)=c*a+s*b                                                    
  260 r(j,n+1)=-s*a+c*b                                                 
      y=x(n,n)                                                          
      z=x(n+1,n)                                                        
      if (dabs(y).lt.dabs(z)) go to 270                                   
      t=z/y                                                             
      c=1.0/dsqrt(1.0+t*t)                                               
      s=c*t                                                             
      go to 280                                                         
  270 t=y/z                                                             
      s=1.0/dsqrt(1.0+t*t)                                               
      c=s*t                                                             
  280 do 290 j=1,3                                                      
      v=x(n,j)                                                          
      w=x(n+1,j)                                                        
      a=q(j,n)                                                          
      b=q(j,n+1)                                                        
      x(n,j)=c*v+s*w                                                    
      x(n+1,j)=-s*v+c*w                                                 
      q(j,n)=c*a+s*b                                                    
  290 q(j,n+1)=-s*a+c*b                                                 
      if (n.ge.nn) go to 300                                            
      y=x(n,n+1)                                                        
      z=x(n,n+2)                                                        
  300 continue                                                          
      go to 150                                                         
  310 do 320 i=1,3                                                      
  320 e(i)=x(i,i)                                                       
      nit=0
  330 n0=0                                                              
      nit=nit+1
      if (nit.gt.10000) then
c        print '(a)','!!!! Over 10000 iterations in SIVADE!!!!!'
        non_conv=.true.
        return
      endif
      do 360 i=1,3                                                      
      if (e(i).ge.0.0d0) go to 350                                        
      e(i)=-e(i)                                                        
      do 340 j=1,3                                                      
  340 q(j,i)=-q(j,i)                                                    
  350 if (i.eq.1) go to 360                                             
      if (dabs(e(i)).lt.dabs(e(i-1))) go to 360                           
      call switch(i,1,q,r,e)                                            
      n0=n0+1                                                           
  360 continue                                                          
      if (n0.ne.0) go to 330                                            
      if (dabs(e(3)).gt.small*xnrm) go to 370                            
      e(3)=0.0d0                                                          
      if (dabs(e(2)).gt.small*xnrm) go to 370                            
      e(2)=0.0d0                                                          
  370 dt=det(q(1,1),q(1,2),q(1,3))*det(r(1,1),r(1,2),r(1,3))            
*     write (1,501) (e(i),i=1,3)                                        
      return                                                            
  501 format (/,5x,'singular values - ',3e15.5)                         
      end                                                               
      subroutine givns(a,b,m,n)                                         
      implicit real*8 (a-h,o-z)
      dimension a(3,3),b(3,3)                                           
      if (dabs(a(m,n)).lt.dabs(a(n,n))) go to 10                          
      t=a(n,n)/a(m,n)                                                   
      s=1.0/dsqrt(1.0+t*t)                                               
      c=s*t                                                             
      go to 20                                                          
   10 t=a(m,n)/a(n,n)                                                   
      c=1.0/dsqrt(1.0+t*t)                                               
      s=c*t                                                             
   20 do 30 j=1,3                                                       
      v=a(m,j)                                                          
      w=a(n,j)                                                          
      x=b(j,m)                                                          
      y=b(j,n)                                                          
      a(m,j)=c*v-s*w                                                    
      a(n,j)=s*v+c*w                                                    
      b(j,m)=c*x-s*y                                                    
   30 b(j,n)=s*x+c*y                                                    
      return                                                            
      end                                                               
      subroutine switch(n,m,u,v,d)                                      
      implicit real*8 (a-h,o-z)
      dimension u(3,3),v(3,3),d(3)                                      
      do 10 i=1,3                                                       
      tem=u(i,n)                                                        
      u(i,n)=u(i,n-1)                                                   
      u(i,n-1)=tem                                                      
      if (m.eq.0) go to 10                                              
      tem=v(i,n)                                                        
      v(i,n)=v(i,n-1)                                                   
      v(i,n-1)=tem                                                      
   10 continue                                                          
      tem=d(n)                                                          
      d(n)=d(n-1)                                                       
      d(n-1)=tem                                                        
      return                                                            
      end                                                               
      subroutine mvvad(b,xav,yav,t)                                     
      implicit real*8 (a-h,o-z)
      dimension b(3,3),xav(3),yav(3),t(3)                               
c     dimension a(3,3),b(3),c(3),d(3)                                   
c     do 10 j=1,3                                                       
c     d(j)=c(j)                                                         
c     do 10 i=1,3                                                       
c  10 d(j)=d(j)+a(j,i)*b(i)                                             
      do 10 j=1,3                                                       
      t(j)=yav(j)                                                       
      do 10 i=1,3                                                       
   10 t(j)=t(j)+b(j,i)*xav(i)                                           
      return                                                            
      end                                                               
      double precision function det (a,b,c)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3),c(3)                                          
      det=a(1)*(b(2)*c(3)-b(3)*c(2))+a(2)*(b(3)*c(1)-b(1)*c(3))         
     1  +a(3)*(b(1)*c(2)-b(2)*c(1))                                     
      return                                                            
      end                                                               
      subroutine mmmul(a,b,c)                                           
      implicit real*8 (a-h,o-z)
      dimension a(3,3),b(3,3),c(3,3)                                    
      do 10 i=1,3                                                       
      do 10 j=1,3                                                       
      c(i,j)=0.0d0                                                        
      do 10 k=1,3                                                       
   10 c(i,j)=c(i,j)+a(i,k)*b(k,j)                                       
      return                                                            
      end                                                               
      subroutine matvec(uvec,tmat,pvec,nback)                           
      implicit real*8 (a-h,o-z)
      real*8 tmat(3,3),uvec(3,nback), pvec(3,nback)                     
c                                                                       
      do 2 j=1,nback                                                    
         do 1 i=1,3                                                     
         uvec(i,j) = 0.0d0                                                
         do 1 k=1,3                                                     
    1    uvec(i,j)=uvec(i,j)+tmat(i,k)*pvec(k,j)                        
    2 continue                                                          
      return                                                            
      end                                                               

