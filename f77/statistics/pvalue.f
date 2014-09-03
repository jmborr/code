c     MORE STUFF IN f77/miscelanea

c     #############################3
c     Inverse of the error function
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
