      program random
      implicit integer(i-z)

      no=1124
#      if(no.gt.0)no=-no
      do 11 n=1,10000
         aff=aranzy(no)
         write(*,*) aff
 11   continue
      stop
      end

****************************************************
*     no=-1 (less than 0)                          *
*     rrr=ranzy(no)                                *
****************************************************
      function aranzy(no)
c   random numbers uniformly distributed between 0 and 1.
c   (code by Y. Zhang, ITP, Acdemia Sincica, 1999.6
c   can be replaced by any other suitable generator for random numbers)
      common/asdfa/ix1,ix2,ix3,rm1,rm2,r(99)
      data ia1,ic1,m1/1279,351762,1664557/
      data ia2,ic2,m2/2011,221592,1048583/
      data ia3,ic3,m3/15551,6150,29101/
      if(no.ge.0) go to 2
      ix1=mod(-no,m1)
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ix1,m2)
      ix1=mod(ia1*ix1+ic1,m1)
      ix3=mod(ix1,m3)
      rm1=1./float(m1)
      rm2=1./float(m2)
      do 1 j=1,99
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
 1    r(j)=(float(ix1)+float(ix2)*rm2)*rm1
 2    ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(99*ix3)/m3
      aranzy=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      no=ix1
      return
      end
      
