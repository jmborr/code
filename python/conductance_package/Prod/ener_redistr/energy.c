/* gcc -O3 -o energy.x energy.c -lm
   This program scales velocities of a sub-set (active-site) atoms 
   and redistributes them into water atom 

   Read: Natoms
   Read: Scaling factor (could be in deg K in ref to 300K)
   Read: An array of Natoms, where
        1 = active-site atoms
        2 = water oxygen atoms
        3 = water hydrogen atoms
        0 = otherwise
   (Note: The 2 water hydrogen atoms should follow oxygen atom)
   Read: mass
   Read: coordinates, velocities

   Put the difference of energy (1/2mv^2) into a bin variable

   Redistribute the energy from the bin equally into the water atoms

   Check if total KE is conserved

   Punch out coordinates with velocities
    
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

main (int argc, char **argv)
{ int i,j,count,count2,n,list[60000];
  double x[60000],y[60000],z[60000];
  double vx[60000],vy[60000],vz[60000],mass[60000];
  double factor,factorw,kex,key,kez; 
  double kdelta,ktemp, ktemp2, ewat, eperwater, bin;

  int count1;
  double Eold, Enew;
  if(argc!=2){
    printf("Usage: energy.x delta-K(Kcal/mol)\n");
    return 1;
  }

   count = 0;
   kex   = 0.000;
   key   = 0.000;
   kez   = 0.000;

   kdelta=atof(argv[1])*1.26 ; //(Kcal/mol) to Amber-Energy-Units
   //printf("%s argv[1]=%lf delta-K=%lf\n",argv[1],atof(argv[1]),kdelta);

   scanf(" %d",&n); //number of atoms in system
   //printf("Atoms %d\n",n);

   for (i=1; i<=n; i++)
       {scanf(" %d",&list[i]);
//       printf(" %d\n",list[i]);
       }

//   printf("Mass\n");
   for (i=1; i<=n; i++)
       {scanf("%lf",&mass[i]);
//        printf("%12.5f\n",mass[i]);
       }
   for (i=1; i<=n; i++)
       {scanf("%lf %lf %lf",&x[i],&y[i],&z[i]);
       }
   for (i=1; i<=n; i++)
       {scanf("%lf %lf %lf",&vx[i],&vy[i],&vz[i]);
       kex=kex + mass[i]*vx[i]*vx[i];
       key=key + mass[i]*vy[i]*vy[i];
       kez=kez + mass[i]*vz[i]*vz[i];
       }
    kex=0.500e+00*kex;
    key=0.500e+00*key;
    kez=0.500e+00*kez;
    Eold=kex+key+kez;
    //printf("Check %12.5f  %12.5f  %12.5f %12.5f\n",kex,key,kez,Eold);

    count1=0;
    ktemp=0.0; //calculate initial energy of hot-spot
    for (i=1; i<=n; i++)
      if (list[i] == 1){
	count1+=1;
	ktemp+=0.5000*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
      }
   

    factor=(ktemp-kdelta)/ktemp; //rescale velocities of hot-spot
    if (factor<0) factor=0.0;     //better be safe!
    factor=sqrt(factor);
    for (i=1; i<=n; i++)
      if (list[i] == 1){
	vx[i] = factor*vx[i]; 
	vy[i] = factor*vy[i]; 
	vz[i] = factor*vz[i]; 
      }

    ktemp=0.0; //calculate initial energy of first 24 water molecules
    for (i=1; i<=n; i++){
      if (list[i] == 2 || list[i] == 3){
	ktemp+=0.5000*mass[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
      }
    }
    
    factor=sqrt((ktemp+kdelta)/ktemp); //rescale velocities of water atoms
    //printf("water factor=%lf\n",factor);
    for (i=1; i<=n; i++){
      if (list[i] == 2 || list[i] == 3){
	vx[i] = factor*vx[i]; 
	vy[i] = factor*vy[i]; 
	vz[i] = factor*vz[i]; 
      }
    }

  count2 = 0 ; 
  for (i=1; i<=n; i++)
      {
       printf("%12.7f%12.7f%12.7f",x[i],y[i],z[i]);
       count2 = count2 + 1;
       if (count2 == 2)
          {count2 = 0;
           printf("\n");
          }
      }
  for (i=1; i<=n; i++)
      {
       printf("%12.7f%12.7f%12.7f",vx[i],vy[i],vz[i]);
       count2 = count2 + 1;
       if (count2 == 2)
          {count2 = 0;
           printf("\n");
          }
      }
       if (count2 == 1)
          {printf("\n");
          }

   kex   = 0.000;
   key   = 0.000;
   kez   = 0.000;

   for (i=1; i<=n; i++)
       {
       kex=kex + mass[i]*vx[i]*vx[i];
       key=key + mass[i]*vy[i]*vy[i];
       kez=kez + mass[i]*vz[i]*vz[i];
       }
    kex=0.500e+00*kex;
    key=0.500e+00*key;
    kez=0.500e+00*kez;
    Enew=kex+key+kez;
    //printf("Check %12.5f  %12.5f  %12.5f %12.5f\n",kex,key,kez,Enew-Eold);
}
