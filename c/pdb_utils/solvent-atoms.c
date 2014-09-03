/* Read: mass [list or residues]
   Read: coordinates, velocities

   Calculate a list of water molecules which are within
      the first solvation shell for the protein [list of residues]
     - 1 to N protein residues
           Any water O within 3 \AA of heavy atom
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main(int argc, char** argv)
{ 
  int i,j,k,m,n,nwat;
  double *x, *y, *z;
  double xx,yy,zz;
  double *mass;
  int *resi;
  double dist;
  
  // n = number of atoms in system 
  n=atoi(argv[2]);
  mass = malloc(n * sizeof(double));
  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));
  z = malloc(n * sizeof(double));

  // m = number of protein atoms
  m=atoi(argv[4]);
  resi = malloc(m * sizeof(int));
  for(i = 1; i <= m; i++)
	  scanf("%d", &resi[i]);
  
  for (i=1; i<=n; i++)
  {
	  scanf("%lf",&mass[i]);
  }

  for (i=1; i<=n; i++)
  {
	  scanf("%lf %lf %lf",&x[i],&y[i],&z[i]);
  }
  
   for (i=1; i<=m; i++)
   {
        if (mass[i] > 1.10e+00)
        {
		for (j=m+1; j<=n; j++)
                {
			if (mass[j] > 1.10e+00)
                    	{
				xx=x[i]-x[j];
                     		yy=y[i]-y[j];
                     		zz=z[i]-z[j];
                     		dist=xx*xx + yy*yy + zz*zz;
                     		dist= sqrt(dist);
                     		if (dist <= 3.00e+00)
                        	{
					printf("%10d%10d\n", i, resi[i]);
                        	}
                    	}
                 }
            }
       }
}
