/*======================================================================*/
/*               INT TYPE MULTIDIMENSIONAL ARRAYS                    */
/*======================================================================*/
int *allocI1(int);
void dumpI1(int*,int,int,int);
int **allocI2(int,int);
void deallocI2(int**);
void dumpI2(int**,int,int,int);
int ***allocI3(int,int,int);
void deallocI3(int***);
int **allocsymI2(int);
void deallocsymI2(int**);
void dumpsymI2(int**,int,int);

/*======================================================================*/
/*               DOUBLE TYPE MULTIDIMENSIONAL ARRAYS                    */
/*======================================================================*/
double *allocD1(int);
double *dumpD1(double*,int,int,int);
double **allocD2(int,int);
void deallocD2(double**);
void dumpD2(double**,int,int,int,int);
double ***allocD3(int,int,int);
void deallocD3(double***);
/*======================================================================*/
/*               CHAR TYPE MULTIDIMENSIONAL ARRAYS                    */
/*======================================================================*/
char *allocC1(int);
char **allocC2(int,int);
void deallocC2(char**);
char ***allocC3(int,int,int);
void deallocC3(char***);
char **allocsymC2(int);/*symmetric matrix. Only (i,j) with i<j is filled*/
void deallocsymC2(char**);
