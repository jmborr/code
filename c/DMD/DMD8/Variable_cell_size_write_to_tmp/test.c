#include <stdio.h>
#include "cell_size.h"

int main(){
  int i;
  double range=7.5;
  dimensions bound[NDIM];

  for(i=0;i<NDIM;i++)
    {
      bound[i].length=100;
      bound[i].dl=1.501;
      bound[i].period=100;
    }

  /*printf("dl=%f, period=%d \n",bound[0].dl,bound[0].period);*/
  set_amso(range,bound);
  allocate_shell_order_scheme(bound);



  return 0;
}
