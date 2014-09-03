#include<stdio.h>
#include<stdlib.h>
#include <time.h>
#include "seldon.h"

extern char tmp_dir[13];

void get_tmp_dir()
{
  time_t seconds;
  char command[40];
  seconds= time(NULL);
  sprintf(tmp_dir,"/tmp/%06d/", seconds%1000000);
  sprintf(command,"mkdir %s", tmp_dir);
  system(command);
  return;
}
