#ifndef _PROSPECTOR_2_LISTPDBATOM_
#define _PROSPECTOR_2_LISTPDBATOM_

using namespace std;
#include "pdbClasses2.h"
#include<iostream>
#include<fstream>
#include<string>

bool read_prosp_line(ifstream &, string &);
int read_prosp_file(const char *, listsPDBatom &);

#endif
