/*Input template coordinates from chain.dat*/
#ifndef _CHAIN_2_LISTPDBATOM_
#define _CHAIN_2_LISTPDBATOM_

using namespace std;
#include "pdbClasses2.h"
#include<iostream>
#include<fstream>
#include<string>

#define _INFc2l_ (10e10)

bool read_chain_line(ifstream &, string &) ;
int read_chain_file(const char *, listsPDBatom &);

#endif
