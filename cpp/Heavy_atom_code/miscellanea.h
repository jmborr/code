#ifndef _MISC_
#define _MISC_

string char2string( const char * ) ;
void string2char( string &, char * ) ;
bool go_to( ifstream &, const string & ) ;
int remove_firsts( const string &, string & ) ;
int remove_lasts( const string &, string & ) ;
string* split( const string &, int & ) ;
double *split_line2doubles( const string &, int & ) ;
int *split_line2ints( const string &, int & ) ;
#endif/*_MISC_*/
