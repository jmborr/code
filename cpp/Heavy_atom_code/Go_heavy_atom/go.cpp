#include "go.h"

/*======================================================================*/
double d2( double *u, double *v ){
  return (u[0]-v[0])*(u[0]-v[0]) + (u[1]-v[1])*(u[1]-v[1]) + 
    (u[2]-v[2])*(u[2]-v[2]) ;
}
/*======================================================================*/
/*backbone atoms are not considered Go atoms*/
void output_nat_cont( PDBchain &native, listPDBatom &nat_p, listPDBatom &nat_q,
		      const int &min_sep, const double &go_f, ifstream &DAT,
		      int *&p, int *&q, double *&go_ranges, int &nnc ){

  /*clean arrays and variables to be filled later*/
  nnc = 0 ;
  if( p ){ delete [ ] p ; } 
  if( q ){ delete [ ] q ; } 
  if( go_ranges ){ delete [ ] go_ranges ; }

  /*check input*/
  /*cout<<nat_p<<"TER\n"<<nat_q<<"min_sep="<<min_sep<<"  go_f="<<go_f ;
  cout<<"  p="<<p<<"  q="<<q<<"  go_ranges="<<go_ranges<<"  nnc="<<nnc<<endl;/*

  /*prepare list of hardcore radii, prepare list of residue sequence numbers*/
  int lp=nat_p.length( ),lq=nat_q.length( ); 
  /*cout<<"lp="<<lp<<" lq="<<lq<<endl;*/
  PDBatom *atom_p, *atom_q ;
  int *resSeq_p = new int[ lp+1 ], *resSeq_q = new int[ lq+1 ] ;
  string amino_and_atom ;
  double *hcr_list_p = new double[ lp+1 ], *hcr_list_q = new double[ lq+1 ] ;
  for( int i=1; i<=lp; i++ ){
    atom_p = nat_p.getAtomAt( i ) ;
    resSeq_p[ i ] = atom_p->getResSeq( ) ;
    amino_and_atom = atom_p->get_resName_nameII( ) ;
    hcr_list_p[ i ] = out_hard_core_radius( DAT, amino_and_atom ) ;
    /*cout<<*atom_p<<endl<<hcr_list_p[ i ]<<endl;*/
    delete atom_p ;
  }
  for( int i=1; i<=lq; i++ ){
    atom_q = nat_q.getAtomAt( i ) ;
    resSeq_q[ i ] = atom_q->getResSeq( ) ;
    amino_and_atom = atom_q->get_resName_nameII( ) ;
    hcr_list_q[ i ] = out_hard_core_radius( DAT, amino_and_atom ) ;
    /*cout<<*atom_p<<endl<<hcr_list_p[ i ]<<endl;*/
    delete atom_q ;
  }

  /*create temporary array*/
  int max, max_nc = 100 ;
  lp>lq ? max=lp : max= lq ;  max *= max_nc ;
  int *tmp_p = new int[max], *tmp_q = new int[max] ;
  double *tmp_go_ranges = new double[max] ;

  /*search for native contacts*/
  int a, b ;
  double hc, hc1, hc2 ;
  PDBatom atom_pII,atom_qII;
  for( int i=1; i<=lp; i++ ){
    atom_p = nat_p.getAtomAt( i ) ;
    if( !atom_p->is_backbone( ) ){
      a = atom_p->getSerialNumber( ) ;
      atom_pII=native.getAtomWithIndex(a);
	hc1 = hcr_list_p[ i ] ;
      for( int j=1; j<=lq; j++ ){
	atom_q = nat_q.getAtomAt( j ) ;
	if( ( !atom_q->is_backbone( ) )           &&
	    ( abs(resSeq_p[i]-resSeq_q[j]) > min_sep) ){
	  hc2 = hcr_list_q[ j ] ; 
	  hc = (1+go_f )*(hc1+hc2) ;
	  b = atom_q->getSerialNumber( ) ;
	  atom_qII=native.getAtomWithIndex(b);
	  if( atom_pII.d2( atom_qII ) < hc*hc ){
	    tmp_p[nnc]=a ;  tmp_q[nnc]=b ; tmp_go_ranges[nnc]=hc ;  nnc++ ;
	  }
	}
	delete atom_q ;
      }    
    }
    delete atom_p ;
  }

  /*write to "p" and "q" from the temporary arrays*/
  p = new int[nnc] ; q = new int[nnc] ; go_ranges = new double[nnc] ;
  for( int i=0; i<nnc; i++){ 
    p[ i ] = tmp_p[ i ] ; q[ i ] = tmp_q[ i ] ; 
    go_ranges[ i ] = tmp_go_ranges[ i ] ;
  }
}
/*======================================================================*/
void shift_to_center( double **r, const int &natoms, const double &boundary ){
  double bH = boundary/2.0 ;
  double dr[3], old_r[3]={r[0][0], r[0][1], r[0][2] } ;
  for( int i=1; i<natoms; i++ ){/*remove periodic boundary conditions*/
    for( int j=0; j<3; j++ ){
      dr[j] = r[i][j]-old_r[j] ;
      if( dr[j] > bH ){ dr[j] -= boundary ; }
      else if( dr[j] < -bH ){ dr[j] += boundary ; }
      old_r[j] = r[i][j] ;
      r[i][j] = r[i][j] + dr[j] ;
    }
  }
  double gc[3]={ 0.0, 0.0, 0.0 } ;/*geometric center*/
  for(int i=0;i<natoms;i++){for(int j=0;j<3;j++){gc[j]+=r[i][j];}}
  for(int j=0;j<3;j++){ gc[j] /= natoms ; }
  for(int i=0;i<natoms;i++){for(int j=0;j<3;j++){r[i][j]+=bH-gc[j];}}
}
/*=================================================================*/
/*given a list of native contacts (p[i],q[i]), each contact with range
  go_range[i], the list has n_nc native contacts, given a list of 
  coordinates r. The program outputs the fraction of native contacts present.
  NOTE: example, for nat cont (3, 56 ) we would look at r[2][xyz] and
  r[55][xyz], because we assume atoms are shifted by one*/
double output_fraction( PDBchain &native, double **r, int *p, int *q,
		 double *go_ranges, const int &n_nc ){
  int n=0 ;
  double range2 ;
  double *s=new double[3];
  for( int i=0; i<n_nc; i++ ){
    native.getCoordOfAtomWithAtomIndex( p[i], s );
    /*cout<<"p="<<p[i]<<" p_rnat="<<s[0]<<" "<<s[1]<<" "<<s[2]<<" p_r="<<r[p[i]-1][0]<<" "<<r[p[i]-1][1]<<" "<<r[p[i]-1][2]<<" || ";
    native.getCoordOfAtomWithAtomIndex( q[i], s );
    cout<<"q="<<q[i]<<" q_rnat="<<s[0]<<" "<<s[1]<<" "<<s[2]<<" q_r="<<r[q[i]-1][0]<<" "<<r[q[i]-1][1]<<" "<<r[q[i]-1][2]<<"  ";
    cout<<"go_ranges="<<go_ranges[i]<<endl;*/
    range2 = go_ranges[i] * go_ranges[i] ;
    if( d2( r[ p[i]-1 ], r[ q[i]-1 ] ) < range2*1.01 ){ n++ ; }
  }
  /*cout<<"n="<<n<<" n_nc="<<n_nc<<endl;*/
  return double(n)/n_nc ;
}
/*=================================================================*/
