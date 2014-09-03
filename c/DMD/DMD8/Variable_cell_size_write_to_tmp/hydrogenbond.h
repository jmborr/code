#ifndef _HB_
#define _HB_
#include "bcp.h"

/* I have considerably simplified the hydrogen bond reaction. If (p,q)
undergo forming/breaking of hydrogen bond, then I only change
collision types for the (p,q) collision and for the four auxiliary
collisions. I do not change the collision type for other (p,r) and
(q,r) collisions. For example, say atom p is a nitrogen, type=1, that
undergoes a reaction with q, a carbon of type=2. Then atom type of
 nitrogen p changes to type=3 and atom type of carbon q changes to
 type=4. Imagine that p was bonded to its CA (type 5). Before (p,q)
 reaction, the link was between an atom of type 1 and an atom of
 type 5. Thus it corresponds some collision type "ct". After the
 reaction, the link is between an atom of type 3 and atom of
 type 5. Thus it corresponds some collision type ct'. We do not
 change ct to ct' because the specifications (distance, change in
 energy, and so on) are the same for ct than for ct'. In fact, the
 whole potential that describes a bond between one atom of type 1
 and one atom of type 5 is the same than the potential describing 
 a bond between one atom of type 3 and one atom of type 5. Thus
 there's no need to change collision types. The same goes for any
 other atom interacting with p, except for q and for the two
 associations of q that are now bonded to p.
*/

void set_amso_for_HB(int value);
int isHBRelated(atom*, atom*);
coll_type process_hbII(atom *a,atom_number i1,atom_number i2,coll_type ct1,
		     double sc,double x,double y,double z,coll_type ct,
		     coll_type rtype, int revers);
coll_type process_hbII2(atom *a,atom_number i1,atom_number i2,coll_type ct1,
		     double sc,double x,double y,double z,coll_type ct,
		     coll_type rtype, int revers);
int checkAssociatedBond(atom*, int, atom*, int);
int checkDessociatedBond(atom*, int, atom*, int);
#endif /*_HB_*/
