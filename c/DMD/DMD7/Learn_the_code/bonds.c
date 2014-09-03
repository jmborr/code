#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "bonds.h"

                           /*Remember that typedef int bond_type;            */
static bond_type *friend;  /*                                                */
static bond_type *next;    /*						     */
static bond_type *freeInd; /*						     */
static bond_type *first;   /*						     */
static bond_type nAtom;    /*Number of atoms                                 */
static bond_type nFriend;  /*Maximum allowable number of links               */
static bond_type nTotal;   /*nFriend+nAtom                                   */
static bond_type freeCount;/*                                                */
static bond_type newBonds; /*                                                */
static bond_type newTypes; /*fix it 6                                        */


/*n=number of atoms; nrt=number of reaction types; lbt=length of bond table*/
int allocBonds(int n, int numbond,int nrt,int lbt)
{
  int i;
  nAtom=n;
  if(!nrt)
    nFriend=2*numbond;
  else
    {
      nFriend=NBONDS*nAtom;/*NBONDS=12 by default                            */
      if(lbt)              /*                                                */
	nFriend=lbt;       /*maximum number of bonds passed in the input     */
    }                      /*configuration file                              */
  nTotal=nFriend+nAtom;
  next =(bond_type *)malloc((nTotal)*sizeof(bond_type));
  if(!(next))return 0; 
  first=next+nFriend;
  for(i = 0; i <nAtom ; i++)
    first[i] = -1;
  freeCount=nFriend;
  if(!freeCount)return 1;
  friend =(bond_type *)malloc(nFriend*sizeof(bond_type));
  if(!(friend))return 0;
  freeInd =(bond_type *)malloc(nFriend*sizeof(bond_type));
  if(!(freeInd))return 0;

  for(i=0; i <freeCount ; i++)
    {
      freeInd[i]=i;
      friend[i]=-1;
    }
  return 1;
}


/*creates some kind of link list with the help of arrays next, friend and    */
/*freeInd                                                                    */
int setBond(int friend1, int friend2)
{ 
  int i;
  bond_type index;
  bond_type newindex;
  bond_type newfriend;
  int found;
  int res=0;
  if (friend1==friend2)return 0;
  for(i=0;i<2;i++)
    {
      if(i)                        /*we count twice the bond,i.e.,we set a*/
	{                          /*bond friend1--friend2 and we set a   */
	  index=friend1+nFriend;   /*bond friend2--friend1                */
	  newfriend=friend2;
	}
      else
	{
	  index=friend2+nFriend;
	  newfriend=friend1;
	}
      found=0;
      while(next[index]>=0) /*if we had previously allocated a bond in for   */
	{                   /*"index", then go through the list of bonded    */
	  index=next[index];/*friends until we finish list (next[index]==-1) */
	  if(friend[index]==newfriend){/*we can stop search through the list */
	    found=1;break;             /*if we actually find the friend we   */
	  }                            /*want to allocate ("newfriend")      */
	}
      if(!found)
	{
	  if(!(freeCount))return -1;  /*freeCount is nFriend minus the twice */
	  (freeCount)--;              /*the number of allocated bonds        */
	  newindex=freeInd[freeCount];/*(remember we allocate twice each bond*/
	  next[index]=newindex;       /*and also that nFriend is twice the   */
	  next[newindex]=-1;          /*total number of bonds                */
	  friend[newindex]=newfriend;
	  res++;
	}
    }
  if(res>0)newBonds=1;
  return res;
}/*Matches int setBond(int friend1, int friend2)*/

int breakBond(int friend1, int friend2)/*39 lines of code*/
{ 
  int i;
  bond_type index;
  bond_type oldindex;
  bond_type oldfriend;
  int found;
  int res=0;
  for(i=0;i<2;i++)
    {
      if(i)
	{
	  index=friend1+nFriend;
	  oldfriend=friend2;
	}
      else
	{
	  index=friend2+nFriend;
	  oldfriend=friend1;
	}
      found=0;
      while(next[index]>=0)
	{
	  oldindex=next[index];
	  if(friend[oldindex]==oldfriend)
	    {
	      res++;
              next[index]=next[oldindex];
              next[oldindex]=-1;
              friend[oldindex]=-1;
              freeInd[freeCount]=oldindex;
              freeCount++;
	      break;
	    }
	  index=oldindex;
	}
    }
  if(res>0)newBonds=1;
  return res;
}/*Matches int breakBond(int friend1, int friend2)*/

	
void setNewBonds(int value)
{
  newBonds = value;
}

void setNewTypes(int value)
{
  newTypes = value;
}

int getNewBonds(void)
{
  return (int)newBonds;
}

int getNewTypes(void)
{
  return (int)newTypes;
}
    
int isFriend(bond_type atomNumber, bond_type friendNumber)
{
  bond_type index=nFriend+atomNumber;
  if(atomNumber == friendNumber)
    return 1;
  while(next[index]>=0)
    {
      index=next[index];
      if(friend[index]==friendNumber)return 1;
    }
  return 0;
}

bond_type nextFriend(int atomNumber, int * index)
{
  bond_type newindex;
  if(index[0]==-1)
    newindex=nFriend+atomNumber;
  else
    newindex=index[0];
  if((newindex<0)||(newindex>=nTotal)){index[0]=-1;return -1;}
  
  newindex=next[newindex];
  if(newindex<0)
    {index[0]=-1;return -1;}
  else
    {index[0]=newindex;return friend[newindex];}
}

int getMaxBonds(void)
{return (int)(nFriend>>1);}






