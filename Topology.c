#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Lyapunov3D.h"

extern int inn[NPMAX][4], degree[NPMAX], cross[NPMAX], crossnn[NPMAX][4], qcore[NPMAX];

extern int qindex[NPMAX][4], qmax[4];

typedef struct{int i, j;} vector_int;
typedef struct{vector_int min,max;} rectangle;

void Topology( int lp, int mp, int sp)
{
  int i,j,k, dir;
  int q,p;

  int cont[4];

  int  ncore;

  int east, north, west, south;

  vector_int incr[4];
  rectangle domain[4];
  int  i_old, j_old;


  /* Inicializing Variables:*/

  incr[0].i = 1;
  incr[0].j = 0;
  for(dir=1; dir<4; dir++)
    {
      incr[dir].i = -incr[dir-1].j;
      incr[dir].j =  incr[dir-1].i;
    }

  for(dir=0; dir<4; dir++)
    {
      domain[dir].min.i = -1*(incr[dir].i<0)*incr[dir].i ;
      domain[dir].min.j = -1*(incr[dir].j<0)*incr[dir].j ;
      domain[dir].max.i = lp - (incr[dir].i>0)*incr[dir].i ;
      domain[dir].max.j = mp - (incr[dir].j>0)*incr[dir].j ;
    }

  cont[0] = cont[1] = cont[2] = cont[3] = 0;
  for(k=0; k < sp; k++)
    {
      for(j = 0; j < mp; j++)
	{
	  for(i = 0; i < lp; i++)
	    {
	      q = i + j * lp + k * lp * mp;
	      degree[q] = 0;	      
	      cross[q] = 0;
	      for(dir=0; dir<4; dir++)
		{
		  crossnn[q][dir] = 0;
		  if(i>=domain[dir].min.i && i<domain[dir].max.i && j>=domain[dir].min.j && j<domain[dir].max.j)
		    {
		      inn[q][dir] = q + incr[dir].i + incr[dir].j * lp;
		      qindex[cont[dir]][dir] = q;
		      cont[dir]++;
		      qmax[dir] = cont[dir];
		    }
		  else
		    inn[q][dir] = -1;
		}
	    }
	}
    }

  p = 0;
  for(k=0; k < sp; k++)
    {
      for(j = 1; j < (mp-1); j++)
	{
	  for(i = 1; i < (lp-1); i++)
	    {
	      q = i + j * lp + k * lp * mp;
	      qcore[p]=q;
	      p++;
	    }
	}
    }
  ncore = (lp - 2) * (mp - 2) * sp;
 
  for(p = 0; p < ncore; p++)
    {
      q = qcore[p];
      CreateCross(q);
    }
}

void CreateCross(int q)
{
  int dir;

  degree[q] += 4;
  cross[q] = 1;
  for(dir = 0; dir < 4; dir++)
    {	  
      degree[inn[q][dir]]++;
      crossnn[inn[q][dir]][(dir+2) % 4] = 1;
    }
}
void DestroyCross(int q)
{
  int dir;

  degree[q] -= 4;
  cross[q] = 0;
  for(dir = 0; dir < 4; dir++)
    {	  
      degree[inn[q][dir]]--;
      crossnn[inn[q][dir]][(dir+2) % 4] = 0;
    }
}
