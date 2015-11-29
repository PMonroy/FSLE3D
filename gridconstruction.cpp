#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#include "vectorXYZ.h"
#include "constants.h"

int gridfsle3d(vector<vectorXYZ> *itracer, int *ni, int *nj, int *nk, vectorXYZ domainmin, vectorXYZ intergrid, vectorXYZ domainmax)
{
  double x,y,z;
  int i,j,k;
  int ntracers;
  double mumin,mumax;

  mumin = log((1.0/cos(rads*domainmin.y))+tan(rads*domainmin.y))*degrees;
  mumax = log((1.0/cos(rads*domainmax.y))+tan(rads*domainmax.y))*degrees;

  ntracers = (int)((domainmax.z-domainmin.z)/intergrid.z)+1;
  ntracers *= (int)((mumax-mumin)/intergrid.y)+1;
  ntracers *= (int)((domainmax.x-domainmin.x)/intergrid.x)+1;

  if(ntracers<=0)
    return 1;
  
  (*itracer).reserve(ntracers);
  for(z=domainmin.z, k=0; z<=domainmax.z; z+=intergrid.z, k++)
    {
      for(y=domainmin.y,j=0; y<=domainmax.y; y+=intergrid.x*cos(rads*y),j++)
	{
	  for(x=domainmin.x,i=0; x<=domainmax.x; x+=intergrid.x,i++)
	    {
	      (*itracer).push_back(vectorXYZ(x,y,z));
	    }
	}
    }

  *ni=i;
  *nj=j;
  *nk=k;

  return 0;
}

vector<int> neighbors(int ni, int nj, int nk)
{

  int i;
  int j;
  int k;
  int ntracers;
  int nlayer;

  int q;
  int qmin[6],qincr[6],qmax[6];
  int qneighbor;
  int dir;

  vector<int> neighbor;

  nlayer=ni*nj;
  ntracers=nk*nlayer;

  neighbor.reserve(6*ntracers);

  qincr[0]=1;
  qincr[1]=ni;
  qincr[2]=nlayer;
  qincr[3]=-1;
  qincr[4]=-ni;
  qincr[5]=-nlayer;

  qmin[2]=qmin[5]=0;
  qmax[2]=qmax[5]=ntracers;
 
  for(k=0; k<nk; k++)
    {
      qmin[1]=qmin[4]=k*nlayer;
      qmax[1]=qmax[4]=(k+1)*nlayer;      
      for(j=0; j<nj; j++)
	{
	  qmin[0]=qmin[3]=j*ni+qmin[1];
	  qmax[0]=qmax[3]=(j+1)*ni+qmin[1];
	  for(i=0; i<ni; i++)
	    {
	      q=i+j*ni+k*nlayer;
	      for(dir=0; dir<6; dir++)
		{
		  qneighbor=q+qincr[dir];
		  if(qneighbor>=qmin[dir] && qneighbor<qmax[dir])
		    neighbor.push_back(qneighbor);
		  else
		    neighbor.push_back(-1);
		}
	    }
	}
    }

  return neighbor;
}
