#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <netcdf.h>
#include "Lyapunov3D.h"

extern point ptm[NPMAX];

int LinearInterpolation(double t, int ipoint, point pt, vector *vint)
{
  point proms[16];
  vector vroms[16];

  int index;
  int tau, i, j, k;
  double alpha, beta;

  int t0;

  /* Locate Box: Calculates the vectors vroms[15] and points proms[15] */
    
  if(LocateBox(t, ipoint, pt, proms, vroms))
    return 1;
  
  /* Depth Interpolation: */

  for(index = 0; index < 16; index=index+2)
    {
      alpha = (proms[index+1].dpt - pt.dpt)/(proms[index+1].dpt - proms[index].dpt);
      beta = 1 - alpha;
      INTERPOLATION(vroms[index/2], vroms, index, alpha, beta);
      proms[index/2] = proms[index]; 
    }

  /* Mu Interpolation: */
  alpha = (proms[1].mu - pt.mu)/(proms[1].mu - proms[0].mu);
  beta = 1.0 - alpha;
  for(index = 0; index < 8; index=index+2)
    {
      INTERPOLATION(vroms[index/2], vroms, index, alpha, beta);
      proms[index/2] = proms[index];
    }

  /* Lon Interpolation */
  alpha = (proms[1].lon - pt.lon)/(proms[1].lon - proms[0].lon);
  beta = 1.0 - alpha;
  for(index = 0; index < 4; index=index+2)
    {
      INTERPOLATION(vroms[index/2], vroms, index, alpha, beta);
      proms[index/2] = proms[index];
    }

  /* Time Interpolation: */
  t0 = (int) (t + 1); 
  alpha = ((double) t0) - t;  
  beta = 1.0 - alpha; 
  INTERPOLATION(vroms[0], vroms, 0, alpha, beta);

  *vint = vroms[0];

  return 0;
}
