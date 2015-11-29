/*
  Program for calculating lagrangian trajectories:

  The data is from ROMS benguela.

  P. Monroy , IFISC, July 2014
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Lyapunov3D.h"

point ptm[NPMAX];

int inn[NPMAX][4], degree[NPMAX], cross[NPMAX], crossnn[NPMAX][4], qcore[NPMAX];

int qindex[NPMAX][4], qmax[4];

int main()
{
 
  static point ptgrid[NPMAX];
  static double delta[NPMAX][4], deltaf[NPMAX];
  static double tau[NPMAX];
  double deltamax;
  double delta0_rads, mu_p_min, mu_p_max;
  
  int lp, mp, sp, np; 
  int ipoint;
 
  int i, j, k;

  char nameout[256]; 
  FILE *output;
 
  double t, tmin, tmax;

  double phi0, phi1, theta0, theta1;
  double cos_Phi;

  /* TEST VARIABLES*/
  int ncore,q;
  int cont;
  int dir,p;
  /* START VALUE : */

  mu_p_min = LAT_TO_MU(LAT_P_MIN);
  mu_p_max = LAT_TO_MU(LAT_P_MAX);  
  delta0_rads = RADS * DELTA0;


  lp = ((int) ((LON_P_MAX-LON_P_MIN)/DELTA0) + 3);       // max index xi of grid points
  mp = ((int) ((mu_p_max-mu_p_min)/delta0_rads) + 3);    // max index eta of grid points
  sp = ((int) ((DPT_P_MAX-DPT_P_MIN)/DELTA0_DEPTH) + 1); //max index s of material points
  np = lp * mp * sp;

  printf("lp = %d\n", lp);
  printf("mp = %d\n", mp);
  printf("sp = %d\n", sp);
  printf("np = %d\n", np);

  if(np>NPMAX)
    {
      printf("I had stopped because np>NPMAX, sorry.\n");
      return 1;
    }

  sprintf(nameout,"./grid.dat");
  output = fopen(nameout,"w");
 
  ipoint = 0; 
  for(k=0; k < sp; k++) 
    {
      for(j = 0; j < mp; j++)
	{
	  for(i = 0; i < lp; i++)
	    {
	      ptgrid[ipoint].lon = i * DELTA0 + LON_P_MIN - DELTA0 ;

	      ptgrid[ipoint].mu =  j * delta0_rads + mu_p_min - delta0_rads;
	      ptgrid[ipoint].lat = MU_TO_LAT(ptgrid[ipoint].mu);

	      ptgrid[ipoint].dpt = k * DELTA0_DEPTH + DPT_P_MIN;
	      ptgrid[ipoint].r = R_EARTH + ptgrid[ipoint].dpt;

	      ptm[ipoint] = ptgrid[ipoint];
	      fprintf(output,"%lf %lf %lf\n", ptgrid[ipoint].lon, ptgrid[ipoint].lat, ptgrid[ipoint].dpt);
	      tau[ipoint] = DURATION;
	      for(dir=0; dir<4; dir++)
		{
		  delta[ipoint][dir]=0.0;
		}
	      ipoint++;
	    }
	  fprintf(output,"\n");
	}
    }

  fclose(output);

  tmin = 360.0 * INIT_YEAR + 30.0 * (INIT_MONTH - 1.0) + INIT_DAY;
  tmax = tmin + DURATION; 

  initializeVariablesROMS();
  InitializeVariablesLocate(np);
  Topology( lp, mp, sp);

  /* TEST: ------------------------------------------------------*/
  for(dir=0; dir<4; dir++)
    {
      sprintf(nameout,"./grid_dir%d.dat", dir);
      output = fopen(nameout,"w");
      printf("qmax[%d]=%d\n", dir, qmax[dir]);
      for(p=0; p<qmax[dir]; p++)
	    {
	      q = qindex[p][dir];
	      fprintf(output,"%lf %lf %lf\n", ptgrid[q].lon, ptgrid[q].lat, ptgrid[q].dpt);
	    }
      fclose(output);
    }

  sprintf(nameout,"./grid_core.dat");
  output = fopen(nameout,"w");
  ncore = (lp - 2) * (mp - 2) * sp;
  for(q=0; q<ncore; q++)
    {
      fprintf(output,"%lf %lf %lf\n", ptgrid[qcore[q]].lon, ptgrid[qcore[q]].lat, ptgrid[qcore[q]].dpt);
    }
  fclose(output);

  for(t=tmin; t<tmax; t=t+0.25) 
    {
      for(ipoint=0; ipoint<np; ipoint++) 
	{
	  if(degree[ipoint]!=0)
	    {
	      if(rk4(t, ipoint, 0.25)==1)
		{
		  //destruimos todos los cross nearest
		  if(cross[ipoint]==1)		    
		    DestroyCross(ipoint); 
		  
		  for(j=0; j<4; j++)
		    {
		      if(crossnn[ipoint][j]==1)
			DestroyCross(inn[ipoint][j]);
		    }
		}
	    }	  
	}

      //Calculation of relative distances:

      for(dir = 0; dir < 2; dir++)
	{
	  for(p=0; p<qmax[dir]; p++)
	    {
	      q = qindex[p][dir];
	      if((degree[q]*degree[inn[q][dir]])!=0)
		{
		  phi1 = RADS * ptm[inn[q][dir]].lon;
		  phi0 = RADS * ptm[q].lon;
		  theta1 = RADS * ptm[inn[q][dir]].lat;
		  theta0 = RADS * ptm[q].lat;
		  cos_Phi = cos(phi1-phi0)*cos(theta1)*cos(theta0) + sin(theta1)*sin(theta0);
		  delta[q][dir] = R_EARTH * sqrt(2.0*(1.0-cos_Phi));
		  delta[inn[q][dir]][dir+2] = delta[q][dir];
		}
	    }
	}

      //Calculation of time escape tau 

      for(q=0; q<np; q++)
	{ 
	  if(cross[q]==1)
	    {
	      //mira las distancias si alguna es mayor que deltamax
	      deltamax = delta[q][0];
	      for(dir = 1; dir < 4; dir++)
		{
		  if(delta[q][dir] > deltamax)
		    deltamax = delta[q][dir];
		}
	      //si es mayor que deltamax: actualiza degrees, tau y delta[]
	      if(deltamax > DELTAMAX)
		{
		  tau[q] = t + 0.25 - tmin;
		  deltaf[q] = deltamax;
		  DestroyCross(q);
		}
	    }
	}
    }


  sprintf(nameout,"./tau.dat");
  output = fopen(nameout,"w");
  for(ipoint=0; ipoint<np; ipoint++) 
    {
      if((ipoint % lp)==0)
	    fprintf(output,"\n");
      fprintf(output,"%lf %lf %lf %lf\n", ptgrid[ipoint].lon, ptgrid[ipoint].lat, ptgrid[ipoint].dpt, 1.0/tau[ipoint]);
    } 
  fclose(output);

  return 0;
}


