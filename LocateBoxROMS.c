#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <netcdf.h>
#include "Lyapunov3D.h"

extern point ptm[NPMAX];

/* Longitude and latitude of rho-points */
static double lon[L]; // Longitude(xi) 0 =< xi < L
static double lat[M]; // Latitude(eta) 0 =< eta < M
static double mu[M];  // mu(eta) 0 =< eta < M 

static unsigned long xib[NPMAX], etab[NPMAX], sb[NPMAX][2][2][2];


int initializeVariablesROMS(int np)//Esto habrá que cambiarlo
{
  char file_name_nc[255]; // Name of the file
  int ncid;                    // The netCDF ID for the file
  int lon_rho_id, lat_rho_id;  // and data variable

  static size_t start[] = { 0, 0};
  static size_t count[] = { M, L};
  static ptrdiff_t stride[] = { 1, 1};  

  int xi, eta, retval; //Loop indexes, and error handling.

  int ipoint, i, j, k;

  for(ipoint=0; ipoint<np; ipoint++)
    {
      xib[ipoint] = 0;
      etab[ipoint] = 0;
      for(i=0; i<2; i++)
	{
	  for(j=0; j<2; j++)
	    {
	      for(k=0; k<2; k++)
		{
		  sb[ipoint][i][j][k] = 0;
		}
	    }
	}
    }

  sprintf(file_name_nc,"%sextract_roms_avg_Y%dM%d.nc.1", PATH, INIT_YEAR, INIT_MONTH);
  
  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
   * to the file.*/
  if ((retval = nc_open(file_name_nc, NC_NOWRITE, &ncid)))
    ERR(retval);
  
  /* Get the varid of the data variable, based on its name. */
  if ((retval = nc_inq_varid(ncid, "lon_rho", &lon_rho_id)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncid, "lat_rho", &lat_rho_id)))
    ERR(retval);
      
  /* Read the data. */

  count[0] = 1;
  count[1] = L;

  if ((retval = nc_get_vars_double(ncid, lon_rho_id, start, count, stride, &lon[0])))
    ERR(retval);

  count[0] = M;
  count[1] = 1;

  if ((retval = nc_get_vars_double(ncid, lat_rho_id, start, count, stride, &lat[0])))
    ERR(retval);
    
  for(eta = 0; eta < M; eta++)
      mu[eta] = LAT_TO_MU(lat[eta]);

  if ((retval = nc_close(ncid)))
	    ERR(retval);

  return 0;
}

int LocateBox(double t, int ipoint, point pt, point proms[16], vector vroms[16])
{

  double *llon, *mmu;
  unsigned long eta, xi, s ;
  unsigned long xib_1offset, etab_1offset;

  static double dpt[2][2][2][S];

  int troms;

  double *ddpt;
  unsigned long sb_1offset;
  int i,j,k,tau;

  int index;
  vector vroms_array[2][2][2][2];

  int contw;
  /* Localizamos xib y etab */

  llon = lon - 1;
  xib_1offset = xib[ipoint] + 1;
  hunt(llon, L, pt.lon, &xib_1offset);

  if(xib_1offset == 0 || xib_1offset == L)
    return 1;
  else
    xib[ipoint] = xib_1offset - 1;

  mmu = mu - 1;
  etab_1offset = etab[ipoint] + 1;
  hunt(mmu, M, pt.mu, &etab_1offset);

  if(etab_1offset == 0 || etab_1offset == M)
    return 1;
  else
    etab[ipoint] = etab_1offset - 1;

  /* Leemos las columnas de profundidad */
  
  troms = (int) t;
  ReadDepth(troms, ipoint, dpt);
 
  /* Localizamos la profundidad */

  for(tau = 0; tau < 2; tau++)
    {
      for(i = 0; i < 2; i++)
	{
	  for(j = 0; j < 2; j++)
	    {
	      ddpt = dpt[tau][i][j] - 1;
	      sb_1offset = sb[ipoint][tau][i][j] + 1;
	      hunt(ddpt, S, pt.dpt, &sb_1offset);
	     
	      if(sb_1offset == 0 || sb_1offset == S)
		return 1;
	      else
		{
		  sb[ipoint][tau][i][j] = sb_1offset - 1;
		}
	    }
	}
    }

  /* Leemos las velocidades*/

  ReadV(troms, ipoint, vroms_array);

  /* Vectors and points with only one int index*/

  index = 0;
  for(tau = 0; tau < 2; tau++)
    {
      for(i=0; i < 2; i++)
	{
	  for(j = 0; j < 2; j++)
	    {
	      for(k = 0; k < 2; k++)
		{
		  xi = xib[ipoint]+i;
		  eta = etab[ipoint]+j;
		  s = sb[ipoint][tau][i][j]+k;
		  proms[index].lon = lon[xi];
		  proms[index].mu = mu[eta];
		  proms[index].dpt = dpt[tau][i][j][s];

		  vroms[index]=vroms_array[tau][i][j][k];
		  if(vroms[index].w==0)
		    contw++;
		  index++; 
		}
	    }
	}
    }

  if(contw==16)
    return 1;

  return 0;
}

void hunt(double *xx, unsigned long n, double x, unsigned long *jlo)
/* Given an array xx[1..n], and given a value x, returns a value jlo such that x is between
 * xx[jlo] and xx[jlo+1]. xx[1..n] must be monotonic, either increasing or decreasing.
 * jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input is taken as the
 * initial guess for jlo on output.
 */
{
  unsigned long jm,jhi,inc;
  int ascnd;

  ascnd=(xx[n] >= xx[1]);
  if (*jlo <= 0 || *jlo > n) {
    jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    if (x >= xx[*jlo] == ascnd) {
      if (*jlo == n) return;
      jhi=(*jlo)+1;
      while (x >= xx[jhi] == ascnd) {
	*jlo=jhi;
	inc += inc;
	jhi=(*jlo)+inc;
	if (jhi > n) {
	  jhi=n+1;
	  break;
	}
      }
    } else {
      if (*jlo == 1) {
	*jlo=0;
	return;
      }
      jhi=(*jlo)--;
      while (x < xx[*jlo] == ascnd) {
	jhi=(*jlo);
	inc <<= 1;
	if (inc >= jhi) {
	  *jlo=0;
	  break;
	}
	else *jlo=jhi-inc;
      }
    }
  }
  while (jhi-(*jlo) != 1) {
      jm=(jhi+(*jlo)) >> 1;
      if (x >= xx[jm] == ascnd)
	*jlo=jm;
      else
	jhi=jm;
  }
  if (x == xx[n]) *jlo=n-1;
  if (x == xx[1]) *jlo=1;
}

int ReadDepth(int t0, int ipoint, double dpt[2][2][2][S])
{
  char file_name_nc[255];
  date tdate[2]; // Date format for time
  
  double depth[2][S][2][2];

  static size_t start[] = { 0, 0, 0, 0};
  static size_t count[] = { 1, S, 2, 2};
  static ptrdiff_t stride[] = { 1, 1, 1, 1};
  
  int ncid, retval;
  int depth_id;
  int t1, i, j, k, s;

  t1 = t0 + 1;

  TIME_TO_DATE( tdate[0], t0);
  TIME_TO_DATE( tdate[1], t1);

  start[2] = etab[ipoint];
  start[3] = xib[ipoint];
 
  count[0] = 1 + (tdate[0].month == tdate[1].month);

  for(i = 0; i < 2; i += count[0])
    {
      sprintf(file_name_nc,"%sextract_roms_avg_Y%dM%d.nc.1", PATH, tdate[i].year, tdate[i].month);
  
      /* Open the file. NC_NOWRITE tells netCDF we want read-only access
       * to the file.*/
      if ((retval = nc_open(file_name_nc, NC_NOWRITE, &ncid)))
	ERR(retval);
      
      if ((retval = nc_inq_varid(ncid, "depth", &depth_id)))
	ERR(retval);

      start[0] = tdate[i].day;

      if ((retval = nc_get_vars_double(ncid, depth_id, start, count, stride, &depth[i][0][0][0])))
	ERR(retval);
  
      /* Close the file, freeing all resources. */
      if ((retval = nc_close(ncid)))
	ERR(retval);
    
    }

  /* Tranpose the depth matrix*/
  
  for(k = 0; k < 2; k++)
    {
      for(i = 0; i < 2; i++)
	{
	  for(j = 0; j < 2; j++)
	    {
	      for(s = 0; s < S; s++)
		{
		  dpt[k][i][j][s] = depth[k][s][j][i];
		}
	    }
	}
    }
  return 0;
}

int ReadV(int t0, int ipoint, vector vf[2][2][2][2])
{
  char file_name_nc[255];
  date tdate[2]; // Date format for time
  
  static size_t index[] = { 0, 0, 0, 0};

  static double u[2][2][2][2], v[2][2][2][2], w[2][2][2][2];
  
  int ncid, retval;
  int u_id, v_id, w_id;
  int t1;
  unsigned long  tau, i, j, k;
  unsigned long xi, eta;

  t1 = t0 + 1;

  TIME_TO_DATE( tdate[0], t0);
  TIME_TO_DATE( tdate[1], t1);

  for(tau = 0; tau < 2; tau++)
    {
      sprintf(file_name_nc,"%sextract_roms_avg_Y%dM%d.nc.1", PATH, tdate[tau].year, tdate[tau].month);


      /* Open the file. NC_NOWRITE tells netCDF we want read-only access
       * to the file.*/
      if ((retval = nc_open(file_name_nc, NC_NOWRITE, &ncid)))
	ERR(retval);
      
      if ((retval = nc_inq_varid(ncid, "u", &u_id)))
	ERR(retval);
      if ((retval = nc_inq_varid(ncid, "v", &v_id)))
	ERR(retval);
      if ((retval = nc_inq_varid(ncid, "w", &w_id)))
	ERR(retval);



      index[0] = tdate[tau].day;

      for(k = 0; k < 2; k++)
	{
	  for(i = 0; i < 2; i++)
	    {
	      xi = xib[ipoint] + i;
	      for(j = 0; j < 2; j++)
		{
		  eta = etab[ipoint] + j;

		  index[1] = sb[ipoint][tau][i][j] + k;
		  /* Velocity u:*/	
		  index[2] = eta;
		  index[3] = xi - (xi!=0);

		  if ((retval = nc_get_var1_double(ncid, u_id, index, &u[tau][k][j][i])))
		    ERR(retval);

		  vf[tau][i][j][k].u = u[tau][k][j][i];

		  index[3] = xi - (xi==(L-1));
		  if ((retval = nc_get_var1_double(ncid, u_id, index, &u[tau][k][j][i])))
		  ERR(retval);

		  vf[tau][i][j][k].u += u[tau][k][j][i];
		  vf[tau][i][j][k].u = vf[tau][i][j][k].u / 2.00;

		  /* Velocity v:*/

		  index[2] = eta - (eta!=0);
		  index[3] = xi;

		  if ((retval = nc_get_var1_double(ncid, v_id, index, &v[tau][k][j][i])))
		    ERR(retval);

		  vf[tau][i][j][k].v = v[tau][k][j][i];

		  index[2] = eta - (eta==(M-1));

		  if ((retval = nc_get_var1_double(ncid, v_id, index, &v[tau][k][j][i])))
		  ERR(retval);

		  vf[tau][i][j][k].v += v[tau][k][j][i];
		  vf[tau][i][j][k].v = vf[tau][i][j][k].v / 2.00;
		  
		   /* Velocity w */

		  index[2] = eta;
		  index[3] = xi;
		  if ((retval = nc_get_var1_double(ncid, w_id, index, &w[tau][k][j][i])))
		    ERR(retval);
		  vf[tau][i][j][k].w = w[tau][k][j][i]; 
		}
	    }
	}
      /* Close the file, freeing all resources. */
      if ((retval = nc_close(ncid)))
	ERR(retval);
    }

  return 0;
}
