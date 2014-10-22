#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Lyapunov3D.h"

extern point ptm[NPMAX];

int rk4(double t, int ipoint, double deltat)
{
	double deltat2,deltat6;
	double deltat_sec, deltat2_sec, deltat6_sec;
	double tau;
	double h;
	point phi, pt;
	vector v1,v2,v3,v4;
	

	deltat = deltat;
	deltat2=deltat*0.5;
	deltat6=deltat/6.0;

	deltat_sec = deltat * SECONDS_DAY;
	deltat2_sec = deltat2 * SECONDS_DAY;
	deltat6_sec = deltat6 * SECONDS_DAY;

	if(LinearInterpolation(t, ipoint, ptm[ipoint], &v1))
	  return 1;

	tau = t + deltat2;
	h = R_EARTH * cos(RADS * ptm[ipoint].lat);
	v1.u = v1.u * (GRADS / h);
	v1.v = v1.v / h;
	SUM_INC(phi, ptm[ipoint], deltat2_sec, v1);
	phi.lat = LAT_TO_MU(phi.mu);
	if(LinearInterpolation(tau, ipoint, phi, &v2))
	  return 1;

	h = R_EARTH * cos(RADS * phi.lat);
	v2.u = v2.u * (GRADS / h);
	v2.v = v2.v / h;
	SUM_INC(phi, ptm[ipoint], deltat2_sec, v2);
	phi.lat = LAT_TO_MU(phi.mu);
	if(LinearInterpolation(tau, ipoint, phi, &v3))
	  return 1;

	h = R_EARTH * cos(RADS * phi.lat);
	v3.u = v3.u * (GRADS / h);
	v3.v = v3.v / h;
	tau = t + deltat;
	SUM_INC(phi, ptm[ipoint], deltat_sec, v3);
	phi.lat = LAT_TO_MU(phi.mu);
	if(LinearInterpolation(tau, ipoint, phi, &v4))
	  return 1;
	h = R_EARTH * cos(RADS * phi.lat);
	v4.u = v4.u * (GRADS / h);
	v4.v = v4.v / h;


	ptm[ipoint].lon = ptm[ipoint].lon + deltat6_sec * (v1.u + v4.u + 2.0 * (v2.u + v3.u));
	ptm[ipoint].mu  = ptm[ipoint].mu  + deltat6_sec * (v1.v + v4.v + 2.0 * (v2.v + v3.v));
	ptm[ipoint].dpt = ptm[ipoint].dpt + deltat6_sec * (v1.w + v4.w + 2.0 * (v2.w + v3.w));
	ptm[ipoint].lat = MU_TO_LAT(ptm[ipoint].mu);

	return 0;
}

