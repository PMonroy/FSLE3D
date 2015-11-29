#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;
#include "readparameters.h"
#include "vectorXYZ.h" 
#include "date.h"
#include "gridconstruction.h" 
#include "integration.h"
#include "vflow.h" 
#include "constants.h"


char velocitydir[] = "/scratch/pmonroy/";
//char velocitydir[] = "/data/geo/escola/roms_benguela/";

date reference_date = {8,  //year
		       1,  //month
		       1,  //day
};

int (*velocity)(double ,vectorXYZ , vectorXYZ* );

int main(int argc, char **argv)
{

  /**********************************************
   * READ COMAND LINE PARAMETERS
   **********************************************/
  string fnameparams;
  if(GetcmdlineParameters(argc, argv, &fnameparams))
    return 1;

  /* VERBOSE */
  if(verbose)
    {
      cout << "READ COMMAND LINE ARGUMENTS:" <<endl;
      cout << " Parameters file: " << fnameparams <<endl;
      cout << endl;
    }

  /**********************************************
   * READ PARAMETERS FROM FILE
   **********************************************/
  if(GetfileParameters(fnameparams))
    return 1;

  /* VERBOSE */
  if(verbose == 1)
    {
      cout << "PARAMETERS FROM FILE: "<< fnameparams <<endl; 
      cout << " vflow = "<<vflow<<endl ;
      cout << " domainmin = "<< domainmin<<endl;
      cout << " intergrid = " <<intergrid<<endl;
      cout << " domainmax = "<< domainmax<<endl;
      cout << " seeddate = "<< seeddate.day<<"-"<<seeddate.month<<"-"<<seeddate.year<<endl ;
      cout << " intstep = "<<intstep<<endl ;
      cout << " tau = "<<tau<<endl ;
      cout << " deltamax = "<<deltamax<<endl ;
      cout << endl;
    }

 

  /**********************************************
   * GRID CONSTRUCTION
   **********************************************/
  vector<vectorXYZ> grid;
  int ni,nj,nk;
  if(gridfsle3d(&grid, &ni, &nj, &nk, domainmin, intergrid, domainmax))
    {
      cout << "Error in contruction grid" << endl;
      return 1;
    }
  if(verbose == 1)
    {
      cout << "GRID CONSTRUCTION: "<<endl; 
      cout << " num. nodes = "<< grid.size()<<endl;
      cout << " ni = "<<ni<<endl;
      cout << " nj = "<<nj<<endl;
      cout << " nk = "<<nk<<endl;
      cout << endl;
    }
  vector<int> neighbor;
  neighbor = neighbors(ni, nj, nk);
 
  /**********************************************
   * INITIAL relative distances
   **********************************************/

  vectorXYZ delta,scalefactor;
  vector<double> ilength;
  unsigned int p;

  ilength.reserve(6*grid.size());

  for(unsigned int q=0; q<grid.size(); q++)
    {
      p=6*q;
      for(int dir=0; dir<6; dir++)
	{
	  if(neighbor[p+dir]>=0)
	    {
	      delta=grid[neighbor[p+dir]]-grid[q];
	      
	      delta.x=rads*delta.x;
	      delta.y=rads*delta.y;
	      
	      scalefactor.x=rearth*cos(rads*grid[q].y); 
	      scalefactor.y=rearth; 
	      scalefactor.z=1.0;

	      delta=delta*scalefactor;
	      delta*=delta;

	      ilength.push_back(sqrt(delta.x+delta.y+delta.z));
	    }
	  else
	    ilength.push_back(-1.0);
	}
    }

  /**********************************************
   * INITIALIZE VARIBLES
   **********************************************/

  vector<double> exit_time;
  vector<double> response;

  exit_time.reserve(grid.size());
  response.reserve(grid.size());

  for(unsigned int q=0; q<grid.size(); q++)
    {
      p=6*q;
      if(neighbor[p]>0 && neighbor[p+1]>0 && neighbor[p+3]>0 && neighbor[p+4]>0)
	{
	  exit_time.push_back(tau+intstep);
	  response.push_back(1.0);
	}
     else
	{
	  exit_time.push_back(0.0);
	  response.push_back(0.0);
	}
    }

  /**********************************************
   * TIME LOOP
   **********************************************/

  // Setup velocity field
  if(loadvflowgrid(reference_date, velocitydir)!=0)
    {
      cout << "Error in reading reference date netcdf"<< endl;
      return 1;
    }

  if(loadvflow(seeddate, (int) tau, velocitydir)!=0)
    {
      cout << "Error in reading velocities"<< endl;
      return 1;
    }

   double tstart;
   double tend;
   double h;
   int ascnd;
   ascnd = tau > 0;
   if(ascnd)
    {
      tstart = 0.0;
      tend = tau-intstep;
      h=intstep;
    }
  else
    {
      tend = 0.0+intstep;
      tstart = abs(tau);
      h=-1.0*intstep;
    }

   double t;
   int count;

   vector<vectorXYZ> tracer;
   tracer = grid;

   vector<double> length;
   length = ilength;

   double lengthmax;
   int dirmax;

   int n0, n1, n3, n4; 
   int ndir;
   for(t=tstart,count=0; ((t<tend)==ascnd) || (t==tend); t+=h,count++)
     {

       /* Compute the position of tracer in time t=t+h*/
       for (unsigned int q = 0; q<tracer.size() ; q++)
	 {
	   //index four neighboirs
	   n0 = neighbor[6*q];
	   n1 = neighbor[6*q+1];
	   n3 = neighbor[6*q+3];
	   n4 = neighbor[6*q+4];
	   
	   if((exit_time[q]>tau) ||
	      (n0>=0 && exit_time[n0]>tau) ||
	      (n1>=0 && exit_time[n1]>tau) ||
	      (n3>=0 && exit_time[n3]>tau) ||
	      (n4>=0 && exit_time[n4]>tau))
	     {
	       if(RK4(t, h, &tracer[q], getvflow))// Semi-implicit 4th order Runge-Kutta
		 {
		   exit_time[q]=0.0;
		   if(n0>=0) exit_time[n0]=0.0;
		   if(n1>=0) exit_time[n1]=0.0;
		   if(n3>=0) exit_time[n3]=0.0;
		   if(n4>=0) exit_time[n4]=0.0;
		 }
	     }
	 }
       /* Compute the relative distances */
       for(unsigned int q=0; q<tracer.size(); q++)
	 {
	   p=6*q;
	   for(int dir=0; dir<2; dir++)
	     {
	       ndir = neighbor[p+dir];
	       if((exit_time[q]<tau && ndir>=0) || (ndir>=0 && exit_time[ndir]))
		 {
		   delta=tracer[ndir]-tracer[q];
		   
		   delta.x=rads*delta.x;
		   delta.y=rads*delta.y;
		   
		   scalefactor.x=rearth*cos(rads*grid[q].y); 
		   scalefactor.y=rearth; 
		   scalefactor.z=1.0;
		   
		   delta=delta*scalefactor;
		   delta*=delta;
		   
		   length[p+dir]=sqrt(delta.x+delta.y+delta.z);
		   length[6*ndir+3]=length[p+dir];
		 }
	     }
	 }

       /* Compute the max length*/
       for(unsigned int q=0; q<tracer.size(); q++)
	 {
	   if(exit_time[q]>tau)
	     {
	       p=6*q;
	       lengthmax=length[p];
	       for(int dir=0; dir<2; dir++)
		 {
		   if(length[p+dir]>lengthmax)
		     {
		       lengthmax=length[p+1];
		       dirmax=dir;
		     }
		   if(length[p+dir+3]>lengthmax)
		     {
		       lengthmax=length[p+dir+3];
		       dirmax=dir+3;
		     }
		 }

	       if(lengthmax>deltamax)
		 {
		   exit_time[q]=t+h;
		   response[q]=lengthmax/ilength[p+dirmax];
		 }
	     }
	 }
     }

  /****************
   * FREE MEMORY
   ****************/

  freevflowgrid();
  freevflow((int) (tau));

  /**********************************************************
   * COMPUTE FINITE SIZE LYAPUNOV EXPONENT
   **********************************************************/

  vector<double> fsle;
  fsle.reserve(grid.size());

  for(unsigned int q=0; q<tracer.size(); q++)
	 {
	   if(exit_time[q]>0.0 && exit_time[q]<=tau)
	     {
	       fsle.push_back(log(response[q])/exit_time[q]);
	     }
	   else
	     fsle.push_back(0.0); 
	 }

 /**********************************************************
   * WRITING RESULT IN VTK FILE
   **********************************************************/
  ofstream ofile("fsle_debug.vtk");

  ofile<<"# vtk DataFile Version 3.0"<<endl;
  ofile<<"Complete vector field of ROMS Benguela"<<endl; 
  ofile<<"ASCII"<<endl;
  ofile<<"DATASET STRUCTURED_GRID"<<endl;
  ofile<<"DIMENSIONS "<<ni<<" "<<nj<<" "<<nk<<endl;
  ofile<<"POINTS "<<ni*nj*nk<<" float"<<endl;
  for(unsigned int q=0; q<grid.size(); q++) 
    {
      ofile<<grid[q].x<<" "<<grid[q].y<<" "<<grid[q].z/1000.0 <<endl;
    }
  ofile<<endl;
  ofile<<"POINT_DATA "<<ni*nj*nk<<endl;
  ofile<<"SCALARS fsle float"<<endl;
  ofile<<"LOOKUP_TABLE default"<<endl;
  for(unsigned int q=0; q<fsle.size(); q++) 
    {
      ofile<<fsle[q]<<endl;
    }
 
  return 0;
}

