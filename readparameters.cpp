#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <string>
#include <sstream>
#include <unistd.h>
#include <fstream>
using namespace std;
#include "readparameters.h"

// GLOBALS PARAMTERS
int verbose=0; //Verbosity is disabled by default

int vflow;
vectorXYZ domainmin;
vectorXYZ intergrid;
vectorXYZ domainmax;
date seeddate;
double  intstep;
double tau;
double deltamax;

/* FSLE  domain:*/

//FUNTIONS
int GetcmdlineParameters(int narg,char ** cmdarg, string *fnameparams)
{ 
  string me="GetcmdlineParameters";
  int opt;

  int verboseflag=0; // Verbose flag option
  int fnameparamsflag=0; // File parameters flag option

  while ((opt = getopt(narg, cmdarg, "f:vh")) != -1) 
    {
      switch(opt)
	{

	case 'f':
	  fnameparamsflag++;
	  if (optarg)
	    *fnameparams = optarg;
	  else
	    return 1;
	  break;

	case 'v':
	  verboseflag++;
	  verbose=1;
	  break;

	case 'h':
	  cout<<"Usage: "<< cmdarg[0] <<" [OPTIONS]"<<endl;
	  cout<<" -f [file]      Where [file] is the input file parmeters (mandatory)" <<endl;
	  cout<<" -h             Print this help and exit (optional)"<<endl;
         return 1;

	case '?':
	  cout << "Try "<<cmdarg[0]<<" -h for more information"<<endl;
	  return 1;
	}
    }

  /* Check mandatory repeated or options*/
  if(fnameparamsflag==0)// File parameters option
    {
      cout << "Error("<<me <<"): Option -f <file> is mandatory" <<endl;
      cout << "Try "<<cmdarg[0]<<" -h for more information"<<endl;
      return 1;
    }
  else if (fnameparamsflag>1)
    cout << "Warning("<<me <<"): Option -f is repeated. Parameters file now is the last one ("<<*fnameparams<<")"<<endl;
  
  if(verboseflag>1)// Verbose option
    {
      cout << "Warning("<<me <<"): Option -v is repeated. Verbose mode continues to be enable. "<<endl;
    }

  return 0;
}

int GetfileParameters(string nfileparameters)
{
  string me="readparams()";
  ifstream fparameters(nfileparameters.c_str());

  enum  enum_parameters { VFLOW,
			  DOMAINMIN,
			  INTERGRID,
			  DOMAINMAX,
			  SEEDDATE,
			  INTSTEP,
			  TAU,
			  DELTAMAX,
			  NPARAMETERS
  };
  string pname[NPARAMETERS];

  pname[VFLOW]="vflow";
  pname[DOMAINMIN]="domainmin";
  pname[INTERGRID]="intergrid";
  pname[DOMAINMAX]="domainmax";
  pname[SEEDDATE]="seeddate";
  pname[INTSTEP]="intstep";
  pname[TAU]="tau";
  pname[DELTAMAX]="deltamax";

  int delimiter,end;
  string name, value;
  stringstream svalue;
  string line;

  int *pflag;
  int parameter;
  
  if (!fparameters.is_open())
    {
      cout << me <<": Skipping unreadable file \"" << nfileparameters.c_str() << "\" "<<endl; 
      return 1;
    }

  pflag = (int*) calloc (NPARAMETERS,sizeof(int));

  while(!fparameters.eof())
    {
      getline(fparameters,line);

      if (line[0] == '#') 
	continue;  /* ignore comment line which starts with #*/

      if (isspace(line[0])) 
	  continue; /* ignore blank line */

      delimiter = line.find("=");
      end = line.length();

      value = line.substr(delimiter+1, end);
      name = line.substr(0, delimiter);

      for(parameter=0; parameter<NPARAMETERS; parameter++)
	{
	  if (name.compare(pname[parameter]) == 0)
	    {
	      pflag[parameter]++;
	      if(pflag[parameter]>1)
		{
		  cout  << me << ": Parameter "<< name << " repeated"<<endl;
		  return 1;
		}
	       if ((delimiter+1) == end) 
		 {
		   cout << me << ": Parameter "<< name << " has no value"<<endl;
		   return 1;
		 }
	       switch (parameter) 
		{
		case VFLOW:
                  vflow = atoi(value.c_str());
		  break;
		case DOMAINMIN: 
		  svalue<<value;
		  if((svalue>>domainmin)==0)
		    {
		      cout  << me << ": Format of domainmin is incorrect" << endl;
		      return 1;
		    }
		  svalue.clear();//clear any bits set
		  svalue.str(string());
		  break;

		case INTERGRID:
		  svalue << value;		  
		  if((svalue>>intergrid)==0)
		    {
		      cout  << me << ": Format of domainmin is incorrect" << endl;
		      return 1;
		    }
		  svalue.clear();//clear any bits set
		  svalue.str(string());
		  break;
		case DOMAINMAX:
		  svalue << value;		  
		  if((svalue>>domainmax)==0)
		    {
		      cout  << me << ": Format of domainmin is incorrect" << endl;
		      return 1;
		    }
		  svalue.clear();//clear any bits set
		  svalue.str(string());
		  break;
		case SEEDDATE:
		  if(sscanf(value.c_str(),"%u-%u-%u",&seeddate.day,&seeddate.month,&seeddate.year)!=3)
		    {
		      cout  << me << ": Date format of seeddate is incorrect" << endl;
		      return 1;
		    }
		  break;
		case INTSTEP:
		  intstep = atof(value.c_str());
		  break;
		case TAU:
		  tau = atof(value.c_str());
		  break;  
		case DELTAMAX:
		  deltamax = atof(value.c_str());
		  break;  
		default:
		  cout  << me << ": Unknown parameter "<< name <<endl;
		  return 1;
		}
	       break;
	    }
	}
      if(parameter==NPARAMETERS)
	{
	  cout  << me << ": Unknown parameter "<< name <<endl ;
	  return 1;
	}
    }

  for(parameter=0; parameter<NPARAMETERS; parameter++)
	{
	  if(pflag[parameter]==0)
	    {
	      cout  << me << ": parameter "<< pname[parameter] << " is not defined"<<endl;
	      return 1;
	    }
	}  
  fparameters.close();

  free(pflag);

  return 0;
}
