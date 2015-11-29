
#include "date.h"
#include "vectorXYZ.h"

extern int verbose;
extern int vflow;
extern vectorXYZ domainmin;
extern vectorXYZ intergrid;
extern vectorXYZ domainmax;
extern date seeddate;
extern double  intstep;
extern double tau;
extern double deltamax;

int GetcmdlineParameters(int narg,char ** cmdarg, string *fnameparams);
int GetfileParameters(string nfileparameters);
