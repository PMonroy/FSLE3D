
#define L 194 // Max index xi of rho points 
#define LU 193 //Max index xi of u momentum = L-1
#define M 314 // Max index eta of rho points
#define MV 313 // Max index eta of v momemtun = M-1
#define T  30 // How many days in every file
#define S  32 // Max index s of rho points 

#define LON_MIN 3.833333 // Westernmost longitude
#define LON_MAX 19.91667 // Easternmost longitude 

#define LAT_MIN -35.64948 // Southernmost latitude
#define LAT_MAX -12.05083 // Northernmost_latitude

#define DELTA_T    0.25 // Time increment
#define INIT_YEAR  8    // Only two years [8,9]
#define INIT_MONTH 9    // 1 <= (num. month) <= 12
#define INIT_DAY   17    // 0 <= (num. days) < 30
#define DURATION 120


#define DELTA0 0.0277 //resolution in degrees
#define DELTA0_DEPTH 20.0
#define DELTAMAX 100000.0

#define LON_P_MIN 8.0 // Westernmost longitude
#define LON_P_MAX 16.0 // Easternmost longitude 

#define LAT_P_MIN -30.0 // Southernmost latitude
#define LAT_P_MAX -20.0 // Northernmost_latitude

#define DPT_P_MIN -120.0 // Southernmost latitude
#define DPT_P_MAX -120.0 // Northernmost_latitude

#define NPMAX 5000000

/* Path to the directoy of files :*/ 
#define PATH "/scratch/pmonroy/"

/* Math and Physic Constants*/
#define PI 3.1415926  
#define RADS PI/180.0
#define GRADS 180.0/PI
#define R_EARTH 6371000.0
#define SECONDS_DAY 86400.0


/* Handle errors in the lecture of .nc files by printing 
 * an error message and exiting with a  non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#define TIME_TO_DATE(u,x)	     \
  u.year = (x) / 360;		     \
  u.month = ((x) % 360) / 30 + 1;    \
  u.day = ((x) % 360) % 30;

#define LAT_TO_MU(x)  log(fabs((1.0/cos((PI/180.0)*(x)))+tan((PI/180.0)*(x)))); 
#define MU_TO_LAT(x) (180.0/PI)*(PI/2.0 - 2.0 * atan(exp(-1.0*(x)))); // lat result in degrees

#define INTERPOLATION(Y, y, j, alpha, beta)	\
  Y.u = alpha * y[j].u + beta * y[j+1].u;	\
  Y.v = alpha * y[j].v + beta * y[j+1].v;	\
  Y.w = alpha * y[j].w + beta * y[j+1].w	

#define SUM_INC(R, r, alpha, vec)		  \
  R.lon = r.lon + (alpha) * vec.u;		  \
  R.mu  = r.mu  + (alpha) * vec.v;		  \
  R.dpt = r.dpt + (alpha) * vec.w

        
typedef struct{int year, month, day;} date;
typedef struct{double lon, lat, mu, r, dpt;} point;
typedef struct{double u, v, w;} vector;



int ReadDepth(int t0, int ipoint, double dpt[2][2][2][S]);
int ReadV(int t0, int ipoint, vector vf[2][2][2][2]);

int initializeVariablesROMS(void );
int initializeVariablesTopology(int n);

void Topology( int lp, int mp, int sp);
void CreateCross(int q);
void DestroyCross(int q);

int LocateBox(double t, int ipoint, point pt, point proms[16], vector vroms[16]);
void hunt(double *xx, unsigned long n, double x, unsigned long *jlo);

int LinearInterpolation(double t, int ipoint, point pt, vector *vint);

int rk4(double t, int ipoint, double deltat);
