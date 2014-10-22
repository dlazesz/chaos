#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <ctime>

#define VERBOSE_DEBUG false
// Restriction to the subset of the Z plane
// where the vector field points upward:
// Cut if x,y,z satisfies: x*y >1-z
#define CUT_THE_EDGES false

using namespace std;

//-----------------------------------------------------------------------------
// Struct for better handling of points & prettyprint, conversion rutines

struct point3{
        double t;
        double x;
        double y;
        double z;
        };

double* point3ToArray(point3 p){

double* ret;
double y[3];
ret=y;
y[0]=p.x;
y[1]=p.y;
y[2]=p.z;

return ret;
} // point3ToArray

point3 ArrayTopoint3(double y[],double t){

point3 p;
p.t=t;
p.x=y[0];
p.y=y[1];
p.z=y[2];
return p;
} // ArrayTopoint3

string PPrintPoint3(string name, point3 p,bool with_t){
stringstream s;

s.setf(ios::fixed,ios::floatfield);
s.precision(20);

if (with_t)
   s << name <<"t("<< p.t <<")\t";

   s << name <<"x("<< p.x <<")\t"
     << name <<"y("<< p.y <<")\t"
     << name <<"z("<< p.z <<")";

return s.str();
}

string PPrintPoint3WOZ(string name, point3 p,bool with_t){
stringstream s;

s.setf(ios::fixed,ios::floatfield);
s.precision(20);

if (with_t)
   s << name <<"t("<< p.t <<")\t";

   s << name <<"x("<< p.x <<")\t"
     << name <<"y("<< p.y <<")";

return s.str();
}

//-----------------------------------------------------------------------------

// General function signature
int (*fun)(double t, const double y[], double f[],void *params);

// Real global functions for the equations
int func (double t, const double y[], double f[],void *params);
int funcP(double t, const double y[], double f[],void *params);

bool (*cond)(point3, point3);

class ODEEngine{

 public:
  const gsl_odeiv2_step_type *T;
        gsl_odeiv2_step      *s;
        gsl_odeiv2_control   *c;
        gsl_odeiv2_evolve    *e;
        gsl_odeiv2_system   sys;

 size_t dim;
 double h;
 double tend;
 double z;
 int status;
 void integrate(point3 &y){
  double t=y.t;
  double *yp = point3ToArray(y);
  double ypp[3]={yp[0],yp[1],yp[2]}; // I DON'T NOW WHY I NEED THIS!
  status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,tend,&h,ypp);
  y=ArrayTopoint3(ypp,t);
 }
 bool success(){ return status == GSL_SUCCESS; };
 bool reset(){ status=gsl_odeiv2_evolve_reset(e); return success();}

ODEEngine(){
   h    = 1e-6;
   tend = 10000.0;
   z    = 0.3;// Lorenz: 0.27;
   dim  = 3;
// Declare and init Solver
   T= gsl_odeiv2_step_rkf45;
   s= gsl_odeiv2_step_alloc(T, dim);
   c= gsl_odeiv2_control_y_new(h, 0.0);
   e= gsl_odeiv2_evolve_alloc(dim);
   sys= {func, NULL, dim, NULL};
   }

ODEEngine(double h, 
		  int (DEfun)(double t, const double y[], double f[],void *params), 
		  double tend, 
		  double z) : h(h), tend(tend), z(z), dim(3) {
// Declare and init Solver
   T= gsl_odeiv2_step_rkf45;
   s= gsl_odeiv2_step_alloc(T, dim);
   c= gsl_odeiv2_control_y_new(h, 0.0);
   e= gsl_odeiv2_evolve_alloc(dim);
   sys= {DEfun, NULL, dim, NULL};
   }

~ODEEngine(){
   gsl_odeiv2_evolve_free (e);
   gsl_odeiv2_control_free(c);
   gsl_odeiv2_step_free   (s);
}
};

//-----------------------------------------------------------------------------
// Point orientations per equations

// Above Z independent from the actual equations
bool aOFZinz(ODEEngine *OE, point3 p){
return (p.z > OE->z);
}

// Under Z independent from the actual equations
bool uOFZinz(ODEEngine *OE, point3 p){
return (p.z < OE->z);
}

// These depend from the actual equations
bool changeinz(ODEEngine *OE, point3 first, point3 second){
return uOFZinz(OE, first) && !uOFZinz(OE, second);
};

bool reversalinz(ODEEngine *OE, point3 p){
return uOFZinz(OE, p);
};

bool passinz(ODEEngine *OE, point3 p){
return aOFZinz(OE, p);
};

int func(double t, const double y[], double f[],void *params){
//ENSO
 double p=0.0; // If this is 0.0 then same as lorenz
 f[0]= 102*y[1]-3*(y[0]+p);
 f[1]=  y[0]*y[2]-y[1];
 f[2]= -y[0]*y[1]-y[2]+1;
return GSL_SUCCESS;
};

// Lower "line" Originally for 0.3 will not work for z=0.3995
bool lowerLine(point3 p){

double m=(0.207+0.212)/(3.202+3.343);
// less then(<) is the upper crescent-shaped rope
return ((m*(p.x-3.202) -p.y+0.207) >0);
}

// Lower "line" For z=0.3995 will not work for z=0.3
bool lowerLine2(point3 p){

double m=(0.203706+0.203207)/(2.81225+2.81443);
// less then(<) is the upper crescent-shaped rope
return ((m*(p.x+2.81443)-p.y-0.203207) >0);
}

// Line to separate the returning fold from the main "line"
bool return_fold(point3 p){

double m=(-0.0145424-0.084999)/(-4.39417+3.41740);
// less then(<) is the upper crescent-shaped rope
return ((m*(p.x+3.41740)-p.y+0.084999) <0);
}

bool upperLine(point3 p){ return !lowerLine(p); }

bool all(point3 first, point3 second){ return true; }

bool Upper_Upper(point3 first, point3 second){
    return !lowerLine(first) && !lowerLine(second);
}

bool Lower_Lower(point3 first, point3 second){
    return lowerLine(first) && lowerLine(second);
}

bool Upper_Lower(point3 first, point3 second){
    return !lowerLine(first) && lowerLine(second);
}

bool Lower_Upper(point3 first, point3 second){
    return lowerLine(first) && !lowerLine(second);
}

//-----------------------------------------------------------------------------

int funcP(double t, const double y[], double f[],void *params){
//ENSO
 double p=0.83; // If this is 0.0 then same as lorenz
 f[0]= 102*y[1]-3*(y[0]+p);
 f[1]=  y[0]*y[2]-y[1];
 f[2]= -y[0]*y[1]-y[2]+1;
return GSL_SUCCESS;
};

// Lower "line" p=0.83
bool lowerLineP(point3 p){

double m=(0.209+0.199)/(2.885+3.531);
// less then(<) is the upper crescent-shaped rope
return ((m*(p.x-2.885) -p.y+0.209) >0);
}

bool upperLineP(point3 p){ return !lowerLineP(p); }

bool allP(point3 first, point3 second){ return true; }

bool Upper_UpperP(point3 first, point3 second){
    return !lowerLineP(first) && !lowerLineP(second);
}

bool Lower_LowerP(point3 first, point3 second){
    return lowerLineP(first) && lowerLineP(second);
}

bool Upper_LowerP(point3 first, point3 second){
    return !lowerLineP(first) && lowerLineP(second);
}

bool Lower_UpperP(point3 first, point3 second){
    return lowerLineP(first) && !lowerLineP(second);
}

//-----------------------------------------------------------------------------

//Lorenz
bool changeinzL(ODEEngine *OE, point3 first, point3 second){
return aOFZinz(OE, first) && !aOFZinz(OE, second);
};

bool reversalinzL(ODEEngine *OE, point3 p){
return aOFZinz(OE, p);
};

bool passinzL(ODEEngine *OE, point3 p){
return uOFZinz(OE, p);
};

int funcL(double t, const double y[], double f[],void *params){
//Lorenz classic
 f[0]= 10.0*(y[1]-y[0]);
 f[1]=  28.0*y[0]-y[1]-y[0]*y[2];
 f[2]= y[0]*y[1]-(8.0/3.0)*y[2];
return GSL_SUCCESS;
};

// Lower lorenz
bool lowerLineL(point3 p){

double m=4.0/3.0;

return ((m*p.x -p.y) >0);
}

bool upperLineL(point3 p){ return !lowerLineL(p); }

bool allL(point3 first, point3 second){ return true; }

bool Upper_UpperL(point3 first, point3 second){
    return !lowerLineL(first) && !lowerLineL(second);
}

bool Lower_LowerL(point3 first, point3 second){
    return lowerLineL(first) && lowerLineL(second);
}

bool Upper_LowerL(point3 first, point3 second){
    return !lowerLineL(first) && lowerLineL(second);
}

bool Lower_UpperL(point3 first, point3 second){
    return lowerLineL(first) && !lowerLineL(second);
}

//-----------------------------------------------------------------------------

point3 iterate(ODEEngine *OE, point3 yp, int returns){

yp.t =0.0;
double m=0.0;
int current_return=0;
point3 pnull,pend,y,out;
point3 passtrough_coord_f,passtrough_coord_s;
bool change=false;

// Saves the last and the current interval for the logarithmic search
// (See the other programs using VNODE)
pnull=yp; // Last
pend=yp;  // Current
y=yp; // This will change

// Reset engine
if (!OE->reset())
 {
  cout << "RESET ERROR(mk)!"
       << endl;
  exit(1);
 }

// Integrate
 while (y.t < OE->tend)
 {
   OE->integrate(y);
   if (!OE->success())
    {
     cout << "FATÁL ERROR(mk)!"
          << yp.x <<"\t"<< yp.y <<endl;
     exit(1);
    }

 if(VERBOSE_DEBUG) cout << PPrintPoint3("",y,true)<<endl;

 // Updates the saved points
 pnull=pend;
 pend=y;

 // The position wrt Z is changed, we've the first coordinate
 if(changeinz(OE, pnull,pend))
  {
   if(VERBOSE_DEBUG) cout << "Change!" <<endl;
   change=true;
   passtrough_coord_f=pnull;
  }

 // The trajectory turned back. No pass, just intersection. False alarm!
 if(change && reversalinz(OE, pend))
  {
   if(VERBOSE_DEBUG) cout << "False alarm!" <<endl;
   change=false;
  }

 // The trajectory passed trough Z, we've the second coordinate AND passtrough
 if(change && passinz(OE, pend))
  {
   if(VERBOSE_DEBUG) cout << "Pass!" <<endl;
   change=false;
   current_return++;
   passtrough_coord_s=pend;
  }

 // The number of Poencaré rerurns is enough
 if (current_return>=returns)
  {
   if(VERBOSE_DEBUG)
    cout<< "Reached the "<<returns<<"th pass!" << endl;

   break; // We're done here
  } // If Pass. There is no else!

} // while

 if (y.t>=OE->tend){
  cout << "THERE WEREN'T "<<returns<<" POENCARÉ RETURNS!" <<endl;
  cerr << "THERE WEREN'T "<<returns<<" POENCARÉ RETURNS!" <<endl;
  double t = y.t;                  // In batch processing we can't exit
  double yp[3] = { 0.0, 0.0, 0.0}; // Just tell the outer program what happened
  point3 p = ArrayTopoint3(yp,t);
  return p;
 }

// Linear interpolation
m=(OE->z-pnull.z)/(pend.z-pnull.z);
out.x=pnull.x +m*(pend.x-pnull.x);
out.y=pnull.y +m*(pend.y-pnull.y);
out.z=OE->z;
out.t=pend.t;

return out; // Not fatal
} // iterate

void cross_section(ODEEngine *OE, point3 p, int ret, long int count,
                    bool (separationCond)(point3 p, point3 p2), ofstream &fs){
 long int i=0,j=0;
 point3 first,second;
 second=first=p;

 // Init
 while(i<50)
 {
  second = iterate(OE,first,ret);
  if (second.x==0.0 and second.y==0.0 and second.z==0.0) {
            cout << "Equation time elapsed(50):" << second.t << endl; return; }
  if (separationCond(first, second)
      #if CUT_THE_EDGES
      && (second.x*second.y)<1-second.z && (first.x*first.y)<1-second.z
      # endif
     ) i++;
  first=second;
 }

 while(j<count)
 {
  second = iterate(OE, first,ret);
  if (second.x==0.0 and second.y==0.0 and second.z==0.0) {
            cout << "Equation time elapsed:" << second.t << endl; return; }
  if (separationCond(first, second)
      #if CUT_THE_EDGES
      && (second.x*second.y)<1-second.z && (first.x*first.y)<1-second.z
      # endif
     )
   {
    j++;
    fs <<  first.x << " " <<  first.y << " " <<  first.z << " " <<  first.t
       << " "
       << second.x << " " << second.y << " " << second.z << " " << second.t
       << endl;
    if (j % 1000 ==0) cerr << j <<endl;
   }
  first=second;
 }
}

// Parameters:
// ul,lr; // upper right lower left
// m,n; The size of the grid. For Lorenz 500 1500
// disc_x,disc_y // Min. jump distance for discontinuity. For Lorenz 1.0 1.0
void discont(ODEEngine *OE, int returns, point3 ul, point3 lr, int m, int n,
                                                double disc_x, double disc_y){

 bool first;
 point3 ret,last,test,test2,test2_last;

 // This will be the positive diameter of the intervals
 double dx=fabs(ul.x-lr.x);
 double dy=fabs(lr.y-ul.y);

 for(int j=0;j<m;j++)
  {
   test.t = 0.0;
   test.z = OE->z;
   test.x=ul.x + (j*dx)/m;
   first=true;

   for(int i=0;i<n;i++)
    {
     test2=test;
     test2.y=lr.y + (i*dy)/n;
     cout <<"("<<j<<","<<i<<")" <<endl;
     if (j % 10 == 0)
     cerr <<"("<<j<<","<<i<<")" <<endl;
     ret=iterate(OE, test2,returns);
     if(!first &&
        (fabs(last.x-ret.x)>disc_x ||
         fabs(last.y-ret.y)>disc_y
        )
       )
     {
    cout <<"Discontinuity: "
         << test2.x      << " " << test2.y       << " "
         << test2_last.x << " " << test2_last.y  << endl
         <<"Delta(x):" << last.x-ret.x <<endl
         <<"Delta(y):" << last.y-ret.y <<endl;
     }
     last=ret;
     test2_last=test2;

     first=false;
    } //fori
  } //forj

}

// Solve till T_end, no linear interpolation to plane Z_0
void solve(ODEEngine *OE, point3 yp){

point3 y=yp; // This will change

// Reset engine
if (!OE->reset())
 {
  cout << "RESET ERROR(mk)!"
       << endl;
  exit(1);
 }

// Integrate
 while (y.t < OE->tend)
 {
   OE->integrate(y);
   if (!OE->success())
    {
     cout << "FATÁL ERROR(mk)!"
          << yp.x <<"\t"<< yp.y <<endl;
     exit(1);
    }

if (y.t > 100) cout << y.x << " " << y.y << " "<< y.z << " " << y.t << endl;
if (fmod(y.t,1000) < 10.0) cerr << y.t<< endl;

} // while


return; // Not fatal
} // solve

void make_plot_data_to_file(
                  int num_of_returns,
                  bool (lineCond)(point3 p, point3 p2),
                  double z,
                  double h,
                  int (f)(double t, const double y[], double f[],void *params),
                  double tend,
                  point3 init_point,
                  long int num_of_points,
                  ofstream &fs)
{

    ODEEngine* OE = new ODEEngine(h, f, tend, z);
    cross_section(OE, init_point, num_of_returns, num_of_points, lineCond, fs);
    delete OE;
}

void make_plot_data_batch(
                  int num_of_returns,
                  bool (lineCond)(point3 p, point3 p2),
                  double z_from,
                  double z_to,
                  double z_step,
                  double h,
                  int (f)(double t, const double y[], double f[],void *params),
                  double tend,
                  point3 init_point,
                  long int num_of_points,
                  string ext)
{

 time_t start,end;
 string name;
 string type;
 string prefix;
 ofstream fs;

      if (lineCond == Lower_Lower  or
          lineCond == Lower_LowerP or
          lineCond == Lower_LowerL) type="LL";
 else if (lineCond == Upper_Upper  or
          lineCond == Upper_UpperP or
          lineCond == Upper_UpperL) type="UU";
 else if (lineCond == Upper_Lower  or
          lineCond == Upper_LowerP or
          lineCond == Upper_LowerL) type="UL";
 else if (lineCond == Lower_Upper  or
          lineCond == Lower_UpperP or
          lineCond == Lower_UpperL) type="LU";
 else if (lineCond == all) type="CC";
 else { cerr << "ERROR: INVALID lineCond SPECIFIED!" << endl; exit(1); }

      if (f == func ) prefix="P00";
 else if (f == funcP) prefix="P83";
 else if (f == funcL) prefix="LCS";
 else { cerr << "ERROR: INVALID f SPECIFIED!" << endl; exit(1); }

 for (double z= z_from; z <= z_to; z = z + z_step){
    ODEEngine* OE = new ODEEngine(h, f, tend, z);
    name = prefix + "_" + to_string(num_of_returns) + "_" + type + "_" +
           to_string(z) + "." + ext;
    cout << "opening file: "<<name <<endl;
    fs.open (name, std::fstream::out);
    if ( !fs.is_open() ) {
       cout<< "The file could not be opened!"<<endl;
       }
    fs.setf(ios::fixed,ios::floatfield);
    fs.precision(20);
    time (&start);
    cross_section(OE, init_point, num_of_returns, num_of_points, lineCond, fs);
    time (&end);
    double dif = difftime(end,start);
    cout  << "Elasped time in seconds: " << dif << endl;
    fs.close();
    delete OE;
 }
}

void make_plot_data_from_cin(
                  int num_of_returns,
                  bool (lineCond)(point3 p, point3 p2),
                  double z_from,
                  double z_to,
                  double z_step,
                  double h,
                  int (f)(double t, const double y[], double f[],void *params),
                  double tend,
                  point3 init_point,
                  long int num_of_points,
                  string name)
{

ofstream fs;
cout << "opening file: "<<name <<endl;
fs.open(name, std::ofstream::out | std::ofstream::app);
if ( !fs.is_open() ) {
 cout<< "The file could not be opened!"<<endl;
 }
fs.setf(ios::fixed,ios::floatfield);
fs.precision(20);

while (cin){
cin >> init_point.x >> init_point.y;
cout << init_point.x << " " << init_point.y << endl;
fs << init_point.x << " " << init_point.y << " ";
make_plot_data_to_file(num_of_returns, cond, z_from, h, fun, 
                                           tend, init_point, num_of_points,fs);
}
fs.close();
}

int main (int argc, char **argv)
{
cout.setf(ios::fixed,ios::floatfield);
cout.precision(20);
cerr.setf(ios::fixed,ios::floatfield);
cerr.precision(20);
cin.setf(ios::fixed,ios::floatfield);
cin.precision(20);

double t = 0.0;
double yp[3] = { 0.5, 1, 0.3};
point3 init_point = ArrayTopoint3(yp,t);

//-----------------------------------------------------------------------------

     int num_of_returns   = 2;
long int number_of_points = 100000;
double h = 1e-12;
double tend = 100000.0;
double z_from = -0.6;
double z_to   = -0.4;
double z_step = 0.0005;
string ext("log");

int section_cond_num = 1;
cond=&all;
fun = &func;

//-----------------------------------------------------------------------------

if (argc!=1 && argc!=4) {cerr << "Argument count must be 0 or 3! Currently: "
                              << argc << endl; exit(1);}
else if (argc==4){
 num_of_returns   = atoi(argv[1]);
 number_of_points = atol(argv[2]);
 section_cond_num = atoi(argv[3]);

 switch (section_cond_num){
  case 1:  cond=&Upper_Upper;  break;
  case 2:  cond=&Lower_Lower;  break;
  case 3:  cond=&Upper_Lower;  break;
  case 4:  cond=&Lower_Upper;  break;
  case 5:  cond=&all;          break;
  case 6:  cond=&Upper_UpperP; break;
  case 7:  cond=&Lower_LowerP; break;
  case 8:  cond=&Upper_LowerP; break;
  case 9:  cond=&Lower_UpperP; break;
  case 10: cond=&allP;         break;
  case 11: cond=&Upper_UpperL; break;
  case 12: cond=&Lower_LowerL; break;
  case 13: cond=&Upper_LowerL; break;
  case 14: cond=&Lower_UpperL; break;
  case 15: cond=&allL;         break;

  default:
   cerr << "condition error: " << section_cond_num << endl
        << "Options:"          << endl
        << "Symmetric ENSO:"   << endl
        << "1  Upper_Upper"    << endl
        << "2  Lower_Lower"    << endl
        << "3  Upper_Lower"    << endl
        << "4  Lower_Upper"    << endl
        << "5  All"            << endl
        << "Asymmetric ENSO:"  << endl
        << "6  Upper_Upper"    << endl
        << "7  Lower_Lower"    << endl
        << "8  Upper_Lower"    << endl
        << "9  Lower_Upper"    << endl
        << "10 All"            << endl
        << "Lorenz Classic:"   << endl
        << "11 Upper_Upper"    << endl
        << "12 Lower_Lower"    << endl
        << "13 Upper_Lower"    << endl
        << "14 Lower_Upper"    << endl
        << "15 All"            << endl;
   exit(1);
 }
      if (section_cond_num>=1  && section_cond_num<=5 ) fun = &func;
 else if (section_cond_num>=6  && section_cond_num<=10) fun = &funcP;
 else if (section_cond_num>=10 && section_cond_num<=15) fun = &funcL;
 cout << "#using values for paramter:" << num_of_returns << " "
      << number_of_points << " " << section_cond_num << endl;
} else cout << "#using default values for paramter:" << num_of_returns <<" "
            << number_of_points << " " << section_cond_num << endl;

ODEEngine* OE = new ODEEngine(h, func, tend, z_from);

//-----------------------------------------------------------------------------

// Batch processing
make_plot_data_batch(num_of_returns, cond, z_from, z_to, z_step, h, fun,
                                      tend, init_point, number_of_points, ext);

// Get Points fromc cin
make_plot_data_from_cin(num_of_returns, cond, z_from, z_to, z_step, h, fun,
                               tend, init_point, number_of_points, "data.txt");

// Solve till T_end -> Butterfly
solve(OE, init_point); // Butterfly tend = 550

// Discontinuity searcher
// Parameters:
// ul,lr; // upper right lower left
// m,n; The size of the grid. For Lorenz 500 1500
// disc_x,disc_y // Min. jump distance for discontinuity. For Lorenz 1.0 1.0
point3 ul,lr;
int m=2000, n=2000;
double disc_x=1.0, disc_y=1.0;
ul.x=-8.0001;
ul.y=8.0001;
lr.x=8.0;
lr.y=-8.0;
discont(OE, num_of_returns, ul, lr, m, n, disc_x, disc_y);

// Divide a vertical/horizontal line given by it's endpoints
// to *max* points and use them as initial values
double dy=0.002;
int max=1000;
point3 test;
test.x=-1.861;
test.y=-0.076;
ofstream fs("data.txt");
for(int j=0;j<(max-1);j++){
   init_point=test;
   init_point.y=init_point.y + (j*dy)/m;
   make_plot_data_to_file(num_of_returns, cond, z_from, h, fun,
                                        tend, init_point, number_of_points,fs);
}
fs.close();

delete OE;
return EXIT_SUCCESS;
}
