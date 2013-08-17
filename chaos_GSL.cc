#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <sstream>
#include <math.h>

#define VERBOSE_DEBUG false

using namespace std;

struct point3{
        double t;
        double x;
        double y;
        double z;
        };

double* point3ToArray(point3 p);
point3 ArrayTopoint3(double y[],double t);
point3 iterate(point3 yp,int returns);
bool aOFZinz(point3 p);
bool uOFZinz(point3 p);
int func (double t, const double y[], double f[],void *params);
bool lowerLine(point3 p);
bool lowerLineP(point3 p);

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

//------------------------------------------------------------------------

class ODEEngine{

 public:
  size_t dim;
  double h;
  const gsl_odeiv2_step_type *T;
        gsl_odeiv2_step      *s;
        gsl_odeiv2_control   *c;
        gsl_odeiv2_evolve    *e;
        gsl_odeiv2_system   sys;

 double tend;
 double z;
 int status;
 void integrate(point3 &y,double tend){
  double t=y.t;
  double *yp = point3ToArray(y);
  double ypp[3]={yp[0],yp[1],yp[2]}; // I DON'T NOW WHY I NEED THIS!
  status = gsl_odeiv2_evolve_apply (e,c,s,&sys,&t,tend,&h,ypp);
  y=ArrayTopoint3(ypp,t);
 }
 bool success(){ return status == GSL_SUCCESS; };
 bool reset(){ status=gsl_odeiv2_evolve_reset(e); return success();}

ODEEngine(){
   dim=3;
   h = 1e-6;
// Declare and init Solver
   T= gsl_odeiv2_step_rkf45;
   s= gsl_odeiv2_step_alloc (T, dim);
   c= gsl_odeiv2_control_y_new (h, 0.0);
   e= gsl_odeiv2_evolve_alloc (dim);
   sys= {func, NULL, dim, NULL};
   tend=30.0;
   z=0.3;// Lorenz: 0.27;
   }

~ODEEngine(){
   gsl_odeiv2_evolve_free (e);
   gsl_odeiv2_control_free(c);
   gsl_odeiv2_step_free   (s);
}
} OE;

// ENSO
bool changeinz(point3 first, point3 second){
return uOFZinz(first) && !uOFZinz(second);
};

bool reversalinz(point3 p){
return uOFZinz(p);
};

bool passinz(point3 p){
return aOFZinz(p);
};

// Lower "line"
bool lowerLine(point3 p){

double m=(0.207+0.212)/(3.202+3.343);
// less then(<) is the upper crescent-shaped rope
return ((m*(p.x-3.202) -p.y+0.207) >0);
}

// Lower "line" p=0.38
bool lowerLineP(point3 p){

double m=(0.209+0.199)/(2.885+3.531);

return ((m*(p.x-2.885) -p.y+0.209) >0);
}

int func (double t, const double y[], double f[],void *params){
//ENSO
 double p=0.83; // If this is 0.0 then same as lorenz
 f[0]= 102*y[1]-3*(y[0]+p);
 f[1]=  y[0]*y[2]-y[1];
 f[2]= -y[0]*y[1]-y[2]+1;
return GSL_SUCCESS;
};

/*
//Lorenz
/*bool changeinz(point3 first, point3 second){
return aOFZinz(first) && !aOFZinz(second);
};

bool reversalinz(point3 p){
return aOFZinz(p);
};

bool passinz(point3 p){
return uOFZinz(p);
};

// Lower lorenz
bool lowerLine(point3 p){

double m=4.0/3.0;

return ((m*p.x -p.y) >0);
}

int func (double t, const double y[], double f[],void *params){
//Lorenz classic
 f[0]= 10.0*(y[1]-y[0]);
 f[1]=  28.0*y[0]-y[1]-y[0]*y[2];
 f[2]= y[0]*y[1]-(8.0/3.0)*y[2];
return GSL_SUCCESS;
};
*/

point3 iterate(point3 yp,int returns){

double m=0.0;
int current_return=0;
point3 pnull,pend,y,out;
point3 passtrough_coord_f,passtrough_coord_s;
bool change=false;

// Saves the last and the current interval for the logarithmic search
pnull=yp; // Last
pend=yp;  // Current
y=yp; // This will change

// Reset engine
if (!OE.reset())
 {
  cout << "RESET ERROR(mk)!"
       << endl;
  exit(1);
 }

// Integrate
 while (y.t < OE.tend)
 {
   OE.integrate(y,OE.tend);
   if (!OE.success())
    {
     cout << "FATÁL ERROR(mk)!"
          << yp.x <<"\t"<< yp.y <<endl;
     exit(1);
    }

 if(VERBOSE_DEBUG==true) cout << PPrintPoint3("",y,true)<<endl;

 // Updates the saved points
 pnull=pend;
 pend=y;

 // The position wrt Z is changed, we've the first coordinate
 if(changeinz(pnull,pend))
  {
   if(VERBOSE_DEBUG) cout << "Change!" <<endl;
   change=true;
   passtrough_coord_f=pnull;
  }

 // The trajectory turned back. No pass, just intersection. False alarm!
 if(change && reversalinz(pend))
  {
   if(VERBOSE_DEBUG) cout << "False alarm!" <<endl;
   change=false;
  }

 // The trajectory passed trough Z, we've the second coordinate AND passtrough
 if(change && passinz(pend))
  {
   if(VERBOSE_DEBUG==true) cout << "Passtrough" <<endl;
   change=true;
   current_return++;
   passtrough_coord_s=pend;
  }

 // The trajectory has passed trough Z
 if(change && passinz(pend))
  {
   if(VERBOSE_DEBUG) cout << "Pass!" <<endl;
   change=false;
   current_return++;
   passtrough_coord_s=pend;
  }

 // The number of Poencaré rerurns is enough
 if (current_return>=returns)
  {
   if(VERBOSE_DEBUG==true)
    cout<< "Reached the "<<returns<<"th pass!" << endl;

   break; // We're done here
  } // If Pass. There is no else!

} // while

 if (y.t>=OE.tend){
  cout << "TERE WEREN'T "<<returns<<" POENCARÉ RETURNS!" <<endl;
  exit(1);
 }

// Linear interpolation
m=(OE.z-pnull.z)/(pend.z-pnull.z);
out.x=pnull.x +m*(pend.x-pnull.x);
out.y=pnull.y +m*(pend.y-pnull.y);
out.z=OE.z;
out.t=0.0;

return out; // Not fatal
} // iterate

// Above Z
bool aOFZinz(point3 p){
return (p.z > OE.z); // 0.3
}

// Under Z
bool uOFZinz(point3 p){
return (p.z < OE.z); // 0.3
}

void cross_section(point3 p,int ret,long int count){
 long int i=0,j=0;
 point3 first,second;
 second=first=p;

 // Init
 while(i<50)
 {
  second = iterate(first,ret);
  if (lowerLine(first) && lowerLine(second)) i++;
  first=second;
 }

 while(j<count)
 {
  second = iterate(first,ret);
  if (lowerLine(first) && lowerLine(second))
   {
    j++;
    cout <<"x->X:"<< first.x <<" "<< second.x <<" % "<< first.y<<" "<<second.y <<endl;
    cout << j <<endl;
    cerr << j <<endl;
   }
  first=second;
  cout << PPrintPoint3("",second,false)<<endl;
 }
}

void discont(void){

 point3 ret,last,test0,test,test2;
 bool first;
 point3 ul,lr; // upper right lower left
 int m=300,n=300; //Lorenz 500 1500
 double dx,dy;

 test0.t = 0.0;
 test0.z = OE.z;

/*
 ul.x=-8.0001;
 ul.y=8.0001;

 lr.x=8.0;
 lr.y=-8.0;

 dx=fabs(ul.x-lr.x); // This will be positive
 dy=fabs(lr.y-ul.y);
*/

 ul.x=-1.0;
 ul.y=0.05;

 lr.x=-0.9;
 lr.y=0.0;

 dx=fabs(ul.x-lr.x); // This will be positive
 dy=fabs(lr.y-ul.y);

 for(int j=0;j<m;j++)
  {
   test=test0;
   test.x=ul.x + (j*dx)/m;
   first=true;

   for(int i=0;i<n;i++)
    {
     test2=test;
     test2.y=lr.y + (i*dy)/n;
     cout <<"("<<j<<","<<i<<")" <<endl;
     cerr <<"("<<j<<","<<i<<")" <<endl;
     //cout << PPrintPoint3("object",test2,false)<<endl;
     ret=iterate(test2,3);

     if(!first &&
        (fabs(last.x-ret.x)>0.05 || // Lorenz 1.0 1.0
         fabs(last.y-ret.y)>0.05
        )
       ) //5.0E-2
     { cout << PPrintPoint3WOZ("Discontinuity:",test2,false)<<endl;

     //cout << PPrintPoint3WOZ("",ret,false)<<endl;
     cout <<"Delta(x):" << last.x-ret.x <<endl;
     cout <<"Delta(y):" << last.y-ret.y <<endl;
     }
     last=ret;

     first=false;
    } //fori
  } //forj

}

int main (void)
{
cout.setf(ios::fixed,ios::floatfield);
cout.precision(20);
cerr.setf(ios::fixed,ios::floatfield);
cerr.precision(20);

 double t = 0.0;
 double yp[3] = { 0.5, 1, 0.3};
 point3 p = ArrayTopoint3(yp,t);

//cross_section(p,2,100000);
/*
point3 test,test0;
double dy;
int m=1000;
test.z=0.3;
test.t=0.0;
test.x=-1.861;
test.y=-0.076;

 dy=0.002;

 for(int j=0;j<(m-1);j++)
  {
   test0=test;
   test0.y=test.y + (j*dy)/m;
   cout << PPrintPoint3WOZ("object",test0,false)<< endl;
   cout << PPrintPoint3WOZ("f",iterate(test0,2),false)<<endl;

  }
 exit(1);
*/
 discont();

return 0;
}
