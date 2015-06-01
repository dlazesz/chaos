#include <ostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "vnode.h"

using namespace std;
using namespace vnodelp;

#define VERBOSE_DEBUG false
#define VERBOSE_DEBUG_IT false
#define VERBOSE_DEBUG2 false
#define VERBOSE_DEBUG2_5 true
#define VERBOSE_DEBUG2_5_2 false
#define VERBOSE_DEBUG3 true
#define VERBOSE_DEBUG_SAME false
#define GNUPLOT_IMG true
#define GNUPLOT_OBJ false
#define GNUPLOT_OBJ2 true
#define PRINT2 false
#define GNUPLOT_ESTIMATE true

#define MIN_INTERVAL_RADIUS 1.0E-10
#define MIN_INT_TIME 1.0E-3

struct point3{
        interval t;
        interval x;
        interval y;
        interval z;
        };

struct point3_fs{
        point3 f;
        point3 s;
        };

// Functions:

void discont(void);
point3 aPosterioriEst(point3_fs pm);
vector<point3> branch(point3 p);
bool bound(point3 p);
vector<point3> BB(point3 p,int l);
bool iterate(point3 yp,point3_fs *pm);
point3_fs logsearch(bool &fatal, point3 null, point3 end);
point3 solve(interval tend,point3 p);
void header();
void print(point3 p, double GlobalExcess);
iVector point3ToiVector(point3 p);
point3 iVectorTopoint3(iVector y,interval t);
double LEnd(interval x);
double UEnd(interval x);
string PPrintPoint3(string name, point3 p,bool with_t);
string PPrintGnuPlot(point3 p,bool with_t);
string PPrintPoint3WOZ(string name, point3 p,bool with_t);

// Tests related to the lines:
bool lOFrinRx(point3 p);
bool rOFlinRx(point3 p);
bool rOFrinRx(point3 p);
bool lOFlinRx(point3 p);
bool inRx(point3 p);
bool NOTinRx(point3 p);
bool lOFtinRy(point3 p);
bool uOFbinRy(point3 p);
bool uOFtinRy(point3 p);
bool lOFbinRy(point3 p);
bool inRy(point3 p);
bool NOTinRy(point3 p);
bool NOTinR(point3 p);
bool NOTinRT(point3 p);
bool lOFrinLx(point3 p);
bool rOFlinLx(point3 p);
bool rOFrinLx(point3 p);
bool lOFlinLx(point3 p);
bool inLx(point3 p);
bool NOTinLx(point3 p);
bool lOFtinLy(point3 p);
bool uOFbinLy(point3 p);
bool uOFtinLy(point3 p);
bool lOFbinLy(point3 p);
bool inLy(point3 p);
bool NOTinLy(point3 p);
bool NOTinL(point3 p);
bool NOTinLT(point3 p);
bool uZ(point3 p);
bool aZ(point3 p);
bool oZ(point3 p);
bool aOFZinz(point3 p);
bool uOFZinz(point3 p);
bool onZinz(point3 p);
bool IsSame(point3 p1,point3 p2);

// Need to set these appropriately for the Poincaré section to work correctly
bool firstBetterZ(point3 first, point3 second);
bool secondBetterZ(point3 first, point3 second);
bool firstGoodZ(point3 p);
bool secondGoodZ(point3 p);
bool changeinz(point3 first, point3 second);
bool reversalinz(point3 p);
bool passinz(point3 p);

template<typename var_type>
void Chaos( int              n,     // SIZE is the problem
            var_type*       yp,     // yp is a pointer to output variables
            const var_type*  y,     // y is a pointer to input variables
            var_type         t,     // t is the time variable
            void*           param)  // param is a pointer to additional
{
    interval a(102.0),b(3.0),c(1.0),p(0.0); // OR p(0.83);

    yp[0]=     a*y[1] -b*(y[0]+p);
    yp[1]=  y[0]*y[2] -   y[1];
    yp[2]= -y[0]*y[1] -   y[2]+c;
}

class vnodeEngine{
 private:
 AD* ad;

 public:
 VNODE* Solver;
 interval tend;
 interval z;
 interval zTol;
 int returns;
 double MaxGlobalExcess;

vnodeEngine(){ // Declare and init Solver
    ad= new FADBAD_AD(4,Chaos,Chaos);
    Solver= new VNODE(ad);
    Solver->setOrder(50);
    MaxGlobalExcess=1000.0;
    tend=100.0;
    returns=1;
    z=0.3995;
    zTol=0.0003;
   }

~vnodeEngine(){
    delete ad;
    delete Solver;
}
} VE;

// The first approximation (in time) IT'S DEPENDENT ON THE DE! (Direction of the trajectory)
bool firstBetterZ(point3 old_version, point3 new_version){
return UEnd(old_version.z)<UEnd(new_version.z) && LEnd(old_version.z)<LEnd(new_version.z);
};

// The second approximation (in time) IT'S DEPENDENT ON THE DE! (Direction of the trajectory)
bool secondBetterZ(point3 old_version, point3 new_version){
return UEnd(old_version.z)>UEnd(new_version.z) && LEnd(old_version.z)>LEnd(new_version.z);
};

// The first approximation (in time)
bool firstGoodZ(point3 p){
return uOFZinz(p);
};

// The second approximation (in time)
bool secondGoodZ(point3 p){
return aOFZinz(p);
};

// The position wrt Z is changed
bool changeinz(point3 first, point3 second){
return uOFZinz(first) && !uOFZinz(second);
};

bool reversalinz(point3 p){
return uOFZinz(p);
};

bool passinz(point3 p){
return aOFZinz(p);
};

//------------------------------------------------------------------------

void discont(void){ // Discontinuity searcher

 point3_fs last,ret;
 point3 test0,test,test2;
 bool fatal,first;
 int m=1500,n=500; // Pieces in x and y
 double dx,dy;

 test0.t = 0.0;
 test0.z = 0.3;
 test0.x = interval(-1.739047,-1.68); // Cover of the area to be checked
 test0.y = interval(-0.08, -0.0750931);

 dx=UEnd(test0.x)-LEnd(test0.x);
 dy=UEnd(test0.y)-LEnd(test0.y);

 for(int j=0;j<m;j++)
  {
   test=test0;
   test.x=interval(LEnd(test0.x) + (j*dx)/m,LEnd(test0.x) + (j*dx)/m);
   first=true;

   for(int i=0;i<n;i++)
    {
     test2=test;
     test2.y=interval(LEnd(test0.y) + (i*dy)/n,LEnd(test.y) + (i*dy)/n);
     cout <<"("<<j<<","<<i<<")" <<endl;
     cerr <<"("<<j<<","<<i<<")" <<endl;
     //cout << PPrintPoint3("",test2,false)<<endl;
     fatal=iterate(test2,&ret);
     if (fatal){ cout << "fatal"<< endl; exit(1);}
     if(!first &&
        (fabs(midpoint(last.f.x)-midpoint(ret.f.x))>1.0E-2 ||
         fabs(midpoint(last.s.x)-midpoint(ret.s.x))>1.0E-2 ||
         fabs(midpoint(last.f.y)-midpoint(ret.f.y))>1.0E-2 ||
         fabs(midpoint(last.s.y)-midpoint(ret.s.y))>1.0E-2
        )
      ){
        cout << PPrintPoint3WOZ("Discontinuity:",test2,false)<<endl;
        //cout << PPrintPoint3WOZ("f",ret.f,false)<<endl
        //     << PPrintPoint3WOZ("s",ret.s,false)<<endl;
        cout <<"Delta(lastfx):" << midpoint(last.f.x)-midpoint(ret.f.x) <<endl
             <<"Delta(lastsx):" << midpoint(last.s.x)-midpoint(ret.s.x) <<endl
             <<"Delta(lastfy):" << midpoint(last.f.y)-midpoint(ret.f.y) <<endl
             <<"Delta(lastsy):" << midpoint(last.s.y)-midpoint(ret.s.y) <<endl;
       }
     last=ret;

     first=false;
    } //fori
  } //forj
} // end discont

//------------------------------------------------------------------------

point3 aPosterioriEst(point3_fs pm){
// From the lower and upper approximation we get estimates for
// the projected image on the plane. The estimation itself:
 interval T;
 point3 estimation;

 T=interval(0 , UEnd(pm.s.t) - LEnd(pm.f.t) );
 estimation.x=pm.f.x + T*(102*pm.f.y - 3*pm.f.x);
 estimation.y=pm.f.y + T*( pm.f.x*pm.f.z - pm.f.y);
 estimation.z=pm.f.z;

 if(VERBOSE_DEBUG2_5)
  cout << "aPosterioriEst(Before growth):"
       << PPrintPoint3("",estimation,false) <<endl;
                                     // Theoretical computations based on
 estimation.x+=interval(-18,7 )*T*T; // Grönwall's inequality, relying on T
 estimation.y+=interval(-4,2)*T*T;   // Generously overestimated
                                     // (Maybe better bounds are needed)
 if(VERBOSE_DEBUG2_5)
  cout << "T: ("<< LEnd(T) <<","<< UEnd(T) <<") "<< rad(T) <<endl
       << "aPosterioriEst:"<< PPrintPoint3("",estimation,false) <<endl;

return estimation;
}

//------------------------------------------------------------------------

bool bound(point3 p){
// If intersects the TABOO set, then inmediately branch
// No idea, if it ever stops... Because there is no depth limit

// Not in the TABOO set of L AND Not in the TABOO set of R
// The working mechanism is unknown. At least for now... ;)
return (NOTinLT(p) && NOTinRT(p));
}

//------------------------------------------------------------------------

vector<point3> branch(point3 p){
 point3 temp;
 vector<point3> out;
 interval nlower,nupper;

 if (rad(p.x) > rad (p.y)) // Note: The priority here affects performance
  {   // Cut by X          // If L and R lies parallel to X or Y
   if (rad(p.x)<MIN_INTERVAL_RADIUS) //Not saves us from wrong definition of L and R
    {
       cout << "BRANCH: FATÁL ERROR, INTERVAL IS ONLY ONE POINT LONG!"
            << endl;
       exit(1);
    }

   if(VERBOSE_DEBUG2) cout << "NEXT: Cut by X" <<endl;

   // [a;b]-> [a;a] és [b;b]
   interval lower( LEnd(p.x), LEnd(p.x));
   interval upper( UEnd(p.x), UEnd(p.x));

   // ([a;a] + [b;b]) /2 = [c;d]
   interval middle = (lower+upper)/2;

   // [a;d] && [c;b]
   // => c in [a;d] & d in [c;b]
   interval nlower (LEnd(p.x   ), UEnd(middle));
   interval nupper (LEnd(middle), UEnd(p.x   ));

   temp.x = nlower;
   temp.y = p.y;
   temp.z = p.z;
   temp.t = p.t;
   out.push_back(temp);

   temp.x = nupper;
   temp.y = p.y;
   temp.z = p.z;
   temp.t = p.t;
   out.push_back(temp);
  }
 else
 {
   // Cut by Y
   if (rad(p.y)<MIN_INTERVAL_RADIUS)
   {
       cout <<"BRANCH: FATÁL ERROR, INTERVAL IS ONLY ONE POINT LONG!" <<endl;
       exit(1);
   }

   if(VERBOSE_DEBUG2) cout << "NEXT: Cut by Y" <<endl;

   // [a;b]-> [a;a] és [b;b]
   interval lower( LEnd(p.y), LEnd(p.y));
   interval upper( UEnd(p.y), UEnd(p.y));

   if(VERBOSE_DEBUG2_5_2)
    cout <<"lower("<< LEnd(p.y) <<","<< LEnd(p.y)<<")"<<rad(lower)<<endl
         <<"upper("<< UEnd(p.y) <<","<< UEnd(p.y)<<")"<<rad(upper)<<endl;

   // ([a;a] + [b;b]) /2 = [c;d]
   interval middle = (lower+upper)/2;

   if(VERBOSE_DEBUG2_5_2)
    cout <<"middle"<< midpoint(middle) <<" "<<rad(upper)<<endl;

   // [a;d] && [c;b]
   // => c in [a;d] & d in [c;b]
   nlower=interval (LEnd(p.y   ), UEnd(middle));
   nupper=interval (LEnd(middle), UEnd(p.y   ));

   if(VERBOSE_DEBUG2_5_2)
    cout << "nlower" << midpoint(nlower) <<" "<<rad(nlower  ) <<endl
         << "p.y"    << midpoint(p.y   ) <<" "<<rad(p.y     ) <<endl
         << "nlower" << LEnd(    nlower) <<" "<< UEnd(nlower) <<endl
         << "p.y"    << LEnd(    p.y   ) <<" "<< UEnd(p.y   ) <<endl
         << "L"      <<((LEnd(nlower))==(LEnd(p.y)))<<endl
         << "U"      <<((UEnd(nlower))- (UEnd(p.y)))<<endl;

   if(VERBOSE_DEBUG2)
    cout << "This was given to Branch: " << PPrintPoint3("p",p,false) << endl;

   temp.x = p.x;
   temp.y = nlower;
   temp.z = p.z;
   temp.t = p.t;
   out.push_back(temp);

   if(VERBOSE_DEBUG2)
    cout << "NEXT:lower: " << PPrintPoint3("temp",temp,false) << endl;

   temp.x = p.x;
   temp.y = nupper;
   temp.z = p.z;
   temp.t = p.t;
   out.push_back(temp);

   if(VERBOSE_DEBUG2)
    cout << "NEXT: upper: " << PPrintPoint3("temp",temp,false) << endl;
 }

 return out;
}

//------------------------------------------------------------------------

// Branch & Bound
vector<point3> BB(point3 p,int l){

vector<point3> b,out,ret;
vector<point3>::iterator i,oi;
point3_fs pm;
point3 estimation;
bool fatal=false;

if(GNUPLOT_OBJ) cerr << "GNUPLOT(OBJECT):" << PPrintGnuPlot(p,false) <<endl;
if(VERBOSE_DEBUG2) cout << "NEXT: if NOTinR(p) && NOTinL(p)" <<endl;

 // This is the object interval:
 // If the object interval does not intersect L U R, accept instantly

 if (NOTinR(p) && NOTinL(p))

 {
     out.push_back(p); // Accepted, push in the output vector

     if(VERBOSE_DEBUG2)
      cout << "BOUND(ObjectNOTin):" <<  PPrintPoint3("",p,false) <<endl;
     if(GNUPLOT_OBJ2)
      cerr << "GNUPLOT(OBJECT-NOTin):" << PPrintGnuPlot(p,false) <<endl;
 }
 else
 {
  if(VERBOSE_DEBUG2)
   cout << "NEXT: else NOTinR(p) && NOTinL(p)" <<endl;

  pm.f.x=pm.f.y=pm.f.z=pm.s.x=pm.s.y=pm.s.z=0; // Just for pretty output
  fatal=iterate(p,&pm);

  if(VERBOSE_DEBUG2)
   cout << "After iterate" <<endl;
  if(VERBOSE_DEBUG2_5)
   cout << "imgF: " << PPrintPoint3("",pm.f,false) << endl
        << "imgS: " << PPrintPoint3("",pm.s,false) << endl;

  // From the lower and upper approximation estimates the image on the plane
  // And use this estimation later
  if(!fatal) estimation= aPosterioriEst(pm); // If fatal, no need for estimation
  cout << "bound conditions(true,"<<(true)<<"):"
       << " bound(pm.f):" << bound(pm.f)
       << " bound(pm.s):" << bound(pm.s)
       << " IsSame(pm.f,pm.s):" << IsSame(pm.f,pm.s)
       << " bound(estimation):" << bound(estimation) << endl;
  if(!fatal && bound(pm.f) && bound(pm.s) &&
      IsSame(pm.f,pm.s) && bound(estimation))
  {
    out.push_back(p); // Accepted, push in the output vector

    cout << "BOUND(Accepted Object):" << PPrintPoint3("",p,false) << endl;

    // This is the image interval:
    if(VERBOSE_DEBUG3)
     cout << "NEXT:BOUND(Accepted Image)F " << PPrintPoint3("",pm.f,true)<<endl
          << "NEXT:BOUND(Accepted Image)S " << PPrintPoint3("",pm.s,true)<<endl;

    // GnuPlot data goes to stderr to have clear debug output on stdout
    if(GNUPLOT_IMG)
     cerr << "GNUPLOT:" << PPrintGnuPlot(pm.f,false)  <<endl
          << "GNUPLOT:" << PPrintGnuPlot(pm.s,false)  <<endl
          << "T:"       << UEnd(pm.s.t) - LEnd(pm.f.t)<<endl
          << "Depth:"   << l                          <<endl;

    if(GNUPLOT_OBJ2)
     cerr << "GNUPLOT(OBJECT):"     << PPrintGnuPlot(p,false)           <<endl;

    if(GNUPLOT_ESTIMATE)
     cerr << "GNUPLOT(ESTIMATION):" << PPrintGnuPlot(estimation,false)  <<endl;

  }
  else
  {
    if(VERBOSE_DEBUG2)
    cout << "NEXT: branchF"<< PPrintPoint3("",pm.f,false) << endl
         << "NEXT: branchS"<< PPrintPoint3("",pm.s,false) << endl;

    b=branch(p);
    // Foreach
    for(i=b.begin();i!=b.end();i++)
    {
    if(VERBOSE_DEBUG3)
     cout << "NEXT:BB recursively(Object)" << PPrintPoint3("",*i,false) <<endl;

    ret = BB(*i,l++);
    // Push the accepted ones to the output vector
    for(oi=ret.begin();oi!=ret.end();oi++) out.push_back(*oi);
    }
  }
 }

return out; // Return the produced good branches
}

//------------------------------------------------------------------------

bool iterate(point3 yp,point3_fs *pm){

VNODE* Solver=VE.Solver;
interval tend=VE.tend;
interval t; // Initialize later to yp.t

iVector y(3);
point3 pnull,pend;
point3_fs passtrough_coord;
bool change=false,fatal=false;
int current_return=0;

if(VERBOSE_DEBUG2)
 cout << "NEXT: iterate(this comes in)" << PPrintPoint3("",yp,false) << endl;

// Saves the last and the current interval for the logarithmic search
pnull=yp; // Last
pend=yp;  // Current

// Integrate init
y=point3ToiVector(yp);
t=yp.t; //tnull

// Integrate
if(VERBOSE_DEBUG2) cout << "NEXT: iterate(this comes in)"
                        <<"t0("  << LEnd(t)    <<","<< UEnd(t)    <<")\t"
                        <<"tend("<< LEnd(tend) <<","<< UEnd(tend) <<")"
                        << endl;

Solver->setFirstEntry();
Solver->setHmin(0);
Solver->setOneStep(on);
while(t!=tend)
{
 Solver->integrate(t,y,tend);
 if(Solver->successful()&& Solver->getGlobalExcess()<=VE.MaxGlobalExcess)
  {
   if(VERBOSE_DEBUG_IT)
    {
     header();
     print(iVectorTopoint3(y,t),Solver->getGlobalExcess());
    }
  }
 else
  {
   cout << "FATÁL ERROR(mk)!"
        << midpoint(yp.x) <<"\t"<< midpoint(yp.y) <<endl;
   return true;
  }

 // Updates the saved points
 pnull=pend;
 pend=iVectorTopoint3(y,t);

 // The position wrt Z is changed, we have the first coordinate
 if(changeinz(pnull,pend))
  {
   if(VERBOSE_DEBUG) cout << "Change!" <<endl;
   change=true;
   passtrough_coord.f=pnull;
  }

 // The trajectory turned back. No pass, just intersection. False alarm!
 if(change && reversalinz(pend))
  {
   if(VERBOSE_DEBUG) cout << "False alarm!" <<endl;
   change=false;
  }

 // The trajectory has passed trough Z
 if(change && passinz(pend))
  {
   if(VERBOSE_DEBUG) cout << "Pass!" <<endl;
   change=false;
   current_return++;
   passtrough_coord.s=pend;
  }

 // The number of Poencaré rerurns is enough
 if (current_return>=VE.returns)
  {

   if(VERBOSE_DEBUG)
    cout<< "Reached the "<<VE.returns<<"th pass!" << endl
        << "NEXT: logsearch(this goes in)"
        << PPrintPoint3("pnull",passtrough_coord.f,false) << endl
        << "NEXT: logsearch(this goes in)"
        << PPrintPoint3("pend" ,passtrough_coord.s,false) << endl;

   *pm=logsearch(fatal,passtrough_coord.f,passtrough_coord.s);
    if (fatal) return true;

   // This is different from the usual output for easy reading
   if(VERBOSE_DEBUG){
                     cout << "pf\t"; print(pm->f,0);
                     cout << "ps\t"; print(pm->s,0);
                    }

   break;
  } // If Pass. There is no else!

} // While

if (UEnd(t)>=LEnd(tend)){
 cout << "TERE WEREN'T "<<VE.returns<<" POENCARÉ RETURNS!" <<endl;
 exit(1);
}
return false; // Not fatal
} // Nth iterate

//------------------------------------------------------------------------

// Logarithmic search
point3_fs logsearch(bool &fatal, point3 null, point3 end){

point3_fs current;
point3 pm=null,abs_middle; // Middle interval
point3_fs out; // Final approximations: first and second
interval time;

current.f=null; // Save to preserve the original
current.s=end;

if(VERBOSE_DEBUG)
 cout <<"Before while"
      << endl
      << "NEXT: logsearch(while)" << PPrintPoint3("current.f",current.f,false)
      << endl
      << "NEXT: logsearch(while)" << PPrintPoint3("current.s",current.s,false)
      << endl
      << "NEXT: logsearch(while)" << PPrintPoint3("pm",pm,false)
      << endl
      << "while"<<endl;

// If we find an intersection, then we're done
while (!onZinz(pm))
{
 // null.t -> (current.f.t+current.s.t)/2
 pm=solve((current.f.t+current.s.t)/2,null);

 // Was a good interval found? -> Update the respective one
 if (firstGoodZ (pm)) current.f=pm;
 if (secondGoodZ(pm)) current.s=pm;

 if(VERBOSE_DEBUG)
  cout << "NEXT: logsearch(while)" << PPrintPoint3("pm",pm,false)
       << endl
       << "NEXT: logsearch(while)" << PPrintPoint3("current.f",current.f,false)
       << endl
       << "NEXT: logsearch(while)" << PPrintPoint3("current.s",current.s,false)
       << endl;
} // While

// Intersects Z
if(VERBOSE_DEBUG)
 cout << "Intersects Z:" << PPrintPoint3("current.f",current.f,false)
      << endl
      << "Intersects Z:" << PPrintPoint3("current.s",current.s,false)
      << endl
      << "Intersects Z:" << PPrintPoint3("pm",pm,false)
      << endl;

abs_middle=pm; // This will be the absolute middle

if (!oZ(abs_middle)) {
 cout << "FATÁL ERROR(abs_middle)!" <<endl;
 fatal=true;
 return current;
}

//--------------------------------

// Triing to put closer the first interval to the middle interval

pm = current.f;
time = (current.f.t+abs_middle.t)/2;

// If we find an intersection, then we're done
// We need to keep in mind three things:
// 1) Time must not be too short, unless the zTol is too low
// 2) The "better" point should not intersect the selected point in Z
// 3) If the point is already good enough, don't make it worse!
while (LEnd(time) >= MIN_INT_TIME && !onZinz(pm) && !(uZ(current.f) || aZ(current.f)))
{
 // null.t -> (current.f.t+abs_middle.t)/2
 time = (current.f.t+abs_middle.t)/2;
 pm=solve(time,null);

 // If we find better, update. If not, leave it.
 // Still good AND better then the original
 if (firstGoodZ(pm) && firstBetterZ(current.f,pm)) current.f=pm;

 if(VERBOSE_DEBUG)
  cout << "NEXT: logsearch(while2)" << PPrintPoint3("pm",pm,false)
       << endl
       << "NEXT: logsearch(while2)" << PPrintPoint3("current.f",current.f,false)
       << endl
       << "NEXT: logsearch(while2):" << PPrintPoint3("abs_middle",abs_middle,false)
       << endl;
} // While

// Intersects Z
if(VERBOSE_DEBUG)
 cout << "Intersects Z_2(while)" << PPrintPoint3("current.f",current.f,false)
      << endl
      << "Intersects Z_2(while)" << PPrintPoint3("current.s",current.s,false)
      << endl
      << "Intersects Z_2(while)" << PPrintPoint3("abs_middle",abs_middle,false)
      << endl;

if (LEnd(time)<MIN_INT_TIME) {
 cout << "MIN_INT_TIME REACHED!" <<endl;
 exit(1);
}

//--------------------------------

// Triing to put closer the first interval to the middle interval

pm = current.s;
time = (current.s.t+abs_middle.t)/2;

// If we find an intersection, then we're done
// We need to keep in mind three things:
// 1) Time must not be too short, unless the zTol is too low
// 2) The "better" point should not intersect the selected point in Z
// 3) If the point is already good enough, don't make it worse!
while (LEnd(time)>=MIN_INT_TIME && !onZinz(pm) && !(uZ(current.s) || aZ(current.s)))
{
 // null.t -> (current.s.t+abs_middle.t)/2
 time = (current.s.t+abs_middle.t)/2;
 pm=solve(time,null);

 // If we find better, update. If not, leave it.
 // Still good AND better then the original
 if (secondGoodZ(pm) && secondBetterZ(current.s,pm)) current.s=pm;

 if(VERBOSE_DEBUG)
  cout << "NEXT: logsearch(while3)" << PPrintPoint3("pm",pm,false)
       << endl
       << "NEXT: logsearch(while3)" << PPrintPoint3("current.s",current.s,false)
       << endl
       << "NEXT: logsearch(while3):" << PPrintPoint3("abs_middle",abs_middle,false)
       << endl;
} // While

// Intersects Z
if(VERBOSE_DEBUG)
 cout << "Intersects Z_3(while)" << PPrintPoint3("current.f",current.f,false)
      << endl
      << "Intersects Z_3(while)" << PPrintPoint3("current.s",current.s,false)
      << endl
      << "Intersects Z_3(while)" << PPrintPoint3("abs_middle",abs_middle,false)
      << endl;

if (LEnd(time)<MIN_INT_TIME) {
 cout << "MIN_INT_TIME REACHED!" <<endl;
 exit(1);
}

return current;
} // LogSearch

//------------------------------------------------------------------------

point3 solve(interval tend,point3 p){ // Solve at a given t from p initial value
VNODE* Solver=VE.Solver;
point3 p2=p;

iVector y(3);
y=point3ToiVector(p);

interval t=p.t;

Solver->setFirstEntry();
Solver->setHmin(0);
Solver->setOneStep(on);
while(t!=tend)
{
  Solver->integrate(t,y,tend);
  if(Solver->successful()&& Solver->getGlobalExcess()<=VE.MaxGlobalExcess)
  {
    p2=iVectorTopoint3(y,t);
    if(VERBOSE_DEBUG) print(p2,Solver->getGlobalExcess());
  }
  else
  {
    cout << "FATÁL ERROR(solve)!" <<endl;
    exit(1);
  }
} // While

return p2;
}

//------------------------------------------------------------------------

void header(){

      cout << "midpoint(t)" <<"\t"
           << "midpoint(x)" <<"\t"
           << "midpoint(y)" <<"\t"
           << "midpoint(z)" <<"\t"
           << "rad(t)"      <<"\t"
           << "rad(x)"      <<"\t"
           << "rad(y])"     <<"\t"
           << "rad(z)"      <<"\t"
           << "Solver->getGlobalExcess()"
           << endl;
} // header

void print(point3 p, double GlobalExcess){

if (PRINT2)
    {
      cout <<  "t("<< LEnd(p.t) <<","<< UEnd(p.t) <<")"
           <<"\tx("<< LEnd(p.x) <<","<< UEnd(p.x) <<")"
           <<"\ty("<< LEnd(p.y) <<","<< UEnd(p.y) <<")"
           <<"\tz("<< LEnd(p.z) <<","<< UEnd(p.z) <<")"
           <<"\t"  << GlobalExcess
           << endl;
    }
else{
      cout << midpoint(p.t) <<"\t"
           << midpoint(p.x) <<"\t"
           << midpoint(p.y) <<"\t"
           << midpoint(p.z) <<"\t"
           << rad(p.t) <<"\t"
           << rad(p.x) <<"\t"
           << rad(p.y) <<"\t"
           << rad(p.z) <<"\t"
           << GlobalExcess
           << endl;
    }
}// print

iVector point3ToiVector(point3 p){

iVector y(3);
y[0]=p.x;
y[1]=p.y;
y[2]=p.z;
return y;
} // point3ToiVector

point3 iVectorTopoint3(iVector y,interval t){

point3 p;
p.t=t;
p.x=y[0];
p.y=y[1];
p.z=y[2];
return p;
} // iVectorTopoint3

double LEnd(interval x){ // Gives the lower end of an interval
return midpoint(x)-rad(x);
}

double UEnd(interval x){ // Gives the upper end of an interval
return midpoint(x)+rad(x);
}

string PPrintPoint3(string name, point3 p,bool with_t){
stringstream s;

s.setf(ios::fixed,ios::floatfield);
s.precision(20);

if (with_t) s << name <<"t("<< LEnd(p.t) <<","<< UEnd(p.t) <<")\t";

s << name <<"x("<< LEnd(p.x) <<","<< UEnd(p.x) <<")\t"
  << name <<"y("<< LEnd(p.y) <<","<< UEnd(p.y) <<")\t"
  << name <<"z("<< LEnd(p.z) <<","<< UEnd(p.z) <<")";

return s.str();
}

string PPrintPoint3WOZ(string name, point3 p,bool with_t){
stringstream s;

s.setf(ios::fixed,ios::floatfield);
s.precision(20);

if (with_t) s << name <<"t("<< LEnd(p.t) <<","<< UEnd(p.t) <<")\t";

s << name <<"x("<< LEnd(p.x) <<","<< UEnd(p.x) <<")\t"
  << name <<"y("<< LEnd(p.y) <<","<< UEnd(p.y) <<")";

return s.str();
}

string PPrintGnuPlot(point3 p,bool with_t){
stringstream s;

s.setf(ios::fixed,ios::floatfield);
s.precision(20);

 s <<  midpoint(p.x) <<" "<< midpoint(p.y) <<" "
   << LEnd(p.x)      <<" "<<     UEnd(p.x) <<" "
   << LEnd(p.y)      <<" "<<     UEnd(p.y) <<" #"
   << LEnd(p.z)      <<" "<<     UEnd(p.z) ;

if (with_t)
 s <<" "<< LEnd(p.t) <<","<<     UEnd(p.t) <<")";

return s.str();
}

//------------------------------------------------------------------------

static const struct T_R{
interval m,br,bl,bt,bb;
T_R(){
    interval m_denom=0.007;
    interval m_nom  =0.097;
    m = m_denom/m_nom;
    br=-0.195;
    bl=-0.26;
    bt=m*0.266 - 0.0335;
    bb=m*0.266 - 0.0365;
   }
} R;

static const struct T_L{
interval m,br,bl,bt,bb;
T_L(){
    interval m_denom= 0.0120;
    interval m_nom  = 0.176;
    m = m_denom/m_nom;
    br=-0.6;
    bl=-0.7;
    bt=m*(0.727) - 0.0645;
    bb=m*(0.727) - 0.0675;
   }
} L;

bool lOFrinRx(point3 p){
return (UEnd(p.x - R.br) < 0);
}

bool rOFlinRx(point3 p){
return (LEnd(p.x - R.bl) > 0);
}

bool rOFrinRx(point3 p){ // The whole interval is right of R
return (LEnd(p.x - R.br) > 0);
}

bool lOFlinRx(point3 p){ // The whole interval is left of R
return (UEnd(p.x - R.bl) < 0);
}

bool inRx(point3 p){
return lOFrinRx(p) && rOFlinRx(p);
}

bool NOTinRx(point3 p){ // The whole interval is right OR left of R
return rOFrinRx(p) || lOFlinRx(p);
}

// top bottom upper lower // 0< equation -> Lower then line
// lowerside OF top       // 0> equation -> Upper then line
bool lOFtinRy(point3 p){
return (LEnd(R.m*p.x - p.y + R.bt) > 0);
}

// top bottom upper lower
// upperside OF bottom
bool uOFbinRy(point3 p){
return (UEnd(R.m*p.x - p.y + R.bb) < 0);
}

// top bottom upper lower
// upperside OF top
bool uOFtinRy(point3 p){ // The whole interval is above R
return (UEnd(R.m*p.x - p.y + R.bt) < 0);
}

// top bottom upper lower
// lowerside OF bottom
bool lOFbinRy(point3 p){ // The whole interval is under R
return (LEnd(R.m*p.x - p.y + R.bb) > 0);
}

bool inRy(point3 p){
return lOFtinRy(p) && uOFbinRy(p);
}

bool NOTinRy(point3 p){ // The whole interval is above OR under R
return uOFtinRy(p) || lOFbinRy(p);
}

bool NOTinR(point3 p){
// If the whole interval is NOT in R wrt x [y], then y [x] doesn't matter
// If the inteval intersects R wrt x then it could intersect Taboo set too!
// This is checked elsewhere
return NOTinRx(p) || NOTinRy(p);
}

bool NOTinRT(point3 p){ // NOT in the Taboo set of R
// If the whole interval is in R wrt y, it's not in the Taboo set OR
// If the whole interval is NOT in R wrt x There's no need to check furhter...
return inRy(p) || NOTinRx(p);
}

bool lOFrinLx(point3 p){
return (UEnd(p.y - L.br) < 0);
}

bool rOFlinLx(point3 p){
return (LEnd(p.x - L.bl) > 0);
}

bool rOFrinLx(point3 p){ // The whole interval is right of L
return (LEnd(p.x - L.br) > 0);
}

bool lOFlinLx(point3 p){ // The whole interval is left of L
return (UEnd(p.x - L.bl) < 0);
}

bool inLx(point3 p){
return lOFrinLx(p) && rOFlinLx(p);
}

bool NOTinLx(point3 p){
return rOFrinLx(p) || lOFlinLx(p);
}

// top bottom upper lower
// lowerside OF top
bool lOFtinLy(point3 p){
return (LEnd(L.m*p.x - p.y + L.bt) > 0);
}

// top bottom upper lower
// upperside of bottom
bool uOFbinLy(point3 p){
return (UEnd(L.m*p.x - p.y + L.bb) < 0);
}

// top bottom upper lower
// upperside OF top
bool uOFtinLy(point3 p){ // The whole interval is above L
return (UEnd(L.m*p.x - p.y + L.bt) < 0);
}

// top bottom upper lower
// lowerside OF bottom
bool lOFbinLy(point3 p){ // The whole interval is under L
return (LEnd(L.m*p.x - p.y + L.bb) > 0);
}

bool inLy(point3 p){
return lOFtinLy(p) && uOFbinLy(p);
}

bool NOTinLy(point3 p){
return uOFtinLy(p) || lOFbinLy(p);
}

bool NOTinL(point3 p){
return NOTinLx(p) || NOTinLy(p);
}

bool NOTinLT(point3 p){ // NOT in the Taboo set of L
// If the whole interval is in L wrt y, it's not in the Taboo set OR
// If the whole interval is NOT in L wrt x There's no need to check furhter...
return inLy(p) || NOTinLx(p)  ;
}

bool uZ(point3 p){ // Under Z AND in Tol
return (LEnd(p.z) > UEnd(VE.z - VE.zTol)) &&
       (UEnd(p.z) < LEnd(VE.z)) ;
}

bool aZ(point3 p){ // Above Z AND in Tol
return (LEnd(p.z) > UEnd(VE.z)) &&
       (UEnd(p.z) < LEnd(VE.z + VE.zTol)) ;
}

bool oZ(point3 p){ // ON Z AND in Tol
return (LEnd(p.z) > UEnd(VE.z - VE.zTol)) &&
       (UEnd(p.z) < LEnd(VE.z + VE.zTol)) ;
}

bool aOFZinz(point3 p){ // Above Z
return (LEnd(p.z) > UEnd(VE.z));
}

bool uOFZinz(point3 p){ // Under Z
return (UEnd(p.z) < LEnd(VE.z));
}

bool onZinz(point3 p){ // ON Z (Intersects Z)
return (!aOFZinz(p) && !uOFZinz(p));
}

bool IsSame(point3 p1,point3 p2){
// If the two approximations are not at the same "side" of the quadrangles
// check the stricter rules for both!

  bool noneInRx = NOTinRx(p1)  &&  NOTinRx(p2); // None intersects R in x
  bool  oneInRx = !noneInRx;                    // (At least) ONE intersects R in x
  bool bothInRy =     inRy(p1) &&     inRy(p2); // Both are in R in y
  bool safeInRy = oneInRx && bothInRy;          // One is in Rx => Both must be in Ry to be sure!

  bool noneInLx = NOTinLx(p1)  &&  NOTinLx(p2); // None intersects L in x
  bool  oneInLx = !noneInLx;                    // (At least) ONE intersects L in x
  bool bothInLy =     inLy(p1) &&     inLy(p2); // Both are in L in y
  bool safeInLy = oneInLx && bothInLy;          // One is in Lx => Both must be in Ly to be sure!

  bool sameInRx = (                             // Both points has the same orientation
   lOFrinRx(p1) == lOFrinRx(p2) &&
   rOFlinRx(p1) == rOFlinRx(p2) &&
   rOFrinRx(p1) == rOFrinRx(p2) &&
   lOFlinRx(p1) == lOFlinRx(p2));

  bool sameInLx = (                             // Both points has the same orientation
   lOFrinLx(p1) == lOFrinLx(p2) &&
   rOFlinLx(p1) == rOFlinLx(p2) &&
   rOFrinLx(p1) == rOFrinLx(p2) &&
   lOFlinLx(p1) == lOFlinLx(p2));

  bool sameInRy = (                             // Both points has the same orientation
   lOFtinRy(p1) == lOFtinRy(p2) &&
   uOFbinRy(p1) == uOFbinRy(p2) &&
   uOFtinRy(p1) == uOFtinRy(p2) &&
   lOFbinRy(p1) == lOFbinRy(p2));

  bool sameInLy = (                             // Both points has the same orientation
   lOFtinLy(p1) == lOFtinLy(p2) &&
   uOFbinLy(p1) == uOFbinLy(p2) &&
   uOFtinLy(p1) == uOFtinLy(p2) &&
   lOFbinLy(p1) == lOFbinLy(p2));

  bool uaZ      =       uZ(p1) && aZ(p2);       // In Z one is upper the other is lower
  bool auZ      =       aZ(p1) && uZ(p2);       // In Z one is lower the other is upper

  if (VERBOSE_DEBUG_SAME) {
    cout << "IsSame(true," << true << "):" << endl
         << PPrintPoint3("p1",p1,false)  << endl
         << PPrintPoint3("p2",p2,false)  << endl
         << "noneInRx: "   << noneInRx   << endl
         << "bothInRy: "   << bothInRy   << endl
         << "noneInLx: "   << noneInLx   << endl
         << "bothInLy: "   << bothInLy   << endl
         << "sameInRx: "   << sameInRx   << endl
         << "sameInRy: "   << sameInRx   << endl
         << "sameInLx: "   << sameInRx   << endl
         << "sameInLy: "   << sameInRx   << endl
         << "uaZ || auZ: " << (uaZ || auZ) << endl;
  }
  return
  (safeInRy || sameInRx) &&
  (safeInLy || sameInLx) &&
// These are the vertical checks
  (noneInRx || sameInRy) &&
  (noneInLx || sameInLy) &&
// Check for Z too: p1 should be under/above Z AND p2 should be above/under Z
       (uaZ || auZ);
}

//------------------------------------------------------------------------

int main (int argc, char **argv){
// Just for precision
cout.setf(ios::fixed,ios::floatfield);
cout.precision(20);
cerr.setf(ios::fixed,ios::floatfield);
cerr.precision(20);
cerr << "zTol:" << VE.zTol<< endl;

interval tnull= 0.0;
interval eps = interval(-1,1)*1.0E-8;

//A1-2: (-0.169,-0.0265),(-0.169,-0.0295) //-0.195
//B1-2: (-0.266,-0.0335),(-0.266,-0.0365) //-0.26
//C1-2: (-0.551,-0.0525),(-0.551,-0.0555) //-0.6
//D1-2: (-0.727,-0.0645),(-0.727,-0.0675) //-0.7

point3 test;
test.t = tnull;
test.z = 0.3995 + eps;

//L
//~ test.x=interval(-0.7,-0.6);
//~ test.y=interval(-0.0675,-0.0525);


//R
//~ test.x=interval(-0.26,-0.195);
//~ test.y=interval(-0.0365,-0.0265);

// A line
//~ test.x = interval(-0.195,-0.195);
//~ test.y = interval(-0.0365,-0.0265);

// B line
//~ test.x = interval(-0.26,-0.26);
//~ test.y = interval(-0.0365,-0.0265);

// C line
//~ test.x = interval(-0.6,-0.6);
//~ test.y = interval(-0.0675,-0.0525);

// D line
//~ test.x = interval(-0.7,-0.7);
//~ test.y = interval(-0.0675,-0.0525);

vector<point3> ret;

ret = BB(test,1);
cout << "Final num. of pieces:"<< ret.size() << endl;
cerr << "Final num. of pieces:"<< ret.size() << endl;

return EXIT_SUCCESS;
}
