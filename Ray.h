#include "TROOT.h"
#include "TMath.h"
#include "TMatrixD.h"
#include <iostream> 

// Always assume the optical system goes along Z.
#ifndef __RAY__
#define __RAY__

class Ray{

 private: 
  double x,y,z;  // origin. 
  
  double xin,yin,zin; 

  double vx,vy,vz; // direction
  double vxin,vyin,vzin;

  double indexrefraction; 

  bool debug; 

 public: 

  Ray(double x0,double y0,double z0, double vx0,double vy0, double vz0, double indx) {

    xin = x = x0; 
    yin = y = y0; 
    zin = z = z0; 

    indexrefraction = indx;

    debug = false; 
    double norm = TMath::Sqrt(vx0*vx0+vy0*vy0+vz0*vz0);

    vxin = vx = vx0/norm;
    vyin = vy = vy0/norm;
    vzin = vz = vz0/norm;
   }

  void Transport(double znew){
    double lambda = (znew-z)/vz;
    z = znew; 
    x+=lambda*vx;
    y+=lambda*vy;

    if( debug ) Print("Transport");
  } 

  double GetX(void) { return x;}
  double GetY(void) { return y;}
  double GetZ(void) { return z;}

  double GetVX(void) { return vx;}
  double GetVY(void) { return vy;}
  double GetVZ(void) { return vz;}

  double GetIdxR(void) { return indexrefraction; } 

  void SetIdxR(double a) { indexrefraction = a; } 

  void SetPos(double ax,double ay,double az) { x = ax; y = ay;z = az;}

  void SetDir(double ax,double ay,double az) { 
    double norm = TMath::Sqrt(ax*ax+ay*ay+az*az);
    vx = ax/norm; vy = ay/norm;vz = az/norm;
  }

  void SetDebug(bool a) { debug = a;} 
 
  void Print(string label = "" ){

    std::cout << " --------------------------------------------------- " << std::endl;
    std::cout << " ---------------- Ray --------------- " << std::endl;
    std::cout << label << std::endl;
    std::cout << " Position   (" << x << " , " << y << " , " << z << " ) " << std::endl;
    std::cout << " Direction  (" << vx << " , " << vy << " , " << vz << " ) " << std::endl;
    std::cout << " Index of refraction " << indexrefraction << std::endl;
    std::cout << " Initial Position   (" << xin << " , " << yin << " , " << zin << " ) " << std::endl;
    std::cout << " Initial Direction  (" << vxin << " , " << vyin << " , " << vzin << " ) " << std::endl;
    std::cout << " --------------------------------------------------- " << std::endl;
    
  }

};

#endif
