#include "TROOT.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "Ray.h"
#include <iostream> 

#ifndef __SURFACE__
#define __SURFACE__ 

class PlanarSurface {
 private: 
  double z; 
  double TR;  // this is the transverse radious
  bool debug; 

 public:
  PlanarSurface(double z0,double TR0) {
    z = z0;
    TR = TR0; 
  }

  bool Transport(Ray &ray){
    if( z < ray.GetZ() ) return false; 
    ray.Transport(z); // Transport the ray to this location.
    if( TR < TMath::Sqrt( ray.GetX()*ray.GetX()+ ray.GetY()*ray.GetY()) ) return false; 

    if( debug ) ray.Print(" Ray at the entrance of planar surface "); 
    return true;
  }

  bool Refraction(Ray &ray,double index){
    double sinangle = TMath::Sqrt(1.-ray.GetVZ()*ray.GetVZ());
    double sinagleout = ray.GetIdxR()*sinangle/index; 
    double phi = TMath::ATan2(ray.GetVY(),ray.GetVX());
    
    if( TMath::Abs(sinagleout) > 1. ) return false; 

    double vx = sinagleout*TMath::Cos(phi);
    double vy = sinagleout*TMath::Sin(phi);
    double vz = TMath::Sqrt(1.-sinagleout*sinagleout);

    ray.SetDir(vx,vy,vz); 

    ray.SetIdxR(index);
  
    if( debug ) ray.Print(" Ray at the exit of planar surface "); 
  
    return true; 
  }
  
  void SetDebug(bool a) { debug = a;} 
};

class SphericalSurface {
 private: 
  double R;
  double z0; // Coordinate of the center of the surface circumference.  
  double TR; // Transverse radious  
  bool debug; 

 public:
  SphericalSurface(double zc,double R0,double TR0) {
    z0 = zc;
    R = R0; 
    TR = TR0;
  }

  void SetDebug(bool a) { debug = a;}  

  bool Transport(Ray &ray){
    // Find the intersection. 
    double a = 1.; 
    double b = 2.*(ray.GetVX()*ray.GetX()+ray.GetVY()*ray.GetY()+ray.GetVZ()*(ray.GetZ()-z0));
    double c = ray.GetX()*ray.GetX()+ray.GetY()*ray.GetY()+(ray.GetZ()-z0)*(ray.GetZ()-z0)-R*R;

    double discr = b*b-4*a*c;

    if( discr <= 0 ) { 
      return false;
    } 
    
    double lambda1 = (-b+TMath::Sqrt(discr))/(2.*a);  
    double lambda2 = (-b-TMath::Sqrt(discr))/(2.*a);  

    double z1 = ray.GetZ()+lambda1*ray.GetVZ();
    double z2 = ray.GetZ()+lambda2*ray.GetVZ();
    double z; 

    if( z1 < ray.GetZ() && z2 < ray.GetZ() ) {       
      return false; 
    }
    else if ( z1 > ray.GetZ() && z2 > ray.GetZ() ) 
      z = TMath::Min(z1,z2);
    else if( z1 < ray.GetZ() ) 
      z = z2; 
    else 
      z = z1;

    ray.Transport(z); // Transport the ray to this location.
    if( TR < TMath::Sqrt( ray.GetX()*ray.GetX()+ ray.GetY()*ray.GetY()) ) return false; 

    if( debug ) ray.Print(" Ray at the entrace of spherical surface "); 

    return true;
  }

  bool Refraction(Ray &ray,double index){

    // find the normal to the surface at contact point. 
    // Local reference system 
    double x[3],y[3],z[3]; 

    z[0] = ray.GetX(); 
    z[1] = ray.GetY(); 
    z[2] = ray.GetZ()-z0;

    double aa = TMath::Sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]); 
    z[0]/=aa;
    z[1]/=aa;
    z[2]/=aa;

    x[0] = 1.; x[1] = 0.; x[2] = 0.;

    // Y = Z x X 

    y[0] = z[1]*x[2]-z[2]*x[1];
    y[1] = z[2]*x[0]-z[0]*x[2];
    y[2] = z[0]*x[1]-z[1]*x[0];

    aa = TMath::Sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]); 
    y[0]/=aa;
    y[1]/=aa;
    y[2]/=aa;

    // X = Y x Z 

    x[0] = y[1]*z[2]-y[2]*z[1];
    x[1] = y[2]*z[0]-y[0]*z[2];
    x[2] = y[0]*z[1]-y[1]*z[0];

    //    std::cout << "X ( " << x[0] << " , " << x[1] << " , " << x[2] << " ) " << std::endl;



    double cosangle = z[0]*ray.GetVX()+z[1]*ray.GetVY()+z[2]*ray.GetVZ();
 
    double sinangle = TMath::Sqrt(1.-cosangle*cosangle);
    double sinangleout = ray.GetIdxR()*sinangle/index; 

    if( TMath::Abs(sinangleout) > 1. ) { return false; }

    double phi = TMath::ATan2(y[0]*ray.GetVX()+y[1]*ray.GetVY()+y[2]*ray.GetVZ(),
			      x[0]*ray.GetVX()+x[1]*ray.GetVY()+x[2]*ray.GetVZ());
    
    double cosangleout = TMath::Sqrt(1.-sinangleout*sinangleout);

    double sign = 1.; 

    if( ray.GetZ() < z0 ) sign = -1.; 

    double vx = sinangleout*(TMath::Cos(phi)*x[0]+TMath::Sin(phi)*y[0])+sign*cosangleout*z[0];
    double vy = sinangleout*(TMath::Cos(phi)*x[1]+TMath::Sin(phi)*y[1])+sign*cosangleout*z[1];
    double vz = sinangleout*(TMath::Cos(phi)*x[2]+TMath::Sin(phi)*y[2])+sign*cosangleout*z[2];

    ray.SetDir(vx,vy,vz); 

    ray.SetIdxR(index);
    
    if( debug ) ray.Print(" Ray at the exit of spherical surface "); 

    return true; 
  }


  bool Reflexion(Ray &ray){

    // find the normal to the surface at contact point. 
    // Local reference system 
    double x[3],y[3],z[3]; 

    z[0] = ray.GetX(); 
    z[1] = ray.GetY(); 
    z[2] = ray.GetZ()-z0;

    double aa = TMath::Sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]); 
    z[0]/=aa;
    z[1]/=aa;
    z[2]/=aa;

    x[0] = 1.; x[1] = 0.; x[2] = 0.;

    // Y = Z x X 

    y[0] = z[1]*x[2]-z[2]*x[1];
    y[1] = z[2]*x[0]-z[0]*x[2];
    y[2] = z[0]*x[1]-z[1]*x[0];

    aa = TMath::Sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]); 
    y[0]/=aa;
    y[1]/=aa;
    y[2]/=aa;

    // X = Y x Z 

    x[0] = y[1]*z[2]-y[2]*z[1];
    x[1] = y[2]*z[0]-y[0]*z[2];
    x[2] = y[0]*z[1]-y[1]*z[0];


    double cosangle = z[0]*ray.GetVX()+z[1]*ray.GetVY()+z[2]*ray.GetVZ();
     
    if( cosangle > 0 ) cosangle *=-1.;

    double vx = ray.GetVX()+2.*cosangle*z[0]; 
    double vy = ray.GetVY()+2.*cosangle*z[1]; 
    double vz = ray.GetVZ()+2.*cosangle*z[2]; 

    ray.SetDir(vx,vy,vz);
   
    if( debug ) ray.Print(" Ray at the exit of spherical surface "); 

    return true; 
  }


};


class ParabolicalSurface {
 private: 
  double R;
  double zc;
  double z0; // Coordinate of the center of the surface circumference.  
  double D; // Transverse radious  
  double A;
  bool debug; 

 public:
  ParabolicalSurface(double z0,double R0,double D0) {
    zc = z0;
    R = R0; 
    D = D0;
    A = D/R/R; 
  }

  double GetFocus(void) { return zc-1./4./A;}

  void SetDebug(bool a) { debug = a;}  

  bool Transport(Ray &ray){

    // Find the intersection 
    // ZC-Z = d/R^2 (X*X + Y*Y ) 
    


    double x0 = ray.GetX(); 
    double y0 = ray.GetY(); 
    double z0 = ray.GetZ(); 
    double vx = ray.GetVX(); 
    double vy = ray.GetVY(); 
    double vz = ray.GetVZ(); 

    if( z0 > zc ) return false; // Beyond the parabola.  

    double tx = vx/vz;
    double ty = vy/vz;


    double a = A*(tx*tx+ty*ty); 
    double b = 2.*A*(tx*x0+ty*y0);
    double c = A*(x0*x0+y0*y0);

    // Intersection at (zc-z) = c + a * (z-z0)^2 + b (z-z0) 
    // Taking r = (z-z0) 
    //  zc+z0-r = c + a *r^2 + b r 
    //  a * r2 + (b+1) r + c - zc - z0 = 0. 
    //

    double disc = (b+1)*(b+1) -4.*a*(c-zc-z0); 

    if( disc < 0 ) {
      std::cout << " No solution found " << disc <<  std::endl; 
      return false; 
    }
      
    double r1 = ( -(b+1) - sqrt(disc) ) /(2.*a) + z0 ;  
    double r2 = ( -(b+1) + sqrt(disc) ) /(2.*a) + z0 ;  
   
    double zint = 0; 

    if( r1 > z0 ) 
      zint = r1; 
    else if ( r2 > z0 ) 
      zint = r2; 
    else {
      std::cout << " No solution found " << r1 << "  "<< r2 << "  "<< z0 <<  std::endl; 
      return false;
    } 

    ray.Transport(zint); 

    if( debug ) ray.Print(" Ray at the entrace of spherical surface "); 

    return true;
  }

  bool Reflexion(Ray &ray){

    // find the normal to the surface at contact point. 
    // Local reference system 
    double ep[3],y[3],z[3]; 

    ep[0] = ray.GetX(); 
    ep[1] = ray.GetY(); 
    ep[2] = ray.GetZ();

    double phi = atan2(ep[1],ep[0])+3.141592; 

    double zx = ep[2]-zc; 

    // Surface direction

    double dx = -2.*A*zx*cos(phi)/sqrt(1+4*A*A*zx*zx); 
    double dy = -2.*A*zx*sin(phi)/sqrt(1+4*A*A*zx*zx); 
    double dz = .1/sqrt(1+4*A*A*zx*zx); 

    double cosang = dx*ray.GetVX()+dy*ray.GetVY()+dz*ray.GetVZ();
    
    if( cosang > 0 ) cosang *=-1.;

    double vx = ray.GetVX()+2.*cosang*dx; 
    double vy = ray.GetVY()+2.*cosang*dy; 
    double vz = ray.GetVZ()+2.*cosang*dz; 

    ray.SetDir(vx,vy,vz);

    if( debug ) ray.Print(" Ray at the exit of spherical surface "); 

    return true; 
  }

};


#endif 
