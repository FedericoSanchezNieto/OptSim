#include "TROOT.h"
#include "TMath.h"
#include "TMatrixD.h"
#include <iostream> 
#include "Ray.h"
#include "Surface.h"

#ifndef __OPTICS__
#define __OPTICS__

class Element {

 public: 
  int m_id; 

 Element(int id): m_id( id ) {;}
  
  int getId() const { return m_id; }
   

}; 


class Disc: public Element {
 private: 
  double TR; 
  double posz; // Center of the window
  double width; // Width of the window 
  PlanarSurface *inputsurface; 
  PlanarSurface *outputsurface; 
  double indexrefraction;
  bool debug; 

 public: 
    Disc(double posz0, double width0, double TR0,  double indexrefraction0, int id = 0 ):Element(id){
    posz = posz0;
    width = width0;
    TR = TR0;
    indexrefraction = indexrefraction0; 
    debug = false;

    inputsurface = new PlanarSurface(posz-width/2.,TR);
    outputsurface = new PlanarSurface(posz+width/2.,TR); 
  }

  void SetDebug(bool a ) {
    debug = a;
    inputsurface->SetDebug(debug);
    outputsurface->SetDebug(debug);
  }

  bool Transport( Ray &ray ){ 

    double outsideindexrefraction = ray.GetIdxR();

    if( ! inputsurface->Transport(ray) ) return false; 

    if( ! inputsurface->Refraction(ray,indexrefraction) ) return false;

    if( ! outputsurface->Transport(ray) ) return false; 

    if( ! outputsurface->Refraction(ray,outsideindexrefraction) ) return false;

    return true;
  }
};

class Lens: public Element {
 private: 
  double indexrefraction;
  double TR;
  double R;
  double posz; // Center of the lens at center
  double width; // Width of the lens at center 
  SphericalSurface *inputsurface; 
  SphericalSurface *outputsurface; 
  bool debug; 

 public: 
    Lens(double posz0, double width0, double R0, double TR0,  double indexrefraction0, int id = 0 ):Element(id) {
   
    posz = posz0;
    width = width0;
    TR = TR0;
    R = R0; 
    indexrefraction = indexrefraction0; 
    debug = false; 

    inputsurface = new SphericalSurface(posz-width/2.+R,R,TR);
    outputsurface = new SphericalSurface(posz+width/2.-R,R,TR); 
  }

  double GetFocus(void){
    return R/2./(indexrefraction-1.);
  }


  void SetDebug(bool a ) {
    debug = a;
    inputsurface->SetDebug(debug);
    outputsurface->SetDebug(debug);
  }

  bool Transport( Ray &ray ){ 

    double outsideindexrefraction = ray.GetIdxR();

    if( ! inputsurface->Transport(ray) ) return false; 

    if( ! inputsurface->Refraction(ray,indexrefraction) ) return false;

    if( ! outputsurface->Transport(ray) ) return false; 

    if( ! outputsurface->Refraction(ray,outsideindexrefraction) ) return false;

    return true;
  }

};


class Aperture: public Element {
 private: 
  double TR; 
  double posz; // Center of the window
  double width; // Width of the window 
  PlanarSurface *inputsurface; 
  PlanarSurface *outputsurface; 
  bool debug; 

 public:
    Aperture(double posz0, double width0, double TR0, int id = 0 ):Element(id){
    posz = posz0;
    width = width0;
    TR = TR0;
    debug = false;
    inputsurface = new PlanarSurface(posz-width/2.,TR);
    outputsurface = new PlanarSurface(posz+width/2.,TR); 
  }

  void SetDebug(bool a ) {
    debug = a;
    inputsurface->SetDebug(debug);
    outputsurface->SetDebug(debug);
  }

  bool Transport( Ray &ray ){ 

    double outsideindexrefraction = ray.GetIdxR();

    if( ! inputsurface->Transport(ray) ) return false; 
    if( ! outputsurface->Transport(ray) ) return false; 

    return true;
  }
};

class Parabola: public Element {
 private: 
  double indexrefraction;
  double D;
  double R;
  double posz; // Center of the lense at center
  ParabolicalSurface *surface; 
  bool debug; 

 public: 
    Parabola(double posz0, double R0, double D0, int id = 0 ):Element(id) {
   
    posz = posz0;
    R = R0; 
    D = D0; 
    debug = false; 

    surface = new ParabolicalSurface(posz,R,D);
  }

  void SetDebug(bool a ) {
    debug = a;
    surface->SetDebug(debug);
  }

  bool Transport( Ray &ray ){ 

    if( ! surface->Transport(ray) ) return false; 

    if( ! surface->Reflexion(ray) ) return false;

    return true;
  }

  double GetFocus(void){
    return surface->GetFocus();
  }

};


class SphericalMirror: public Element {
 private: 
  double TR;
  double R;
  double posz; // Center of the lense at center
  SphericalSurface *surface; 
  bool debug; 

 public: 
    SphericalMirror(double posz0, double R0, double TR0, int id = 0 ):Element(id) {
    posz = posz0;
    TR = TR0;
    R = R0; 
    debug = false; 
    surface = new SphericalSurface(posz-R,R,TR); 
  }

  void SetDebug(bool a ) {
    debug = a;
    surface->SetDebug(debug);
  }

  bool Transport( Ray &ray ){ 

    if( ! surface->Transport(ray) ) return false; 

    if( ! surface->Reflexion(ray) ) return false;

    return true;
  }

  double GetFocus(void){
    return posz-R/2.;
  }

};



#endif 











    
  

