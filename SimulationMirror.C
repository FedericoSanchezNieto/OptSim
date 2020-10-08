#include "OpticSystemTransport.h" 
#include "TRandom.h"
#include "TF1.h"  
#include "TH1D.h" 
#include "TH2D.h" 
#include "TH3D.h" 
#include "TProfile.h"
#include "TCanvas.h" 

double SimulationMirror(double pos = 500.,double radious = 400., double off = 0., bool debug = false){
  SphericalMirror parab(pos,radious,radious);

  parab.SetDebug(debug); 

  // You can recompute the matrices for different distances by (i.e. if you want to scan parameters): 
  //
  // opt.setdistances(d_source_silica,d_silica_lens,d_lens_sensor)
  // opt.ComputeMatrices(); 

  TRandom r;


  TH2D *rz = new TH2D("rz","",1000,0.,500,1000,0,500); 

  int ntot = 1e+6;
  double size = 20.;
  double anglimit = TMath::ATan(radious/pos);

  double ss = 8.0;

  TH2D *map = new TH2D("map","",16,-ss,ss,16,-ss,ss); 

  for(int i = 0; i < ntot;i++) {

    double x = r.Uniform(-size,size);
    double y = x;

    double z = 0.;
    double costheta = r.Uniform(TMath::Cos(anglimit),1.); 
    double phi = r.Uniform(0.,2.*3.141592); 
    double sintheta = TMath::Sqrt(1.-costheta*costheta);

    double vx = sintheta*TMath::Cos(phi); 
    double vy = sintheta*TMath::Sin(phi); 
    double vz = costheta; 

    Ray ray(x,y,z,vx,vy,vz,1.); 
    
    if( !parab.Transport(ray) ) continue;

    rz->Fill(ray.GetZ(),sqrt(ray.GetY()*ray.GetY()+ray.GetX()*ray.GetX()));

    ray.Transport(parab.GetFocus()+off); 

    if( TMath::Abs(ray.GetX()) < 8. && TMath::Abs(ray.GetY()) < 8. )
      map->Fill(ray.GetX(),ray.GetY()); 

  }

  map->Draw("colz"); 

  //rz->Draw("colz"); 

  return 0.;
}

