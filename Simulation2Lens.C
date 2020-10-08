#include "OpticSystemTransport.h" 
#include "TRandom.h"
#include "TF1.h"  
#include "TH1D.h" 
#include "TH2D.h" 
#include "TH3D.h" 
#include "TProfile.h"
#include "TCanvas.h" 

double Simulation2Lens(double size = 1.,double lensposition = 100., double sensorpos = 126.,bool debug = false){

  Disc  w(30.,8.,50.,1.458,1); 
  Aperture  a(lensposition+17.,1.,15.,2); 
  //  Lens lens(lensposition,16.,74.78,49./2.,1.78472,3);
  // Lens lens2(lensposition+16.+0.1,16.,74.78,49./2.,1.78472,3);

  double d = 0;
  
  #if 0
  Lens lens(lensposition,6.5,45.61,15.,1.833,3);
  #else
  Lens lens(lensposition,16.,74.78,50./2.,1.833,3);
  #endif

  #if 0
  Lens lens2(lensposition+6.5+d,6.5,45.61,15.,1.833,3);
  #else 
  Lens lens2(lensposition+16.+d,16.,45.61,15.,1.833,3);
  #endif

  lens.SetDebug(debug);
  lens2.SetDebug(debug);
  a.SetDebug(debug);  
  w.SetDebug(debug); 

  // You can recompute the matrices for different distances by (i.e. if you want to scan parameters): 
  //
  // opt.setdistances(d_source_silica,d_silica_lens,d_lens_sensor)
  // opt.ComputeMatrices(); 

  TRandom r;

  TH2D     *h2x = new TH2D("h2x","  ",100,-size,size,100,-8.,8.); 
  TH2D     *h2y = new TH2D("h2y","  ",100,-size,size,100,-8.,8.); 
  TProfile *hpx = new TProfile("hpx","  ",100,-size,size,-8.,8.,"S"); 
  TProfile *hpy = new TProfile("hpy","  ",100,-size,size,-8.,8.,"S"); 

  TH2D     *hsensornolens = new TH2D("hsensornolens","  ",16,-8.,8.,16,-8.,8.); 
  TH2D     *hsensorlens = new TH2D("hsensorlens","  ",16,-8.,8.,16,-8.,8.); 

  TH2D     *hwindow = new TH2D("hwindow","  ",100,-50.,50.,100,-50.,50.); 

  TH1D *haberr = new TH1D("haberr","",1000,0.,size/10.);  // RMS is always positive

  TH3D *testspherical = new TH3D("testspherical","",100,-1.,1.,100,-1.,1.,100,-1.,1.); 
  TH3D *hdxr = new TH3D("hdxr","",100,0.,25.,100,-size,size,100,-size/5.,size/5.); 

  double D = 100.;
  double d1 = 30.;

  double anglimit = TMath::ATan((D/2.)/d1);
 
  std::cout << anglimit << "  " << TMath::Cos(anglimit) << std::endl;

  int ntot = 1e+7;

  int sensornolens = 0;

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

    testspherical->Fill(vx,vy,vz);

    Ray ray(x,y,z,vx,vy,vz,1.); 
    
    // The angle is wrt the lens surface vector (perpendicular to surface). 

    double thx = (vx/vz);
    double thy = (vy/vz);
   
    if( TMath::Sqrt((x+thx*d1)*(x+thx*d1)+(y+thy*d1)*(y+thy*d1)) < D/2. ) 
      hwindow->Fill(x+thx*d1,y+thy*d1);

    if( TMath::Abs(x+thx*d1)<8. &&  TMath::Abs(y+thy*d1)< 8. ) {
      hsensornolens->Fill(x+thx*d1,y+thy*d1);
    } 

    double xf; 
    double yf; 
    double thxf;
    double thyf;
                
    if( !w.Transport(ray) ) continue;
    if( !lens.Transport(ray) ) continue;
    if( !lens2.Transport(ray) ) continue;
    //if( !a.Transport(ray) ) continue;
    double Ratlens = TMath::Sqrt(ray.GetX()*ray.GetX()+ray.GetY()*ray.GetY());

    // Image 
    ray.Transport(sensorpos);
    //ray.Print();
      
    if( TMath::Abs(ray.GetX())<8. &&  TMath::Abs(ray.GetY())< 8. ) {
      hsensorlens->Fill(ray.GetX(),ray.GetY());

    // scaled image. 

      h2x->Fill(x,ray.GetX()); 
      hpx->Fill(x,ray.GetX()); 
      h2y->Fill(y,ray.GetY()); 
      hpy->Fill(y,ray.GetY()); 
    }

    hdxr->Fill(Ratlens,x,ray.GetX()); 

  }

  h2x->FitSlicesX(); // This fits gaussians in slices of Y
  h2y->FitSlicesX(); // This fits gaussians in slices of Y 

  TH1D *h2xmean = (TH1D*) gROOT->FindObject("h2x_1"); // This gets the histogram of  mean of the fitted gaussian as function of reco points.
  TH1D *h2ymean = (TH1D*) gROOT->FindObject("h2y_1");

  h2xmean->Fit("pol1","","",-7.,7.);
  h2ymean->Fit("pol1","","",-7.,7.);

  double scale = ((h2xmean->GetFunction("pol1")->GetParameter(1))+(h2ymean->GetFunction("pol1")->GetParameter(1))) /2.;
  
  // I define the aberration as the sigma of points in the true position for a given position in the detector. 
  
  TH1D *h2xsigma = (TH1D*) gROOT->FindObject("h2x_2"); // This gets the histogram of  sigmas of the fitted gaussian as function of reco points. 


  h2xsigma->Fit("pol0","","",-3.,3.);   // To obtain the aberration, this is mainly the spread of true points 

  double aberr = (h2xsigma->GetFunction("pol0")->GetParameter(0)); 

  scale = TMath::Abs(scale); 
  
  TCanvas *c = new TCanvas("c","",800,900);
  c->Divide(2,4); 
  c->cd(1); 
  h2x->Draw("colz"); 
  hpx->Fit("pol1","","same"); 
 
  for(int i = 0; i < 100; i++ ) { 
    haberr->Fill(hpx->GetBinError(i+1)); 
  }
  c->cd(2); 
  h2y->Draw("colz"); 
  hpy->Fit("pol1","","same"); 
 
  for(int i = 0; i < 100; i++ ) { 
    haberr->Fill(hpy->GetBinError(i+1)); 
  }
  
  c->cd(3); 

  hsensorlens->Draw("colz"); 

  c->cd(4);
  hsensornolens->Draw("colz");

  c->cd(5); 
  hwindow->Draw("colz");

  c->cd(6);
  hdxr->Draw("box");

  c->cd(7);
  h2xmean->Draw();
  h2xmean->GetYaxis()->SetRangeUser(-7.*scale*1.1,7.*scale*1.1); 

  c->cd(8);
  h2ymean->Draw();
  h2ymean->GetYaxis()->SetRangeUser(-7.*scale*1.1,7.*scale*1.1); 
  c->Update(); 

  TCanvas *c1 = new TCanvas("c1","",400,400);
  testspherical->Draw("colz");
  c1->Update();

  std::cout <<std::endl<<std::endl; 

 
  
  std::cout << " ------------------------------------------------------------ " << std::endl;
  std::cout << " Scale factor  " <<  scale << std::endl; 
  std::cout << " Mean dispersion at sensor " <<  aberr/scale << std::endl; 
  std::cout << " Mean dispersion at image " <<  aberr << std::endl; 
  std::cout << " Fraction of photons reaching the sensor " << std::endl; 
  std::cout << "       With lens " << (double)hsensorlens->GetEntries()/(double)ntot<< std::endl;
  std::cout << "       With no lense " << (double)hsensornolens->GetEntries()/(double)ntot << std::endl;
  std::cout << " Maximum number for " <<  hsensorlens->GetMaximum() << std::endl; 
  std::cout << " ------------------------------------------------------------------- " << std::endl;

  TCanvas *ccc = new TCanvas("ccc");
  h2xsigma->Draw();
  ccc->Update(); 
  
  return aberr;
}

