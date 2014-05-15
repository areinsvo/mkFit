#include <cmath>

#include "TMath.h"

#include "Simulation.h"

void setupTrackByToyMC(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, std::vector<Hit>& hits, int& charge, float pt) {

  unsigned int nTotHit = 10;

  //assume beam spot width 1mm in xy and 1cm in z
  pos=SVector3(gRandom->Gaus(0,0.1),gRandom->Gaus(0,0.1),gRandom->Gaus(0,1.0));

  if (charge==0) {
    if (gRandom->Rndm()>0.5) charge = -1;
    else charge = 1;
  }

  float phi = 0.5*TMath::Pi()*gRandom->Rndm(); // make an angle between 0 and pi/2
  float px = pt * cos(phi);
  float py = pt * sin(phi);
  // float px = pt*gRandom->Rndm();
  // float py = sqrt(pt*pt-px*px);
  if (gRandom->Rndm()>0.5) px*=-1.;
  if (gRandom->Rndm()>0.5) py*=-1.;
  float pz = pt*(2.3*(gRandom->Rndm()-0.5));//pz flat between -2*pt and +2*pt
  mom=SVector3(px,py,pz);
  covtrk=ROOT::Math::SMatrixIdentity();
  //initial covariance can be tricky
  for (unsigned int r=0;r<6;++r)
    for (unsigned int c=0;c<6;++c) {
      if (r==c) covtrk(r,c)=1;
      else covtrk(r,c)=0.5;
    }

  //std::cout << "track with p=" << px << " " << py << " " << pz << " pt=" << sqrt(px*px+py*py) << " p=" << sqrt(px*px+py*py+pz*pz) << std::endl;

  float hitposerr = 0.01;//assume 100mum uncertainty in each coordinate
  float k=charge*100./(-0.299792458*3.8);
  float curvature = pt*k;
  float ctgTheta=mom.At(2)/pt;

  TrackState initState;
  initState.parameters=SVector6(pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
  initState.errors=covtrk;
  initState.charge=charge;

  TrackState tmpState = initState;

  //do 4 cm in radius using propagation.h
  for (unsigned int nhit=1;nhit<=nTotHit;++nhit) {
    TrackState propState = propagateHelixToR(tmpState,4.*float(nhit));//radius of 4*nhit
    float hitx = gRandom->Gaus(0,hitposerr)+propState.parameters.At(0);
    float hity = gRandom->Gaus(0,hitposerr)+propState.parameters.At(1);
    //float hity = sqrt((pos.At(0) + k*(px*sinAP-py*(1-cosAP)))*(pos.At(0) + k*(px*sinAP-py*(1-cosAP)))+
    //          	(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))*(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))-
    //	   	        hitx*hitx);//try to get the fixed radius
    float hitz = gRandom->Gaus(0,hitposerr)+propState.parameters.At(2);
    //std::cout << "hit#" << nhit << " " << hitx << " " << hity << " " << hitz << std::endl;
    SVector3 x1(hitx,hity,hitz);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=hitposerr*hitposerr; 
    covx1(1,1)=hitposerr*hitposerr;
    covx1(2,2)=hitposerr*hitposerr;
    Hit hit1(x1,covx1);    
    hits.push_back(hit1);  
    tmpState = propState;
  }
  
  /*
  //do 4 cm along path
  for (unsigned int nhit=1;nhit<=nTotHit;++nhit) {
    float distance = 4.*float(nhit);//~4 cm distance along curvature between each hit
    float angPath = distance/curvature;
    float cosAP=cos(angPath);
    float sinAP=sin(angPath);
    float hitx = gRandom->Gaus(0,hitposerr)+(pos.At(0) + k*(px*sinAP-py*(1-cosAP)));
    float hity = gRandom->Gaus(0,hitposerr)+(pos.At(1) + k*(py*sinAP+px*(1-cosAP)));
    //float hity = sqrt((pos.At(0) + k*(px*sinAP-py*(1-cosAP)))*(pos.At(0) + k*(px*sinAP-py*(1-cosAP)))+
    //          	(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))*(pos.At(1) + k*(py*sinAP+px*(1-cosAP)))-
    //	   	        hitx*hitx);//try to get the fixed radius
    float hitz = gRandom->Gaus(0,hitposerr)+(pos.At(2) + distance*ctgTheta);    
    //std::cout << "hit#" << nhit << " " << hitx << " " << hity << " " << hitz << std::endl;
    SVector3 x1(hitx,hity,hitz);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=hitposerr*hitposerr; 
    covx1(1,1)=hitposerr*hitposerr;
    covx1(2,2)=hitposerr*hitposerr;
    Hit hit1(x1,covx1);    
    hits.push_back(hit1);
  }
  */

}

//temporary... from a dump of cmssw events
void setupTrackByHand(SVector3& pos, SVector3& mom, SMatrixSym66& covtrk, std::vector<Hit>& hits, int& charge, int pt) {
  if (pt==1) {
    charge=-1;
    pos=SVector3(-4.44741,-1.34288,-0.716961);
    mom=SVector3(-0.9199,-0.373767,-0.154805);
    covtrk = ROOT::Math::SMatrixIdentity();
    covtrk(0,0)=0.00653892 ;
    covtrk(0,1)=-0.0104617 ;
    covtrk(0,2)=-0.0135972 ;
    covtrk(0,3)=0.000270623;
    covtrk(0,4)=0.00428196 ;
    covtrk(0,5)=0.00374236 ;
    covtrk(1,1)=0.0284513  ;
    covtrk(1,2)=-0.00652737;
    covtrk(1,3)=-0.00205512;
    covtrk(1,4)=-0.0111767 ;
    covtrk(1,5)=0.000726528;
    covtrk(2,2)=0.096559   ;
    covtrk(2,3)=0.00335384 ;
    covtrk(2,4)=0.00154082 ;
    covtrk(2,5)=-0.0239924 ;
    covtrk(3,3)=0.392652   ;
    covtrk(3,4)=0.151538   ;
    covtrk(3,5)=0.0585285  ;
    covtrk(4,4)=0.106961   ;
    covtrk(4,5)=0.0252197  ;
    covtrk(5,5)=0.0590578  ; 
    SVector3 x1(-4.44738,-1.34304,-0.717461);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=9.8227e-08; 
    covx1(0,1)=-5.57455e-07;
    covx1(0,2)=8.29384e-12;
    covx1(1,1)=3.16365e-06;
    covx1(1,2)=-6.35349e-11;
    covx1(2,2)=1.1557e-05;
    Hit hit1(x1,covx1);
    SVector3 x2(-7.10575,-2.47484,-1.16735);
    SMatrixSym33 covx2 = ROOT::Math::SMatrixIdentity();
    covx2(0,0)=7.94892e-08; 
    covx2(0,1)=-2.44671e-07;
    covx2(0,2)=1.6297e-11;
    covx2(1,1)=7.5311e-07;
    covx2(1,2)=-3.84282e-11;
    covx2(2,2)=1.16604e-05;
    Hit hit2(x2,covx2);
    SVector3 x3(-9.71576,-3.69152,-1.61728);
    SMatrixSym33 covx3 = ROOT::Math::SMatrixIdentity();
    covx3(0,0)=6.39644e-08; 
    covx3(0,1)=-1.63007e-07;
    covx3(0,2)=-2.92856e-11;
    covx3(1,1)=4.15408e-07;
    covx3(1,2)=-3.41337e-11;
    covx3(2,2)=1.16777e-05;
    Hit hit3(x3,covx3);
    hits.push_back(hit1);
    hits.push_back(hit2);
    hits.push_back(hit3);
  } else if (pt==10) {
    charge = 1;
    pos=SVector3(2.94053,3.5543,3.05683);
    mom=SVector3(6.49551,7.58125,-11.6539);
    covtrk = ROOT::Math::SMatrixIdentity();
    covtrk(0,0)=0.00829206  ;
    covtrk(0,1)=-0.000988223;
    covtrk(0,2)=0.00397884  ;
    covtrk(0,3)=-0.174796   ;
    covtrk(0,4)=-0.184843   ;
    covtrk(0,5)=0.277992    ;
    covtrk(1,1)=0.0077582   ;
    covtrk(1,2)=0.00449615  ;
    covtrk(1,3)=0.141826    ;
    covtrk(1,4)=0.148706    ;
    covtrk(1,5)=-0.261275   ;
    covtrk(2,2)=0.00514256  ;
    covtrk(2,3)=-0.00516294 ;
    covtrk(2,4)=-0.00628755 ;
    covtrk(2,5)=-0.0150243  ;
    covtrk(3,3)=77.1274     ;
    covtrk(3,4)=88.3183     ;
    covtrk(3,5)=-136.784    ;
    covtrk(4,4)=101.395     ;
    covtrk(4,5)=-156.858    ;
    covtrk(5,5)=243.153     ;
    SVector3 x1(2.94069,3.55417,3.05572);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=5.78482e-07 ; 
    covx1(0,1)=-4.85419e-07;
    covx1(0,2)=-2.09565e-11;
    covx1(1,1)=4.07328e-07 ;
    covx1(1,2)=-3.08018e-11;
    covx1(2,2)=4.39754e-06 ;
    Hit hit1(x1,covx1);
    SVector3 x2(4.53878,5.41477,0.193781);
    SMatrixSym33 covx2 = ROOT::Math::SMatrixIdentity();
    covx2(0,0)=5.50755e-07 ; 
    covx2(0,1)=-4.00035e-07;
    covx2(0,2)=-1.61696e-11;
    covx2(1,1)=2.90562e-07 ;
    covx2(1,2)=3.45213e-13 ;
    covx2(2,2)=2.74192e-06 ;
    Hit hit2(x2,covx2);
    SVector3 x3(6.41319,7.58286,-3.1512);
    SMatrixSym33 covx3 = ROOT::Math::SMatrixIdentity();
    covx3(0,0)=4.57591e-07 ; 
    covx3(0,1)=-4.24647e-07;
    covx3(0,2)=-2.75564e-11;
    covx3(1,1)=3.94076e-07 ;
    covx3(1,2)=-2.13072e-11;
    covx3(2,2)=4.39821e-06 ;
    Hit hit3(x3,covx3);
    hits.push_back(hit1);
    hits.push_back(hit2);
    hits.push_back(hit3); 
  } else if (pt==100) {
    charge=-1;
    pos=SVector3(0.556433,-4.13765,1.17947);
    mom=SVector3(6.93565,-99.8987,106.122);
    covtrk = ROOT::Math::SMatrixIdentity();
    covtrk(0,0)=0.00668929   ;
    covtrk(0,1)=-0.000405335 ;
    covtrk(0,2)=-0.00081875  ;
    covtrk(0,3)=-0.654305    ;
    covtrk(0,4)=8.77271      ;
    covtrk(0,5)=-9.28359     ;
    covtrk(1,1)=0.00585989   ;
    covtrk(1,2)=0.00554276   ;
    covtrk(1,3)=-0.0221493   ;
    covtrk(1,4)=0.336286     ;
    covtrk(1,5)=-0.533675    ;
    covtrk(2,2)=0.00527125   ;
    covtrk(2,3)=0.0219121    ;
    covtrk(2,4)=-0.25678     ;
    covtrk(2,5)=0.104355     ;
    covtrk(3,3)=365.307      ;
    covtrk(3,4)=-5081.05     ;
    covtrk(3,5)=5397.81      ;
    covtrk(4,4)=70741.4      ;
    covtrk(4,5)=-75154.2     ;
    covtrk(5,5)=79853.4      ; 
    SVector3 x1(0.555814,-4.13765,1.18106);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=7.89828e-07; 
    covx1(0,1)=1.71343e-10;
    covx1(0,2)=-3.65348e-11;
    covx1(1,1)=3.75336e-14;
    covx1(1,2)=3.71905e-11;
    covx1(2,2)=3.81485e-06;
    Hit hit1(x1,covx1);
    SVector3 x2(0.759732,-7.0482,4.26977);
    SMatrixSym33 covx2 = ROOT::Math::SMatrixIdentity();
    covx2(0,0)=1.16538e-06; 
    covx2(0,1)=3.15313e-10;
    covx2(0,2)=-9.24943e-12;
    covx2(1,1)=8.59631e-14;
    covx2(1,2)=-4.24538e-11;
    covx2(2,2)=2.7723e-06;
    Hit hit2(x2,covx2);
    SVector3 x3(0.790956,-7.52219,4.77523);
    SMatrixSym33 covx3 = ROOT::Math::SMatrixIdentity();
    covx3(0,0)=5.06206e-06; 
    covx3(0,1)=1.07417e-06;
    covx3(0,2)=-8.4071e-11;
    covx3(1,1)=2.27938e-07;
    covx3(1,2)=1.37342e-10;
    covx3(2,2)=1.521e-05  ;
    Hit hit3(x3,covx3);
    SVector3 x4(0.990696,-10.3606,7.79052);
    SMatrixSym33 covx4 = ROOT::Math::SMatrixIdentity();
    covx4(0,0)=1.32161e-06 ; 
    covx4(0,1)=1.98982e-07 ;
    covx4(0,2)=-1.69703e-12;
    covx4(1,1)=2.99588e-08 ;
    covx4(1,2)=-1.4976e-11 ;
    covx4(2,2)=3.81449e-06 ;
    Hit hit4(x4,covx4);
    hits.push_back(hit1);
    hits.push_back(hit2);
    hits.push_back(hit3);
    hits.push_back(hit4);
  } else {
    charge = 1;
    pos=SVector3(-3.97545,-1.39023,-3.20927);
    mom=SVector3(-822.213,-348.139,-765.718);
    covtrk = ROOT::Math::SMatrixIdentity();
    covtrk(0,0)=0.00339944  ;
    covtrk(0,1)=-0.000129504;
    covtrk(0,2)=-0.00359138 ;
    covtrk(0,3)=-86.3911    ;
    covtrk(0,4)=-36.5052    ;
    covtrk(0,5)=-79.7691    ;
    covtrk(1,1)=0.00386787 ;
    covtrk(1,2)=-0.0016195 ;
    covtrk(1,3)=154.307    ;
    covtrk(1,4)=65.1894    ;
    covtrk(1,5)=143.913    ;
    covtrk(2,2)=0.00459267;
    covtrk(2,2)=22.6083   ;
    covtrk(2,2)=9.55969   ;
    covtrk(2,2)=20.2235   ;
    covtrk(3,3)=4.59221e+07;
    covtrk(3,3)=1.94154e+07;
    covtrk(3,3)=4.27286e+07;
    covtrk(4,4)=8.20865e+06;
    covtrk(4,4)=1.80653e+07;
    covtrk(5,5)=3.9758e+07;
    SVector3 x1(-3.97592,-1.38942,-3.21347);
    SMatrixSym33 covx1 = ROOT::Math::SMatrixIdentity();
    covx1(0,0)=2.8478e-07  ;
    covx1(0,1)=-4.93169e-07;
    covx1(0,2)=-1.68421e-11;
    covx1(1,1)=8.54048e-07 ;
    covx1(1,2)=-2.08965e-11;
    covx1(2,2)=2.82947e-06 ;
    Hit hit1(x1,covx1);
    SVector3 x2(-4.40679,-1.57342,-3.60635);
    SMatrixSym33 covx2 = ROOT::Math::SMatrixIdentity();
    covx2(0,0)=4.19391e-09 ;
    covx2(0,1)=-2.38012e-08;
    covx2(0,2)=3.16473e-12 ;
    covx2(1,1)=1.35076e-07 ;
    covx2(1,2)=-2.23185e-11;
    covx2(2,2)=3.05883e-06 ;
    Hit hit2(x2,covx2);
    SVector3 x3(-7.03704,-2.68632,-6.06063);
    SMatrixSym33 covx3 = ROOT::Math::SMatrixIdentity();
    covx3(0,0)=1.26699e-07 ;
    covx3(0,1)=-3.89985e-07;
    covx3(0,2)=2.77271e-12 ;
    covx3(1,1)=1.20039e-06 ;
    covx3(1,2)=-5.58307e-12;
    covx3(2,2)=2.93281e-06 ;
    Hit hit3(x3,covx3);
    SVector3 x4(-9.67263,-3.80139,-8.51372);
    SMatrixSym33 covx4 = ROOT::Math::SMatrixIdentity();
    covx4(0,0)=2.05148e-06 ;
    covx4(0,1)=-5.22726e-06;
    covx4(0,2)=-1.94228e-10;
    covx4(1,1)=1.33193e-05 ;
    covx4(1,2)=-1.02983e-10;
    covx4(2,2)=1.88032e-05 ;
    Hit hit4(x4,covx4);
    hits.push_back(hit1);
    hits.push_back(hit2);
    hits.push_back(hit3);
    hits.push_back(hit4);
  }  
}

