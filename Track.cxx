#include "TAxis.h"
#include "Riostream.h"
#include "TMath.h"
#include "Track.h"
#include "TRandom3.h"
#include "TClonesArray.h"

ClassImp(Track)

    Track::Track(double thetaValue, double phiValue)
{
  fTheta = thetaValue;
  fPhi = phiValue;
}
void Track::SetFixedPhi(double phiValue)
{
  fPhi = phiValue;
}
void Track::SetFixedTheta(double thetaValue)
{
  fTheta = thetaValue;
}

//Set functions for randomic direction
void Track::SetRandomPhi()
{
  double a = gRandom->Rndm();
  fPhi = 2 * TMath::Pi() * a;
}
void Track::SetRandomPhiP()
{
  double a = gRandom->Rndm();
  fPhiP = 2 * TMath::Pi() * a;
};
void Track::SetGaussThetaP()
{
  fThetaP = gRandom->Gaus(0., 0.001); //rad
}

//Get functions
double Track::GetPhi() const
{
  return fPhi;
}
double Track::GetTheta() const
{
  return fTheta;
}
double Track::GetPhiP() const
{
  return fPhiP;
}
double Track::GetThetaP() const
{
  return fThetaP;
}

//pseudorapidity extraction from a distribution
void Track::SetRandomTheta(TH1F *heta2)
{

  double rndmheta = heta2->GetRandom();
  double n = exp(-rndmheta);
  fTheta = 2 * TMath::ATan(n);
}

//Print functions
void Track::PrintPhi() const
{
  cout << "phi (radianti) =" << fPhi << endl;
}
void Track::PrintTheta() const
{
  cout << "theta (radianti) =" << fTheta << endl;
}
void Track::PrintPhiP() const
{
  cout << "phi (radianti) =" << fPhiP << endl;
}
void Track::PrintThetaP() const
{
  cout << "theta (radianti) =" << fThetaP << endl;
}

//---------------------------------------------------
void Track::Rotate(double *cd)
{
  double mr[3][3];

  mr[0][0] = -cd[1] / TMath::Sqrt(1 - cd[2] * cd[2]);
  mr[1][0] = cd[0] / TMath::Sqrt(1 - cd[2] * cd[2]);
  mr[2][0] = 0.;
  mr[0][1] = -(cd[0] * cd[2]) / TMath::Sqrt(1 - cd[2] * cd[2]);
  mr[1][1] = -(cd[1] * cd[2]) / TMath::Sqrt(1 - cd[2] * cd[2]);
  mr[2][1] = TMath::Sqrt(1 - cd[2] * cd[2]);
  mr[0][2] = cd[0];
  mr[1][2] = cd[1];
  mr[2][2] = cd[2];

  double cdp[3];
  cdp[0] = TMath::Sin(fThetaP) * TMath::Cos(fPhiP);
  cdp[1] = TMath::Sin(fThetaP) * TMath::Sin(fPhiP);
  cdp[2] = TMath::Cos(fThetaP);
  for (int i = 0; i < 3; i++)
  {
    cd[i] = 0.;
    for (int j = 0; j < 3; j++)
    {
      cd[i] += mr[i][j] * cdp[j];
    }
  }
}
//--------------------------------------------------------

void Track::SetDirection(Vertex &vertice, TH1F *heta2, TClonesArray *ptrhits[], const int NLayer, double *R, double H, bool varMS) //varMS=1 (multiple scattering on) varMS=0 (multiple scattering off)
{
  double noise[3]; //vector for noise coordinates

  //tools for particle propagation
  double cVect[3]; //coefficient vector
  double disc, t;
  bool success; //variable to manage the hit saving per layer

  //TClonesArray to save hits on layer 1,2 after smearing
  TClonesArray &hitsSmearing1 = *ptrhits[0];
  TClonesArray &hitsSmearing2 = *ptrhits[1];
  int index[2] = {0, 0}; //to count the elements in hitsSmearing1,hitsSmearing2

  for (int j = 0; j < vertice.GetMult(); j++)
  {
    SetRandomPhi();
    SetRandomTheta(heta2);
    cVect[0] = TMath::Sin(GetTheta()) * TMath::Cos(GetPhi());
    cVect[1] = TMath::Sin(GetTheta()) * TMath::Sin(GetPhi());
    cVect[2] = TMath::Cos(GetTheta());
    success = true;

    for (int iLayer = 0; iLayer < NLayer; iLayer++)
    {
      if (success)
      {
        disc = vertice.computeDisc(cVect, R[iLayer]);
        t = vertice.computeT(cVect, disc);
        success = (vertice.check(disc, t, cVect[2], H));
        if (success)
        {
          //hit coordinates
          double xH = vertice.GetX() + cVect[0] * t;
          double yH = vertice.GetY() + cVect[1] * t;
          double zH = vertice.GetZ() + cVect[2] * t;

          //------------------MULTIPLE SCATTERING------
          // varMs=1: multiple scattering on
          if (varMS && iLayer != 2)
            {
              SetRandomPhiP();
              SetGaussThetaP();
              Rotate(cVect);
              disc = vertice.computeDisc(cVect, R[iLayer]);
              t = vertice.computeT(cVect, disc);
            }
          if (vertice.check(disc, t, cVect[2], H) && iLayer != 0) //smearing
            {
              double phiS = (gRandom->Gaus(0, 0.003)) / R[iLayer];
              double zS = zH;
              do //check on z smeared
              {
                zS = zH + gRandom->Gaus(0, 0.012);
              } while (zS < -H / 2 || zS > H / 2);
              if (iLayer == 1)
              {
                new (hitsSmearing1[index[iLayer - 1]])
                    Point(xH * TMath::Cos(phiS) - yH * TMath::Sin(phiS), xH * TMath::Sin(phiS) + yH * TMath::Cos(phiS), zS);
                //cout<<"X"<<xH * TMath::Cos(phiS) - yH * TMath::Sin(phiS);
              }
              else
              {
                new (hitsSmearing2[index[iLayer - 1]])
                    Point(xH * TMath::Cos(phiS) - yH * TMath::Sin(phiS), xH * TMath::Sin(phiS) + yH * TMath::Cos(phiS), zS);
                //cout<<"X"<<xH * TMath::Cos(phiS) - yH * TMath::Sin(phiS);
              }
              index[iLayer - 1]++;
             }
        }
        success = (vertice.check(disc, t, cVect[2], H));
      }
      //noise added
      if (j == vertice.GetMult() - 1)
      {
        for (int iLayer = 1; iLayer < 3; iLayer++)
        {
          for (int iNoise = 0; iNoise < nNoise; iNoise++)
          {
            Noise(R[iLayer], H, noise);
            if (iLayer == 1)
              new (hitsSmearing1[index[iLayer - 1] + iNoise]) Point(noise[0], noise[1], noise[2]);
            if (iLayer == 2)
              new (hitsSmearing2[index[iLayer - 1] + iNoise]) Point(noise[0], noise[1], noise[2]);
          }
        }
      }
    }
    //saving in case of acceptance study
    /*if (index[0]!=index[1]){ 
       hitsSmearing1.RemoveAt(index[0]);
      index[0]--;
   }*/
    // cout<<"index lay1"<< index[0]<<"  index lay2 "<<index[1]<<endl;

    /*if (j == vertice.GetMult() - 1&&!success){
      hitsSmearing1.RemoveAt(index[0]);
   }*/
  }
}

void Track::Noise(double R, double H, double *noise)
{
  double phi = gRandom->Uniform(0, 2 * TMath::Pi());
  noise[0] = R * TMath::Cos(phi);
  noise[1] = R * TMath::Sin(phi);
  noise[2] = gRandom->Uniform(-H/2 , H/2 );
}
