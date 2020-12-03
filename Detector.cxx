#include "Detector.h"
#include <string.h>
#include "Riostream.h"
#include "TFile.h"
#include "Track.h"
#include "Vertex.h"

ClassImp(Detector)

    Detector::Detector()
{
  fH = length;
  for (int i = 0; i < nLayer; i++)
    fR[i] = radius[i];
  for (int i = 0; i < nLayer; i++)
    fNames[i] = names[i];
  for (int i = 0; i < 2; i++)
    fPtrhits[i] = new TClonesArray("Point", 100); //declaration of type
}

void Detector::SetH(double H)
{
  fH = H;
}

void Detector::SetRLayer(double *R)
{
  for (int i = 0; i < nLayer; i++)
    fR[i] = R[i];
}

void Detector::SetNames(string *names)
{
  for (int i = 0; i < nLayer; i++)
    fNames[i] = names[i];
}

void Detector::PrintDet() const
{
  cout << "----Detector features----" << endl;
  cout << "Length= " << fH << " cm" << endl;
  cout << "Radius: { ";
  for (int iLayer = 0; iLayer < nLayer; iLayer++)
  {
    cout << fNames[iLayer] << ":" << fR[iLayer] << " ";
  }
  cout << "}"
       << " cm" << endl;
}

void Detector::Simulation(bool varMS, Vertex &vertice, Track &traccia, TH1F *distMult, TH1F *heta2){
//void Detector::Simulation(bool varMS, Vertex &vertice, Track &traccia, TH1F *heta2){
  vertice.SetMult(distMult);
 // vertice.SetMult();
  vertice.SetRandomX();
  vertice.SetRandomY();
  vertice.SetRandomZ();
  traccia.SetDirection(vertice, heta2, fPtrhits, nLayer, fR, fH, varMS);
}

Detector::~Detector()
{
}
