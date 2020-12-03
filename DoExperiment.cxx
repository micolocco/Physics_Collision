
#include "Detector.h"
#include "Riostream.h"
#include "string.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TStopwatch.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"

TH1F *ExtractDist(TH1F *); //function to extract the pseudorapidity distribution in [-2,2]

//varMS=variable for multiple scattering. 1=multiple scattering on, 0=multiple scattering off
//ncoll=number of collisions

void DoExperiment(int ncoll, bool varMS)
{
  TStopwatch timer;
  timer.Start();

  Detector dect;
  Vertex vertice;
  Track traccia;
  //extraction of multiplicity and pseudorapidity distributions
  TFile F("kinem.root");
  TH1F *distMult = (TH1F *)F.Get("hmul");
  TH1F *distHeta = (TH1F *)F.Get("heta");
  distHeta->SetDirectory(0);
  TH1F *heta2 = ExtractDist(distHeta);

  // opening of output file
  TFile hfile("htree.root", "RECREATE");
  // TTree declaration
  TTree *tree = new TTree("T", "T");
  TClonesArray *ptrhits[2]; //array of TClonesArray to save the hits on each layer
  for (int iLayer = 0; iLayer < 2; iLayer++)
  {
    ptrhits[iLayer] = new TClonesArray("Point", 100);
    ptrhits[iLayer] = dect.GetPtrhits(iLayer);
  }

  tree->Branch("VertMult", "Vertex", &vertice);
  tree->Branch("HitsLayer1", &ptrhits[0]);
  tree->Branch("HitsLayer2", &ptrhits[1]);

  for (int iColl = 0; iColl < ncoll; iColl++)
  {
    dect.Simulation(varMS, vertice, traccia, distMult, heta2);
    //dect.Simulation(varMS, vertice, traccia, heta2); //in case of fixed or uniform multiplicity
    if (iColl % 10000 == 0)
      cout << "Collision n° " << iColl << endl;
    /*
    printf("COLLISIONE N° %d - moltepl: %d\n", iColl, vertice.GetMult());
    printf("VERTICE: x= %f ; y= %f; z= %f \n", vertice.GetX(), vertice.GetY(), vertice.GetZ());
    cout << "Punti di noise: " << nNoise << endl;
    
    // Debug
    for (int iLayer = 0; iLayer < 2; iLayer++)
    {      
      cout << "-----------" << dect.GetLName(iLayer + 1) << "----------" << endl;
      int tot = vertice.GetMult() + nNoise;
      for (int j = 0; j < dect.GetPtrhits(iLayer)->GetEntries(); j++)
      {
        Point *tst = (Point *)dect.GetPtrhits(iLayer)->At(j); //cast di ptrhits per avere un array di punti; sta riempiendo l'array
        cout << "Hit " << j << ") x, y, z = " << tst->GetX() << "; " << tst->GetY() << "; " << tst->GetZ() << endl;
      }
      cout << endl;
    }
    // fine del debug*/

    tree->Fill();
    for (int iLayer = 0; iLayer < 2; iLayer++)
    {
      dect.GetPtrhits(iLayer)->Clear();
    }
  }

  hfile.Write();

  delete tree;
  delete ptrhits[0];
  delete ptrhits[1];
  delete distMult;
  delete distHeta;
  F.Close();
  hfile.Close();
  timer.Stop();
  timer.Print();
}

TH1F *ExtractDist(TH1F *disteta)
{

  TAxis *xa = disteta->GetXaxis();
  int b1 = xa->FindBin(-2.);
  int b2 = xa->FindBin(2.);
  double xlow = xa->GetBinLowEdge(b1);
  double xhig = xa->GetBinUpEdge(b2);
  int nobins = b2 - b1 + 1;
  TH1F *heta2 = new TH1F("heta2", "#eta distribution 2", nobins, xlow, xhig);
  int j = 1;
  for (int i = b1; i <= b2; i++)
    heta2->SetBinContent(j++, disteta->GetBinContent(i));
  return heta2;
}
