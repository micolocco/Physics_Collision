#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Detector.h"
#include <TCanvas.h>
#include <TH1.h>
#include <vector>
#include <bits/stdc++.h>
#include "TNtuple.h"
#include "TStopwatch.h"
#include "TLeaf.h"

int NBin = 250;
void ReconstructVertex()
{
  TStopwatch timer;
  timer.Start();
  Vertex *vertice = new Vertex();
  Detector dect;
  double sum, count, rapp, zMax, zUp, zLow, m, q, zTry, zMean = 0;

  // declaration of TClonesArray type
  TClonesArray *hitsLayer[2];
  hitsLayer[0] = new TClonesArray("Point", 100);
  hitsLayer[1] = new TClonesArray("Point", 100);

  //opening of input file
  TFile inputFile("htree.root");
  //reading of TTree
  TTree *tree = (TTree *)inputFile.Get("T");
  TBranch *b1 = tree->GetBranch("VertMult");
  TBranch *b2 = tree->GetBranch("HitsLayer1");
  TBranch *b3 = tree->GetBranch("HitsLayer2");
  // Definition of branch address
  b1->SetAddress(&vertice);
  b2->SetAddress(&hitsLayer[0]);
  b3->SetAddress(&hitsLayer[1]);
  //opening of output file
  TFile outputFile("zVertez.root", "recreate");

  vector<double> storage;                      //vector to store the attempts for reconstructing z vertex by association
  TNtuple zVect("zVect", "Zvect", "zVertice"); //TNtupla to save the reconstructed zeta
  //TH1D *hPhi=new TH1D("Sen", "Sin(#delta#phi)", 100, -0.07,0.07); //histogram for acceptance study

  //loop on tree
  for (int ev = 0; ev < tree->GetEntries(); ev++)
  {
    tree->GetEvent(ev);
    if (ev % 10000 == 0)
      cout << "Reconstruction of event " << ev << endl;
    TH1D *h1 = new TH1D("zVertex", Form("zVertex-%i_%f", vertice->GetMult(), vertice->GetZ()), NBin, -25, 25); //histogram to evaluate the most probable z
    h1->GetXaxis()->SetTitle("z [cm]");

    /*//debug 
    //NOTE: iLayer=0 is layer1, iLayer=1 is layer2
    for (int iLayer = 0; iLayer < 2; iLayer++)
    {
      int num = hitsLayer[iLayer]->GetEntries();
      cout << "Numero di elementi nel TClonesArray (hits Layer" << iLayer + 1 << ") " << num << endl;
      for (int j = 0; j < num; j++)
      {
        Point*tst = (Point*)hitsLayer[iLayer]->At(j);
        cout << "hit" << j << ") x, y, z = " << tst->GetX() << "; " << tst->GetY() << "; " << tst->GetZ() << endl;
      }
    }
    //end debug*/

    //small acceptance study
    //int molt=vertice->GetMult();
    /*for (int iPart=0; iPart<hitsLayer[1]->GetEntries();iPart++){
      Point*tst1 = (Point*)hitsLayer[0]->At(iPart);
      Point*tst2 = (Point*)hitsLayer[1]->At(iPart);
      double rapp=(tst1->GetY()*tst2->GetX()-tst1->GetX()*tst2->GetY())/(dect.GetRLayer(1)*dect.GetRLayer(2));
      hPhi->Fill(rapp);
          }*/

    for (int iNum = 0; iNum < hitsLayer[0]->GetEntries(); iNum++)
    {
      Point *p1 = (Point *)hitsLayer[0]->At(iNum);
      for (int jNum = 0; jNum < hitsLayer[1]->GetEntries(); jNum++)
      {
        Point *p2 = (Point *)hitsLayer[1]->At(jNum);
        rapp = (p1->GetY() * p2->GetX() - p1->GetX() * p2->GetY()) / (dect.GetRLayer(1) * dect.GetRLayer(2)); // azimut difference computing

        if (TMath::Abs(rapp) <= 0.0052) //small azimut difference condition
        {
          m = (dect.GetRLayer(1) - dect.GetRLayer(2)) / (p1->GetZ() - p2->GetZ());
          q = dect.GetRLayer(1) - m * p1->GetZ();
          zTry = -q / m; //computing of z
          if (TMath::Abs(zTry) < 24.9)
          {
            h1->Fill(zTry);
            storage.push_back(zTry);
           }
        }
      }
    }

    sort(storage.begin(), storage.end());
    //opening of a window around the zTry distribution peak satisfying some conditions
    // if ((h1->GetBinContent(h1->GetMaximumBin() - 1) != 0 || h1->GetBinContent(h1->GetMaximumBin() + 1) != 0) || h1->GetBinContent(h1->GetMaximumBin()) > 1)
    if ((h1->GetBinContent(h1->GetMaximumBin()) > 1) || (h1->GetBinContent(h1->GetMaximumBin()) == 1 && (h1->GetBinContent(h1->GetMaximumBin() - 1) != 0 || h1->GetBinContent(h1->GetMaximumBin() + 1) != 0 || h1->GetEntries() == 1)))

    {
      zMax = h1->GetBinCenter(h1->GetMaximumBin());
      zUp = zMax + 0.5;
      zLow = zMax - 0.5;
      //mean in the window [zMax-0.5,zMax+0.5]cm
      sum = 0.;
      count = 0.;
      for (int iSearch = 0; iSearch < (int)storage.size(); iSearch++)
      {
        if (storage.at(iSearch) <= zUp && storage.at(iSearch) >= zLow)
        {
          sum += storage.at(iSearch);
          count++;
        }
      }
      zMean = sum / count;
      zVect.Fill(zMean);
    }
    else
    {
      //if(ev<10000)h1->Write();
      zVect.Fill(50); //TNTupla is filled in z=50cm if conditions on the peak are not satisfied
    }
    storage.clear();
    h1->Delete();
  }
  // cout << "c1  " << count1 << "  c2 " << count2 << endl;
  // cout << "c1g  " << count1g << "  c2g " << count2g << endl;
  //canvas for azimut study
  //TCanvas *t=new TCanvas();
  //hPhi->DrawCopy();
  //t->SaveAs("SmallAcceptanceStudy.png");
  //hPhi->Write();
  //hPhi->Delete();

  outputFile.Write();
  vertice->Delete();
  hitsLayer[0]->Delete();
  hitsLayer[1]->Delete();
  inputFile.Close();
  outputFile.Close();
  timer.Stop();
  timer.Print();
}