#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Detector.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"

void ResultsStudy()
{
  Vertex *vertice = new Vertex();
  //opening input file
  TFile inputFile1("htree.root");
  TFile inputFile2("zVertez.root");
  TFile outputFile("Results.root", "recreate");
  //reading
  TTree *tree = (TTree *)inputFile1.Get("T");
  //tree->SetDirectory(0);
  TBranch *b1 = tree->GetBranch("VertMult");
  b1->SetAddress(&vertice);
  //definition of variable bins for TH3D
  const int bin = 10;
  double limitMult[bin + 1] = {2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 11.5, 20.5, 30.5, 40.5, 50.5};
  // double limitMult[bin + 4] = {2, 3,4,5, 6,7, 8,10,15, 20, 30, 40, 60, 80}; //in case of mult from uniform distribution
  double limitZ[bin + 1] = {-20, -15, -10, -5, -2.5, 0, 2.5, 5, 10, 15, 20};
  //double limitZ[bin +7 ] = {-24,-21,-18,-15, -12,-8,-5, -2.5,0, 2.5, 5,8,12, 15,18, 21,24}; //in case of z from uniform distribution
  const int binRes = 200;
  double limitRes[binRes + 1];
  int var = -700;
  for (int i = 0; i <= binRes; i++)
  {
    limitRes[i] = var;
    var += 7;
  }

  TH1D *effic[2]; //array of histogram: effic[0]=efficiency(multiplicity), effic[1]=efficiency(zGen)
  TH1D *MultGen = new TH1D("MultGen", "MultGen", bin, limitMult);
  TH1D *MultRec = new TH1D("MultRec", "MultRec", bin, limitMult);
  // TH1D *MultGen = new TH1D("MultGen", "MultGen", 13, limitMult); //in case of mult from uniform distribution
  // TH1D *MultRec = new TH1D("MultRec", "MultRec", 13, limitMult);
  TH1D *zetaGen = new TH1D("zetaGen", "zetaGen", bin, limitZ);
  //TH1D *zetaGen = new TH1D("zetaGen", "zetaGen", 16, limitZ); //in case of z from uniform distribution
  zetaGen->GetXaxis()->SetTitle("zGenerated [cm]");
  TH1D *zetaRec = new TH1D("zetaRec", "zetaRec", bin, limitZ);
  //TH1D *zetaRec = new TH1D("zetaRec", "zetaRec", 16, limitZ); //in case of z from uniform distribution

  TH3D *resol = new TH3D("Resolution", "Resolution", bin, limitMult, binRes, limitRes, bin, limitZ);
  //TH3D *resol = new TH3D("Resolution", "Resolution", bin, limitMult, binRes, limitRes, 16, limitZ); //in case of z from uniform distribution
  //TH3D *resol = new TH3D("Resolution", "Resolution", 13, limitMult, binRes, limitRes, 16, limitZ); //in case of mult from uniform distribution
  TH1D *resHist = new TH1D("Ris", "Resolution (All)", 100, -1500, 1500);
  resHist->GetXaxis()->SetTitle("zGenerated-zReconstructed [#mum]");
  resol->GetYaxis()->SetTitle("zGenerated-zReconstructed [#mum]");
  resol->GetXaxis()->SetTitle("Multeplicity");
  resol->GetZaxis()->SetTitle("zVertice [cm]");
  TF1 *g1 = new TF1("m1", "gaus", -250, 250);
  TF1 *g2 = new TF1("m2", "gaus", -700, 700);
  TNtuple *zVect = (TNtuple *)inputFile2.Get("zVect");
  //zVect->SetDirectory(0);
  double zRec[zVect->GetEntries()];
  for (int iEvent = 0; iEvent < zVect->GetEntries(); iEvent++)
  {
    zVect->GetEntry(iEvent);
    zRec[iEvent] = (zVect->GetLeaf("zVertice")->GetValue(0)) * 10000;
  }

  //Efficiency study and TH3D resol filling
  double res;
  for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++)
  {
    tree->GetEvent(iEvent);
    if (TMath::Abs(vertice->GetZ()) <= 3 * 5.3)
    {
      // if (TMath::Abs(vertice->GetZ()) <= 5.3){
      MultGen->Fill(vertice->GetMult());
      zetaGen->Fill(vertice->GetZ());
    }
    if (zRec[iEvent] != 500000)
    {
      res = 10000 * vertice->GetZ() - zRec[iEvent];
      resHist->Fill(res);
      // if (TMath::Abs(vertice->GetZ()) <= 5.3){
      if (TMath::Abs(vertice->GetZ()) <= 3 * 5.3)
      {
        resol->Fill(vertice->GetMult(), res, vertice->GetZ());
        MultRec->Fill(vertice->GetMult());
        zetaRec->Fill(vertice->GetZ());
      }
    }
  }

  effic[0] = (TH1D *)MultGen->Clone("efficiency VS Mult");
  effic[1] = (TH1D *)zetaGen->Clone("efficiency VS zGen");
  TCanvas canvas1("eff", "eff");
  canvas1.cd();
  for (int i = 0; i < 2; i++)
  {
    effic[i]->Sumw2();
    if (i == 0)
    {
      effic[i]->SetStats(0);
      effic[i]->Divide(MultRec, MultGen, 1, 1, "B");
    }
    else
      effic[i]->Divide(zetaRec, zetaGen, 1, 1, "B");
    effic[i]->UseCurrentStyle();
    effic[i]->SetMarkerSize(0.8);
    effic[i]->SetMarkerColor(2);
    effic[i]->SetMarkerStyle(20);
    effic[i]->SetLineColor(1);
    effic[i]->SetStats(0);
    if (i == 0)
    {
      effic[i]->SetTitle("efficiency VS multeplicity");
      effic[i]->SetName("efficiency VS multeplicity");
      canvas1.SetName("efficiency VS multeplicity");
    }
    else
    {
      effic[i]->SetTitle("efficiency VS zVertex");
      effic[i]->SetName("efficiency VS zVertex");
      canvas1.SetName("efficiency VS zVertex");
    }
    effic[i]->Draw("E1");
    effic[i]->Draw("same LHist");
    canvas1.Write();
  }

  //ProjectionY of TH3D (resol) to study resolution VS multiplicity, resolution VS zVert
  //resolution obtained as the std deviation of a gaussian fit
  TH1D *fProj;
  double resVectMult[bin], resVectMultError[bin], Mult[bin], MultError[bin]; //vector to store the resolution and the mult or zVert bin associated
  //double resVectMult[13], resVectMultError[13], Mult[13], MultError[13]; //in case of uniform distribution
  double resVectZ[bin], resVectZError[bin], zVert[bin], zVertError[bin];
  //double resVectZ[16], resVectZError[16], zVert[16], zVertError[16]; //in case of uniform distribution

  int nBins[] = {resol->GetNbinsX(), resol->GetNbinsZ()};

  for (int iProj = 0; iProj < 2; iProj++)
  {
    for (int iBin = 1; iBin <= nBins[iProj]; ++iBin)
    {
      if (iProj == 0)
      {
        fProj =
            resol->ProjectionY(Form("Multeplicity= [%.2f,%.2f] ", resol->GetXaxis()->GetBinLowEdge(iBin), resol->GetXaxis()->GetBinLowEdge(iBin) + resol->GetXaxis()->GetBinWidth(iBin)), iBin, iBin); // prjection of iBin
        fProj->SetTitle(Form("Multeplicity= [%.2f,%.2f] ", resol->GetXaxis()->GetBinLowEdge(iBin), resol->GetXaxis()->GetBinLowEdge(iBin) + resol->GetXaxis()->GetBinWidth(iBin)));
      }
      else
      {
        fProj =
            resol->ProjectionY(Form("zVert= [%.2f,%.2f] ", resol->GetZaxis()->GetBinLowEdge(iBin), resol->GetZaxis()->GetBinLowEdge(iBin) + resol->GetZaxis()->GetBinWidth(iBin)), 1, resol->GetNbinsX(), iBin, iBin);
        fProj->SetTitle(Form("zVert= [%.2f,%.2f] ", resol->GetZaxis()->GetBinLowEdge(iBin), resol->GetZaxis()->GetBinLowEdge(iBin) + resol->GetZaxis()->GetBinWidth(iBin)));
      }
      fProj->GetYaxis()->SetTitle("counts");

      g1->SetParameters(1, 0);
      g1->SetParameters(2, 150);
      g2->SetParameters(1, 0);
      //g2->SetParameters(2, 300);
      if ((iProj == 1 && TMath::Abs(resol->GetZaxis()->GetBinCenter(iBin)) > 15)||iProj==0)
        //if (iProj == 1 ||(resol->GetXaxis()->GetBinCenter(iBin) < 20 && iProj == 0)) //in case of uniform distribution (multiplicity, vertex)
        fProj->Fit(g2, "RLQ");
      else
        fProj->Fit(g1, "RLQ");
      gStyle->SetOptFit(1111);
      gStyle->SetOptStat(00000);

      if (iProj == 0)
      {
        // if (resol->GetXaxis()->GetBinCenter(iBin) > 20) //in case of uniform distribution (multiplicity, vertex)
        // {
        resVectMult[iBin - 1] = g2->GetParameter(2);
        resVectMultError[iBin - 1] = g2->GetParError(2);
        Mult[iBin - 1] = resol->GetXaxis()->GetBinCenter(iBin);
        MultError[iBin - 1] = resol->GetXaxis()->GetBinWidth(iBin) / 2;
        // }
        // else //in case of uniform distribution (multiplicity, vertex)
        // {
        //   resVectMult[iBin - 1] = g2->GetParameter(2);
        //   resVectMultError[iBin - 1] = g2->GetParError(2);
        //   Mult[iBin - 1] = resol->GetXaxis()->GetBinCenter(iBin);
        //   MultError[iBin - 1] = resol->GetXaxis()->GetBinWidth(iBin) / 2;
        // }
      }

      if (iProj == 1)
      {
        if (TMath::Abs(resol->GetZaxis()->GetBinCenter(iBin)) > 15)
        {
          resVectZ[iBin - 1] = g2->GetParameter(2);
          resVectZError[iBin - 1] = g2->GetParError(2);
          zVert[iBin - 1] = resol->GetZaxis()->GetBinCenter(iBin);
          zVertError[iBin - 1] = resol->GetZaxis()->GetBinWidth(iBin) / 2;
        }
        else
        {
          resVectZ[iBin - 1] = g1->GetParameter(2);
          resVectZError[iBin - 1] = g1->GetParError(2);
          zVert[iBin - 1] = resol->GetZaxis()->GetBinCenter(iBin);
          zVertError[iBin - 1] = resol->GetZaxis()->GetBinWidth(iBin) / 2;
        }
      }
    }
  }
  //Resolution VS Multeplicity, Resolution VS zVertice
  auto resMult = new TGraphErrors(bin, Mult, resVectMult, MultError, resVectMultError);
  //auto resMult = new TGraphErrors(13, Mult, resVectMult, MultError, resVectMultError); //in case of mult from uniform distribution
  resMult->GetXaxis()->SetTitle("Multeplicity");
  resMult->SetName("Resolution VS Multeplicity");
  resMult->SetTitle("Resolution VS Multeplicity");
  resMult->Draw("ALP");
  resMult->GetYaxis()->SetTitle("Resolution [#mum]");
  resMult->SetMarkerColor(2);
  resMult->SetMarkerStyle(20);
  resMult->Write();

  auto resZvert = new TGraphErrors(bin, zVert, resVectZ, zVertError, resVectZError);
  //auto resZvert = new TGraphErrors(16, zVert, resVectZ, zVertError, resVectZError); //in case of z from uniform distribution
  resZvert->GetXaxis()->SetTitle("zVert [cm]");
  resZvert->SetName("Resolution VS zVert");
  resZvert->SetTitle("Resolution VS zVert");
  resZvert->Draw("ALP");
  resZvert->GetYaxis()->SetTitle("Resolution [#mum]");
  resZvert->SetMarkerColor(2);
  resZvert->SetMarkerStyle(20);
  resZvert->Write();

  MultGen->Delete();
  MultRec->Delete();
  zetaGen->Delete();
  zetaRec->Delete();
  TCanvas canvas("Resolution(All)", "Resolution(All)");
  canvas.cd();
  resHist->Draw("pe");

  //effic[1]->Draw("L");
  resHist->SetMarkerColor(2);
  resHist->SetMarkerStyle(20);
  gStyle->SetOptStat(1111);
  canvas.SaveAs("ResAll");
  canvas.Write();
  g1->Delete();
  g2->Delete();
  vertice->Delete();
  outputFile.Write();
  resol->Delete();
  resHist->Delete();
  resZvert->Delete();
  resMult->Delete();
  inputFile1.Close();
  inputFile2.Close();
  outputFile.Close();
}