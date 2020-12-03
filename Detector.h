#ifndef DETECTOR_h
#define DETECTOR_h

#include <string.h>
using namespace std;
#include "TClonesArray.h"
#include "TH1F.h"
#include "Vertex.h"
#include "Track.h"
#include "TTree.h"

const int nLayer = 3;                                     // number of layers (Beam Pipe included)
double radius[nLayer] = {3, 4, 7};                        //cm radius of each layer
double length = 27;                                       //cm detector length
std::string names[nLayer] = {"Beam", "Layer1", "Layer2"}; // names of each layer

class Detector : public TObject
{
public:
    Detector();
    virtual ~Detector();
    void SetRLayer(double *); //setting radius
    void SetNames(string *);
    int GetRLayer(int i) const { return (i >= 0 && i < nLayer) ? fR[i] : -1; }
    string GetLName(int i) const { return (i >= 0 && i < nLayer) ? fNames[i] : "error"; }
    void PrintDet() const; //print detector features
    double GetH() const { return fH; }
    TClonesArray *GetPtrhits(int i) const { return fPtrhits[i]; }
    void SetH(double);
    //Simulation: function to start the simulation
    void Simulation(bool, Vertex&, Track &, TH1F *, TH1F *);
    //void Simulation(bool, Vertex &, Track &, TH1F *); //in case of fixed or uniform multiplicity

private:
    double fR[nLayer]; // radius of each layer
    string fNames[nLayer];
    double fH;                 // detector length
    TClonesArray *fPtrhits[2]; //array of TClonesArray to save the hits on each layer

    ClassDef(Detector, 1)
};
#endif
