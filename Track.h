#ifndef TRACK_H
#define TRACK_H

#include "TObject.h"
#include "TH1F.h"
#include "Vertex.h"
#include "TClonesArray.h"

int nNoise = 10; //number of noise points
class Track : public TObject
{
public:
    Track(double = 0, double = 0);
    double GetPhi() const;
    double GetPhiP() const;
    void SetFixedPhi(double);
    void SetFixedTheta(double);
    double GetTheta() const;
    double GetThetaP() const;
    void SetRandomPhi();
    void SetRandomPhiP();
    void SetRandomTheta(TH1F *);
    void SetGaussThetaP();
    void PrintPhi() const;
    void PrintTheta() const;
    void PrintPhiP() const;
    void PrintThetaP() const;
    void Rotate(double *);                                                                   //function to change frame reference (from particle RF to lab RF)
    void Noise(double, double, double *);                                                    //function to create noise coordinates
    void SetDirection(Vertex &, TH1F *, TClonesArray **, const int, double *, double, bool); //function to determinate the propagation direction

private:
    double fTheta;
    double fPhi;
    double fThetaP;//Theta for multiple scattering
    double fPhiP; //Phi for multiple scattering

    ClassDef(Track, 1)
};
#endif