#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Vertex.h"
#include "TRandom3.h"
#include "TH1F.h"

ClassImp(Vertex)

Vertex::Vertex(double xValue, double yValue, double zValue, int multValue) : Point(xValue, yValue, zValue)
{
  fMult = multValue;
}

//funzioni Set per direzione randomica
void Vertex::SetRandomX()
{
  double a = gRandom->Gaus(0., 0.01);
  Point::SetX(a); //cm
}
void Vertex::SetRandomY()
{
  double b = gRandom->Gaus(0., 0.01);
  Point::SetY(b); //cm
}
void Vertex::SetRandomZ()
{
double c = gRandom->Gaus(0., 5.3);//cm
//double c = gRandom->Uniform(-25, 25);
  Point::SetZ(c); //cm
}

void Vertex::Stamp() const
{
  cout << "---vertex---" << endl;
  Point::Stamp();
  cout << "multiplicicty " << fMult << endl;
}

 void Vertex::SetMult(TH1F *distMult){
  fMult = (int)distMult->GetRandom();
}

// void Vertex::SetMult()
// {
//  fMult = 10;  
//  fMult= gRandom->Integer(78)+2;
// }

int Vertex::GetMult() const
{
  return fMult;
}

double Vertex::computeDisc(double *c, double R)
{
  double disc = (Point::GetX() * c[0] + Point::GetY() * c[1]) * (Point::GetX() * c[0] + Point::GetY() * c[1]) - (c[0] * c[0] + c[1] * c[1]) * (Point::GetX() * Point::GetX() + Point::GetY() * Point::GetY() - R * R);
  return disc;
}

double Vertex::computeT(double *c, double disc)
{
  double t = (-(Point::GetX() * c[0] + Point::GetY() * c[1]) + TMath::Sqrt(disc)) / (c[0] * c[0] + c[1] * c[1]);
  return t;
}
bool Vertex::check(double disc, double t, double c3, double H)
{
  if (disc > 0 && t > 0 && (Point::GetZ() + c3 * t) >= -H / 2 && (Point::GetZ() + c3 * t) <= H / 2)
    return true;
  else
    return false;
}
/*else
      {
        //cout << "condizione su z0 non valida" << endl;
        return false;
      }
    }
    else
    {
      //cout << "t è negativo " << endl;
      return false;
    }
  }
  else
  {
    //cout << "ERRORE: il discriminante è negativo!" << endl;
    return false;
  }
}*/