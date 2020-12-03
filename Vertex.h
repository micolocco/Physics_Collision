#ifndef VERTEX_H
#define VERTEX_H
#include "Point.h"
#include "TH1F.h"

class Vertex : public Point
{
public:
  Vertex(double = 0, double = 0, double = 0, int = 0);
  void SetRandomX();
  void SetRandomY();
  void SetRandomZ();
  virtual void Stamp() const;
  void SetMult(TH1F *);//multiplicity distribution
  //void SetMult(); //fixed or uniform multiplicity
  int GetMult() const;
  double computeDisc(double *, double);
  double computeT(double *, double);
  bool check(double, double, double, double);

private:
  int fMult;
  ClassDef(Vertex, 1)
};

#endif
