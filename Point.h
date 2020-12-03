#ifndef POINT_H
#define POINT_H

#include "TObject.h"

class Point : public TObject
{

public:

Point();
Point(double X, double Y, double Z);
virtual void Stamp()const; //dynamic binding: object (not pointer) determines function called
void SetX(double);
void SetY(double);
void SetZ(double);

virtual ~Point();

 double GetX() const {return fX;} 
 double GetY() const {return fY;}
 double GetZ() const {return fZ;}

private:


double fX;
double fY;
double fZ;

ClassDef(Point,1)
};


#endif 


