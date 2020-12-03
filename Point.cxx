
#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Point.h"



ClassImp(Point)

//________________________________________________________________________
Point::Point():TObject(),
 fX(0.),
 fY(0.),
 fZ(0.){
   // default constructor
 }


//___________________________________________________________________________
Point::Point(Double_t X, Double_t Y, Double_t Z):TObject(),
 fX(X),
 fY(Y),
 fZ(Z){
	//standard constructor 
}	     
//___________________________________________________________________________
void Point::SetX(double xValue)
{
 fX = xValue;
}
void Point::SetY(double yValue)
{
  fY = yValue;
}
void Point::SetZ(double zValue)
{
  fZ = zValue;
}
//____________________________________________________________________________
void Point::Stamp()const{
  cout<<"Point coordinate [x,y,z]= ["<<fX<<","<<fY<<","<<fZ<<"]"<<endl;
}
//___________________________________________________________________________
Point::~Point()	 {
  // destructor
}

