#include "Observables.h"

// Get Primitive Energy Estimator
// See eq. 12.6.40, Ref. 2
double Paths::getPE()
{
  rowvec dR(nD);
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      dR = Displacement( R(iPart,iBead) , R(iPart,bL(iBead+1)) );
      total += 0.5 * mwp2 * dot(dR,dR) - oneOvernBead * getV(iPart, iBead);
    }
  }
  
  return nDnPartnBeadOver2Beta - total;
}

// Get Virial Energy Estimator
// See eq. 12.6.40, Ref. 2
double Paths::getVE()
{
  rowvec rc(nD), dR(nD);
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    rc.zeros();
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) rc += R(iPart,iBead);
    rc = rc * oneOvernBead;
    PutInBox(rc);
    
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      dR = Displacement( R(iPart,iBead) , rc );
      total += 0.5 * dot(dR,dR) * getdV(iPart, iBead) + getV(iPart, iBead);
    }
  }
  
  return nDnPartOver2Beta + oneOvernBead*total;
}

// Position Estimator
double Paths::getR()
{
  rowvec Rsum(nD);
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {  
    Rsum.zeros();
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1)
      Rsum += R(iPart,iBead);

    total += Rsum(0);
  }
  
  return oneOvernPartnBead * total;
}

// Position Squared Estimator
double Paths::getR2()
{
  double total = 0.0;

  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iBead = 0; iBead < nBead; iBead += 1) {
      total += dot( R(iPart,iBead) , R(iPart,iBead) );
    }
  }
  
  return oneOvernPartnBead * total;
}
