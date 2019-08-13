// -*- C++ -*-
#if !defined(thrust_H)
#define thrust_H

#include <vector>
#include "TVector2.h"


namespace nTupleHelperTools {

  TVector2 calcThrust(const std::vector<TVector2>& inputPts, bool debug = false);
  
  void calcT(const std::vector<TVector2>& momenta, double& t, TVector2& taxis);

}
#endif // thrust_H
