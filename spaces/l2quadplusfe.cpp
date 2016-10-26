
#include <fem.hpp>
#include "l2quadplusfe.hpp"

namespace dpg { 


  L2EnrichedQuad::L2EnrichedQuad (int k)
    : ScalarFiniteElement<2> ((k+1)*(k+1)+1, k+1),
    l2quad(k), _k(k) {}

  

  void L2EnrichedQuad::CalcShape (const IntegrationPoint & ip, 
				  BareSliceVector<> shape) const {
    double x = ip(0), y = ip(1);
    l2quad.CalcShape(ip, shape);
    shape( (_k+1)*(_k+1) ) = T_CalcLastShape(x,y);
  }

  void L2EnrichedQuad::CalcDShape (const IntegrationPoint & ip, 
				   SliceMatrix<> dshape) const {
    AutoDiff<2> x(ip(0),0),  y(ip(1),1);
    l2quad.CalcDShape(ip, dshape);
    int last = (_k+1)*(_k+1);
    AutoDiff<2> s = T_CalcLastShape(x,y);
    dshape(last, 0) = s.DValue(0);
    dshape(last, 1) = s.DValue(1);
  }

  template <class T>
  T L2EnrichedQuad::T_CalcLastShape(const T & x, const T & y) const  {
    
    if (_k%2==0)  {
      
      int n = (_k - 2) / 2;
      return (x*(1-x)-y*(1-y))*(2*x-1)*(2*y-1)*( pow(x*(1-x), n) +
						 pow(y*(1-y), n) );      
    }
    else  {
      
      int n = (_k - 1) / 2;
      return (x*(1-x)-y*(1-y)) * ( pow(x*(1-x), n) +
				   pow(y*(1-y), n) );
    }
  }
}
