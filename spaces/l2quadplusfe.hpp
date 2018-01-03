#ifndef FILE_LENRICHEDQUADELEM_HPP 
#define FILE_LENRICHEDQUADELEM_HPP

/*
  This is an implementation of a finite element that equals NGsolve's
  Quadrilateral L2HighOrderFE plus one function of higher degree. 
 */


using namespace ngfem;

namespace dpg {

  ///  L2EnrichedQuad(k) = Q_{k,k} + one_degree_k+1_function
  
  class L2EnrichedQuad : public ScalarFiniteElement<2>   {
    
    int vnums[4];
    L2HighOrderFE<ET_QUAD> l2quad;  // Q_{k,k} 
    int _k;                         // degree of Q_{k,k}
    
  public:

    L2EnrichedQuad (int k);

    int Order() {return _k+1 ;}     // highest degree of shapes
    
    virtual ELEMENT_TYPE ElementType() const { return ET_QUAD; }
    
    void SetVertexNumber (int i, int v) {
      vnums[i] = v;
      l2quad.SetVertexNumber(i,v);
    }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceVector<> shape) const;  
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceMatrix<> dshape) const;
    
  private:

    template <class T> T T_CalcLastShape(const T & x, const T & y) const;
  };
  
}

namespace ngstd  {

  /// Integer powers of AutoDiff variables (used in L2EnrichedQuad shape)
  
  template<int D, typename SCAL>
  AutoDiff<D,SCAL> pow (const AutoDiff<D,SCAL> & f, int n)   {

    switch (n) {
      
    case 0: return AutoDiff<D,SCAL> (1.0);
      
    case 1: return f;
      
    default:

      AutoDiff<D,SCAL> rslt(0);
      rslt.Value() = pow(f.Value(),n);
      for (int i=0; i<D; i++)
	rslt.DValue(i) = n*pow(f.Value(),n-1)*f.DValue(i);

      return rslt;
      
    }
  }
}

#endif  // FILE_LENRICHEDQUADELEM_HPP

