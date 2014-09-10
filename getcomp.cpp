#include <solve.hpp>

using namespace ngsolve;
using namespace ngfem;

namespace ngfem  {


class NumProcGetComponent : public NumProc  {
  /*
    Numproc GetComponent, a trivial NumProc!
    --------------------
    
    Returns a component of a compound FE function
    Required flags:
      -comp=<n>
          returns component #n of the compound FE function
      -re 
          return only real part of component #n 
      -im
          return only imag part of component #n 
      -compoundgf=<fname>
          name of the gridfunction in compound FE space
      -componentgf=<fname>
          name of the gridfunction for storing comp #n.
      
   */
  
  GridFunction * gf1;
  GridFunction * gf2;
  int ind;
  bool re;
  bool im;

public:
  
  NumProcGetComponent ( PDE & apde, const Flags & flags) : NumProc(apde) {
    
    gf1 = pde.GetGridFunction(flags.GetStringFlag("compoundgf",NULL));
    ind = flags.GetNumFlag("comp",1) - 1 ;
    gf2 = pde.GetGridFunction(flags.GetStringFlag("componentgf",NULL));
    re = flags.GetDefineFlag("re");
    im = flags.GetDefineFlag("im");
  }

  virtual void Do(LocalHeap & lh)  {

    
    cout << "GetComponent " << ind+1 << ", of type " 
	 << gf2 -> GetFESpace(). GetClassName() << flush ;
    int n = gf2->GetVector().FV<double>(). Size();

    if(re) {
      
      cout<<", REAL part"<<endl;
      
      for (int i=0; i< n; i++ ) 
	gf2->GetVector().FV<double>()[i]  = 
	  (gf1->GetComponent(ind)->GetVector().FV<Complex>()[i] ). real();
    }
    else if(im) {

      cout<<", IMAG part"<<endl;
      for (int i=0; i< n; i++ ) 
	gf2->GetVector().FV<double>()[i]  = 
	  (gf1->GetComponent(ind)->GetVector().FV<Complex>()[i] ). imag();
    }
    else {
      cout << endl ;
      gf2->GetVector()=gf1->GetComponent(ind)->GetVector();    
    }
  }


};

  static RegisterNumProc<NumProcGetComponent> npinitgetcomp("getcomp");

};
