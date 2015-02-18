#include <solve.hpp>


using namespace ngsolve;
using namespace ngfem;

namespace dpg  {

  template<typename SCAL>
  class NumProcEnorms : public NumProc   {

  /*

    Numproc Enorms
    ----------------

    DPG error representation functions are components of the compund
    DPG solution. Its Y-norms are used to provide an error esimator
    for marking elements for refinement in an adaptive scheme. This
    numproc computes element-wise Y-norms of the computed error
    representation functions.

    Required flags:
    
    -solution=<compoundsol>
        computed solution in a compound grid function, out of which 
	the error representation components are extracted using 
	indices provided in flag "yspaces".

    -estimator=<eestim>
        element-wise constant function, containing error estimator values
	for each element

    -fespace=<fs>
        the compound FE space

    -bilinearform=<dpg>
       the DPG bilinear form in mixed form on the compound FE space

    -yintegrators=[i1,i2,...iN]
       the Y-inner product is based on the integrators numbered
       i1,..,iN in the dpg bilinearform. 

       Example: If the pde file contains 

            define bilinearform dpg -fespace=fs
	    gradgrad    lam
	    laplace     one  -comp=3
	    mass        one  -comp=3
	    <many other integrators...>

       then -yintegrators=[2,3] indicates that Y-norm is the H1-norm.

    -yspaces=[j1,j2,...jM]
       If the compound FE space is  Y[1] x Y[2] x ..., then this 
       flag indicates that the test spaces consist of spaces  
       Y[j1] x Y[j2] x ... x Y[jM].  

       In the above example, -yspaces=[3] indicate the test 
       space is Y[3], the third component in the compound space.

   */

  protected:
    
    shared_ptr<FESpace>       fes;
    shared_ptr<GridFunction>  sol;
    shared_ptr<GridFunction>  est;
    shared_ptr<BilinearForm>  bfa;
    Array<int> Yintegrators;
    Array<int> Yspaceind;

  public:
    
    NumProcEnorms(shared_ptr<PDE>  apde, const Flags & flags) : NumProc(apde) {

      cout << "Constructor of NumProcEnorms" 
	   << "with flags:" << endl << flags;

      fes = GetPDE()->GetFESpace(flags.GetStringFlag("fespace",NULL));
      bfa = GetPDE()->GetBilinearForm (flags.GetStringFlag("bilinearform",NULL));
      sol = GetPDE()->GetGridFunction(flags.GetStringFlag("solution",NULL));
      est = GetPDE()->GetGridFunction(flags.GetStringFlag("estimator",NULL));

      Array<double> yintegrators (flags.GetNumListFlag("yintegrators"));
      for (int i=0; i<yintegrators.Size(); i++) {
	Yintegrators.Append(int(yintegrators[i]-1));
      }

      Array<double> yspaceind (flags.GetNumListFlag("yspaces"));
      for (int i=0; i<yspaceind.Size(); i++) {
	Yspaceind.Append( int( yspaceind[i]-1)  );
      }
      
    }

    virtual string GetClassName () const {

      return "Enorms";
    }

    
    virtual void Do (LocalHeap & lh)  {

      BaseVector& estvec = est->GetVector();   
      BaseVector& solvec = sol->GetVector();   
      estvec.FVDouble() = 0.0;

      for(int k=0; k<ma->GetNE(); k++)  {
	
	const FiniteElement & fel = fes->GetFE(k,lh);
	const CompoundFiniteElement & cfel = 
	  dynamic_cast<const CompoundFiniteElement&>(fel);


	int ndofel   = cfel.GetNDof();
	Matrix<SCAL> elmat(ndofel), xelmat(ndofel);
	elmat = SCAL(0.0);

	// Compute the Y-norm
	for (int ii=0; ii<bfa->NumIntegrators(); ii++) {
	  
	  if (Yintegrators.Contains(ii)) {
	    bfa->GetIntegrator(ii)->
	      CalcElementMatrix(cfel,ma->GetTrafo(k,0,lh),xelmat,lh);
	    elmat += xelmat;
	  }
	} 


	for (int ii=0; ii<Yspaceind.Size(); ii++) {
	  const FiniteElement & efel = cfel[Yspaceind[ii]];
	  IntRange ri = cfel.GetRange(Yspaceind[ii]);
	  int ni = ri.Size();

	  Matrix<SCAL> Aii(ni,ni);
	  Aii = elmat.Rows(ri).Cols(ri);

	  Array<int> Gdofnrs, Gedofnrs;
	  fes->GetDofNrs(k,Gdofnrs);	
	  Gedofnrs = Gdofnrs[ri];
	  Vector<SCAL> e(ni);
	  e = solvec.FV<SCAL>()(Gedofnrs);

	  estvec.FVDouble()[k] += fabs(InnerProduct( e, Aii * e));
	}
      }
      cout << "Error estimator total norm = "
	   << L2Norm(estvec.FVDouble()) << endl;
    }
  
  };

  static RegisterNumProc<NumProcEnorms<double> >  npdpgest ("enorms");
  static RegisterNumProc<NumProcEnorms<Complex> > npdpgestc("enormsc");


} // namespace

