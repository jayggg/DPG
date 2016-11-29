#include <solve.hpp>
#include <iostream>


/* Calculate error in numerical fluxes by polynomial extension,      */


using namespace ngsolve;
using namespace ngfem;


namespace dpg {

template<typename SCAL>
class NumProcFluxError : public NumProc  {
  /*

    Numproc FluxError
    -----------------
    
    Given computed numerical fluxes q on the element interfaces 
    and the exact flux Q, this numproc computes an approximation
    to the H^(-1/2) norm of q - Q on the element interfaces (sum
    of the squares of the H^(-1/2) norms of q - Q on the boundary
    of each element). The approximation to the fractional norm is
    computed using the best polynomial extension (of the difference 
    between q and an interpolant of Q) in the H(div) space.

    Required flags:
    
    -exactq=<Q>
       exact flux vector Q
    -discreteq=<q>
       a grid function q with no inner dofs
    -extensionspace=<RT>
       H(div) conforming fe space of NGSolve (the space into which the
       numerical fluxes are polynomially extended)
    -fespace=<fs>
       the space in which q lies 
    -errorsquareq=<qerrsqr>
       piecewise constant grid function whose value on each element 
       equals the square of the error norm
    -hdivproduct=<hdivipe>
       bilinear form representing H(div) inner product


   */

  shared_ptr<BilinearForm> hdivip;
   
  shared_ptr<GridFunction> Q;
  shared_ptr<GridFunction> q;
  shared_ptr<GridFunction> err;

  shared_ptr<FESpace> ext;
  shared_ptr<FESpace> fes;

public:

  NumProcFluxError ( shared_ptr<PDE>  apde, const Flags & flags) : NumProc(apde) {

    fes    = GetPDE()->GetFESpace(flags.GetStringFlag("fespace",NULL));
    ext    = GetPDE()->GetFESpace(flags.GetStringFlag("extensionspace",NULL));
    hdivip = GetPDE()->GetBilinearForm(flags.GetStringFlag("hdivproduct",NULL));
    q      = GetPDE()->GetGridFunction(flags.GetStringFlag("discreteq",NULL));
    Q      = GetPDE()->GetGridFunction(flags.GetStringFlag("exactq",NULL));
    err    = GetPDE()->GetGridFunction(flags.GetStringFlag("errorsquareq",NULL));
  }
  
  void Do(LocalHeap & lh) {    
    // We proceed in three steps:
    // 1.  Compute the difference between Q and q
    // 2.  Compute the H(div) Schur complement 
    // 3.  Apply Schur complement to the difference


    // grid function with (interpolated) exact flux, grad(u) 
    BaseVector& vecQ = Q->GetVector();    
    // numerical flux q
    BaseVector& vecq = q->GetVector(); 
    // p.w. constant gridfunction to store element-wise error
    BaseVector& errvec = err->GetVector();   
    errvec.FV<double>() = 0.0;
    double sqer =0.0;   // this will contain the total error square
    
    for(int k=0; k<ma->GetNE(); k++)  {
      
      ElementId ei (VOL, k);

      size_t elndof = ext->GetFE(ei,lh).GetNDof();
      
      Vector<SCAL> diff(elndof);           
      // dof nrs: global, global inner, local inner, local Schur
      Array<int>  Gn,     Ginn,         Linn,        Lsn;

      // compute the difference between Q and q
      ext->GetDofNrs(ei,Gn);        // Global# of all dofs on element k
      diff = SCAL(0.0);
      for(int j=0; j<elndof; j++)
	diff[j] = vecQ.FV<SCAL>()[Gn[j]] - vecq.FV<SCAL>()[Gn[j]];
      
      // H(div) Gram matrix (given in two parts in pde file)
      Matrix<double> elmat(elndof), elmat2(elndof);
      elmat = 0.0; elmat2 = 0.0;
      hdivip->GetIntegrator(0)->
	CalcElementMatrix(ext->GetFE(ei,lh),ma->GetTrafo(ei,lh),elmat,lh);
      hdivip->GetIntegrator(1)->
	CalcElementMatrix(ext->GetFE(ei,lh),ma->GetTrafo(ei,lh),elmat2,lh);
      elmat += elmat2;
    
      // compute the H(div) Schur complement 
      ext->GetInnerDofNrs(k,Ginn); // Global# of inner dofs on element k
      for(int j=0; j<elndof; j++)
	if (Ginn.Contains( Gn[j] ))
	  Linn.Append(j);          // Local#  of inner dofs on element k
	else
	  Lsn.Append(j);           // Local#  of Schur dofs on element k

      int ielndof = Linn.Size();
      Matrix<double> elmati(ielndof),elmatiinv(ielndof);
      elmati = elmat.Rows(Linn).Cols(Linn);
      CalcInverse(elmati,elmatiinv);
            
      // apply Schur complement to the difference
      int selndof = elndof - ielndof;
      Vector<SCAL> diffs(selndof);
      Matrix<double> S(selndof), Asi(selndof,ielndof);
      diffs = diff(Lsn);

      //      S  =  A_ss  -  A_si  * inv(A_ii) *  A_is
      Asi = elmat.Rows(Lsn).Cols(Linn);
      S   = elmat.Rows(Lsn).Cols(Lsn);
      S  -= Asi  * elmatiinv * Trans(Asi);
      //      error  = (S * diffs, diffs)
      errvec.FVDouble()[k] = fabs(InnerProduct(diffs,  S * diffs));
      sqer += errvec.FVDouble()[k];
    }
    
    cout<<"Discrete H^(-1/2) norm of error in q = "<<sqrt(sqer)<<endl;
    
    // write file (don't know what the last argument of AddVariable 
    // does, but 6 seems to be the value everywhere! It seems to intializes 
    // an object  of class IM (important message).
    GetPDE()->AddVariable (string("fluxerr.")+GetName()+".value", sqrt(sqer), 6);  

  }
   
  virtual string GetClassName () const {
      return "FluxError";
  }

    
};//class

static RegisterNumProc<NumProcFluxError<double>> npinitfluxerr("fluxerr");
static RegisterNumProc<NumProcFluxError<Complex>> npinitfluxerrc("fluxerrc");

};//namespace
