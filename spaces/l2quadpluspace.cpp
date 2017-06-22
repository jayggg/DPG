#include <comp.hpp>    
#include "l2quadplusfe.hpp"
#include "l2quadpluspace.hpp"

namespace dpg {
  
  L2EnrichedQuadFESpace::L2EnrichedQuadFESpace (shared_ptr<MeshAccess> ama,
				  const Flags & flags)
    : FESpace (ama, flags)   {
    
    _k = int(flags.GetNumFlag ("order", 2));
    evaluator[VOL] =
      make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] =
      make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
    evaluator[BND] =
      make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
    integrator[VOL] = GetIntegrators() .
      CreateBFI("mass", ma->GetDimension(), 
		make_shared<ConstantCoefficientFunction>(1));
  }

  
  void L2EnrichedQuadFESpace::Update(LocalHeap & lh)   {
   
    int n_cell = ma->GetNE();  
    int ii = 0;

    first_cell_dof.SetSize (n_cell+1);
    for (int i = 0; i < n_cell; i++, ii+=(_k+1)*(_k+1)+1)
      first_cell_dof[i] = ii;
    first_cell_dof[n_cell] = ii;
    ndof = ii;

    ctofdof.SetSize(ndof);
    ctofdof = LOCAL_DOF;
  }

  void L2EnrichedQuadFESpace::GetDofNrs (ElementId ei, Array<int> & dnums) const {

    int elnr = ei.Nr();
    dnums.SetSize(0);

    if ( ei.VB() == VOL ) {
      int first = first_cell_dof[elnr];
      int next  = first_cell_dof[elnr+1];
      for (int j = first; j < next; j++)  dnums.Append (j);
    }
  }
  
//  void L2EnrichedQuadFESpace::GetSDofNrs (int elnr, Array<int> & dnums) const 
//  { }


  //const FiniteElement & 
  //L2EnrichedQuadFESpace::GetFE (int elnr, LocalHeap & lh) const  {
  FiniteElement & 
  L2EnrichedQuadFESpace::GetFE (ElementId ei, Allocator & lh) const  {

    L2EnrichedQuad * quad = new (lh) L2EnrichedQuad(_k);
    //dd: Ngs_Element ngel = ma->GetElement (elnr);
    Ngs_Element ngel = ma->GetElement (ei);

    for (int i = 0; i < 4; i++)
      quad->SetVertexNumber (i, ngel.vertices[i]);
    
    return *quad;
  }

  /* 
  const FiniteElement &
  L2EnrichedQuadFESpace::GetSFE(int selnr, LocalHeap & lh) const {

    throw Exception("Called GetSFE! Not implemented yet!");
  }
  */

  
  static RegisterFESpace<L2EnrichedQuadFESpace> initifes ("l2quadplus");
}
