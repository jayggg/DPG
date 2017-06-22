#ifndef FILE_ENRICHQUADFESPACE_HPP
#define FILE_ENRICHQUADFESPACE_HPP

/* 
   Implementation of a DG space on quadrilateral meshes, where the
   function space on each element equals 
          Q_{k,k} +   an_extra_degree_k+1_function.
   This is the smallest known test space for the DPG method on rectangles,
   which provably provides optimal H1 convergence rates. 
 */

using namespace ngcomp;

namespace dpg {

  class L2EnrichedQuadFESpace : public FESpace  {

    int _k;
    int ndof;    
    Array<int> first_cell_dof;
    
  public:

    
    L2EnrichedQuadFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    virtual ~L2EnrichedQuadFESpace () {;}

    virtual string GetClassName () const override 
    { 
        return "L2EnrichedQuadFESpace"; 
    }

    virtual void Update(LocalHeap & lh) override;
    virtual size_t GetNDof () const override { return ndof; }

    virtual void GetDofNrs (ElementId ei, Array<int> & dnums) const override;

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
  };

}

#endif //  FILE_ENRICHQUADFESPACE_HPP
