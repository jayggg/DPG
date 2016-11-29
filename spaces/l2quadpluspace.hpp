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

    virtual string GetClassName () const { return "L2EnrichedQuadFESpace"; }

    virtual void Update(LocalHeap & lh);
    virtual size_t GetNDof () const { return ndof; }

    virtual void GetDofNrs (ElementId ei,  Array<int> & dnums) const;
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
  };


}

#endif //  FILE_ENRICHQUADFESPACE_HPP
