#include <comp.hpp>
using namespace ngcomp;

class PeriodicHCurlSpace : public HCurlHighOrderFESpace
{
private:

  Array<int> vertmapx;
  Array<int> vertmapy;

  Array<int> dofmapx;
  Array<int> dofmapy;

  // In 2D: if segment i is copied to segment j, then idx = j.
  // In 3D: index is the number of the "identify periodic " statement 
  //        in geo file.
  
  int idx; 
  int idy;  

public:

  PeriodicHCurlSpace (shared_ptr<MeshAccess> ama, const Flags & flags);

  virtual ~PeriodicHCurlSpace ();

  virtual string GetClassName () const
    {
      return "PeriodicHCurlSpace";
    }

  virtual void Update (LocalHeap & lh);

  virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
  virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

  virtual FiniteElement & GetFE (int enr, LocalHeap & lh) const;
  virtual FiniteElement & GetSFE (int enr, LocalHeap & lh) const;
  virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const;
};

PeriodicHCurlSpace :: PeriodicHCurlSpace (shared_ptr<MeshAccess> ama, const Flags & flags)
  : HCurlHighOrderFESpace (ama, flags)
{
  idx = 1;
  idy = 2; 
}


PeriodicHCurlSpace :: ~PeriodicHCurlSpace ()
{
  ;
}


void PeriodicHCurlSpace :: Update (LocalHeap & lh)
{

  dofmapx.SetSize(0); // Sentinels: dofmaps are not yet set
  dofmapy.SetSize(0); 

  HCurlHighOrderFESpace::Update(lh);

  // set dof maps to identity 
  dofmapx.SetSize (GetNDof());
  for (int i = 0; i < dofmapx.Size(); i++)
    dofmapx[i] = i;

  dofmapy.SetSize (GetNDof());
  for (int i = 0; i < dofmapy.Size(); i++)
    dofmapy[i] = i;

  // vertex-pair to edge hashtable
  HashTable<INT<2>, int> vp2e(ma->GetNEdges());

  for (int enr = 0; enr < ma->GetNEdges(); enr++)
  {
    int v1, v2;
    ma->GetEdgePNums (enr, v1, v2);
    if (v1 > v2) Swap (v1, v2);
    vp2e[INT<2>(v1,v2)] = enr;
  }

  // vertex-triple-or-quartet to face hashtable
  // (triangular faces get a dummy vertex number equal to -1)
  HashTable<INT<4>, int> v2f(ma->GetNFaces());

  Array<int> pnums;
  for (int fnr = 0; fnr < ma->GetNFaces(); fnr++)
  {
    ma->GetFacePNums (fnr, pnums);
    INT<4> i4;
    if (pnums.Size() == 3) 
      i4 = {-1, pnums[0], pnums[1], pnums[2]};
    if (pnums.Size() == 4) 
      i4 = {pnums[0], pnums[1], pnums[2], pnums[3]};
    i4.Sort();
    v2f[i4] = fnr;
  }

  // idx

  // vertex slave -> master array
  vertmapx.SetSize(ma->GetNV());
  for (int i = 0; i < vertmapx.Size(); i++)
    vertmapx[i] = i;

  Array<INT<2> > periodic_verts; 
  ma->GetPeriodicVertices(idx, periodic_verts);

  for (auto pair : periodic_verts)
    vertmapx[pair[1]] = pair[0];

  // find periodic edges (using vertex-pair to edge hashtable)
  for (int enr = 0; enr < ma->GetNEdges(); enr++)
  {
    int v1, v2;
    ma->GetEdgePNums (enr, v1, v2);
    // number of master-vertices
    int mv1 = vertmapx[v1];   // 
    int mv2 = vertmapx[v2];
    if (v1 != mv1 && v2 != mv2) // edge shall be mapped
    {      
      if (mv1 > mv2) Swap (mv1, mv2);
      int menr = vp2e[INT<2>(mv1,mv2)];  // the master edge-nr

      dofmapx[enr] = menr;

      IntRange edofs = GetEdgeDofs (enr);   // dofs on slave edge
      IntRange medofs = GetEdgeDofs (menr); // dofs on master edge
      for (int i = 0; i < edofs.Size(); i++)
	dofmapx[edofs[i]] = medofs[i];
    }
  }

  // find periodic faces (using vertex-triple to face hashtable)
  for (int fnr = 0; fnr < ma->GetNFaces(); fnr++)
  {
    ma->GetFacePNums (fnr, pnums);
    INT<4> i4;

    if (pnums.Size() == 3) 
    {
      i4 = {-1, vertmapx[pnums[0]], vertmapx[pnums[1]], vertmapx[pnums[2]]};
      if (i4[1] != pnums[0] && i4[2] != pnums[1] && i4[3] != pnums[2])
      {
	i4.Sort();
	int mfnr = v2f[i4];	
	IntRange fdofs = GetFaceDofs (fnr);
	IntRange mfdofs = GetFaceDofs (mfnr);
	for (int i = 0; i < fdofs.Size(); i++)
	  dofmapx[fdofs[i]] = mfdofs[i];
      }
    }

    if (pnums.Size() == 4) 
    {
      i4 = {vertmapx[pnums[0]], vertmapx[pnums[1]], vertmapx[pnums[2]], vertmapx[pnums[3]]};  
      if (i4[0] != pnums[0] && i4[1] != pnums[1] && i4[2] != pnums[2] && i4[3] != pnums[3])
      {
	i4.Sort();
	int mfnr = v2f[i4];

	IntRange fdofs = GetFaceDofs (fnr);
	IntRange mfdofs = GetFaceDofs (mfnr);
	for (int i = 0; i < fdofs.Size(); i++)
	  dofmapx[fdofs[i]] = mfdofs[i];
      }
    }
  }

  
  // idy
  vertmapy.SetSize(ma->GetNV());
  for (int i = 0; i < vertmapy.Size(); i++)
    vertmapy[i] = i;

  ma->GetPeriodicVertices(idy, periodic_verts);

  for (auto pair : periodic_verts)
    vertmapy[pair[1]] = pair[0];

  // find periodic edges
  for (int enr = 0; enr < ma->GetNEdges(); enr++)
  {
    int v1, v2;
    ma->GetEdgePNums (enr, v1, v2);
    // number of master-vertices
    int mv1 = vertmapy[v1];   // 
    int mv2 = vertmapy[v2];
    if (v1 != mv1 && v2 != mv2) // edge shall be mapped
    {
      if (mv1 > mv2) Swap (mv1, mv2);
      int menr = vp2e[INT<2>(mv1,mv2)];  // the master edge-nr

      dofmapy[enr] = menr;

      IntRange edofs = GetEdgeDofs (enr);   // dofs on slave edge
      IntRange medofs = GetEdgeDofs (menr); // dofs on master edge
      for (int i = 0; i < edofs.Size(); i++)
	dofmapy[edofs[i]] = medofs[i];
    }
  }

  // find periodic faces
  for (int fnr = 0; fnr < ma->GetNFaces(); fnr++)
  {
    ma->GetFacePNums (fnr, pnums);
    INT<4> i4;

    if (pnums.Size() == 3) 
    {
      i4 = {-1, vertmapy[pnums[0]], vertmapy[pnums[1]], vertmapy[pnums[2]]};
      if (i4[1] != pnums[0] && i4[2] != pnums[1] && i4[3] != pnums[2])
      {
	i4.Sort();
	int mfnr = v2f[i4];	
	IntRange fdofs = GetFaceDofs (fnr);
	IntRange mfdofs = GetFaceDofs (mfnr);
	for (int i = 0; i < fdofs.Size(); i++)
	  dofmapy[fdofs[i]] = mfdofs[i];
      }
    }

    if (pnums.Size() == 4) 
    {
      i4 = {vertmapy[pnums[0]], vertmapy[pnums[1]], vertmapy[pnums[2]], vertmapy[pnums[3]]};  
      if (i4[0] != pnums[0] && i4[1] != pnums[1] && i4[2] != pnums[2] && i4[3] != pnums[3])
      {
	i4.Sort();
	int mfnr = v2f[i4];
	IntRange fdofs = GetFaceDofs (fnr);
	IntRange mfdofs = GetFaceDofs (mfnr);
	for (int i = 0; i < fdofs.Size(); i++)
	  dofmapy[fdofs[i]] = mfdofs[i];
      }
    }
  }

  for (int i = 0; i < dofmapx.Size(); i++)
    if (dofmapx[i] != i)
      ctofdof[i] = UNUSED_DOF;

  for (int i = 0; i < dofmapy.Size(); i++)
    if (dofmapy[i] != i)
      ctofdof[i] = UNUSED_DOF;
}



void PeriodicHCurlSpace::GetDofNrs (int elnr, Array<int> & dnums) const
{
  HCurlHighOrderFESpace::GetDofNrs(elnr,dnums);

  if (dofmapx.Size())
    for(int i=0; i<dnums.Size(); i++)
      dnums[i] = dofmapy[dofmapx[dnums[i]]];
}

void PeriodicHCurlSpace::GetSDofNrs (int selnr, Array<int> & dnums) const
{
  HCurlHighOrderFESpace::GetSDofNrs(selnr,dnums);

  if (dofmapx.Size()) //
    for(int i=0; i<dnums.Size(); i++)
      dnums[i] = dofmapy[dofmapx[dnums[i]]];
}


FiniteElement & PeriodicHCurlSpace::GetFE (int enr, LocalHeap & lh) const
{
  GetFE (ElementId(VOL, enr), lh);
}

FiniteElement & PeriodicHCurlSpace::GetSFE (int enr, LocalHeap & lh) const
{
  GetFE (ElementId(BND, enr), lh);
}

FiniteElement & PeriodicHCurlSpace::GetFE (ElementId ei, Allocator & alloc) const
{
  Ngs_Element ngel = ma->GetElement (ei);
  switch (ngel.GetType())
  {
  case ET_TRIG:
  {
    auto fe = new(alloc) HCurlHighOrderFE<ET_TRIG> (order);
    fe->SetVertexNumbers (vertmapy[vertmapx[ngel.Vertices()]]);
    if(ma->GetDimension()==2){
      fe->SetOrderCell (order_inner[ei.Nr()]);
      fe->ComputeNDof();
    }
    return *fe;
  }
  case ET_TET:
  {
    auto fe = new(alloc) HCurlHighOrderFE<ET_TET> (order);
    fe->SetVertexNumbers (vertmapy[vertmapx[ngel.Vertices()]]);
    fe->SetOrderCell (order_inner[ei.Nr()]);
    fe->ComputeNDof();
    return *fe;
  }
  case ET_PRISM:
  {
    auto fe = new(alloc) HCurlHighOrderFE<ET_PRISM> (order);
    fe->SetVertexNumbers (vertmapy[vertmapx[ngel.Vertices()]]);
    fe->SetOrderCell (order_inner[ei.Nr()]);
    fe->ComputeNDof();
    return *fe;
  }
  case ET_QUAD:
  {
    auto fe = new(alloc) HCurlHighOrderFE<ET_QUAD> (order);
    fe->SetVertexNumbers (vertmapy[vertmapx[ngel.Vertices()]]);
    if(ma->GetDimension()==2){
      fe->SetOrderCell (order_inner[ei.Nr()]);
      fe->ComputeNDof();
    }
    return *fe;
  }
  }
  throw Exception("GetFE: undefined element");
}
  


static RegisterFESpace<PeriodicHCurlSpace> init ("hcurlho_periodic");



