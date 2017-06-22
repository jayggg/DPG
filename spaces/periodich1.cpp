/*
  H(grad) high order FE Space with periodicity in two directions.

  (See H(curl) periodic space for more comments on periodic dof maps!)

*/

#include <comp.hpp>
using namespace ngcomp;

class PeriodicH1Space : public H1HighOrderFESpace  {

private:

  Array<int> dofmapx;    
  Array<int> dofmapy; 

  Array<int> xid, yid;
  Array<double> xends, yends;

  void SetPeriodicIds();
  
  template <ELEMENT_TYPE ET> FiniteElement &
  T_GetFE (int elnr, Allocator & lh) const;

public:

  PeriodicH1Space (shared_ptr<MeshAccess> ama, const Flags & flags);
  virtual ~PeriodicH1Space () {;} 
  virtual string GetClassName () const { return "PeriodicH1Space"; }

  virtual void Update (LocalHeap & lh);
 
  virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
  virtual void GetSDofNrs (int elnr, Array<int> & dnums) const;

  virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const;
};


PeriodicH1Space :: PeriodicH1Space (shared_ptr<MeshAccess> ama, const Flags & flags)
  : H1HighOrderFESpace (ama, flags) {


  if ( flags.NumListFlagDefined("xends") &&
       flags.NumListFlagDefined("yends") )  {
    xends = flags.GetNumListFlag("xends");
    yends = flags.GetNumListFlag("yends");
    QuickSort(xends);
    QuickSort(yends);
  }
  else {
    cerr << "**** Provide flags xends and yends to use xy periodic spaces!"
	 << endl;
    exit(1);
  }
    
  cout << " Initializing H1 space " << endl;
  if (xends.Size()>0) 
    cout << "   with periodic identification of x="
	 << xends[0] << " and x="<< xends[1] << endl;
  if (yends.Size()>0) 
    cout << "   with periodic identification of y="
	 << yends[0] << " and y="<< yends[1] << endl;
  
  SetPeriodicIds();

  cout << "   Ids of x-periodic identifications: " << endl << xid << endl;
  cout << "   Ids of y-periodic identifications: " << endl << yid << endl;
}


void PeriodicH1Space::SetPeriodicIds () {

  
  Array<INT<2> > pairs;
  ma->GetPeriodicVertices(0, pairs);
  int totpairs = pairs.Size();
  int countpairs = 0 ;
  int id = 1;
  Vec<3> pt0, pt1;
    
  while (countpairs < totpairs) {

    ma->GetPeriodicVertices(id, pairs);
    ma->GetPoint(pairs[0][0], pt0);
    ma->GetPoint(pairs[0][1], pt1);


    bool xpair = false;
    if (pt0[0] < pt1[0]) {
      if ( (abs(xends[0]-pt0[0])<1.e-15) && (abs(xends[1]-pt1[0])<1.e-15) )
	xpair = true;
    }
    else if ( (abs(xends[1]-pt0[0])<1.e-15) && (abs(xends[0]-pt1[0])<1.e-15) )
      xpair = true;

    bool ypair = false;
    if (pt0[1] < pt1[1]) {
      if ( (abs(yends[0]-pt0[1])<1.e-15) && (abs(yends[1]-pt1[1])<1.e-15) )
	ypair = true;
    }
    else if ( (abs(yends[1]-pt0[1])<1.e-15) && (abs(yends[0]-pt1[1])<1.e-15) )
      ypair = true;

    if (xpair) xid.Append(id);
    if (ypair) yid.Append(id);

    // if (xpair || ypair) {
    //   cout << "id = " << id << endl;
    //   // << " pairs:" << endl
    //   // << pairs << endl;
    //   cout << "first pair:" << pt0 <<"," << pt1 << endl;

    // }
    
    countpairs += pairs.Size();
    id +=1;    
  }

}
  
void PeriodicH1Space::Update (LocalHeap & lh)  {

  dofmapx.SetSize(0); // Sentinels: dofmaps are not yet set
  dofmapy.SetSize(0); 

  H1HighOrderFESpace::Update (lh);
    
  // set dof maps to identity
  dofmapx.SetSize (GetNDof());
  for (int i = 0; i < dofmapx.Size(); i++)
    dofmapx[i] = i;
    
  dofmapy.SetSize (GetNDof());
  for (int i = 0; i < dofmapy.Size(); i++)
    dofmapy[i] = i;

  // build vertex-pair to edge hashtable:
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

  
  // idx ////////////////////////////////////////////////////////

  Array<INT<2> > pcbPairs;

  for (int idx : xid) {
    
    ma->GetPeriodicVertices(idx, pcbPairs);
  
    // first dofs are vertex dofs    
    for (auto pair : pcbPairs) {

      if (pair[1] < pair[0]) Swap(pair[1], pair[0]);
    
      dofmapx[pair[1]] = pair[0];
    }

    // find periodic edges (using vertex-pair to edge hashtable)
    for (int enr = 0; enr < ma->GetNEdges(); enr++)
      {
	int v1, v2;
	ma->GetEdgePNums (enr, v1, v2);
	// number of master-vertices
	// use that dofmap[0:nv] is exactly the vertex-map
	int mv1 = dofmapx[v1];   // 
	int mv2 = dofmapx[v2];
	if (v1 != mv1 && v2 != mv2) // edge shall be mapped
	  {
	    if (mv1 > mv2) Swap (mv1, mv2);
	    int menr = vp2e[INT<2>(mv1,mv2)];  // the master edge-nr

	    IntRange edofs = GetEdgeDofs (enr);   // dofs on slave edge
	    IntRange medofs = GetEdgeDofs (menr); // dofs on master edge

	    if ( edofs.First() < medofs.First() ) Swap(edofs, medofs);

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
	    i4 = {-1, dofmapx[pnums[0]], dofmapx[pnums[1]], dofmapx[pnums[2]]};
	    if (i4[1] != pnums[0] && i4[2] != pnums[1] && i4[3] != pnums[2])
	      {
		i4.Sort();
		int mfnr = v2f[i4];

		IntRange fdofs = GetFaceDofs (fnr);
		IntRange mfdofs = GetFaceDofs (mfnr);

		if ( fdofs.First() < mfdofs.First() ) Swap(fdofs, mfdofs);

		for (int i = 0; i < fdofs.Size(); i++)
		  dofmapx[fdofs[i]] = mfdofs[i];
	      }
	  }

	if (pnums.Size() == 4) 
	  {
	    i4 = {dofmapx[pnums[0]], dofmapx[pnums[1]], dofmapx[pnums[2]], dofmapx[pnums[3]]};  
	    if (i4[0] != pnums[0] && i4[1] != pnums[1] && i4[2] != pnums[2] && i4[3] != pnums[3])
	      {
		i4.Sort();
		int mfnr = v2f[i4];

		IntRange fdofs = GetFaceDofs (fnr);
		IntRange mfdofs = GetFaceDofs (mfnr);

		if ( fdofs.First() < mfdofs.First() ) Swap(fdofs, mfdofs);

		for (int i = 0; i < fdofs.Size(); i++)
		  dofmapx[fdofs[i]] = mfdofs[i];
	      }
	  }
      }
  }
    
  // idy  /////////////////////////////////////////////////////////

  for (int idy : yid) {
    ma->GetPeriodicVertices(idy, pcbPairs);

    // first dofs are vertex dofs    
    for (auto pair : pcbPairs) {

      if (pair[1] < pair[0]) Swap(pair[1], pair[0]);

      dofmapy[pair[1]] = pair[0];
    }

    // find periodic edges (using vertex-pair to edge hashtable)
    for (int enr = 0; enr < ma->GetNEdges(); enr++)
      {
	int v1, v2;
	ma->GetEdgePNums (enr, v1, v2);

	// number of master-vertices
	// use that dofmap[0:nv] is exactly the vertex-map
	int mv1 = dofmapy[v1];   // 
	int mv2 = dofmapy[v2];
	if (v1 != mv1 && v2 != mv2) // edge shall be mapped
	  {
	    if (mv1 > mv2) Swap (mv1, mv2);
	    int menr = vp2e[INT<2>(mv1,mv2)];  // the master edge-nr

	    IntRange edofs = GetEdgeDofs (enr);   // dofs on slave edge
	    IntRange medofs = GetEdgeDofs (menr); // dofs on master edge

	    if ( edofs.First() < medofs.First() ) Swap(edofs, medofs);
	    
	    for (int i = 0; i < edofs.Size(); i++)
	      dofmapy[edofs[i]] = medofs[i];
	  }
      }


    // find periodic faces (using vertex-triple to face hashtable)
    for (int fnr = 0; fnr < ma->GetNFaces(); fnr++)
      {
	ma->GetFacePNums (fnr, pnums);
	INT<4> i4;
	
	if (pnums.Size() == 3) 
	  {
	    i4 = {-1, dofmapy[pnums[0]],
		  dofmapy[pnums[1]], dofmapy[pnums[2]]};
	    
	    if (i4[1] != pnums[0] && i4[2] != pnums[1] && i4[3] != pnums[2])
	      {
		i4.Sort();
		int mfnr = v2f[i4];

		IntRange fdofs = GetFaceDofs (fnr);
		IntRange mfdofs = GetFaceDofs (mfnr);

		if ( fdofs.First() < mfdofs.First() ) Swap(fdofs, mfdofs);

		for (int i = 0; i < fdofs.Size(); i++)
		  dofmapy[fdofs[i]] = mfdofs[i];
	      }
	  }

	if (pnums.Size() == 4) 
	  {
	    i4 = {dofmapy[pnums[0]], dofmapy[pnums[1]],
		  dofmapy[pnums[2]], dofmapy[pnums[3]]};  
	    if (i4[0] != pnums[0] && i4[1] != pnums[1] &&
		i4[2] != pnums[2] && i4[3] != pnums[3])
	      {
		i4.Sort();
		int mfnr = v2f[i4];

		IntRange fdofs = GetFaceDofs (fnr);
		IntRange mfdofs = GetFaceDofs (mfnr);

		if ( fdofs.First() < mfdofs.First() ) Swap(fdofs, mfdofs);

		for (int i = 0; i < fdofs.Size(); i++)
		  dofmapy[fdofs[i]] = mfdofs[i];
	      }
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


void PeriodicH1Space :: GetDofNrs (int elnr, Array<int> & dnums) const
{
  H1HighOrderFESpace::GetDofNrs (elnr, dnums);

  if (dofmapx.Size())
    for(int i=0; i<dnums.Size(); i++)
      dnums[i] = dofmapy[dofmapx[dnums[i]]];
}

void PeriodicH1Space :: GetSDofNrs (int elnr, Array<int> & dnums) const
{
  H1HighOrderFESpace::GetDofNrs (ngfem::ElementId(BND,elnr), dnums);

  if (dofmapx.Size()) //
    for(int i=0; i<dnums.Size(); i++)
      dnums[i] = dofmapy[dofmapx[dnums[i]]];
}


FiniteElement &
PeriodicH1Space::GetFE (ElementId ei, Allocator & alloc) const {

  int elnr = ei.Nr();
  Ngs_Element ngel = ma->GetElement(ei);
  ELEMENT_TYPE eltype = ngel.GetType();

  switch (eltype) {

    case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, alloc);
                
    case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, alloc);
    case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, alloc);
                
    case ET_TET:     return T_GetFE<ET_TET> (elnr, alloc);
    case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, alloc);
    case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, alloc);
    case ET_HEX:     return T_GetFE<ET_HEX> (elnr, alloc);
                
    default:
      throw Exception ("illegal element in H1HoFeSpace::GetFE");
  }
}


template <ELEMENT_TYPE ET> FiniteElement &
PeriodicH1Space::T_GetFE (int elnr, Allocator & lh) const  {
  
  Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (elnr);
  auto hofe =  new (lh) H1HighOrderFE<ET> ();

  hofe->SetVertexNumbers( dofmapy[dofmapx[ngel.Vertices()]] );  
  
  switch (int(ET_trait<ET>::DIM))  {
    
    case 1:
      {
	hofe -> SetOrderEdge (0, order_inner[elnr][0]);
	break;
      }

    case 2:
      {
	hofe -> SetOrderEdge (order_edge[ngel.Edges()] );
	// hofe -> SetOrderFace (0, order_inner[elnr]);
	hofe -> SetOrderFace (0, order_face[ma->GetSElFace(elnr)]);
	break;
      }

    case 3: default:  
      {
	hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
	hofe -> SetOrderFace (order_face[ngel.Faces()]);
	hofe -> SetOrderCell (order_inner[elnr]);
	break;
      }
  }

  hofe -> ComputeNDof();
  return *hofe;
}



// FiniteElement & PeriodicH1Space :: GetFE (ElementId ei, Allocator & alloc) const
// {
//   Ngs_Element ngel = ma->GetElement (ei);
//   switch (ngel.GetType())
//   {
//   case ET_TRIG:
//   {
//     auto fe = new(alloc) H1HighOrderFE<ET_TRIG> (order);
//     fe->SetVertexNumbers (dofmapy[dofmapx[ngel.Vertices()]]);
//     if(ma->GetDimension()==2){
//       fe->SetOrderCell (order_inner[ei.Nr()]);
//       fe->ComputeNDof();
//     }
//     return *fe;
//   }
//   case ET_TET:
//   {
//     auto fe = new(alloc) H1HighOrderFE<ET_TET> (order);
//     fe->SetVertexNumbers (dofmapy[dofmapx[ngel.Vertices()]]);
//     fe->SetOrderCell (order_inner[ei.Nr()]);
//     fe->ComputeNDof();
//     return *fe;
//   }
//   case ET_HEX:
//   {
//     auto fe = new(alloc) H1HighOrderFE<ET_HEX> (order);
//     fe->SetVertexNumbers (dofmapy[dofmapx[ngel.Vertices()]]);
//     fe->SetOrderCell (order_inner[ei.Nr()]);
//     fe->ComputeNDof();
//     return *fe;
//   }
//  case ET_PRISM:
//   {
//     auto fe = new(alloc) H1HighOrderFE<ET_PRISM> (order);
//     fe->SetVertexNumbers (dofmapy[dofmapx[ngel.Vertices()]]);
//     fe->SetOrderCell (order_inner[ei.Nr()]);
//     fe->ComputeNDof();
//     return *fe;
//   }
//   case ET_PYRAMID:
//   {
//     auto fe = new(alloc) H1HighOrderFE<ET_PYRAMID> (order);
//     fe->SetVertexNumbers (dofmapy[dofmapx[ngel.Vertices()]]);
//     fe->SetOrderCell (order_inner[ei.Nr()]);
//     fe->ComputeNDof();
//     return *fe;
//   }
//   case ET_QUAD:
//   {
//     auto fe = new(alloc) H1HighOrderFE<ET_QUAD> (order);
//     fe->SetVertexNumbers (dofmapy[dofmapx[ngel.Vertices()]]);
//     if(ma->GetDimension()==2){
//       fe->SetOrderCell (order_inner[ei.Nr()]);
//       fe->ComputeNDof();
//     }
//     return *fe;
//   }
//   }
//   throw Exception("GetFE: undefined element");
// }
  

static RegisterFESpace<PeriodicH1Space> myinitifes ("h1ho_periodic");
 
