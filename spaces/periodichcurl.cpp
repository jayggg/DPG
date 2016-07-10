/*
  H(curl) high order FE Space with periodicity in two directions
*/

#include <comp.hpp>
using namespace ngcomp;

class PeriodicHCurlSpace : public HCurlHighOrderFESpace {

private:

  Array<int> vertmapx; // vertmapx[ slave_vertex# ] = master_vertex#
  Array<int> vertmapy; // (similarly for y)

  Array<int> dofmapx;  // dofmap[ slave_dof# ] = master_dof#
  Array<int> dofmapy;

  // Some more notes on Fx=dofmapx and Fy=dofmapy:

  // On non-identified dofs (not on periodic surfaces) Fx and Fy
  // behave like identity
  //       Fx(i) = i   and Fy(i) = i.
  //
  // When x=x0 and x=x1 faces are identified to be the same, and in
  // addition y=y0 and y=y1 faces are also identified, then the dofs
  // on the 4 edges where both x and y are constant are all identified
  // to be in the same equivalence class. In order to have a unique
  // master_dof# for elements of this class, we adopt the policy that
  //      Fx(i) is not greater than i, and 
  //      Fy(j) is not greater that j,
  // so master_dof# is always the lowest of the identified dof numbers.
  // Since Fx( Fy( i ) ) and  Fy( Fx( i ) ) are less than or equal to i,
  // any composition of these maps return the master dof number.
  
  
  Array<int> xid, yid;
  Array<double> xends, yends;

  void SetPeriodicIds();

  // int idx;             // periodic surface identification from mesh 
  // int idy;
  
  // Some more notes on idx and idy: 
  // 
  // In 2D: if segment i is copied to segment j, then idx = j.
  //
  // In 3D: index is the number of the "identify periodic " statement 
  //        in geo file.

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


PeriodicHCurlSpace :: PeriodicHCurlSpace (shared_ptr<MeshAccess> ama,
					  const Flags & flags)
  : HCurlHighOrderFESpace (ama, flags) {

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

  cout << " Initializing H(curl) space " << endl;
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


PeriodicHCurlSpace :: ~PeriodicHCurlSpace () {;}


void PeriodicHCurlSpace::SetPeriodicIds () {

  
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


void PeriodicHCurlSpace :: Update (LocalHeap & lh) {

  // Here we make dofmap, so that later when GetDofNrs is called by
  // periodic FEspace, dofmap will give the revised dof numbers
  // after the periodic identifications.
    
    
  // Warning: Note that GetDofNrs is called by the base class
  // HCurlHighOrderFESpace::Update, even before this derived class's
  // Update is called. Until dofmap is set, GetDofNrs must use the
  // base class GetDofNrs. This sentinel shows when dofmap is set:

  dofmapx.SetSize(0); // Sentinels: dofmaps are not yet set
  dofmapy.SetSize(0); 

  
  HCurlHighOrderFESpace::Update(lh);

  // Setting dofmap:
  
  // First set dofmaps to identity maps
  dofmapx.SetSize (GetNDof());
  for (int i = 0; i < dofmapx.Size(); i++)
    dofmapx[i] = i;
  
  dofmapy.SetSize (GetNDof());
  for (int i = 0; i < dofmapy.Size(); i++)
    dofmapy[i] = i;


  // Make a vertex-pair to edge hashtable
  HashTable<INT<2>, int> vp2e(ma->GetNEdges());

  for (int enr = 0; enr < ma->GetNEdges(); enr++)  {
    int v1, v2;
    ma->GetEdgePNums (enr, v1, v2);
    if (v1 > v2) Swap (v1, v2);
    vp2e[INT<2>(v1,v2)] = enr;
  }

  // Make a vertex-triple-or-quartet to face hashtable
  // (triangular faces get a dummy vertex number equal to -1)
  HashTable<INT<4>, int> v2f(ma->GetNFaces());

  Array<int> pnums;
  for (int fnr = 0; fnr < ma->GetNFaces(); fnr++)   {
    ma->GetFacePNums (fnr, pnums);
    INT<4> i4;
    if (pnums.Size() == 3) 
      i4 = {-1, pnums[0], pnums[1], pnums[2]};
    if (pnums.Size() == 4) 
      i4 = {pnums[0], pnums[1], pnums[2], pnums[3]};
    i4.Sort();
    v2f[i4] = fnr;
  }

  /* idx  **************************************************** */

  // Make a vertex slave -> master  array
  
  vertmapx.SetSize(ma->GetNV());
  for (int i = 0; i < vertmapx.Size(); i++)
    vertmapx[i] = i;

  Array<INT<2> > periodic_verts; // array of identified integer pairs 

  for (int idx : xid) {          // loop over all peroidic surface id nums

    /*  We are not sure if NGSolve policy is to give more than one
	periodic surface id in each direction, so for now we are
	looping over potential multiple x-periodic ids.
     */
    
    
    ma->GetPeriodicVertices(idx, periodic_verts);

    for (auto pair : periodic_verts) {
      int p1 = pair[1];
      int p0 = pair[0];
      if (p1<p0) Swap(p1,p0);
      vertmapx[p1] = p0;         // p0 is declared master here 
    }
    
    // Find periodic edges (using vertex-pair to edge hashtable)
  
    for (int enr = 0; enr < ma->GetNEdges(); enr++) {

      int v1, v2;
      ma->GetEdgePNums (enr, v1, v2);
    
      int mv1 = vertmapx[v1];                // master vertex numbers
      int mv2 = vertmapx[v2];
      if (v1 != mv1 && v2 != mv2) {          // edge shall be mapped
	if (mv1 > mv2) Swap (mv1, mv2);
	int menr = vp2e[INT<2>(mv1,mv2)];    // the master edge-nr
      
	// Edge numbers are the lowest dof numbers in H(curl) spaces,
	// hence we set dofmap for lowest order dofs  now:
	if (enr < menr) Swap (enr, menr);
	dofmapx[enr] = menr;
      
	// Note how the F(i) =< i policy is enforced above and below:
	// after the above swap, dofmapx[i] is less than or equal to i

	IntRange edofs = GetEdgeDofs (enr);   // dofs on slave edge
	IntRange medofs = GetEdgeDofs (menr); // dofs on master edge

	if ( edofs.First() < medofs.First() ) Swap(edofs, medofs);
      
	for (int i = 0; i < edofs.Size(); i++)
	  dofmapx[edofs[i]] = medofs[i];
      }
    }

    // Find periodic faces (using vertex-triple to face hashtable)
  
    for (int fnr = 0; fnr < ma->GetNFaces(); fnr++)  {

      ma->GetFacePNums (fnr, pnums);
      INT<4> i4;

      if (pnums.Size() == 3) 
	{
	  i4 = {-1, vertmapx[pnums[0]], vertmapx[pnums[1]], 
		vertmapx[pnums[2]]};
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
	  i4 = {vertmapx[pnums[0]], vertmapx[pnums[1]], vertmapx[pnums[2]], 
		vertmapx[pnums[3]]};  
	  if (i4[0] != pnums[0] && i4[1] != pnums[1] && i4[2] != pnums[2] && 
	      i4[3] != pnums[3])
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
  
  /* idy ********************************************************** */
  
  vertmapy.SetSize(ma->GetNV());
  for (int i = 0; i < vertmapy.Size(); i++)
    vertmapy[i] = i;

  for (int idy : yid) {
    ma->GetPeriodicVertices(idy, periodic_verts);

    for (auto pair : periodic_verts) {
      int p1 = pair[1];
      int p0 = pair[0];
      if (p1<p0) Swap(p1,p0);
      vertmapy[pair[1]] = pair[0];
    }

    // Find periodic edges
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

	    if (enr < menr) Swap (enr, menr);
	    dofmapy[enr] = menr;

	    IntRange edofs = GetEdgeDofs (enr);   // dofs on slave edge
	    IntRange medofs = GetEdgeDofs (menr); // dofs on master edge

	    if ( edofs.First() < medofs.First() ) Swap(edofs, medofs);
      
	    for (int i = 0; i < edofs.Size(); i++)
	      dofmapy[edofs[i]] = medofs[i];
	  }
      }

    // Find periodic faces
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

		if ( fdofs.First() < mfdofs.First() ) Swap(fdofs, mfdofs);

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



void PeriodicHCurlSpace::GetDofNrs (int elnr, Array<int> & dnums) const {

  // If dofmap is not yet set, then only do this:

  HCurlHighOrderFESpace::GetDofNrs(elnr,dnums);

  // If dofmap is set, then make the periodic adjustment:

  if (dofmapx.Size())
    for(int i=0; i<dnums.Size(); i++)
      dnums[i] = dofmapy[dofmapx[dnums[i]]];
}

void PeriodicHCurlSpace::GetSDofNrs (int selnr, Array<int> & dnums) const {

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

      INT<2> p(order_inner[ei.Nr()][0], order_inner[ei.Nr()][1]);
      FlatArray<INT<2> > of(1, &p);
      fe -> SetOrderFace (of);
      
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
  case ET_HEX:
  {
    auto fe = new(alloc) HCurlHighOrderFE<ET_HEX> (order);
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
  case ET_PYRAMID:
  {
    auto fe = new(alloc) HCurlHighOrderFE<ET_PYRAMID> (order);
    fe->SetVertexNumbers (vertmapy[vertmapx[ngel.Vertices()]]);
    fe->SetOrderCell (order_inner[ei.Nr()]);
    fe->ComputeNDof();
    return *fe;
  }
  case ET_SEGM:
  {
    auto fe = new(alloc) HCurlHighOrderFE<ET_SEGM> (order);
    fe->SetVertexNumbers (vertmapy[vertmapx[ngel.Vertices()]]);
    // fe->SetOrderCell (order_inner[ei.Nr()]);
    // fe->ComputeNDof();
    return *fe;
  }
  case ET_QUAD:
  {
    auto fe = new(alloc) HCurlHighOrderFE<ET_QUAD> (order);
    fe->SetVertexNumbers (vertmapy[vertmapx[ngel.Vertices()]]);
    if(ma->GetDimension()==2){
      fe->SetOrderCell (order_inner[ei.Nr()]);

      INT<2> p(order_inner[ei.Nr()][0], order_inner[ei.Nr()][1]);
      FlatArray<INT<2> > of(1, &p);
      fe -> SetOrderFace (of);

      fe->ComputeNDof();
    }
    return *fe;
  }
  }
  throw Exception("GetFE: undefined element");
}

  
// FiniteElement &
// HCurlHighOrderFESpace::GetFE (int elnr, LocalHeap & lh) const  {

//   Ngs_Element ngel = ma->GetElement(elnr);
//   ELEMENT_TYPE eltype = ngel.GetType();
    
//   switch (eltype)
//     {
//     case ET_SEGM:    return T_GetFE<ET_SEGM> (elnr, lh);
        
//     case ET_TRIG:    return T_GetFE<ET_TRIG> (elnr, lh);
//     case ET_QUAD:    return T_GetFE<ET_QUAD> (elnr, lh);
        
//     case ET_TET:     return T_GetFE<ET_TET> (elnr, lh);
//     case ET_PRISM:   return T_GetFE<ET_PRISM> (elnr, lh);
//     case ET_PYRAMID: return T_GetFE<ET_PYRAMID> (elnr, lh);
//     case ET_HEX:     return T_GetFE<ET_HEX> (elnr, lh);

//     default:
//       throw Exception ("illegal element in HCurlHoFeSpace::GetFE");
//     }
// }


// template <ELEMENT_TYPE ET>
// FiniteElement &
// HCurlHighOrderFESpace::T_GetFE (int elnr, LocalHeap & lh) const {
  
//   Ngs_Element ngel = ma->GetElement<ET_trait<ET>::DIM,VOL> (elnr);
//   fe->SetVertexNumbers (vertmapy[vertmapx[ngel.Vertices()]]);

//   if (!DefinedOn (ngel))
//     return * new (lh) HCurlDummyFE<ET>();

//   HCurlHighOrderFE<ET> * hofe =  new (lh) HCurlHighOrderFE<ET> ();
    
//   hofe -> SetVertexNumbers (ngel.vertices);
//   hofe -> SetOrderEdge (order_edge[ngel.Edges()]);
//   hofe -> SetUseGradEdge (usegrad_edge[ngel.Edges()]);

//   switch (int(ET_trait<ET>::DIM))
//     {
//     case 1:
//       throw Exception("no 1D elements in H(curl)");
//     case 2:
//       {
// 	hofe -> SetOrderCell (order_inner[elnr]);   // old style
// 	INT<2,TORDER> p(order_inner[elnr][0], order_inner[elnr][1]);
// 	FlatArray<INT<2,TORDER> > of(1, &p);
// 	hofe -> SetOrderFace (of);
          
// 	hofe -> SetUseGradCell (usegrad_cell[elnr]);  // old style
// 	FlatArray<bool> augf(1,&usegrad_cell[elnr]);
// 	hofe -> SetUseGradFace (augf); 
// 	break;
//       }
//     case 3:
//       {
// 	hofe -> SetOrderFace (order_face[ngel.Faces()]);
// 	hofe -> SetUseGradFace (usegrad_face[ngel.Faces()]);
          
// 	hofe -> SetOrderCell (order_inner[elnr]);
// 	hofe -> SetUseGradCell (usegrad_cell[elnr]); 
// 	break;
//       }
//     }
//   hofe -> SetType1 (type1);          
//   hofe -> ComputeNDof();
//   // hofe -> SetDiscontinuous(discontinuous);
  
//   return *hofe;
// }




static RegisterFESpace<PeriodicHCurlSpace> init ("hcurlho_periodic");



