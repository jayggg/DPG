#include <comp.hpp>
using namespace ngcomp;


///////////////////////////////////////////////////////////////
// The PROBLEM addressed by this code:
//
//   NGSolve's "L2HighOrderFESpace" does not have the ability to refer
//   to traces of its functions over the boundary of a domain. This
//   design was intentional, because theoretically one should not be
//   taking traces of L2 functions.
//
//   However, the code makes no distinction between methods that are
//   conforming in broken H^1, or conforming in L^2. Finite element
//   subspaces in either case are DG spaces. Methods conforming in
//   broken H^1 often need to integrate traces on global boundary, 
//   and this was not possible with the facilities in the code.
//
// The WORKAROUND for the problem implemented in this code:
// 
//   First define a new trace finite element. Here orientations 
//   must be carefully tracked.
//
//   Next, define an "L^2 high order FE space with traces", which 
//   we call "L2HighOrderFESpaceTrace" below, inheriting it from the
//   standard "L2HighOrderFESpace". The L2HighOrderFESpace class 
//   returns a dummy finite element whenever GetSFE() is called, 
//   but the  newly derived class's GetSFE() makes a new trace 
//   element and returns it.

// ------------------------------------------------------------


///////////////////////////////////////////////////////////////
// The first step: Define the "trace finite element" 
// of a volume element in DIM dimensions:

template <int DIM> 
class TraceElement : public ScalarFiniteElement<DIM-1> {

  // Element from which this trace element is obtained
  const ScalarFiniteElement<DIM> & vol_element;

  // Local facet number 
  int facnr;   

  ELEMENT_TYPE et;

  // Vertices of the vol_element
  Array<int> vnums;


  // To undertand the next data members, we need some background:
  // There are two possible orientations on any trace element:
  //
  // Facet orientation: The lowest vertex number is listed first, the
  // next highest is counted next, etc.
  //
  // Surface orientation: The vertices are listed in the order they
  // appear in the mesh.
  //
  // The following three data members should now make sense:


  // Trace element's vertices in surface orientation
  Array<int> svnums; 

  // trafo: maps points with DIM-1 coordinates in facet orientation to
  // give points in DIM coordinates of vol_element.
  Facet2ElementTrafo trafo; 

  // strafo: maps points with DIM-1 coordinates in facet orientation
  // to points with DIM-1 coordinates in surface orientation.
  Facet2SurfaceElementTrafo strafo; 


  // To understand the next data members, we need more concepts:

  // When constructor of trace element is called, we input vertex
  // numbers in the surface orientation (because at that point we have
  // these vertices in that orientation from mesh access). Hence
  // "svnums" has the surface orientation.
  
  // The vertices in svnums are mapped to a reference trace element.
  // Let's fix the reference trace element as follows:
  // In DIM=3 case: Triangle with points (0,0), (1,0), (0,1).
  // In DIM=2 case: Segment with points 0, 1.

  // All finite elements must have a CalcShape(ip,shape) routine which
  // evaluates its shape functions on a reference element point "ip". In
  // the trace element's case, "ip" is in reference element coordinates
  // with the surface orientation. We would like to have CalcShape
  // call the vol_element's CalcShape(). But to do this, we must make
  // these transformations:
  //    (DIM-1) coordinate S in surface orientation 
  //         --> (DIM-1) coordinate F facet orientation 
  //            ---> DIM-coordinate V volume element point.
  // ie, we need the maps S -> F -> V. 

  // The maps between S and F are made easy by matrices "a", "ainv"
  // and vector "b", which are the next data members.

  // Map S --> F  is  F = ainv*(S - b)
  // Map F --> S  is  S = a*F + b,   or simply "strafo".  

  Mat<DIM-1,DIM-1> a, ainv; 
  Vec<DIM-1> b;


public:
  


  TraceElement (const FiniteElement & ave,
                int afacnr,           // local facet # 
		ELEMENT_TYPE aet,     // elt type of surface elt
		Array<int> & avnums,  // vertex #s of volume elt
		Array<int> & asvnums); // vertex #s of  surface elt


  virtual ELEMENT_TYPE ElementType () const { return et; }


  virtual void CalcShape (const IntegrationPoint & ip, 
                          BareSliceVector<> shape) const ;

  virtual void CalcDShape (const IntegrationPoint & ip, 
			   BareSliceMatrix<> dshape) const  {

    throw Exception("CalcDShape not available for trace element");
    
  }

  virtual void Print(ostream & ost) const {

    ost << "A " << DIM-1 << "- dimensional trace element with "
	<< endl << "vertex numbers "<< endl<<svnums << endl;
  }
};



///////////////////////////////////////////////////////////////
// The second step: Make an L^2 high order finite element space
// with traces on the domain boundary:

class L2HighOrderFESpaceTrace : public L2HighOrderFESpace  {

public:
  L2HighOrderFESpaceTrace (shared_ptr<MeshAccess> ama, 
			   const Flags & flags, 
                           bool parseflags=false)
    : L2HighOrderFESpace (ama, flags, parseflags)  { 

    static ConstantCoefficientFunction one(1);

    integrator[BND] = // this integrator needed for visualization
      GetIntegrators().CreateBFI("robin", ma->GetDimension(), &one);
  }


  // Return surface element (needed for boundary integrals)
  //virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
  virtual const FiniteElement & GetSFE (ElementId sei, LocalHeap & lh) const;

  // Return dofs of a surface element
  virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

};





///////////////////////////////////////////////////////////////
// Member function definitions of "trace finite element" 

template <int DIM> 
TraceElement<DIM>::TraceElement(const FiniteElement & ave,
                int afacnr,           // local facet # 
		ELEMENT_TYPE aet,     // elt type of surface elt
		Array<int> & avnums,  // vertex #s of volume elt
		Array<int> & asvnums) // vertex #s of  surface elt
    : ScalarFiniteElement<DIM-1> (ave.GetNDof(), ave.Order()), 
    // // Note : For older versions use this constructor instead: 
    // ScalarFiniteElement<DIM-1> (aet,ave.GetNDof(),ave.Order()),
    vol_element(dynamic_cast<const ScalarFiniteElement<DIM>&>(ave)), 
    facnr(afacnr), et(aet), vnums(avnums), svnums(asvnums),
    trafo(vol_element.ElementType(), vnums),
    strafo(et, svnums)  { 

    // To make matrix "a" representing the map F --> S, namely 
    // S = a*F + b, we simply use the same map given by "strafo".
    //
    // DIM=3 case:    
    //     b =  strafo(0,0)
    //     a = [strafo(1,0)-strafo(0,0); strafo(1,0)-strafo(0,0)]
    // DIM=2 case:
    //     b =  strafo(0)
    //     a = [strafo(1)-strafo(0)]

    if (DIM == 3) {

      IntegrationPoint ip0(0.0,0.0), ip1(1.0,0.0), ip2(0.0,1.0);

      ip0 = strafo (ip0); 
      ip1 = strafo (ip1); 
      ip2 = strafo (ip2); 

      for (int j = 0; j < DIM-1; j++) {
	b(j)   = ip0(j);       
	a(j,0) = ip1(j)-ip0(j); 
	a(j,1) = ip2(j)-ip0(j);
      }
    }
    else if (DIM == 2)  {
      
      IntegrationPoint ip0(0.0), ip1(1.0);
      
      ip0 = strafo (ip0); 
      ip1 = strafo (ip1); 

      for (int j = 0; j < DIM-1; j++) {
	b(j)   = ip0(j);       
	a(j,0) = ip1(j)-ip0(j); 
      }
    }
    
    ainv = Inv(a);

}

template <int DIM> 
void TraceElement<DIM>::CalcShape (const IntegrationPoint & ip, 
				   BareSliceVector<> shape) const  {
    IntegrationPoint facet_ip;
    
    // We set up things to use the vol_element's CalcShape().
    // Input "ip" is in surface orientation, so we need to  
    // map  S (surface otn) -> F (facet otn)  -> V (volume coords).

    // First perform S --> F map using F = ainv*(S - b):
    Vec<DIM-1> vip_se = ip.Point(); 
    Vec<DIM-1> vip_facet = ainv * (vip_se - b); // F = vip_facet
    facet_ip = IntegrationPoint(vip_facet);

    // Next perform F --> V map using "trafo":
    IntegrationPoint vol_ip = trafo(facnr, facet_ip);

    // Compute vol_element shape function values and return:
    vol_element.CalcShape (vol_ip, shape);
    
}



///////////////////////////////////////////////////////////////
// Member function definitions of L^2 high order finite element 
// space with traces on the domain boundary:

const FiniteElement &  L2HighOrderFESpaceTrace::
GetSFE (ElementId sei, LocalHeap & lh) const {

	ArrayMem<int,10> fnums, elnums, vnums, svnums;
	ELEMENT_TYPE et = ma->GetElType (sei);

	fnums = ma->GetElFacets(sei);  /* fnums = facet numbers of 
				     surface elt number selnr */
	int fac = fnums[0];
	ma->GetFacetElements(fac,elnums);/* elnums = elt numbers of the elt
				     sharing facet number fac */
	int el = elnums[0];
    auto ei = ngfem::ElementId(VOL, el);
	fnums = ma->GetElFacets(ei);  /* fnums = facet numbers of 
				     elt number el */
	const FiniteElement & fel = GetFE (ei, lh);
	int facnr = 0;                  /* facnr = local facet number of
				     facet numbered fac globally */
	for (int k=0; k<fnums.Size(); k++)
	if(fac==fnums[k]) facnr = k;

	vnums = ma->GetElVertices (ei);     
	svnums = ma->GetElVertices (sei);     

	if (ma->GetDimension() == 2)
	return *new (lh) TraceElement<2> (fel, facnr, et, vnums, svnums);
	else {
		return *new (lh) TraceElement<3> (fel, facnr, et, vnums, svnums);
	}
}



void  L2HighOrderFESpaceTrace::
GetSDofNrs(int selnr, Array<int> & dnums) const {

    Array<int> fnums, elnums;
    fnums = ma->GetElFacets(ngfem::ElementId(BND, selnr));
    int fac = fnums[0];
    ma->GetFacetElements(fac,elnums);

    int el = elnums[0];
    if (elnums.Size() != 1)
      cerr << "surface element not at outer boundary" << endl;
    GetDofNrs (el, dnums);
}



static RegisterFESpace<L2HighOrderFESpaceTrace> init ("l2ho_trace");


