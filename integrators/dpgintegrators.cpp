#include <fem.hpp>
#include "dpgintegrators.hpp"

// See end of file for all integrators provided


using namespace ngsolve;

namespace dpg {

  //////////////////////////////////////////////////////////////
  // Integrate a(x)*grad u . grad e, where u and e are in different spaces

  template<int D> template <class SCAL>
  void GradGrad<D>::T_CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<SCAL> elmat,
			    LocalHeap & lh) const {
    
    const CompoundFiniteElement &  cfel  // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const ScalarFiniteElement<D> & fel_u =  // u space
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[GetInd1()]);
    const ScalarFiniteElement<D> & fel_e =  // e space
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[GetInd2()]);

    elmat = SCAL(0.0);

    // u dofs [ru.First() : ru.Next()-1],  e dofs [re.First() : re.Next()-1]
    IntRange ru = cfel.GetRange(GetInd1()); 
    IntRange re = cfel.GetRange(GetInd2()); 
    int ndofe = re.Size();
    int ndofu = ru.Size();

    FlatMatrixFixWidth<D> dum(ndofu,lh); // to store grad(u-basis)  
    FlatMatrixFixWidth<D> dem(ndofe,lh); // to store grad(e-basis)

    ELEMENT_TYPE eltype                  // get the type of element: 
      = fel_u.ElementType();             // ET_TRIG in 2d, ET_TET in 3d.

    const IntegrationRule &              // Note: p = fel_u.Order()-1
      ir = SelectIntegrationRule(eltype, fel_u.Order()+fel_e.Order()-2);
    
    FlatMatrix<SCAL> submat(ndofe,ndofu,lh);
    submat = 0.0;

    for(int k=0; k<ir.GetNIP(); k++) {	
      
      MappedIntegrationPoint<D,D> mip (ir[k],eltrans);
      // set grad(u-basis) and grad(e-basis) at mapped pts in dum and dem.
      fel_u.CalcMappedDShape( mip, dum ); 
      fel_e.CalcMappedDShape( mip, dem );

      // evaluate coefficient
      SCAL fac = coeff_a -> T_Evaluate<SCAL>(mip);
      fac *= mip.GetWeight() ;
      
      //             [ndofe x D] * [D x ndofu]
      submat +=  fac *  dem     * Trans(dum) ;
    }
    
    elmat.Rows(re).Cols(ru) += submat;
    if (GetInd1() != GetInd2())
      elmat.Rows(ru).Cols(re) += Conj(Trans(submat));
    
  }


  //////////////////////////////////////////////////////////////
  // Integrate d(x) * q.n * e on all element boundaries
  // where q is an Hdiv space and e is in a scalar space,
  // and d is a complex or real coefficient.

  template<int D> template <class SCAL>
  void FluxTrace<D>::T_CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<SCAL> elmat,
			    LocalHeap & lh) const {
    
    const CompoundFiniteElement &  cfel  // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const HDivFiniteElement<D>   & fel_q =  // q space
      dynamic_cast<const HDivFiniteElement<D>&  > (cfel[GetInd1()]);
    const ScalarFiniteElement<D> & fel_e =  // e space
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[GetInd2()]);

    elmat = SCAL(0.0);

    IntRange rq = cfel.GetRange(GetInd1()); 
    IntRange re = cfel.GetRange(GetInd2());
    int ndofq = rq.Size();
    int ndofe = re.Size();

    FlatMatrix<SCAL> submat(ndofe,ndofq,lh);
    FlatMatrixFixWidth<D> shapeq(ndofq,lh);  // q-basis (vec) values
    FlatVector<>          shapee(ndofe,lh);  // e-basis basis 
    
    ELEMENT_TYPE eltype                      // get the type of element: 
      = fel_q.ElementType();                 // ET_TRIG in 2d, ET_TET in 3d.

    // transform facet integration points to volume integration points
    Facet2ElementTrafo transform(eltype);

    int nfa = ElementTopology::GetNFacets(eltype); /* nfa = number of facets
						      of an element */    
    submat = 0.0;

    for(int k = 0; k<nfa; k++) {
      
      // type of geometry of k-th facet
      ELEMENT_TYPE eltype_facet = ElementTopology::GetFacetType(eltype, k); 
      
      const IntegrationRule & facet_ir =
	SelectIntegrationRule (eltype_facet, fel_q.Order()+fel_e.Order()); 

      // reference element normal vector
      FlatVec<D> normal_ref = ElementTopology::GetNormals(eltype) [k]; 

      for (int l = 0; l < facet_ir.GetNIP(); l++) {

	// map 1D facet points to volume integration points
	IntegrationPoint volume_ip = transform(k, facet_ir[l]);
	MappedIntegrationPoint<D,D> mip (volume_ip, eltrans);
	
	// compute normal on physcial element
	Mat<D> inv_jac = mip.GetJacobianInverse();
	double det = mip.GetJacobiDet();
	Vec<D> normal = fabs(det) * Trans(inv_jac) * normal_ref;       
	double len = L2Norm(normal);
	normal /= len;
	double weight = facet_ir[l].Weight()*len;
	
	// mapped H(div) basis fn values and DG fn (no need to map) values
	fel_q.CalcMappedShape(mip,shapeq); 
	fel_e.CalcShape(volume_ip,shapee); 
	
	// evaluate coefficient
	SCAL dd = coeff_d -> T_Evaluate<SCAL>(mip);

	//                   [ndofe x 1]      [ndofq x D] *  [D x 1] 	
	submat += (dd*weight) * shapee * Trans( shapeq    *  normal ) ;
      }
    }
    elmat.Rows(re).Cols(rq) += submat;
    elmat.Rows(rq).Cols(re) += Conj(Trans(submat));
  }


  //////////////////////////////////////////////////////////////
  // Integrate a(x)* u * e, where u and e are in different spaces

  template<int D> template <class SCAL>
  void EyeEye<D>::T_CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<SCAL> elmat,
		     LocalHeap & lh) const {
   

    const CompoundFiniteElement &  cfel  // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);
    const ScalarFiniteElement<D> & fel_u =  // u space
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[GetInd1()]);
    const ScalarFiniteElement<D> & fel_e =  // e space
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[GetInd2()]);

    elmat = SCAL(0.0);

    IntRange ru = cfel.GetRange(GetInd1()); 
    IntRange re = cfel.GetRange(GetInd2()); 
    int ndofe = re.Size();
    int ndofu = ru.Size();

    Vector<> ushape(ndofu);
    Vector<> eshape(ndofe);

    ELEMENT_TYPE eltype = fel_u.ElementType();      
    const IntegrationRule &         
      ir = SelectIntegrationRule(eltype, fel_u.Order()+fel_e.Order());
    FlatMatrix<SCAL> submat(ndofe,ndofu,lh);
    submat = SCAL(0.0);

    for(int k=0; k<ir.GetNIP(); k++) {	
      
      MappedIntegrationPoint<D,D> mip (ir[k],eltrans);

      fel_u.CalcShape( ir[k], ushape ); 
      fel_e.CalcShape( ir[k], eshape );

      SCAL fac = (coeff_a -> T_Evaluate<SCAL>(mip))* mip.GetWeight() ;
      //               [ndofe x D] * [D x ndofu]
      submat +=  fac *  eshape     * Trans(ushape) ;
    }
    
    elmat.Rows(re).Cols(ru) += submat;
    if (GetInd1() != GetInd2())
      elmat.Rows(ru).Cols(re) += Conj(Trans(submat));
  }
   

  //////////////////////////////////////////////////////////////
  // Integrate c(x) * u * v over ELEMENT BOUNDARIES, where u and v 
  // are in different spaces

  template<int D> template <class SCAL>
  void TraceTrace<D>::T_CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<SCAL> elmat,
			    LocalHeap & lh) const {
    
    const CompoundFiniteElement &  cfel  // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const ScalarFiniteElement<D> & fel_u =  // u space
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[GetInd1()]);
    const ScalarFiniteElement<D> & fel_e =  // e space
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[GetInd2()]);

    elmat = SCAL(0.0);

    IntRange ru = cfel.GetRange(GetInd1()); 
    IntRange re = cfel.GetRange(GetInd2()); 
    int ndofe = re.Size();
    int ndofu = ru.Size();
    FlatVector<>      shapee(ndofe,lh);  
    FlatVector<>      shapeu(ndofu,lh);  
    FlatMatrix<SCAL>  submat(ndofe,ndofu, lh);  
    submat = SCAL(0.0);

    ELEMENT_TYPE eltype = fel_u.ElementType();         
    Facet2ElementTrafo transform(eltype);
    int nfa = ElementTopology :: GetNFacets(eltype); 

    for(int k = 0; k<nfa; k++) {
      
      // type of geometry of k-th facet
      ELEMENT_TYPE eltype_facet = ElementTopology::GetFacetType(eltype, k); 
      
      const IntegrationRule & facet_ir =
	SelectIntegrationRule (eltype_facet, fel_u.Order()+fel_e.Order()); 

      // reference element normal vector
      FlatVec<D> normal_ref = ElementTopology::GetNormals(eltype) [k]; 

      for (int l = 0; l < facet_ir.GetNIP(); l++) {

	// map 1D facet points to volume integration points
	IntegrationPoint volume_ip = transform(k, facet_ir[l]);
	MappedIntegrationPoint<D,D> mip (volume_ip, eltrans);
	
	// compute normal on physcial element
	Mat<D> inv_jac = mip.GetJacobianInverse();
	double det = mip.GetJacobiDet();
	Vec<D> normal = fabs(det) * Trans(inv_jac) * normal_ref;       
	double len = L2Norm(normal);
	normal /= len;
	double weight = facet_ir[l].Weight()*len;
	
	fel_e.CalcShape(volume_ip,shapee); 
	fel_u.CalcShape(volume_ip,shapeu); 

	SCAL cc = coeff_c  -> T_Evaluate<SCAL>(mip);

	//                     [ndofe x 1]  [1 x ndofu] 
	submat +=  (cc*weight) * shapee * Trans(shapeu);
      }
    }
    
    elmat.Rows(re).Cols(ru) += submat;
    if (GetInd1() != GetInd2())
      elmat.Rows(ru).Cols(re) += Conj(Trans(submat));    
  }


  //////////////////////////////////////////////////////////////
  // Integrate c(x) * q.n * r.n over the GLOBAL BOUNDARY, where 
  // q and r may be from different H(div)-type spaces.

  template<int D> template <class SCAL>
  void FluxFluxBoundary<D> ::
  T_CalcElementMatrix (const FiniteElement & base_fel,
		       const ElementTransformation & eltrans, 
		       FlatMatrix<SCAL> elmat,
		       LocalHeap & lh) const {

    const CompoundFiniteElement &  cfel      // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);
    
    // This FE is already multiplied by normal:
    const HDivNormalFiniteElement<D-1> & fel_q = // q.n space
      dynamic_cast<const HDivNormalFiniteElement<D-1>&> (cfel[GetInd1()]);

    const HDivNormalFiniteElement<D-1> & fel_r = // r.n space
      dynamic_cast<const HDivNormalFiniteElement<D-1>&> (cfel[GetInd2()]);
    
    elmat = SCAL(0.0);

    IntRange rq = cfel.GetRange(GetInd1());
    IntRange rr = cfel.GetRange(GetInd2());
    int ndofq = rq.Size();
    int ndofr = rr.Size();
 
    FlatMatrix<SCAL> submat(ndofr, ndofq, lh);  
    submat = SCAL(0.0);
    
    FlatVector<> qshape(fel_q.GetNDof(), lh);
    FlatVector<> rshape(fel_r.GetNDof(), lh);

    const IntegrationRule ir(fel_q.ElementType(), 
			     fel_q.Order() + fel_r.Order());

    for (int i = 0 ; i < ir.GetNIP(); i++) {

      MappedIntegrationPoint<D-1,D> mip(ir[i], eltrans);

      SCAL cc = coeff_c -> T_Evaluate<SCAL>(mip);

      fel_r.CalcShape (ir[i], rshape);
      fel_q.CalcShape (ir[i], qshape);
      // mapped q.n-shape is simply reference q.n-shape / measure
      qshape *= 1.0/mip.GetMeasure();
      rshape *= 1.0/mip.GetMeasure();
      //                              [ndofr x 1]  [1 x ndofq]
      submat += (cc*mip.GetWeight()) * rshape  * Trans(qshape);
    }       

    elmat.Rows(rr).Cols(rq) += submat;
    if (GetInd1() != GetInd2())
      elmat.Rows(rq).Cols(rr) += Conj(Trans(submat));
  }

  /////////////////////////////////////////////////////////////////
  // Integrate c(x) * u * e over the GLOBAL BOUNDARY, where u and e
  // are from different spaces in general. This works whenever the u
  // and the e spaces have "surface elements", i.e., elements on the
  // boundary facets that represent traces of u and e on the global
  // boundary. 
  //
  // E.g., NGSolve's L2 finite element spaces do not have such surface
  // elements, but H1 spaces do and so does the modified
  // "L2HighOrderFESpaceTrace" space in this library.
  //
  // If one of the spaces is L2HighOrderFESpaceTrace, then 
  //     -eliminate_internal
  // flag should NOT be used!
  // 
  // (See also next integrator for a different approach.)
  //
  template<int D> template <class SCAL>
  void TraceTraceBoundary<D> ::
  T_CalcElementMatrix (const FiniteElement & base_fel,
		       const ElementTransformation & eltrans, 
		       FlatMatrix<SCAL> elmat,
		       LocalHeap & lh) const {

    const CompoundFiniteElement &  cfel      // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);

    // get surface elements
    const ScalarFiniteElement<D-1> & fel_u = // u space
      dynamic_cast<const ScalarFiniteElement<D-1>&> (cfel[GetInd1()]);
    const ScalarFiniteElement<D-1> & fel_e = // u space
      dynamic_cast<const ScalarFiniteElement<D-1>&> (cfel[GetInd2()]);

    elmat = SCAL(0.0);

    IntRange ru = cfel.GetRange(GetInd1());
    IntRange re = cfel.GetRange(GetInd2());
    int ndofu = ru.Size();
    int ndofe = re.Size();
 
    FlatMatrix<SCAL> submat(ndofe, ndofu, lh);  
    submat = SCAL(0.0);
    
    FlatVector<> ushape(fel_u.GetNDof(), lh);
    FlatVector<> eshape(fel_e.GetNDof(), lh);

    const IntegrationRule ir(fel_u.ElementType(), 
			     fel_u.Order() + fel_e.Order());

    for (int i = 0 ; i < ir.GetNIP(); i++) {

      MappedIntegrationPoint<D-1,D> mip(ir[i], eltrans);

      SCAL cc = coeff_c -> T_Evaluate<SCAL>(mip);

      fel_u.CalcShape (ir[i], ushape);
      fel_e.CalcShape (ir[i], eshape);
      //                             [ndofe x 1]  [1 x ndofu]
      submat += (cc*mip.GetWeight()) * eshape  * Trans(ushape);
    }       

    elmat.Rows(re).Cols(ru) += submat;
    if (GetInd1() != GetInd2())
      elmat.Rows(ru).Cols(re) += Conj(Trans(submat));
  }


  /////////////////////////////////////////////////////////////////
  // The next integrator also integrates c(x) * u * e * ds over the
  // GLOBAL BOUNDARY. The difference is that now u and e are in finite
  // element spaces which may not have "surface elements" that
  // represent their global boundary traces. The idea is to integrate
  // on those facet elements of a volume element that lie on the
  // global boundary.
  //
  // This integrator works correctly with "-eliminate_internal" flag.
  //
  // This integrator currently does not work for general coefficients
  // (because we haven't implemented a mechanism to take a volume
  // coefficient and make it into a boundary coefficient) but it
  // works, for example, for constant coefficients.
  //
  
  template<int D> template <class SCAL>
  void RobinVolume<D> ::
  T_CalcElementMatrix (const FiniteElement & base_fel,
                       const ElementTransformation & eltrans, 
                       FlatMatrix<SCAL> elmat,
                       LocalHeap & lh) const {
    
    ELEMENT_TYPE eltype                
      = base_fel.ElementType();        
    const CompoundFiniteElement &  cfel     // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);

    // note how we do NOT refer to D-1 elements here:
    const ScalarFiniteElement<D> & fel_u =  // u space
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[GetInd1()]);
    const ScalarFiniteElement<D> & fel_e =  // e space
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[GetInd2()]);
    
    elmat = SCAL(0);
    IntRange ru = cfel.GetRange(GetInd1());
    IntRange re = cfel.GetRange(GetInd2());
    int ndofe = re.Size();
    int ndofu = ru.Size();
            
    FlatVector<> ushape(fel_u.GetNDof(), lh);
    FlatVector<> eshape(fel_e.GetNDof(), lh);
    FlatMatrix<SCAL> submat(ndofe,ndofu,lh);
    submat = SCAL(0);

    int nfacet = ElementTopology::GetNFacets(eltype);
    Facet2ElementTrafo transform(eltype); 
    FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);    
    const MeshAccess & ma = *(const MeshAccess*)eltrans.GetMesh();

    Array<int> fnums, sels;
    fnums = ma.GetElFacets (eltrans.GetElementId());
      
    for (int k = 0; k < nfacet; k++)    {

      ma.GetFacetSurfaceElements (fnums[k], sels);

      // if interior element, then do nothing:
      if (sels.Size() == 0) continue; 

      // else: 

      Vec<D> normal_ref = normals[k];

      ELEMENT_TYPE etfacet=ElementTopology::GetFacetType(eltype, k);

      IntegrationRule ir_facet(etfacet, fel_e.Order()+fel_u.Order());
      
      // map the facet integration points to volume reference elt ipts
      IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
      // ... and further to the physical element 
      MappedIntegrationRule<D,D> mir(ir_facet_vol, eltrans, lh);
        
      for (int i = 0 ; i < ir_facet_vol.GetNIP(); i++) {
	
	SCAL val = coeff_c->T_Evaluate<SCAL> (mir[i]);

	// this is contrived to get the surface measure in "len"
	Mat<D> inv_jac = mir[i].GetJacobianInverse();
	double det = mir[i].GetMeasure();
	Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
	double len = L2Norm (normal);    

	val *= len * ir_facet[i].Weight();
	
	fel_u.CalcShape (ir_facet_vol[i], ushape);
	fel_e.CalcShape (ir_facet_vol[i], eshape);
        
	submat += val * eshape * Trans(ushape);
      }    
    }
    elmat.Rows(re).Cols(ru) += submat;
    if (GetInd1() != GetInd2())
      elmat.Rows(ru).Cols(re) += Conj(Trans(submat));
  }


  /////////////////////////////////////////////////////////////////
  // A source integrator (to accompany the last integrator) that
  // integrates (G.n + g) * e * ds over the GLOBAL BOUNDARY, where e
  // may be in a finite element space without "surface elements", and
  // g is a volume coefficient (which will be evaluated only on the
  // boundary). The integrator is called in this form in 2D
  //    neumannvol <ind> <g> <Gx> <Gy>  
  // and 
  //    neumannvol <ind> <g> <Gx> <Gy> <Gz>
  // in 3D. 
  //
  // This integrator currently does not work for coefficients on
  // multiple boundary parts (because we haven't implemented a
  // mechanism to take a volume coefficient and make it into a
  // boundary coefficient).
  //

  template<int D> template <class SCAL>
  void NeumannVolume<D> ::
  T_CalcElementVector (const FiniteElement & base_fel,
		       const ElementTransformation & eltrans, 
		       FlatVector<SCAL> elvec,
		       LocalHeap & lh) const {

    const CompoundFiniteElement &  cfel  
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const ScalarFiniteElement<D> & fel = 
      dynamic_cast<const ScalarFiniteElement<D>&> (cfel[indx]);

    FlatVector<> ushape(fel.GetNDof(), lh);
    elvec = SCAL(0);    
    IntRange re = cfel.GetRange(indx);
    int ndofe = re.Size();
    FlatVector<SCAL> subvec(ndofe,lh);
    subvec = SCAL(0);

    const IntegrationRule ir(fel.ElementType(), 2*fel.Order());
    ELEMENT_TYPE eltype = base_fel.ElementType();        
    int nfacet = ElementTopology::GetNFacets(eltype);
    Facet2ElementTrafo transform(eltype); 
    FlatVector< Vec<D> > normals = ElementTopology::GetNormals<D>(eltype);

    const MeshAccess & ma = *(const MeshAccess*)eltrans.GetMesh();

    Array<int> fnums, sels;
    fnums = ma.GetElFacets (eltrans.GetElementId());

    for (int k = 0; k < nfacet; k++)    {

      ma.GetFacetSurfaceElements (fnums[k], sels);

      // if interior element, then do nothing:
      if (sels.Size() == 0) continue; 

      // else: 

      Vec<D> normal_ref = normals[k];

      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);

      IntegrationRule ir_facet(etfacet, 2*fel.Order());
      
      // map the facet integration points to volume reference elt ipts
      IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
      // ... and further to the physical element 
      MappedIntegrationRule<D,D> mir(ir_facet_vol, eltrans, lh);

      for (int i = 0 ; i < ir_facet_vol.GetNIP(); i++) {
	
       	SCAL G[3] ;
	G[0] = coeff_Gx -> T_Evaluate<SCAL>(mir[i]);
	G[1] = coeff_Gy -> T_Evaluate<SCAL>(mir[i]);
	if (D==3)  G[2] = coeff_Gz -> T_Evaluate<SCAL>(mir[i]);
	FlatVector<SCAL> Gval(D,lh);	
	for (int dd=0; dd<D; dd++)  Gval[dd] = G[dd];
	SCAL g = coeff_g -> T_Evaluate<SCAL>(mir[i]);

	// this is contrived to get the surface measure in "len"
	Mat<D> inv_jac = mir[i].GetJacobianInverse();
	double det = mir[i].GetMeasure();
	Vec<D> normal = det * Trans (inv_jac) * normal_ref;       
	double len = L2Norm (normal);    
	
	SCAL gg = (InnerProduct(Gval,normal) + g*len)
	          * ir_facet[i].Weight();
		
	fel.CalcShape (ir_facet_vol[i], ushape);
	        
	subvec += gg * ushape;
      }   
    }
    elvec.Rows(re) += subvec;
  }


  /////////////////////////////////////////////////////////////////
  // Integrate c(x) * q.n * w over the GLOBAL BOUNDARY, were q is in
  // any H(div)-like space and w is in a scalar finite element space.
  //
  template<int D> template <class SCAL>
  void FluxTraceBoundary<D> ::
  T_CalcElementMatrix (const FiniteElement & base_fel,
		       const ElementTransformation & eltrans, 
		       FlatMatrix<SCAL> elmat,
		       LocalHeap & lh) const {

    const CompoundFiniteElement &  cfel      // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);

    // This FE is already multiplied by normal:
    const HDivNormalFiniteElement<D-1> & fel_q = // q.n space
      dynamic_cast<const HDivNormalFiniteElement<D-1>&> (cfel[GetInd1()]);

    const ScalarFiniteElement<D-1> & fel_w =     // w space
      dynamic_cast<const ScalarFiniteElement<D-1>&> (cfel[GetInd2()]);

    elmat = SCAL(0.0);

    IntRange rq = cfel.GetRange(GetInd1());
    IntRange rw = cfel.GetRange(GetInd2());
    int ndofq = rq.Size();
    int ndofw = rw.Size();
 
    FlatMatrix<SCAL> submat(ndofw, ndofq, lh);  
    submat = SCAL(0.0);
    
    FlatVector<> qshape(fel_q.GetNDof(), lh);
    FlatVector<> wshape(fel_w.GetNDof(), lh);

    const IntegrationRule ir(fel_q.ElementType(), 
			     fel_q.Order() + fel_w.Order());

    for (int i = 0 ; i < ir.GetNIP(); i++) {

      MappedIntegrationPoint<D-1,D> mip(ir[i], eltrans);

      SCAL cc = coeff_c -> T_Evaluate<SCAL>(mip);

      fel_q.CalcShape (ir[i], qshape);
      // mapped q.n-shape is simply reference q.n-shape / measure
      qshape *= 1.0/mip.GetMeasure();
      fel_w.CalcShape (ir[i], wshape);

      submat += (cc*mip.GetWeight()) * wshape  * Trans(qshape);
    }       

    elmat.Rows(rw).Cols(rq) += submat;
    elmat.Rows(rq).Cols(rw) += Conj(Trans(submat));
  }


  //////////////////////////////////////////////////////////////


  static RegisterBilinearFormIntegrator<GradGrad<2>> 
  initlaplacedpg2d("gradgrad", 2, 3);
  static RegisterBilinearFormIntegrator<GradGrad<3>> 
  initlaplacedpg3d("gradgrad", 3, 3);

  static RegisterBilinearFormIntegrator<FluxTrace<2>> 
  initflxtrc2d("flxtrc", 2, 3);
  static RegisterBilinearFormIntegrator<FluxTrace<3>> 
  initflxtrc3d("flxtrc", 3, 3);

  static RegisterBilinearFormIntegrator<EyeEye<2>> 
  initmassdpg2d("eyeeye", 2, 3);
  static RegisterBilinearFormIntegrator<EyeEye<3>> 
  initmassdpg3d("eyeeye", 3, 3);

  static RegisterBilinearFormIntegrator<TraceTrace<2>> 
  inittrcdpg2d("trctrc", 2, 3);
  static RegisterBilinearFormIntegrator<TraceTrace<3>> 
  inittrcdpg3d("trctrc", 3, 3);

  static RegisterLinearFormIntegrator<NeumannVolume<2>> 
  initneumanndg2d("neumannvol", 2, 4);
  static RegisterLinearFormIntegrator<NeumannVolume<3>> 
  initneumanndg3d("neumannvol", 3, 5);

  static RegisterBilinearFormIntegrator<RobinVolume<2>> 
  initrobinvol2d("robinvol", 2, 3);
  static RegisterBilinearFormIntegrator<RobinVolume<2>> 
  initrobinvol3d("robinvol", 3, 3);

  static RegisterBilinearFormIntegrator<FluxFluxBoundary<2>> 
  initflxflxb2d("flxflxbdry", 2, 3);
  static RegisterBilinearFormIntegrator<FluxFluxBoundary<3>> 
  initflxflxb3d("flxflxbdry", 3, 3);

  static RegisterBilinearFormIntegrator<TraceTraceBoundary<2>> 
  inittrctrcb2d("trctrcbdry", 2, 3);
  static RegisterBilinearFormIntegrator<TraceTraceBoundary<3>> 
  inittrctrcb3d("trctrcbdry", 3, 3);

  static RegisterBilinearFormIntegrator<FluxTraceBoundary<2>> 
  initflxtrcb2d("flxtrcbdry", 2, 3);
  static RegisterBilinearFormIntegrator<FluxTraceBoundary<3>> 
  initflxtrcb3d("flxtrcbdry", 3, 3);


} // closing namespace dpg

