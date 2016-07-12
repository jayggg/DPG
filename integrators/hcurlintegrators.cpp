#include <fem.hpp>
#include "dpgintegrators.hpp"


// See end of file for all integrators provided

using namespace ngsolve;

namespace dpg {

  //////////////////////////////////////////////////////////////
  // Integrate a(x) * curl U . curl V, where U and V are in 
  // Hcurl spaces, and a(x) is a complex or real coefficient


  template<int D> class CurlCurlPG : public DPGintegrator {

    shared_ptr<CoefficientFunction> coeff_a;

    template<class SCAL>
    void T_CalcElementMatrix (const FiniteElement & base_fel,
			      const ElementTransformation & eltrans, 
			      FlatMatrix<SCAL> elmat,
			      LocalHeap & lh)  const ;        
  public:
    
    CurlCurlPG(const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : DPGintegrator(coeffs), coeff_a(coeffs[2])  {

      cout << "Using DPG integrator " << Name() << " with components "
	   << GetInd1()+1 << " and " << GetInd2()+1 << endl ; 
    }

    virtual bool IsSymmetric() const { return !coeff_a->IsComplex() ; }
    
    virtual string Name () const { return "CurlCurlPG"; }

    virtual bool BoundaryForm () const { return false; }

    void CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<double> elmat,
			    LocalHeap & lh) const {
      T_CalcElementMatrix<double>(base_fel,eltrans,elmat,lh);
					       
    }
    void CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<Complex> elmat,
			    LocalHeap & lh) const {
     T_CalcElementMatrix<Complex>(base_fel,eltrans,elmat, lh);    
    }
  };


  template<int D> template <class SCAL>
  void CurlCurlPG<D>::T_CalcElementMatrix (const FiniteElement & base_fel,
					const ElementTransformation & eltrans, 
					FlatMatrix<SCAL> elmat,
					LocalHeap & lh) const {

    const CompoundFiniteElement &  cfel  // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const HCurlFiniteElement<D> & fel_u =  // U space
      dynamic_cast<const HCurlFiniteElement<D>&> (cfel[GetInd1()]);
    const HCurlFiniteElement<D> & fel_v =  // V space
      dynamic_cast<const HCurlFiniteElement<D>&> (cfel[GetInd2()]);

    elmat = SCAL(0.0);

    // U dofs [ru.First() : ru.Next()-1],  v dofs [rv.First() : rv.Next()-1]
    IntRange ru = cfel.GetRange(GetInd1()); 
    IntRange rv = cfel.GetRange(GetInd2()); 
    int ndofu = ru.Size();
    int ndofv = rv.Size();

    FlatMatrixFixWidth<D> curl_um(ndofu,lh); // to store curl(U-basis)
    FlatMatrixFixWidth<D> curl_vm(ndofv,lh); // to store curl(V-basis)

    ELEMENT_TYPE eltype                  // get the type of element: 
      = fel_u.ElementType();             // ET_TET in 3d.

    const IntegrationRule &              // Note: p = fel_u.Order()-1
      ir = SelectIntegrationRule(eltype, fel_u.Order()+fel_v.Order()-2);
    
    FlatMatrix<SCAL> submat(ndofv,ndofu,lh);
    submat = SCAL(0.0);

    for(int k=0; k<ir.GetNIP(); k++) {	
      
      MappedIntegrationPoint<D,D> mip (ir[k],eltrans);
      // set curl(U-basis) and curl(V-basis) at mapped pts in curl_um and curl_vm.
      fel_u.CalcMappedCurlShape( mip, curl_um ); 
      fel_v.CalcMappedCurlShape( mip, curl_vm );

      // evaluate coefficient
      SCAL fac = coeff_a -> T_Evaluate<SCAL>(mip);
      fac *= mip.GetWeight() ;
      
      //             [ndofv x D] * [D x ndofu]
      submat +=  fac * curl_vm  * Trans(curl_um) ;
    }
    
    elmat.Rows(rv).Cols(ru) += submat;

    if (GetInd1() != GetInd2())
      elmat.Rows(ru).Cols(rv) += Conj(Trans(submat));
  }


  /////////////////////////////////////////////////////////////////
  // Integrate d(x) * H . (F x n)  on all ELEMENT BOUNDARIES
  // where H and F are in Hcurl spaces and d is a complex or real
  // coefficient.


  template<int D> class TraceTraceXn : public DPGintegrator {
    
    shared_ptr<CoefficientFunction>  coeff_d;
    
    template<class SCAL>
    void T_CalcElementMatrix (const FiniteElement & base_fel,
			      const ElementTransformation & eltrans, 
			      FlatMatrix<SCAL> elmat,
			      LocalHeap & lh)  const ;
    
  public:
    
    TraceTraceXn(const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : DPGintegrator(coeffs), coeff_d(coeffs[2])  {

      cout << "Using DPG integrator " << Name() << " with components "
	   << GetInd1()+1 << " and " << GetInd2()+1 << endl ; 
    }

    virtual bool IsSymmetric() const { return !coeff_d->IsComplex() ; }
    
    virtual string Name () const { return "TraceTraceXn"; }

    virtual bool BoundaryForm () const { return false; }

    void CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<double> elmat,
			    LocalHeap & lh) const {
      T_CalcElementMatrix<double>(base_fel,eltrans,elmat,lh);
					       
    }
    void CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<Complex> elmat,
			    LocalHeap & lh) const {
      T_CalcElementMatrix<Complex>(base_fel,eltrans,elmat, lh);    
    }
  };



  template<int D> template <class SCAL>
  void TraceTraceXn<D>::T_CalcElementMatrix (const FiniteElement & base_fel,
					  const ElementTransformation & eltrans, 
					  FlatMatrix<SCAL> elmat,
					  LocalHeap & lh) const {
    
    const CompoundFiniteElement &  cfel  // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);

    const HCurlFiniteElement<D>   & fel_h =  // H space
      dynamic_cast<const HCurlFiniteElement<D>&  > (cfel[GetInd1()]);
    const HCurlFiniteElement<D> & fel_f =  // F space
      dynamic_cast<const HCurlFiniteElement<D>&> (cfel[GetInd2()]);

    elmat = SCAL(0.0);

    IntRange rh = cfel.GetRange(GetInd1()); 
    IntRange rf = cfel.GetRange(GetInd2());
    int ndofh = rh.Size();
    int ndoff = rf.Size();

    FlatMatrix<SCAL> submat(ndoff,ndofh,lh);
    FlatMatrixFixWidth<D> shapeh(ndofh,lh);  // H-basis (vec) values
    FlatMatrixFixWidth<D> shapef(ndoff,lh);  // F-basis (vec) values 
    FlatMatrixFixWidth<D> cp(ndoff,lh); // Cross products F x n

    ELEMENT_TYPE eltype                   // get the type of element: 
      = fel_h.ElementType();              // ET_TET in 3d.

    // transform facet integration points to volume integration points
    Facet2ElementTrafo transform(eltype);

    int nfa = ElementTopology::GetNFacets(eltype); /* nfa = number of 
                                                      facets of an
                                                      element */    
    submat = SCAL(0.0);

    for(int k = 0; k<nfa; k++) {
      
      // type of geometry of k-th facet
      ELEMENT_TYPE eltype_facet = ElementTopology::GetFacetType(eltype, k); 
      
      const IntegrationRule & facet_ir =
	SelectIntegrationRule (eltype_facet, fel_h.Order()+fel_f.Order()); 

      // reference element normal vector
      FlatVec<D> normal_ref = ElementTopology::GetNormals(eltype) [k]; 

      for (int l = 0; l < facet_ir.GetNIP(); l++) {

	// map facet points to volume integration points
	IntegrationPoint volume_ip = transform(k, facet_ir[l]);
	MappedIntegrationPoint<D,D> mip (volume_ip, eltrans);
	
	// compute normal on physcial element
	Mat<D> inv_jac = mip.GetJacobianInverse();
	double det = mip.GetJacobiDet();
	Vec<D> normal = fabs(det) * Trans(inv_jac) * normal_ref;       
	double len = L2Norm(normal);
	normal /= len;
	double weight = facet_ir[l].Weight()*len;
	
	// mapped H(curl) basis fn values 
 	fel_h.CalcMappedShape(mip,shapeh); 
	fel_f.CalcMappedShape(mip,shapef); 

	// F x n
	cp.Col(0) = normal(2)*shapef.Col(1)-normal(1)*shapef.Col(2);
	cp.Col(1) = normal(0)*shapef.Col(2)-normal(2)*shapef.Col(0);
	cp.Col(2) = normal(1)*shapef.Col(0)-normal(0)*shapef.Col(1);

	// evaluate coefficient
	SCAL dd = coeff_d -> T_Evaluate<SCAL>(mip);

	//                  [ndoff x D]    [ndofh x D] 	
	submat += (dd*weight) * cp * Trans( shapeh ) ;
      }
    }

    elmat.Rows(rf).Cols(rh) += submat;

    if (GetInd1() != GetInd2())
      elmat.Rows(rh).Cols(rf) += Conj(Trans(submat));
  }


  //////////////////////////////////////////////////////////////
  // Integrate a(x)* u . e, where u and e are in H(curl) spaces

  template<int D> class EyeEyeEdge : public DPGintegrator  {
    
    shared_ptr<CoefficientFunction> coeff_a;
    
    template<class SCAL>
    void T_CalcElementMatrix (const FiniteElement & base_fel,
			      const ElementTransformation & eltrans, 
			      FlatMatrix<SCAL> elmat,
			      LocalHeap & lh)  const ;

  public:
    
    EyeEyeEdge(const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : DPGintegrator(coeffs), coeff_a(coeffs[2])  {
    
      cout << "Using DPG integrator " << Name() << " with components "
	   << GetInd1()+1 << " and " << GetInd2()+1 << endl ; 
    }

    virtual bool IsSymmetric() const { return !coeff_a->IsComplex() ; }
    
    virtual string Name () const { return "EyeEyeEdge"; }

    virtual bool BoundaryForm () const { return false; }

    void CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<double> elmat,
			    LocalHeap & lh) const {
      T_CalcElementMatrix<double>(base_fel,eltrans,elmat,lh);
					       
    }
    void CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<Complex> elmat,
			    LocalHeap & lh) const {
      T_CalcElementMatrix<Complex>(base_fel,eltrans,elmat, lh);    
    }
  };


  template<int D> template <class SCAL>
  void EyeEyeEdge<D>::T_CalcElementMatrix (const FiniteElement & base_fel,
  				       const ElementTransformation & eltrans, 
  				       FlatMatrix<SCAL> elmat,
  				       LocalHeap & lh) const {

    const CompoundFiniteElement &  cfel  // product space 
      =  dynamic_cast<const CompoundFiniteElement&> (base_fel);
    const HCurlFiniteElement<D> & fel_u =  // u space
      dynamic_cast<const HCurlFiniteElement<D>&> (cfel[GetInd1()]);
    const HCurlFiniteElement<D> & fel_e =  // e space
      dynamic_cast<const HCurlFiniteElement<D>&> (cfel[GetInd2()]);

    elmat = SCAL(0.0);

    IntRange ru = cfel.GetRange(GetInd1()); 
    IntRange re = cfel.GetRange(GetInd2()); 
    int ndofu = ru.Size();
    int ndofe = re.Size();

    FlatMatrixFixWidth<D> ushape(ndofu,lh);
    FlatMatrixFixWidth<D> eshape(ndofe,lh);

    ELEMENT_TYPE eltype = fel_u.ElementType();      
    const IntegrationRule &         
      ir = SelectIntegrationRule(eltype, fel_u.Order()+fel_e.Order());
    FlatMatrix<SCAL> submat(ndofe,ndofu,lh);
    submat = SCAL(0.0);

    for(int k=0; k<ir.GetNIP(); k++) {	
      
      MappedIntegrationPoint<D,D> mip (ir[k],eltrans);

      fel_u.CalcMappedShape( mip, ushape ); 
      fel_e.CalcMappedShape( mip, eshape );     

      SCAL fac = (coeff_a -> T_Evaluate<SCAL>(mip))* mip.GetWeight() ;
      //               [ndofe x D] * [D x ndofu]
      submat +=  fac *   eshape    * Trans(ushape) ;
    }
   
    elmat.Rows(re).Cols(ru) += submat;

    if (GetInd1() != GetInd2())
      elmat.Rows(ru).Cols(re) += Conj(Trans(submat));
  }


//////////////////////////////////////////////////////////////
// Integrate c(x) * H * (W x n) over the GLOBAL BOUNDARY, 
// where H and W may be from different H(curl) spaces.

template<int D> 
class XnBoundary : public DPGintegrator   {
    
  shared_ptr<CoefficientFunction> coeff_c;
    
  template<class SCAL>
  void T_CalcElementMatrix (const FiniteElement & base_fel,
			    const ElementTransformation & eltrans, 
			    FlatMatrix<SCAL> elmat,
			    LocalHeap & lh)  const ; 
public:

  XnBoundary(const Array<shared_ptr<CoefficientFunction>> & coeffs)
    : DPGintegrator(coeffs), coeff_c(coeffs[2]) {

    cout << "Using DPG integrator " << Name() << " with components "
	 << GetInd1()+1 << " and " << GetInd2()+1 << endl ; 
  }

  virtual bool IsSymmetric() const { return !coeff_c->IsComplex() ; }
    
  virtual string Name () const { return "XnBoundary"; }
  virtual int DimElement () const { return D-1; }
  virtual int DimSpace () const { return D; }		
  virtual bool BoundaryForm () const { return 1; }

  void CalcElementMatrix (const FiniteElement & base_fel,
			  const ElementTransformation & eltrans, 
			  FlatMatrix<double> elmat,
			  LocalHeap & lh) const {
    T_CalcElementMatrix<double>(base_fel,eltrans,elmat,lh);
					       
  }
  void CalcElementMatrix (const FiniteElement & base_fel,
			  const ElementTransformation & eltrans, 
			  FlatMatrix<Complex> elmat,
			  LocalHeap & lh) const {
    T_CalcElementMatrix<Complex>(base_fel,eltrans,elmat, lh);    
  }
};


template<int D> template <class SCAL>
void XnBoundary<D> ::
T_CalcElementMatrix (const FiniteElement & base_fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<SCAL> elmat,
		     LocalHeap & lh) const {

  const CompoundFiniteElement &  cfel      // product space 
    =  dynamic_cast<const CompoundFiniteElement&> (base_fel);
    
  const HCurlFiniteElement<D-1> & fel_h = // H space
    dynamic_cast<const HCurlFiniteElement<D-1>&> (cfel[GetInd1()]);

  const HCurlFiniteElement<D-1> & fel_w = // W space
    dynamic_cast<const HCurlFiniteElement<D-1>&> (cfel[GetInd2()]);

  elmat = SCAL(0.0);

  IntRange rh = cfel.GetRange(GetInd1());
  IntRange rw = cfel.GetRange(GetInd2());
  int ndofh = rh.Size();
  int ndofw = rw.Size();
 
  FlatMatrix<SCAL> submat(ndofw, ndofh, lh);  
  submat = SCAL(0.0);
    
  FlatMatrixFixWidth<D> shapeh_ref(ndofh, lh); // H-basis on reference
                                               // surface element
  FlatMatrixFixWidth<D> shapeh(ndofh, lh);     // H-basis on physical
                                               // element
  FlatMatrixFixWidth<D> shapew_ref(ndofw, lh); // W-basis on reference
                                               // surface element
  FlatMatrixFixWidth<D> shapew(ndofw, lh);     // W-basis on physical
                                               // element
  FlatMatrixFixWidth<D> cpw(ndofw,lh);         // Cross products W x n

  const IntegrationRule ir(fel_h.ElementType(), 
			   fel_h.Order() + fel_w.Order());

  Mat<D> shape_map; 
  Vec<D> normal;

  for (int i = 0 ; i < ir.GetNIP(); i++) {

    MappedIntegrationPoint<D-1,D> mip(ir[i], eltrans);

    SCAL cc = coeff_c -> T_Evaluate<SCAL>(mip);

    fel_w.CalcShape (ir[i], shapew_ref);
    fel_h.CalcShape (ir[i], shapeh_ref);
    normal = mip.GetNV();

    for (int ii = 0; ii < (D - 1); ii++)
      for (int jj = 0; jj < D; jj++)
	shape_map(ii, jj) = mip.GetJacobian()(jj, ii);
    for (int jj = 0; jj < D; jj++)
      shape_map(D - 1, jj) = normal(jj);

    shapew = shapew_ref * Inv( Trans(shape_map) );
    shapeh = shapeh_ref * Inv( Trans(shape_map) );

    // W x n
    cpw.Col(0) = normal(2)*shapew.Col(1)-normal(1)*shapew.Col(2);
    cpw.Col(1) = normal(0)*shapew.Col(2)-normal(2)*shapew.Col(0);
    cpw.Col(2) = normal(1)*shapew.Col(0)-normal(0)*shapew.Col(1);

    //                              [ndofw x D]  [D x ndofh]
    submat += (cc*mip.GetWeight()) * cpw  * Trans(shapeh);
  }       

  elmat.Rows(rw).Cols(rh) += submat;

  if (GetInd1() != GetInd2())
    elmat.Rows(rh).Cols(rw) += Conj(Trans(submat));
}
  
  //////////////////////////////////////////////////////////////

  static RegisterBilinearFormIntegrator<CurlCurlPG<3>> 
  initcurlcurlpg("curlcurlpg", 3, 3);

  static RegisterBilinearFormIntegrator<TraceTraceXn<3>> 
  inittrctrcxn("trctrcxn", 3, 3);

  static RegisterBilinearFormIntegrator<EyeEyeEdge<2>> 
  initmassedge2d("eyeeyeedge", 2, 3);
  static RegisterBilinearFormIntegrator<EyeEyeEdge<3>> 
  initmassedge3d("eyeeyeedge", 3, 3);

  static RegisterBilinearFormIntegrator<XnBoundary<3>> 
  initxnb3d("xnbdry", 3, 3);

} // closing namespace dpg

