/*

  Provides Schwarz subspace correction preconditioners using vertex
  patches.

  All dofs, after condensation, that are on facets connected to a
  vertex, define the subspaces. Corrections on these subspaces are
  multiplicatively combined. 

  An optional coarse solve using wirebasket dofs is turned off by
  default.

*/


#include <solve.hpp>


namespace ngcomp  {
  
  class VertexPatchSchwarz : public Preconditioner  {

    shared_ptr<BilinearForm> bfa;
    shared_ptr<BaseBlockJacobiPrecond> jacobi;
    shared_ptr<BaseMatrix>  coarseinv;
    bool                    addcoarse; 

  public:

    VertexPatchSchwarz (const PDE & pde, const Flags & flags, 
			const string & aname);
    VertexPatchSchwarz (shared_ptr<BilinearForm> abfa, const Flags & aflags,
                        const string aname = "mgprecond");

    ~VertexPatchSchwarz ();

    virtual void Update();

    virtual int VHeight() const { return jacobi->VHeight(); }
    
    virtual int VWidth() const { return jacobi->VWidth(); }
	
    virtual void Mult (const BaseVector & f, BaseVector & u) const  {

      // jacobi -> Mult (f, u);

      u = 0.0;
      jacobi -> GSSmooth (u, f);
      jacobi -> GSSmoothBack (u, f);

      if (addcoarse)
	coarseinv->MultAdd( 1, f, u );  // u = u + 1 * inv(A0) * f
    }

    virtual const BaseMatrix & GetAMatrix() const   {

      return bfa -> GetMatrix();
    }

  };


  VertexPatchSchwarz::VertexPatchSchwarz (const PDE & pde, 
		      const Flags & flags, const string & aname)
    : Preconditioner (&pde, flags, aname), jacobi(NULL)  {

    addcoarse = flags.GetDefineFlag("addcoarse");

    cout << endl << "Constructor of VertexPatchSchwarz" ;
    if (addcoarse) cout << "with coarse solve" ;      
	
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
  }
    
  VertexPatchSchwarz::
  VertexPatchSchwarz (shared_ptr<BilinearForm> abfa, 
                      const Flags & aflags, const string aname)
    : Preconditioner (abfa, aflags, aname)
  {
    addcoarse = flags.GetDefineFlag("addcoarse");

    cout << endl << "Constructor of VertexPatchSchwarz" ;
    if (addcoarse) cout << "with coarse solve" ;      
	
    bfa = abfa;
  }


  VertexPatchSchwarz :: ~VertexPatchSchwarz ()  
  { 
    // noting to delete with shared ptr
    // delete jacobi; 
  }

  
  void VertexPatchSchwarz :: Update()   {

    // delete jacobi;

    const BaseSparseMatrix & mat 
      = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix());
    const BitArray * freedofs = bfa->GetFESpace()->GetFreeDofs();
    shared_ptr<FESpace> fes = bfa -> GetFESpace();
    

    BitArray used (fes->GetNDof());
    used.Clear();
    for (int i = 0; i < ma->GetNE(); i++)
      {
        Array<int> dofs;
        fes->GetDofNrs (i, dofs);
        for (int j = 0; j < dofs.Size(); j++)
          used.Set (dofs[j]);
      }
    for (int i = 0; i < fes->GetNDof(); i++)
      if (freedofs->Test(i) && !used.Test(i))
        cerr << "freedof, but never used: " << i << endl;
        
    

    /**** First attempt at making the required sparse matrices: */

    // int nv = ma->GetNV();

    // // cnt[p] = number of all dofs connected to vertex number p
    // Array<int> cnt(nv);
    // cnt = 1;   // vertex counted

    // for (int i=0; i < ma->GetNEdges(); i++) {
      
    //   Array<int> dofs, pnums;
    //   ma->GetEdgePNums(i,pnums);  // pnums = vertices of the edge
    //   fes->GetEdgeDofNrs(i,dofs);   // 
      
    //   for (int j=0; j<2; j++) 
	
    // 	cnt[pnums[j]] += dofs.Size();

     
    // }  
      
    // Table<int> & blocks = *( new Table<int>(cnt) );

    // for (int p=0; p<nv; p++ ) 
    //   blocks[p][0] = p;

    // cnt = 1;

    // for (int i=0; i < ma->GetNEdges(); i++) {
      
    //   Array<int> dofs, pnums;
    //   ma->GetEdgePNums(i,pnums);  // pnums = vertices of the edge
    //   fes->GetEdgeDofNrs(i,dofs);   // 
      
    //   for (int j=0; j<2; j++) 
    // 	for (int l=0; l<dofs.Size(); l++) 

    // 	  blocks[pnums[j]][cnt[pnums[j]]++] = dofs[l];

    // }
    // cout << "Blocks: "<< endl << blocks << endl;



    /****  Second attempt (better): Uses the NGSolve table creator  */
    /*
        A table-creator creates a table in compressed form.  
	A filtered table-creator also filters out Dirichlet-dofs and
	eliminated internal dofs.
     */

    FilteredTableCreator creator(freedofs);

    while (!creator.Done())  {
      
      Array<int> vdofs, dofs, pnums;

      for (int i=0; i<ma->GetNV(); i++) {
	
	// vdofs = dof# of vertex i (many for compound space)
	fes->GetVertexDofNrs(i,vdofs);

	// to the block of vertex i, add that vertex dof
	creator.Add(i, vdofs);  // Add(block#, array of dof#s)
      }

      for (int i=0; i < ma->GetNEdges(); i++) {
      	
	// pnums = endpoint vertices of edge i
	  ma->GetEdgePNums(i,pnums); 
	  
	  // collect dof nums interior to edge i in dofs.
	  fes->GetEdgeDofNrs(i,dofs);
	
	for (int j=0; j<2; j++) 
	  creator.Add(pnums[j],dofs);
      }
      
      if (ma->GetDimension()==3) {   // 3D case 

	for (int i=0; i < ma->GetNFaces(); i++) {      

	  // pnums = vertices of face i
	  ma->GetFacePNums(i,pnums); 
	  
	  // collect dof nums interior to face i in dofs.
	  fes->GetFaceDofNrs(i,dofs);
	  
	  for (int j=0; j<pnums.Size(); j++) 
	    creator.Add(pnums[j],dofs);
	}
      }
      
      creator++;
    } 
    
    //cout << "Blocks: "<< endl << *creator.GetTable() << endl;
    
    jacobi = mat.CreateBlockJacobiPrecond (*creator.GetTable());

    if (addcoarse) {
      int ndof = fes->GetNDof();
      BitArray * coarsedofs = new BitArray(ndof);
      coarsedofs->Clear();

      for (int i=0; i<ndof; i++) 
	if (fes->GetDofCouplingType(i) == WIREBASKET_DOF)
	  coarsedofs->Set(i);
      coarsedofs->And(*fes->GetFreeDofs());

      coarseinv = mat.InverseMatrix(coarsedofs);
    }
    
    if (test) Test();
  }
  

  static RegisterPreconditioner<VertexPatchSchwarz> 
  initschwarz ("vertexschwarz");

}

