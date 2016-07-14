### Periodic mesh

Periodic identifications between surfaces can be made in geometry. Netgen uses it to generate meshes suitable for imposition of periodic boundary conditions. 

- Example: [periodicbox.py](../projects/nanogap/learn2mesh/periodicbox.py). You can ask Netgen to mark vertices that are identified to be the same due to periodicity from the menu option `View->Mesh->Show identified points` and you should see a mesh  with purple markings like this:

![Image of mesh](figs/mesh2periodic.png)


### Periodic finite element space 

Once you have a periodic mesh, the next step in implementing a periodic boundary condition is to implement the periodic finite element space you need. This can be done by **identifying** degrees of freedoms. An example with one pair of x-periodic surfaces is available  within the  NGSolve developer's tutorial [MyLittleNGSolve](https://sourceforge.net/p/ngsolve/ml-ngs/ci/master/tree/periodic.cpp). Here, we provide  two additional xy-periodic spaces  in `libDPG`: 

- [Periodic H1 space](../spaces/periodich1.cpp)
- [Periodic H(curl) space](../spaces/periodichcurl.cpp)

Usage examples: 

- [magnet.pde](../pde/magnet.pde) - uses periodic space and standard FEM (not DPG). 
- [periodicmaxwell.py](../python/periodicmaxwell.py) - uses periodic space and DPG forms. 
