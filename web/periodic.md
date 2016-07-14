### Periodic mesh

Periodic identifications between surfaces can be made in geometry. Netgen uses it to generate meshes suitable for imposition of periodic boundary conditions. 

- Example: [periodicbox.py](../projects/nanogap/learn2mesh/periodicbox.py). You can ask Netgen to mark (in purple) vertices that are identified to be the same due to periodicity from the menu option `View->Mesh->Show identified points` and you should get a picture like this:

![Image of mesh](figs/mesh2periodic.png)


### Periodic finite element space 

Once you have a periodic mesh, the next step in implementing a periodic boundary condition is to implement the periodic finite element space you need. This can be done by **identifying** degrees of freedoms. Two examples are available:

- [Periodic H1 space](../spaces/periodich1.cpp)
- [Periodic H(curl) space](../spaces/periodichcurl.cpp)