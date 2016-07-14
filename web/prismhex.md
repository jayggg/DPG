### Meshes with Prismatic elements

Netgen can be alerted to close surfaces which will automatically trigger prismatic element generation. Here are two examples: 

- [thincyl.py](../projects/nanogap/learn2mesh/thincyl.py).
- [thinlayer.py](../projects/nanogap/learn2mesh/thinlayer.py).

The first example  shows how thin annular gaps can be meshed by curved prisms.

![Image of mesh](figs/mesh6cyl.png)

The second example gives prisms, which can be seen after zooming into the thin layer in the middle:

![Image of mesh](figs/mesh5prism.png)


### Meshes with Hexahedral elements

Sometimes we have a situation where two pairs of close surfaces intersect, such as in this figure:

![Image of mesh](figs/mesh7cylayer.png)

In such cases, Netgen will attempt hexahedral element generation. Here is an example: 

-  [thinlayercyl.py](../projects/nanogap/learn2mesh/thinlayercyl.py).

The generated mesh, when zoomed in near the niddle thin layer shows a ring of hexahderal elements:

![Image of mesh](figs/mesh8hexring.png)