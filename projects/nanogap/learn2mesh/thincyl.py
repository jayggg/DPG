from netgen.csg import *

X=50
Y=50
Z=5
r0=16
w = 1
d = 1

rinn = r0
rout = r0+w
outercyls = []
innercyls = []
outercyls.append(Cylinder( Pnt(0,0,-1), Pnt(0,0,0), rout))
innercyls.append(Cylinder( Pnt(0,0,-1), Pnt(0,0,0), rinn))

rings = []           # layer part in between outer & inner cylinders
rings.append( (outercyls[0] - innercyls[0]) *
              OrthoBrick(Pnt(-X/2,-Y/2,-Z/2),Pnt(X/2,Y/2,Z/2)) )

geo = CSGeometry()
geo.Add(rings[0])

geo.CloseSurfaces(outercyls[0],innercyls[0]) # declare close cylinders


ngmesh = geo.GenerateMesh()
m = Mesh(ngmesh)
m.Curve(6)
Draw(m)
