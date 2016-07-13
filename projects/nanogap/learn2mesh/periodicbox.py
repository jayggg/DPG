from netgen.csg import *

X=50
Y=50
Z=200
r0=16

# periodic enclosure

xneg = Plane(Pnt(-X/2,  0,   0), Vec(-1,  0, 0)).bc("x-")
xpos = Plane(Pnt( X/2,  0,   0), Vec( 1,  0, 0)).bc("x+")
yneg = Plane(Pnt(  0, -Y/2,  0), Vec( 0, -1, 0)).bc("y-")
ypos = Plane(Pnt(  0,  Y/2,  0), Vec( 0,  1, 0)).bc("y+")
enclperiodic = xneg * xpos * yneg * ypos


airbelow = OrthoBrick(Pnt(-X/2,-Y/2,-Z/2),Pnt(X/2,Y/2,0)) * enclperiodic
airbelow.mat('airbot').bc('airbelow').maxh(20)

airabove = OrthoBrick(Pnt(-X/2,-Y/2,0),Pnt(X/2,Y/2,Z/2)) * enclperiodic
airabove.mat('airtop').bc('airabove').maxh(20)                       

geo = CSGeometry()

geo.Add(airabove)
geo.Add(airbelow)

geo.PeriodicSurfaces(xneg, xpos)      # declare x periodicity 
geo.PeriodicSurfaces(yneg, ypos)      # declare y periodicity 


ngmesh = geo.GenerateMesh()
m = Mesh(ngmesh)
Draw(m)
