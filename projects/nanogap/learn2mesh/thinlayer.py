from netgen.csg import *

X=50
Y=50
Z=200
r0=16
d = 1

# periodic enclosure

xneg = Plane(Pnt(-X/2,  0,   0), Vec(-1,  0, 0)).bc("x-")
xpos = Plane(Pnt( X/2,  0,   0), Vec( 1,  0, 0)).bc("x+")
yneg = Plane(Pnt(  0, -Y/2,  0), Vec( 0, -1, 0)).bc("y-")
ypos = Plane(Pnt(  0,  Y/2,  0), Vec( 0,  1, 0)).bc("y+")
enclperiodic = xneg * xpos * yneg * ypos


halfspaces = []      # layers are in between halfspaces
halfspaces.append( Plane(Pnt(0,0, 0), Vec(0,0,1)) )
halfspaces.append( Plane(Pnt(0,0,-d), Vec(0,0,1)) )
layer =  (halfspaces[0] - halfspaces[1]) * enclperiodic
#halfspaces[0].maxh(3)
            
airabove = OrthoBrick(Pnt(-X/2,-Y/2,0),Pnt(X/2,Y/2,Z/2)) * enclperiodic
airabove.mat('airtop').bc('airabove').maxh(20)                       

airbelow = OrthoBrick(Pnt(-X/2,-Y/2,-Z/2),Pnt(X/2,Y/2,-d)) * enclperiodic
airbelow.mat('airbot').bc('airbelow').maxh(20)


geo = CSGeometry()

geo.Add(layer)
geo.Add(airabove)
geo.Add(airbelow)

geo.CloseSurfaces(halfspaces[0], halfspaces[1])

geo.PeriodicSurfaces(xneg, xpos)      # declare x periodicity 
geo.PeriodicSurfaces(yneg, ypos)      # declare y periodicity 


ngmesh = geo.GenerateMesh()
m = Mesh(ngmesh)
Draw(m)
