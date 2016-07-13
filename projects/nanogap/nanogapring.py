from ngsolve import *
from netgen.csg import *
from ctypes import CDLL
from cmath import pi, sqrt

import sys
sys.path.append('../pyutils')
from pcg import pcg

    
def ringgeom(w, nlayers,
             ncyl = 1,
             fine = False,
             thickness={'alox'  : 0.010,   # alox layer thickness
                        'gold'  : 0.150,   # gold layer thickness
                        'glass' : 0.500 }, # glass layer thickness
             X=50, Y=50, Z=200, r0=16 ):

    """ Create and return Netgen geometry. 

    INPUT GEOMETRY PARAMETERS  in length units of micrometers:

    w:  Annular nanogap width
    r0: Inner radius of annular gap
    nlayers: Each material slab is split into number of layers given
    thickness: Thickness of each material slab
    X:  Length of the enclosure along x-axis
    Y:  Length of the enclosure along y-axis
    Z:  Length of the enclosure along z-axis
    fine: If false, make a coarse mesh """

    # split materials in the interesting region to layers
    
    r1 = r0 + w                 # outer  radius of annular gap
    thickness['airtop'] = max(10*w,0.5)  # add top air layer
    thickness['airbot'] = max(10*w,0.5)  # add air layer below

    assert nlayers['gold'], 'There must be a gold layer!'
    hlayer={                    # thickness of layers in each material 
        'gold' : thickness['gold'] / nlayers['gold'] }
    for key in nlayers:         # revise material slab thickness 
        if nlayers[key] == 0:   # to 0 if user specified numlayers=0   
            thickness[key] = 0.0
        else:
            hlayer[key] = thickness[key] / nlayers[key]     
    smallest_airbot_layer_h = hlayer['gold']
    smallest_airtop_layer_h = 2*hlayer['gold']

    zs={'airtop': [ 0 ],        # preparing to set layer z-coordinates
        'alox'  : [ 0 ],        # alox-air interface is z=0
        'gold'  : [-thickness['alox']],  
        'glass' : [-thickness['alox']-thickness['gold']],
        'airbot': [-thickness['alox']-thickness['gold']-thickness['glass']] }

    for key in hlayer:          # loop over alox, gold, glass to set z
        for i in range(nlayers[key]):
            zs[key].append( zs[key][-1] - hlayer[key] )
            
    zs['airtop'].insert(0, smallest_airtop_layer_h)
    while zs['airtop'][0] < thickness['airtop']:
        # always add air layers of doubling thickness 
        zs['airtop'].insert(0, 2*zs['airtop'][0])
    depth0 = -thickness['alox']-thickness['gold']-thickness['glass']
    i = 0 
    while  depth0 - zs['airbot'][-1] < thickness['airbot']:
        i = i+1
        zs['airbot'].append(  zs['airbot'][-1] - smallest_airbot_layer_h * 2**i)
        
    zs['airtop'][-1] = zs['alox'][0]  # a layer ends where another begins
    zs['alox']  [-1] = zs['gold'][0]
    zs['gold']  [-1] = zs['glass'][0]
    zs['glass'] [-1] = zs['airbot'][0]

    print('\n\nMesh layer z-coordinates: zs = ', zs)
    
    # make the ring with nanogap

    origin = Pnt(0,0,0)            # alox-air interface is z=0
    outercyls = []
    innercyls = []
    rinn = r0
    rout = r1
    outercyls.append(Cylinder( Pnt(0,0,-1), Pnt(0,0,0), rout))
    innercyls.append(Cylinder( Pnt(0,0,-1), Pnt(0,0,0), rinn))
    
    for i in range(1,ncyl):
        rout += w * 2**i
        rinn -= w * 2**i 
        outercyls.append(Cylinder( Pnt(0,0,-1), Pnt(0,0,0), rout))
        innercyls.append(Cylinder( Pnt(0,0,-1), Pnt(0,0,0), rinn))
            
    # periodic enclosure
    
    xneg = Plane(Pnt(-X/2,  0,   0), Vec(-1,  0, 0)).bc("x-")
    xpos = Plane(Pnt( X/2,  0,   0), Vec( 1,  0, 0)).bc("x+")
    yneg = Plane(Pnt(  0, -Y/2,  0), Vec( 0, -1, 0)).bc("y-")
    ypos = Plane(Pnt(  0,  Y/2,  0), Vec( 0,  1, 0)).bc("y+")
    enclperiodic = xneg * xpos * yneg * ypos

    # objects we will now make:
    
    olayers = []         # layer parts outside the outer cylinder
    ilayers = []         # layer parts inside the inside cylinder
    rings = []           # layer part in between outer & inner cylinders
    orings = []
    irings = []
    layers = []
    halfspaces = []      # layers are in between halfspaces
    
    def add_material(nl, Hl, materialname):
        for i in range(nl+1,len(Hl)+nl):
            print('Added plane at ', Hl[i-nl], ' making layer ', i)
            halfspaces.append( Plane(Pnt(0,0,Hl[i-nl]), Vec(0,0,1)) )
            
            rings.append( (outercyls[0] - innercyls[0]) * halfspaces[i-1]
                          - halfspaces[i] )                          
            if materialname == 'gold' : 
                rings[-1].mat('alox')
            else:       # only gold layer is embedded with AlOx
                rings[-1].mat(materialname)
                          
            for j in range(1,ncyl):
                orings.append( (outercyls[j]-outercyls[j-1]) * halfspaces[i-1]
                               - halfspaces[i] )
                irings.append( (innercyls[j-1]-innercyls[j]) * halfspaces[i-1]
                               - halfspaces[i] )
                orings[-1].mat(materialname)
                irings[-1].mat(materialname)    

                
            layers.append ( halfspaces[i-1] - halfspaces[i] )
            olayers.append( (layers[i-1] - outercyls[-1]) * enclperiodic )
            ilayers.append( layers[i-1] * innercyls[-1] )            
                
            olayers[i-1].mat(materialname)
            ilayers[i-1].mat(materialname)        

            
    Hl = zs['airtop']   # start with first airtop plane & add layers below
    nl = 0              # running number of layers
    halfspaces.append( Plane(Pnt(0,0,Hl[0]), Vec(0,0,1)) )
    print('Added plane at ', Hl[0])

    add_material(nl, Hl, 'airtop')
    nl += len(Hl) - 1 

    if nlayers['alox'] > 0: 
        Hl = zs['alox']
        add_material(nl, Hl, 'alox')
        nl += len(Hl) - 1 

    Hl = zs['gold']
    add_material(nl, Hl, 'gold')
    nl += len(Hl) - 1 

    if nlayers['glass'] > 0: 
        Hl = zs['glass']
        add_material(nl, Hl, 'glass')
        nl += len(Hl) - 1 

    Hl = zs['airbot']
    add_material(nl, Hl, 'airbot')
    nl += len(Hl) - 1 

    # fill with more air below
    
    d = zs['airbot'][-1]
    airbelow = OrthoBrick(Pnt(-X/2,-Y/2,-Z/2),Pnt(X/2,Y/2,d)) * enclperiodic
    airbelow.mat('airbot').bc('airbelow').maxh(20)

    # fill with more air above
    
    d = zs['airtop'][0]
    airabove = OrthoBrick(Pnt(-X/2,-Y/2,d),Pnt(X/2,Y/2,Z/2)) * enclperiodic
    airabove.mat('airtop').bc('airabove').maxh(20)                       

    # add objects to make the total  geometry
    
    geo = CSGeometry()

    for i in range(len(rings)):
        geo.Add(rings[i])

    for i in range(len(olayers)):
        
        geo.Add(olayers[i])
        geo.Add(ilayers[i])
        
        # alert close planar layers to mesher: 
        geo.CloseSurfaces(halfspaces[i], halfspaces[i+1])

    
    for j in range(len(orings)):
        
        geo.Add(orings[j])
        geo.Add(irings[j])

    geo.CloseSurfaces(outercyls[0],innercyls[0]) # declare close cylinders
    
    for j in range(1,ncyl):

        geo.CloseSurfaces(outercyls[j], outercyls[j-1])
        geo.CloseSurfaces(innercyls[j-1], innercyls[j])

    geo.Add(airabove)
    geo.Add(airbelow)

    geo.PeriodicSurfaces(xneg, xpos)      # declare x periodicity 
    geo.PeriodicSurfaces(yneg, ypos)      # declare y periodicity 

    return geo



def genmesh(w, nlayers, ncyl, savemeshfile='',
             X=50, Y=50, Z=200 ):

    geom = ringgeom(w, nlayers, ncyl, X=X, Y=Y, Z=Z)
    ngmesh = geom.GenerateMesh()
    if len(savemeshfile):
        ngmesh.Save(savemeshfile)
    m = Mesh(ngmesh)
    return m



def solve(meshfile,
          p=1,
          freq=0.625e12,
          localprec=False,
          cgiterations=10000,
          dpglib='../../libDPG.so',          
          X=50, Y=50, Z=200 ):
    """
    Solve using the DPG method. INPUTS: 
    
    p :    polynomial degree 
    freq:  incident wave frequency 
    localprec: If true, use local preconditioner, else use direct solve
    
    """
    
    # load DPG C++ lib  & load (or make) mesh 
    
    libDPG = CDLL(dpglib)

    mesh = Mesh(meshfile)
    
    Draw(mesh)
    mesh.Curve(max(3,p))

    # material properties 

    omega= 2 * pi * freq
    mu0  = pi * 4e-7             # kg * m * s^(-2) * A^(-2)
    ep0  = 8.854e-12             # A^2 * s^4 * kg^(-1) * m^(-3) 
    mu0s = mu0 * 1e6             # scaled mu0 in micrometers 
    ep0s = ep0 * 1e-18           # scaled ep0 in micrometers 
    omp =  1.37e16               # s^(-1), used in Druid model
    gam =  40.7e12               # s^(-1)
    k0 = omega * sqrt(ep0s*mu0s) # final scaled wavenumber used for air
    materialk = \
    {'air'   : k0,     'airtop': k0,     'airbot': k0,
     'alox'  : k0 * sqrt(2.34),
     'gold'  : k0 * sqrt( 1-omp*omp/(omega*omega+gam*gam) +
                          1j*gam/(omega*(omega*omega+gam*gam)) ),
     'glass' : k0 * sqrt(1.95)  }
    
    klist = [ materialk[mat] for mat in mesh.GetMaterials() ]
    k = CoefficientFunction( klist )

    # boundary parts
    
    bcs = mesh.GetBoundaries()
    bdrylist = [0] * len(bcs)
    airabovebcs = [i for i, bc in enumerate(bcs) if bc=='airabove']
    airbelowbcs = [i for i, bc in enumerate(bcs) if bc=='airbelow']
    for index in airabovebcs:  
        bdrylist[index] = 1   # indicator fn of top/bot air boundaries
    for index in airbelowbcs:
        bdrylist[index] = 1
        
    bdry  = CoefficientFunction( bdrylist )
    kbdry = materialk['air'] * bdry

    # spaces
    
    S0 = FESpace("hcurlho", mesh, order=p+2, complex=True,
                 flags={"discontinuous":True})
    S1 = FESpace("hcurlho_periodic", mesh, order=p, complex=True,
                 flags={'xends':[-X/2,X/2], 'yends':[-Y/2,Y/2]})
    S2 = FESpace("hcurlho_periodic", mesh, order=p+1, complex=True,
                 flags={"orderinner": 0,
                        'xends':[-X/2,X/2], 'yends':[-Y/2,Y/2]})
    S = FESpace( [S0,S1,S2], flags={"complex":True})

    e,E,M = S.TrialFunction()
    v,F,W = S.TestFunction()

    Einc = GridFunction(S1, 'Incident')
    einc = CoefficientFunction( (exp(1j * k0 * z), 0,0) )
    with TaskManager():
        Einc.Set(einc)
    Draw(Einc)
    k0 = CoefficientFunction(materialk['air'])
    
    f = (k*k - k0*k0) * einc     # this is the rhs for the pde

    # forms
    
    b = LinearForm(S)
    b+= SymbolicLFI(f * v)
    b.Assemble()

    a = BilinearForm(S, symmetric=False, flags={"eliminate_internal" : True})
    a+= BFI("curlcurlpg", coef=[2,1,1])       # (curl E, curl v)
    a+= BFI("eyeeyeedge", coef=[2,1,-k*k])    # -(k*k E, v)
    a+= BFI("trctrcxn",   coef=[3,1,1j] )     # i<<M, v x n>>
    a+= BFI("xnbdry", coef=[2,3,kbdry])       # <k E, W x n>
    a.components[1]+= BFI("robinedge",        # -<k*kbar E x n, F x n>
                          coef=-kbdry * Conj(kbdry))                        
    a.components[2]+= BFI("robinedge",
                          coef=-bdry)         # -<M x n, W x n>
    a.components[0]+= BFI("massedge",         
                          coef=1.0)           # (e, v)
    a.components[0]+= BFI("curlcurledge",     # (curl e, curl v)
                          coef=1.0)        

    # assemble
    
    eEM = GridFunction(S, 'Scattered')
    
    if localprec:
        c = Preconditioner(a, type="local")
    else:
        c = Preconditioner(a, type="direct")

    with TaskManager():    
        a.Assemble(heapsize=int(5e8))

    iterates_to_save = [i*cgiterations//10 for i in range(1,11)]
    def save_pcg_iterate(x, iter):
        if iter in iterates_to_save:
            solfilename='solpcg%d'%iter+'.sol'
            eEM.Save(solfilename)
    
    b.vec.data += a.harmonic_extension_trans * b.vec
    eEM.vec[:] = 0.0

    # solve
    
    with TaskManager():    
        eEM.vec.data = pcg(a.mat, c.mat, b.vec, x=eEM.vec,
                           maxits=cgiterations,
                           saveitfn=save_pcg_iterate)
        
        eEM.vec.data += a.harmonic_extension * eEM.vec
        eEM.vec.data += a.inner_solve * b.vec
        
    Esct = eEM.components[1]
    Draw(Esct)

    Etot = GridFunction(S1, 'Total')
    Etot.vec.data = Einc.vec + Esct.vec
    Draw(Etot)

    return Etot

def loadsol(solfileEtot, meshfile, p,
            dpglib='../../libDPG.so',
            X=50, Y=50, Z=200):

    libDPG = CDLL('../../libDPG.so')
    mesh = Mesh(meshfile)
    mesh.Curve(max(3,p))
    Draw(mesh)

    S0 = FESpace("hcurlho", mesh, order=p+2, complex=True,
                 flags={"discontinuous":True})
    S1 = FESpace("hcurlho_periodic", mesh, order=p, complex=True,
                 flags={'xends':[-X/2,X/2], 'yends':[-Y/2,Y/2]})
    S2 = FESpace("hcurlho_periodic", mesh, order=p+1, complex=True,
                 flags={"orderinner": 0,
                        'xends':[-X/2,X/2], 'yends':[-Y/2,Y/2]})
    S = FESpace( [S0,S1,S2], flags={"complex":True})


    Etot = GridFunction(S1, 'Total')
    Etot.Load(solfileEtot)
    Draw(Etot)
    
    return Etot
    
 
