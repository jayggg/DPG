"""Automatic adaptivity for spacetime wave equation using DPG method."""


import ngsolve as ngs
from ngsolve import VTKOutput, sqrt, Mesh, exp, x, GridFunction
from netgen.geom2d import unit_square
from wave import vec, waveA, makeforms


# PARAMETERS:

p = 3          # polynomial degree
h0 = 1         # coarse mesh size for unit square domain
markprm = 0.5  # percentage of max total error for marking

# SET UP:

mesh = Mesh(unit_square.GenerateMesh(maxh=h0))
q_zero = 'bottom'              # Mesh boundary parts where q and
mu_zero = 'bottom|right|left'  # mu has essential b.c
u00 = exp(-1000 * ((x - 0.5) * (x - 0.5)))  # Nonzero initial condition
F = ngs.CoefficientFunction((0, 0))         # Zero source

a, f, X, sep = makeforms(mesh, p, F, q_zero, mu_zero, epsil=1.e-10)

euz = GridFunction(X)           # Contains solution at each adaptive step
q = euz.components[sep[0]]      # Volume (L2) components
mu = euz.components[sep[0]+1]
zq = euz.components[sep[1]]     # Interface components
zmu = euz.components[sep[1]+1]

zq.Set(u00, definedon='bottom')
zmu.Set(-u00, definedon='bottom')

ngs.Draw(mu, autoscale=False,   # draw only one of the solution components
         min=-1.0, max=1.0)


# ADAPTIVE LOOP:

globalerr = 1      # estimated total error
itcount = 0        # iteration count


def solve_on_current_mesh():

    # assemble the problem on current mesh:
    X.Update()
    euz.Update()
    euz.vec[:] = 0
    zq.Set(u00, definedon='bottom')
    zmu.Set(-u00, definedon='bottom')
    a.Assemble()
    f.Assemble()

    # solve the problem on the current mesh:
    r = f.vec.CreateVector()
    r.data = f.vec - a.mat * euz.vec
    r.data += a.harmonic_extension_trans * r
    euz.vec.data += a.mat.Inverse(X.FreeDofs(True)) * r
    euz.vec.data += a.harmonic_extension * euz.vec
    euz.vec.data += a.inner_solve * r

    # save solution for display later
    meshfilename = 'outputs/pressuremesh' + str(int(itcount))
    solfilename = 'outputs/pressure' + str(int(itcount))
    VTKOutput(ma=mesh, coefs=[mu],
              names=["mesh"], filename=meshfilename,
              subdivision=0).Do()
    VTKOutput(ma=mesh, coefs=[mu, mu],
              names=["pressure"], filename=solfilename,
              subdivision=p).Do()


def mark_for_refinement():

    # extract estimator component from solution:
    eh = [euz.components[i] for i in range(sep[0])]
    elerr = ngs.Integrate(vec(eh)*vec(eh) + waveA(eh)*waveA(eh),
                          mesh, ngs.VOL, element_wise=True)
    maxerr = max(abs(elerr.NumPy()))

    # mark elements
    for el in mesh.Elements():
        mesh.SetRefinementFlag(el, bool(abs(elerr[el.nr]) > markprm * maxerr))
    globalerr = sqrt(abs(sum(elerr)))
    print("Adaptive step %d: Estimated error=%g, Ndofs=%d"
          % (itcount, globalerr, X.ndof))


while X.ndof < 40000 and globalerr > 1.e-4:

    itcount += 1
    mesh.Refine()
    solve_on_current_mesh()
    ngs.Redraw(blocking=True)
    mark_for_refinement()
