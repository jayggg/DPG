"""
Spacetime wave simulation using DPG in 2 and 3 dimensional spacetime.

Test norm: (e,w)_W = (waveA(e, cwave), waveA(w, cwave)) + (e,w)

DPG form: b((u,z), w) = -(u, waveA(w, cwave)) + <waveD(n, z, cwave, w>

Here waveA is the first order wave operator and waveD is the
corresponding boundary operator, both defined below as python
functions.
"""

from ngsolve import CoefficientFunction, L2, H1, FESpace, specialcf
from ngsolve import BilinearForm, SymbolicBFI, GridFunction, Preconditioner
from ngsolve import LinearForm, BVP, SymbolicLFI, sin, cos
from netgen.geom2d import unit_square
from netgen.csg import unit_cube
import ngsolve as ngs
from math import log, pi
from generate_cubic_mesh import GenerateCubeMesh


def vec(listfun):
    return CoefficientFunction(tuple(f for f in listfun))


def waveA(u,cwave):  # Return action of the first order wave operator on u

    if len(u) == 2:
        q, mu = u
        dxq, dtq = ngs.grad(q)
        dxmu, dtmu = ngs.grad(mu)
        return CoefficientFunction((dtq - cwave*dxmu, dtmu - cwave*dxq))
    elif len(u) == 3:
        q1, q2, mu = u
        dxq1, dyq1, dtq1 = ngs.grad(q1)
        dxq2, dyq2, dtq2 = ngs.grad(q2)
        dxmu, dymu, dtmu = ngs.grad(mu)
        return CoefficientFunction((dtq1 - cwave*dxmu,
                                    dtq2 - cwave*dymu,
                                    dtmu - cwave*dxq1 - cwave*dyq2))
    else:
        raise ValueError('Wave operator cannot act on vectors of length %d'
                         % len(u))


def waveD(n, z, cwave):  # Return action of boundary matrix Dn on z

    if len(z) == 2:
        return CoefficientFunction((n[1]*z[0] - cwave*n[0]*z[1],
                                    n[1]*z[1] - cwave*n[0]*z[0]))
    elif len(z) == 3:
        return CoefficientFunction((n[2]*z[0] - cwave*n[0]*z[2],
                                    n[2]*z[1] - cwave*n[1]*z[2],
                                    n[2]*z[2] - cwave*n[0]*z[0] - cwave*n[1]*z[1]))
    else:
        raise ValueError('Boundary operator given vector of length %d'
                         % len(z))


def makeforms(mesh, p, F, q_zero, mu_zero, cwave, epsil=0):

    d = mesh.dim
    W = L2(mesh, order=p+d)
    U = L2(mesh, order=p)
    Zq = H1(mesh, order=p+1, dirichlet=q_zero, orderinner=0)
    Zmu = H1(mesh, order=p+1, dirichlet=mu_zero, orderinner=0)

    spacelist = [W]*d + [U]*d + [Zq]*(d-1) + [Zmu]
    X = FESpace(spacelist)
    # separate    W...W, U...U,  Zq...Zq,  Zmu
    separators = [d,     2*d,    2*d+(d-1)]

    Xtrials = X.TrialFunction()
    e = Xtrials[0:separators[0]]
    u = Xtrials[separators[0]:separators[1]]
    zu = Xtrials[separators[1]:]

    Xtests = X.TestFunction()
    w = Xtests[0:separators[0]]
    v = Xtests[separators[0]:separators[1]]
    zv = Xtests[separators[1]:]

    n = specialcf.normal(d)

    a = BilinearForm(X, symmetric=True, eliminate_internal=True)
    a += SymbolicBFI(vec(e) * vec(w))
    a += SymbolicBFI(waveA(e, cwave) * waveA(w, cwave))
    a += SymbolicBFI(-vec(u) * waveA(w, cwave))
    a += SymbolicBFI(-waveA(e, cwave) * vec(v))
    a += SymbolicBFI(vec(e) * waveD(n, zv, cwave), element_boundary=True)
    a += SymbolicBFI(waveD(n, zu, cwave) * vec(w), element_boundary=True)
    a += SymbolicBFI(-epsil * vec(zu) * vec(zv))

    f = LinearForm(X)
    f += SymbolicLFI(F*vec(w))

    return (a, f, X, separators)


def compute_error(euz, sep, exactu):

    uh = vec([euz.components[i] for i in range(sep[0], sep[1])])
    if exactu is None:
        er = None
    else:
        er = ngs.sqrt(ngs.Integrate((uh-exactu)*(uh-exactu), mesh))
        print("|| u - uh || = ", er)
    return er


def solvewave(mesh, p, F, q_zero, mu_zero, cwave, exactu=None):

    a, f, X, sep = makeforms(mesh, p, F, q_zero, mu_zero, cwave)
    euz = GridFunction(X)
    c = Preconditioner(type="local", bf=a)
    a.Assemble()
    f.Assemble()

    BVP(bf=a, lf=f, gf=euz, pre=c, maxsteps=10000, prec=1.e-10).Do()
    er = compute_error(euz, sep, exactu)

    return (er, euz, sep, X)


def solvewavedirect(mesh, p, F, q_zero, mu_zero, cwave,
                    exactu=None, epsil=1.e-9):

    a, f, X, sep = makeforms(mesh, p, F, q_zero, mu_zero, cwave, epsil=epsil)

    euz = GridFunction(X)
    a.Assemble()
    f.Assemble()

    f.vec.data += a.harmonic_extension_trans * f.vec
    euz.vec.data = a.mat.Inverse(X.FreeDofs(True)) * f.vec
    euz.vec.data += a.harmonic_extension * euz.vec
    euz.vec.data += a.inner_solve * f.vec

    er = compute_error(euz, sep, exactu)

    return (er, euz, sep, X)


def print_rates(err, hs):
    """
    Tabulate eigenvalue errors (ewerr) with their rates assuming that
    meshsize (hs) is decreased successively by half.
    """

    print("\nCONVERGENCE RATES:---------+")
    print("---------------------------+")
    for i in range(len(err)):
        if i > 0:
            line = 'h/%-5d ' % 2**i
        else:
            line = 'h=%4.3f ' % hs[i]

        line += ' %7.4e  ' % abs(err[i])

        if abs(err[i]) > 1.e-15 and i > 0:
            r = log(err[i-1]/err[i])/log(2)
            line += format(r, '+5.2f')
        else:
            line += '  -  '
        print(line + ' |')
    print("---------------------------+")


###########################################

if __name__ == '__main__':

    ngs.ngsglobals.msg_level = 1
    demos = ['2d_triangular',
             '2d_rectangular',
             '3d_tetrahedral',
             '3d_hexahedral']

    example = demos[3]     # Choose the demo to run

    if example[:3] == '2d_':

        maxr = 5       # max refinement
        h0 = 0.25      # coarsest mesh size
        p = 1          # polynomial degree

        if example[3:] == 'triangular':
            mesh = ngs.Mesh(unit_square.GenerateMesh(maxh=h0))
        elif example[3:] == 'rectangular':
            mesh = ngs.Mesh(unit_square.GenerateMesh(maxh=h0,
                                                     quad_dominated=True))

        q_zero = 'bottom'  # mesh boundary parts where q = 0,
        mu_zero = 'bottom|right|left'         # where mu = 0.

        x = ngs.x
        t = ngs.y
        cwave = 1
        f = 2*pi*pi*sin(pi*x)*cos(2*pi*t) + pi*pi*sin(pi*x)*sin(pi*t)*sin(pi*t)
        F = CoefficientFunction((0, f))
        exactu = CoefficientFunction((pi*cos(pi*x)*sin(pi*t)*sin(pi*t),
                                      pi*sin(pi*x)*sin(2*pi*t)))

        hs = []
        er = []
        for l in range(maxr):
            h = h0 * 2**(-l)
            hs.append(h)
            e, *rest = solvewavedirect(mesh, p, F, q_zero, mu_zero, cwave, exactu)
            # e, *rest = solvewave(mesh, p, F, q_zero, mu_zero, cwave, exactu)
            er.append(e)
            mesh.Refine()

        print_rates(er, hs)

    elif example == '3d_tetrahedral':

        maxr = 4       # max refinement
        h0 = 1.0       # coarsest mesh size
        p = 0          # polynomial degree

        mesh = ngs.Mesh(unit_cube.GenerateMesh(maxh=h0))
        q_zero = 'bottom'      # mesh boundary parts where q = 0,
        mu_zero = 'bottom|right|left|front|back'  # where mu = 0.
        x = ngs.x
        y = ngs.y
        t = ngs.z
        cwave = 1
        F = CoefficientFunction((0, 0, sin(pi*y)*sin(pi*x)*(2+2*pi*pi*t*t)))
        exactu = CoefficientFunction((pi*cos(pi*x)*sin(pi*y)*t*t,
                                      pi*cos(pi*y)*sin(pi*x)*t*t,
                                      2*sin(pi*x)*sin(pi*y)*t))

        hs = []
        er = []
        for l in range(maxr):
            h = h0 * 2**(-l)
            hs.append(h)
            e, *rest = solvewavedirect(mesh, p, F, q_zero, mu_zero, cwave, exactu)
            er.append(e)
            mesh.Refine()

        print_rates(er, hs)

    elif example == '3d_hexahedral':

        maxr = 4       # max refinement
        h0 = 0.5       # coarsest mesh size
        p = 0          # polynomial degree

        q_zero = 'Z0'               # mesh boundary parts where q = 0,
        mu_zero = 'Z0|X0|X1|Y0|Y1'  # where mu = 0.
        x = ngs.x
        y = ngs.y
        t = ngs.z
        cwave = 1
        F = CoefficientFunction((0, 0, sin(pi*y)*sin(pi*x)*(2+2*pi*pi*t*t)))
        exactu = CoefficientFunction((pi*cos(pi*x)*sin(pi*y)*t*t,
                                      pi*cos(pi*y)*sin(pi*x)*t*t,
                                      2*sin(pi*x)*sin(pi*y)*t))

        hs = []
        er = []
        for l in range(maxr):
            n = 2**l
            h = h0 / n
            hs.append(h)
            mesh = ngs.Mesh(GenerateCubeMesh(n))
            e, *rest = solvewavedirect(mesh, p, F, q_zero, mu_zero, cwave, exactu)
            er.append(e)

        print_rates(er, hs)
