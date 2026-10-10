"""Polarization evolution through a homogeneous brane under S10's supplied MAIN action
and S10's XFORM_ANISO inertia control (FORM control: changes the kinetic quadratic form).

Action (S10_two_transverse_photons.md:187-193):
    L = (rho/2) sum_j (d_t u_j)^2 - (mu/2) S_curl,  S_curl = (1/2) sum_ij (d_i u_j - d_j u_i)^2
ANISO (S10:88-92): kinetic form s*(d_t u_1)^2 + (d_t u_2)^2 + (d_t u_3)^2, s>0, s!=1.

For a plane wave u = a exp(i(k.x - w t)) the Euler-Lagrange equation is
    w^2 R a = K(k) a,
with R the inertia matrix and K the stiffness matrix built here from the action by
differentiation (not typed).  We print the generalized eigenvalues w^2, then the relative
transverse phase accumulated over a path length Lp at fixed frequency and the Stokes vector
of an arbitrary input polarization before and after.  Operands and residuals are printed;
no conclusion is stated.
"""
import sympy as sp

rho, mu, s, k, w, Lp = sp.symbols('rho mu s k omega L_p', positive=True)
t = sp.symbols('t', real=True)
x = sp.symbols('x1:4', real=True)
a = sp.symbols('a1:4')
nx, ny, nz = sp.Rational(2, 7), sp.Rational(3, 7), sp.Rational(6, 7)   # generic unit direction
n = [nx, ny, nz]
print('unit check |n|^2 =', sum(c**2 for c in n))

def stiffness_matrix():
    # quadratic form of S_curl for u_j = a_j * phase, gradient d_i u_j -> i k n_i a_j
    G = sp.Matrix(3, 3, lambda i, j: k * n[i] * a[j])
    Scurl = sp.Rational(1, 2) * sum((G[i, j] - G[j, i])**2 for i in range(3) for j in range(3))
    V = (mu / 2) * Scurl
    return sp.hessian(V, a)

def inertia_matrix(kind):
    T = (rho / 2) * (s * a[0]**2 + a[1]**2 + a[2]**2) if kind == 'ANISO' else (rho / 2) * sum(ai**2 for ai in a)
    return sp.hessian(T, a)

def stokes(Ex, Ey):
    S0 = sp.Abs(Ex)**2 + sp.Abs(Ey)**2
    S1 = sp.Abs(Ex)**2 - sp.Abs(Ey)**2
    S2 = 2 * sp.re(sp.conjugate(Ex) * Ey)
    S3 = 2 * sp.im(sp.conjugate(Ex) * Ey)
    return [sp.nsimplify(sp.simplify(v)) for v in (S0, S1, S2, S3)]

for kind in ('MAIN', 'ANISO'):
    K = stiffness_matrix()
    R = inertia_matrix(kind)
    lam = sp.symbols('lam')
    charpoly = sp.factor((K - lam * R).det())
    roots = sp.solve(charpoly, lam)
    print('\n===', kind, '===')
    print('char poly det(K - lam R) =', charpoly)
    print('w^2 roots =', [sp.simplify(r) for r in roots])
    pos = [r for r in roots if sp.simplify(r) != 0]
    # multiplicities
    mult = sp.roots(sp.Poly(charpoly, lam))
    print('root multiplicities =', {sp.simplify(r): m for r, m in mult.items()})
    nonzero = [r for r, m in mult.items() for _ in range(m) if sp.simplify(r) != 0]
    print('nonzero w^2 operands (with multiplicity) =', [sp.simplify(r) for r in nonzero])
    if len(nonzero) == 2:
        resid = sp.simplify(nonzero[0] - nonzero[1])
        print('residual w^2_(1) - w^2_(2) =', resid)
        # phase velocities v_i = w / k_i at fixed w: w^2 = c_i^2 k^2 => k_i = w / c_i
        c1 = sp.sqrt(sp.simplify(nonzero[0] / k**2)); c2 = sp.sqrt(sp.simplify(nonzero[1] / k**2))
        dphi = sp.simplify(w * Lp * (1 / c1 - 1 / c2))
        print('relative transverse phase over L_p at fixed omega, dphi =', dphi)
        # Stokes vector of arbitrary input (Ex, Ey) = (cos th, e^{i d} sin th) in the eigenbasis
        th, d = sp.Rational(1, 5), sp.Rational(1, 3)
        Ein = (sp.cos(th), sp.exp(sp.I * d) * sp.sin(th))
        # sample numeric values for the ANISO coefficient, path and frequency
        subs = {s: sp.Rational(3, 2), rho: 1, mu: 1, w: 1, Lp: 10}
        dphi_num = sp.N(dphi.subs(subs))
        Eout = (Ein[0], Ein[1] * sp.exp(sp.I * dphi_num))
        Sin = [sp.N(v) for v in stokes(*Ein)]
        Sout = [sp.N(v) for v in stokes(*Eout)]
        print('  (sample: s=3/2, rho=mu=omega=1, L_p=10) dphi =', dphi_num)
        print('  Stokes in  =', Sin)
        print('  Stokes out =', Sout)
        print('  Stokes residual out-in =', [sp.N(o - i_) for o, i_ in zip(Sout, Sin)])
