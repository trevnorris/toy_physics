"""Static response to a depth twist on a flat ordered slab with shear-free faces."""

import numpy as np
import sympy as sp
from scipy.linalg import solve_banded


def emit(label, value):
    print(f"{label} = {value}", flush=True)


def uniform_slab():
    a, b, k, q, mu, eta, z = sp.symbols("a b k q mu eta z", positive=True)
    # Half slab. Unit source at b, difference or sum with its mirror on the
    # negative half. z=0 is Dirichlet for odd, Neumann for even.
    def unit_load(left_shape):
        right_shape = sp.cosh(q*(a-z))
        wronskian = (sp.diff(left_shape,z)*right_shape
                     -left_shape*sp.diff(right_shape,z)).subs(z,b)
        scale = mu*sp.trigsimp(wronskian)
        left = left_shape*right_shape.subs(z,b)/scale
        right = left_shape.subs(z,b)*right_shape/scale
        residuals = [sp.diff(left,z,2)-q*q*left, sp.diff(right,z,2)-q*q*right,
                     sp.diff(right,z).subs(z,a), (left-right).subs(z,b),
                     (mu*sp.diff(left-right,z)).subs(z,b)-1]
        assert all(sp.trigsimp(v)==0 for v in residuals)
        return left, 2*left.subs(z,b)

    odd, Kodd = unit_load(sp.sinh(q*z))
    even, Keven = unit_load(sp.cosh(q*z))
    assert odd.subs(z,0)==sp.diff(even,z).subs(z,0)==0
    qk = k*sp.sqrt(eta/mu)
    pole = lambda K: sp.simplify(sp.limit(k*k*K.subs(q,qk),k,0))
    emit("twist-to-twist kernel", Kodd)
    emit("same-sign displacement kernel", Keven)
    emit("zero-wave-number twist compliance", sp.limit(Kodd,q,0))
    emit("Coulomb residues, difference and sum", [pole(Kodd),pole(Keven)])

    H, r = sp.symbols("H r", positive=True)
    m1 = sp.pi/H*sp.sqrt(mu/eta)
    emit("inverse range of first mixed-shear mode", m1)
    # b=H/4: sin^2(n*pi/4)=1/2 for every odd n. Divide the odd response by
    # the 1/r contribution of the same-sign control of equal source strength.
    m = sp.symbols("m", integer=True, nonnegative=True)
    rel = sp.summation(sp.exp(-(2*m+1)*sp.pi*r/H),(m,0,sp.oo))
    ratios = [float(rel.subs(r,i*H)) for i in (1,2,4)]
    emit("odd response relative to even Coulomb term at r/H=1,2,4; isotropic", ratios)
    isotropic = {a:sp.Rational(1,2), b:sp.Rational(1,4), mu:1, eta:1}
    return [float(f.subs(q,qk).subs(isotropic).subs(k,1)) for f in (Kodd,Keven)] + ratios


def longitudinal_control():
    # Curl plus positive 4D divergence energy. For u_x=U sin(kx),
    # u_w=W cos(kx), P=mu*(U'+k W), N=B*(W'+k U).
    z, k, a = sp.symbols("z k a", positive=True)
    U, W, fx, fw = [sp.Function(n)(z) for n in ("U","W","fx","fw")]
    mu, B = sp.symbols("mu B",positive=True)
    energy = (mu*(sp.diff(U,z)+k*W)**2+B*(sp.diff(W,z)+k*U)**2)/2
    P = sp.diff(energy,sp.diff(U,z))
    N = sp.diff(energy,sp.diff(W,z))
    ex = sp.diff(energy,U)-sp.diff(P,z)-fx
    ew = sp.diff(energy,W)-sp.diff(N,z)-fw
    elim = -sp.diff(P,z,2)+k*k*P-sp.diff(fx,z)-k*fw
    identity = sp.simplify(elim-sp.diff(ex,z)-k*ew)
    # Unit midplane normal force, no tangential body force. P is even and
    # vanishes at both faces. The derivative jump fixes its normalization.
    shape = sp.sinh(k*(a-z))
    half = shape*k/(-2*sp.diff(shape,z).subs(z,0))
    fraction = sp.simplify(2*k*sp.integrate(half,(z,0,a)))
    jump = sp.simplify(2*sp.diff(half,z).subs(z,a)/k)
    emit("coupled stress equation residual", identity)
    emit("normal unit force, in-plane mixed stress fraction and face-stress jump",
         [fraction,jump])
    emit("zero-wave-number fractions", [sp.limit(fraction,k,0),sp.limit(jump,k,0)])
    assert identity==0 and sp.simplify(fraction-jump-1)==0
    return [float(expr.subs({a:sp.Rational(1,2),k:1})) for expr in (fraction,jump)]


def variable_slab(elements, k, odd):
    # Linear elements on 0 <= w <= 1. The reflected probe is at w=1/2.
    nodes = np.linspace(0,1,elements+1)
    h = 1/elements
    diag = np.zeros(elements+1)
    upper = np.zeros(elements)
    for node,left in enumerate(nodes[:-1]):
        for t in (0.5-np.sqrt(3)/6,0.5+np.sqrt(3)/6):
            w=left+h*t
            grad=(1+w*w)/(2*h)
            mass=k*k*(2+w*w)*h/2
            diag[node] += grad+mass*(1-t)**2
            diag[node+1] += grad+mass*t*t
            upper[node] += -grad+mass*t*(1-t)
    rhs=np.zeros_like(nodes)
    rhs[elements//2]=1
    if odd:
        upper[0]=0
        diag[0]=1
    matrix=np.zeros((3,elements+1))
    matrix[0,1:]=matrix[2,:-1]=upper
    matrix[1]=diag
    profile=solve_banded((1,1),matrix,rhs)
    return 2*profile[elements//2]


if __name__=="__main__":
    common = uniform_slab()+longitudinal_control()
    # Cross-check the extended, even moduli against independent continuum
    # identities. k=0 is nonsingular on the odd half-domain.
    target=2*np.arctan(.5)
    errors=[]
    for cells in (32,64,128):
        compliance=variable_slab(cells,0,True)
        residue=.01**2*variable_slab(cells,.01,False)
        errors.append(abs(compliance-target))
        emit(f"variable moduli, elements {cells}, odd at k=0 and even k^2 G at k=.01",
             [compliance,residue])
    emit("variable exact odd compliance and even residue", [target,6/7])
    print("comparison",*(f"{x:.14g}" for x in common+[target,6/7]),flush=True)
    assert errors[2]<errors[1]<errors[0] and errors[2]<5e-6
    assert abs(residue-6/7)<1e-4
