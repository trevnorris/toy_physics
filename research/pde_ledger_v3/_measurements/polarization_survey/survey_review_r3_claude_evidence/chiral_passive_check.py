#!/usr/bin/env python3
# Question: can a parity-odd (chiral) quadratic stiffness term split the two transverse helicities
# while being lossless (Hermitian stiffness, real omega^2) and reciprocal (K_ij(k) = K_ji(-k)),
# i.e. without any non-passive coupling or reservoir?
# Stiffness tensor for an isotropic in-plane displacement u with curl stiffness mu and a
# parity-odd third-derivative term eta: K_ij(k) = mu (k^2 d_ij - k_i k_j) + i eta k^2 eps_ijl k_l.
import sympy as sp
mu, eta, rho = sp.symbols('mu eta rho', positive=True)
kx, ky, kz = sp.symbols('k_x k_y k_z', real=True)
kv = sp.Matrix([kx, ky, kz]); k2 = kx**2 + ky**2 + kz**2
def K(kvec):
    kk = kvec.dot(kvec)
    M = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            M[i, j] = mu*(kk*sp.KroneckerDelta(i, j) - kvec[i]*kvec[j]) \
                + sp.I*eta*kk*sum(sp.LeviCivita(i, j, l)*kvec[l] for l in range(3))
    return M
Kk = K(kv)
print("K(k) - K(k)^dagger (Hermiticity residual):", sp.simplify(Kk - Kk.H))
print("K(k) - K(-k)^T (Onsager reciprocity residual):", sp.simplify(Kk - K(-kv).T))
# Parity: under u -> -u(-x), k -> -k, a parity-even tensor satisfies K(-k) = K(k).
print("K(-k) - K(k) (parity-odd part):", sp.simplify(K(-kv) - Kk))
kz_ = sp.symbols('k', positive=True)
Kz = K(sp.Matrix([0, 0, kz_]))
print("K for k along z:", Kz)
for val, mult, vecs in Kz.eigenvects():
    print("eigenvalue rho*omega^2 =", sp.simplify(val), " multiplicity", mult, " eigenvector", [list(v) for v in vecs])
# Control: eta = 0 restores the degenerate transverse pair.
print("eta=0 eigenvalues:", Kz.subs(eta, 0).eigenvals())
