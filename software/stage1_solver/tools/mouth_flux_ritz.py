"""Finite-basis 3D Maxwell/coherence collar, in units B3, at matched radius.

Use radial gauge and fix the residual exact part of the angular connection
to zero at the inner mouth.  Both types of tangential vector spherical
harmonics and both complex monopole sections vary throughout the collar.
Matching the tangential trace to the least-energy exterior continuation
gives l/Lambda times the squared coexact trace and zero exact-trace cost.
The minimization is nonconvex; outputs describe attained trial fields.
"""

import argparse
from math import factorial, pi
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre, legfit, leggauss
from scipy.optimize import minimize
from scipy.special import eval_jacobi, lpmv


ASPECT = 1.25
GAMMA = .2
B3 = .225  # In units 4*pi*K**2/(u*R).


def quadrature(order, lo, hi):
    nodes, weights = leggauss(order)
    return lo + (1+nodes)*(hi-lo)/2, weights*(hi-lo)/2


class Collar:
    def __init__(self, n, gains, shells, flux_order, matter_radial, flux_radial, extra=0):
        self.n, self.gains = n, gains
        self.ms = (n+2, n-2)
        span = max(abs(m) for m in self.ms) + 2*shells
        total_degree = max(2*span, 2*flux_order + span)
        z, wz = quadrature(total_degree//2+3+extra, 0, 1)
        nphi = total_degree+3+2*extra
        self.z = np.repeat(z, nphi)
        self.x = 2*self.z-1
        self.phi = np.tile(2*pi*np.arange(nphi)/nphi, len(z))
        self.sin = 2*np.sqrt(self.z*(1-self.z))
        self.wa = np.repeat(wz/nphi, nphi)
        self.r, wr = quadrature(max(7, 2*matter_radial+flux_radial+3)+extra,
                                1, ASPECT)
        s = 2*(self.r-1)/(ASPECT-1)-1
        self.w = wr[:, None]*self.wa
        self.r2 = self.r[:, None]**2
        self.wp = wr
        order = max(matter_radial, flux_radial)
        self.V = np.stack([Legendre.basis(j)(s) for j in range(order+1)], axis=1)
        self.dV = np.stack([Legendre.basis(j).deriv()(s)*2/(ASPECT-1)
                           for j in range(order+1)], axis=1)
        self.nr = matter_radial+1
        self.sections = [self._section(m, shells) for m in self.ms]
        self.gauge = self._gauge(flux_order, flux_radial)
        self.keys = []
        for h, (labels, _, _, _, _) in enumerate(self.sections):
            for part in ("re", "im"):
                self.keys.extend((part, h, j, p, k) for j in range(self.nr)
                                 for p, k in labels)
        for typ, labels, V, _, _, _, _, _ in self.gauge:
            self.keys.extend((typ, j+int(typ=="grad"), l, q)
                             for j in range(V.shape[1]) for l, q in labels)

    def _section(self, m, shells):
        z, x = self.z, self.x
        labels, basis, polar, azimuth, scale = [], [], [], [], []
        for p in range(min(0, m)-shells, max(0, m)+shells+1):
            a, b = abs(p), abs(p-m)
            edge = (a+b-abs(m))//2
            for j in range(shells-edge+1):
                labels.append((p, j))
                norm = (factorial(j+a)*factorial(j+b) /
                        ((2*j+a+b+1)*factorial(j)*factorial(j+a+b)))
                pref = z**(b/2)*(1-z)**(a/2)*np.exp(1j*p*self.phi)/np.sqrt(norm)
                poly = eval_jacobi(j, a, b, x)
                deriv = (j+a+b+1)*eval_jacobi(j-1, a+1, b+1, x) if j else 0
                f = pref*poly
                basis.append(f)
                polar.append(-self.sin/2*pref*(deriv + poly*(b/(2*z)-a/(2*(1-z)))))
                azimuth.append(1j*(p-m*(1-z))/self.sin*f)
                ell = j+(a+b)/2
                lam = ell*(ell+1)-m*m/4
                scale.append(1/np.sqrt(self.wp @ (
                    self.r2*self.dV[:, :self.nr]**2 + (1+lam)*self.V[:, :self.nr]**2)))
        B, Dt, Dp = map(np.array, (basis, polar, azimuth))
        normerr = np.max(abs((B*self.wa) @ B.conj().T-np.eye(len(labels))))
        lap = (Dt*self.wa) @ Dt.conj().T + (Dp*self.wa) @ Dp.conj().T
        expected = [(j+(abs(p)+abs(p-m))/2)*(j+(abs(p)+abs(p-m))/2+1)-m*m/4
                    for p, j in labels]
        if normerr > 1e-11 or np.max(abs(lap-np.diag(expected))) > 1e-10:
            raise ArithmeticError("Monopole basis failed angular metric check.")
        return labels, B, Dt, Dp, np.asarray(scale).T

    def _gauge(self, flux_order, flux_radial):
        labels, dt, dp = [], [], []
        for l in range(1, flux_order+1):
            for q in range(-l, l+1):
                labels.append((l, q))
                k = abs(q)
                coeff = np.sqrt((2*l+1)*(1+(k>0))*factorial(l-k)/factorial(l+k)
                                /(l*(l+1)))
                ang = np.sin(k*self.phi) if q<0 else np.cos(k*self.phi)
                dang = k*np.cos(k*self.phi) if q<0 else -k*np.sin(k*self.phi)
                P = lpmv(k, l, self.x)
                previous = lpmv(k, l-1, self.x) if k<l else 0
                dt.append(coeff*ang*(l*self.x*P-(l+k)*previous)/self.sin)
                dp.append(coeff*dang*P/self.sin)
        dt, dp = map(np.asarray, (dt, dp))
        identity = (dt*self.wa) @ dt.T + (dp*self.wa) @ dp.T
        cross = (dt*self.wa) @ dp.T - (dp*self.wa) @ dt.T
        if max(np.max(abs(identity-np.eye(len(labels)))), np.max(abs(cross))) > 1e-10:
            raise ArithmeticError("Vector harmonic basis failed angular metric check.")
        result = []
        for typ, Bt, Bp in (("curl", -dp, dt), ("grad", dt, dp)):
            j0 = int(typ=="grad")
            V = self.V[:, j0:flux_radial+1].copy()
            dV = self.dV[:, j0:flux_radial+1]
            if j0:
                V -= (-1.)**np.arange(j0, flux_radial+1)
            derivative = (dV.T*self.wp) @ dV
            stiffness = (V.T*(self.wp/self.r**2)) @ V
            H = np.array([GAMMA*(derivative +
                         (l*(l+1)*stiffness + l/ASPECT if not j0 else 0))
                          for l, q in labels])
            scale = 1/np.sqrt(np.diagonal(H, axis1=1, axis2=2).T)
            result.append((typ, labels, V, H, scale, Bt, Bp, dV))
        return result

    def unpack(self, vector):
        parts, position = [], 0
        for labels, _, _, _, scale in self.sections:
            d = self.nr*len(labels)
            real, imag = vector[position:position+2*d].reshape(2, self.nr, len(labels))
            parts.append((real + 1j*imag)*scale)
            position += 2*d
        for _, _, _, _, scale, _, _, _ in self.gauge:
            parts.append(vector[position:position+scale.size].reshape(scale.shape)*scale)
            position += scale.size
        return parts

    def physical(self, vector):
        parts = self.unpack(vector)
        arrays = [v for a in parts[:2] for v in (a.real.ravel(), a.imag.ravel())]
        arrays.extend(a.ravel() for a in parts[2:])
        return dict(zip(self.keys, np.concatenate(arrays)))

    def vector(self, physical):
        scales = ([k[-1].ravel() for k in self.sections for _ in range(2)]
                  + [k[4].ravel() for k in self.gauge])
        return np.array([physical.get(k, 0.) for k in self.keys])/np.concatenate(scales)

    def energy(self, vector, detail=False):
        parts = self.unpack(vector)
        At, Ap = np.zeros_like(self.w), np.zeros_like(self.w)
        field_cost, field_deriv = 0., []
        for a, (_, _, V, H, _, Bt, Bp, _) in zip(parts[2:], self.gauge):
            amps = V @ a
            At += amps @ Bt
            Ap += amps @ Bp
            Ha = np.einsum("qij,jq->iq", H, a)
            field_cost += np.sum(a*Ha)
            field_deriv.append(Ha)
        fields = []
        V, dV = self.V[:, :self.nr], self.dV[:, :self.nr]
        for c, (_, B, Dt, Dp, _) in zip(parts, self.sections):
            coeff = V @ c
            f = coeff @ B
            fields.append((f, (dV @ c) @ B, coeff @ Dt-1j*At*f, coeff @ Dp-1j*Ap*f))
        S = sum(abs(f[0])**2 for f in fields)
        material = np.sum(self.w*self.r2*S**2)/2
        jt, jp = np.zeros_like(S), np.zeros_like(S)
        grads, radial_energy = [], 0.
        for gain, (f, dr, polar, az), (_, B, Dt, Dp, scale) in zip(
                self.gains, fields, self.sections):
            rad_h = np.sum(self.w*self.r2*abs(dr)**2)
            radial_energy += rad_h
            material += rad_h + np.sum(self.w*(abs(polar)**2 + abs(az)**2 - gain*abs(f)**2))
            d0 = (self.r2*S-gain)*f + 1j*(At*polar+Ap*az)
            grad = V.T @ ((self.w*d0) @ B.conj().T +
                          (self.w*polar) @ Dt.conj().T + (self.w*az) @ Dp.conj().T)
            grad += dV.T @ ((self.w*self.r2*dr) @ B.conj().T)
            grads.extend(((grad.real*scale).ravel(), (grad.imag*scale).ravel()))
            jt += np.imag(polar.conj()*f)
            jp += np.imag(az.conj()*f)
        for Ha, (_, _, Vg, _, scale, Bt, Bp, _) in zip(field_deriv, self.gauge):
            grad = Ha + Vg.T @ ((self.w*jt) @ Bt.T + (self.w*jp) @ Bp.T)
            grads.append((grad*scale).ravel())
        bg = GAMMA*self.n**2/4
        total = (bg+field_cost+material)/B3
        if detail:
            return np.array([total, bg/B3, field_cost/B3, material/B3, radial_energy/B3])
        return total, 2*np.concatenate(grads)/B3

    def seed(self, simple=False):
        n = self.n
        if n == 0 and simple:
            k = (0, -2)
            amplitudes = [np.sqrt(.8)*self.z, np.sqrt(.8)*(1-self.z)]
            b = -.25 + 0*self.r
        elif n == 0 and self.gains[0] == self.gains[1]:
            k = (2, 0)
            amplitudes = [(1-self.z)*(1.051+.192*self.x+.002*self.x**2),
                          self.z*(1.051-.192*self.x+.002*self.x**2)]
            b = .355-.056*(self.r**2-1)
        elif n == 0:
            k = (0, -2)
            amplitudes = [self.z*(1.146-.278*self.x+.024*self.x**2),
                          (1-self.z)*(.631+.042*self.x-.009*self.x**2)]
            b = -.302+.048*(self.r**2-1)
        elif n == -1:
            k = (0, 0)
            amplitudes = [np.sqrt(self.z)*(1.308-.352*self.x+.093*self.x**2), 0*self.z]
            b = -.270+.042*(self.r**2-1)
        else:
            k = (0, 0)
            amplitudes = [0*self.z, np.sqrt(self.z)*(1.087-.233*self.x+.049*self.x**2)]
            b = .213-.034*(self.r**2-1)
        values = {}
        for h, (p, a, (labels, B, _, _, _)) in enumerate(zip(k, amplitudes, self.sections)):
            projections = B.conj() @ (self.wa*a*np.exp(1j*p*self.phi))
            for (q, j), c in zip(labels, projections):
                values["re",h,0,q,j] = c.real
                values["im",h,0,q,j] = c.imag
        for j, c in enumerate(legfit(2*(self.r-1)/(ASPECT-1)-1, -np.sqrt(2/3)*b, 2)):
            values["curl",j,1,0] = c
        return self.vector(values)


def run_case(n, gains, levels, directory):
    rng = np.random.default_rng(37851+n)
    previous = None
    name = f"n{n}_{gains[1]}"
    for li, level in enumerate(levels):
        c = Collar(n, gains, *level)
        fine = Collar(n, gains, *level, extra=4)
        seed = c.seed()
        if n==0 and gains[1]==1.5:
            control = c.energy(c.seed(True), True)[0]
            if abs(control + 338/2025)>2e-12:
                raise ArithmeticError(f"Flux patch normalization or monopole coupling: {control}")
        test = seed + .015*rng.normal(size=len(seed))
        d = rng.normal(size=len(seed))
        d /= np.linalg.norm(d)
        step = 2e-6
        numerical = (c.energy(test+step*d)[0]-c.energy(test-step*d)[0])/(2*step)
        grad_error = abs(numerical - c.energy(test)[1] @ d)
        if grad_error > 3e-8:
            raise ArithmeticError(f"Complex matter/flux gradient mismatch: {grad_error}")
        if previous is None:
            # One exact-axis seed and two unrelated complex, all-m starts.
            starts = [seed, seed+.10*rng.normal(size=len(seed)),
                      .13*rng.normal(size=len(seed))]
        else:
            seed = c.vector(previous)
            starts = [seed, seed+.07*rng.normal(size=len(seed))]
        results = []
        for si, x0 in enumerate(starts):
            res = minimize(c.energy, x0, jac=True, method="L-BFGS-B",
                           options={"ftol": 4e-15, "gtol": 1e-8, "maxcor": 16,
                                    "maxiter": 2000, "maxls": 40})
            grad = np.linalg.norm(res.jac, np.inf)
            energy_fine = fine.energy(res.x, True)
            qerror = abs(energy_fine[0]-res.fun)
            print(f"{name} level={level} start={si} E={res.fun:.12f} "
                  f"dof={len(seed)} it={res.nit} ok={res.success} grad={grad:.2e} "
                  f"quadrature={qerror:.2e} gradient_check={grad_error:.2e}", flush=True)
            if not res.success or grad>8e-6 or qerror>2e-10:
                raise ArithmeticError(f"Optimizer or quadrature unresolved: {res.message}")
            results.append((res.fun, res.x, energy_fine))
        energy, best, values = min(results, key=lambda o: o[0])
        previous = c.physical(best)
        if directory:
            np.savez(directory/f"{name}_l{li}.npz", vector=best, level=level,
                     n=n, gains=gains, energies=values)
        print(f"{name} best {level} E={energy:.12f} mono={values[1]:.10f} "
              f"patch={values[2]:.10f} matter={values[3]:.10f} "
              f"matter_radial={values[4]:.10f}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fields", type=Path, help="directory for attained Ritz coefficients")
    args = parser.parse_args()
    if args.fields:
        args.fields.mkdir(parents=True, exist_ok=True)
    levels = [(1, 2, 0, 2), (2, 3, 0, 2), (2, 3, 1, 2), (3, 4, 2, 3)]
    for n, gains in ((0, (1.5,1.5)), (-1, (1.5,1.5)), (0, (1.5,1.2)),
                     (-1, (1.5,1.2)), (1, (1.5,1.2))):
        run_case(n, gains, levels, args.fields)


if __name__ == "__main__":
    main()
