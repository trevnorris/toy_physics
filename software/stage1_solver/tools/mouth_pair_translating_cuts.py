"""Linear outgoing torsional field of two prescribed translating spherical cuts.

Lengths in cut radii. Each surface has a fixed primitive sector, homogeneous
attachment in its own frame, and free site rolls. Axisymmetric motion along
the line of centers. The static field and the outgoing waves are solved jointly
on the two separate surfaces, retaining every harmonic through the stated L.
"""
from dataclasses import dataclass
from math import pi

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_legendre, spherical_jn, spherical_yn


FLUX = 2*pi
MODULUS, INERTIA = 4., 8.
SPEED = np.sqrt(MODULUS/INERTIA)
AMPLITUDE = .01


def polynomial(l, u):
    p = eval_legendre(l,u)
    if l == 0:
        return p, np.zeros_like(u)
    # The grid uses interior Gauss nodes (no derivatives at a pole).
    sin = np.sqrt(np.maximum(0., (1-u)*(1+u)))
    return p, l*(eval_legendre(l-1,u)-u*p)/sin


def hankel(l, x, derivative=False):
    return spherical_jn(l,x,derivative)+1j*spherical_yn(l,x,derivative)


@dataclass
class Pair:
    centers: np.ndarray
    lmax: int
    binding: float  # G R / C; infinity pins the co-moving normal flux profile

    def __post_init__(self):
        self.centers = np.asarray(self.centers)
        self.charges = np.array([FLUX,-FLUX])[:len(self.centers)]
        self.u, self.weights = leggauss(max(48,3*self.lmax))
        self.sin = np.sqrt(1-self.u**2)
        radial, tangential = [], []
        for l in range(1,self.lmax+1):
            pl, tl = polynomial(l,self.u)
            norm = self.weights*(2*l+1)/2
            radial.append(norm*pl)
            # The tangential projection norm has an extra l(l+1);
            # this cancels the same factor in the surface response.
            tangential.append(norm*tl/self.binding)
        self.pr, self.pt = np.array(radial), np.array(tangential)
        self.monopoles = self.charges/(4*pi)
        self.static = self.solve_static()

    def field(self, origin, radius, k=None, monopole=False, shift=0.):
        """Radial/tangential D and azimuthal angle a on a sphere about origin.
        shift translates the sample points along z while keeping the test
        normals fixed, as required by translation of an entire prepared cut.
        Each column is a unit source at its own center, no self-field removal.
        """
        cols = []
        for center in self.centers:
            dz = (origin-center)+shift+radius*self.u
            dr = radius*self.sin
            r = np.hypot(dr,dz)
            cos, sin = dz/r, dr/r
            normal_dot = self.u*cos+self.sin*sin
            theta_dot = self.sin*cos-self.u*sin
            for l in ((0,) if monopole else range(1,self.lmax+1)):
                p, t = polynomial(l,cos)
                if k is None:
                    scale = r**(-l-2)
                    radial, theta, angle = (l+1)*scale*p, scale*t, np.zeros_like(r)
                else:
                    h0 = hankel(l,k)
                    f = hankel(l,k*r)/h0
                    fp = k*hankel(l,k*r,True)/h0
                    radial, theta, angle = l*(l+1)*f*p/r, -(fp+f/r)*t, f*t
                cols.append([normal_dot*radial+theta_dot*theta,
                             -theta_dot*radial+normal_dot*theta, angle])
        return np.stack(cols,axis=-1)

    def project(self, field):
        return self.pr@field[0]+self.pt@field[1]

    def solve_static(self):
        mat, rhs = [], []
        for z in self.centers:
            mat.append(self.project(self.field(z,1)))
            rhs.append(-self.project(self.field(z,1,monopole=True))@self.monopoles)
        return np.linalg.solve(np.vstack(mat),np.concatenate(rhs))

    def dipole(self):
        return self.charges@self.centers+4*pi*sum(self.static[::self.lmax])

    def background(self,z,shift):
        return (self.field(z,1,shift=shift)@self.static
                +self.field(z,1,monopole=True,shift=shift)@self.monopoles)

    def emit(self,k,motion):
        mat, rhs = [], []
        step=1e-4
        for z, x in zip(self.centers,motion):
            mat.append(self.project(self.field(z,1,k)))
            derivative = ((self.background(z,-2*step)-self.background(z,2*step))
                          +8*(self.background(z,step)-self.background(z,-step)))/(12*step)
            rhs.append(-x*self.project(derivative))
        source = np.concatenate(rhs)
        coeff = np.linalg.solve(np.vstack(mat),source)
        return coeff

    def flux(self, coeff, origin, radius, k):
        _, theta, angle = self.field(origin,radius,k)@coeff
        angular_velocity = -1j*k*SPEED*angle
        radial_work = MODULUS*(theta*angular_velocity.conj()).real/2
        return 2*pi*radius**2*(self.weights@radial_work)

    def outgoing(self, coeff, k):
        # Outgoing amplitude of r a_phi exp(-ikr) on an enclosing sphere.
        amplitude = np.zeros_like(self.u,dtype=complex)
        for z, block in zip(self.centers,coeff.reshape(-1,self.lmax)):
            phase = np.exp(-1j*k*z*self.u)
            for l,c_l in enumerate(block,1):
                t = polynomial(l,self.u)[1]
                amplitude += c_l*phase*(-1j)**(l+1)*t/(k*hankel(l,k))
        return pi*MODULUS*SPEED*k**2*(self.weights@abs(amplitude)**2)


def dipole_power(p,k):
    return MODULUS*SPEED*k**4*abs(p)**2/(12*pi)


def main():
    for binding in (1., 2., 8., np.inf):
        single = Pair(np.array([0.]),8,binding)
        k = .2
        power = single.outgoing(single.emit(k,[AMPLITUDE]),k)
        ratio = power/dipole_power(FLUX*AMPLITUDE,k)
        expected = abs(np.exp(-1j*k)/(1-1j*k-k**2/(1+binding)))**2
        print(f"isolated g={binding:g} kR={k:g} power_ratio={ratio:.10f} analytic={expected:.10f}")
        assert abs(ratio-expected)<1e-8
    motion = AMPLITUDE*np.array([1.,-1.])
    point_dipole = 2*FLUX*AMPLITUDE
    for d,binding in ((4.,1.),(4.,2.),(4.,8.),(4.,np.inf),(6.,8.),(10.,8.)):
        centers = np.array([-d/2,d/2])
        pair = Pair(centers,10,binding)
        step = .001
        response = (Pair(centers+step*motion,10,binding).dipole()
                    - Pair(centers-step*motion,10,binding).dipole())/(2*step)
        print(f"d/R={d:g} g={binding:g} static_dipole={pair.dipole():.10f} "
              f"dp_per_relative_translation_over_q={response/point_dipole:.9f}")
        for k in (.008,.04,.12):
            coeff = pair.emit(k,motion)
            work_out = pair.outgoing(coeff,k)
            near = sum(pair.flux(coeff,z,1.22,k) for z in centers)
            far = pair.flux(coeff,0,d+3/k,k)
            common = pair.outgoing(pair.emit(k,[AMPLITUDE]*2),k)
            # Second partial-wave truncation is an independent resolution check.
            coarse = Pair(centers,6,binding)
            coarse_work = coarse.outgoing(coarse.emit(k,motion),k)
            balance = max(abs(far/work_out-1),abs(near/work_out-1))
            print(f"d/R={d:g} g={binding:g} kR={k:g} outward_work_rate={work_out:.10g} "
                  f"power_over_static_dipole={work_out/dipole_power(response,k):.9f} "
                  f"common_over_opposed={common/work_out:.9g} control_sphere_error={balance:.3g} "
                  f"L6_vs_L10={abs(coarse_work/work_out-1):.3g}")
            assert work_out>0 and balance < 2e-5 and abs(coarse_work/work_out-1)<1e-4
            if k==.008:
                assert abs(work_out/dipole_power(response,k)-1)<.002
        assert not np.any(pair.emit(.12,[0.,0.]))
        if binding==2:
            assert abs(response/point_dipole-1)<1e-8


if __name__ == "__main__":
    main()
