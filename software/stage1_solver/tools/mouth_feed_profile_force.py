"""Small-acceleration response of the translated amplitude-restoring surface energy."""

import json

import numpy as np

from mouth_prescribed_acceleration import MovingFeed
from signed_mouth_collar import RUNS, load_drives


class FeedTangent:
    """Linearize the semidiscrete fixed-intake feed at its saved stationary field."""

    def __init__(self, feed):
        self.p = p = feed.patch
        self.psi = feed.base
        self.face = self.psi[-1]
        self.n = abs(self.face)**2
        self.u = np.sqrt(self.n)
        h = np.maximum(1-p.x**2-feed.transverse_squared,0)
        self.g = h**3
        self.g1 = 6*p.x*h**2             # derivatives with respect to the rim position X
        self.g2 = -6*h**2+24*p.x**2*h
        self.ref_grad = np.fft.irfft(1j*feed.k*feed.reference_spectrum, n=p.n, axis=0)
        self.scale = feed.rate/(p.dx**3*self.g.sum())
        self.j = self.scale*self.g
        self.normal_log = -2*self.g*(1-1/self.u)
        self.c = self.normal_log/2 + 1j*self.j/self.n
        nbulk = abs(self.psi[1:])**2
        self.linear_potential = p.mu0*(nbulk**4-1)
        self.density_potential = 4*p.mu0*nbulk**3*self.psi[1:]

    def boundary(self, f, x):
        dn = 2*(self.face.conj()*f).real
        du = dn/(2*self.u)
        raw = self.g1*x + self.g*(dn+self.ref_grad*x)/self.n
        dj = self.scale*(raw-self.g*raw.sum()/self.g.sum())
        dc_real = -self.g1*x*(1-1/self.u)-self.g*du/self.n
        dc = dc_real+1j*(dj/self.n-self.j*dn/self.n**2)
        return self.c*f+self.face*dc, dj, 2*dc_real

    def generator(self, field):
        p = self.p
        lap = np.zeros_like(field)
        lap[1:-1] = (field[2:]-2*field[1:-1]+field[:-2])/p.dz**2
        normal = self.boundary(field[-1],0)[0]
        lap[-1] = 2*(field[-2]-field[-1]+p.dz*normal)/p.dz**2
        f = field[1:]
        for axis in (1,2,3):
            lap[1:] += (np.roll(f,1,axis)-2*f+np.roll(f,-1,axis))/p.dx**2
        dn = 2*(self.psi[1:].conj()*f).real
        lap[1:] = 1j*(lap[1:]/2-self.linear_potential*f-self.density_potential*dn)
        return lap

    def measurements(self, f, x):
        p = self.p
        dn = 2*(self.face.conj()*f).real
        du = dn/(2*self.u)
        _, dj, dlog = self.boundary(f,x)
        grad = p.x_derivative
        vx = (self.face.conj()*grad(self.face)).imag/self.n
        dvx = ((f.conj()*grad(self.face)+self.face.conj()*grad(f)).imag-vx*dn)/self.n
        carried = p.dx**3*(dj*vx+self.j*dvx)[p.mouth].sum()
        stress_terms = -p.dx**3/4*(dn*grad(self.normal_log)+self.n*grad(dlog))
        stress = stress_terms[p.mouth].sum()
        fixed = -p.dx**3/2*(self.g2*(self.u-1)**2).sum()*x
        response = -p.dx**3*(self.g1*(self.u-1)*du).sum()
        return np.array([carried,stress,fixed,response,stress_terms.sum()])


def main():
    _, drives = load_drives()
    times = drives["positive"][1]
    feed = MovingFeed()
    model = FeedTangent(feed)
    # For X=A*t^2/2, the tangent per unit A is
    # sum L^j B t^(j+3)/(j+3)! starting from zero, B the unit shift source.
    term = np.zeros_like(feed.base)
    term[-1] = 1j*model.boundary(np.zeros_like(model.face),1)[0]/model.p.dz * times[-1]**3/6
    face = np.zeros((len(times),*model.face.shape),dtype=complex)
    for order in range(3,62):
        face += (times/times[-1])[:,None,None,None]**order * term[-1]
        if order == 42:
            coarse = face.copy()
        if order < 61:
            term = model.generator(term)*times[-1]/(order+1)
    readings = np.array([model.measurements(f,t*t/2) for f,t in zip(face,times)])
    shorter = np.array([model.measurements(f,t*t/2) for f,t in zip(coarse,times)])
    init = model.measurements(np.zeros_like(model.face),1)
    print("initial force per displacement (carried, mask stress, virtual work) =",init[:3],flush=True)
    dt_error = np.max(abs(readings-shorter),axis=0)
    print("series truncation change per acceleration =",dt_error,flush=True)
    for i,t in enumerate(times):
        carried,stress,fixed,response,whole = readings[i]
        print(json.dumps(dict(t=t,A=0.02,carried=0.02*carried,mask_stress=0.02*stress,
                              frozen_profile_force=0.02*fixed,field_response_force=0.02*response,
                              profile_force=0.02*(fixed+response),whole_plane_stress=0.02*whole)),flush=True)
    # The nonlinear archived outputs independently test the tangent's two face currents.
    def observed(name):
        item = json.loads((RUNS/f"mouth_prescribed_accel_{name}02/stdout").read_text())
        return np.array([[r["carried"],r["stress"]] for r in item["32"]])
    for suffix in ("","_half"):
        name = "positive"+suffix
        amp = drives[name][0]
        odd = (observed(name)-observed("negative"+suffix))/2
        predicted = readings[:,:2]*amp
        error = abs(odd-predicted).max(axis=0)
        print(f"archive A={amp:g} peak absolute error (carried, mask stress) = {error}",flush=True)
        np.testing.assert_allclose(predicted,odd,rtol=3e-4,atol=4e-18)
    np.testing.assert_allclose(readings,shorter,rtol=1e-7,atol=1e-18)


if __name__ == "__main__":
    main()
