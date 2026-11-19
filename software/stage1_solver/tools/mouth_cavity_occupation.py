"""Which bound standing modes and occupations can two cavity loads distinguish?

Compare the archived sharp spherical side with a finite radial step of the
same inside/outside gaps. Both auxiliary cylinders retain the mouth node and
the reflecting end in every transverse section. Apply the stricter frequency
cut for an unconfined receiving medium separately from radial binding.
"""
from functools import lru_cache
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import spherical_jn, spherical_kn


ARCHIVE = Path(__file__).resolve().parents[1]/"reports/pathA_26_derrick_results.yaml"
S = yaml.safe_load(ARCHIVE.read_text())
R0, L0 = S["phase_A"]["q_star"]
PARAMS = S["parameters"]
C, MI, MO = (PARAMS[key] for key in ("c_w","mu_in","mu_out"))


@lru_cache(None)
def zero(l,n):
    if l == 0:
        return n*np.pi
    # Positive zeros interlace between successive orders.
    return brentq(lambda x: spherical_jn(l,x),zero(l-1,n),zero(l-1,n+1),xtol=1e-14)


def radial(l,n,r,soft):
    hard = zero(l,n)
    if not soft:
        return (C*hard/r)**2, 2*C*C*hard**2/r**3
    start = (n-.5)*np.pi if l==0 else zero(l-1,n)
    depth = MO**2-MI**2
    V = r*np.sqrt(depth)/C
    if V <= start:
        return None
    def mismatch(x):
        y = np.sqrt(max(0.,(V-x)*(V+x)))
        slope = -(l+1) if y < 1e-9 else y*spherical_kn(l,y,True)/spherical_kn(l,y)
        return x*spherical_jn(l,x,True)-slope*spherical_jn(l,x)
    x = brentq(mismatch,start,min(hard,V),xtol=1e-14)
    y = np.sqrt((V-x)*(V+x))
    j = spherical_jn(l,x)
    kn = spherical_kn(l,y)
    inside = quad(lambda s: (s*spherical_jn(l,x*s)/j)**2,0,1)[0]
    outside = quad(lambda s: (s*spherical_kn(l,y*s)/kn)**2,1,np.inf)[0]
    # First variation of the finite radial step at fixed geometry and state.
    minus_squared_slope = depth/(r*(inside+outside))
    return (C*x/r)**2, minus_squared_slope


def modes(r,L,soft):
    found=[]
    for l in range(24):
        if l and zero(l-1,1)>=r*np.sqrt(MO**2-MI**2)/C:
            break
        for n in range(1,24):
            ray = radial(l,n,r,soft)
            if ray is None:
                break
            transverse,radial_slope = ray
            if transverse+MI**2>=MO**2:
                break
            for j in range(24):
                kz = C*np.pi*(j+.5)/L
                w = np.sqrt(MI**2+transverse+kz*kz)
                if w>=MO:
                    break
                found.append(((l,n,j),w,radial_slope/(2*w),kz*kz/(w*L)))
    return sorted(found,key=lambda row:row[1])


def frequency(l,n,j,r,L,soft):
    ray = radial(l,n,r,soft)
    return np.sqrt(MI**2+ray[0]+(C*np.pi*(j+.5)/L)**2)


def threshold(soft):
    def det(r):
        kz=C*np.pi/(2*L0)
        x=r*np.sqrt(MO**2-MI**2-kz**2)/C
        if not soft:
            return x-zero(1,1)
        y=r*kz/C
        return x*spherical_jn(1,x,True)-y*spherical_kn(1,y,True)/spherical_kn(1,y)*spherical_jn(1,x)
    return brentq(det,1.35,2.1)


def run():
    print(f"reference R={R0:.12g} L={L0:.12g} c={C:g} gaps={MI:g},{MO:g}")
    for soft in (False,True):
        print(f"side={'finite' if soft else 'hard'} first_dipole_crosses_bulk_cut_R={threshold(soft):.12g}")
        for f in (.9,1.,1.125,1.25,1.5,2.):
            r = f*R0
            rows = modes(r,L0,soft)
            print(f"levels soft={soft} R/R0={f:g} groups=" +
                  ",".join(f"{l}/{n}/{j}@{w:.8f}" for (l,n,j),w,_,_ in rows))
            assert all(w<MO for key,w,b,c in rows)
        r=(1. if soft else 1.125)*R0
        pair=modes(r,L0,soft)
        assert [a[0] for a in pair]==[(0,1,0),(1,1,0)]
        loads=np.array([a[2:] for a in pair]).T
        checks=[]
        for key,w,fr,fl in pair:
            h=1e-5
            checks.extend([abs(fr+(frequency(*key,r+h,L0,soft)-frequency(*key,r-h,L0,soft))/(2*h)),
                           abs(fl+(frequency(*key,r,L0+h,soft)-frequency(*key,r,L0-h,soft))/(2*h))])
        # Two independent synthetic wave force channels. Actions sum over m.
        truth=np.array([.7,.3])
        measured=loads@truth
        fitted=np.linalg.solve(loads,measured)
        print(f"loads soft={soft} R={r:.10f} radial,axial={measured.tolist()} "
              f"action_fit={fitted.tolist()} det={np.linalg.det(loads):.8g} "
              f"max_force_derivative_error={max(checks):.3g}")
        assert max(checks)<1e-8 and np.max(abs(fitted-truth))<1e-11
        # Complex amplitudes from two orthogonal l=1 traction patterns;
        # amplitudes normalized to action, with the uniform component removed.
        for label,xy in (('cw',[1,1j]),('ccw',[1,-1j]),('linear',[np.sqrt(2),0])):
            ax,ay=np.sqrt(truth[1]/2)*np.asarray(xy)
            signed=2*np.imag(np.conj(ax)*ay)
            print(f"quadrature soft={soft} state={label} action_xy={abs(ax)**2+abs(ay)**2:.8g} "
                  f"twice_imag_axstar_ay={signed:.8g}")
        # At doubled radius, two forces can admit distinct nonnegative mixtures.
        first=modes(2*R0,L0,soft)[:3]
        H=np.array([m[2:] for m in first]).T
        null=np.cross(H[0],H[1])
        null/=np.max(abs(null))
        a=.5+null*.4
        b=.5-null*.4
        print(f"unresolved soft={soft} keys={[row[0] for row in first]} "
              f"action_a={a.tolist()} action_b={b.tolist()} load_a={(H@a).tolist()} "
              f"load_b={(H@b).tolist()}")
        assert min(np.r_[a,b])>0 and np.max(abs(H@(a-b)))<1e-12


if __name__=="__main__":
    run()
