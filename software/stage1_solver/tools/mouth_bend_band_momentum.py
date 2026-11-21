"""What do two wall-momentum bands on the draining arm tell us about the bend?"""
import json

import numpy as np
from scipy.integrate import quad

from mouth_bend_depth_participation import FOLDER, spectral_readout
from mouth_projected_bend_spectrum import fields


ENTRANCE = (0., .25)
MIDDLE = (.375, .625)


def f(x):
    return np.sin(np.pi*x)**2


def integrate(func, band, points=()):
    lo,hi = band
    return quad(func,lo,hi,points=[x for x in points if lo<x<hi],
                epsabs=1e-12,limit=200)[0]


def sharp_bounds(e,c):
    # e,c are the dimensionless first f moments in their respective windows.
    # B is the maximum in E and a the minimum in C (supremum outside E,C).
    b = f(ENTRANCE[1])
    a = f(MIDDLE[0])
    assert e>0 and c>0 and e/b+c<1
    mc = min(c/a,1-e/b)
    # Jensen in each band. The minimum assigns as much mass to C as can
    # reduce its mean f, then puts the rest in E at its required mean.
    lower = c*c/mc+e*e/(1-mc)
    # For the supremum, concentrate at f=b in E, at 1 in C and at a just
    # outside C. Sharp endpoint bounds can be approached by smooth layers.
    upper = b*e+c+a*a*(1-e/b-c)
    return lower,upper


def run():
    prior = json.loads((FOLDER/"mouth_bend_coefficients.json").read_text())["runs"]
    geo = fields(prior["bend_projection_wave015"]["stdout_lines"][0])
    radial = float(geo["breathing"])
    half = float(geo["half"])
    fold = float(fields(prior["bend_projection_fold64"]["stdout_lines"][0])["wall_mass_at_fold"])
    depth = json.loads((FOLDER/"mouth_bend_depth_participation.json").read_text())["runs"]
    lines = depth["bend_depth_pressure015"]["stdout_lines"]
    shares = np.array([float(fields(s)["radial_mass_share"]) for s in lines if s.startswith("attached ")])
    edges = np.minimum(np.arange(len(shares)+1)*.015/half,1.)
    assert edges[-1]==1 and np.isclose(shares.sum(),1)
    density = shares/shares.sum()/np.diff(edges)
    def pressure(x):
        return density[min(len(density)-1,np.searchsorted(edges,x,side="right")-1)]
    central_norm = integrate(lambda x:.02+f(x)**4,(0,1))
    def central(x):
        return (.02+f(x)**4)/central_norm
    def report(name,profile,points):
        e,c = [integrate(lambda x:profile(x)*f(x),band,points) for band in (ENTRANCE,MIDDLE)]
        q = integrate(lambda x:profile(x)*f(x)**2,(0,1),points)
        assert np.isclose(integrate(profile,(0,1),points),1)
        lo,hi = sharp_bounds(e,c)
        assert lo<=q<=hi
        print(f"profile={name} one_arm_entrance_P_over_bdot={radial*e/6:.12g} "
              f"one_arm_middle_P_over_bdot={radial*c/6:.12g} middle_over_entrance={c/e:.12g} "
              f"both_arm_bend_M={radial*q/3:.12g} two_band_M_min={radial*lo/3:.12g} "
              f"two_band_M_sup={radial*hi/3:.12g}",flush=True)
        return e,c,radial*q/3
    e,c,mp = report("pressure",pressure,edges)
    report("centrally_weighted",central,[])
    # Smooth normalized bumps within the sampled bands and in an unobserved
    # shoulder. Match both pressure-law readings exactly by construction.
    width = 1/64
    intervals = [(ENTRANCE[1]-width,ENTRANCE[1]),(.5-width/2,.5+width/2),
                 (MIDDLE[1],MIDDLE[1]+width)]
    def bump(x,span):
        lo,hi = span
        return 2*np.sin(np.pi*(x-lo)/width)**2/width if lo<x<hi else 0.
    weights = [moment/integrate(lambda x:bump(x,span)*f(x),span)
               for moment,span in zip((e,c),intervals[:2])]
    weights.append(1-sum(weights))
    assert min(weights)>0
    def matched(x):
        patch = sum(w*bump(x,span) for w,span in zip(weights,intervals))
        return .05*pressure(x)+.95*patch
    points = sorted(set(edges)|{x for span in intervals for x in span})
    e2,c2,m2 = report("matched_readings",matched,points)
    assert max(abs(e2-e),abs(c2-c))<1e-11
    print(f"fixed_stiffness_fold={fold:.12g} one_arm_weighted_P_over_bdot_at_fold={fold/2:.12g} "
          f"constructed_patch_weights={','.join(f'{w:.10g}' for w in weights)}",flush=True)
    spectral_readout(mp,"pressure_bands")
    spectral_readout(m2,"matched_bands")


if __name__=="__main__":
    run()
