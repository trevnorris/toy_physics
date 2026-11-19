"""Check the relaxed small-tilt stiffness and the continuous paths between fitted states."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular

from relaxed_neck import Neck


def check(path):
    model=Neck(*path["shape"],path["free_angle"],*path["wave"],path["quad"])
    states=path["states"]
    zero=states[0]
    a0=np.array(zero["coeff"])
    z0=np.r_[0.,model.pre.T@a0]
    def gradient(z):
        a=solve_triangular(model.pre.T,z[1:],lower=False)
        _,g=model.value(z[0],a)
        return np.r_[g[0],solve_triangular(model.pre,g[1:],lower=True)]
    for step in (0.002,0.001):
        hess=np.column_stack([(gradient(z0+d)-gradient(z0-d))/(2*step)
                              for d in step*np.eye(len(z0))])
        hess=(hess+hess.T)/2
        stiff=hess[0,0]-hess[0,1:]@np.linalg.solve(hess[1:,1:],hess[1:,0])
        print(f"free_angle={path['free_angle']} hessian_step={step} "
              f"relaxed_tilt_stiffness={stiff:+.8f} frozen_increment_stiffness={hess[0,0]:+.8f}",
              flush=True)
    # Check the derivative along the actual interpolant. Disconnected minima
    # could look like a barrier-free path if only optimized endpoints were read.
    greatest_slope=-np.inf
    excess=-np.inf
    for left,right in zip(states,states[1:]):
        lo,hi=left["beta"],right["beta"]
        a=np.array(left["coeff"]); slope=(np.array(right["coeff"])-a)/(hi-lo)
        for f in np.linspace(0,1,9):
            beta=lo+f*(hi-lo)
            e,g=model.value(beta,a+(beta-lo)*slope)
            if beta>0:
                excess=max(excess,e-zero["energy"])
                greatest_slope=max(greatest_slope,g[0]+g[1:]@slope)
    print(f"free_angle={path['free_angle']} largest_positive_tilt_excess={excess:+.9f} "
          f"greatest_path_slope_away_from_zero={greatest_slope:+.9f}",flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("paths",type=Path)
    for path in json.loads(ap.parse_args().paths.read_text()).values():
        check(path)
