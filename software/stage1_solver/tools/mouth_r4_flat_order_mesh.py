"""Signed axial mesh work for one unchanged, planar material interface.

Standard Q1 order and the two Gauss nodes used by the joined carrier. The
continuum reference depends on physical w only. A hypothetical terminal
section changes no domain volume and thus has no physical depth force here.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from mouth_breathing_mobile_patch import KAPPA, WELL, WIDTH
from mouth_exchange_carrier import mapped_nodes
from mouth_last_rings_fluid import GRIDS


def native_depths(length, refinement):
    nr, _, nb = GRIDS['double']
    feed = mapped_nodes(nr, length, 28., 22.)[1][0]
    near = .12
    far = (np.sinh(2*np.linspace(0., 1., nr+1))/np.sinh(2))**2
    blocks = (feed, np.linspace(0., length-near, nr+1), np.linspace(length-near, length, nb+1),
              np.linspace(length, length+near, nb+1), length+near+(22.-length-near)*far)
    return np.r_[np.concatenate([(b[:-1,None]+np.diff(b)[:,None]*np.arange(refinement)/refinement).ravel()
                                  for b in blocks]), blocks[-1][-1]]


def energy(length, refinement, radius=28.):
    w = native_depths(length, refinement)
    h = np.diff(w)
    chi = expit(-w/WIDTH)
    q = (1.+np.array([-1.,1.])/np.sqrt(3.))/2.
    sampled = chi[:-1,None]*(1-q)+chi[1:,None]*q
    gradient = np.sum(.5*KAPPA*np.diff(chi)**2/h)
    potential = h @ np.mean(WELL*(sampled*(1.-sampled))**2, axis=1)
    volume = 4*np.pi*radius**3/3
    return np.array([gradient, potential])*volume


def test(endpoint):
    data = json.loads(endpoint.read_text())
    length = data['depth']
    if data['pole_radius'][0] != 4.:
        raise ValueError('Expected held radius-four preparation')
    trials = []
    for refinement in (1,2,4):
        for step in (1.e-4, 5.e-5):
            a = energy(length*np.exp(-step), refinement)
            b = energy(length*np.exp(step), refinement)
            force = (a-b)/(2*step)
            trials.append(dict(axial_subdivision=refinement, log_depth_step=step,
                               gradient_force=force[0], double_well_force=force[1], total_force=force.sum()))
    return dict(depth=length, gd_terms_at_original=energy(length,1).tolist(), trials=trials)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coarse', type=Path)
    parser.add_argument('fine', type=Path)
    args=parser.parse_args()
    print(json.dumps(dict(scope='Nodal interpolation of the same neutral planar profile; carrier and lattice absent. '
                                'Uniformly subdivided native axial mesh. Exact radial/angular volume. '
                                'No-mouth continuum force on this change of mesh is zero.',
                          coarse=test(args.coarse),fine=test(args.fine)), allow_nan=False))
