"""Linear coherent current through one joined, neutral cell with two uniform phase plates.

Space is (x,y,z,w) with rotational symmetry in the first three directions.
The reference carrier is at rest. A hypothetical order cap occupies r < a;
the surrounding flat face is internal to the carrier. No reflecting wall or
separate carrier phase is imposed on that face. The material mobility affects
the order equation only. This tests a necessary stationary phase condition.
"""
import argparse

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from mouth_flat_phase_control import saved_bath
from mouth_free_front_wave_head import emit


def measure(radial_cells, layered, phase):
    length, radius = 32., 16.
    nz = 2*radial_cells
    dz, dr = length/nz, radius/radial_cells
    r = np.linspace(0., radius, radial_cells+1)
    z = -length/2 + (np.arange(nz)+.5)*dz
    area = 4*np.pi*np.diff(r**3)/3
    # The nonuniform option is a prescribed positive linear-response coefficient.
    # Its sustaining static potential is outside the carrier phase test.
    density = 1./(1.+.4*np.cos(2*np.pi*z/length)) if layered else np.ones(nz)
    ids = np.arange(radial_cells*nz).reshape(radial_cells,nz)
    z_conductance = area[:,None]/(.5*dz*(1/density[:-1]+1/density[1:]))
    r_conductance = 4*np.pi*r[1:-1,None]**2*dz*density[None,:]/dr
    left = np.r_[ids[:,:-1].ravel(),ids[:-1,:].ravel()]
    right = np.r_[ids[:,1:].ravel(),ids[1:,:].ravel()]
    weight = np.r_[z_conductance.ravel(), r_conductance.ravel()]
    plates = np.r_[ids[:,0],ids[:,-1]]
    boundary_weight = np.r_[area*density[0],area*density[-1]]/(dz/2)
    rows = np.r_[left,right,left,right,plates]
    cols = np.r_[left,right,right,left,plates]
    values = np.r_[weight,weight,-weight,-weight,boundary_weight]
    matrix = coo_matrix((values,(rows,cols)),shape=(ids.size,ids.size)).tocsr()
    rhs = np.zeros(ids.size)
    rhs[ids[:,-1]] = phase*boundary_weight[radial_cells:]
    theta = spsolve(matrix,rhs).reshape(ids.shape)
    gap = nz//2-1
    current = z_conductance[:,gap]*(theta[:,gap+1]-theta[:,gap])
    face = r[:-1] >= radius/2
    expected = phase/(dz*np.sum(1./density))
    feed = area*2*density[0]*theta[:,0]/dz
    bulk = area*2*density[-1]*(phase-theta[:,-1])/dz
    emit(test='neutral_connected_cell', nr=radial_cells, nw=nz,
         density_law='reciprocal_cosine' if layered else 'uniform', delta_phase=phase,
         cap_projected_radius=radius/2, outer_radius=radius, plate_distance=length,
         predicted_current_per_projected_area=expected, mean_normal_current_on_flat=sum(current[face])/sum(area[face]),
         central_projected_current=sum(current[~face]), flat_annular_current=sum(current[face]),
         flat_share=sum(current[face])/sum(current), max_current_density_error=max(abs(current/area-expected)),
         feeding_flux=sum(feed), bulk_flux=sum(bulk), half_level_speed_if_locked_and_flat=expected/(1/(1+.4) if layered else 1.))
    if np.max(abs(current/area-expected))>2e-10 or max(abs(sum(current)-sum(feed)),abs(sum(current)-sum(bulk)))>2e-8:
        raise ArithmeticError('Axial flux or global number balance failed')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=(24,48),required=True)
    args=parser.parse_args()
    for nonuniform, multiplier in ((False,1.),(True,1.),(False,-.5)):
        measure(args.mesh,nonuniform,multiplier*saved_bath())
