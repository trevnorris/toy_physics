"""Read neighbor torque and local turning after the archived thermal strain hold."""
import json
from pathlib import Path

import numpy as np

from neutral_patch_resting_order_hold import digest, emit
from neutral_patch_sites import COUNT, GAP, LAWS, PITCH_UNIT, Patch, morse


RECORD = Path('research/pde_ledger_v3/_measurements/neutral_patch_strain_hold.json')


def cross(a, b):
    return a[...,0]*b[...,1]-a[...,1]*b[...,0]


def turn(before, after):
    a = before-before.mean(axis=-2, keepdims=True)
    b = after-after.mean(axis=-2, keepdims=True)
    return np.arctan2(cross(a,b).sum(axis=-1), (a*b).sum(axis=(-2,-1)))


def gauge(points, patch, selected, order):
    """Moment of the outside sites' pair forces on their inside neighbors."""
    p = patch
    delta = points[...,p.j,:]-points[...,p.i,:]
    s = (delta*delta).sum(axis=-1)
    r2 = s.copy()
    geometric = np.ones_like(s)
    if p.law=='layered':
        r2[...,p.cross] = s[...,p.cross]**2/PITCH_UNIT**2
        geometric[...,p.cross] = 2*s[...,p.cross]/PITCH_UNIT**2
    r2[...,p.cross] += GAP**2
    r = np.sqrt(r2)
    well = np.where(p.cross,p.depth*order/.95,p.depth)
    energy, slope = morse(r, well, p.alpha, p.length)
    force_i = (geometric*slope/r)[...,None]*delta
    # Positive x force exerted on the lower plane, checked against the saved gauge.
    lower = force_i[...,p.cross,0].sum(axis=-1)
    cut = selected[p.i]!=selected[p.j]
    inside = np.where(selected[p.i], p.i, p.j)[cut]
    sign = np.where(selected[p.i[cut]],1.,-1.)
    center = points[...,selected,:].mean(axis=-2,keepdims=True)
    moment = cross(points[...,inside,:]-center, force_i[...,cut,:]*sign[:,None]).sum(axis=-1)
    return moment,lower,energy[...,cut].sum(axis=-1)


def compare(law):
    path = Path('_scratch/s11c/strain_hold')/f'hold_{law}.npz'
    pin = json.loads(RECORD.read_text())['primary_output_sha256'][str(path)]
    if digest(path)!=pin:
        raise ValueError('Changed thermal end-state input')
    p = Patch(law)
    with np.load(path) as data:
        xy = data['final_xy']
        gate = data['gate_xy'][0,0].reshape(2,COUNT,2)
        dist = ((gate-gate.mean(axis=1,keepdims=True))**2).sum(axis=2)
        np.testing.assert_allclose(data['amplitudes'], [0.,.001,-.001,.0005,-.0005])
        for count in (4,6):
            ids = np.sort((np.argsort(dist,axis=1)[:,:count]+np.array([[0],[COUNT]])).ravel())
            mask = np.isin(np.arange(2*COUNT),ids)
            for index, chi in enumerate(data['orders']):
                state = xy[:,index]  # step size, bath history, strain sign, site and coordinate
                torque,force,_ = gauge(state,p,mask,chi)
                reference_error = abs(force-data['force'][:,index,:,:,-1]).max()
                control = state[:,:,0:1]
                angle = turn(control[...,mask,:],state[...,mask,:])-turn(control[...,~mask,:],state[...,~mask,:])
                # A small virtual turn tests signed torque as -dU_cut/dangle.
                frame = state[-1,0,1]
                midpoint = frame[mask].mean(axis=0)
                centered = frame[mask]-midpoint
                energies = []
                for h in (1e-6,-1e-6):
                    test = frame.copy()
                    test[mask] = midpoint + centered*np.cos(h)+np.stack((-centered[:,1],centered[:,0]),axis=1)*np.sin(h)
                    energies.append(gauge(test,p,mask,chi)[2])
                work_error = abs((energies[0]-energies[1])/2e-6+torque[-1,0,1])
                row = {}
                for name, array in [('neighbor_torque',torque),('relative_angle',angle)]:
                    odd = (array[:,:,1]-array[:,:,2])/.002
                    half = (array[:,:,3]-array[:,:,4])/.001
                    row[name] = dict(mean=float(odd[-1].mean()),per_history=odd[-1].tolist(),
                                     half_mean=float(half[-1].mean()),coarse_mean=float(odd[0].mean()))
                t_odd = (torque[-1,:,1]-torque[-1,:,2])/.002
                a_odd = (angle[-1,:,1]-angle[-1,:,2])/.002
                emit('warm_endpoint',law=law,order=chi,inner_per_plane=count,ids=ids.tolist(),
                     time=float(data['time'][-1]),input_sha256=pin,torque_virtual_work_abs_error=work_error,
                     archived_layer_force_error=reference_error,opposite_increment_signs=int(np.sum(t_odd*a_odd<0)),**row)
                if reference_error>1e-10 or work_error>2e-7:
                    raise ArithmeticError('Warm force or signed virtual-torque check failed')


if __name__=='__main__':
    for law in LAWS:
        compare(law)
