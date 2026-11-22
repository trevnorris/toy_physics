"""Check the free-exit states with independent order sites and all six collar coordinates."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_exit_free_order import relax_order
from mouth_exit_witness import inputs
from mouth_poly_shallow_release import Trial, emit
from mouth_sealed_shape_relax import Meridian


def run(radius, mesh, origin, folder):
    trial = Trial(radius,mesh,bonded_core=True)
    shifts,chi,seed,digest,status,_,_ = inputs(origin,radius,mesh,trial)
    if status != 'balanced':
        raise ValueError('The four-coordinate state has no validated local balance')
    folder.mkdir(exist_ok=True,parents=True)
    shape0 = Meridian(radius,trial.x+shifts[:4])
    initial_shift = shifts[1]
    x = shifts.copy()
    x[1] = 1.

    def material_shifts(x):
        if x[1]<=0.:
            raise ValueError('Reached boundary of the positive-radius material chart')
        new = x.copy()
        new[1] = initial_shift + np.log(x[1])
        return new

    def fields(x,chi,seed):
        shifts = material_shifts(x)
        shape = Meridian(radius,trial.x+shifts[:4])
        if shape.wall()['minimum_radius']<=0:
            raise ValueError('Material chart reached the symmetry axis')
        data = trial.evaluate(shifts,chi,seed)
        info = data[4]
        if (not np.isfinite(data[0]) or info['current']<=0 or info['mode_overlap']<.94
                or info['pole_residual']>1e-7 or abs(info['current']-info['outflow'])>1e-6*info['current']):
            raise ArithmeticError('Unresolved wave or carrier balance on released material')
        info['exit_radius'] = shape.exit
        return data

    def checkpoint(x,chi,seed,name):
        temp = folder/(name+'.tmp.npz')
        file = folder/(name+'.npz')
        np.savez_compressed(temp, radius=radius,mesh=mesh,shifts=material_shifts(x),
                            chi=chi,pole=seed[0],vector=seed[1],parent_sha256=digest)
        temp.replace(file)
        return hashlib.sha256(file.read_bytes()).hexdigest()

    def slopes(x,chi,seed,h):
        grad = []
        for i in range(6):
            dx = np.zeros(6)
            dx[i] = min(.05*x[1],h) if i==1 else h
            ep = fields(x+dx,chi,seed)[0]
            em = fields(x-dx,chi,seed)[0]
            grad.append((ep-em)/(2*dx[i]))
        return np.array(grad)

    def actuate_exit(x,chi,seed,base):
        for h in (.001,.0005):
            dx = np.zeros(6)
            dx[1] = min(h,.05*x[1])
            a = fields(x-dx,chi,seed)
            b = fields(x+dx,chi,seed)
            emit(kind='independent_exit',step=dx[1],minus_work=a[0]-base,plus_work=b[0]-base,
                 exit_energy_slope=(b[0]-a[0])/(2*dx[1]))

    initial = fields(x,chi,seed)
    initial_energy = initial[0]
    emit(kind='from_exit_balance',radius=radius,mesh=mesh,parent_sha256=digest,**initial[4])
    metric = np.diag([1.,1.,.15,.1,1.,.3])/(16.*radius**2)
    previous = None
    iteration = 0
    while True:
        base,_,owner,loads,base_info,elastic = fields(x,chi,seed)
        fixed_part = base_info['energy']['outer_order']+base_info['energy']['static_outside']
        def progress(chi,seed,**readout):
            checkpoint(x,chi,seed,'latest')
            emit(kind='order_progress',geometry_iteration=iteration,**readout)
        chi,order_info = relax_order(chi,owner,loads,elastic,progress)
        seed = loads.wave.seed
        state = fields(x,chi,seed)
        value,seed,info = state[0],state[3].wave.seed,state[4]
        if (value>base+3e-7 or abs(value-fixed_part-order_info['inside_energy'])>3e-7):
            raise ArithmeticError('Free order differs from the independently reassembled total energy')
        emit(kind='order_state',geometry_iteration=iteration,change=value-initial_energy,
             stationary=order_info['stationary'],reason=order_info['reason'],**info)
        if not order_info['stationary']:
            status = 'order_stationarity_unresolved'
            break
        grad = slopes(x,chi,seed,.0004)
        emit(kind='all_shape_forces',geometry_iteration=iteration,gradient=grad.tolist(),
             coordinates=x.tolist(),change=value-initial_energy)
        if max(abs(grad))<.001:
            finer = slopes(x,chi,seed,.0002)
            status = 'stationary_in_represented_coordinates' if max(abs(finer))<.002 else 'shape_precision_unresolved'
            break
        if previous is not None:
            old_x,old_grad = previous
            dx = x-old_x
            dg = grad-old_grad
            curvature = dg@dx
            if curvature>1e-10*np.linalg.norm(dx)*np.linalg.norm(dg):
                left = np.eye(6)-np.outer(dx,dg)/curvature
                metric = left@metric@left.T+np.outer(dx,dx)/curvature
        step = -metric@grad
        if step@grad>=0 or not np.isfinite(step).all():
            raise ArithmeticError('Invalid geometric search direction')
        step /= max(1.,np.max(abs(step))/.16,-step[1]/(.3*x[1]))
        factor = 1.
        reason = 'No resolved decrease along geometry direction'
        while max(abs(factor*step))>5e-7:
            try:
                proposed = fields(x+factor*step,chi,seed)
                if proposed[0]<=value+1e-4*factor*(step@grad):
                    break
            except (ValueError,ArithmeticError,RuntimeError) as err:
                reason = str(err)
            factor /= 2.
        else:
            emit(kind='unresolved_shape_step',reason=reason)
            status = 'geometric_step_unresolved'
            break
        previous = x.copy(),grad
        x = x+factor*step
        seed = proposed[3].wave.seed
        iteration += 1
    final_sha = checkpoint(x,chi,seed,'endpoint')
    emit(kind='material_endpoint',status=status,snapshot_sha256=final_sha,change=value-initial_energy,
         pole=[seed[0].real,seed[0].imag],**info)
    actuate_exit(x,chi,seed,value)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',choices=(1,2),type=int,required=True)
    p.add_argument('--mesh',choices=(48,72),type=int,required=True)
    p.add_argument('--origin',type=Path,required=True)
    p.add_argument('--folder',type=Path,required=True)
    a=p.parse_args()
    run(a.radius,a.mesh,a.origin,a.folder)
