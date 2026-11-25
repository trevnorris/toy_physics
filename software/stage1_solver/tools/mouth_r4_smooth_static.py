"""Frozen old-to-smooth R4 gauge work, shape variations and orientation controls.

Order is Q1-interpolated before the smoothstep. The harmonic part uses the
two-point physical-volume rule. The barrier retains the eight corner gradients
of each cell, and its mean is integrated with the same weighted volume.
"""
import argparse
import json
import sys
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from mouth_r4_projected_fields import EXPECTED, MANIFEST, digest, read

BITS = np.array(list(product((0,1), repeat=3)))
RULE = tuple(product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)), repeat=3))
LD = np.longdouble


class Gauge:
    def __init__(self, trial, reference, order):
        from mouth_r4_eulerian_release import block_axes
        self.axes = block_axes(trial.b)
        self.shape, self.polar, self.reference = trial.shape, trial.polar_shift, reference
        self.cells = np.indices(tuple(n-1 for n in order.shape)).reshape(3,-1).T
        self.ids = np.array([np.ravel_multi_index(tuple((self.cells+b).T),order.shape) for b in BITS]).T
        self.order = order.ravel()[self.ids]
        self.base = np.stack(np.meshgrid(*self.axes,indexing='ij'),axis=-1)

    def cell_barrier(self, offsets, widths, gradient):
        """Cell mean and optional nodal derivative on the transported chart."""
        from mouth_cp10_label_gauge import gradient_basis
        dtype = offsets.dtype
        barrier = np.zeros(len(offsets),dtype=dtype)
        barrier_grad = np.zeros_like(offsets) if gradient else None
        identity = np.eye(3,dtype=dtype)
        for vertex in BITS:
            basis = gradient_basis(vertex).astype(dtype)
            D = (offsets.transpose(0,2,1) @ basis)/widths[:,None,:]+identity
            col = D.transpose(0,2,1)
            cofactor = np.stack([np.cross(col[:,(k+1)%3],col[:,(k+2)%3]) for k in range(3)],axis=2)
            det = np.sum(D[:,:,0]*cofactor[:,:,0],axis=1)
            if np.min(det.real)<=0:
                raise ValueError('Nonpositive corner Jacobian')
            barrier += (det-1-np.log(det))/8
            if gradient:
                barrier_grad += ((1-1/det)/8)[:,None,None]*(basis @ (cofactor/widths[:,None,:]).transpose(0,2,1))
        return barrier, barrier_grad

    def evaluate(self, labels, dq, *, enabled=None, gradient=False):
        from mouth_cp10_label_gauge import depth_grid, gradient_basis, profile
        # Promote the determinants, volume, field and signed sums upstream.
        dtype = np.result_type(dq,LD)
        x,t,a = [np.asarray(axis,dtype=dtype) for axis in self.axes]
        length = LD(self.shape.length)*np.exp(dq[0])
        moved = depth_grid(length,self.reference).astype(dtype)
        cells = self.cells if enabled is None else self.cells[enabled]
        ids = self.ids if enabled is None else self.ids[enabled]
        order = self.order if enabled is None else self.order[enabled]
        axes = (x,moved,a)
        lo = np.stack([ax[cells[:,k]] for k,ax in enumerate(axes)],axis=1)
        widths = np.stack([np.diff(ax)[cells[:,k]] for k,ax in enumerate(axes)],axis=1)
        # Shape motion fixes xi minus the normalized background coordinate.
        off = labels.astype(dtype).reshape(-1,3)[ids]-self.base.reshape(-1,3)[ids].astype(dtype)
        offsets = off-off[:,:1]
        barrier, barrier_grad = self.cell_barrier(offsets, widths, gradient)
        energy = np.zeros(len(cells),dtype=dtype)
        grad = np.zeros_like(off) if gradient else None
        for q in RULE:
            X,T,A = (lo+widths*q).T
            R,Rt=profile(self.shape,T,dq)
            y=(T-.06)/.56
            inside=(y.real>0)&(y.real<1)
            bump=np.where(inside,16*y*y*(1-y)**2,0)
            slope=np.where(inside,32*y*(1-y)*(1-2*y)/.56,0)
            polar=self.polar+dq[5]
            P=(3*np.cos(A)**2-1)/2
            r=X*R+X**3*polar*bump*P
            rx=R+3*X*X*polar*bump*P
            rt=X*Rt+X**3*polar*slope*P
            ra=-3*X**3*polar*bump*np.cos(A)*np.sin(A)
            vol=.5*LD(np.pi)*widths.prod(axis=1)*rx*r*r*length*np.sin(A)
            N=np.prod(np.where(BITS,q,1-np.array(q)),axis=1)
            chi=order.astype(dtype)@N.astype(dtype)
            s=(1-chi)**2*(1+2*chi)
            basis=gradient_basis(q).astype(dtype)
            D=(offsets.transpose(0,2,1) @ basis)/widths[:,None,:]
            gr=D[:,:,0]/rx[:,None]
            gw=(D[:,:,1]-gr*rt[:,None])/length
            ga=(D[:,:,2]-gr*ra[:,None])/r[:,None]
            w=s*vol
            energy += w*(barrier+.5*np.sum(gr*gr+gw*gw+ga*ga,axis=1))
            if gradient:
                c1=gw/length
                c2=ga/r[:,None]
                c0=(gr-rt[:,None]*c1-ra[:,None]*c2)/rx[:,None]
                dual=np.stack([c0,c1,c2],axis=-1)/widths[:,None,:]
                grad += w[:,None,None]*(basis@dual.transpose(0,2,1)+barrier_grad)
        if not gradient:
            return energy
        assembled=np.zeros(labels.shape,dtype=dtype)
        np.add.at(assembled.reshape(-1,3),ids.ravel(),grad.reshape(-1,3))
        return energy,assembled


def frozen(area, tmp):
    pin = area/'pin/software/stage1_solver/tools'
    sys.path.insert(0,str(pin))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_r4_eulerian_release import Driver,order_grid
    from mouth_r4_stage2_outcome import verify_guard
    if digest(area/'sources.json')!=MANIFEST or any(digest(pin/k)!=v for k,v in read(area/'sources.json').items()):
        raise ValueError('Held source manifest differs')
    verify_guard(area/'primary_guard',1)
    state=area/'primary/held_stop.npz'
    if digest(state)!=EXPECTED['held_stop']:
        raise ValueError('Wrong frozen state')
    with np.load(state,allow_pickle=False) as f:
        fields, labels, x = (f[k] for k in ('field','labels','x'))
        if (int(f['radius']),int(f['winding']))!=(4,1):
            raise ValueError('Changed preparation')
    d=Driver(area.parent,tmp,1e-8)
    if not np.array_equal(x,d.x):
        raise ValueError('Mouth shape changed')
    trial=ShapeTrial(d.base,d.prep,x,8,24,6,3,2,precondition=False,
                     axes=axes_for_depth(d.shape(x).length,d.reference_depth))
    return d, trial, fields, labels, order_grid(trial,fields)


def negative_controls(opt,labels,order):
    from mouth_r4_eulerian_relax_labels import InfeasibleLabels, chart_orientation
    from mouth_cp10_label_gauge import gradient_basis
    cells = np.indices(opt.active.shape).reshape(3,-1).T
    vertices=cells[:,None,:]+BITS
    values=order[tuple(vertices.transpose(2,0,1))]
    out={}
    # Swap the opposing radial faces of one target cell. Its determinant
    # reverses sign throughout the cell; all its labels stay in the domain.
    choose=(cells[:,0]>1)&(cells[:,1]>0)&(cells[:,2]>0)
    q = np.array([.5-.5/np.sqrt(3)]*3)
    N=np.prod(np.where(BITS,q,1-q),axis=1)
    qchi=values@N
    for key,mask in [('ordered',np.min(values,axis=1)>0),('fringe',(np.max(values,axis=1)>0)&(np.min(values,axis=1)==0))]:
        relevant=np.flatnonzero(choose&mask)
        idx=relevant[np.argmin(qchi[relevant]) if key=='fringe' else np.argmax(qchi[relevant])]
        ix=vertices[idx]
        coords=tuple(ix.T)
        trial=labels.copy()
        trial[coords]=labels[coords][np.arange(8)^4]
        width=np.array([np.diff(ax)[i] for ax,i in zip(opt.axes,cells[idx])])
        def determinant(array):
            local=array[coords]
            return np.linalg.det(((local-local[0]).T@gradient_basis(q))/width)
        before,after=determinant(labels),determinant(trial)
        try:
            chart_orientation(opt.axes,trial,opt.active)
        except InfeasibleLabels as error:
            if not (before>0 and after<0 and qchi[idx]>0):
                raise ArithmeticError('Orientation control missed target')
            out[key]=dict(cell=cells[idx].tolist(),gauss_fraction=q.tolist(),order=float(qchi[idx]),
                          before=float(before),after=float(after),rejection=error.witness)
        else:
            raise ArithmeticError('Guard admitted an inverted cell')
    return out


def run(area,record):
    with TemporaryDirectory() as directory:
        tmp=Path(directory)
        d,t,fields,labels,order=frozen(area,tmp)
        from mouth_r4_eulerian_relax_labels import LabelRelaxation
        from mouth_r4_eulerian_solid import SpatialCrystal
        from mouth_r4_melt_barrier import MeltBarrier
        from mouth_r4_barrier_shape import unit_work
        from mouth_r4_slaved_melt import FLOW
        crystal=SpatialCrystal(t.b)
        elastic,forces,source,_=crystal.evaluate(labels,order)
        inner=np.bincount(t.b.ids.ravel(),source.transpose(0,2,1).ravel(),minlength=t.model.size)
        currents=t.measure(fields,t.ext+inner)
        opt=LabelRelaxation(t.b,labels,order,1e-8)
        old=MeltBarrier(t.b,opt,tmp)
        old_e,old_g,_,_=old.sample(labels)
        smooth=Gauge(t,d.reference_depth,order)
        zero=np.zeros(6,LD)
        val,grad=smooth.evaluate(labels,zero,gradient=True)
        new_e=val.sum()
        difference=old.extract(np.asarray(grad,float)-old_g,'free')
        physical=old.extract(forces,'free')
        old_shape=[];new_shape=[];fd=[]
        for offset in np.eye(6,dtype=LD):
            dq=offset.astype(np.clongdouble)*np.clongdouble(1.e-24j)
            new_shape.append(float(-smooth.evaluate(labels,dq).sum().imag/1.e-24))
            old_shape.append(float(-unit_work(t.shape,t.polar_shift,d.reference_depth,opt.axes,labels,
                                             ~opt.active,np.asarray(dq,np.complex128)).imag/1.e-24))
            banks=[]
            for h in [LD(2)**-18,LD(2)**-20]:
                minus=smooth.evaluate(labels,-h*offset)
                plus=smooth.evaluate(labels,h*offset)
                banks.append(float(np.sum(minus-plus)/(2*h)))
            fd.append(banks)
        signed=np.asarray(new_shape)-np.asarray(old_shape)
        control=negative_controls(opt,labels,order)
        previous=read(area/'primary/held_stop.json')['current']
        flow_error={k:float(currents[k]-previous[k]) for k in FLOW}
        if any(flow_error.values()):
            raise ArithmeticError(f'Physical current on fixed state changed: {flow_error}')
        result=dict(source_sha256=digest(Path(__file__)), state_sha256=EXPECTED['held_stop'],
                    command=sys.orig_argv, gauge_quadrature='2 point Gauss s(Q1 chi), mean 8 corner log barriers',
                    unit_energy=dict(boolean=float(old_e),smooth=float(new_e)),
                    frozen_physical_crystal_energy=elastic,max_abs_physical_label_derivative=float(max(abs(physical))),
                    new_unit_six_shape_forces=new_shape,old_unit_six_shape_forces=old_shape,
                    new_unit_six_shape_forces_extended_fd=fd,negative_controls=control,
                    direct_current_readings={k:float(currents[k]) for k in FLOW},current_readback_delta=flow_error,
                    static_at_epsilon={})
        for eps in (0.,1e-10,1e-8):
            result['static_at_epsilon'][str(eps)]=dict(delta_physical_crystal_energy=0.,
                  delta_gauge_energy=float(eps*(new_e-old_e)),max_abs_free_active_force_change=float(eps*max(abs(difference))),
                  six_shape_force_changes=(eps*signed).tolist(),physical_label_selector_unchanged=True,
                  fixed_field_currents_delta=flow_error)
        record.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(json.dumps({k:result[k] for k in ('unit_energy','static_at_epsilon')}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('area',type=Path)
    parser.add_argument('record',type=Path)
    args=parser.parse_args()
    run(args.area.resolve(),args.record.resolve())
