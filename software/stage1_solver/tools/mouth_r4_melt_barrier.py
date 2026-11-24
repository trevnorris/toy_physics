"""Fixed-trace Eulerian melted-label barrier test on the two held R4 states.

The normalized gauge uses the prior physical-volume harmonic term, plus
d-1-log(d), with d=det(dxi/dp) relative to the identity chart. Each strictly
melted Q1 cell distributes its two-point physical four-volume equally across
its eight one-sided corner gradients. Whole-cell positivity is separately
certified. This gauge has no coupling to the registered reference metrics.
"""
import argparse
import ctypes as ct
import json
from itertools import product
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_label_condition import physical_frame
from mouth_cp10_label_grid import BITS, get_cp10
from mouth_exact_checkpoint_control import digest, read
from mouth_r4_eulerian_relax_labels import InfeasibleLabels, LabelRelaxation, absmax, chart_orientation
from mouth_r4_eulerian_solid import SpatialCrystal, InvalidSpatialCrystal
from mouth_whole_r4_relax import atomic_json


class MeltBarrier:
    def __init__(self, body, opt, output):
        self.opt=opt
        self.components=opt._components()
        self.sizes=[int(c.dead.sum()) for c in self.components]
        self.cuts=np.cumsum(self.sizes)[:-1]
        cells=np.argwhere(~opt.active)
        dims=opt.stencil.shape
        self.ids=np.ascontiguousarray(np.stack([np.ravel_multi_index(tuple((cells+b).T),dims)
                                               for b in BITS],axis=1),dtype=np.int32)
        self.width=np.ascontiguousarray(np.column_stack([np.diff(ax)[cells[:,i]] for i,ax in enumerate(opt.axes)]))
        lo=np.column_stack([ax[cells[:,i]] for i,ax in enumerate(opt.axes)])
        vol=np.zeros(len(cells))
        for q in product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)),repeat=3):
            basis,hoop=physical_frame(body,lo+self.width*q)
            vol+=4*np.pi*hoop*np.linalg.det(basis)*self.width.prod(axis=1)/8.
        self.weight=np.ascontiguousarray(vol/8.)
        if np.any(self.weight<=0):
            raise ValueError('A strictly melted cell has zero quadrature weight')
        source=Path(__file__).with_suffix('.cpp')
        library=output/'barrier.so'
        subprocess.run(['g++','-std=c++17','-O3','-fPIC','-shared',str(source),'-o',str(library)],check=True)
        lib=ct.CDLL(str(library.resolve()))
        self.function=lib.melt_barrier
        self.function.argtypes=[ct.c_int]+[ct.c_void_p]*7
        self.function.restype=ct.c_int
        self.work_function=lib.melt_barrier_work
        self.work_function.argtypes=[ct.c_int]+[ct.c_void_p]*6
        self.work_function.restype=ct.c_int
        self.extended_K=None

    def sample(self,labels):
        if any(labels[...,c].min()<ax[0]-2e-13 or labels[...,c].max()>ax[-1]+2e-13 for c,ax in enumerate(self.opt.axes)):
            raise InfeasibleLabels('Label outside registered domain')
        g=np.zeros_like(labels); det=np.empty(self.ids.shape); value=np.zeros(1,np.longdouble)
        rc=self.function(len(self.ids),*[a.ctypes.data for a in (self.ids,self.width,self.weight,labels,g,det,value)])
        if rc: raise InfeasibleLabels('Melted-corner determinant nonpositive')
        d=(labels-self.opt.base).reshape(-1,3)
        har=(self.opt.K@d).reshape(labels.shape)
        harmonic=.5*np.sum(d.ravel()*har.ravel(),dtype=np.longdouble)
        return harmonic+value[0],g+har,dict(harmonic=float(harmonic),barrier=float(value[0]),corner_min=float(det.min())),det

    def extract(self,grad,kind):
        return np.concatenate([cm.sum(grad[...,c])[getattr(cm,kind)] for c,cm in enumerate(self.components)])

    def melt_direction(self,flat):
        step=np.zeros_like(self.opt.base)
        for c,(cm,dc) in enumerate(zip(self.components,np.split(flat,self.cuts))):
            v=np.zeros(len(cm.held));v[cm.dead]=dc
            step[...,c]=v[cm.groups].reshape(self.opt.stencil.shape)
        return step

    def pre(self,grad):
        return np.concatenate([cm.extension.solve(part) for cm,part in zip(self.components,np.split(grad,self.cuts))])

    def local_work(self,old,new):
        """Signed unit-gauge energy difference at two float64 label maps."""
        if not (old.flags.c_contiguous and new.flags.c_contiguous
                and old.dtype==new.dtype==np.float64):
            raise ValueError('Local work needs contiguous float64 label maps')
        value=np.zeros(1,np.longdouble)
        args=(self.ids,self.width,self.weight,old,new,value)
        if self.work_function(len(self.ids),*[a.ctypes.data for a in args]):
            raise InfeasibleLabels('Invalid melted determinant in gauge work')
        if self.extended_K is None:
            self.extended_K=self.opt.K.astype(np.longdouble)
        x=old.reshape(-1,3).astype(np.longdouble)
        delta=new.reshape(-1,3).astype(np.longdouble)-x
        mid=(x-self.opt.base.reshape(-1,3))+.5*delta
        # Use the symmetric part for the exact stored quadratic form even
        # when element assembly leaves the last bit of K slightly asymmetric.
        harmonic=.5*np.sum(delta*((self.extended_K@mid)+(self.extended_K.T@mid)))
        return harmonic+value[0]

    def solve(self,start,name,output,*,outer_force=None):
        labels=start.copy()
        value,gr,parts,det=self.sample(labels)
        gradient=self.extract(gr,'dead')
        history=[]; alphas=[]; last=None
        tolerance=3.e-8
        if outer_force is not None:
            if not np.isfinite(outer_force) or outer_force<=0:
                raise ValueError('Use a positive finite outer force scale for adaptive gauge solves')
            tolerance=min(1.e-5,1.e-3*outer_force)
        with (output/(name+'.jsonl')).open('x') as log:
            while True:
                residual=absmax(gradient)
                log.write(json.dumps(dict(i=len(alphas),unit_energy=float(value),melt_residual=residual,
                    last=last,**parts))+'\n');log.flush()
                if residual<=tolerance:break
                y=gradient.copy(); weights=[]
                for s,v,inv in reversed(history):
                    a=inv*np.dot(s,y);y-=a*v;weights.append(a)
                y=self.pre(y)
                for (s,v,inv),a in zip(history,reversed(weights)):y+=s*(a-inv*np.dot(v,y))
                direction=-y
                if np.dot(direction,gradient)>=0:
                    history=[];direction=-self.pre(gradient)
                step=self.melt_direction(direction)
                # Arbitrarily large entry slopes at nearly singular corners are
                # scaled only in the search metric. The gauge energy is exact.
                scale=max(1.,*(absmax(step[...,k])/(.25*np.median(np.diff(ax))) for k,ax in enumerate(self.opt.axes)))
                step/=scale;direction/=scale
                slope=np.dot(gradient,direction)
                if not slope<0:raise ArithmeticError('Melt search is not descending')
                alpha=1.; reject=None
                while True:
                    trial=labels+alpha*step
                    try:
                        nv,ng,nparts,ndet=self.sample(trial)
                        if outer_force is not None:
                            dv=self.local_work(labels,trial)
                        else:
                            dstep=trial-labels
                            rel=(ndet-det)/det
                            if np.any(rel<=-1):raise InfeasibleLabels('Unresolved relative corner volume')
                            dv=np.sum(self.weight[:,None]*((ndet-det)-np.log1p(rel)),dtype=np.longdouble)
                            mid=((labels-self.opt.base)+.5*dstep).reshape(-1,3)
                            dv+=np.sum(dstep.reshape(-1,3)*(self.opt.K@mid),dtype=np.longdouble)
                        if not dv<1.e-4*alpha*slope:raise InfeasibleLabels('No resolved Armijo descent')
                        orient=chart_orientation(self.opt.axes,trial,self.opt.active)
                    except InfeasibleLabels as e:
                        reject=str(e);alpha*=.5
                        if np.array_equal(labels+alpha*step,labels) or (-alpha*slope<1e-28 and alpha<1e-4):
                            return labels,dict(stationary=False,reason=reject,accepted_steps=len(alphas),
                                alpha=alpha,residual=residual,target=tolerance,**parts)
                        continue
                    break
                next_gradient=self.extract(ng,'dead')
                s=alpha*direction;v=next_gradient-gradient;sy=np.dot(s,v)
                if sy>1.e-12*np.linalg.norm(s)*np.linalg.norm(v):
                    history.append((s,v,1./sy));history=history[-10:]
                last=dict(alpha=alpha,search_normalization=scale,delta=float(dv),chart_lower=orient['chart_det_global_lower'])
                labels,value,gr,parts,det,gradient=trial,nv,ng,nparts,ndet,next_gradient
                alphas.append(alpha)
        return labels,dict(stationary=True,accepted_steps=len(alphas),min_alpha=min(alphas,default=None),last=last,
                           residual=absmax(gradient),target=tolerance,**parts)


def run(scratch,case,output):
    output.mkdir(parents=True,exist_ok=False)
    base=scratch/'r4_eulerian_stage2_01'
    doc=read(base/case/'label_stall.json');fn=base/case/doc['file']
    if digest(fn)!=doc['sha256']:raise ValueError('Stall state hash differs')
    b,cp,axes=get_cp10(scratch)
    with np.load(fn) as f:
        if int(f['radius'])!=4 or int(f['winding'])!=1 or not np.array_equal(f['x'],cp['x']):
            raise ValueError('Incorrect baseline')
        labels=f['labels'];order=np.ascontiguousarray(f['field'][2,b.ids].transpose(0,2,1));eps=float(f['epsilon'])
    opt=LabelRelaxation(b,labels,order,eps)
    barrier=MeltBarrier(b,opt,output)
    solid=SpatialCrystal(b)
    E,grad,_,initial=solid.evaluate(labels,order)
    initial_unit,bg,start_parts,_=barrier.sample(labels)
    result=dict(case=case,epsilon=eps,sha256=doc['sha256'],start=start_parts,initial_chart=initial['min_chart_Jinv_bound'])
    def save():
        atomic_json(output/'readout.json',result)
    save()
    relaxed,info=barrier.solve(labels,'fixed_active_trace',output)
    result['inner']=info
    E2,g2,_,final=solid.evaluate(relaxed,order)
    v,g,parts,_=barrier.sample(relaxed)
    result['static']=dict(energy_before=E,delta_physical=E2-E,active_label_delta=absmax((labels-relaxed)[opt.stencil]),
        physical_force_delta=absmax((grad-g2)[opt.stencil]),total_field_held=True,chart=chart_orientation(axes,relaxed,opt.active),
        sampled_min_Jinv=final['min_Jinv_global'],sampled_min_active_Jinv=final['min_Jinv_active'],unit_gauge=float(v),
        gauge_delta_from_old_harmonic=eps*(float(v)-start_parts['harmonic']),new_gauge_delta=eps*float(v-initial_unit),
        gauge_force_on_free_active=eps*absmax(barrier.extract(g,'free')),gauge_to_elastic=eps*float(v)/E2,
        max_active_force=absmax(barrier.extract(g2+eps*g,'free')))
    np.savez_compressed(output/'fixed_active.npz',labels=relaxed)
    save()
    if not info['stationary']:return
    # A reduced gradient validation at fixed shape uses a normalized,
    # localized, two-sided perturbation on the independent melted labels.
    rng=np.random.default_rng(815)
    direction=barrier.melt_direction(rng.normal(size=sum(barrier.sizes)))
    direction/=absmax(direction)
    result['melt_gradient_check']=[]
    # Move slightly off stationarity for a nontrivial signed control.
    center=relaxed+2e-6*direction
    _,cg,_,_=barrier.sample(center)
    analytic=np.sum(cg*direction,dtype=np.longdouble)
    for h in (2e-8,1e-8):
        pair=[]
        for sign in (-1,1):
            w=center+sign*h*direction
            chart_orientation(axes,w,opt.active)
            pair.append(barrier.sample(w)[0])
        fd=(pair[1]-pair[0])/(2*h)
        result['melt_gradient_check'].append(dict(h=h,analytic=float(analytic),central=float(fd),difference=float(fd-analytic)))
    save()
    # Take exactly one trial in the production preconditioned active direction.
    # Preserve the nonlinear extension offset while predicting the new trace,
    # then re-settle the unit barrier separately at the first admissible size.
    active_gradient=barrier.extract(g2+eps*g,'free')
    direction=-opt._precondition(order)(active_gradient)
    x=opt._pack(relaxed)
    original=opt._lift(x)
    r0=E2+eps*float(v)
    result['outer_trials']=[]
    alpha=1.
    while alpha>=1e-4:
        trial=opt._lift(x+alpha*direction)+(relaxed-original)
        row=dict(alpha=alpha)
        result['outer_trials'].append(row)
        try:
            chart_orientation(axes,trial,opt.active)
            ep=solid.evaluate(trial,order,False)[0]; ug=barrier.sample(trial)[0]
            if ep+eps*ug >= r0+alpha*1e-4*np.dot(direction,active_gradient):
                raise InfeasibleLabels('Active trial fails reduced Armijo sufficient-decrease test')
        except (InfeasibleLabels,InvalidSpatialCrystal) as error:
            row['issue']=str(error);row['witness']=getattr(error,'witness',getattr(error,'diagnostics',{}))
            alpha*=.5;save();continue
        result['outer_predictor_energy_change']=float((ep-E2)+eps*(ug-v))
        save()
        settled,check=barrier.solve(trial,'single_active_probe',output)
        ep,gp,_,ph=solid.evaluate(settled,order)
        uf,gg,_,_=barrier.sample(settled)
        result['outer_probe']=dict(alpha=alpha,inner=check,reduced_energy_change=float((ep-E2)+eps*(uf-v)),
                 active_force=absmax(barrier.extract(gp+eps*gg,'free')),sampled_min_active_Jinv=ph['min_Jinv_active'],
                 orientation=chart_orientation(axes,settled,opt.active))
        save();return
    result['outer_probe']=dict(healthy=False,reason='No admissible sufficient-decrease predictor at alpha >= 1e-4')
    save()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path);p.add_argument('case',choices=('primary','secondary'));p.add_argument('output',type=Path)
    a=p.parse_args();run(a.scratch,a.case,a.output)
