"""Matrix-free Newton–CG on the joint free crystal labels at fixed order and shape.

The accepted objective uses the existing SpatialCrystal and pointwise Gauge
evaluators. Physical H·v differentiates the native moving-cut gradient on
both sides of the weld; the harmonic and determinant terms are differentiated
at fixed spatial quadrature, without a nested melted-region optimization.
"""
import ctypes as ct
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
import os
from pathlib import Path
import subprocess
import tempfile
from time import perf_counter

import numpy as np
from scipy import sparse

from mouth_cp10_label_gauge import VERTICES, gradient_basis
from mouth_r4_corner_projection import ProjectionError, project_halfspaces
from mouth_r4_eulerian_relax_labels import FORCE_TOL, InfeasibleLabels, LabelRelaxation, absmax, chart_orientation
from mouth_r4_eulerian_solid import InvalidSpatialCrystal
from mouth_r4_full_label_work import coloring
from mouth_r4_projected_fields import location
from mouth_r4_smooth_fields import GaugeGeometry
from mouth_r4_smooth_static import Gauge, RULE

RELATIVE_STEP = 2.**-20
UNIT_MELT_TOL = 3e-8


def feasible_steepest(gradient, normals, diagonal):
    """Find strict descent or prove that the gradient belongs to the normal cone.

    Only used if direct KKT verification and its high-precision retry could
    not produce a descending projection. The determinant tangents are homogeneous.
    """
    A = np.asarray(normals)
    touched = np.any(A != 0,axis=0)
    if np.any(gradient[~touched]):
        return np.where(touched,0.,-diagonal*gradient),dict(method='outside corner support')
    # All remaining gradient components are supported on a few corner nodes.
    # Use exact binary-input rationals to distinguish zero from tiny descent.
    from sympy import Matrix, Rational, zeros
    from sympy.matrices.exceptions import NonInvertibleMatrixError
    rows = Matrix(A[:,touched].tolist()).applyfunc(Rational)
    g = Matrix(gradient[touched].tolist()).applyfunc(Rational)
    for count in range(rows.rows+1):
        for bank in combinations(range(rows.rows),count):
            lam = zeros(rows.rows,1)
            if bank:
                active = rows[list(bank),:]
                try:
                    mu = (active*active.T).inv()*(active*g)
                except NonInvertibleMatrixError:
                    continue
                if min(mu)<0:
                    continue
                for i,y in zip(bank,mu):
                    lam[i] = y
            candidate = rows.T*lam-g
            if min(rows*candidate)<0:
                continue
            proof = dict(method='exact rational normal cone', rows=list(bank),
                         original_row_multipliers=[str(t) for t in lam])
            if candidate.is_zero_matrix:
                return None,proof
            out = np.zeros_like(gradient)
            out[touched] = [float(t) for t in candidate]
            return out,proof
    raise ArithmeticError('Exact positive-semidefinite cone projection has no KKT representation')


@lru_cache(maxsize=1)
def corner_kernel():
    source=Path(__file__).with_suffix('.cpp')
    key=sha256(source.read_bytes()+source.with_name('mouth_r4_melt_barrier.cpp').read_bytes()).hexdigest()
    path=source.parents[3]/'_scratch/s11c/r4_newton_build'
    path.mkdir(exist_ok=True,parents=True)
    lib=path/(key+'.so')
    if not lib.is_file():
        fd,tmp=tempfile.mkstemp(dir=path,suffix='.so')
        os.close(fd)
        try:
            subprocess.run(['g++','-std=c++17','-O3','-fPIC','-shared',str(source),'-o',tmp],check=True)
            os.replace(tmp,lib)
        finally:
            Path(tmp).unlink(missing_ok=True)
    handle=ct.CDLL(str(lib))
    for name,n in [('corner_tangent',7),('melt_barrier_work',6)]:
        f=getattr(handle,name)
        f.argtypes=[ct.c_int]+[ct.c_void_p]*n
        f.restype=ct.c_int
    return handle


class WeightedGauge:
    """Fixed-order harmonic stiffness and corner barrier quadrature with s(Q1 χ)."""
    def __init__(self,trial,order):
        if not hasattr(trial,'smooth_order_geometry'):
            trial.smooth_order_geometry=GaugeGeometry(trial)
        geometry=trial.smooth_order_geometry
        self.ids=np.ascontiguousarray(geometry.ids,dtype=np.int32)
        self.width=np.ascontiguousarray(geometry.widths)
        self.base=geometry.base
        cell_order=order.ravel()[self.ids]
        weights=np.zeros(len(self.ids))
        stiffness=np.zeros((len(self.ids),8,8))
        elastic_diagonal=np.zeros_like(order).ravel()
        for q,vol,grads in zip(RULE,geometry.volume,geometry.physical_gradients):
            N=np.prod(np.where(VERTICES,q,1-np.array(q)),axis=1)
            chi=cell_order@N
            s=(1-chi)**2*(1+2*chi)
            weights+=vol*s
            gram=grads@grads.transpose(0,2,1)
            stiffness+=(vol*s)[:,None,None]*gram
            np.add.at(elastic_diagonal,self.ids.ravel(),
                      ((vol*chi)[:,None]*np.diagonal(gram,axis1=1,axis2=2)).ravel())
        self.elastic_diagonal=elastic_diagonal.reshape(order.shape)
        self.weight=np.ascontiguousarray(weights/8)
        ni=order.size
        self.K=sparse.coo_matrix((stiffness.ravel(),(np.repeat(self.ids,8,axis=1).ravel(),
                    np.tile(self.ids,(1,8)).ravel())),shape=(ni,ni)).tocsr()
        self.lib=corner_kernel()

    def tangent(self,labels,step=None):
        out=np.zeros_like(labels)
        args=[self.ids,self.width,self.weight,labels,step,out if step is not None else None,
              out if step is None else None]
        rc=self.lib.corner_tangent(len(self.ids),*[a.ctypes.data if a is not None else None for a in args])
        if rc:
            raise InfeasibleLabels('Nonpositive chart corner in gauge tangent')
        if step is None:
            return out+self.K.diagonal().reshape(out.shape[:-1])[...,None]
        return out+(self.K@step.reshape(-1,3)).reshape(step.shape)

    def work(self,old,new):
        signed=np.zeros(1,np.longdouble)
        args=(self.ids,self.width,self.weight,old,new,signed)
        if self.lib.melt_barrier_work(len(self.ids),*[a.ctypes.data for a in args]):
            raise InfeasibleLabels('Nonpositive corner in signed gauge work')
        # Difference first, then apply the midpoint identity to the quadratic.
        a=old.reshape(-1,3).astype(np.longdouble)
        dv=new.reshape(-1,3).astype(np.longdouble)-a
        midpoint=(a-self.base.reshape(-1,3))+.5*dv
        return signed[0]+np.sum(dv*(self.K@midpoint),dtype=np.longdouble)


class JointCrystal:
    gauge_class = Gauge
    weighted_class = WeightedGauge

    def __init__(self,driver,trial,solid,labels,order):
        self.trial,self.solid,self.order=trial,solid,order
        self.epsilon=driver.epsilon
        self.opt=LabelRelaxation(trial.b,labels,order,self.epsilon)
        self.components=self.opt._components()
        self.cuts=np.cumsum([np.count_nonzero(~cm.held) for cm in self.components])[:-1]
        self.active=np.concatenate([cm.free[~cm.held] for cm in self.components])
        self.gauge=self.gauge_class(trial,driver.reference_depth,order)
        self.weighted=self.weighted_class(trial,order)
        self.initial_force=None
        self.previous_norm=None
        self.tight=[]
        widths=[]
        for c,cm in enumerate(self.components):
            w=np.diff(self.opt.axes[c])
            node_width=np.r_[w[0],np.minimum(w[:-1],w[1:]),w[-1]]
            ix=np.unravel_index(cm.representative[~cm.held],order.shape)[c]
            widths.append(node_width[ix])
        self.width=np.concatenate(widths)
        self.last_hvp=None
        self.jacobi=None
        self.jacobi_anchor=None
        self.jacobi_stats={}

    def lift(self,packed):
        result=np.zeros_like(self.opt.base)
        for c,(cm,block) in enumerate(zip(self.components,np.split(packed,self.cuts))):
            values=np.zeros(len(cm.held))
            values[~cm.held]=block
            result[...,c]=values[cm.groups].reshape(self.order.shape)
        return result

    def reduce(self,nodal):
        return np.concatenate([cm.sum(nodal[...,c])[~cm.held] for c,cm in enumerate(self.components)])

    def compatible(self,labels):
        """Honor a new axis registry after a field rematch, moving strictly melted copies only."""
        corrected=labels.copy()
        for c,cm in enumerate(self.components):
            delta=(labels-self.opt.base)[...,c].ravel()
            values=delta[cm.representative].copy()
            values[cm.free]=delta[cm.source[cm.free]]
            values[cm.held]=0.
            mismatch=(values[cm.groups]-delta).reshape(self.order.shape)
            if absmax(mismatch[self.opt.stencil])>2e-13:
                raise InfeasibleLabels('Incompatible registered labels on the active trace')
            dead=~self.opt.stencil
            corrected[...,c][dead]=self.opt.base[...,c][dead]+values[cm.groups].reshape(self.order.shape)[dead]
        return corrected

    def evaluate(self,labels):
        chart=chart_orientation(self.opt.axes,labels,self.opt.active)
        elastic,ge,_,physics=self.solid.evaluate(labels,self.order)
        if (physics['min_Jinv_active']<=np.sqrt(np.finfo(float).eps)
                or chart['chart_det_active_lower']<=128*np.finfo(float).eps):
            raise InfeasibleLabels('Ordered label orientation unresolved',**physics['worst_active'])
        cells,unit=self.gauge.evaluate(labels,np.zeros(6),gradient=True)
        gu=self.reduce(np.asarray(unit,float))
        gradient=self.reduce(ge+self.epsilon*np.asarray(unit,float))
        ge=self.reduce(ge)
        if absmax(ge[~self.active])>1e-12 or not np.isfinite(gradient).all():
            raise ArithmeticError('Unresolved crystal force on strictly melted labels or invalid gradient')
        largest=absmax(gradient[self.active])
        melt=absmax(gu[~self.active])
        energy=float(cells.sum())
        m=dict(**physics,chart=chart,elastic_energy=elastic,gauge_energy=self.epsilon*energy,
               unit_gauge_energy=energy,gauge_to_elastic=self.epsilon*energy/elastic,
               max_free_label_force=largest,max_free_elastic_force=absmax(ge[self.active]),
               max_unit_melt_residual=melt,max_dead_elastic_force=absmax(ge[~self.active]),
               gradient_l2=float(np.linalg.norm(gradient)),
               stationary=largest<=FORCE_TOL and melt<=UNIT_MELT_TOL)
        worst_index=int(np.argmax(np.where(self.active,abs(gradient),-np.inf)))
        idx=worst_index
        for c,cm in enumerate(self.components):
            mobile=np.flatnonzero(~cm.held)
            if idx < len(mobile):
                group=mobile[idx]
                node=np.unravel_index(cm.source[group],self.order.shape)
                m['worst_free_active']=dict(component=c,group=int(group),chi=float(self.order[node]),
                    signed_force=-float(gradient[worst_index]),
                    **location(self.trial,self.opt.axes,node))
                break
            idx-=len(mobile)
        return elastic+m['gauge_energy'],gradient,m

    def hvp(self,labels,v,relative_step=None):
        """Differentiate the cubic crystal and the gauge at the same base state."""
        if self.solid.body.data.collar_coeff:
            u=self.lift(v)
            result=self.solid.hessian_action(labels,self.order,u)
            result+=self.epsilon*self.weighted.tangent(labels,u)
            self.last_hvp=dict(h=None,relative_step=0.,method='forward_tangent')
            return self.reduce(result)
        # The original triangular registry has separate internal face switches.
        rel=RELATIVE_STEP if relative_step is None else float(relative_step)
        if rel<=0 or not np.isfinite(rel):
            raise ValueError('Expected positive finite cell-relative derivative width')
        u=self.lift(v)
        result=self.epsilon*self.weighted.tangent(labels,u)
        va=np.where(self.active,v,0.)
        scale=absmax(va/self.width)
        if scale:
            step=self.lift(va)
            h=rel/scale
            for _ in range(32):
                try:
                    plus=labels+h*step
                    minus=labels-h*step
                    if np.array_equal(plus,minus):
                        raise ArithmeticError('Unresolved physical Hessian-vector displacement')
                    gp=self.solid.evaluate(plus,self.order)[1]
                    gm=self.solid.evaluate(minus,self.order)[1]
                    break
                except InvalidSpatialCrystal:
                    h*=.5
            else:
                raise ArithmeticError('Unable to obtain feasible symmetric physical Hessian-vector samples')
            result+=(gp-gm)/(2*h)
            self.last_hvp=dict(h=float(h),relative_step=float(h*scale))
        else:
            self.last_hvp=dict(h=None,relative_step=0.)
        return self.reduce(result)

    def inverse_diagonal(self,labels,chart_lower):
        """Simultaneously read true independent-coordinate curvatures on disjoint Q1 supports."""
        # One frozen Jacobi matrix per fixed-order burst. Rebuild sooner if
        # the determinant margin changes sufficiently to change barrier curvature.
        if self.jacobi is not None and .5*self.jacobi_anchor <= chart_lower <= 2*self.jacobi_anchor:
            return self.jacobi,dict(**self.jacobi_stats,rebuilt=False)
        diagonal=np.empty_like(self.width)
        cursor=0
        colors=[]
        hmin=0. if self.solid.body.data.collar_coeff else RELATIVE_STEP
        for cm in self.components:
            groups=np.flatnonzero(~cm.held)
            color=coloring(cm,self.weighted.ids)
            cell_colors=color[cm.groups[self.weighted.ids]]
            cell_groups=cm.groups[self.weighted.ids]
            count=int(max(color))+1
            colors.append(count)
            for k in range(count):
                # Multiple appearances of a single registered group at an axis
                # are permitted; distinct varied groups must never share a cell.
                selected=cell_colors==k
                first=np.max(np.where(selected,cell_groups,-1),axis=1)
                if np.any(selected & (first[:,None]!=cell_groups)):
                    raise ArithmeticError('Invalid Jacobi support coloring')
                ids=cursor+np.flatnonzero(color[groups]==k)
                direction=np.zeros_like(self.width)
                direction[ids]=self.width[ids]
                column_sum=self.hvp(labels,direction)
                diagonal[ids]=column_sum[ids]/direction[ids]
                if self.last_hvp['h'] is not None:
                    hmin=min(hmin,self.last_hvp['relative_step'])
            cursor+=len(groups)
        if cursor!=len(self.width) or not np.all(np.isfinite(diagonal)):
            raise ArithmeticError('Incomplete or nonfinite packed Hessian diagonal')
        # Signs belong to the true Hessian. Use magnitudes only in the
        # positive preconditioner; truncated CG tests the unmodified H itself.
        # Homogeneous lower bound at each node is a roundoff safeguard for a
        # zero diagonal, with the same epsilon scaling as its gauge block.
        unit=self.reduce(self.weighted.tangent(labels))
        physical=self.reduce(np.broadcast_to(self.weighted.elastic_diagonal[...,None],labels.shape))
        floor=np.finfo(float).eps*(self.epsilon*unit+physical)
        positive=np.maximum(abs(diagonal),floor)
        if np.any(positive<=0):
            raise ArithmeticError('Joint label lacks positive curvature scale')
        self.jacobi=1./positive
        self.jacobi_anchor=chart_lower
        self.jacobi_stats=dict(colors_by_component=colors,hvp_calls=sum(colors),
            negative_diagonal=int((diagonal<0).sum()),zero_diagonal=int((diagonal==0).sum()),
            floored=int((positive>abs(diagonal)).sum()),minimum_relative_step=hmin,
            active_range=[float(diagonal[self.active].min()),float(diagonal[self.active].max())],
            melted_range=[float(diagonal[~self.active].min()),float(diagonal[~self.active].max())],
            chart_at_assembly=chart_lower)
        return self.jacobi,dict(**self.jacobi_stats,rebuilt=True)

    def corner_normal(self,labels,cell,local):
        cell=np.asarray(cell,int)
        width=np.array([np.diff(ax)[i] for ax,i in zip(self.opt.axes,cell)])
        ix=tuple((VERTICES+cell).T)
        Q=gradient_basis(np.asarray(local))/width
        points=labels[ix]
        D=(points-points[:1]).T@Q
        cof=np.column_stack([np.cross(D[:,(k+1)%3],D[:,(k+2)%3]) for k in range(3)])
        normal=np.zeros_like(labels)
        normal[ix]=Q@cof.T
        return self.reduce(normal),float(D[:,0]@cof[:,0])

    def step(self,labels,value,g,material):
        start=perf_counter()
        ticks=dict(hvp=0.,preconditioner=0.,projection=0.,trial_evaluation=0.)
        diag,jacobi=self.inverse_diagonal(labels,material['chart']['chart_det_global_lower'])
        ticks['preconditioner']=perf_counter()-start
        max_force=material['max_free_label_force']
        active_polish=max_force<=FORCE_TOL and material['max_unit_melt_residual']>UNIT_MELT_TOL
        moving=~self.active if active_polish else np.ones(len(g),bool)
        rhs=-np.where(moving,g,0.)
        norm0=float(np.linalg.norm(rhs))
        if not norm0:
            raise ArithmeticError('Newton direction requested at zero mobile gradient')
        if self.initial_force is None:
            self.initial_force=max(max_force,FORCE_TOL)
        # Eisenstat–Walker forcing (choice 2), capped by a residual-dependent
        # envelope so the inner system tightens at the prescribed outer target.
        ew=.1 if self.previous_norm is None else .9*(norm0/self.previous_norm)**1.5
        target_ratio=max(0.,(max_force-FORCE_TOL)/self.initial_force)
        eta=max(.001,min(.1,ew,.001+.099*min(1.,target_ratio)))
        self.previous_norm=norm0
        p=np.zeros_like(g)
        residual=rhs.copy()
        z=diag*residual
        direction=z.copy()
        rz=float(residual@z)
        rel=1.
        reason='CG forcing attained'
        cg=0
        replacements=[]
        actual=None
        def check_direction(raw):
            at=perf_counter()
            hp=self.hvp(labels,raw)
            ticks['hvp']+=perf_counter()-at
            full=g+hp
            restricted=np.where(moving,-full,0.)
            return restricted,float(np.linalg.norm(restricted)/norm0),float(np.linalg.norm(full)/np.linalg.norm(g))
        for cg in range(1,int(moving.sum())+1):
            at=perf_counter()
            hd=self.hvp(labels,direction)
            ticks['hvp']+=perf_counter()-at
            hd[~moving]=0.
            curvature=float(direction@hd)
            if not np.isfinite(curvature):
                raise ArithmeticError('Nonfinite Newton curvature')
            if curvature<=0:
                reason='nonpositive curvature'
                if cg==1 or g@p>=0:
                    p=diag*rhs
                break
            beta=rz/curvature
            p+=beta*direction
            residual-=beta*hd
            rel=float(np.linalg.norm(residual)/norm0)
            if rel<=eta:
                true_residual,true_rel,full_rel=check_direction(p)
                if true_rel<=eta:
                    actual=(true_rel,full_rel)
                    break
                replacements.append(dict(iteration=cg,recursive=rel,explicit=true_rel))
                # Nonlinear secant noise may make a very small recurrence
                # residual unattainable for the actual matrix-free operator.
                # Restart CG with the measured residual and detect stagnation.
                if len(replacements)>1 and true_rel>=replacements[-2]['explicit']*(1-1e-8):
                    raise InfeasibleLabels('CG residual replacements stagnated above forcing tolerance',
                        replacements=replacements,eta=eta)
                residual=true_residual
                z=diag*residual
                direction=z.copy()
                rz=float(residual@z)
                continue
            z=diag*residual
            new_rz=float(residual@z)
            direction=z+(new_rz/rz)*direction
            rz=new_rz
        else:
            raise ArithmeticError('CG exhausted its system dimension without attaining the forcing tolerance')
        if g@p>=0:
            p=diag*rhs
            reason+='; used descending Jacobi direction'
            actual=None
        if actual is None:
            _,true_rel,full_rel=check_direction(p)
            actual=(true_rel,full_rel)
        tests=[]
        projections=0
        fallback=False
        projection_events=[]
        # Constrain previously identified near-inverting sites in the same
        # scaled metric used by CG. Certify the whole chart at every trial.
        def project(raw, *, steepest):
            nonlocal projections
            at=perf_counter()
            normals=[]
            sites=[]
            for cell,local in self.tight:
                normal,_=self.corner_normal(labels,cell,local)
                normal[~moving]=0.
                if np.any(normal):
                    normals.append(normal)
                    sites.append(dict(cell=list(cell),local=list(local)))
            if normals:
                projections+=1
                try:
                    out,info=project_halfspaces(raw,normals,gradient=g,diagonal=diag,return_info=True)
                except ProjectionError as error:
                    corner_data=[]
                    for site in sites:
                        nodes=VERTICES+site['cell']
                        corner_data.append(dict(**site,vertex_chi=self.order[tuple(nodes.T)].tolist()))
                    details=dict(sites=corner_data,projection_error=str(error),projection=error.diagnostics)
                    projection_events.append(dict(sites=sites,method='extended',failed=True,reason=str(error)))
                    print(json.dumps(dict(kind='corner_projection_fallback',**details)),flush=True)
                    if not steepest:
                        return None
                    alternate,proof=feasible_steepest(g*moving,normals,diag)
                    if alternate is None:
                        raise InfeasibleLabels('Certified absence of strict descent on the mobile corner cone',
                                               **details,certificate=proof) from error
                    out,info=project_halfspaces(alternate,normals,gradient=g,diagonal=diag,return_info=True)
                    info['alternative']=proof
                event=dict(**info,sites=sites)
                projection_events.append(event)
                if info['method']=='extended':
                    print(json.dumps(dict(kind='corner_projection_fallback',**event)),flush=True)
            else:
                out=raw
            ticks['projection']+=perf_counter()-at
            return out
        for steepest,raw in ((False,p),(True,diag*rhs)):
            plan=project(raw,steepest=steepest)
            alpha=1.
            while plan is not None and alpha>1e-4:
                slope=float(g@plan)
                if slope>=0:
                    break
                trial=labels+alpha*self.lift(plan)
                at=perf_counter()
                try:
                    nv,ng,after=self.evaluate(trial)
                    # No O(10^5) reservoir energy enters the fixed-field work.
                    # For melt polish physical work is identically zero and
                    # corner determinant increments are evaluated before rounding.
                    de=after['elastic_energy']-material['elastic_energy']
                    dg=self.epsilon*self.weighted.work(labels,trial)
                    delta=float(de+dg)
                    if (delta<0 and delta<=1e-4*alpha*slope and
                            (de!=0 or np.array_equal(trial[self.opt.stencil],labels[self.opt.stencil]))):
                        ticks['trial_evaluation']+=perf_counter()-at
                        ticks['total']=perf_counter()-start
                        return trial,nv,ng,after,dict(alpha=alpha,energy_change=delta,physical_delta=de,
                            gauge_delta=float(dg),slope=slope,gradient_l2_before=float(np.linalg.norm(g)),
                            gradient_l2_after=float(np.linalg.norm(ng)),cg_iterations=cg,cg_reason=reason,
                            cg_relative_residual=actual[0],cg_full_relative_residual=actual[1],cg_recursive_relative_residual=rel,
                            forcing_attained=actual[0]<=eta,cg_residual_replacements=replacements,
                            jacobi=jacobi,forcing_eta=eta,polished_melt_only=active_polish,
                            jacobi_fallback=fallback,projections=projections,backtracks=len(tests),
                            projection_events=projection_events,
                            projection_fallbacks=sum(e['method']=='extended' for e in projection_events),
                            rejections=tests,timing_seconds=ticks,min_Jinv=after['min_Jinv_global'],
                            min_ordered_Jinv=after['min_Jinv_active'])
                    exc=InfeasibleLabels('Insufficient resolved signed work',delta=delta)
                except (InfeasibleLabels,InvalidSpatialCrystal) as e:
                    exc=e
                ticks['trial_evaluation']+=perf_counter()-at
                issue=getattr(exc,'witness',getattr(exc,'diagnostics',{}))
                tests.append(dict(alpha=alpha,reason=str(exc),site=issue))
                where=issue.get('chart_fractions',issue.get('local'))
                if where is not None and 'cell' in issue:
                    pair=(tuple(issue['cell']),tuple(where))
                    if pair not in self.tight:
                        self.tight.append(pair)
                        plan=project(raw,steepest=steepest)
                        # This is a different ray. A step fraction rejected on
                        # the previous ray says nothing about its feasible length.
                        alpha=1.
                        continue
                alpha*=.5
            fallback=True
        raise InfeasibleLabels('Joint Newton and projected Jacobi both failed signed Armijo at healthy fractions',
                               tests=tests,material=material,cg_iterations=cg,cg_reason=reason,
                               projection_events=projection_events)
