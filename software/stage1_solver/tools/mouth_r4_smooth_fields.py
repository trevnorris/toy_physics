"""Held R4 order/carrier equilibration with the smooth pointwise gauge in the order work."""
import argparse
import json
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from scipy import sparse

from mouth_continuous_order import ContinuousTrial, I, KAPPA, WELL, norm
from mouth_cp10_label_condition import physical_frame
from mouth_cp10_label_gauge import VERTICES, gradient_basis
from mouth_exchange_carrier import BIAS
from mouth_r4_eulerian_relax_labels import FORCE_TOL
from mouth_r4_eulerian_release import block_axes, order_grid
from mouth_whole_r4_relax import FIELD_TOL


class GaugeGeometry:
    """Held physical quadrature shared by every label/order update."""
    def __init__(self, trial, *, aligned=False):
        self.aligned = aligned
        axes = block_axes(trial.b)
        nodes = trial.number[trial.b.ids.transpose(0,2,1)]
        self.base = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1)
        cells = np.indices(tuple(n-1 for n in nodes.shape)).reshape(3,-1).T
        self.ids = np.stack([np.ravel_multi_index(tuple((cells+v).T),nodes.shape) for v in VERTICES], axis=1)
        self.widths = np.column_stack([np.diff(axis)[cells[:,k]] for k,axis in enumerate(axes)])
        lo = self.base.reshape(-1,3)[self.ids[:,0]]
        columns = nodes.ravel()[self.ids].ravel()
        rows = np.repeat(np.arange(len(cells)),8)
        quadrature = (.5-.5/np.sqrt(3), .5+.5/np.sqrt(3))
        self.volume, self.physical_gradients, interpolation = [], [], []
        for q in product(quadrature, repeat=3):
            frame, hoop = physical_frame(trial.b, lo+self.widths*q)
            self.volume.append(.5*np.pi*self.widths.prod(axis=1)*hoop*np.linalg.det(frame))
            deriv = gradient_basis(q)[None,:,:]/self.widths[:,None,:]
            self.physical_gradients.append(deriv@np.linalg.inv(frame))
            basis = np.where(VERTICES,q,1-np.asarray(q)).prod(axis=1)
            interpolation.append(sparse.coo_matrix((np.tile(basis,len(cells)),(rows,columns)),
                                 shape=(len(cells),len(trial.first))))
        self.V = sparse.vstack(interpolation, format='csr')

    def work_density(self, labels):
        """Per-Gauss-point physical four-volume times the unweighted unit gauge."""
        points = labels.reshape(-1,3)[self.ids]
        offsets = (labels-self.base).reshape(-1,3)[self.ids]
        if self.aligned:
            from mouth_r4_aligned_stencil import barrier_means
            _, mean = barrier_means(points, self.widths)
            barrier = np.asarray(mean, dtype=float)
        else:
            barrier = np.zeros(len(points))
            for v in VERTICES:
                deriv = (points-points[:,:1]).transpose(0,2,1)@gradient_basis(v)/self.widths[:,None,:]
                det = np.linalg.det(deriv)
                if np.any(det <= 0.) or not np.isfinite(det).all():
                    raise ValueError('Nonpositive or nonfinite corner in the held gauge')
                barrier += (det-1-np.log(det))/8
        work = []
        for volume,basis in zip(self.volume,self.physical_gradients):
            grad = offsets.transpose(0,2,1)@basis
            work.append(volume*(barrier+.5*np.sum(grad**2,axis=(1,2))))
        return np.concatenate(work)


class SmoothOrder:
    """Use the native constrained Newton search with exact smooth order energy/gradient.

    The native trial Hessian already acts as a search metric. All accepted
    steps and the smooth gradient at the midpoint use the full functional.
    """
    def __init__(self, trial, labels, epsilon):
        self.trial, self.epsilon = trial, epsilon
        if not hasattr(trial,'smooth_order_geometry'):
            trial.smooth_order_geometry = GaugeGeometry(trial)
        geometry = trial.smooth_order_geometry
        self.Vg = geometry.V
        self.work = geometry.work_density(labels)

    def __getattr__(self,name):
        # ContinuousTrial.solve_order uses the held trial's physical mesh and
        # preconditioner, while dispatching order_terms on this object.
        return getattr(self.trial,name)

    def gauge(self, c):
        x = np.clip(self.Vg@c,0.,1.)
        value = np.sum(self.work*((1-x)**2*(1+2*x)),dtype=np.longdouble)
        grad = self.Vg.T@(self.work*(-6*x*(1-x)))
        return value, grad

    def order_terms(self, c, n, gp, load, rest):
        # Accumulate the physical weak derivative directly as in the native
        # field evaluator. K @ c cancels large off-diagonal terms even when c
        # is locally constant. Keep K as the sparse Newton search metric.
        t,m=self.trial,self.model
        physical_nodes=c[t.number]
        physical_order=m.V@physical_nodes
        gc=[G@physical_nodes for G in m.G]
        weights=m.hq
        grad_sq=sum(b*b for b in gc)
        local=.5*KAPPA*grad_sq + WELL*(physical_order*(1-physical_order))**2
        potential=2*WELL*physical_order*(1-physical_order)*(1-2*physical_order)
        nodal=m.V.T@(weights*potential)
        if rest:
            flux=[KAPPA*b for b in gc]
        else:
            den=n+I*grad_sq
            dot=sum(a*b for a,b in zip(gc,gp))
            pi=I*n*dot/den
            local-=.5*pi*dot
            flux=[KAPPA*b-pi*(a-I*dot/den*b) for a,b in zip(gp,gc)]
        for G,current in zip(m.G,flux):
            nodal+=G.T@(weights*current)
        physical=float(np.sum(np.multiply(weights,local,dtype=np.longdouble),dtype=np.longdouble))+load@c
        gradient=t.pull.T@nodal+load
        unit, g = self.gauge(c)
        return physical+float(self.epsilon*unit), gradient+self.epsilon*g, physical_order


def material_load(t, solid, labels, u):
    elastic,_,source,_ = solid.evaluate(labels,order_grid(t,u),False)
    distributed = np.bincount(t.b.ids.ravel(),weights=source.transpose(0,2,1).ravel(),
                              minlength=t.model.size)
    difference = float(elastic-distributed@u[2])
    if abs(difference)>2e-10*max(1.,abs(elastic)):
        raise ArithmeticError('Crystal order work differs from physical nodal deposition')
    return t.ext+distributed, difference


def state(driver, t, order_solver, u, exchange):
    m=t.model
    physical,g=m.evaluate(u)
    chi=u[2,t.first]
    unit,unit_gradient = order_solver.gauge(chi)
    physical_order=t.pull.T@g[2]
    joined=physical_order+driver.epsilon*unit_gradient
    projected=norm(t.projected(chi,joined))
    kkt=np.where(chi==0,np.minimum(joined,0),np.where(chi==1,np.maximum(joined,0),joined))
    amplitude,phase=(g[c,m.masks[c]] for c in (0,1))
    metric=max(norm(m.solvers[0](amplitude)),norm(m.solvers[1](phase)))
    g[2,t.first] += driver.epsilon*unit_gradient
    return float(physical),dict(physical_energy=float(physical),unit_gauge=float(unit),
                gauge_energy=float(driver.epsilon*unit),phase_raw=norm(phase),amplitude_raw=norm(amplitude),
                order_raw=norm(joined[t.free]),order_kkt=norm(kkt[t.free]),order_projected=projected,
                carrier_metric=metric,shared_work_error=exchange,
                order_knot_transport=driver.forces_at_order_knots(t,u,g))


def read_fields(driver, t, solid, labels, fields):
    """Simultaneous residuals at precisely the supplied labels and physical fields."""
    load,exchange=material_load(t,solid,labels,fields)
    t.model.material_load=load
    order_solver=SmoothOrder(t,labels,driver.epsilon)
    energy,info=state(driver,t,order_solver,fields,exchange)
    return energy,load,info


def coherent(driver,t,solid,labels,fields,*,field_target=FIELD_TOL,raw_target=FORCE_TOL):
    """Rematch chi, amplitude and phase with the active and melted labels fixed."""
    m=t.model
    u=t.transported(fields)
    u[0,~m.masks[0]]=1.
    u[1,m.feed],u[1,m.bulk]=0.,BIAS
    load,exchange=material_load(t,solid,labels,u)
    m.material_load=load
    # The field model also allows an independently supplied quadrature-linear
    # elastic density, separate from the fixed crystal load on the order nodes.
    joined=t.pull.T@load+t.V.T@(m.hq*m.D[:,-1])
    order_solver=SmoothOrder(t,labels,driver.epsilon)
    physical,info=state(driver,t,order_solver,u,exchange)
    before,old_gauge=physical,info['gauge_energy']
    order_tol,carrier_tol=min(2e-9,field_target/2),min(8e-8,field_target/2)
    iterations=[]
    while max(info[k] for k in ('order_kkt','order_projected','carrier_metric'))>=field_target or max(
            info[k] for k in ('phase_raw','amplitude_raw'))>raw_target:
        if info['order_kkt']>=field_target and info['order_projected']<order_tol:
            order_tol=min(order_tol,.5*info['order_projected'])
        if max(info[k] for k in ('phase_raw','amplitude_raw'))>raw_target:
            carrier_tol=min(carrier_tol,.5*info['carrier_metric'])
        os,_=ContinuousTrial.solve_order(order_solver,u,joined,False,order_tol)
        cs,_=t.carrier(u,False,carrier_tol)
        iterations.append([os,cs])
        physical,info=state(driver,t,order_solver,u,exchange)
        if not os and not cs:
            raise ArithmeticError(f'Smooth order/carrier rematch made no field step: {info}')
    delta=(physical-before)+(info['gauge_energy']-old_gauge)
    if delta>64*abs(np.spacing(before)):
        raise ArithmeticError(f'Smooth order/carrier rematch raised total energy by {delta}')
    info['order_step_energy']=dict(physical_change=physical-before,old_gauge=old_gauge,
                                      new_gauge=info['gauge_energy'],delta=delta)
    info['field_iterations']=iterations
    return u,physical,load,info


def check(area,output):
    from mouth_r4_smooth_static import frozen
    from mouth_r4_smooth_order_boundary import gauge_coefficients
    from mouth_r4_projected_fields import EXPECTED, digest
    from mouth_r4_eulerian_solid import SpatialCrystal
    parent=Path(__file__).resolve().parents[3]
    with TemporaryDirectory() as directory:
        driver,thin,fields,labels,order=frozen(area,Path(directory))
        # Use the actual production preconditioner, including raw carrier work.
        del thin
        t=driver.trial()
        solid=SpatialCrystal(t.b)
        physical,load,info=read_fields(driver,t,solid,labels,fields)
        c=fields[2,t.first]
        solver=SmoothOrder(t,labels,driver.epsilon)
        unit,gradient=solver.gauge(c)
        coeff,_,gauge=gauge_coefficients(t,labels,fields)
        j=t.pull.T@load+t.V.T@(t.model.hq*t.model.D[:,-1])
        n=(t.model.D[:,4]*(t.model.V@fields[0]))**2
        gp=[G@fields[1] for G in t.model.G]
        _,joined,_=solver.order_terms(c,n,gp,j,False)
        _,g=t.model.evaluate(fields)
        physicalg=t.pull.T@g[2]
        _,assembled,_=t.order_terms(c,n,gp,j,False)
        assembly_delta=assembled-physicalg
        worst=int(np.argmax(abs(assembly_delta)))
        flux=t.pull.T@sum(G.T@(t.model.hq*(KAPPA*(G@fields[2]))) for G in t.model.G)
        diffusion_delta=t.K@c-flux
        row=t.K.getrow(worst)
        # Use explicit CSR terms for the sign budget, retaining the solver's
        # native K @ c for its actual accumulated residual.
        terms=row.data*c[row.indices]
        assembly=dict(node=worst,free=bool(t.free[worst]),order=float(c[worst]),depth=float(t.depth[worst]),
            signed_total=float(assembly_delta[worst]),signed_diffusion=float(diffusion_delta[worst]),
            assembled_positive_terms=float(terms[terms>0].sum()),assembled_negative_terms=float(terms[terms<0].sum()),
            max_free=norm(assembly_delta[t.free]),max_fixed=norm(assembly_delta[~t.free]),
            max_after_subtracting_diffusion=norm(assembly_delta-diffusion_delta))
        record=json.loads((parent/'research/pde_ledger_v3/_measurements/mouth_r4_projected_fields.json').read_text())['held_stop']
        deltas=dict(unit_energy=float(unit-gauge['new_unit']),unit_all_order_gradients=norm(gradient-coeff[0]),
            full_order_solver_vs_independent=norm(joined-(physicalg+driver.epsilon*coeff[0])),
            physical_energy=float(physical-record['physical_energy']),
            amplitude_raw=info['amplitude_raw']-record['amplitude_raw'],phase_raw=info['phase_raw']-record['phase_raw'])
        result=dict(radius=4,state_sha256=EXPECTED['held_stop'],script_sha256=digest(Path(__file__)),
            full_order_coordinate_count=len(c),free_order_coordinate_count=int(t.free.sum()),
            gauge=info['gauge_energy'],unit_gauge=float(unit),deltas={k:float(v) for k,v in deltas.items()},
            previous_assembled_diffusion=assembly,diagnostics=info)
        result['passed']=(abs(deltas['unit_energy'])<1e-12 and deltas['unit_all_order_gradients']<1e-12
            and deltas['full_order_solver_vs_independent']<1e-10 and max(abs(deltas[k]) for k in
            ('physical_energy','amplitude_raw','phase_raw'))==0)
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(json.dumps(result,allow_nan=False),flush=True)
        if not result['passed']:
            raise ArithmeticError('Held field implementation comparison failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    check(a.area.resolve(),a.output.resolve())
