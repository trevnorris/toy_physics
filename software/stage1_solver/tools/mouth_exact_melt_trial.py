"""R4 melted-collar trial with exact support for the bounded order field.

A zero-cost continuation of W on invalid labels gives the existing lower
bound on reordering work. A positive order coefficient outside the actual
fold carries its original, fully defined quadrature load. The certificate
does not enter the action or its finite differences.
"""
from itertools import product

import numpy as np

from mouth_active_prism_orientation import active_orientation
from mouth_continuous_order import norm
from mouth_melted_collar_survival import MeltTrial, brief
from mouth_prism_support_decision import decide


class VerifiedOverlap(ArithmeticError):
    def __init__(self, record):
        self.record = record
        super().__init__(f"Verified positive order on nonpositive J: {record['witness']}")


class ExactMeltTrial(MeltTrial):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.b.data.exact_melt_support = 1
        self.overlap_sink = None

    def certificate(self,q,order,collect):
        return active_orientation(self.b,self.b.unpack(q),order,
                                  collect_activation=collect,exact_fallback=decide)

    def check_material_domain(self,q,u):
        result=self.certificate(q,u[2],False)
        if result["status"]=="active_fold" and self.overlap_sink is not None:
            self.stop_overlap("material or shape displacement at fixed transported order",q,u,result)
        if result["status"]=="unresolved":
            raise ArithmeticError(f"Geometry certification unresolved: {result['witness']}")
        return brief(result)

    def stop_overlap(self,kind,q,u,certificate,gradient=None):
        witness=certificate["witness"]
        record=dict(kind=kind,witness=witness,x=self.x.tolist(),cell_index_order="depth, radial, polar",
                    radius=self.shape.neck,exit=self.shape.exit,depth=self.shape.length)
        if self.overlap_sink is not None:
            # Save the triggering state before doing point diagnostics. This
            # also retains the exact fields of an infeasible search proposal.
            self.overlap_sink(self,record,q,u,gradient)
        raise VerifiedOverlap(record)

    def crystal_load_for_match(self,q,u):
        result=self.certificate(q,u[2],True)
        if result["status"]=="active_fold":
            self.stop_overlap("incoming transported order",q,u,result)
        if not result["all_active_cells_oriented"]:
            raise ArithmeticError(f"Invalid chart at entry: {brief(result)}")
        _,_,load,details=self.b.lower_bound_load(q,u[2])
        ids=np.union1d(result["potential_activation_nodes"],details["undefined_order_nodes"]).astype(int)
        self.affected=np.unique(self.number[ids])
        self.invalid_points=details["nonpositive_J_points"]
        return load

    def check_settled_fields(self,q,u,gradient):
        actual=self.certificate(q,u[2],True)
        if actual["status"]=="active_fold":
            self.stop_overlap("rematched order",q,u,actual,gradient)
        if not actual["all_active_cells_oriented"]:
            raise ArithmeticError(f"Invalid rematched chart: {brief(actual)}")
        candidates=np.union1d(self.affected,self.number[np.asarray(actual["potential_activation_nodes"],int)])
        c=u[2,self.first]
        g=self.pull.T@gradient[2]
        free=candidates[self.free[candidates]]
        # c>0 has already been certified everywhere, jointly over ALL positive
        # coefficients. Only zero coefficients with an adverse one-sided slope
        # need a second query about their ability to activate on the fold.
        adverse=free[(c[free]==0.) & (g[free]<0.)]
        if adverse.size:
            test=np.isin(self.number,adverse).astype(float)
            connection=self.certificate(q,test,False)
            if connection["status"]=="active_fold":
                raise ArithmeticError("Unresolved one-sided reordering work at zero order; "
                                      f"indicator witness={connection['witness']}, "
                                      f"adverse nodes={adverse.tolist()}, lower derivatives={g[adverse].tolist()}")
            if not connection["all_active_cells_oriented"]:
                raise ArithmeticError("Unable to certify one-sided order support")
        kkt=np.where(c==0,np.minimum(g,0.),np.where(c==1,np.maximum(g,0.),g))
        raw={name:norm(gradient[i,self.model.masks[i]]) for i,name in enumerate(("amplitude","phase"))}
        raw.update(order_free=norm(g[self.free]),order_kkt=norm(kkt[self.free]),
                   order_interior=norm(g[self.free & (c>0)&(c<1)]))
        return dict(raw_residual=raw, material_domain=brief(actual),
                    zero_order_check=dict(enclosure_nodes=len(candidates),invalid_quadrature_points=self.invalid_points,
                           positive_disjoint_nodes=candidates[c[candidates]>0].tolist(),
                           zero_adverse_disjoint_nodes=adverse.tolist()))


def local_order_forces(trial,u,gradient,point):
    """Report signed lower-bound nodal derivatives whose spatial hats are positive here."""
    b=trial.b
    axes=[b.buffers["field_"+s] for s in ("x","t","a")]
    basis=[]
    for axis,value in zip(axes,point):
        i=np.clip(np.searchsorted(axis,value,side="right")-1,0,len(axis)-2)
        f=(value-axis[i])/(axis[i+1]-axis[i])
        basis.append([(i,1.-f),(i+1,f)])
    indices={}
    for (i,wx),(j,wt),(k,wa) in product(*basis):
        hat=wx*wt*wa
        if hat>0.:
            idx=int(trial.number[b.ids[i,k,j]])
            indices[idx]=indices.get(idx,0.)+hat
    joined=trial.pull.T@gradient[2]
    return [dict(joined_node=idx,basis_weight_at_witness=float(weight),order=float(u[2,trial.first[idx]]),
                 lower_bound_dE_dchi=float(joined[idx]),minus_lower_bound_dE_dchi=float(-joined[idx]))
            for idx,weight in indices.items()]
