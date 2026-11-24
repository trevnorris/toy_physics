"""Joint R4 shape and material descent with one absolute, spatial order field.

Start with the accepted crystal of the requested solid mesh. The rematched held
field is an initial guess only. Every accepted geometry is rematched with that
geometry's crystal load and the same physical bath trace.
"""
import argparse
from pathlib import Path
from weakref import proxy

import numpy as np

from mouth_continuous_depth_mesh import axes_for_depth, neutral_force
from mouth_continuous_evidence import held_receipt
from mouth_continuous_order import ContinuousTrial
from mouth_exchange_carrier import BIAS
from mouth_joint_endpoint_audit import solid_readout
from mouth_joined_front_audit import Front
from mouth_last_rings_fluid import GRIDS
from mouth_rim_pinned_comparison import digest, read, validate_guard
from mouth_sealed_shape_relax import Meridian
from mouth_whole_r4_relax import Release, STEP, atomic_json, log


class ShapeTrial(ContinuousTrial):
    def measure(self, u, load):
        m = self.model
        m.material_load = load
        result=m.measure(u)
        result["disk_minus_feed_current"]=float(result["weak_disk_current"]-result["feed_current"])
        order=u[2,m.ids[1,0]]
        result["absolute_order_range"]=[float(u[2].min()),float(u[2].max())]
        result["axis_depth_order"]=order[0,0,:].tolist()
        result["wall_pole_depth_order"]=order[-1,0,:].tolist()
        result["wall_equator_depth_order"]=order[-1,-1,:].tolist()
        result["depth_knots"]=m.z[1].tolist()
        return result


class ContinuousRelease(Release):
    trial_type = ShapeTrial

    def __init__(self, args):
        self.args=args
        case=args.original/args.level
        validate_guard(case)
        held_receipt(args.held/"held")
        validate_guard(args.held/"readback")
        proof=read(args.edges)["results"][args.level]
        self.endpoint=read(case/"endpoint.json")
        prior=case/self.endpoint["file"]
        args.state=Path(read(case/"invocation.json")["command"][5])
        for path, key in ((case/"endpoint.json","endpoint_sha256"),(prior,"state_sha256"),
                          (args.state,"preparation_sha256")):
            if digest(path)!=proof[key]:
                raise ValueError(f"Changed preparation {path}")
        held=read(args.held/"saved/record.json")
        seed=args.held/"saved/wound.npz"
        if digest(seed)!=held["solutions"]["wound"]["sha256"] or held["state"]!=proof["preparation_sha256"]:
            raise ValueError("Mismatched continuous field or material preparation")
        self.base=Front(args.state,GRIDS["double"],1)
        self.discretization=proof["solid_discretization"]
        self.reference_depth=self.endpoint["depth"]
        with np.load(prior) as f, np.load(seed) as fields:
            self.q=f["q"].copy()
            self.x=f["x"].copy()
            self.u=fields["field"].copy()
        if BIAS!=.4510084282857874 or self.endpoint["pole_radius"][0]!=4.:
            raise ValueError("Mismatched radius or bath bias")
        self.trial=self.build(self.x,True)
        if self.u.shape!=self.trial.model.seed.shape:
            raise ValueError("Field charts have different topology")
        self.history=[]
        self.iteration=0
        self.shape_metric=None
        self.last=None
        log("input",level=args.level,initial_crystal=digest(prior),held_field=digest(seed),
            state_sha=digest(args.state),radius=4,winding=1,bias=BIAS,solid=self.discretization,
            reference_depth=self.reference_depth,absolute_order=True)

    def build(self,x,solve_fields=False):
        shape=Meridian(4.,self.base.sh.params+x[:4]*[1.,.25,.25,.25],
                       extra=self.base.sh.extra+.25*x[4])
        axes=axes_for_depth(shape.length,self.reference_depth)
        trial=getattr(self, "trial_type", ShapeTrial)(self.base,self.args.state,x,**self.discretization,
                         axes=axes,precondition=solve_fields)
        # Signed shapes include the full bath energy on both sides of the difference.
        trial.model.extended_energy_sum=True
        # Release uses trial.f for the field object. Here the two are identical.
        # Drop each signed probe's large sparse matrices as soon as it finishes.
        trial.f=proxy(trial)
        return trial

    def work(self,x,q,u,e,h):
        length=self.base.sh.length*np.exp(x[0])
        background=neutral_force(length,self.reference_depth,h,24)
        if max(abs(case["total"]) for case in background.values())>1.e-5:
            raise ArithmeticError(f"Neutral bath depth-force control exceeds 1e-5 at {length}: {background}")
        return super().work(x,q,u,e,h)

    def finish(self,t,q,u,load,measured):
        solid=solid_readout(self.args.state,t.shape,t.polar_shift,t.b,t.b.unpack(q))
        report=dict(self.last,solid=solid,current=measured,
                    neutral_control=neutral_force(t.shape.length,self.reference_depth,STEP/2,24))
        atomic_json(self.args.folder/"endpoint.json",report)
        if (measured["feed_current"]<=0 or measured["port_error"]>3e-5 or
                abs(measured["disk_minus_feed_current"])>3e-5 or
                not solid["orientation"]["all_open_cells_oriented"] or solid["all"]["minimum_J"]<=0):
            raise ArithmeticError("Shape-force endpoint failed orientation, current or volume check")
        log("mesh_endpoint",endpoint=digest(self.args.folder/"endpoint.json"),
            qualification="Review spatial order, loading, drainage and solid continuity before admission of another mesh")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("original","held","edges","folder"):
        parser.add_argument(key,type=Path)
    parser.add_argument("level",choices=("coarse","fine"))
    args=parser.parse_args()
    run=ContinuousRelease(args)
    try:
        run.run()
    except (ArithmeticError,ValueError,RuntimeError) as exc:
        run.failure(exc)
        raise
