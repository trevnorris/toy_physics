"""Continue the unpinned 24x24 R4 run through the verified checkpoint-10 shape proposal."""
import argparse
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import numpy as np

from mouth_continuous_shape_relax import ContinuousRelease
from mouth_exact_melt_trial import ExactMeltTrial, local_order_forces
from mouth_exchange_carrier import BIAS
from mouth_free_shape_solid import pointer
from mouth_joined_front_audit import Front
from mouth_last_rings_fluid import GRIDS
from mouth_melted_collar_survival import Survival
from mouth_r4_fold_constitutive import Prism, native_callback, point
from mouth_rim_pinned_comparison import digest, read, validate_guard
from mouth_whole_r4_relax import atomic_json, save_fields


CHECKPOINT="51727c2d48eea199c29bacbfd4bf07b7b11a3e997a737ebcc8c2f69925cd112b"


class ExactRelease(Survival):
    trial_type=ExactMeltTrial

    def __init__(self,source,replay,controls,folder):
        self.args=SimpleNamespace(folder=folder)
        self.purpose="reconstructed full shape proposal"
        area=source/"coarse"
        old=read(area/"latest.json")
        if old["iteration"]!=10 or digest(area/old["file"])!=CHECKPOINT or old["sha256"]!=CHECKPOINT:
            raise ValueError("Expected unchanged checkpoint 10")
        previous=Path(read(area/"invocation.json")["command"][-3])/"coarse"
        cmd=read(previous/"invocation.json")["command"]
        original=Path(cmd[-5])/"coarse"
        proof=read(Path(cmd[-3]))["results"]["coarse"]
        self.args.state=Path(read(original/"invocation.json")["command"][5])
        if (digest(self.args.state)!=proof["preparation_sha256"] or
                digest(original/"endpoint.json")!=proof["endpoint_sha256"]):
            raise ValueError("Reference preparation changed")
        check=read(controls/"new.json")
        import mouth_exact_melt_trial, mouth_active_prism_orientation, mouth_free_shape_solid
        for key,module,suffix in (("candidate_sha",mouth_exact_melt_trial,".py"),
                                  ("native_interface_sha",mouth_free_shape_solid,".py"),
                                  ("native_solid_sha",mouth_free_shape_solid,".cpp"),
                                  ("exact_orientation_sha",mouth_active_prism_orientation,".py")):
            if check[key]!=digest(Path(module.__file__).with_suffix(suffix)):
                raise ValueError(f"Runtime differs from controlled source {key}")
        import mouth_prism_support_decision
        overlay=Path(mouth_exact_melt_trial.__file__).resolve().parent
        control_cmd=read(controls/"new_guard/invocation.json")["command"]
        if (Path(control_cmd[control_cmd.index("--candidate")+1]).resolve()!=overlay or
                Path(mouth_prism_support_decision.__file__).resolve().parent!=overlay or
                digest(overlay/"mouth_prism_support_decision.py")!=read(overlay.parent/"inputs.json")["patch/mouth_prism_support_decision.py"]):
            raise ValueError("Exact fallback differs from the verified overlay")
        for root in (controls/"old_guard",controls/"new_guard",replay/"serial_guard"):
            validate_guard(root)
        if (check["checkpoint"]!=CHECKPOINT or not all(check["bitwise_equal"].values()) or
                any(check["exact_differences"].values()) or
                check["negative_control"]["status"]!="active_fold" or BIAS!=.4510084282857874):
            raise ValueError("Checkpoint controls failed")
        self.discretization=proof["solid_discretization"]
        if self.discretization!=dict(n=8,k=24,angles=6,degree=3,subdivisions=2):
            raise ValueError("Wrong crystal mesh")
        self.base=Front(self.args.state,GRIDS["double"],1)
        self.reference_depth=read(original/"endpoint.json")["depth"]
        meta=read(replay/"replay.json")
        prop=read(replay/"proposal.json")
        if (not meta["identical_rejection"] or meta["replayed"]!=meta["archived"] or
                digest(replay/"trial.npz")!=meta["snapshot_sha256"] or
                digest(replay/"proposal.json")!=meta["proposal_sha256"] or prop["input_sha"]!=CHECKPOINT):
            raise ValueError("Full trial did not reproduce the archived rejection")
        with np.load(replay/"trial.npz",allow_pickle=False) as f:
            if int(f["radius"])!=4 or int(f["winding"])!=1 or str(f["state_sha"])!=digest(self.args.state):
                raise ValueError("Unexpected replayed material")
            self.x,self.q,self.u=(f[key].copy() for key in ("x","q","field"))
            self.energy=float(f["energy"])
            trial=self.trial=self.build(self.x,True)
            if not np.array_equal(self.x,prop["x"]) or not np.array_equal(trial.b.unpack(self.q),f["solid"]):
                raise ValueError("Trial geometry/material mismatch")
            if self.energy>old["energy"]+1.e-9:
                raise ValueError("Rematch would increase checkpoint energy")
            # Replay used exactly this lower-bound deposition. Recheck the
            # actual fields with the new certificate, then compare every load.
            cert=trial.certificate(self.q,self.u[2],False)
            if cert["status"]=="active_fold":
                trial.stop_overlap("reconstructed full-step rematched order",self.q,self.u,cert,f["gradient"])
            if not cert["all_active_cells_oriented"]:
                raise ArithmeticError("Unresolved geometry on replay")
            load=trial.crystal_load_for_match(self.q,self.u)+trial.ext
            if not np.array_equal(load,f["load"]):
                raise ArithmeticError("Refreezing work changed in the full-step replay")
            trial.model.material_load=load
            energy,grad=trial.model.evaluate(self.u)
            if energy!=self.energy or not np.array_equal(grad,f["gradient"]):
                raise ArithmeticError("Replay energy or field gradient differs")
            checks=trial.check_settled_fields(self.q,self.u,grad)
            flow=trial.measure(self.u,load)
            if (flow["feed_current"]!=meta["current"]["feed_current"] or
                    abs(flow["disk_minus_feed_current"])>3.e-5 or flow["port_error"]>3.e-5):
                raise ArithmeticError("Unsatisfied conserved-current gate on the full step")
            target=51603
            basis=trial.certificate(self.q,(trial.number==target).astype(float),False)
            if not basis["all_active_cells_oriented"] or self.u[2,trial.first[target]]!=meta["order"]:
                raise ArithmeticError("Guarded full-step comparison for node 51603 disagrees")
            result=dict(checkpoint_sha=CHECKPOINT,replayed_sha=meta["snapshot_sha256"],
                        archived_rejection=meta["replayed"],order=meta["order"],
                        joined_dE_dchi=meta["lower_bound_order_derivative"], node_basis=basis,
                        depth=trial.shape.length,exit=trial.shape.exit,energy=self.energy,
                        rematched_field=meta["field"],guard_checks=checks,
                        current=flow,forces_at_checkpoint=old["shape_force"])
            atomic_json(folder/"full_trial_gate.json",result)
        self.iteration=11
        self.history=[]
        self.shape_metric=np.asarray(prop["metric"])
        self.last=old
        self.parent_x=np.array(old["x"])
        self.parent_g=-np.array(old["shape_force"])

    def build(self,x,solve_fields=False):
        trial=ContinuousRelease.build(self,x,solve_fields)
        trial.overlap_sink=self.save_overlap
        return trial

    def material_update(self,*args):
        self.purpose="internal material line search and rematch"
        return super().material_update(*args)

    def shape_step(self,*args):
        self.purpose="global shape line search at fixed transported fields"
        return super().shape_step(*args)

    def shape_update(self,*args):
        self.purpose="global shape step and rematch"
        return super().shape_update(*args)

    def work(self,*args):
        self.purpose="signed virtual-work probe"
        try:
            return super().work(*args)
        finally:
            self.purpose="relaxation proposal"

    def save_overlap(self,trial,record,q,u,gradient):
        folder=self.args.folder
        path=folder/"overlap_state.npz"
        if path.exists():
            raise FileExistsError(path)
        record.update(proposal_type=self.purpose,checkpoint_sha=CHECKPOINT)
        save_fields(path,q=q,x=trial.x,field=u,shape_params=trial.shape.params,
                    shape_extra=trial.shape.extra,state_sha=digest(self.args.state),radius=4,winding=1)
        record["saved_state_sha256"]=digest(path)
        atomic_json(folder/"overlap.json",record)
        if gradient is None:
            # A true W is undefined for this attempted overlap. Calculate the
            # original best-case load W=0 only at invalid quadrature points;
            # all positive-J terms use the real W, independent of current χ.
            _,_,inner,_=trial.b.lower_bound_load(q,np.zeros_like(u[2]))
            trial.model.material_load=inner+trial.ext
            _,gradient=trial.model.evaluate(u)
        w=record["witness"]
        joined=trial.pull.T@gradient[2]
        save_fields(folder/"overlap_order_variation.npz",joined_dE_dchi=joined,number=trial.number)
        record["order_variations"]=local_order_forces(trial,u,gradient,w["placement_x_t_angle"])
        record["order_force_scope"]="Signed joined nodal dE/dchi and its negative with zero lower bound for undefined W. Basis weights are evaluated at the displayed witness. All nodal derivatives are saved. Order drive at an inverted point has no unique constitutive value."
        from mouth_r4_fold_order import physical
        record["physical"]=dict(r_w_theta=physical(trial.shape,trial.polar_shift,w["placement_x_t_angle"]).tolist(),
                                chi=w["order"])
        atomic_json(folder/"overlap.json",record)
        try:
            solid=trial.b.unpack(q)
            correction=np.ascontiguousarray(u[2,trial.b.ids])
            trial.b.data.offset=pointer(solid)
            trial.b.data.correction=pointer(correction)
            callback=native_callback(Path(__import__("mouth_free_shape_solid").__file__).parent,folder/"point_reader")
            with mp.workdps(70):
                pos=point(Prism(trial.b,solid,w["cell_t_x_angle"],w["triangle"]),w["barycentric"],
                          w["angular_fraction"],trial,u,callback)
            record["physical"].update(reference_x_t_theta=pos["reference_x_t_theta"],
                                      J=pos["chart_to_physical_J_scale"]*w["chart_determinant"],
                                      density=pos["n"],principal_stretches=pos["principal_stretches"])
        except (ArithmeticError,ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
            record["point_readout_issue"]=str(exc)
        atomic_json(folder/"overlap.json",record)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("source","replay","controls","output"):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    if any((a.output/f).exists() for f in ("full_trial_gate.json","overlap.json","endpoint.json","issue.json")):
        raise FileExistsError("Preserve existing gated run")
    a.output.mkdir(parents=True,exist_ok=True)
    try:
        run=ExactRelease(a.source.resolve(),a.replay.resolve(),a.controls.resolve(),a.output.resolve())
        run.purpose="relaxation proposal"
        run.run(matched_energy=run.energy,last_x=run.parent_x,last_g=run.parent_g)
    except Exception as exc:
        record=dict(error=type(exc).__name__,detail=str(exc))
        if (a.output/"latest.json").is_file():
            record["last"]=read(a.output/"latest.json")["file"]
        else:
            record["last"]="parent checkpoint 10"
        atomic_json(a.output/"issue.json",record)
        raise
