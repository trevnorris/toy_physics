"""Move the fixed-density lateral reservoir while keeping the two intake laws fixed."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from lateral_companion_load import CompanionTrace
from lateral_mouth_radius import Target, held_reading, steady
from lateral_paired_acceleration import PairedTranslation
from lateral_reservoir_field import LateralSlab
from lateral_source_rate import verified_packet


def preparation(slab):
    archive = json.loads(Path("docs/_measurements/lateral_companion_load.json").read_text())
    run, = [row for row in archive["runs"] if row["grid"]=="e48_c4"]
    packet, = [Path(name) for name in run["reference_fields"] if Path(name).name=="d0.npz"]
    with verified_packet(packet,{packet:run["reference_fields"][str(packet)]}) as data:
        x,r,dz = data["x"],data["r"],float(data["dz"])
        if int(data["d"])!=0 or dz!=slab.dz or 2*r[0]!=slab.dx or r[-1]+r[0]!=48:
            raise ValueError("Wrong fixed reference calibration")
        density = abs(data["field"][-1])**2
        dv = 2*np.pi*r*slab.dx**2
        mean = np.sum(dv*data["currents"][0]/density)/(4*np.pi/3)
    # Only the two mouth footprints sample this calibration. Embed the archived
    # profile so that changing the reservoir cannot retune either opening.
    crop = slice(int(np.searchsorted(slab.x,x[0])),int(np.searchsorted(slab.x,x[-1]))+1)
    if not np.array_equal(slab.x[crop],x) or not np.array_equal(slab.r[:len(r)],r):
        raise ValueError("Larger domain does not contain the original cell centers")
    reference = np.ones((len(slab.x),len(slab.r)))
    reference[crop,:len(r)] = density
    target = Target(slab,1.5,reference,mean)
    original = Target(SimpleNamespace(x=x,r=r,dv=dv),1.5,density,mean)
    error = max(np.max(abs(new[crop,:len(r)]-old)) for new,old in
                ((target.g,original.g),(target.speed,original.speed),(target.dspeed,original.dspeed)))
    if not np.isfinite(error) or error>1e-12:
        raise ArithmeticError("Target intake profile changed")
    print(f"calibration_L=48 mean_speed={mean:.12g} max_profile_error={error:.6g}",flush=True)
    return target,reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("extent",type=int)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    if args.extent<=48:
        parser.error("Choose a reservoir farther away than the archived L=48 case")
    slab = LateralSlab(args.extent,4,16)
    target,reference = preparation(slab)
    distance = 4
    base,currents = steady(slab,target,distance,reference)
    held = held_reading(slab,target,base,currents,distance)
    motion = PairedTranslation(target,base,currents,distance,reference)
    p,rp = motion.static_tangent()
    bp = motion.boundary_check(p,1)
    ix = slice(np.searchsorted(slab.x,-10.),np.searchsorted(slab.x,10.))
    ir = slice(0,np.searchsorted(slab.r,4.))
    trace_grid = SimpleNamespace(dx=slab.dx,dz=slab.dz,dv=slab.dv[ir])
    trace = CompanionTrace(trace_grid,base[-4:,ix,ir],currents[0,ix,ir],
                           motion.source_g[ix,ir],motion.source_shape[ix,ir])
    print(f"L={args.extent} d={distance} target_radius={target.radius:g} dx={slab.dx:g} dz={slab.dz:g} "
          f"companion_rate={np.sum(slab.dv*currents[0]):.12g} target_static={sum(held[:2]):.12g} "
          f"companion_static={trace.nonlinear_load(trace.face):.12g} rp={rp:.6g} bp={bp:.6g}",flush=True)
    args.output.mkdir(exist_ok=False)
    with (args.output/"base.npz").open("xb") as stream:
        np.savez_compressed(stream,x=slab.x,r=slab.r,dz=slab.dz,field=base,currents=currents,
                            local_x=slab.x[ix],local_r=slab.r[ir],p=p[:,ix,ir],d=distance,target_radius=target.radius)
    for omega in (0.,0.01):
        h,residual,steps = motion.dynamic_tangents(p,omega)
        bh = motion.boundary_check(h,0)
        local_h = h[:,ix,ir]
        carry,stress,direct = trace.load(local_h[-4:])
        cross = carry+stress
        step = 1e-4/max(1.,np.max(abs(local_h[-1])))
        fd,fd_half = [trace.acceleration_difference(local_h[-1],a) for a in (step,step/2)]
        target_carry,target_stress,_ = motion.load(h)
        print(f"omega={omega:g} carry_per_A={carry:.12g} stress_per_A={stress:.12g} cross_per_A={cross:.12g} "
              f"direct={direct:.12g} fd={fd:.12g} fd_half={fd_half:.12g} "
              f"target_reaction={-(target_carry+target_stress):.12g} residual={residual:.6g} "
              f"bh={bh:.6g} steps={steps}",flush=True)
        if (not np.isfinite([bh,bp,cross,direct,residual]).all() or max(bh,bp)>3e-8
                or not np.allclose([fd,fd_half],cross,rtol=2e-6,atol=2e-11)):
            raise ArithmeticError("Companion load derivative or paired boundary tangent failed")
        if omega:
            shift = 0.002*slab.dx
            matched = [trace.matched_difference(p[-1,ix,ir],local_h[-1],omega,s) for s in (shift,shift/2)]
            print(f"shift={shift:g} matched={matched[0]:.12g} matched_half={matched[1]:.12g}",flush=True)
            if not np.allclose(matched,cross,rtol=5e-5,atol=4e-10):
                raise ArithmeticError("Matched-position subtraction failed")
        with (args.output/f"h_{omega:g}.npz").open("xb") as stream:
            np.savez_compressed(stream,h=local_h,omega=omega)
        del local_h,h


if __name__=="__main__":
    main()
