"""Local axial momentum budgets for the slow translation of an isolated draining mouth."""

import argparse
import json
from pathlib import Path

import numpy as np

from lateral_mouth_radius import SpeedTranslation, Target, held_reading, steady
from lateral_reservoir_field import LateralSlab
from lateral_source_rate import verified_packet


def slow_tangent(experiment, rhs):
    """Solve one of L u = i p and L h = i u, retaining both full fields."""
    tangent = np.zeros_like(rhs)
    for _ in range(500):
        update = experiment.step(tangent, rhs)
        change = np.max(abs(update-tangent))
        tangent = update
        if not np.isfinite(change):
            raise ArithmeticError("Nonfinite translation tangent")
        if change < 2e-13:
            residual = experiment.residual(tangent, rhs)
            if residual > 3e-9:
                raise ArithmeticError(f"Translation equation residual {residual}")
            return tangent, residual
    raise ArithmeticError("Translation tangent did not converge")


def cross(a, b, da, db):
    """Variation of the complex bilinear a* b."""
    return da.conj()*b + a.conj()*db


def link_budget(f, u, h, normal, dnormal, x, r, dx, dz, reach, stride):
    # Sum x-link momentum on the discrete carrier grid. The side and cap
    # currents are read from face fields; the storage comes from the velocity
    # tangent u, and the stress and all face readouts from acceleration tangent h.
    m = stride
    scale = m*dx
    link_x = (x[:-m]+x[m:])/2
    # Use integer half-cell positions to put a half weight on an endcap link.
    twice_x = np.rint(2*link_x/dx).astype(int)
    edge = round(2*reach/dx)
    wx = (abs(twice_x)<edge).astype(float) + 0.5*(abs(twice_x)==edge)
    nr = round(reach/dx)
    if not np.isclose(nr*dx,reach) or not 0 < nr < len(r) or wx[0] or wx[-1]:
        raise ValueError("Cylinder does not fit the local field")
    wz = np.full(len(f),dz)
    wz[[0,-1]] *= 0.5
    disk = 2*np.pi*r[:nr]*dx
    vol = wz[:,None,None]*wx[None,:,None]*dx*disk
    # Difference of adjacent control weights selects the two outward cap fluxes.
    caps = wz[:,None,None]*(wx[:-1]-wx[1:])[None,:,None]*disk
    left, right = f[:,:-m], f[:,m:]
    hl, hr = h[:,:-m], h[:,m:]
    stored = np.sum(vol*cross(left,right,u[:,:-m],u[:,m:]).imag[:,:,:nr])/scale
    n = abs(f)**2
    dn = 2*(f.conj()*h).real
    dp = n**4*dn            # P(n) = n**5/5
    dv = n**3*dn            # V(n) = (n**4-1)/4
    v = (n**4-1)/4
    links = len(link_x)
    pressure = sum(dp[:,k:k+links-1,:nr] for k in range(1,m+1))/m
    cap_pressure = np.sum(caps*pressure)
    kinetic = (cross(f[:,1:-m],right[:,:-1],h[:,1:-m],hr[:,:-1])
               - cross(left[:,:-1],right[:,1:],hl[:,:-1],hr[:,1:])).real/(2*dx*scale)
    cap_gradient = np.sum(caps*kinetic[:,:,:nr])
    mixed = (cross(left[:,:,nr],right[:,:,nr-1],hl[:,:,nr],hr[:,:,nr-1])
             - cross(left[:,:,nr-1],right[:,:,nr],hl[:,:,nr-1],hr[:,:,nr])).real/(2*dx*scale)
    side = np.sum(wz[:,None]*wx[None,:]*mixed)*(2*np.pi*reach*dx)
    mouth = (cross(normal[:-m],right[-1],dnormal[:-m],hr[-1])
             - cross(left[-1],normal[m:],hl[-1],dnormal[m:])).real/(2*scale)
    mouth_out = np.sum(wx[:,None]*dx*disk*mouth[:,:nr])
    # Finite differences of V and P do not satisfy n dV=dP exactly.
    # Evaluate the resulting body term directly on every link.
    defect = ((dp[:,m:]-dp[:,:-m]) -
              cross(left,right,hl,hr).real*(v[:,m:]-v[:,:-m]) -
              (left.conj()*right).real*(dv[:,m:]-dv[:,:-m]))/scale
    body = np.sum(vol*defect[:,:,:nr])
    return np.array([stored,cap_pressure,cap_gradient,side,mouth_out,body])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet",type=Path)
    parser.add_argument("radius",type=float,choices=(0.75,1.,1.5))
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    prior = json.loads(Path("docs/_measurements/lateral_mouth_radius.json").read_text())
    record, = [row for row in prior["runs"] if row["radius"]==args.radius
               and str(args.packet) in row["reference_fields"]]
    hashes = {Path(k):v for k,v in record["reference_fields"].items()}
    with verified_packet(args.packet,hashes) as packet:
        x,r,dz = packet["x"],packet["r"],float(packet["dz"])
        s = LateralSlab(round(r[-1]+r[0]),round(1/(2*r[0])),round(1/dz))
        if int(packet["d"]) or not np.array_equal(s.x,x) or not np.array_equal(s.r,r) or dz!=s.dz:
            raise ValueError("Mismatched radius-one calibration")
        unit,unit_j = packet["field"],packet["currents"]
    ref = abs(unit[-1])**2
    mean = np.sum(s.dv*unit_j[0]/ref)/(4*np.pi/3)
    target = Target(s,args.radius,ref,mean)
    if args.radius==1:
        base,currents = unit,unit_j
    else:
        del unit
        base,currents = steady(s,target,0,ref)
    isolated = held_reading(s,target,base,currents,0)
    if max(abs(isolated))>1e-11:
        raise ArithmeticError("Isolated static load failed symmetry check")
    motion = SpeedTranslation(target,base,currents[0])
    p,rp = motion.static_tangent()
    u,ru = slow_tangent(motion,1j*p)
    h,rh = slow_tangent(motion,1j*u)
    carry,stress,direct = motion.load(h)
    reaction = -(carry+stress)
    archived = next(l for l in record["stdout"] if l.startswith("omega=0 "))
    old = float(dict(w.split("=") for w in archived.split())["reaction"])
    print(f"radius={args.radius:g} extent={r[-1]+r[0]:g} dx={s.dx:.12g} dz={dz:.12g} "
          f"reaction={reaction:.12g} prior={old:.12g} direct={-direct:.12g} "
          f"rp={rp:.5g} ru={ru:.5g} rh={rh:.5g}",flush=True)
    if not np.isfinite(reaction) or not np.isclose(reaction,old,rtol=2e-8,atol=1e-13):
        raise ArithmeticError("Isolated preparation differs from the radius experiment")
    # Store the local fields needed for a new cut or a separate flux reading.
    cut_x = np.flatnonzero(abs(x)<10)
    cut = slice(cut_x[0],cut_x[-1]+1)
    radial = slice(0,np.count_nonzero(r<10))
    f = base[:,cut,radial]
    velocity = u[:,cut,radial]
    acceleration = h[:,cut,radial]
    normal = (motion.normal_factor*motion.face)[cut,radial]
    dnormal = motion.boundary(h)[cut,radial]
    with args.output.open("xb") as stream:
        np.savez_compressed(stream,x=x[cut],r=r[radial],dz=dz,field=f,u=velocity,h=acceleration,
                            normal=normal,dnormal=dnormal,radius=args.radius)
    for reach in (2,4,8):
        terms = [link_budget(f,velocity,acceleration,normal,dnormal,x[cut],r[radial],s.dx,dz,reach,m)
                 for m in (1,2)]
        for order,budget in ((2,terms[0]),(4,(4*terms[0]-terms[1])/3)):
            stored,pressure,kinetic,side,mouth,body = budget
            lateral = pressure+kinetic+side
            error = stored+lateral+mouth-body
            print(f"reach={reach} order={order} storage={stored:.12g} cap_pressure={pressure:.12g} "
                  f"cap_gradient={kinetic:.12g} side={side:.12g} lateral={lateral:.12g} "
                  f"mouth_out={mouth:.12g} lattice_body={body:.12g} closure={error:.7g} "
                  f"storage_over_reaction={stored/reaction:.9g} lateral_over_reaction={lateral/reaction:.9g} "
                  f"physical_gap_over_reaction={(stored+lateral-reaction)/reaction:.9g}",flush=True)
            if not np.isfinite(budget).all() or abs(error)>2e-8*max(abs(stored),abs(lateral),reaction)+2e-10:
                raise ArithmeticError(f"Discrete axial-momentum balance failed at reach {reach}, order {order}")


if __name__=="__main__":
    main()
