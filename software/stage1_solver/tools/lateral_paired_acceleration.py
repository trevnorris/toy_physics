"""Target-only axial acceleration with a fixed draining companion in the lateral reservoir."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from lateral_mouth_acceleration import Translation
from lateral_mouth_radius import SpeedTranslation, Target, held_reading, steady
from lateral_reservoir_field import LateralSlab
from lateral_source_rate import verified_packet
from mouth_feed_weak_check import dx4


class PairedTranslation(SpeedTranslation):
    def __init__(self, target, field, currents, distance, unit_density):
        super().__init__(target, field, currents[-1])
        self.target = target
        slab = self.slab
        self.source_g = np.maximum(1-(slab.x[:, None]+distance)**2-slab.r**2, 0)**3
        self.source_shape = self.source_g/np.roll(unit_density, round(-distance/slab.dx), axis=0)
        self.source_rate = 0.25
        self.source_speed = currents[0]/self.n[-1]
        self.total_g = self.g+self.source_g
        self.normal_factor = -self.total_g*(1-1/self.amp) + 1j*(self.speed+self.source_speed)
        # displaced_normal from the parent translates only the target. The fixed
        # source can respond to density changes locally at constant total intake.

    def boundary(self, h):
        dh = h[-1]
        dn = 2*(self.face.conj()*dh).real
        mean = np.sum(self.slab.dv*self.source_speed*dn)/self.source_rate
        return (self.normal_factor*dh + self.face*(-self.total_g*dn/(2*self.amp**3)
                                                    - 1j*self.source_speed*mean))

    def boundary_check(self, tangent, displacement):
        # Differentiate the full nonlinear law at a displaced target and fixed source.
        def evaluate(eps):
            face = self.face+eps*tangent[-1]
            tg, ts, _, _ = self.target.profile(displacement*eps)
            ss = self.source_rate*self.source_shape/np.sum(self.slab.dv*self.source_shape*abs(face)**2)
            return (-(tg+self.source_g)*(1-1/abs(face)) + 1j*(ts+ss))*face
        eps = 1e-5/max(1,np.max(abs(tangent[-1])))
        observed = (evaluate(eps)-evaluate(-eps))/(2*eps)
        expected = self.boundary(tangent)+displacement*self.displaced_normal
        return np.max(abs(observed-expected))

    def load(self, h):
        # Translation.load reads the target footprint alone. Undo its fixed-
        # rate correction because the target is run at fixed mean normal speed.
        carry, stress, direct = Translation.load(self, h)
        dn = 2*(self.face.conj()*h[-1]).real
        dj = np.sum(self.slab.dv*self.speed*dn)
        vx = (self.face.conj()*dx4(self.face,self.slab.dx)).imag/self.n[-1]
        carry += dj/self.rate*np.sum(self.slab.dv*self.j*vx)
        return carry, stress, direct


def target_load(face, weight, speed, slab):
    """Signed nonlinear load in the +x direction on this target footprint."""
    n = abs(face)**2
    advected = speed*(face.conj()*dx4(face,slab.dx)).imag
    stress = -weight*(1-1/abs(face))*dx4(n,slab.dx)/2
    return np.sum((advected+stress)*slab.dv)


def turning_point_check(motion, p, h, omega, displacement):
    # At a turning point A=-omega^2 X and V=0. Reconstruct static and dynamic
    # target loads at the identical shifted opening, once for each sign of X.
    differences, static_loads = [], []
    for x in (-displacement,displacement):
        w, speed, _, _ = motion.target.profile(x)
        static = motion.face+x*p[-1]
        f_static = target_load(static,w,speed,motion.slab)
        f_dynamic = target_load(static-omega**2*x*h[-1],w,speed,motion.slab)
        differences.append((f_dynamic-f_static)/(omega**2*x))
        static_loads.append(f_static)
    return np.mean(differences), (static_loads[1]-static_loads[0])/(2*displacement)


def check_displaced_steadies(motion, p, distance, reference, epsilon):
    faces, loads = [], []
    for x in (-epsilon,epsilon):
        weight,speed,_,_ = motion.target.profile(x)
        field, _ = steady(motion.slab,SimpleNamespace(g=weight,speed=speed),distance,reference)
        faces.append(field[-1].copy())
        loads.append(target_load(faces[-1],weight,speed,motion.slab))
    derivative = (faces[1]-faces[0])/(2*epsilon)
    return (np.max(abs(derivative-p[-1]))/np.max(abs(p[-1])),
            (loads[1]-loads[0])/(2*epsilon))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet",type=Path)
    parser.add_argument("radius",type=float,choices=(0.75,1.,1.5))
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    archive = json.loads(Path("docs/_measurements/lateral_mouth_radius.json").read_text())
    record, = [r for r in archive["runs"] if r["radius"]==args.radius and str(args.packet) in r["reference_fields"]]
    hashes = {Path(k):v for k,v in record["reference_fields"].items()}
    with verified_packet(args.packet, hashes) as packet:
        x,r,dz = packet["x"],packet["r"],float(packet["dz"])
        slab = LateralSlab(round(r[-1]+r[0]), round(1/(2*r[0])), round(1/dz))
        if int(packet["d"]) or not np.array_equal(x,slab.x) or not np.array_equal(r,slab.r) or dz!=slab.dz:
            raise ValueError("Calibration geometry mismatch")
        ref = abs(packet["field"][-1])**2
        mean = np.sum(slab.dv*packet["currents"][0]/ref)/(4*np.pi/3)
    target = Target(slab,args.radius,ref,mean)
    args.output.mkdir(exist_ok=False)
    print(f"radius={args.radius:g} extent={r[-1]+r[0]:g} dx={slab.dx:.12g} dz={dz:.12g}",flush=True)
    for d in (4,8):
        field, currents = steady(slab,target,d,ref)
        parts = held_reading(slab,target,field,currents,d)
        attraction = parts[0]+parts[1]
        old_static = dict(w.split("=") for w in next(l for l in record["stdout"] if l.startswith(f"d={d} ")).split())
        if not np.isclose(attraction,float(old_static["load"]),rtol=2e-9,atol=1e-13):
            raise ArithmeticError("Static paired load differs from the radius preparation")
        motion = PairedTranslation(target,field,currents,d,ref)
        p,rp = motion.static_tangent()
        bp = motion.boundary_check(p,1)
        print(f"static_d={d} attraction={attraction:.12g} rp={rp:.4g} boundary_error={bp:.4g}",flush=True)
        local_x, local_r = abs(x)<10, r<4
        # Preserve the full paired base for future perturbations, and current
        # and tangent fields across both mouths for repeatable local readouts.
        fields = dict(field=field,currents=currents,x=x,r=r,dz=dz,local_x=x[local_x],local_r=r[local_r],
                      p=p[:,local_x,:][:,:,local_r],d=d,target_radius=args.radius)
        for omega in (0.,0.01):
            h,residual,steps = motion.dynamic_tangents(p,omega)
            carry,stress,direct = motion.load(h)
            reaction = -(carry+stress)
            old = dict(w.split("=") for w in next(l for l in record["stdout"] if l.startswith(f"omega={omega:g} ")).split())
            isolated, iso_direct = float(old["reaction"]),float(old["direct"])
            bh = motion.boundary_check(h,0)
            print(f"d={d} omega={omega:g} carry_per_A={carry:.12g} stress_per_A={stress:.12g} "
                  f"paired={reaction:.12g} isolated={isolated:.12g} change_pct={100*(reaction-isolated)/isolated:.9g} "
                  f"direct={-direct:.12g} isolated_direct={iso_direct:.12g} "
                  f"direct_change_pct={100*(-direct-iso_direct)/iso_direct:.9g} residual={residual:.4g} "
                  f"boundary_error={bh:.4g} steps={steps}",flush=True)
            if not np.isfinite([reaction,direct,bp,bh]).all() or max(bp,bh)>3e-8:
                raise ArithmeticError("Paired boundary tangent mismatch")
            if omega:
                eps = 0.002*slab.dx
                signed, slope = turning_point_check(motion,p,h,omega,eps)
                signed_half, slope_half = turning_point_check(motion,p,h,omega,eps/2)
                print(f"subtract_d={d} omega={omega:g} eps={eps:.7g} reaction_fd={signed:.12g} "
                      f"reaction_fd_half={signed_half:.12g} static_slope={slope_half:.12g}",flush=True)
                if max(abs(signed-reaction),abs(signed_half-reaction))>3e-5*abs(reaction)+1e-10:
                    raise ArithmeticError("Matched-position subtraction failed")
                if record["grid"]=="e24_c4" and args.radius==1 and d==4:
                    # The small static force gradient comes from differencing
                    # two larger loads. Resolve its centered-difference error.
                    samples = []
                    for delta in (eps/2,eps/4,eps/8):
                        error, slope_solve = check_displaced_steadies(motion,p,d,ref,delta)
                        _, slope_tangent = turning_point_check(motion,p,h,omega,delta)
                        print(f"displaced_d={d} shift={delta:.9g} tangent_relative_error={error:.6g} "
                              f"static_slope_solve={slope_solve:.12g} static_slope_tangent={slope_tangent:.12g}",
                              flush=True)
                        samples.append((error,slope_solve,slope_tangent))
                    extrapolated = (4*np.array(samples[-1][1:])-samples[-2][1:])/3
                    print(f"extrapolated_d={d} static_slope_solve={extrapolated[0]:.12g} "
                          f"static_slope_tangent={extrapolated[1]:.12g}",flush=True)
                    if (samples[-1][0]>3e-5 or
                            not np.isclose(*samples[-1][1:],rtol=2e-5,atol=1e-12) or
                            not np.isclose(*extrapolated,rtol=2e-5,atol=1e-12)):
                        raise ArithmeticError("Displaced steady limit disagrees with translation tangent")
            fields["h_slow" if omega==0 else "h_001"] = h[:,local_x,:][:,:,local_r]
        with (args.output/f"d{d}.npz").open("xb") as stream:
            np.savez_compressed(stream,**fields)
        del fields,motion,p,h,field


if __name__=="__main__":
    main()
