"""Signed axial load on the fixed companion per prescribed target acceleration."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from lateral_mouth_acceleration import Translation
from lateral_paired_acceleration import target_load
from lateral_source_rate import verified_packet


class CompanionTrace(Translation):
    """Read the fixed-rate mouth on the saved patch. Positive x points toward the target."""

    def __init__(self, slab, layers, current, weight, intake_shape):
        self.slab = slab
        self.f = layers
        self.face = layers[-1]
        self.amp = abs(self.face)
        self.n = abs(layers)**2
        self.j = current
        self.g = weight
        self.rate = 0.25
        self.speed = current/self.n[-1]
        self.intake_shape = intake_shape
        if (np.any(current[weight == 0] != 0) or abs(np.sum(current*slab.dv)-self.rate)>3e-14
                or np.max(abs(self.speed-self.normal_speed(self.face)))>2e-13):
            raise ValueError("Archived current disagrees with the companion intake law")

    def normal_speed(self, face):
        return self.rate*self.intake_shape/np.sum(self.slab.dv*self.intake_shape*abs(face)**2)

    def nonlinear_load(self, face):
        return target_load(face,self.g,self.normal_speed(face),self.slab)

    def acceleration_difference(self, h, step):
        return (self.nonlinear_load(self.face+step*h)-self.nonlinear_load(self.face-step*h))/(2*step)

    def matched_difference(self, p, h, omega, shift):
        # At each turning point subtract the static load with the target at the
        # identical position. The companion footprint itself never moves.
        rates = []
        for displacement in (-shift,shift):
            a = -omega**2*displacement
            static = self.face+displacement*p
            rates.append((self.nonlinear_load(static+a*h)-self.nonlinear_load(static))/a)
        return np.mean(rates)


def read_pair(record, reference, slab_x, slab_r, distance, old_rows):
    saved, = [Path(name) for name in record["fields"] if Path(name).name==f"d{distance}.npz"]
    hashes = {Path(k):v for k,v in record["fields"].items()}
    with verified_packet(saved,hashes) as data:
        full_x,full_r = data["x"],data["r"]
        x,r = data["local_x"],data["local_r"]
        dz = float(data["dz"])
        ix = np.flatnonzero(abs(full_x)<10)
        ir = np.flatnonzero(full_r<4)
        if (not np.array_equal(full_x,slab_x) or not np.array_equal(full_r,slab_r)
                or not np.array_equal(x,full_x[ix]) or not np.array_equal(r,full_r[ir])
                or int(data["d"])!=distance or float(data["target_radius"])!=record["radius"]):
            raise ValueError("Paired field geometry mismatch")
        crop = slice(ix[0],ix[-1]+1)
        radial = slice(0,len(ir))
        layers = data["field"][-4:,crop,radial].copy()
        current = data["currents"][0,crop,radial].copy()
        p = data["p"][-1]
        tangents = {0.:data["h_slow"][-4:], 0.01:data["h_001"][-4:]}
    dx = 2*r[0]
    slab = SimpleNamespace(dx=dx,dz=dz,dv=2*np.pi*r*dx**2)
    weight = np.maximum(1-(x[:,None]+distance)**2-r**2,0)**3
    if (not np.any(weight) or np.any(weight[:2]) or np.any(weight[-2:]) or np.any(weight[:,-1])
            or any(h.shape!=layers.shape for h in tangents.values())):
        raise ValueError("Insufficient companion footprint or derivative stencil")
    offset = round(-distance/dx)
    if not np.isclose(offset*dx,-distance):
        raise ValueError("Unaligned source calibration")
    shape = weight/np.roll(reference,offset,axis=0)[crop,radial]
    trace = CompanionTrace(slab,layers,current,weight,shape)
    static_load = trace.nonlinear_load(trace.face)
    print(f"radius={record['radius']:g} d={distance} dx={dx:.12g} dz={dz:.12g} "
          f"static_companion_load={static_load:.12g} companion_rate={np.sum(current*slab.dv):.12g}",flush=True)
    for omega,h in tangents.items():
        carry,stress,direct = trace.load(h)
        total = carry+stress
        # Verify the load derivative from the nonlinear fixed-rate law for both
        # signs of the drive, on the same stationary footprint.
        step = 1e-4/max(1.,np.max(abs(h[-1])))
        fd = trace.acceleration_difference(h[-1],step)
        fd_half = trace.acceleration_difference(h[-1],step/2)
        self_row, = [row for row in old_rows if row["dist"]==distance and row["omega"]==omega]
        self_reaction = self_row["paired"]
        print(f"d={distance} radius={record['radius']:g} omega={omega:g} carry_per_A={carry:.12g} "
              f"stress_per_A={stress:.12g} cross_per_A={total:.12g} direct={direct:.12g} "
              f"fd={fd:.12g} fd_half={fd_half:.12g} target_reaction={self_reaction:.12g} "
              f"signed_cross_over_self={total/self_reaction:.12g}",flush=True)
        if not np.isfinite([total,direct,fd,fd_half]).all() or not np.allclose([fd,fd_half],total,rtol=2e-6,atol=2e-11):
            raise ArithmeticError("Companion fixed-rate load derivative failed")
        if omega:
            shift = 0.002*dx
            matched = [trace.matched_difference(p,h[-1],omega,delta) for delta in (shift,shift/2)]
            print(f"subtract_d={distance} radius={record['radius']:g} shift={shift:.8g} "
                  f"matched={matched[0]:.12g} matched_half={matched[1]:.12g}",flush=True)
            if not np.allclose(matched,total,rtol=5e-5,atol=4e-10):
                raise ArithmeticError("Companion turning-point subtraction failed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("grid",choices=("e24_c4","e24_c6","e48_c4"))
    args = parser.parse_args()
    archive = json.loads(Path("docs/_measurements/lateral_paired_acceleration.json").read_text())
    records = [row for row in archive["runs"] if row["grid"]==args.grid]
    if sorted(row["radius"] for row in records)!=[0.75,1,1.5]:
        raise ValueError("Incomplete or repeated radius preparation")
    cal, = {Path(path) for row in records for path in row["reference_fields"]}
    hashes = {cal:records[0]["reference_fields"][str(cal)]}
    with verified_packet(cal,hashes) as data:
        if int(data["d"])!=0:
            raise ValueError("Wrong calibration field")
        reference = abs(data["field"][-1])**2
        x,r = data["x"],data["r"]
    for row in records:
        if (row["reference_fields"][str(cal)]!=hashes[cal] or not row["limits_verified"]
                or row["worker_exit"]!=0 or row["guard_exit"]!=0):
            raise ValueError("Unverified paired preparation")
        old_rows = [v for v in archive["comparisons"] if v["grid"]==args.grid and v["radius"]==row["radius"]]
        for distance in (4,8):
            read_pair(row,reference,x,r,distance,old_rows)


if __name__=="__main__":
    main()
