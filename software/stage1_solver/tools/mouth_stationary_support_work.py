"""Pair registered-collar work with the same shape displacements in the field records."""
import hashlib
import json
from math import fsum
from pathlib import Path


DATA = Path(__file__).resolve().parents[3] / "research/pde_ledger_v3/_measurements"
INTERVALS = ((-.002, -.001), (.001, .002), (-.001, .001))


def readings(source, case):
    r = source["runs"][case]
    raw = r["literal_stdout"]
    digest = r["stdout_sha256"] if "stdout_sha256" in r else r["sha256_stdout"]
    if isinstance(raw, list):
        raw = "\n".join(raw)+"\n"
    if hashlib.sha256(raw.encode()).hexdigest() != digest:
        raise ValueError(f"Source checksum failed in {case}")
    return [json.loads(line) for line in raw.splitlines()]


def main():
    crystal = json.loads((DATA/"mouth_finite_bore_support.json").read_text())
    simple = json.loads((DATA/"mouth_shallow_puncture_fields.json").read_text())
    deep = json.loads((DATA/"mouth_unified_strip_hold.json").read_text())
    old = readings(deep, "shape")
    for radius in (1,2):
        collar = readings(crystal, f"r{radius}g2")
        shallow = readings(simple, f"r{radius}_m72")
        static_energy = {depth: fsum(r["energy"] for r in collar if r["kind"]=="central"
                                     and r["depth"]==depth) for depth in ("long","unit")}
        endpoint, = [r for r in shallow if r["kind"]=="center" and abs(r["depth"]-1)<1e-10]
        field_delta = fsum(endpoint["energy_minus_deep"].values())
        print(json.dumps(dict(radius=radius, static_elastic_energies=static_energy,
                              fixed_chart_noncrystal_change=field_delta, comparison_energy_change=
                              static_energy["unit"]-static_energy["long"]+field_delta)))
        for axis in ("radial","axial"):
            base = [r for r in old if r["case"]==f"r{radius}_broad_72" and r["axis"]==axis]
            if len(base) != 2:
                raise ValueError("Expected two signed step sizes")
            states = {}
            for r in base:
                for sign, leg in zip((-1,1), r["legs"], strict=True):
                    offset = sign*r["offset"]
                    states[offset] = (leg["depth"] if axis=="axial" else offset,
                                      [leg[k] for k in ("order_inside","order_outside",
                                                       "carrier_two_sided","stored_wave")])
            for depth in ("long","unit"):
                c, = [r for r in collar if r["kind"]=="forces" and r["axis"]==axis and r["depth"]==depth]
                endpoints = [r[0] for r in c["rows"]]
                for i,(lo,hi) in enumerate(INTERVALS):
                    if depth=="long":
                        xl,el = states[lo]
                        xh,eh = states[hi]
                        passive = fsum(a-b for a,b in zip(el,eh))/(xh-xl)
                    else:
                        phase, = [r for r in shallow if r["kind"]=="secant" and r["axis"]==axis
                                  and r["interval"]==[lo,hi]]
                        xl,xh = phase["coordinate"]
                        passive = phase["force_before_static_crystal"]
                    j,k = ((0,1),(2,3),(1,2))[i]
                    if abs(endpoints[j]-xl)+abs(endpoints[k]-xh)>1e-11:
                        raise ValueError("Work refers to different physical displacements")
                    print(json.dumps(dict(radius=radius,depth=depth,axis=axis,interval=[lo,hi],
                                          collar=c["total"][i], fields=passive,
                                          comparison_sum=passive+c["total"][i]),allow_nan=False))


if __name__=="__main__":
    main()
