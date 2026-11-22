"""Read the zero-drive common-face homogeneous preparations before attempting a junction."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.stats import t

from mouth_receiver_junction_geometry import local_order, shell_vectors
from mouth_receiver_material_screen import signature


PRESSURE_MARGIN = .25  # One previously probed native pressure step. Operational only.


def spatial_order(state, midpoint, ordered):
    """One-length slabs use the actual common face normal and eight shifted window centers."""
    box = np.asarray(state["box"])
    pos = np.asarray(state["pos"])
    height = box[3, 3]
    phase = (pos @ np.linalg.inv(box))[:,3] % 1.
    q = local_order(state)
    slabs = []
    for center in np.arange(8)/8:
        mask = np.abs((phase-center+.5)%1.-.5)*height <= .5
        values = q[mask]
        slabs.append(float(values.mean()) if len(values) and np.isfinite(values).all() else None)
    correct = sum(v is not None and (v>midpoint)==ordered for v in slabs)
    return slabs, correct


def contacts(state):
    n = len(state["pos"])
    make, retain = (np.zeros((n,n),bool) for _ in range(2))
    for i, ids, dx in shell_vectors(state):
        retain[i,ids] = True
        make[i,ids[np.sum(dx*dx,axis=1)<1.15**2]] = True
    np.fill_diagonal(make, False)
    np.fill_diagonal(retain, False)
    return make, retain


def signed_interval(block_vectors):
    """Describe finite-path scatter using four successive two-block means from the late half."""
    pairs=np.array(block_vectors).reshape(4,2,4).mean(axis=1)
    center=pairs.mean(axis=0)
    radius=t.ppf(.975,3)*pairs.std(axis=0,ddof=1)/2
    return dict(mean=center.tolist(), descriptive_t95=np.array([center-radius,center+radius]).T.tolist(),
                late_two_block_means=pairs.tolist(),
                within_operational_margin=bool(np.all(np.abs(center)+radius<=PRESSURE_MARGIN)))


def history(r, threshold, ordered):
    assert r["T"]==1.088 and r["Pext"]==0 and len(r["blocks"])==16
    first=r["prepared"]
    n=len(first["pos"])
    assert n==128 and len(r["origin"])==n and set(r["origin"])=={int(not ordered)}
    first_box=np.asarray(first["box"])
    assert np.max(np.abs(first_box[:3,3]))<1.e-12 and first_box[3,3]>0
    start_density=float(n/np.linalg.det(first_box))
    start_signature=signature(first) if start_density>=.3 else None
    start_slabs,start_correct=spatial_order(first,threshold,ordered) if start_signature else ([],0)
    initial_links,_=contacts(first)
    links=initial_links.copy()
    rows=[]
    tractions=[]
    for b in r["blocks"]:
        s=b["state"]
        bx=np.array(s["box"])
        np.testing.assert_allclose(bx[:3],first_box[:3], atol=1.e-10)
        rho=n/np.linalg.det(bx)
        assert rho>0
        new,retain=contacts(s)
        links=(links&retain)|new
        _,groups=connected_components(links,directed=False)
        cluster=int(max(np.bincount(groups)))
        p=-np.array(b["mean"]["stress_tensile"])[:,3]
        np.testing.assert_allclose(np.mean([m["pressure_w"] for m in b["profile"]],axis=0),p,atol=1.e-8)
        np.testing.assert_allclose(np.mean([m["rho"] for m in b["profile"]]), b["mean"]["density"],atol=1.e-8)
        tractions.append(p)
        # Dilute cells cannot pass the material test, and their reciprocal search is unbounded.
        fourier=signature(s) if rho>=.3 else None
        if fourier is None:
            slabs, correct = [], 0
        else:
            slabs, correct = spatial_order(s, threshold, ordered)
        expected=bool(fourier and ((min(fourier["four"])>=.2) if ordered else (fourier["peak"]<=.2)))
        rows.append(dict(time=b["t_end"],steps=b["steps"], density=b["mean"]["density"],endpoint_density=float(rho), cluster=cluster,
                         pressure_w=p.tolist(),global_four=fourier["four"] if fourier else None,
                         global_expected=expected, local_slabs=slabs, correct_slabs=correct))
    final=r["final"]
    assert final==r["blocks"][-1]["state"]
    dp=np.asarray(final["pos"])@np.linalg.inv(final["box"])@first_box-np.array(first["pos"])
    dp-=dp.mean(axis=0)
    rms=float(np.sqrt(np.mean(np.sum(dp*dp,axis=1))))
    lost=float(1-(initial_links & links).sum()/initial_links.sum()) if initial_links.any() else None
    late=rows[-8:]
    def density(blocks):
        return np.average([x["density"] for x in blocks],weights=[x["steps"] for x in blocks])
    density_drift=float(density(late[4:])/density(late[:4])-1)
    structural=all(x["density"]>=1.2 and x["cluster"]>=116 and x["global_expected"] and x["correct_slabs"]==8 for x in late)
    structural &= abs(density_drift)<=.10 and (ordered or (lost is not None and lost>=.5 and rms>=1.))
    return dict(material_pass=bool(structural),native_elapsed=r["measure"],warm_elapsed=r["warm"], final_rms=rms,
                post_warm_prepared=dict(density=start_density,global_four=start_signature["four"] if start_signature else None,
                                        local_slabs=start_slabs,correct_slabs=start_correct),
                absent_at_end=lost,density_drift=density_drift,pressure=signed_interval(tractions[-8:]),samples=rows)


def read_preflight(root, output):
    threshold=json.loads((root/"reference.json").read_bytes())
    assert threshold["ordered_increasing"] and threshold["separates_all_homogeneous_snapshot_means"]
    cases={}
    group_result={}
    pair_result={}
    names_by_group={}
    nominal=48.88171243122866
    for case in ("training","heldout"):
        for material in ("ordered","receiver"):
            group=case+"_"+material
            names_by_group[group]=[group+suffix for suffix in ("", "_rep2", "_rep3")]
            for name in names_by_group[group]:
                raw=(root/"runs"/(name+".json")).read_bytes()
                r=json.loads(raw)
                np.testing.assert_allclose([r["warm"],r["measure"]],[nominal,nominal],rtol=0,atol=1.e-8)
                cases[name]=history(r,threshold["training_midpoint"],material=="ordered")
                cases[name].update(sha256=sha256(raw).hexdigest(), seed=r["seed"], prepared_family=r["input_hash"])
            reps=[cases[k] for k in names_by_group[group]]
            expected_seeds={"training":[51001,51011,51023],"heldout":[61111,61121,61129]}
            assert [rep["seed"] for rep in reps]==expected_seeds[case] and len({rep["prepared_family"] for rep in reps})==1
            group_result[group]=replicate_interval([rep["pressure"]["mean"] for rep in reps])
        a,b=[names_by_group[case+"_"+k] for k in ("ordered","receiver")]
        # The ordered and receiver trajectories with the same replicate index form each contrast.
        deltas=[np.subtract(cases[i]["pressure"]["mean"],cases[j]["pressure"]["mean"]) for i,j in zip(a,b,strict=True)]
        pair_result[case]=replicate_interval(deltas)
    material=all(v["material_pass"] for v in cases.values())
    pressure=all(v["within_operational_margin"] for v in (*group_result.values(),*pair_result.values()))
    result=dict(cases=cases,conditional_traction_by_preparation=group_result,ordered_minus_receiver=pair_result,material_signatures_pass=material,
                all_conditional_pressure_intervals_within_margin=pressure,operational_pressure_margin=PRESSURE_MARGIN,
                interval_scope="Each group and paired difference uses three independent thermal realizations conditional on its initial geometry, with marginal two-sided 95% t intervals. The separate four time-block intervals describe temporal scatter. Neither set establishes general behavior over new founding geometries.",
                physical_scope="Finite homogeneous common-plane preparations. Joining displacements, coexistence and crossing remain unmeasured.")
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for key, record in cases.items():
        print(key,"material",record["material_pass"],"P_xw,yw,zw,ww",record["pressure"]["mean"],"descriptive intervals",record["pressure"]["descriptive_t95"])
    print("material criteria",material,"all finite-path component intervals in pressure margin",pressure)


def replicate_interval(readings):
    values=np.asarray(readings)
    assert values.shape==(3,4)
    mean=values.mean(axis=0)
    half=t.ppf(.975,2)*values.std(axis=0,ddof=1)/np.sqrt(3)
    return dict(mean=mean.tolist(),t95=np.array([mean-half,mean+half]).T.tolist(),replicate_readings=values.tolist(),
                within_operational_margin=bool(np.all(np.abs(mean)+half<=PRESSURE_MARGIN)))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("output",type=Path)
    cli=parser.parse_args()
    read_preflight(cli.root,cli.output)
