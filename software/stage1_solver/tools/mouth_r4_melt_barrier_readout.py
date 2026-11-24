"""Validate the R4 static melt-only barrier and its conditional active probe."""
import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from mouth_cp10_label_grid import get_cp10
from mouth_exact_checkpoint_control import digest, read
from mouth_r4_eulerian_corner_audit import point
from mouth_r4_eulerian_relax_labels import cell_order
from mouth_r4_refinement_review import _collar_reference


def measure(scratch, output):
    root=scratch/'r4_slaved_melt_barrier01'
    prior=scratch/'r4_eulerian_stage2_01'
    body,_,axes=get_cp10(scratch)
    result=dict(gauge='Physical-volume harmonic Dirichlet energy plus the cell volume divided by eight times the sum of d-1-log(d) at eight one-sided Q1 corners. Only exactly melted cells contribute; d=det(dxi/dp).',
        scope='Fixed active trace, followed by one fixed-field preconditioned active-trace trial if its inner residual passes. No shape update.',
        inner_tolerance_unscaled=3.e-8, cases={}, sources={})
    for src in ('mouth_r4_melt_barrier.cpp','mouth_r4_melt_barrier.py','mouth_r4_melt_barrier_readout.py'):
        result['sources'][src]=digest(Path(__file__).with_name(src))
    for case in ('primary','secondary'):
        guard=root/(case+'_guard')
        outcome=read(guard/'outcome.json')
        if not (outcome['limitsVerified'] and outcome['exitCode']==outcome['childOutcome']['exitCode']==0
                and outcome['childOutcome']['guardReason'] is None and int(read(guard/'effective-limits.json')['memory.swap.max'])==0):
            raise RuntimeError(f'Unverified guard {case}')
        path=root/case
        row=read(path/'readout.json'); baseline=read(prior/case/'label_stall.json')
        state=prior/case/baseline['file']
        if digest(state)!=row['sha256']:raise ValueError('Stall baseline changed')
        with np.load(state) as f, np.load(path/'fixed_active.npz') as update:
            order=f['field'][2,body.ids].transpose(0,2,1);labels=update['labels']
            active=cell_order(order); st=np.zeros(order.shape,bool)
            for v in np.ndindex(2,2,2):
                st[tuple(slice(k,k+n) for k,n in zip(v,active.shape))] |= active
            if not np.array_equal(labels[st],f['labels'][st]):raise ValueError('Active trace changed')
            witness=baseline['label']['last_rejection'];cell=witness['cell']
            local=witness.get('local',witness.get('chart_fractions'))
            with mp.workdps(65):
                before=point(body,axes,f['labels'],order,cell,local,_collar_reference)
                after=point(body,axes,labels,order,cell,local,_collar_reference)
        journal=[json.loads(line) for line in (path/'fixed_active_trace.jsonl').read_bytes().splitlines()]
        if any(line['last']['delta']>=0 for line in journal[1:]):raise ArithmeticError('Recorded gauge step was not descending')
        if after['max_vertex_chi']!=0 or before['max_vertex_chi']!=0:raise ValueError('Wrong support')
        result['cases'][case]=dict(row,initial_limiting_physical_Jinv=before['Jinv'],
            fixed_trace_limiting_physical_Jinv=after['Jinv'],fixed_trace_limiting_chart_det=after['chart_det'],
            fixed_trace_limiting_gradient=after['grad_xi_frobenius'],same_active_nodes=int(st.sum()),
            strictly_melted_nodes=int((~st).sum()),fixed_trace_state_sha256=digest(path/'fixed_active.npz'),
            accepted_fixed_trace_changes=len(journal)-1,after_corner=after)
    with output.open('x') as doc:
        json.dump(result,doc,allow_nan=False,indent=2);doc.write('\n')
    for name,case in result['cases'].items():
        print(name,case['initial_limiting_physical_Jinv'],case['fixed_trace_limiting_physical_Jinv'],
              case['inner']['residual'],case.get('outer_probe',{}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();measure(a.scratch.resolve(),a.output)
