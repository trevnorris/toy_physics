#!/usr/bin/env python3
"""Compare all regenerated b cases with independent action and saved sources."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import sys
import time
STARTED=time.monotonic()
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'_measurements')]
import sympy as sp
import S11c_b_brane_operator_sympy_audit as b
import S11c_d_mixing_scattering_sympy_audit as d
from ledger_fold import _restore
from S11c_inertia_artifact_audit import export_data,locate
from S11c_d_end_pairing_emit import small_literal


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def assoc(value):return {str(k):v for k,v in value}
def mechanical(body):
    return sp.ImmutableMatrix([*b.named_tuple_row(b.named_tuple_row(body,'U_BODY_BALANCE'),'EXPANDED'),
        b.named_tuple_row(b.named_tuple_row(body,'E_W_BALANCE'),'EXPANDED')])


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--baseline',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True)
    args=parser.parse_args();base=args.run_directory;base.mkdir(parents=True,exist_ok=False)
    paths=[ROOT/'scripts/S11c_b_exports.py',ROOT/'scripts/S11c_b_brane_operator_sympy_audit.py',
        ROOT/'scripts/S11c_a_exports.py',ROOT/'directives/S11c_b_SHARED_PHYSICS.md',Path(__file__),
        args.baseline/'scripts/S11c_b_exports.py',args.baseline/'scripts/S11c_b_brane_operator_sympy_audit.py']
    baseline_export,baseline_script=paths[-2:]
    paths += [ROOT/name for name in ('scripts/ledger_fold.py','scripts/S11c_d_mixing_scattering_sympy_audit.py',
        'scripts/S11c_d_output_codec.py','_measurements/S11c_inertia_artifact_audit.py',
        '_measurements/S11c_d_end_pairing_emit.py','_measurements/S11c_d_end_pairing_check.py',
        '_measurements/S11c_d_modal_current_check.py','_measurements/S11c_d_joint_sheet_check.py')]
    pins={str(p):sha(p) for p in paths}
    snapshots={}
    for i,path in enumerate(paths):
        relative=Path('source')/str(i)/path.name;target=base/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
        snapshots[str(path)]=str(relative)
    old_values,_,_=export_data(baseline_export);new_values,export_pins,_=export_data(paths[0])
    before=dict(_restore(old_values['slab_operator']));after=dict(_restore(new_values['slab_operator']))
    tree=ast.parse(baseline_script.read_text())
    helper=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=='kinetic_balance_from_energy')
    namespace=dict(vars(b));exec(compile(ast.Module(body=[helper],type_ignores=[]),str(baseline_script),'exec'),namespace)
    baseline_inertia=namespace['kinetic_balance_from_energy']
    lam=sp.Symbol('s11cThicknessCoordinateExportLambda');records=[];residuals={};keys=set()
    def output(name,value,units=None):
        body=d.cas(value);key='s11cThicknessCoordinateB'+name
        if key in keys or key in b.INCOMING_LEDGER:raise ValueError('write-key collision')
        keys.add(key);metadata=[]
        for path,leaf in d.leaves(body):
            if isinstance(leaf,d.Str):continue
            unit=units(path) if callable(units) else units
            if unit is None:unit=tuple(b.dimension_of(leaf))
            expression=leaf.xreplace(b.PROFILE_GRADE_SUBS)
            polynomial=sp.Poly(sp.expand(expression),b.epsilon,b.eta_bg,b.sigma_W)
            support=sp.Poly(sp.expand(expression.subs({b.eta_bg:lam,b.sigma_W:lam*b.W0/b.L_W})),b.epsilon,lam)
            metadata.append({'path':path,'dimensionLTM':unit,
                'multigrade':sorted(g for g,c in polynomial.terms() if c!=0),
                'epsilonLambdaSupport':sorted(g for g,c in support.terms() if c!=0)})
        heavy=not small_literal(body)
        record={'writeKey':key,'value':d.carrier_fingerprint(body) if heavy else body,
                'representation':'CARRIER_PIT_SHA' if heavy else 'LITERAL','metadata':metadata}
        print(json.dumps({k:sp.srepr(d.cas(v)) for k,v in record.items()}),flush=True)
        records.append({'key':key,'body':body,'record':record})
        if name.endswith('Residual'):
            values=[v for _,v in d.leaves(body) if not isinstance(v,d.Str)]
            residuals[name]={'scalars':len(values),'nonzero':sum(v!=0 for v in values)}
    output('SourcePins',pins)
    output('ExportPinResidual',{name:sp.Integer(sha(locate(name))!=digest) for name,digest in export_pins.items()},tuple(b.DIM_ZERO))
    changed=sorted(k for k in new_values.keys()&old_values.keys() if new_values[k]!=old_values[k])
    output('ChangedValueSerializations',changed)
    output('KeySetDifference',(sorted(new_values.keys()-old_values.keys()),sorted(old_values.keys()-new_values.keys())))
    t=sp.Symbol('s11cThicknessCoordinateActionTime',real=True)
    fields=tuple(sp.Function('s11cThicknessCoordinateActionField'+str(i))(t) for i in range(4))
    accelerations=(*b.u_tt,b.e_tt);reverse={sp.diff(f,t,2):a for f,a in zip(fields,accelerations)}
    zero=dict.fromkeys(accelerations,sp.S.Zero)
    units=lambda path:tuple(b.DIM_BODY_U if path[0]<3 else b.DIM_ENERGY)
    preserved={}
    for axes,payload in after.items():
        case=tuple(map(str,axes));label=''.join(part.title().replace('_','') for part in case)
        prior=b.named_tuple_row(before[axes],'VALUE');current=b.named_tuple_row(payload,'VALUE')
        rho=b.density_pair(case[1])[1]
        # Independent supplied action in the defining physical coordinate.
        action=b.epsilon**2*(rho*sum(sp.diff(f,t)**2 for f in fields[:3])+b.mu_W*sp.diff(b.W0*fields[3],t)**2)/2
        action_rows=b.retained_grade(sp.ImmutableMatrix([sp.diff(sp.diff(action,sp.diff(f,t)),t).xreplace(reverse)/b.epsilon for f in fields]))
        old_source=b.retained_grade(sp.ImmutableMatrix([*baseline_inertia(rho)[0],baseline_inertia(rho)[1]]))
        old_rows,new_rows=mechanical(prior),mechanical(current)
        extracted=new_rows.applyfunc(lambda row:sp.expand(sum(sp.diff(row,a)*a for a in accelerations)))
        old_extracted=old_rows.applyfunc(lambda row:sp.expand(sum(sp.diff(row,a)*a for a in accelerations)))
        expected_delta=action_rows-old_source
        for name,value in [('BeforeRows',old_rows),('AfterRows',new_rows),('NativeInertia',extracted),
            ('BaselineInertia',old_extracted),('ActionInertia',action_rows),('BaselineSourceInertia',old_source),
            ('RowDelta',new_rows-old_rows),('ActionDerivedDelta',expected_delta),
            ('ActionResidual',extracted-action_rows),('BaselineSourceResidual',old_extracted-old_source),
            ('NonkineticPreservationResidual',new_rows.xreplace(zero)-old_rows.xreplace(zero)),
            ('FullDeltaAccountingResidual',new_rows-old_rows-expected_delta)]:
            output(label+name,value.applyfunc(sp.expand),units)
        slots=[]
        for key,value in current:
            other=b.named_tuple_row(prior,str(key))
            members=((str(slot),obj,b.named_tuple_row(other,str(slot))) for slot,obj in value if str(slot)!='EXPANDED') if str(key) in ('U_BODY_BALANCE','E_W_BALANCE') else (('',value,other),)
            for slot,obj,old in members:
                slots.append((str(key),slot,hashlib.sha256(sp.srepr(old).encode()).hexdigest(),
                              hashlib.sha256(sp.srepr(obj).encode()).hexdigest(),sp.Integer(obj!=old)))
        output(label+'OtherSlotIdentityOperands',slots,tuple(b.DIM_ZERO))
        preserved['/'.join(case)]=len(slots)
        output(label+'OtherSlotPreservationResidual',tuple(item[-1] for item in slots),tuple(b.DIM_ZERO))
    output('SourcePinStabilityResidual',{name:sp.Integer(sha(Path(name))!=digest) for name,digest in pins.items()},tuple(b.DIM_ZERO))
    elapsed=time.monotonic()-STARTED
    output('WallSeconds',sp.Float(elapsed),tuple(b.DIM_T))
    output('PeakRssKiB',sp.Integer(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),tuple(b.DIM_ZERO))
    target=base/'objects.pickle';target.write_bytes(pickle.dumps(records,protocol=5))
    summary={'stage':'b','newExportSha256':sha(paths[0]),'sourcePins':pins,
        'sourcePinsAfter':{p:sha(Path(p)) for p in pins},'sourceSnapshots':snapshots,'checks':residuals,'objects':len(records),
        'objectsSha256':sha(target),'cases':[list(map(str,c)) for c in after],
        'changedValueSerializations':changed,'addedKeys':sorted(new_values.keys()-old_values.keys()),
        'removedKeys':sorted(old_values.keys()-new_values.keys()),'otherSlotsPerCase':preserved,
        'nonzeroResidualScalars':sum(v['nonzero'] for v in residuals.values()),
        'wallSeconds':elapsed,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (base/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    if summary['nonzeroResidualScalars']:raise ValueError('b export residual; inspect emitted operands')


if __name__=='__main__':run()
