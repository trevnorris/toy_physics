#!/usr/bin/env python3
"""NP1–NP4 compact, read-only identification of existing synthetic pole operands.

No imports/execution of S11c scripts; no production or diagnostic regeneration.
This checks source ASTs and selected recorded operators/maps, not every output.
"""
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sympy as sp

BASE = Path(__file__).resolve().parents[1]
OUT = BASE/'_measurements/S11_lean_pole_source_checks.json'
PROBE = BASE/'_measurements/S11c_d_nonlinear_pole_contract_probe.py'
CHECK = BASE/'_measurements/S11c_d_nonlinear_pole_contract_check.py'
z = sp.Symbol('z')

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(v): return sp.Matrix([[sp.sympify(x,locals={'z':z}) for x in row] for row in v])
def expr(v): return sp.sympify(v,locals={'z':z})
def same(a,b):
    residual=a-b
    return all(sp.cancel(x)==0 for x in residual) if isinstance(residual,sp.MatrixBase) else sp.cancel(residual)==0

def assignment(tree, function, target):
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==function)
    return next(n.value for n in fn.body if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])==target)

def main():
    pre=json.loads((BASE/'_measurements/S11_lean_pole_preserved_inputs.json').read_text())
    inputs={rel:sha(BASE/rel) for rel in pre['read_only_inputs']}
    assert inputs==pre['read_only_inputs'],'Shared input drift: reconcile, do not overwrite.'
    probe=json.loads(PROBE.with_suffix('.json').read_text())
    record=json.loads((BASE/'_measurements/S11c_d_nonlinear_pole_contract_checks.json').read_text())
    assert record['schema']=='nonlinearPoleV2ExactControls'
    pt=ast.parse(PROBE.read_text());ct=ast.parse(CHECK.read_text())
    cases=assignment(pt,'main','cases')
    native_ast={n.elts[0].value:n.elts[1] for n in cases.elts}
    expected={'scalarDouble':'sp.Matrix([[z**2]])','scalarBothRoots':'sp.Matrix([[z**2 - 1]])',
              'affineJordan':'sp.Matrix([[z, -1], [0, z]])'}
    checks=[]
    def check(name,passed):
        checks.append({'name':name,'passed':bool(passed)})
        assert passed,name
    for name,code in expected.items():
        check(name+'_source_AST',ast.dump(native_ast[name])==ast.dump(ast.parse(code,mode='eval').body))
    check('realization_scalar_source_AST',ast.dump(assignment(ct,'realization','polynomial'))==ast.dump(ast.parse('z**2',mode='eval').body))
    check('forcing_observation_source_AST',ast.dump(assignment(ct,'realization','(source, observe)'))==ast.dump(ast.parse('(1+3*z, 2+5*z)',mode='eval').body))
    cases={c['name']:c for c in probe['cases']}
    N=sp.Matrix([[0,1],[0,0]])
    operators={'scalarDouble':sp.Matrix([[z**2]]),'scalarBothRoots':sp.Matrix([[z**2-1]]),'affineJordan':z*sp.eye(2)-N}
    for name,operator in operators.items():check(name+'_recorded_operator',same(mat(cases[name]['pencil']),operator))
    realization=record['realization']
    check('same_state_operator',same(mat(realization['stateOperator']),N))
    check('same_forcing_coordinate',same(mat(realization['stateInjection']),sp.Matrix([0,1])))
    check('same_observed_coordinate',same(mat(realization['fieldRecovery']),sp.Matrix([[1,0]])))
    check('same_affine_forcing',same(expr(realization['forcing']),1+3*z))
    check('same_affine_observation',same(expr(realization['observation']),2+5*z))
    # Read existing evidence; this is not a claim that its computation was rerun.
    wanted=['RealizationTransfer','AnalyticMapsFullLaurentReconstruction',
            'OmittedPairingBasisChangeRejected','HigherOrderCouplingSurvivesZeroResidue',
            'FrozenForcingObservationAtDoublePoleRejected','doubleUnrestrictedProjectorRejected']
    evidence={c['name']:c for c in record['checks'] if c['name'] in wanted}
    assert set(evidence)==set(wanted)
    for name,c in evidence.items():
        check('historical_evidence:'+name,c['satisfied'])
        if 'mutationResiduals' in c:assert any(expr(v)!=0 for v in c['mutationResiduals'])
        if 'residuals' in c:assert all(expr(v)==0 for v in c['residuals'])
    controls={
        'same_count_wrong_Jordan_sign_rejected':not same(mat(cases['affineJordan']['pencil']),z*sp.eye(2)+N),
        'same_kernel_scaled_scalar_rejected':not same(mat(cases['scalarDouble']['pencil']),sp.Matrix([[2*z**2]])),
        'wrong_forcing_coordinate_rejected':not same(mat(realization['stateInjection']),sp.Matrix([1,0])),
        'frozen_observation_rejected':not same(expr(realization['observation']),sp.Integer(2)),
    }
    assert all(controls.values())
    assert inputs=={rel:sha(BASE/rel) for rel in inputs}
    report={'status':'PASS','scope':'NP compact read-only synthetic object/source identification; no physical poles or production runs.',
        'instrument_sha256':sha(Path(__file__)),'source_sha256':inputs,'sympy_version':sp.__version__,
        'checks':checks,'translation_controls':controls,'historical_evidence':evidence,
        'selected_native_operator_names':list(operators),'completed_utc':datetime.now(timezone.utc).isoformat(),
        'interpretation':'The original diagnostic reports remain historical. Selected operands are identified to the Lean definitions by this instrument and source review; translation software is not kernel-certified.'}
    OUT.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','checks':len(checks),'controls':len(controls)}))
if __name__=='__main__':main()
