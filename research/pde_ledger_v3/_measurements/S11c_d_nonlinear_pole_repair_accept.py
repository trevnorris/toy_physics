#!/usr/bin/env python3
"""Accept the additive pole contract repair without loading any physics operands."""
import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
M = ROOT/'_measurements'
RUN = REPO/'_scratch/s11c/s11c-wide-three-20260916/production/complete'
SCRATCH = REPO/'_scratch/s11c/s11c-nonlinear-pole-repair-20260916'
ADDENDUM = 'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md'
BASELINE = 'd9d4787e'
SUFFIX = 'f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    checks_path = M/'S11c_d_nonlinear_pole_contract_checks.json'
    exact = json.loads(checks_path.read_text())
    residuals = [v for c in exact['checks'] for v in c.get('residuals', [])]
    nonzero = [c for c in exact['checks'] if 'mutationResiduals' in c]
    controls = [c for c in nonzero if 'Rejected' in c['name'] or 'FailsFullRank' in c['name']]
    failures = []

    def require(condition, name):
        if not condition:
            failures.append(name)

    require(all(c['satisfied'] for c in exact['checks']), 'exact control guards')
    require(all(v == '0' for v in residuals), 'literal zero residuals')
    require(all(any(Fraction(v) != 0 for v in c['mutationResiduals']) for c in nonzero), 'nonzero controls/witnesses')
    require(checks_path.read_bytes() == (SCRATCH/'check.stdout').read_bytes(), 'exact checks/stdout identity')
    require((SCRATCH/'check.stderr').stat().st_size == 0, 'empty exact stderr')
    for name, expected in exact['sourceFiles'].items():
        require(digest(ROOT/name) == expected, 'exact source:'+name)
    pins = json.loads((RUN/'preflight.json').read_text())['sourceFiles']
    pin_records = {}
    for name, expected in pins.items():
        current, frozen = digest(ROOT/name), digest(RUN/'source'/name)
        pin_records[name] = {'expected': expected, 'current': current, 'frozen': frozen}
        require(current == frozen == expected, 'active production source:'+name)
    builder = (M/'S11c_d_sympy_builder_report.md').read_bytes()
    suffix = builder[builder.index(b'## Retained user-approved solver/export contract'):]
    suffix_digest = hashlib.sha256(suffix).hexdigest()
    require(suffix_digest == SUFFIX, 'retained builder suffix')
    documents = [
        'directives/S11c_d_sympy_build_PROGRAM_BRIEF.md',
        'directives/S11c_d_sympy_build_directive.md',
        'directives/_measurements/S11c_d_sympy_build_completion_addendum.md',
    ]
    document_records = {}
    for name in documents:
        path = ROOT/name
        old = subprocess.check_output(['git', 'show', BASELINE+':'+str(path.relative_to(REPO))], cwd=REPO)
        text = path.read_text()
        require(name not in pins, 'amended entry point is not pinned:'+name)
        require('S11c_d_NONLINEAR_POLE_CONTRACT.md' in text and 'nonlinearPoleV2' in text, 'corrected authority pointer:'+name)
        require('normalized Riesz residues/projectors' not in text and 'bound pole set + Riesz +' not in text,
                'unrestricted old shorthand removed:'+name)
        document_records[name] = {'baselineCommit': BASELINE, 'baselineSha256': hashlib.sha256(old).hexdigest(),
                                  'correctedSha256': digest(path)}
    tree = ast.parse((ROOT/'scripts/S11c_d_mixing_scattering_sympy_audit.py').read_text())
    gate_records = []
    for class_name, method_name, condition in [('EndModeFrequencyData', 'construct', 'rank == n'),
                                                ('EndResolventAudit', 'pole', 'rank == nullity')]:
        cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == class_name)
        method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == method_name)
        expected = ast.dump(ast.parse(condition, mode='eval').body)
        gates = [n for n in ast.walk(method) if isinstance(n, ast.If) and ast.dump(n.test) == expected]
        assignments = [n for gate in gates for n in ast.walk(gate)
                       if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'projector' for t in n.targets)]
        require(len(gates) == 1 and len(assignments) == 1, 'existing guarded projector:'+class_name)
        gate_records.append({'class': class_name, 'method': method_name,
                             'methodAstSha256': hashlib.sha256(ast.dump(method).encode()).hexdigest(),
                             'pairingGate': condition, 'projectorAssignment': [ast.unparse(n) for n in assignments],
                             'scope': 'Source inspection of unchanged gate; not a new numerical revalidation.'})
    require('POLES_RIESZ_OVERLAP' in {n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)},
            'physical pole construction remains outstanding')
    record = {
        'status': 'ACCEPTED_ADDITIVE_CONTRACT_REPAIR' if not failures else 'UNRESOLVED',
        'approval': 'User approved the proposed nonlinear-pole repair after d9d4787e.',
        'authority': ADDENDUM, 'authoritySha256': digest(ROOT/ADDENDUM), 'schema': 'nonlinearPoleV2',
        'baselineAuthority': 'directives/S11c_d_SHARED_PHYSICS.md',
        'baselineAuthoritySha256': digest(ROOT/'directives/S11c_d_SHARED_PHYSICS.md'),
        'exactValidation': {'path': str(checks_path.relative_to(ROOT)), 'sha256': digest(checks_path),
                            'checks': len(exact['checks']), 'zeroResidualScalars': len(residuals),
                            'rejectionControls': len(controls), 'nonzeroOrderCouplingWitnesses': len(nonzero)-len(controls),
                            'stdoutSha256': digest(SCRATCH/'check.stdout'), 'stderrBytes': (SCRATCH/'check.stderr').stat().st_size},
        'entryPointJoins': document_records, 'existingEndModeGateInspection': gate_records,
        'activeProductionRun': str(RUN), 'activePreflightSha256': digest(RUN/'preflight.json'),
        'pinnedSourceCount': len(pins), 'pinnedSources': pin_records,
        'builderRetainedSuffixSha256': suffix_digest,
        'acceptanceInstrumentSha256': digest(Path(__file__).resolve()),
        'historicalDiagnosticSha256': digest(M/'S11c_d_nonlinear_pole_contract_probe.json'),
        'physicalResultDisposition': 'No physical pole solve, numerical rerun, source rebase or output relabeling performed.',
        'nextStageProvenance': 'Retain accepted baseline hashes and explicitly add the corrected authority and this dependency disposition.',
        'failures': failures,
    }
    payload = json.dumps(record, indent=2)+'\n'
    (M/'S11c_d_nonlinear_pole_repair_checkpoint.json').write_text(payload)
    print(payload, end='')
    if failures:
        raise RuntimeError('Pole repair acceptance failed; see emitted checkpoint')


if __name__ == '__main__':
    main()
