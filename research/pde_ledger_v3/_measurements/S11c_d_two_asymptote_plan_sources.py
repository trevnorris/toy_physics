#!/usr/bin/env python3
"""Pin selected source/JSON routes for two-end method design; no science imports."""
import ast
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path('/var/projects/toy_physics')
PROJECT = ROOT / 'research/pde_ledger_v3'
M = PROJECT / '_measurements'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def route(path):
    return {'path': str(path), 'canonicalPath': str(path.resolve(strict=True)),
            'bytes': path.stat().st_size, 'sha256': sha(path)}

def read(path):
    return json.loads(path.read_text())

def source(path, definitions):
    text = path.read_text()
    tree = ast.parse(text)
    records = []
    for qualified in definitions:
        body = tree.body
        for name in qualified.split('.'):
            node = next(v for v in body if getattr(v, 'name', None) == name)
            body = node.body
        excerpt = ''.join(text.splitlines(keepends=True)[node.lineno-1:node.end_lineno])
        records.append({'definition': qualified, 'firstLine': node.lineno,
                        'lastLine': node.end_lineno,
                        'sourceSha256': hashlib.sha256(excerpt.encode()).hexdigest()})
    return {**route(path), 'definitions': records}

def packets(checkpoint_name, names, role):
    checkpoint = M / checkpoint_name
    data = read(checkpoint)
    base = Path(data['runDirectory'])
    artifacts = []
    for name in names:
        actual = route(base / name)
        expected = data['artifacts'][name]
        assert actual['sha256'] == expected['sha256'] and actual['bytes'] == expected['bytes']
        artifacts.append({'manifestKey': name, **actual, 'restored': False})
    return {'checkpoint': route(checkpoint), 'scopeAsStored': data.get('scope'),
            'proposedRole': role, 'artifacts': artifacts}

def main():
    sources = [
        source(PROJECT/'scripts/S11c_d_mixing_scattering_sympy_audit.py', [
            'ReducedPencil.__init__', 'ReducedActionAssembly.construct',
            'BoundedSourceFourierAssembly.bounded', 'BoundedSourceFourierAssembly.construct',
            'EdgeReduction.subtraction', 'EdgeReduction.prescribe', 'EdgeReduction.weak_limit',
            'ConstantEndPencil.background', 'ConstantEndPencil.profile_limit_operands',
            'RectangularModeJets.__init__', 'RectangularModeJets.pair']),
        source(M/'S11c_d_continuum_grades.py', ['split', 'term_joins', 'main']),
        source(M/'S11c_d_continuum_boundary.py', ['phase', 'construct_end']),
        source(M/'S11c_d_continuum_response.py', ['systems', 'solve', 'channels']),
        source(M/'S11c_d_uniform_source.py', ['construct']),
    ]
    selected = [
        packets('S11c_d_reduced_action_source_checkpoint.json',
                ['reduced-action.pickle','actions.pickle'], 'Native infinite-domain source rows and probe actions; restore only under containment.'),
        packets('S11c_d_reduced_action_assembly_checkpoint.json', ['assembly.pickle'],
                'Reuse original local/nonlocal action assembly, without reconstructing it.'),
        packets('S11c_d_continuum_grade_checkpoint.json', ['continuum-grades.pickle'],
                'Reuse symbolic coefficient grades. Finite-domain term-factor integrals are not whole-line action operands.'),
        packets('S11c_d_continuum_boundary_checkpoint.json',
                ['left-pencil.pickle','right-pencil.pickle','left-mode-16.pickle','left-mode-17.pickle',
                 'right-mode-16.pickle','right-mode-17.pickle','continuum-boundary.pickle'],
                'Reuse symbolic end pencil/Taylor tables; numerical invariant pairs are comparison/gauge evidence, not exact tail cancellation.'),
        packets('S11c_d_source_fourier_factorization_checkpoint.json', ['source-factorization.pickle'],
                'Boundary evidence only: stored reconstruction is finite-domain and does not license an infinite-limit interchange.'),
    ]
    ref = read(M/'S11c_d_reference_kernel_inputs.json')
    own_routes = []
    for name in ('uniformPacket','reductionPacket','physicalInput'):
        p = Path(ref[name]); actual = route(p)
        expected = next(v for v in ref['files'] if v['path'] == str(p))
        assert all(actual[k] == expected[k] for k in ('bytes','sha256'))
        own_routes.append({'role': name, **actual, 'restored': False})
    candidate_checkpoint = M/'S11c_d_outgoing_prescription_checkpoint.json'
    candidate = read(candidate_checkpoint)
    assert candidate['status'] == 'VALIDATED_FIXED_INPUT_SEPARATED_POINT_OUTGOING_PRESCRIPTION_CANDIDATE'
    accepted = []
    for item in [candidate['productionCheckpoint'], candidate['validationCheckpoint'],
                 *candidate['primaryArtifacts'].values()]:
        actual = route(Path(item['path']))
        assert all(actual[k] == item[k] for k in ('bytes','sha256'))
        accepted.append(actual)
    source_root = Path(candidate['primaryArtifacts']['checks.json']['path']).parent
    tails = {p.name: read(p) for p in sorted(source_root.glob('inverse-*-tail-*.json'))}
    assert len(tails) == 50 and all(v['decays'] for v in tails.values())
    suffix_bytes = (M/'S11c_d_sympy_builder_report.md').read_bytes()
    suffix_bytes = suffix_bytes[suffix_bytes.index(b'## Retained user-approved solver/export contract'):]
    suffix_sha = hashlib.sha256(suffix_bytes).hexdigest()
    assert suffix_sha == 'f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2'
    record = {'status': 'TWO_ASYMPTOTE_METHOD_PLAN_SOURCE_IDENTITY_ONLY',
        'recordedUtc': datetime.now(timezone.utc).isoformat(),
        'authorizationBasis': 'User said Continue after the completed bounded prescription and the stated next step of two-asymptote response planning.',
        'method': 'Standard-library AST, JSON, filesystem metadata and byte hashes only; no scientific pickle restoration or execution.',
        'sourceDefinitions': sources, 'selectedPackets': selected, 'ownReferenceRoutes': own_routes,
        'validatedCandidateCheckpoint': route(candidate_checkpoint), 'verifiedAcceptedArtifacts': accepted,
        'savedTailPowerCensus': {'files':50, 'powers':sorted({int(v['powerOfReciprocalMomentum']) for v in tails.values()}),
                               'inference':'At least one entry has a 1/k leading tail. Decay alone does not certify absolute Fourier integrability at coincident positions.'},
        'context': [route(PROJECT/n) for n in [
            'directives/S11c_d_FORM_constructor_plan.md','directives/S11c_d_FORM_build_directive.md',
            'directives/S11c_d_SCATTERING_FORM_AMENDMENT.md','directives/S11c_d_SHARED_PHYSICS.md',
            '_measurements/S11c_d_FORM_implementation_coverage.md','_measurements/S11c_d_A9_dependency_report.md']],
        'protectedBuilderSuffixSha256':suffix_sha, 'sharedGuard':route(ROOT/'scripts/s11c_guarded_run.py'),
        'scientificOperations':0, 'scientificPacketRestorations':0, 'externalReviewsLaunched':0,
        'implementationExists':False, 'methodClearance':False, 'productionAuthorized':False,
        'referenceLiterature':{'url':'https://dlmf.nist.gov/1.16','purpose':'Definitions of distribution/test-function action and Fourier transforms, including the Heaviside transform. Not authority for this coupled physical method; native signs/measures govern.'},
        'outputs':[route(PROJECT/'directives/S11c_d_two_asymptote_response_plan.md'),
                   route(PROJECT/'directives/S11c_d_two_asymptote_first_stage_contract.md')],
        'generator':route(Path(__file__).resolve())}
    destination=M/'S11c_d_two_asymptote_plan_sources.json'
    with destination.open('x') as f:
        json.dump(record,f,indent=2);f.write('\n')
    print(json.dumps({'status':record['status'],'selectedPacketFiles':sum(len(v['artifacts']) for v in selected),
        'sourceFiles':len(sources),'receipt':route(destination),'scientificPacketRestorations':0},indent=2))

if __name__ == '__main__':
    main()
