#!/usr/bin/env python3
"""Finish saved factor output; retain native factor and certificate workers."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import shutil
import sympy as sp
import S11c_d_remaining_case_factor_recover as recovery

c, f, engine = recovery.c, recovery.f, recovery.engine
ORIGIN = f.STORE/'s11c-remaining-case-factors-20260919/preflight-recovery-02/complete'
PLAN = f.M/'S11c_d_remaining_case_factor_output_plan.md'
PROVENANCE = 'PY_S11CD_'+c.PREFIX+'_NEW_0_PROVENANCE'
NATIVE_OUTPUT = c.output_replay
NATIVE_EMIT = c.emit


def decoded(path):
    values = {}
    for line in recovery.grades.decoded_lines(path):
        tag, separator, body = line.rstrip('\n').partition(': ')
        f.require(separator and tag not in values, 'unique complete decoded tag')
        values[tag] = recovery.grades._restore(body)
    return values


def provenance_units(name, value, physical, **kwargs):
    if name.startswith(c.PREFIX+'_NEW_') and name.endswith('_PROVENANCE'):
        f.require(isinstance(value, dict) and set(value) == {'source', 'unionIndex'}
                  and type(value['unionIndex']) is int and value['unionIndex'] >= 0,
                  'actual factor provenance address')
        f.require(not kwargs, 'unmodified original provenance call')
        # Association paths use their literal key, not a positional integer.
        return physical(name, value, zero_dimensions={('unionIndex',): (0, 0, 0)})
    return physical(name, value, **kwargs)


def fixed_emit(*args):
    original = engine.physical
    def physical(name, value, **kwargs):
        return provenance_units(name, value, original, **kwargs)
    engine.physical = physical
    try:
        return NATIVE_EMIT(*args)
    finally:
        engine.physical = original


def load_with_output_sources(base):
    states, accepted, catalogue, manifest = recovery.load_saved(base)
    old = json.loads((ORIGIN/'inputs.json').read_text())
    f.require(old == manifest, 'exact completed recovery source and input manifest')
    outcome = json.loads((ORIGIN.parent/'factors_recover.invocation.json').read_text())
    f.require(outcome == json.loads((ORIGIN.parent/'active.json').read_text()) and outcome['exitCode'] == 1,
              'actual final metadata failure retained')
    f.require((ORIGIN.parent/'factors_recover.stderr').read_text().rstrip().endswith(
        'ValueError: complete units and independent grades'), 'recorded final unit guard')
    for path in ORIGIN.rglob('*'):
        if path.is_file() and 'source' not in path.relative_to(ORIGIN).parts:
            manifest['inputPackets'][str(path)] = f.digest(path)
    for name in ('active.json', 'factors_recover.invocation.json', 'factors_recover.stdout', 'factors_recover.stderr'):
        path = ORIGIN.parent/name; manifest['inputPackets'][str(path)] = f.digest(path)
    for path in (Path(__file__).resolve(), PLAN):
        name = str(path.relative_to(f.ROOT)); manifest['sourceFiles'][name] = f.digest(path)
        target = base/'source'/name; target.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(path, target)
    f.save(base/'inputs.json', manifest)
    f.save(base/'output-reuse.json', {'originalOutcome': outcome, 'completedRows': [0, 3, 9, 28],
        'nativeFactorComputationsRepeated': 0, 'certificatesRecomputed': 0,
        'sourceManifestBeforeOutputAdapter': old, 'coordinatorJoin': coordinator()[1]})
    return states, accepted, catalogue, manifest


def resume_checks(previous, manifest):
    if previous == ORIGIN:
        catalogue = f.unpickle(previous/'integral-catalogue.pickle')
        selected = json.loads((recovery.ORIGIN/'preflight.json').read_text())['selected']
        f.require(selected == [0, 3, 9, 28], 'all completed preflight rows')
        for index in selected:
            c.verify_row(previous/'rows'/str(index).zfill(3), catalogue['unique'][index])
        # This permits saved-row reuse only. The failed whole run is not accepted.
        return {'mode': 'preflight', 'savedRowsOnly': True}
    checks = json.loads((previous/'checks.json').read_text())
    outcome = json.loads((previous.parent/'factors_finish.invocation.json').read_text())
    f.require(outcome == json.loads((previous.parent/'active.json').read_text()) and outcome['exitCode'] == 0
              and (previous.parent/'factors_finish.stderr').read_bytes() == b''
              and json.loads((previous.parent/'factors_finish.stdout').read_text()) == checks,
              'clean fully completed output preflight')
    for name, record in checks['artifacts'].items():
        f.require(c.artifact(previous/name) == record, 'accepted preflight artifact hash')
    return checks


def resume_manifest_matches(previous, manifest):
    old = json.loads((previous/'inputs.json').read_text())
    if previous != ORIGIN:
        return old == manifest
    expected = copy.deepcopy(manifest)
    # Only new provenance entries may differ; every old entry must be identical.
    for key in ('sourceFiles', 'inputPackets'):
        f.require(all(expected[key].get(name) == value for name, value in old[key].items()),
                  'unchanged prior source and input entry')
        expected[key] = old[key]
    return expected == old


def corrected_provenance_metadata(value):
    changed = 0; rows = []
    for path, descriptor in value:
        entries = []
        for key, operand in descriptor:
            if tuple(str(v) for v in path) == ('unionIndex',) and str(key) == 'DIMENSION_L_T_M':
                f.require(str(operand) == 'ZERO_MAP', 'the actual zero-index missing unit')
                operand = sp.Tuple(0, 0, 0); changed += 1
            entries.append(sp.Tuple(key, operand))
        rows.append(sp.Tuple(path, sp.Tuple(*entries)))
    f.require(changed == 1, 'one precise provenance unit correction')
    return sp.Tuple(*rows)


def output_replay(base, states, unique, catalogue, packets, manifest):
    result = NATIVE_OUTPUT(base, states, unique, catalogue, packets, manifest)
    if set(packets) == {0, 3, 9, 28}:
        old, new = decoded(ORIGIN/'full.out'), decoded(base/'full.out')
        f.require(set(old) == set(new), 'complete original and final tag census')
        metadata_tag = PROVENANCE.replace('PY_S11CD_', 'PY_S11CD_METADATA_', 1)
        expected = corrected_provenance_metadata(old[metadata_tag])
        f.require(new[metadata_tag] == expected and old[PROVENANCE] == new[PROVENANCE],
                  'unchanged actual provenance and precise dimensionless index metadata')
        input_tag = 'PY_S11CD_'+c.PREFIX+'_INPUTS'
        f.require(old[input_tag] == engine.cas(json.loads((ORIGIN/'inputs.json').read_text()))
                  and new[input_tag] == engine.cas(manifest), 'actual old and new input manifests')
        exceptions = {metadata_tag, input_tag, input_tag.replace('PY_S11CD_', 'PY_S11CD_METADATA_', 1)}
        differences = [tag for tag in old if old[tag] != new[tag]]
        f.require(set(differences) <= exceptions, ('only declared metadata/provenance differences', differences))
        # Changed physical units and misplaced path annotations do not satisfy
        # the actual saved metadata comparison.
        controls = []
        for wrong in (sp.Tuple(1, 0, 0), sp.Tuple(0, 1, 0), sp.Tuple(0, 0, 1)):
            mutated = expected.xreplace({sp.Tuple(0, 0, 0): wrong})
            f.require(mutated != new[metadata_tag], 'changed-index-unit control rejected')
            controls.append(str(wrong))
        evidence = {'identicalTags': len(old)-len(differences), 'changedTags': differences,
                    'unitPath': ['unionIndex'], 'beforeUnit': 'ZERO_MAP', 'afterUnit': [0, 0, 0],
                    'rejectedUnitControls': controls, 'allPhysicalPayloadsIdentical': True,
                    'originalTranscriptSha256': f.digest(ORIGIN/'full.out'),
                    'finalTranscriptSha256': f.digest(base/'full.out')}
        f.save(base/'emission-differences.json', evidence); result['outputRepair'] = evidence
    return result


def coordinator():
    original = ast.parse(inspect.getsource(recovery.main)).body[0]
    changed = copy.deepcopy(original); counts = {'load': 0, 'checks': 0, 'manifest': 0, 'emit': 0}
    class Forward(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if ast.unparse(node.func) == 'load_saved':
                node.func.id = 'load_with_output_sources'; counts['load'] += 1
            return node
        def visit_Assign(self, node):
            self.generic_visit(node)
            if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
                if node.targets[0].id == 'previous_checks':
                    node.value = ast.parse('resume_checks(previous, manifest)', mode='eval').body; counts['checks'] += 1
                if node.targets[0].id == 'previous_emit':
                    node.value = ast.Name(id='fixed_emit', ctx=ast.Load()); counts['emit'] += 1
            return node
        def visit_Compare(self, node):
            self.generic_visit(node)
            if ast.unparse(node) == "json.loads((previous / 'inputs.json').read_text()) == manifest":
                counts['manifest'] += 1
                return ast.parse('resume_manifest_matches(previous, manifest)', mode='eval').body
            return node
    Forward().visit(changed)
    class Reverse(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if ast.unparse(node.func) == 'load_with_output_sources': node.func.id = 'load_saved'
            if ast.unparse(node.func) == 'resume_checks':
                return ast.parse("json.loads((previous/'checks.json').read_text())", mode='eval').body
            if ast.unparse(node.func) == 'resume_manifest_matches':
                return ast.parse("json.loads((previous/'inputs.json').read_text()) == manifest", mode='eval').body
            return node
        def visit_Assign(self, node):
            self.generic_visit(node)
            if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id == 'previous_emit':
                node.value = ast.parse('c.emit', mode='eval').body
            return node
    f.require(counts == {'load': 1, 'checks': 1, 'manifest': 1, 'emit': 1}
              and ast.dump(Reverse().visit(copy.deepcopy(changed))) == ast.dump(original),
              'whole coordinator reverse AST join: output and saved-row routing only')
    namespace = dict(vars(recovery), load_with_output_sources=load_with_output_sources,
                     resume_checks=resume_checks, resume_manifest_matches=resume_manifest_matches, fixed_emit=fixed_emit)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[changed], type_ignores=[])), str(Path(__file__)), 'exec'), namespace)
    return namespace['main'], {'wholeCoordinatorReverseAstJoin': True, 'changes': counts,
        'originalAstSha256': hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'derivedAstSha256': hashlib.sha256(ast.dump(changed).encode()).hexdigest()}


if __name__ == '__main__':
    c.output_replay = output_replay
    coordinator()[0]()
