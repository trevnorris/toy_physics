#!/usr/bin/env python3
"""Build the current Lean sources and record the CAS bridge's verified scope.

Run the mutation instrument first. Its canonical hashes must still match;
this command refreshes the full build and checks every selected axiom audit.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys

BASE = Path(__file__).resolve().parents[1]
LEAN = BASE/'lean'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    subprocess.run([sys.executable,str(BASE/'scripts/S10_lean_cas_bridge.py'),'--check'],check=True)
    checks_path = BASE/'_measurements/S10_lean_cas_bridge_checks.json'
    checks = json.loads(checks_path.read_text())
    assert checks['status'] == 'PASS'
    assert checks['canonical_sha256_before'] == checks['canonical_sha256_after']
    for rel,sha in checks['canonical_sha256_after'].items():
        assert digest(BASE/rel) == sha, 'rerun mutation checks after changing '+rel
    assert digest(BASE/'_measurements/S10_lean_cas_bridge_check.py') == checks['instrument_sha256']
    sources = sorted(p for step in ('s9','s10') for p in (LEAN/step).rglob('*.lean') if '_scratch' not in p.parts)
    before = {p:digest(p) for p in sources}
    log_path = LEAN/'verification.log'
    with log_path.open('w') as output:
        subprocess.run(['lake','build'],cwd=LEAN,stdout=output,stderr=subprocess.STDOUT,check=True,
                       env={**os.environ,'LAKE_CACHE_DIR':'.lake/cache'})
    assert before == {p:digest(p) for p in sources}, 'sources changed during build'
    log = log_path.read_text()
    audit_pattern = r"'([^']+)' (?:depends on axioms: \[([^\]]*)\]|does not depend on any axioms)"
    observed = re.findall(audit_pattern,log)
    expected = set()
    for path in sources:
        source = path.read_text()
        assert not re.search(r'\b(sorry|admit|native_decide)\b|^\s*axiom\s|linter\.[\w.]+\s+false',source,re.M),path
        expected.update(re.findall(r'^#print axioms\s+([\w.]+)',source,re.M))
    assert len(observed) == len(dict(observed)) and set(dict(observed)) == expected
    assert all({a.strip() for a in ax.split(',') if a.strip()} <=
               {'propext','Classical.choice','Quot.sound'} for _,ax in observed)
    assert not re.search(r'\b(error|warning):',log)
    completion = re.search(r'Build completed successfully \(\d+ jobs\)\.',log)
    assert completion
    manifest_path = BASE/'_measurements/S10_lean_cas_bridge_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    arithmetic = [r for key in ('engines','minor_engines','rerun_engines','count_engines','generic_count_engines','root_engines','coincidence_engines','record_engines') for d in manifest[key].values() for r in d['records']]
    metadata_records = manifest['record_metadata']
    record_arithmetic = [r for d in manifest['record_engines'].values() for r in d['records']]
    assert len(metadata_records) == 48 and len(record_arithmetic) == 9
    assert checks['independently_parsed_metadata_records'] == len(metadata_records)
    assert all(r['lean_semantic_claims'] for r in metadata_records)
    assert all('S10Audit.CAS.'+name in expected for r in metadata_records for name in r['lean_semantic_claims'])
    assert all('S10Audit.CAS.'+name in expected for r in record_arithmetic for name in r['record_semantics']['lean_semantic_claims'])
    signs = [r for r in metadata_records if r['kind']=='root_sign']
    assert len(signs) == 16 and sum(r['reported_sign']=='undecided' for r in signs) == 2
    assert all(any(n.endswith('_computed') for n in r['lean_semantic_claims']) for r in signs)
    root_records = [r for d in manifest['root_engines'].values() for r in d['records']]
    root_filters = manifest['root_filter_records']
    coincidence_records = manifest['coincidence_records']
    coincidence_arithmetic = [r for d in manifest['coincidence_engines'].values() for r in d['records']]
    assert len(coincidence_records) == 56 and len(coincidence_arithmetic) == 13
    assert len(coincidence_records) == checks['independently_parsed_coincidence_records']
    assert all(r['lean_semantic_claims'] for r in coincidence_records)
    assert all('S10Audit.CAS.'+name in expected for r in coincidence_records for name in r['lean_semantic_claims'])
    assert all('S10Audit.CAS.'+name in expected for r in coincidence_arithmetic
               for name in r['coincidence_arithmetic'].get('lean_semantic_claims',[]))
    assert len(root_records) == 27 and len(root_filters) == 3
    assert all('S10Audit.CAS.'+name in expected for r in root_records
               for name in r['root_semantics']['lean_semantic_claims'])
    assert all('S10Audit.CAS.'+name in expected for r in root_filters for name in r['lean_semantic_claims'])
    assert all('lean_raw_tree' in c for r in root_records if r['root_semantics']['kind'] == 'candidates'
               for c in r['cells'])
    counts = [r for key in ('count_engines','generic_count_engines') for d in manifest[key].values() for r in d['records']]
    generic_counts = [r for d in manifest['generic_count_engines'].values() for r in d['records']]
    assert len(counts) == 112 and len(generic_counts) == 42
    assert all(r['count_semantics']['lean_chart'] == 'GenericChart sigma k' and
               r['count_semantics']['coordinate_domain'] == ['k0 != 0','k1 != 0','k2 != 0'] and
               r['count_semantics']['coefficient_domain'] == ['rho != 0','mu != 0','0 < sigma','sigma != 1']
               for r in generic_counts)
    assert all('S10Audit.CAS.'+r['count_semantics']['lean_semantic_claim'] in expected for r in counts)
    assert all(len(r['cells']) == 1 and r['cells'][0]['dimension'] == [0,0,0] for r in counts)
    cells = sum(len(r['cells']) for r in arithmetic)
    predicates = sum(len(d['loci']) for d in manifest['locus_engines'].values())
    points = sum(len(d['points']) for d in manifest['locus_engines'].values())
    assert cells == checks['independent_parser_cells']
    assert predicates*27 == checks['coordinate_sign_pattern_comparisons']
    assert points == checks['independently_parsed_points']
    assert len(root_filters) == checks['independently_parsed_root_filters']
    rejected = sum(c['outcome']=='REJECTED' for c in checks['lean_checks'])
    controls = sum(c['outcome']=='PASS' for c in checks['lean_checks'])
    assert rejected == 36 and controls == 14
    pins = [LEAN/name for name in ('lakefile.toml','lake-manifest.json','lean-toolchain')]
    instruments = [BASE/'scripts'/name for name in ('S10_lean_cas_bridge.py','S10_lean_cas_minors.py',
                   'S10_lean_cas_loci.py','S10_lean_cas_reruns.py','S10_lean_cas_counts.py',
                   'S10_lean_cas_generic_counts.py','S10_lean_cas_roots.py','S10_lean_cas_coincidence.py','S10_lean_cas_records.py',
                   'S10_cross_engine_comparator.py','S10_exports.py')]
    inputs = [BASE/'scripts/out/S10_anisotropic_strata_sympy_audit.out',
              BASE/'mathematica/out/S10_anisotropic_strata_mathematica_audit.out']
    additional = [manifest_path,BASE/'_measurements/S10_lean_cas_bridge_check.py',checks_path,Path(__file__)]
    assert digest(BASE/'scripts/S10_exports.py') == 'bc8de16bae05dcf6caa71d82184f5aa95e2a9d6fd157fdbb674b88f185ed34c9'
    lines = ['S10 CAS expression, minor, locus, rerun and generic/exceptional count verification',
             'Command: LAKE_CACHE_DIR=.lake/cache lake build',
             'Working directory: research/pde_ledger_v3/lean','Exit status: 0',
             'Library targets: S9Pilot, S10Pilot, S10Controls, S10Anisotropic, S10Audit',
             f'Canonical modules: warningAsError = true; {len(sources)} Lean source files',
             f'Total selected axiom audits: {len(expected)}',
             f'CAS bridge selected audits: {sum(n.startswith("S10Audit.CAS.") for n in expected)}',
             'Observed axioms: propext, Classical.choice, Quot.sound',
             f'Axiom-free selected declarations: {sum(not ax for _,ax in observed)}',
             'No proof admissions, custom axioms, native_decide, or disabled linters.',
             completion[0],'Full build log SHA-256: '+digest(log_path),'',
             f'Arithmetic: {len(arithmetic)} tagged records, {cells} scalar expressions, 2 engines.',
             f'Generic and exceptional count records: {len(counts)}; every semantic binding is included in the axiom audit.',
             f'Root-list/count records: {len(root_records)}; raw filter-list records: {len(root_filters)}; all semantic claims audited.',
             f'Coincidence records: {len(coincidence_arithmetic)} arithmetic, {len(coincidence_records)} logical/container records; semantic bindings audited.',
             f'Aggregate/metadata records: {len(record_arithmetic)} arithmetic, {len(metadata_records)} full metadata payloads; all semantic bindings audited.',
             f'Coordinate predicates: {predicates} locus records; {points} targeted point records.',
             'Provenance command: python3 ../scripts/S10_lean_cas_bridge.py --check',
             'Regression command: python3 ../_measurements/S10_lean_cas_bridge_check.py',
             f'Independent parser comparison: PASS, {cells} scalar expressions.',
             f'Independent coincidence parser comparison: PASS, {len(coincidence_records)} complete logical/container payloads.',
             f'Independent metadata parser comparison: PASS, {len(metadata_records)} complete payloads, including {len(signs)} reported root signs.',
             f'Coordinate parser comparison: PASS, {predicates*27} locus sign patterns and {points} points.',
             f'Parser rejection checks: PASS, {len(checks["strict_parser_rejections"])} arithmetic/transcript cases, '+
             f'{len(checks["coordinate_parser_rejections"])} coordinate cases.',
             f'Lean mutation checks: PASS, {rejected} rejected mutations and {controls} passing positive controls.',
             'Canonical source, input and frozen export hashes unchanged by mutation checks.',
             'Parsers/transcript translation remain tested software; generated proofs are kernel checked.','',
             'CAS bridge axiom audits from the full build:']
    lines += [line for line in log.splitlines() if "'S10Audit.CAS." in line and re.search(audit_pattern,line)]
    lines += ['','All verified source hashes; paths relative to research/pde_ledger_v3/:']
    lines += [digest(p)+'  '+str(p.relative_to(BASE)) for p in sources+pins+instruments+inputs+additional]
    (LEAN/'s10/CAS_BRIDGE_VERIFICATION.txt').write_text('\n'.join(lines)+'\n')
    print(f'PASS: {len(expected)} standard-axiom audits, {len(sources)} sources; verification record written.')


if __name__ == '__main__':
    main()
