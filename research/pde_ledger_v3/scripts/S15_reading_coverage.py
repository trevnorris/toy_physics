"""Publish item-3 exponent coverage and item-6 source-backed slots from a saved reading build.

Run under scripts/s11c_guarded_run.py. This reads native operands and checkpoints;
it does not evaluate any physical field or rerun the producing engines.
"""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import re
import tarfile

import sympy as s

import S15_reading_math as math


def member(archive, name):
    with archive.extractfile(name) as stream:
        data = stream.read()
    return json.loads(gzip.decompress(data) if name.endswith('.gz') else data)


@lru_cache(maxsize=None)
def feasible_selections(planes):
    """A feasible selection is the set of all equalities holding on an affine flat."""
    def flat(rows):
        reduced, pivots = s.Matrix(rows).rref()
        if 4 in pivots:  # coefficients of a1,a2,b1,b2, then the affine constant
            return None
        return tuple(tuple(reduced.row(i)) for i in range(len(pivots)))

    flats = {()}
    for plane in planes:
        intersections = {flat(rows + (plane,)) for rows in flats}
        flats.update(intersections - {None})
    # Over the reals finitely many proper subflats cannot cover a flat.
    # All planes outside this closure can therefore be false simultaneously.
    return {tuple(i for i,p in enumerate(planes) if flat(rows + (p,)) == rows): len(rows)
            for rows in flats}


def exponent_coverage(archive, config):
    entries = []
    members = set(archive.getnames())
    for index, name in enumerate(config['required_objects']):
        suffix = f'{index:03d}.json.gz'
        obj = member(archive, 'objects/' + suffix)
        if (obj['object'], obj['index'], obj['build_id']) != (name, index, config['build_id']):
            raise ValueError('object provenance mismatch: ' + name)
        ref_name = 'refinement/' + suffix
        refinement = member(archive, ref_name) if ref_name in members else {}
        if refinement and (refinement['object'], refinement['build_id']) != (name, config['build_id']):
            raise ValueError('refinement provenance mismatch: ' + name)
        branch = refinement.get('branch')
        native_planes = branch['single_planes'] if branch else obj['exponent_loci']
        if not native_planes:
            continue
        planes = tuple(tuple(map(Fraction, p)) for p in native_planes)
        if any(len(p)!=5 or not any(p[:4]) for p in planes):
            raise ValueError('invalid affine exponent condition')
        row, = (r for r in obj['rows'] if r['role']=='cross' and r['field']=='item3/object')
        original = obj['rows'][row['same_comparison_as']] if 'same_comparison_as' in row else row
        if original['status']=='exact_identity':
            continue
        samples = list(original['samples'])
        for field in refinement.get('fields', []):
            if field['field']=='item3/object':
                samples.extend(field['additional_samples'])
        feasible = feasible_selections(planes)
        points = {tuple(map(Fraction, q['exponents'])) for q in samples}
        if any(len(point)!=4 for point in points):
            raise ValueError('invalid exponent tuple')
        hit = {tuple(i for i,p in enumerate(planes) if sum(x*y for x,y in zip(p[:4], point)) + p[4]==0)
               for point in points}
        if not hit.issubset(feasible):
            raise ValueError('archived point outside the computed selection roster: '+name)
        absent = sorted(feasible.keys()-hit)
        entries.append(dict(object=name, archive_index=index, operand_ids=row['left']+row['right'],
                            planes=native_planes, sampled=sorted(hit), untested=absent,
                            untested_by_codimension=dict(sorted(Counter(feasible[x] for x in absent).items()))))
    return entries


def source_arguments(source):
    """Decode the accepted Wolfram registration loop and its dimensional basis."""
    pattern = (r'Do\[register\[centreVelocityChange\[j\],([^,]+),\{(\w+)\}\],'
               r'\{j,(\d+)\}\];')
    declarations = list(re.finditer(pattern, source))
    if len(declarations)!=1:
        raise ValueError('unrecognized centreVelocityChange declaration')
    declaration = declarations[0]
    dim = declaration.group(2)
    vector, = re.findall(r'\b'+re.escape(dim)+r'=\{(-?\d+,-?\d+,-?\d+)\};', source)
    return declaration.group(), source[:declaration.start()].count('\n')+1, {
        f'centreVelocityChange[{i}]':list(map(int,vector.split(',')))
        for i in range(1,int(declaration.group(3))+1)}


def argument_slots(archive, source, source_info):
    declaration, line, signatures = source_arguments(source)
    rows = member(archive,'context_checks.json.gz')['argument_slots']
    family, = (q for q in member(archive,'context.json.gz') if q['field']=='family/centre_velocity_changes')
    if family['status']!='same_values':
        raise ValueError('centre velocity binding unresolved')
    rid, = family['right']
    record, = (r for r in member(archive,'coverage.json.gz') if r['id']==rid)
    with gzip.open(archive.extractfile('operands.jsonl.gz'),'rt') as f:
        raw, = (v['native_operand'] for v in map(json.loads,f) if v['sha256']==record['payload'])
    mapping = {}
    for rule in math.original.wl_parse(raw).args:
        lhs, rhs = rule.args
        head = math.original.key_text(rhs.args[0])
        if head not in signatures or len(rhs.args)!=2 or rhs.args[1].value!='t':
            raise ValueError('unknown live velocity response binding')
        py = math.translate(rhs,'wl',())
        j = re.fullmatch(r'centreVelocityChange\[(\d+)\]',head)[1]
        if math.translate(lhs,'wl',()) != math.parse(f'Derivative(X_1_{j}(t),t)','py'):
            raise ValueError('source centre or velocity slot differs from declared binding')
        mapping[(py.func.__name__,0,'t')] = head
    if set(mapping.values())!=set(signatures):
        raise ValueError('velocity rule roster differs from registrations')
    found = set()
    for row in rows:
        key=tuple(row['argument'])
        if key not in mapping:
            continue
        if key in found:
            raise ValueError('duplicate velocity argument')
        found.add(key)
        head=mapping[key]
        dim=signatures[head]
        row.update(wl_head=head, wl_slot=1, wl_expected=dim, wl_velocity_rules=rid,
                   wl_source=dict(**source_info,line=line,declaration=declaration),
                   status='matches_source_declared_slot' if row['expected']==row['observed']==[dim]
                   else 'slot_dimension_difference')
    if found!=set(mapping):
        raise ValueError('missing SymPy observed velocity argument')
    return rows


def publish(manifest, archive_path, repo, output):
    archive_info=manifest['detail_archive']
    with archive_path.open('rb') as f:
        digest=hashlib.file_digest(f,'sha256').hexdigest()
    if digest!=archive_info['sha256'] or archive_path.stat().st_size!=archive_info['bytes']:
        raise ValueError('detail archive differs from its manifest')
    wl='research/pde_ledger_v3/mathematica/S15_separation_departure.wl'
    source=(repo/wl).read_bytes()
    source_sha=hashlib.sha256(source).hexdigest()
    if source_sha!=manifest['engine_sources'][wl]:
        raise ValueError('Wolfram source differs from the archived engine')
    with tarfile.open(archive_path,'r:') as archive:
        config=member(archive,'configuration.json')
        if config['build_id']!=manifest['implementation'] or config['inputs']!=manifest['repository_inputs']:
            raise ValueError('configuration differs from manifest')
        entries=exponent_coverage(archive,config)
        slots=argument_slots(archive,source.decode(), dict(path=wl,revision=manifest['engine_commits']['reading'],sha256=source_sha))
    counts=dict(objects_with_planes=len(entries),objects_with_gaps=sum(bool(e['untested']) for e in entries),
                sampled=sum(len(e['sampled']) for e in entries),untested=sum(len(e['untested']) for e in entries))
    data=dict(build_id=config['build_id'],archive_sha256=digest,field='item3/object',
              plane_convention='Each [A,B,C,D,E] represents A*a1+B*a2+C*b1+D*b2+E = 0 over real exponents.',
              selection_convention='Each list gives exactly the zero-based plane indices that hold; all other planes are false. Counts are per object.',
              untested_reason='No exponent point in the archived item3/object grid realizes this selection. The grid targeted generic points, individual planes and the explicitly printed Wolfram conjunctions.',
              scope='Selections count the printed affine conditions on sampled emitted family results. They do not classify the conditions for inverse-square scaling.',
              counts=counts,objects=entries,argument_counts=dict(Counter(r['status'] for r in slots)),
              argument_slots=slots)
    output.parent.mkdir(parents=True,exist_ok=True)
    tmp=output.with_suffix('.tmp')
    with tmp.open('w') as f:
        f.write('{\n')
        for k,v in data.items():
            f.write('  '+json.dumps(k)+': ')
            if k in ('objects','argument_slots'):
                f.write('[\n    '+',\n    '.join(json.dumps(x,ensure_ascii=False) for x in v)+'\n  ]')
            else:
                f.write(json.dumps(v,ensure_ascii=False))
            f.write('\n' if k=='argument_slots' else ',\n')
        f.write('}\n')
    tmp.replace(output)
    return counts,data['argument_counts']


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('manifest',type=Path)
    parser.add_argument('archive',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    repository=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
    counts, slots=publish(json.loads(args.manifest.read_text()),args.archive,repository,args.output)
    print(json.dumps(dict(exponent_coverage=counts,argument_counts=slots),sort_keys=True))
