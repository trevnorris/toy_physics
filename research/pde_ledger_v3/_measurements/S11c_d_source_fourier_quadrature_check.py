#!/usr/bin/env python3
"""Evaluate every accepted bounded source factor on the approved test fields."""
import argparse
import ast
import contextlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import time

import numpy as np
import sympy as sp
from sympy.core.function import AppliedUndef

import S11c_d_numerical_action_check as native
from S11c_d_numerical_action_check import ROOT, STORE, engine, digest, save, atomic_pickle, number, progress
from S11c_d_output_codec import decoded_lines, restore_emission_index
from ledger_fold import _restore

M = ROOT/'_measurements'
CHECKPOINT = M/'S11c_d_source_fourier_factorization_checkpoint.json'
NUMERICAL = M/'S11c_d_numerical_action_checkpoint.json'
SOURCE = M/'S11c_d_reduced_action_source_checkpoint.json'
PLAN = M/'S11c_d_source_fourier_quadrature_plan.md'
PREFIX = 'SOURCE_FOURIER_QUADRATURE_LAB_HELD_RHO4_CONSTANT'
SOURCES = tuple(dict.fromkeys((*native.SOURCES, CHECKPOINT, NUMERICAL, SOURCE, PLAN, Path(__file__).resolve())))


def accepted(path):
    record = json.loads(path.read_text()); base = Path(record['runDirectory'])
    for name, item in record['artifacts'].items():
        if digest(base/name) != item['sha256']:
            raise ValueError(('accepted operand hash changed', str(path), name))
    if digest(ROOT/record['publication']['path']) != record['publication']['sha256']:
        raise ValueError(('accepted publication changed', str(path)))
    return record, base


def unpickle(path):
    with path.open('rb') as stream:
        return pickle.load(stream)


def load():
    fourier, base = accepted(CHECKPOINT)
    numerical, nb = accepted(NUMERICAL)
    source, sb = accepted(SOURCE)
    assembly, ab = accepted(native.ASSEMBLY)
    if fourier['nonzeroNormalizedResidualScalars'] or fourier['nonzeroProofResidualScalars']:
        raise ValueError('source factorization is not certified')
    original_base = Path(fourier['sourceRunDirectory'])
    if digest(original_base/'checks.json') != fourier['sourceChecksSha256']:
        raise ValueError('original factorization inventory changed')
    original_checks = json.loads((original_base/'checks.json').read_text())
    for item in original_checks['savedRows']:
        if digest(original_base/item['path']) != item['sha256']:
            raise ValueError('original factorization row hash changed')
    for item in fourier['certificateArtifacts']:
        if digest(base/item['path']) != item['sha256']:
            raise ValueError('saved source reconstruction certificate changed')
    engine_name = str(engine.HERE.relative_to(ROOT))
    for name, sha in fourier['sourceFiles'].items():
        if digest(base/'source'/name) != sha:
            raise ValueError(('accepted Fourier snapshot changed', name))
        if name != engine_name and digest(ROOT/name) != sha:
            raise ValueError(('consumed Fourier helper changed', name))
    before = ast.parse((base/'source'/engine_name).read_text())
    after = ast.parse(engine.HERE.read_text())
    additions = [n for n in after.body if getattr(n, 'name', None) == 'BoundedSourceFourierQuadrature']
    if len(additions) != 1:
        raise ValueError('expected exactly one source-quadrature class addition')
    after.body.remove(additions[0])
    if ast.dump(before) != ast.dump(after):
        raise ValueError('engine changed beyond source quadrature addition')
    packet = unpickle(Path(fourier['sourcePacket']['path']))
    if digest(Path(fourier['sourcePacket']['path'])) != fourier['sourcePacket']['sha256']:
        raise ValueError('accepted source factor packet changed')
    reduced = unpickle(sb/'reduced-action.pickle'); actions = unpickle(sb/'actions.pickle')
    assembled = unpickle(ab/'assembly.pickle')
    if (actions['reducedActionSha256'] != digest(sb/'reduced-action.pickle') or
            assembled['provenance'] != assembly['provenance'] or
            assembled['provenance']['ACTIONS_SHA256'] != digest(sb/'actions.pickle') or
            packet['provenance'] != fourier['provenance'] or
            packet['provenance']['ASSEMBLY_CHECKPOINT_SHA256'] != digest(native.ASSEMBLY) or
            packet['provenance']['SOURCE_ACTION_PACKET_SHA256'] != digest(sb/'actions.pickle')):
        raise ValueError('source/assembly/factorization packet provenance join')
    r, dimensions = native.source.restore_context(reduced)
    dimensions.__dict__.update(packet['dimensionState'])
    values = {k: engine.named(body, 'VALUE') for (k,_),body in reduced['payloads'].items()}
    pencil = engine.ReducedPencil(*(values[k] for k in engine.CLOSED_KEYS), r)
    if (pencil.strong != actions['strong'] or pencil.kernel != actions['kernel'] or pencil.fields != actions['fields']):
        raise ValueError('full source action/field/kernel join')
    projection = native.input_check.check()
    expected = json.loads(native.INPUT_JOIN.read_text())
    if projection != {k:v for k,v in expected.items() if k != 'instrument'}:
        raise ValueError('approved numerical input projection changed')
    if packet['provenance']['APPROVED_INPUT_SHA256'] != digest(native.INPUT):
        raise ValueError('factorization numerical input differs')
    builder = engine.BoundedSourceFourierAssembly(r, assembled['result']['NONLOCAL_INTEGRALS'])
    joins = []
    for i,row in enumerate(packet['result']['ROWS']):
        bounded, index, limit, remaining = builder.limit_layout(assembled['result']['NONLOCAL_INTEGRALS'][i])
        joins.append(row['INDEX'] == i and row['ORIGINAL'] == assembled['result']['NONLOCAL_INTEGRALS'][i]
            and row['BOUNDED'] == bounded and row['SOURCE_LIMIT_INDEX'] == index
            and row['SOURCE_LIMIT'] == limit and row['REMAINING_LIMITS'] == remaining
            and all(f['SOURCE_INTEGRAL'] == sp.Integral(f['SOURCE'], limit) for f in row['FACTORS']))
    if len(joins) != len(assembled['result']['NONLOCAL_INTEGRALS']) or not all(joins):
        raise ValueError('full source factor/limit/ordered momentum joins')
    tests = [unpickle(nb/f'test-{i}-operands.pickle') for i in range(2)]
    saved_actions = unpickle(nb/'actions.pickle')['result']
    if tuple((t['width'],t['momentum']) for t in tests) != tuple(saved_actions['tests']):
        raise ValueError('saved test-field parameter join')
    provenance = {'FOURIER_CHECKPOINT_SHA256': digest(CHECKPOINT),
        'FACTOR_PACKET_SHA256': fourier['sourcePacket']['sha256'],
        'CERTIFICATE_PACKET_SHA256': fourier['artifacts']['reconstruction-certificates.pickle']['sha256'],
        'NUMERICAL_CHECKPOINT_SHA256': digest(NUMERICAL), 'ASSEMBLY_CHECKPOINT_SHA256': digest(native.ASSEMBLY),
        'SOURCE_ACTION_CHECKPOINT_SHA256': digest(SOURCE), 'APPROVED_INPUT_SHA256': digest(native.INPUT),
        'INPUT_PROJECTION_SHA256': digest(native.INPUT_JOIN),
        'TEST_PACKET_SHA256': [digest(nb/f'test-{i}-operands.pickle') for i in range(2)]}
    return pencil, assembled['result'], packet['result'], tests, saved_actions, provenance, joins


def bind_sources(pencil, assembly, factors, tests, numerical):
    adapter = engine.NumericalReducedAction(pencil, assembly, json.loads(native.INPUT.read_text()))
    r = pencil.r; momenta = tuple(r.normal_map[g[2]] for g in r.momentum_groups)
    settings = numerical['results'][0]['settings']
    momentum_bound = sp.Rational(str(settings['momentumBound']))
    domains = {k: (-momentum_bound, momentum_bound) for k in momenta}
    bounds = (float(settings['sourceBound']), 1.5*float(settings['sourceBound']))
    records, joins = [], []
    for ti,test in enumerate(tests):
        field = lambda z: engine.memo_xreplace(test['field'], {r.z: z})
        if field(r.z) != sp.exp(-(r.z/test['width'])**2+sp.I*test['momentum']*r.z):
            raise ValueError('saved Gaussian ansatz differs')
        mapping = {p: field for p in pencil.probes}
        cached_integrals = {g for row in test['assembled'] for _,g in row['NONLOCAL']}
        for row in factors['ROWS']:
            current = adapter.bind(engine.dag_substitute(row['ORIGINAL'], mapping))
            joins.append(current in cached_integrals)
        for si,source_integral in enumerate(factors['SOURCE_INTEGRALS']):
            uses = [(row['INDEX'],fi,f) for row in factors['ROWS'] for fi,f in enumerate(row['FACTORS'])
                    if f['SOURCE_INTEGRAL'] == source_integral]
            f = uses[0][2]
            if any(u[2]['SOURCE'] != f['SOURCE'] or u[2]['AMPLITUDE'] != f['AMPLITUDE']
                   or u[2]['FREQUENCY'] != f['FREQUENCY'] for u in uses):
                raise ValueError('source deduplication changed a factor operand')
            amplitude = adapter.bind(engine.dag_substitute(f['AMPLITUDE'], mapping))
            bound_source = adapter.bind(engine.dag_substitute(f['SOURCE'], mapping))
            frequency = adapter.bind(f['FREQUENCY'])
            character = adapter.bind(f['CHARACTER'])
            info = engine.BoundedSourceFourierQuadrature.affine_range(frequency, domains)
            if engine.dag_free_symbols(amplitude)-{r.zp} or engine.dag_free_symbols(bound_source)-{r.zp,*momenta}:
                raise ValueError('source binding leaves unresolved numerical parameters')
            lo,hi = map(float, info['bounds'])
            sample = list(np.linspace(lo, hi, 65))
            sample.extend(x for x in (0.,float(test['momentum'])) if lo <= x <= hi)
            frequencies = np.asarray(sorted(set(sample)), dtype=float)
            # A line between the two extremal box corners attains every value
            # of this scalar affine map. Retain the actual assignments.
            assignments = []
            for nu in frequencies:
                t = 0. if hi == lo else (nu-lo)/(hi-lo)
                point = []
                for k in momenta:
                    a,b = map(float, domains[k]); c = info['coefficients'][k]
                    if c < 0: a,b = b,a
                    point.append(a+t*(b-a) if c != 0 else (a+b)/2)
                assignments.append(point)
            evaluate_frequency = sp.lambdify(momenta, frequency, 'numpy')
            range_residual = np.asarray([float(evaluate_frequency(*p))-nu for p,nu in zip(assignments,frequencies)])
            records.append({'test': ti, 'sourceIndex': si, 'uses': [(i,j) for i,j,_ in uses],
                'originalSourceIntegral': source_integral, 'symbolicAmplitude': f['AMPLITUDE'],
                'symbolicFrequency': f['FREQUENCY'], 'boundAmplitude': amplitude,
                'boundSource': bound_source, 'boundCharacter': character, 'frequency': frequency,
                'range': info, 'frequencies': frequencies, 'assignments': np.asarray(assignments),
                'assignmentResidual': range_residual,
                'amplitudeUnit': engine.PHYSICAL_METADATA.dimensions.measure(f['AMPLITUDE']),
                'integralUnit': engine.PHYSICAL_METADATA.dimensions.measure(source_integral)})
    if not all(joins):
        raise ValueError('bound original integrals do not join accepted numerical test packets')
    return {'records': records, 'momenta': momenta, 'momentumDomains': domains, 'sourceBounds': bounds,
        'orders': (32,64,128), 'workspaceBytes': 16*1024*1024, 'nativeTestIntegralJoins': joins,
        'testWidthsMomenta': [(t['width'],t['momentum']) for t in tests],
        'profileWidth': float(adapter.bind(r.ell))}


def evaluate_record(record, bound, limit_index, pencil):
    r = pencil.r; nparts = bound['orders']; limit = bound['sourceBounds'][limit_index]
    width = float(bound['testWidthsMomenta'][record['test']][0]); ell = bound['profileWidth']
    points = sorted({-limit,limit,0.,*(x for x in (-ell,ell,-width,width) if -limit < x < limit)})
    q = engine.BoundedSourceFourierQuadrature(r.zp, record['boundAmplitude'], workspace_bytes=bound['workspaceBytes'])
    values = [q.gauss(record['frequencies'], points, n) for n in nparts]
    adaptive, errors = q.adaptive(record['frequencies'], points)
    nodes, weights = q.rule(points, nparts[-1])
    # Independently evaluate the literal original bound source operand at the
    # retained momentum assignments. Do not synthesize this operand from its
    # amplitude or substitute a scalar frequency into the answer.
    direct_fn = sp.lambdify((r.zp,*bound['momenta']), record['boundSource'], 'numpy', cse=True)
    direct = np.asarray([np.dot(np.broadcast_to(np.asarray(direct_fn(nodes,*p),dtype=complex),nodes.shape),weights)
                         for p in record['assignments']])
    mutated = q.gauss(record['frequencies'], points, nparts[-1], weight_scale=1.001)
    return {'test': record['test'], 'sourceIndex': record['sourceIndex'], 'limitIndex': limit_index,
        'points': points, 'orders': nparts, 'frequencies': record['frequencies'],
        'gaussValues': np.asarray(values), 'adaptiveValues': adaptive, 'adaptiveEstimates': errors,
        'originalSourceValues': direct, 'sourceResidual': values[-1]-direct,
        'adaptiveResidual': values[-1]-adaptive, 'refinements': np.diff(np.asarray(values),axis=0),
        'measureMutationValues': mutated, 'measureMutationResidual': mutated-values[-1],
        'peakWorkspaceBytes': q.peak_workspace_bytes}


def emit_result(result, bound, pencil, provenance):
    dimensions = engine.PHYSICAL_METADATA.dimensions
    zero = dimensions.zero; zunit = dimensions.measure(pencil.r.zp); kunit = tuple(-v for v in zunit)
    metadata = engine.FullPencilModes.__new__(engine.FullPencilModes); metadata.r = pencil.r
    def numeric(suffix, value, unit, literal=False):
        body = number(value)
        engine.emit(PREFIX+'_'+suffix, body if literal else metadata.compact_fingerprint(body))
        engine.emit('METADATA_'+PREFIX+'_'+suffix, metadata.numeric_metadata(body, unit))
    engine.physical(PREFIX+'_PROVENANCE', provenance)
    engine.physical(PREFIX+'_MOMENTUM_VARIABLES', bound['momenta'])
    numeric('MOMENTUM_BOUNDS', list(bound['momentumDomains'].values()), lambda p: kunit)
    numeric('SOURCE_BOUNDS', bound['sourceBounds'], lambda p: zunit)
    numeric('TEST_WIDTHS_MOMENTA', bound['testWidthsMomenta'], lambda p: zunit if p[1] == 0 else kunit)
    numeric('PROFILE_WIDTH', bound['profileWidth'], lambda p: zunit)
    numeric('WORKSPACE_BYTES', bound['workspaceBytes'], lambda p: zero, True)
    numeric('NATIVE_TEST_INTEGRAL_JOINS', bound['nativeTestIntegralJoins'], lambda p: zero, True)
    for record in bound['records']:
        suffix = str(record['test'])+'_'+str(record['sourceIndex'])
        engine.fingerprinted(PREFIX+'_SYMBOLIC_SOURCE_'+suffix, record['originalSourceIntegral'])
        engine.physical(PREFIX+'_SYMBOLIC_FREQUENCY_'+suffix, record['symbolicFrequency'])
        numeric('ROW_FACTOR_USES_'+suffix, record['uses'], lambda p: zero, True)
        numeric('FREQUENCIES_'+suffix, record['frequencies'], lambda p: kunit)
        numeric('MOMENTUM_ASSIGNMENTS_'+suffix, record['assignments'], lambda p: kunit)
        numeric('ASSIGNMENT_RESIDUAL_'+suffix, record['assignmentResidual'], lambda p: kunit, True)
        numeric('AFFINE_RANGE_RESIDUAL_'+suffix, record['range']['residual'], lambda p: kunit, True)
        numeric('AFFINE_RANGE_BOUNDS_'+suffix, record['range']['bounds'], lambda p: kunit)
    for item in result['evaluations']:
        suffix = '_'.join(str(item[k]) for k in ('test','sourceIndex','limitIndex'))
        record = next(r for r in bound['records'] if (r['test'],r['sourceIndex']) == (item['test'],item['sourceIndex']))
        unit = record['integralUnit']
        numeric('PANELS_'+suffix, item['points'], lambda p: zunit)
        numeric('ORDERS_'+suffix, item['orders'], lambda p: zero, True)
        numeric('PEAK_WORKSPACE_BYTES_'+suffix, item['peakWorkspaceBytes'], lambda p: zero, True)
        for name in ('gaussValues','adaptiveValues','adaptiveEstimates','originalSourceValues','sourceResidual',
                     'adaptiveResidual','refinements','measureMutationValues','measureMutationResidual'):
            numeric(name.upper()+'_'+suffix, item[name], lambda p: unit,
                    literal=name in ('sourceResidual','adaptiveResidual','refinements','measureMutationResidual'))
    for item in result['domainChanges']:
        record = next(r for r in bound['records'] if (r['test'],r['sourceIndex']) == (item['test'],item['sourceIndex']))
        numeric('SOURCE_INTERVAL_CHANGE_'+str(item['test'])+'_'+str(item['sourceIndex']), item['difference'],
                lambda p: record['integralUnit'], True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--resume-from',type=Path)
    args = parser.parse_args(); base=args.run_directory.resolve(); base.relative_to(STORE)
    base.mkdir(parents=True,exist_ok=False); started=time.monotonic()
    pins={str(p.relative_to(ROOT)):digest(p) for p in SOURCES}
    for p in SOURCES:
        dest=base/'source'/p.relative_to(ROOT); dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(p,dest)
    pencil,assembly,factors,tests,numerical,provenance,joins=load()
    # Include all accepted helper snapshots in the new source closure.
    prior=json.loads(CHECKPOINT.read_text())
    for name in prior['sourceFiles']:
        path=ROOT/name
        if name not in pins:
            pins[name]=digest(path); dest=base/'source'/name; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(path,dest)
    preflight={'sourceFiles':pins,'provenance':provenance,'nativeEngineAstJoin':True,'originalLimitJoins':joins}
    save(base/'preflight.json',preflight)
    bound=bind_sources(pencil,assembly,factors,tests,numerical)
    atomic_pickle(base/'bound-sources.pickle',{'result':bound,'provenance':provenance,
        'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    bound_hash=digest(base/'bound-sources.pickle'); progress(base,'sources_bound',sources=len(bound['records']))
    destination=base/'integrals'; destination.mkdir(); inventory=[]; completed={}
    if args.resume_from:
        previous=args.resume_from.resolve(); previous.relative_to(STORE)
        if json.loads((previous/'preflight.json').read_text()) != preflight:
            raise ValueError('resume source/provenance differs; explicit instrument repair join required')
        old=unpickle(previous/'bound-sources.pickle')
        if old['provenance'] != provenance or len(old['result']['records']) != len(bound['records']) or any(
                a[k] != b[k] for a,b in zip(old['result']['records'],bound['records'])
                for k in ('originalSourceIntegral','boundAmplitude','boundSource','frequency')):
            raise ValueError('saved numerical source operand differs')
        for item in json.loads((previous/'integral-inventory.json').read_text()):
            path=previous/item['path']
            if digest(path) != item['sha256']: raise ValueError('saved quadrature record hash changed')
            record=unpickle(path); key=tuple(record[k] for k in ('test','sourceIndex','limitIndex'))
            expected=next((r for r in bound['records'] if (r['test'],r['sourceIndex'])==key[:2]),None)
            if (expected is None or key in completed or key[2] not in range(len(bound['sourceBounds']))
                    or not np.array_equal(record['frequencies'],expected['frequencies'])
                    or tuple(record['orders']) != tuple(bound['orders'])):
                raise ValueError('saved quadrature record identity/grid changed')
            completed[key]=record; shutil.copyfile(path,base/item['path']); inventory.append(item)
        save(base/'resume.json',{'runDirectory':str(previous),'records':len(completed)})
    save(base/'integral-inventory.json',inventory)
    evaluations=[]
    for record in bound['records']:
        for bi in range(len(bound['sourceBounds'])):
            key=(record['test'],record['sourceIndex'],bi)
            if key in completed:
                value=completed[key]
            else:
                value=evaluate_record(record,bound,bi,pencil)
                path=destination/('-'.join(map(str,key))+'.pickle'); atomic_pickle(path,value)
                inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)})
                save(base/'integral-inventory.json',inventory)
                progress(base,'integral_saved',test=key[0],source=key[1],interval=key[2])
            evaluations.append(value)
    changes=[{'test':record['test'],'sourceIndex':record['sourceIndex'],
        'difference':evaluations[2*i+1]['gaussValues'][-1]-evaluations[2*i]['gaussValues'][-1]}
        for i,record in enumerate(bound['records'])]
    result={'evaluations':evaluations,'domainChanges':changes}
    atomic_pickle(base/'quadrature.pickle',{'result':result,'provenance':provenance,'boundPacketSha256':bound_hash,
        'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    before=digest(base/'quadrature.pickle')
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result,bound,pencil,provenance)
        keys={tag:'s11cd'+''.join(w.title() for w in tag.removeprefix('PY_S11CD_').split('_'))
              for tag in engine.EMISSION_LINES if not tag.startswith('PY_S11CD_METADATA_')}
        engine.physical(PREFIX+'_WRITE_KEYS',keys)
        index=engine.emission_index(engine.EMISSION_LINES)
        zero_units={p:(0,0,0) for p,_ in engine.leaves(engine.cas(index))}
        engine.physical(PREFIX+'_EMISSION_LINES',index,zero_dimensions=zero_units)
    entries={}
    for line in decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ')
        if tag in entries: raise ValueError('duplicate source quadrature emission')
        entries[tag]=_restore(body)
    seen=set(); original=engine.emit
    def compare(name,value):
        tag='PY_S11CD_'+name
        if tag in seen or entries.get(tag) != engine.cas(value): raise ValueError(('source quadrature replay mismatch',tag))
        seen.add(tag)
    engine.emit=compare
    try:
        emit_result(result,bound,pencil,provenance)
        engine.physical(PREFIX+'_WRITE_KEYS',keys)
        engine.physical(PREFIX+'_EMISSION_LINES',index,zero_dimensions=zero_units)
    finally:
        engine.emit=original
    if seen != set(entries) or len(keys) != len(set(keys.values())) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('source quadrature emission/write-key census')
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    metadata_paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'): continue
        for record in body:
            if len(record)==2 and isinstance(record[0],sp.Tuple):
                fields={str(k):v for k,v in record[1]}; count=1
            else:
                fields={str(k):v for k,v in record}; count=len(fields['PATHS'])
            unit=fields['DIMENSION_L_T_M']
            if not isinstance(unit,sp.Tuple) or len(unit)!=3 or any(v.free_symbols for v in unit):
                raise ValueError(('unresolved source quadrature units',tag))
            if any(k not in fields for k in ('MULTIGRADE','EPSILON_LAMBDA_SUPPORT')):
                raise ValueError(('missing source quadrature grades',tag))
            metadata_paths+=count
    norms=[]
    for value in evaluations:
        scale=1+np.abs(value['adaptiveValues'])
        norms.append({k:value[k] for k in ('test','sourceIndex','limitIndex','peakWorkspaceBytes')} | {
            'maxSourceResidual':float(np.max(np.abs(value['sourceResidual']))),
            'maxScaledSourceResidual':float(np.max(np.abs(value['sourceResidual'])/scale)),
            'maxAdaptiveResidual':float(np.max(np.abs(value['adaptiveResidual']))),
            'maxScaledAdaptiveResidual':float(np.max(np.abs(value['adaptiveResidual'])/scale)),
            'maxAdaptiveEstimate':float(np.max(value['adaptiveEstimates'])),
            'maxRefinements':[float(np.max(np.abs(row))) for row in value['refinements']],
            'maxMutationDifference':float(np.max(np.abs(value['measureMutationResidual'])))})
    atomic_pickle(base/'dimensions-after-emission.pickle',dict(vars(engine.PHYSICAL_METADATA.dimensions)))
    summary={'runDirectory':str(base),'sourceFiles':pins,'provenance':provenance,'nativeEngineAstJoin':True,
        'originalLimitJoins':joins,'nativeTestIntegralJoins':bound['nativeTestIntegralJoins'],
        'distinctSymbolicSourceIntegrals':len(factors['SOURCE_INTEGRALS']),'boundSourceCount':len(bound['records']),
        'evaluationRecords':len(evaluations),'frequencyEvaluations':sum(len(v['frequencies']) for v in evaluations),
        'maxAssignmentResidual':max(float(np.max(np.abs(r['assignmentResidual']))) for r in bound['records']),
        'nonzeroAffineResiduals':sum(r['range']['residual'] != 0 for r in bound['records']),
        'norms':norms,'sourceIntervalChanges':[float(np.max(np.abs(r['difference']))) for r in changes],
        'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':metadata_paths,
        'boundPacketSha256BeforeEmission':bound_hash,'boundPacketSha256AfterEmission':digest(base/'bound-sources.pickle'),
        'packetSha256BeforeEmission':before,'packetSha256AfterEmission':digest(base/'quadrature.pickle'),
        'integralArtifacts':inventory,'wallSeconds':time.monotonic()-started,
        'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(base.iterdir()) if p.suffix in ('.pickle','.out')},
        'scope':'All distinct accepted bounded source integrals evaluated on both approved Gaussian tests at recorded finite source intervals and derived frequency samples. '
                'Original-source and independent adaptive comparisons, sampled refinement and finite interval changes only. '
                'No uniform interpolation bound, infinite tail, full action convergence, Abel weak limit, scattering or pole solve.'}
    save(base/'checks.json',summary)
    if pins != {name:digest(ROOT/name) for name in pins} or before != digest(base/'quadrature.pickle') or bound_hash != digest(base/'bound-sources.pickle'):
        raise ValueError('source quadrature source/packet changed during execution')
    if (summary['maxAssignmentResidual']>1e-12 or summary['nonzeroAffineResiduals'] or
            any(n['maxScaledSourceResidual']>1e-10 or n['maxScaledAdaptiveResidual']>1e-9 for n in norms) or
            max(n['maxMutationDifference'] for n in norms)<=1e-12):
        raise ValueError('source quadrature residual or sensitivity needs inspection')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
