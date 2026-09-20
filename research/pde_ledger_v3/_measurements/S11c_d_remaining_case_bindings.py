#!/usr/bin/env python3
"""Bind all actual case operands and reuse only complete source identities."""
import argparse
import copy
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import sys
import time
import numpy as np
import sympy as sp
from sympy.core.function import AppliedUndef
import S11c_d_remaining_case_factors as factors
import S11c_d_continuum_grades as grades
import S11c_d_continuum_matrices as interior

f, engine = factors.f, factors.engine
PLAN = f.M/'S11c_d_remaining_case_bindings_plan.md'
FACTOR = f.M/'S11c_d_remaining_case_factors_checkpoint.json'
BASELINE = 'LAB_HELD__RHO4_CONSTANT'


def same(left, right, memo=None):
    """Literal type/content equality, visiting each shared pair once."""
    if left is right:
        return True
    if memo is None:
        memo = {}
    key = (id(left), id(right))
    if key in memo:
        return memo[key][2]
    if isinstance(left, sp.Basic) and isinstance(right, sp.Basic):
        if type(left) is not type(right):
            return False
        a, b = left._hashable_content(), right._hashable_content()
        answer = len(a) == len(b) and all(same(x, y, memo) for x, y in zip(a, b))
    elif isinstance(left, dict) and isinstance(right, dict):
        answer = left.keys() == right.keys() and all(same(left[k], right[k], memo) for k in left)
    elif isinstance(left, (tuple, list)) and isinstance(right, (tuple, list)):
        answer = type(left) is type(right) and len(left) == len(right) and all(same(x, y, memo) for x, y in zip(left, right))
    else:
        answer = bool(left == right)
    # Keep references alive: temporary hashable-content tuples must not recycle
    # an id while this comparison's memo is still in use.
    memo[key] = (left, right, answer)
    return answer


def copy_input(path, target, expected, inputs):
    f.require(f.digest(path) == expected, ('accepted input hash', str(path)))
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
    f.require(f.digest(target) == expected, 'byte-identical accepted operand copy')
    inputs[str(path)] = expected
    return f.unpickle(target)


def load(base):
    cp = json.loads(FACTOR.read_text()); origin = Path(cp['runDirectory'])
    f.require(cp['status'] == 'PUBLISHED_ANNEX_VERIFIED', 'accepted all-case factors')
    publication = f.ROOT/cp['publication']['path']
    f.require(publication.is_symlink() and f.digest(publication) == cp['publication']['sha256'], 'actual factor publication')
    f.require(f.digest(origin/'checks.json') == cp['checksSha256'], 'completed factor checks')
    for name, digest in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(origin/'source'/name) == digest, ('current/frozen factors source', name))
    inputs = {str(origin/'checks.json'): f.digest(origin/'checks.json')}
    cases = {}
    for label in cp['caseFactorPackets']:
        packets = {}
        for kind in ('reduced-action', 'actions', 'assembly'):
            name = 'accepted-cases/'+label+'/'+kind+'.pickle'
            packets[kind] = copy_input(origin/name, base/name, cp['artifacts'][name]['sha256'], inputs)
        name = 'cases/'+label+'/factorization.pickle'
        packets['factorization'] = copy_input(origin/name, base/name, cp['artifacts'][name]['sha256'], inputs)
        cases[label] = packets
    catalogue = copy_input(origin/'integral-catalogue.pickle', base/'integral-catalogue.pickle',
                           cp['artifacts']['integral-catalogue.pickle']['sha256'], inputs)
    grade, gc, gp, gcp = grades.packet('S11c_d_continuum_grade_checkpoint.json', 'continuum-grades.pickle')
    copy_input(gp, base/'accepted-continuum-grades.pickle', gc['artifacts'][gp.name]['sha256'], inputs)
    matrix_cp = f.M/'S11c_d_continuum_matrix_checkpoint.json'
    mc = json.loads(matrix_cp.read_text())
    f.require(mc['status'] == 'PUBLISHED_ANNEX_VERIFIED', 'accepted finite matrix reuse provenance')
    prior_inputs = mc['checks']['inputPackets']
    directory = next(Path(name).parent for name in prior_inputs if name.endswith('/regulator/complete/domain-binding.pickle'))
    old = {}
    for name in ('domain-binding.pickle', 'source-binding.pickle', 'finite-system.pickle'):
        old[name] = copy_input(directory/name, base/('accepted-'+name), prior_inputs[str(directory/name)], inputs)
    old_inputs = json.loads((directory/'inputs.json').read_text())
    inputs[str(directory/'inputs.json')] = f.digest(directory/'inputs.json')
    specification = json.loads((f.M/'S11c_d_variable_profile_development_input.json').read_text())
    f.require(specification == old_inputs['input'], 'same approved physical input')
    settings = old['finite-system.pickle']['settings']
    f.require(settings['sourceBound'] == 64 and settings['momentumBound'] == 4
              and settings['profileBound'] == 14 and settings['regulator'] == 0.1,
              'retained finite settings')
    engine_name = str(engine.HERE.relative_to(f.ROOT))
    joins = factors.cases.definition_joins(Path(gc['runDirectory'])/'source'/engine_name,
                                           {'NumericalReducedAction', 'ChannelInput', 'PhysicalMetadata'})
    paths = {Path(module.__file__).resolve() for module in tuple(sys.modules.values())
             if getattr(module, '__file__', None) and Path(module.__file__).resolve().is_relative_to(f.ROOT)
             and Path(module.__file__).suffix == '.py'}
    paths.update((Path(__file__).resolve(), PLAN, FACTOR, gcp, matrix_cp, f.ACCEPTANCE,
                  f.M/'S11c_d_variable_profile_development_input.json', f.M/'S11c_d_focused_completion_plan.md'))
    paths.update(f.ROOT/name for name in cp['sourceFiles'] if name.startswith('directives/') or name.endswith('_exports.py'))
    pins = {str(path.relative_to(f.ROOT)): f.digest(path) for path in paths}
    for name in pins:
        target = base/'source'/name; target.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(f.ROOT/name, target)
    manifest = {'sourceFiles': pins, 'inputPackets': inputs, 'settings': settings, 'input': specification,
                'nativeBindingDefinitionJoins': joins, 'acceptedFiniteDirectory': str(directory),
                'scope': 'Actual four-case source bindings and independent operator grades; no new quadrature, end modes or responses.'}
    manifest['copiedInputs'] = {str(path.relative_to(base)): f.digest(path) for path in base.rglob('*.pickle')}
    f.save(base/'inputs.json', manifest)
    return cases, catalogue, grade, old, manifest


def bind_case(base, packets, old, manifest):
    r, dimensions, pencil = factors.context(packets['reduced-action'], packets['actions'], packets['assembly'])
    dimensions.__dict__.update(packets['factorization']['dimensionState'])
    assembly = packets['assembly']['result']; fourier = packets['factorization']['result']
    adapter = engine.NumericalReducedAction(pencil, assembly, manifest['input'])
    settings = manifest['settings']; prior = old['domain-binding.pickle']['bound']
    cuts = {symbol: sp.Rational(str(settings['profileBound'] if variable == r.xi else
                    settings['sourceBound'] if variable == r.zp else settings['momentumBound']))
            for variable, symbol in fourier['CUTOFFS'].items()}
    source_indices = {value: i for i, value in enumerate(fourier['SOURCE_INTEGRALS'])}
    rows = []; profiles = {}; sources = {}; jets = {}; gaussian = []; inventory = {}
    moments = tuple(r.normal_map[group[2]] for group in r.momentum_groups)
    domains = {k: (-sp.Rational(str(settings['momentumBound'])), sp.Rational(str(settings['momentumBound']))) for k in moments}
    for row in fourier['ROWS']:
        limits = tuple(engine.memo_xreplace(limit, cuts) for limit in row['REMAINING_LIMITS'])
        record = {'index': row['INDEX'], 'original': row['ORIGINAL'], 'symbolicLimits': row['REMAINING_LIMITS'],
                  'limits': limits, 'sourceLimit': engine.memo_xreplace(row['SOURCE_LIMIT'], cuts),
                  'unit': dimensions.measure(row['ORIGINAL']), 'factors': []}
        for factor in row['FACTORS']:
            value = engine.memo_xreplace(adapter.bind(factor['COEFFICIENT']), cuts)
            unit = dimensions.measure(factor['COEFFICIENT'])
            f.require(not (engine.dag_free_symbols(value)-{r.z, r.regulator, *(v[0] for v in limits)})
                      and not value.has(AppliedUndef, sp.Derivative, sp.Subs), 'fully bound actual momentum coefficient')
            record['factors'].append({'sourceIndex': source_indices[factor['SOURCE_INTEGRAL']],
                'symbolicCoefficient': factor['COEFFICIENT'], 'coefficient': value, 'unit': unit})
            for integral in factor['COEFFICIENT'].atoms(sp.Integral):
                bound = engine.memo_xreplace(adapter.bind(integral), cuts); profile_unit = dimensions.measure(integral)
                f.require(bound not in profiles or profiles[bound] == profile_unit, 'complete nested-profile units')
                profiles[bound] = profile_unit
        path = base/'rows'/f'{row["INDEX"]:03}.pickle'; path.parent.mkdir(exist_ok=True)
        f.atomic_pickle(path, record); inventory[str(row['INDEX'])] = factors.artifact(path)
        f.save(base/'row-inventory.json', inventory)
        f.require(all(tuple(limit[1:]) == domains[limit[0]] for limit in limits)
                  and tuple(record['sourceLimit'][1:]) == (-64, 64), 'actual source and ordered momentum limits')
        rows.append(record)
    for si, integral in enumerate(fourier['SOURCE_INTEGRALS']):
        uses = [(row['INDEX'], fi, factor) for row in fourier['ROWS'] for fi, factor in enumerate(row['FACTORS'])
                if same(factor['SOURCE_INTEGRAL'], integral)]
        first = uses[0][2]
        f.require(all(same((factor['AMPLITUDE'], factor['CHARACTER'], factor['FREQUENCY']),
                           (first['AMPLITUDE'], first['CHARACTER'], first['FREQUENCY'])) for _, _, factor in uses),
                  'complete shared source amplitude, character and frequency')
        source = {'sourceIndex': si, 'originalSourceIntegral': integral, 'symbolicAmplitude': first['AMPLITUDE'],
                  'symbolicFrequency': first['FREQUENCY'], 'frequency': adapter.bind(first['FREQUENCY']),
                  'amplitudeUnit': dimensions.measure(first['AMPLITUDE']), 'integralUnit': dimensions.measure(integral)}
        jet = f.source_jets(source, adapter, r); jets[si] = jet
        info = engine.BoundedSourceFourierQuadrature.affine_range(source['frequency'], domains)
        for ti, (width, k) in enumerate(prior['tests']):
            field = lambda z: sp.exp(-(z/width)**2+sp.I*k*z)
            mapping = {probe: field for probe in pencil.probes}
            direct = adapter.bind(engine.dag_substitute(first['AMPLITUDE'], mapping))
            reconstructed = sum(value*sp.diff(field(r.zp), r.zp, n) for n, value in enumerate(jet['coefficients']))
            residual = sp.expand(direct-reconstructed)
            comparison = {'sourceIndex': si, 'test': ti, 'direct': direct, 'jet': reconstructed,
                          'rawResidual': direct-reconstructed, 'expandedResidual': residual}
            path = base/'source-pairs'/f'{si:03}-{ti}.pickle'; path.parent.mkdir(exist_ok=True)
            f.atomic_pickle(path, comparison)
            f.require(residual == 0, 'actual generic derivative source reproduces Gaussian binding')
            gaussian.append(comparison)
            sources[ti, si] = dict(source, test=ti, uses=[(ri, fi) for ri, fi, _ in uses],
                boundAmplitude=direct, boundCharacter=adapter.bind(first['CHARACTER']), range=info,
                testWidth=width, profileWidth=adapter.bind(r.ell))
    local = adapter.local_matrices(); cells = []
    lookup = integral_addresses(assembly, fourier)
    for ti, (width, k) in enumerate(prior['tests']):
        field = sp.exp(-(r.z/width)**2+sp.I*k*r.z)
        for cell in assembly['ROWS']:
            i, j = cell['ROW'], cell['COLUMN']
            cells.append({'test': ti, 'row': i, 'column': j,
                'local': sum(value[i, j]*sp.diff(field, r.z, n) for n, value in local.items()),
                'terms': [(lookup[id(integral)], adapter.bind(value)) for integral, value in cell['NONLOCAL']]})
    actual_profiles = set().union(*(factor['coefficient'].atoms(sp.Integral) for row in rows for factor in row['factors']))
    f.require(actual_profiles == set(profiles), 'every actual nested profile present')
    denominator_bases = {}; visited = set()
    def denominator_visit(node):
        if id(node) in visited:
            return
        visited.add(id(node))
        if node.is_Pow and node.exp.is_negative and node.base.has(r.regulator):
            pair = tuple(sorted(engine.dag_free_symbols(node.base) & set(moments), key=sp.default_sort_key))
            if pair:
                denominator_bases[node.base] = pair
        for child in node.args:
            denominator_visit(child)
    for row in rows:
        for factor in row['factors']:
            denominator_visit(factor['coefficient'])
    actual_pairs = set(denominator_bases.values())
    f.require(actual_pairs <= set(prior['pairs']), 'accepted concentration panels cover every actual transfer pair')
    source_integral = sp.Integral(r.abel_even, (r.transfer, -sp.oo, sp.oo))
    abel = prior['abel']; k, center = abel['momentum'], abel['center']
    transformed = source_integral.transform(r.transfer, (r.ell*(k-center), k))
    density = adapter.bind(transformed.function)
    f.require(same(source_integral, abel['sourceIntegral']) and same(transformed, abel['transformed'])
              and same(density, abel['density']), 'actual new-case Abel source, measure and density reuse')
    half_height = sp.expand(sp.together(density.subs(k, center+abel['width'])-density.subs(k, center)/2).as_numer_denom()[0])
    f.require(half_height == 0, 'retained width satisfies actual density half-height equation')
    abel_join = {'sourceIntegral': source_integral, 'transformed': transformed, 'density': density,
                 'width': abel['width'], 'halfHeightResidual': half_height,
                 'actualDenominatorPairs': denominator_bases, 'actualPairs': tuple(sorted(actual_pairs, key=str)),
                 'samplingPairs': prior['pairs']}
    f.atomic_pickle(base/'abel-source-join.pickle', abel_join)
    equation_units = [packets['actions']['columnUnits'][0, i] for i in range(5)]
    field_units = [dimensions.known[value] for value in packets['actions']['fields']]
    bound = {'rows': rows, 'sources': sources, 'cells': cells, 'positions': prior['positions'], 'tests': prior['tests'],
             'pairs': prior['pairs'], 'abel': prior['abel'], 'profiles': [{'bound': v, 'unit': u} for v, u in profiles.items()],
             'profileUnits': profiles, 'equationUnits': equation_units, 'momentumUnit': dimensions.measure(moments[0]),
             'cutoffBindings': cuts}
    result = {'bound': bound, 'jets': jets, 'local': local, 'fieldUnits': field_units,
              'equationUnits': equation_units, 'settings': settings, 'gaussianComparisons': gaussian,
              'abelSourceJoin': abel_join, 'dimensionState': dict(vars(dimensions))}
    f.atomic_pickle(base/'binding.pickle', result)
    for row in rows:
        f.require(len({jets[value['sourceIndex']]['column'] for value in row['factors']}) == 1, 'actual input field per row')
    return result, (r, dimensions, pencil, adapter)


def integral_addresses(assembly, fourier):
    buckets = {}
    for row in fourier['ROWS']:
        buckets.setdefault(hash(row['ORIGINAL']), []).append(row)
    result = {}
    for cell in assembly['ROWS']:
        for integral, _ in cell['NONLOCAL']:
            if id(integral) not in result:
                matches = [row['INDEX'] for row in buckets.get(hash(integral), ()) if same(integral, row['ORIGINAL'])]
                f.require(len(matches) == 1, 'unique exact whole-integral case address')
                result[id(integral)] = matches[0]
    return result


def signature(row, binding):
    bound, jets = binding['bound'], binding['jets']
    values = []
    for factor in row['factors']:
        source = bound['sources'][0, factor['sourceIndex']]; jet = jets[factor['sourceIndex']]
        values.append((factor['symbolicCoefficient'], factor['coefficient'], factor['unit'],
            source['originalSourceIntegral'], source['symbolicAmplitude'], source['frequency'],
            jet['column'], jet['probe'], tuple(jet['coefficients']), jet['amplitudeUnit'], jet['integralUnit']))
    return (row['original'], row['symbolicLimits'], row['limits'], row['sourceLimit'], row['unit'], tuple(values))


def operator_grades(base, packets, binding, context, cache):
    r, dimensions, pencil, adapter = context
    records = {}; reused = {}; controls = []; inventory = {}
    generators = engine.PHYSICAL_METADATA.generators
    chosen = grades.specs(packets['assembly']['result'], packets['factorization']['result'],
                         binding['fieldUnits'], binding['equationUnits'], r)
    for key, value, unit, address in chosen:
        f.require(unit is not None and len(unit) == 3, 'actual coefficient unit')
        lookup = (hash(value), tuple(unit))
        matches = [(origin, record) for origin, record in cache.get(lookup, ()) if same(record['ORIGINAL'], value)]
        if matches:
            origin, record = matches[0]; reused[key] = origin
            f.require(tuple(record['UNIT']) == tuple(unit), 'grade proof reuse unit')
        else:
            record = grades.split(value, generators, unit); grades.check(record)
            cache.setdefault(lookup, []).append((str(base/'records'/(key+'.pickle')), record))
        item = {'address': address, 'record': record}; records[key] = item
        path = base/'records'/(key+'.pickle'); path.parent.mkdir(exist_ok=True); f.atomic_pickle(path, item)
        inventory[key] = {'address': address, **factors.artifact(path), 'reusedFrom': reused.get(key)}
        f.save(base/'record-inventory.json', inventory)
        if record['COEFFICIENTS']:
            grade, coefficient = next(iter(record['COEFFICIENTS'].items()))
            mutation = coefficient*sp.prod(v**n for v, n in zip(generators, grade))
            f.require(mutation != 0, 'actual coefficient omission responds'); controls.append((key, grade, mutation))
    assembly = packets['assembly']['result']; fourier = packets['factorization']['result']
    addresses = integral_addresses(assembly, fourier)
    # Supply literally equal, already checked objects from the same graph to
    # the unchanged native dictionary lookup; avoid exponential re-traversal.
    canonical = dict(assembly, ROWS=[dict(cell, NONLOCAL=tuple(
        (fourier['ROWS'][addresses[id(integral)]]['ORIGINAL'], value)
        for integral, value in cell['NONLOCAL'])) for cell in assembly['ROWS']])
    joins = grades.term_joins(canonical, fourier, records)
    result = {'records': records, 'termJoins': joins, 'generators': generators, 'reusedRecords': reused,
              'controls': controls, 'fieldUnits': binding['fieldUnits'], 'equationUnits': binding['equationUnits'],
              'dimensionState': dict(vars(dimensions)), 'finiteCutoffs': packets['factorization']['result']['CUTOFFS']}
    f.atomic_pickle(base/'operator-grades.pickle', result)
    return result


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args(); base = args.run_directory.resolve(); base.relative_to(f.STORE); base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); started = time.monotonic()
    def timeout(*_): raise TimeoutError('case binding budget; preserve all completed operands')
    signal.signal(signal.SIGALRM, timeout); signal.alarm(900)
    cases, catalogue, accepted_grades, old, manifest = load(base)
    prior = {'bound': old['domain-binding.pickle']['bound'], 'jets': old['source-binding.pickle']['jets']}
    signatures = {}
    for row in prior['bound']['rows']:
        signatures.setdefault(hash(row['original']), []).append((BASELINE, row['index'], signature(row, prior)))
    cache = {}
    for key, item in accepted_grades['records'].items():
        record = item['record']; cache.setdefault((hash(record['ORIGINAL']), tuple(record['UNIT'])), []).append(('accepted:'+key, record))
    summaries = {}; output = {}
    for label in (BASELINE, *(v for v in cases if v != BASELINE)):
        target = base/'cases'/label; target.mkdir(parents=True, exist_ok=True)
        binding, context = bind_case(target, cases[label], old, manifest)
        grade = operator_grades(target, cases[label], binding, context, cache)
        reused = []; fresh = []; joins = []
        for row in binding['bound']['rows']:
            actual = signature(row, binding)
            bucket = signatures.setdefault(hash(row['original']), [])
            matches = [(case, index) for case, index, expected in bucket if same(actual, expected)]
            if matches: reused.append({'row': row['index'], 'fromCase': matches[0][0], 'fromRow': matches[0][1]})
            else:
                fresh.append(row['index']); bucket.append((label, row['index'], actual))
            # An actual ordered-limit mutation must not reuse the original row.
            altered = copy.copy(row); limit = altered['limits'][0]
            altered['limits'] = (sp.Tuple(limit[0], limit[1], limit[2]+1), *altered['limits'][1:])
            f.require(not same(signature(altered, binding), actual), 'wrong physical limit reuse rejected')
            joins.append({'row': row['index'], 'reused': bool(matches), 'sources': [v['sourceIndex'] for v in row['factors']]})
        if label == BASELINE:
            f.require(not fresh and len(reused) == 80 and len(binding['jets']) == 35, 'full baseline numerical source replay')
            system = old['finite-system.pickle']; size = len(system['nodes']); local = np.zeros_like(system['localMatrix'])
            for n, matrix in binding['local'].items():
                for i in range(5):
                    for j in range(5):
                        local[i*size:(i+1)*size, j*size:(j+1)*size] += interior.values(matrix[i,j], context[0], system)[:,None]*system['derivativeMatrices'][n]
            residual = local-system['localMatrix']; f.atomic_pickle(target/'baseline-local-comparison.pickle', residual)
            f.require(np.max(np.abs(residual)) == 0, 'independent full accepted baseline local matrix')
        result = {'binding': binding, 'grades': grade, 'reusedRows': reused, 'newRows': fresh, 'rowJoins': joins,
                  'sourceFiles': manifest['sourceFiles'], 'inputPackets': manifest['inputPackets']}
        f.atomic_pickle(target/'case-binding.pickle', result); output[label] = str(target/'case-binding.pickle')
        summaries[label] = {'rows': len(binding['bound']['rows']), 'sources': len(binding['jets']),
            'profiles': len(binding['bound']['profileUnits']), 'nativeTerms': len(grade['termJoins']),
            'gradeRecords': len(grade['records']), 'reusedGradeRecords': len(grade['reusedRecords']),
            'newNumericalRows': fresh, 'reusableRows': reused, 'gaussianChecks': len(binding['gaussianComparisons']),
            'gradeSupport': sorted({str(g) for item in grade['records'].values() for g in item['record']['COEFFICIENTS']})}
        f.save(base/'case-inventory.json', summaries)
    for name, digest in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == digest, 'unchanged current/frozen source')
    for name, digest in manifest['inputPackets'].items(): f.require(f.digest(Path(name)) == digest, 'unchanged accepted operands')
    for name, digest in manifest['copiedInputs'].items(): f.require(f.digest(base/name) == digest, 'unchanged copied pre/post operands')
    artifacts = {str(p.relative_to(base)): factors.artifact(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
    checks = {**manifest, 'runDirectory': str(base), 'cases': summaries, 'casePackets': output, 'artifacts': artifacts,
              'newNumericalIntegrations': 0, 'newModeSolves': 0, 'wallSeconds': time.monotonic()-started,
              'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__': main()
