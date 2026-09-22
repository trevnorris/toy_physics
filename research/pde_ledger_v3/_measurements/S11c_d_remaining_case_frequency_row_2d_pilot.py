#!/usr/bin/env python3
"""One finite 2D row by new numerical factor contraction and two small grids.

The saved coefficient's literal Mul factors are partitioned, never symbolically
rewritten. The finite profile transform, source Fourier calls and contractions
are new numerical operations. Existing source actions and exact Fourier returns
are read from accepted packets. No old compiler, basis, rule or row is replayed.
"""
import argparse
import ast
import gc
import hashlib
import inspect
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import time

import S11c_d_remaining_case_frequency_rows_1d as prior
p, recovery = prior.p, prior.recovery
np, sp, io = p.np, p.sp, p.io
M, F, require, same = p.M, p.F, p.require, p.same
CP = M/'S11c_d_remaining_case_frequency_rows_1d_checkpoint.json'
CP_SHA = '2907d7da91fb048b345e0c4e4aa238a182b90e021672f81ecb453f7cf1fac2d5'
PLAN = M/'S11c_d_remaining_case_frequency_row_2d_pilot_plan.md'
CASE, ROW = 'LAB_HELD__RHOBR_CONSTANT', 46


def release(stream):
    os.posix_fadvise(stream.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''): h.update(block)
        release(stream)
    return h.hexdigest()


class Reader(io.Reader):
    def packet(self, path, expected):
        route = self.retain(path, expected)
        with Path(route['canonical']).open('rb') as stream:
            value = pickle.load(stream)
            release(stream)
        return value


class Journal(p.storage.Journal):
    def write(self, name, value):
        require(not Path(name).is_absolute() and '..' not in Path(name).parts, 'relative new output')
        record = super().write(name, value)
        # Flush only the newly written packet before advising away clean pages.
        with Path(record['path']).open('rb') as stream:
            os.fsync(stream.fileno())
            release(stream)
        return record


def partition(coefficient, k, q, z, xi):
    """Read the actual factor tree; no factor(), expand(), subs() or new Expr."""
    require(isinstance(coefficient, sp.Mul), 'literal multiplicative coefficient')
    groups = {'constant': [], 'output': [], 'input': [], 'profile': []}
    for factor in coefficient.args:
        free = factor.free_symbols
        if isinstance(factor, sp.Integral): groups['profile'].append(factor)
        elif not free: groups['constant'].append(factor)
        elif free <= {k, z}: groups['output'].append(factor)
        elif free <= {q}: groups['input'].append(factor)
        else: raise TypeError(('unseparated actual coefficient factor', repr(factor)))
    require(len(groups['profile']) == 1 and groups['output'] and groups['input'], 'single inspected profile and both momentum legs')
    integral = groups['profile'][0]
    require(len(integral.limits) == 1 and integral.limits[0][0] == xi and
            tuple(map(float, integral.limits[0][1:])) == (-14., 14.), 'actual finite profile limits')
    candidates = []
    for node in sp.preorder_traversal(integral.function):
        if not isinstance(node, sp.Add) or len(node.args) != 2 or k not in node.args: continue
        other = next(a for a in node.args if a != k)
        if isinstance(other, sp.Mul) and len(other.args) == 2 and set(other.args) == {sp.S.NegativeOne, q}:
            candidates.append(node)
    require(len(set(candidates)) == 1, 'literal output-minus-input profile subtree')
    difference = candidates[0]
    def covered(node):
        if node == difference: return
        require(node not in (k, q), 'all profile momentum occurrences inside exact difference subtree')
        for child in node.args: covered(child)
    covered(integral.function)
    require(integral.function.free_symbols <= {xi, k, q}, 'no unbound profile parameters')
    return groups, integral, difference


def trapezoid(lower, upper, panels):
    nodes = np.linspace(lower, upper, panels+1)
    weights = np.full(panels+1, (upper-lower)/panels)
    weights[[0, -1]] *= .5
    return nodes, weights


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory', type=Path, required=True)
    base = ap.parse_args().run_directory.resolve(); base.relative_to(p.REPO/'_scratch/s11c')
    base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900)
    started = time.monotonic(); io.digest = digest
    reader, journal, packets = Reader(), Journal(base), {}
    def packet(record):
        path = record.get('path', record.get('logical'))
        if path not in packets: packets[path] = reader.packet(path, record['sha256'])
        return packets[path]
    def address(route):
        value = packet(route['packet'])
        for key in route['keys']: value = value[key]
        return value
    cp = reader.json(CP, CP_SHA)
    require(cp['status'] == 'ACCEPTED_BOUNDED_CASE_FREQUENCY_1D_ROWS', 'accepted saved 1D rows')
    next_view = reader.json(cp['nextSavedInput']['path'], cp['nextSavedInput']['sha256'])
    selected = next_view['route']
    require(selected['firstOwner'] == [CASE, ROW] and selected['layoutDimension'] == 2 and
            selected['status'] == 'UNMATCHED_FULL_ROW_INPUT' and not selected['savedCompleteMatches'], 'one genuinely missing full 2D row')
    reader.retain(p.READY/'complete/checks.json', p.READY_SHA)
    inventory = reader.json(p.READY/'completed-input-artifact-inventory.json', p.INVENTORY_SHA)
    def metadata(name): return reader.json(p.READY/'complete'/name, inventory[name]['sha256'])
    require(selected == metadata('cases/'+CASE+'/row-'+str(ROW)+'.json'), 'actual next row input route')
    native = metadata('native-callers.json')
    for item in native.values():
        route = reader.retain(item['file']['logical'], item['file']['sha256'])
        tree = ast.parse(Path(route['canonical']).read_text())
        for name, body in item['bodies'].items():
            node = next(n for n in tree.body if getattr(n, 'name', None) == name)
            require(ast.dump(node) == ast.dump(ast.parse(body).body[0]), 'whole native row/source caller')
    for logical, rec in metadata('inputs.json')['consumedRoutes'].items():
        if logical.endswith(('.py', '.md')): reader.retain(logical, rec['sha256'])
    raw = packet(selected['sourceInputs']['packet']); scalars = packet(selected['scalarInputs']); system = packet(selected['basis'])
    physical = reader.json(selected['ownPhysicalRoutes']['packet']['logical'], selected['ownPhysicalRoutes']['packet']['sha256'])
    common = packet(physical['context']); context = common['contextPair'][0]
    require(same(*common['contextPair']) and same(*common['basisPair']), 'full own source context and unit basis')
    for record in physical.values():
        if isinstance(record, dict) and 'logical' in record: reader.retain(record['logical'], record['sha256'])
    row = raw['bound']['rows'][ROW]; require(row['index'] == ROW and len(row['factors']) == 1, 'entire selected single-factor row')
    factor = row['factors'][0]; si = factor['sourceIndex']; source = raw['bound']['sources'][0, si]
    jet = address(selected['jets'][str(si)]); coefficient = scalars['actual']['factor', ROW, 0]
    k, q = (lim[0] for lim in row['limits']); z, xi = context['z'], context['xi']
    positions, settings = system['nodes'], raw['settings']
    require(si == 20 and same(source['frequency'], q), 'actual source momentum identity without old lambda recreation')
    require(same(jet['originalBoundAmplitude'], scalars['actual']['source', si]) and
            same(jet['amplitudeUnit'], source['amplitudeUnit']) and same(jet['integralUnit'], source['integralUnit']), 'whole source amplitude and physical units')
    require(same(settings, system['settings']) and same(settings, context['settings']) and len(positions) == 129 and
            complex(scalars['frequency']) == 1-.01j and settings['momentumBound'] == 4 and settings['profileBound'] == 14,
            'actual fixed frequency, full trial grid and finite domains')
    source_meta = metadata('saved-prepared-basis/'+CASE+'.json'); prepared = packet(source_meta['packet']); saved = prepared['original']
    nodes, weights = saved['source_nodes'], saved['source_weights']
    basis_item = next(v for v in metadata('cases/'+CASE+'/source-basis-inputs.json') if v['sourceIndex'] == si)
    candidates = [v for v in basis_item['savedCoefficientBasisCandidates'] if v['owner'][0] == CASE]
    require(candidates and same(basis_item['coefficientRoute'], selected['jets'][str(si)]) and
            saved['size'] == 129 and nodes.shape == weights.shape == (1024,) and
            same(source_meta['sourceSettings'], json.loads(json.dumps(settings))), 'whole saved native source-rule caller and arrays')
    key = (jet['probe'], tuple(jet['coefficients']), jet['amplitudeUnit'], jet['integralUnit'])
    actions = []
    for candidate in candidates:
        oldjet = address(candidate['input'])
        require(same(key, (oldjet['probe'], tuple(oldjet['coefficients']), oldjet['amplitudeUnit'], oldjet['integralUnit'])) and
                same(address(candidate['nodeRoute']), nodes) and same(address(candidate['weightRoute']), weights), 'full saved source-basis arguments')
        actions.append(address(candidate['result']))
    action = actions[0]; action_route = candidates[0]['result']; node_route = candidates[0]['nodeRoute']
    require(action.shape == (1024, 129) and action.dtype == np.dtype(complex) and np.isfinite(action).all() and
            all(same(action, v) for v in actions), 'complete matching actual source-action returns')
    old_summary = next(v for v in cp['rows'] if v['rowIndex'] == 29); old_row = packet(old_summary['rowInput'])
    require(same(old_row['sourceAction'], action_route) and same(old_row['context'], context) and
            same(old_row['settings'], settings) and same(old_row['fieldUnits'], raw['fieldUnits']) and
            same(old_row['equationUnits'], raw['equationUnits']) and
            same(key, (old_row['jet']['probe'], tuple(old_row['jet']['coefficients']), old_row['jet']['amplitudeUnit'], old_row['jet']['integralUnit'])),
            'accepted row29 Fourier owner has complete same physical source action')
    input_record = journal.write('row-input.pickle', {'row': row, 'coefficient': coefficient, 'source': source, 'jet': jet,
        'context': context, 'settings': settings, 'positions': positions, 'physicalRoutes': physical, 'scalarInputs': selected['scalarInputs'],
        'fieldUnits': raw['fieldUnits'], 'equationUnits': raw['equationUnits'], 'abel': raw['bound']['abel'],
        'pairs': raw['bound']['pairs'], 'profileUnits': raw['bound']['profileUnits'], 'frequency': scalars['frequency'], 'sourceAction': action_route})
    journal.json('source-basis-route.json', {'matches': candidates, 'chosen': action_route, 'nodeRoute': node_route,
        'weightRoute': candidates[0]['weightRoute'], 'actualPriorFourierOwner': old_summary['rowInput'], 'newSourceActionCalls': 0})
    groups, integral, difference = partition(coefficient, k, q, z, xi)
    require(integral in raw['bound']['profileUnits'], 'actual nested profile unit declaration')
    partition_record = journal.write('literal-factor-partition.pickle', {'coefficient': coefficient, 'groups': groups,
        'profile': integral, 'differenceSubtree': difference, 'variables': (k, q, z, xi),
        'profileUnit': raw['bound']['profileUnits'][integral], 'coefficientUnit': factor['unit'], 'rowInput': input_record})
    # Compile the single accepted numerical evaluator; old scientific mains do not run.
    reader.retain(Path(p.__file__), recovery.HELPER_SHA)
    reader.retain(Path(recovery.__file__), '7d9d446165f6e3301dabd5dc225f1585930f5879a1feb65fa58e96ba24f82a40')
    module, adapter_join = recovery.adapted(Path(p.__file__).read_text())
    numerical = ast.Module(body=[module.body[0]], type_ignores=[]); namespace = dict(vars(p))
    exec(compile(numerical, '<accepted-numerical-principal-power-evaluator>', 'exec'), namespace); evaluate = namespace['evaluate']
    baseline = reader.json(M/'S11c_d_frequency_matrix_checkpoint.json', recovery.BASELINE_SHA)
    engine = M.parent/'scripts/S11c_d_mixing_scattering_sympy_audit.py'
    require(baseline['sourceFiles']['scripts/S11c_d_mixing_scattering_sympy_audit.py'] == recovery.ENGINE_SHA, 'actual accepted native compiler source')
    engine_current = reader.retain(engine, recovery.ENGINE_SHA)
    engine_frozen = reader.retain(Path(baseline['runDirectory'])/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py', recovery.ENGINE_SHA)
    for path in (Path(__file__).resolve(), PLAN, Path(prior.__file__), Path(io.__file__), Path(p.storage.__file__),
                 M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'): reader.retain(path)
    journal.json('native-and-new-numerical-callers.json', {'native': native, 'sourceBasisCaller': source_meta['nativeCaller'],
        'engineCurrent': engine_current, 'engineFrozen': engine_frozen, 'evaluatorJoin': adapter_join,
        'evaluator': ast.unparse(numerical), 'partition': inspect.getsource(partition), 'newRule': inspect.getsource(trapezoid),
        'newContractionHelper': reader.retain(Path(__file__).resolve()), 'nativeCompilerBasisRuleCalled': False})
    for path in (Path(__file__).resolve(), PLAN):
        target = base/'source'/path.name; target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as out: out.write(path.read_bytes())
        reader.retain(target, digest(path))
    scope = {'case': CASE, 'rowIndex': ROW, 'sourceIndex': si, 'frequency': {'real': 1., 'imag': -.01},
        'chartAndSheet': 'Exact accepted principal complex coefficient expression branches, finite real momentum square [-4,4]^2; own end at1-.01i unchanged.',
        'trialSize': 129, 'momentumPanels': [512, 1024], 'momentumRule': 'new uniform composite trapezoid',
        'profilePanels': 4096, 'profileComparisonPanels': 2048, 'profileBounds': [-14., 14.], 'profileBatchSize': 64,
        'profileComparisonDeltas': [0., .731, -.731, 8., -8.], 'profileAbsoluteTolerance': 1e-9,
        'sourceRule': 'exact saved1024 source nodes/weights and weighted source action, native sourceOrder256/bound64',
        'originalSettings': settings, 'coefficientCheckTolerance': 2e-12, 'selectedContractionTolerance': 2e-12,
        'rowComparisonTarget': {'relative': .01, 'absolute': 1e-4, 'purpose': 'pilot row resolution indication only, not observable accuracy'},
        'costGate': 'first full coarse grid measured; fine grid only if3x first cost plus60s fits remaining900s',
        'scope': 'One practical analog toy-model finite2D row pilot. New explicit alternative grid, not native16/4/4/profile512 rule completion. No scattering error bound, tiny-effect, pole or domain result; no automatic grid extension.'}
    journal.json('pilot-scope.json', scope)
    def forbidden(*args, **kwargs): raise RuntimeError('completed native science disabled in new2D pilot')
    prior.main = p.source_matrix = p.independent_source = forbidden
    io.native.f.source_jets = io.native.f.polynomial_basis = io.native.f.BasisMomentum.prepare_basis = forbidden
    io.native.Pair.__init__ = io.native.maps = io.native.continue_pair = forbidden
    for name in ('diff', 'lambdify', 'cancel', 'expand', 'factor', 'solve', 'gcd', 'resultant', 'integrate'): setattr(sp, name, forbidden)
    rules = {}
    def rule(name, lower, upper, panels):
        arg = journal.write(name+'/input.pickle', {'lower': lower, 'upper': upper, 'panels': panels, 'method': 'uniform composite trapezoid'})
        x, w = trapezoid(lower, upper, panels)
        rec = journal.write(name+'/value.pickle', {'nodes': x, 'weights': w}); journal.json(name+'/completed.json', {'input': arg, 'value': rec})
        mass = float(np.sum(w)); first = float(w@x)
        journal.json(name+'/measure-check.json', {'mass': mass, 'firstMoment': first, 'expectedMass': upper-lower})
        require(np.isfinite(x).all() and np.all(w>0) and abs(mass-(upper-lower)) < 1e-12 and abs(first-(upper*upper-lower*lower)/2) < 1e-12, 'new finite rule measure and orientation')
        return x, w, rec
    for n in (2048, 4096): rules[n] = rule('profile-rules/'+str(n), -14., 14., n)
    profile_cache = {}; profile_batches = []; profile_serial = 0
    def profiles(deltas, panels=4096):
        nonlocal profile_serial
        ds = [float(v) for v in deltas]; missing = sorted(set(v for v in ds if (panels, v.hex()) not in profile_cache))
        x, w, rule_ref = rules[panels]
        for offset in range(0, len(missing), 64):
            chunk = np.asarray(missing[offset:offset+64]); folder = 'profiles/'+str(profile_serial); profile_serial += 1
            arg = journal.write(folder+'/input.pickle', {'integral': integral, 'differenceSubtree': difference,
                'deltas': chunk, 'xiRule': rule_ref, 'profileUnit': raw['bound']['profileUnits'][integral], 'rowInput': input_record})
            environment = {difference: chunk[:, None], xi: x[None, :]}
            with np.errstate(over='raise', invalid='raise', divide='raise', under='ignore'):
                evaluated = np.broadcast_to(np.asarray(evaluate(integral.function, environment), complex), (len(chunk), len(x)))
            er = journal.write(folder+'/integrand-value.pickle', evaluated)
            cr = journal.write(folder+'/contraction-input.pickle', {'integrand': er, 'rule': rule_ref, 'axis': 1, 'measureIncludedOnce': True})
            value = evaluated@w
            out = journal.write(folder+'/value.pickle', value); journal.json(folder+'/completed.json', {'input': arg, 'integrand': er, 'contractionInput': cr, 'value': out})
            require(np.isfinite(value).all(), 'finite complete profile transform')
            for index, delta in enumerate(chunk): profile_cache[panels, float(delta).hex()] = (value[index], {'packet': out, 'keys': [index]})
            profile_batches.append({'input': arg, 'value': out, 'panels': panels, 'count': len(chunk)})
        return np.asarray([profile_cache[panels, v.hex()][0] for v in ds]), [profile_cache[panels, v.hex()][1] for v in ds]
    probe_deltas = scope['profileComparisonDeltas']
    coarse_profiles, ca = profiles(probe_deltas, 2048); fine_profiles, fa = profiles(probe_deltas)
    journal.json('profile-comparison/input.json', {'coarse': ca, 'fine': fa})
    profile_spread = float(np.max(abs(fine_profiles-coarse_profiles)))
    journal.json('profile-comparison/value.json', {'absoluteSpread': profile_spread, 'tolerance': 1e-9})
    require(profile_spread < 1e-9, 'focused finite profile rule comparison')
    # Read accepted Fourier inputs only where the numeric frequency and complete
    # source/unit arguments match. Other row29 values are not loaded or recomputed.
    catref = cp['savedReview']['pointCatalogue']; catalogue = reader.json(catref['path'], catref['sha256'])
    old_points = {v['momentumHex']: v for v in catalogue if v['rowIndex'] == 29}
    old_checks = reader.json(Path(cp['runDirectory'])/'checks.json', cp['checksSha256'])
    old_artifacts = old_checks['artifacts']; fourier_cache = {}; fourier_routes = []; fourier_serial = 0
    def fouriers(grid):
        nonlocal fourier_serial
        for value in grid:
            qv = float(value); h = qv.hex()
            if h in fourier_cache: continue
            request = {'frequency': complex(qv), 'sourceAction': action_route, 'sourceNodes': node_route,
                'unit': jet['integralUnit'], 'method': 'literal minus-phase vector times accepted weighted source action'}
            if h in old_points:
                owner = old_points[h]; prefix = 'rows/29/points/'+str(owner['pointIndex'])
                arg_ref = old_artifacts[prefix+'/fourier-input.pickle']; receipt_ref = old_artifacts[prefix+'/fourier-completed.json']
                old_arg = packet(arg_ref); receipt = reader.json(receipt_ref['path'], receipt_ref['sha256'])
                require(same(request, old_arg) and receipt == {'input': arg_ref, 'value': owner['operands']['sourceFourier']}, 'actual complete accepted Fourier input/return/receipt')
                out = owner['operands']['sourceFourier']; result = packet(out); disposition = 'ACCEPTED_SAVED'
            else:
                prefix = 'source-fourier/'+str(fourier_serial); fourier_serial += 1
                arg_ref = journal.write(prefix+'/input.pickle', request)
                with np.errstate(over='raise', invalid='raise', divide='raise', under='ignore'):
                    result = np.exp(-1j*qv*nodes)@action
                out = journal.write(prefix+'/value.pickle', result); journal.json(prefix+'/completed.json', {'input': arg_ref, 'value': out})
                disposition = 'NEW'
            require(result.shape == (129,) and result.dtype == np.dtype(complex) and np.isfinite(result).all(), 'full source Fourier return')
            fourier_cache[h] = (result, out)
            fourier_routes.append({'frequencyHex': h, 'input': arg_ref, 'value': out, 'disposition': disposition})
        return np.asarray([fourier_cache[float(v).hex()][0] for v in grid]), [fourier_cache[float(v).hex()][1] for v in grid]
    # Each literal coefficient factor is evaluated only at new grid coordinates.
    # Fine-grid nested entries use the actual coarse saved values and addresses.
    factor_cache = {}; factor_routes = []; constant_product = None
    def coefficient_factors(grid, label):
        nonlocal constant_product
        outputs = {}
        for group in ('constant', 'input', 'output'):
            if group == 'constant' and constant_product is not None:
                outputs[group] = constant_product
                factor_routes.append({'label': label, 'group': group, 'value': constant_product[1], 'disposition': 'EXACT_COMPLETED_PRODUCT'})
                continue
            lists = []
            for fi, expression in enumerate(groups[group]):
                if group == 'constant': wanted = [(None, None)]
                else: wanted = [(float(v).hex(), float(v)) for v in grid]
                missing = [(h,v) for h,v in wanted if (group, fi, h) not in factor_cache]
                if missing:
                    values = np.asarray([v for h,v in missing]) if group != 'constant' else None
                    environment = {} if group == 'constant' else ({q: values} if group == 'input' else {k: values[:,None], z: positions[None,:]})
                    folder = 'coefficient-factors/'+label+'/'+group+'-'+str(fi)
                    arg = journal.write(folder+'/input.pickle', {'expression': expression, 'environment': environment,
                        'group': group, 'factorIndex': fi, 'partition': partition_record, 'coefficientUnit': factor['unit']})
                    with np.errstate(over='raise', invalid='raise', divide='raise', under='ignore'):
                        value = np.asarray(evaluate(expression, environment), complex)
                        value = value.reshape(()) if group == 'constant' else np.broadcast_to(value, (len(missing),129) if group == 'output' else (len(missing),))
                    rec = journal.write(folder+'/value.pickle', value); journal.json(folder+'/completed.json', {'input': arg, 'value': rec})
                    require(np.isfinite(value).all(), 'finite actual factor values')
                    for index,(h,v) in enumerate(missing): factor_cache[group,fi,h] = (value if group == 'constant' else value[index], {'packet':rec,'keys':[] if group == 'constant' else [index]})
                lists.append([factor_cache[group,fi,h] for h,v in wanted])
            route = [[pair[1] for pair in values] for values in lists]
            folder = 'factor-products/'+label+'/'+group
            arg = journal.write(folder+'/input.pickle', {'group': group, 'factors': route, 'partition': partition_record})
            shape = () if group == 'constant' else ((len(grid),129) if group == 'output' else (len(grid),))
            value = np.ones(shape, complex)
            for values in lists: value = value * (values[0][0] if group == 'constant' else np.asarray([v[0] for v in values]))
            rec = journal.write(folder+'/value.pickle', value); journal.json(folder+'/completed.json', {'input': arg, 'value': rec})
            outputs[group] = (value, rec); factor_routes.append({'label': label, 'group': group, 'input': arg, 'value': rec})
            if group == 'constant': constant_product = (value, rec)
        return outputs
    costs = []; results = []; grids = []; summaries = []
    for panels in (512, 1024):
        label = str(panels); folder = 'momentum-grids/'+label
        if costs:
            reserve = 3*costs[0]+60; remaining = 900-(time.monotonic()-started)
            journal.json(folder+'/cost-decision.json', {'firstCompleteGridSeconds': costs[0], 'requiredReserveSeconds': reserve, 'remainingSeconds': remaining})
            require(reserve < remaining, 'measured second small grid fits guarded remaining budget')
        tick = time.monotonic(); grid, mass, grid_ref = rule(folder+'/rule', -4., 4., panels)
        ds = np.arange(-panels, panels+1)*(8./panels)
        profile, profile_refs = profiles(ds)
        source_values, source_refs = fouriers(grid)
        parts = coefficient_factors(grid, label)
        output, output_ref = parts['output']; input_value, input_ref = parts['input']; constant, constant_ref = parts['constant']
        arg = journal.write(folder+'/kernel-input.pickle', {'grid': grid_ref, 'deltas': ds, 'profileValues': profile_refs,
            'inputFactors': input_ref, 'constantFactors': constant_ref, 'differenceIndex': 'outputIndex-inputIndex+panels',
            'rowInput': input_record, 'bothMomentumMeasuresIncludedOnce': True})
        index = np.arange(panels+1)[:,None]-np.arange(panels+1)[None,:]+panels
        kernel = profile[index] * mass[:,None] * mass[None,:] * input_value[None,:] * constant
        kernel_ref = journal.write(folder+'/kernel-value.pickle', kernel); journal.json(folder+'/kernel-completed.json', {'input': arg, 'value': kernel_ref})
        source_ref = journal.write(folder+'/source-value-routes.pickle', {'grid': grid_ref, 'returns': source_refs})
        arg = journal.write(folder+'/inner-input.pickle', {'kernel': kernel_ref, 'sourceValues': source_ref, 'operation': 'kernel @ sourceValues'})
        inner = kernel@source_values
        inner_ref = journal.write(folder+'/inner-value.pickle', inner); journal.json(folder+'/inner-completed.json', {'input': arg, 'value': inner_ref})
        arg = journal.write(folder+'/row-input.pickle', {'outputFactors': output_ref, 'inner': inner_ref, 'operation': 'outputFactors.T @ inner', 'transposeConjugates': False})
        value = output.T@inner
        value_ref = journal.write(folder+'/row-value.pickle', value); journal.json(folder+'/row-completed.json', {'input': arg, 'value': value_ref})
        require(value.shape == (129,129) and np.isfinite(value).all(), 'complete finite full2D row')
        # Focused new checks consume these actual arrays; never rerun coefficient
        # evaluations, profile or source calls for contraction plumbing.
        if panels == 512:
            pairs = ((113,271),(384,63),(71,429)); direct_checks=[]
            for ci,(ki,qi) in enumerate(pairs):
                env={k:float(grid[ki]),q:float(grid[qi]),z:positions,integral:profile[ki-qi+panels]}
                prefix='coefficient-checks/'+str(ci)
                ar=journal.write(prefix+'/input.pickle',{'expression':coefficient,'environment':env,'rowInput':input_record})
                direct=np.broadcast_to(np.asarray(evaluate(coefficient,env),complex),positions.shape)
                dr=journal.write(prefix+'/direct-value.pickle',direct)
                ar2=journal.write(prefix+'/comparison-input.pickle',{'direct':dr,'output':output_ref,'input':input_ref,'constant':constant_ref,
                    'profile':profile_refs[ki-qi+panels],'indices':(ki,qi)})
                grouped=output[ki]*input_value[qi]*constant*profile[ki-qi+panels]
                gr=journal.write(prefix+'/grouped-value.pickle',grouped)
                spread=float(np.max(abs(direct-grouped))/(1+np.max(abs(direct))))
                journal.json(prefix+'/completed.json',{'input':ar,'direct':dr,'comparisonInput':ar2,'grouped':gr,'scaledDifference':spread})
                require(spread<2e-12,'literal whole coefficient versus factor contraction')
                direct_checks.append(spread)
            ri=(0,64,128); cj=(0,11,128)
            ar=journal.write('contraction-check/input.pickle',{'kernel':kernel_ref,'sourceValues':source_ref,'output':output_ref,'row':value_ref,
                'positionIndices':ri,'columnIndices':cj,'method':'literal double summation of nine selected entries'})
            direct=np.asarray([[np.sum(output[:,r,None]*kernel*source_values[None,:,c]) for c in cj] for r in ri])
            dr=journal.write('contraction-check/direct-value.pickle',direct)
            chosen=value[np.ix_(ri,cj)]; spread=float(np.max(abs(direct-chosen))/(1+np.max(abs(chosen))))
            # Changed measure and transposed-kernel controls use saved factors.
            mutated_measure=direct*1.001
            mutated_index=np.asarray([[np.sum(output[:,r,None]*kernel.T*source_values[None,:,c]) for c in cj] for r in ri])
            mr=journal.write('contraction-check/mutation-values.pickle',{'changedMeasure':mutated_measure,'transposedKernel':mutated_index})
            responses={'measure':float(np.max(abs(mutated_measure-direct))), 'kernelTranspose':float(np.max(abs(mutated_index-direct)))}
            journal.json('contraction-check/completed.json',{'input':ar,'direct':dr,'mutations':mr,'scaledDifference':spread,'responses':responses})
            require(spread<2e-12 and all(v>1e-14 for v in responses.values()),'selected contraction, measure and momentum-index controls')
        cost=time.monotonic()-tick;costs.append(cost);results.append(value);grids.append(grid_ref)
        summary={'panels':panels,'rowInput':input_record,'grid':grid_ref,'rowReturn':value_ref,'wallSeconds':cost,
            'rowMaxNorm':float(np.max(abs(value))),'newFourierCallsCumulative':fourier_serial,'profileBatchesCumulative':profile_serial}
        journal.json(folder+'/summary.json',summary);summaries.append(summary)
        del kernel,inner,index;gc.collect()
    journal.json('row-comparison/input.json',{'coarse':summaries[0]['rowReturn'],'fine':summaries[1]['rowReturn']})
    delta=results[1]-results[0];norm=float(np.max(abs(results[1])));absolute=float(np.max(abs(delta)))
    comparison={'rowMaxNorm':norm,'absoluteDifference':absolute,'relativeDifference':absolute/norm if norm else None,
        'pilotTargetMet':absolute<=max(1e-4,.01*norm),'targetAbsolute':1e-4,'targetRelative':.01,
        'scope':'Fixed finite-setting row refinement, not scattering/source-rule error bound. No automatic extension if unresolved.'}
    journal.write('row-comparison/difference.pickle',delta);journal.json('row-comparison/value.json',comparison)
    journal.json('operation-catalogue.json',{'sourceFourier':fourier_routes,'profileBatches':profile_batches,'factorProducts':factor_routes})
    reader.postcheck()
    for record in journal.artifacts.values():require(digest(record['path'])==record['sha256'],'unchanged new input/value/receipt bytes')
    journal.json('inputs.json',{'acceptedRows1D':reader.retain(CP,CP_SHA),'sourceAction':action_route,'consumedRoutes':reader.routes})
    checks={'status':'COMPLETED_BOUNDED_2D_ROW_PILOT','case':CASE,'rowIndex':ROW,'sourceIndex':si,'rows':summaries,
        'comparison':comparison,'profileComparisonAbsolute':profile_spread,'newSourceActions':0,'newNativeBasisCompilerRuleEndMapCalls':0,
        'newSourceFourierCalls':fourier_serial,'acceptedSavedFourierCalls':sum(v['disposition']=='ACCEPTED_SAVED' for v in fourier_routes),
        'uniqueProfileCalls':len(profile_cache),'profileBatches':profile_serial,'allConsumedHashesUnchanged':True,
        'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-started,'scope':scope['scope']}
    journal.json('checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__ == '__main__': main()
