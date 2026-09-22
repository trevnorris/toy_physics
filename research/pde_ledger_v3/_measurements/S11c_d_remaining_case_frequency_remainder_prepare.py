#!/usr/bin/env python3
"""New numerical source actions and finite remainder-profile preparation."""
import argparse
import ast
import gc
import inspect
import json
from pathlib import Path
import resource
import signal
import time

import scipy.special as special
import S11c_d_remaining_case_frequency_remainder_inputs_finish as inputs

old = inputs.original
saved, p, sp, io = old.saved, old.p, old.sp, old.io
np, M, F, require, same = saved.np, old.M, old.F, old.require, old.same
NAME = 'S11c_d_remaining_case_frequency_remainder_prepare'
READY = F/'remainder-inputs-recovery-03/complete'
READY_SHA = 'ecd3be27aedb8f6fa99744b1da51ed4e80a0a555175ee0464092cc8d4a15fbde'


def evaluator():
    module, join = saved.recovery.adapted(Path(p.__file__).read_text())
    namespace = dict(vars(p))
    exec(compile(ast.Module(body=[module.body[0]], type_ignores=[]), '<accepted-principal-numerical-evaluator>', 'exec'), namespace)
    prior = namespace['evaluate']
    def evaluate(expression, environment, memo=None):
        memo = {} if memo is None else memo
        if expression in memo: return memo[expression]
        if expression.func is sp.Heaviside:
            require(len(expression.args) == 2, 'actual native Heaviside origin value')
            argument = np.asarray(evaluate(expression.args[0], environment, memo))
            origin = np.asarray(evaluate(expression.args[1], environment, memo))
            require(np.all(argument.imag == 0) and np.all(origin.imag == 0), 'actual real Heaviside arguments')
            value = np.heaviside(argument.real, origin.real)
            memo[expression] = value
            return value
        return prior(expression, environment, memo)
    namespace['evaluate'] = evaluate
    return evaluate, {'acceptedEvaluator': ast.unparse(module.body[0]), 'acceptedReverseAST': join,
                      'onlyNewNode': 'Heaviside actual argument and origin value -> numpy.heaviside'}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory', required=True, type=Path)
    base = ap.parse_args().run_directory.resolve(); base.relative_to(p.REPO/'_scratch/s11c'); base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900)
    started = time.monotonic(); io.digest = saved.digest
    reader, journal, packets = saved.Reader(), inputs.MetadataJournal(base), {}
    def packet(rec):
        path = rec.get('path', rec.get('logical'))
        if path not in packets: packets[path] = reader.packet(path, rec['sha256'])
        return packets[path]
    def address(route):
        result = packet(route['packet'])
        for key in route['keys']: result = result[key]
        return result
    ready = reader.json(READY/'checks.json', READY_SHA)
    require(ready['status'] == 'COMPLETED_SAVED_REMAINDER_INPUT_INSPECTION' and ready['newScientificCalls'] == 0, 'completed saved operand inspection')
    def metadata(name):
        rec = ready['artifacts'][name]; return reader.json(rec['path'], rec['sha256'])
    caller = metadata('native-caller-and-unsaved-internal-history.json')
    for entry in caller['native'].values():
        rec = reader.retain(entry['file']['logical'], entry['file']['sha256']); tree = ast.parse(Path(rec['canonical']).read_text())
        for name, body in entry['bodies'].items():
            require(ast.dump(next(n for n in tree.body if getattr(n, 'name', None) == name)) == ast.dump(ast.parse(body).body[0]), 'whole native caller identity')
    for rec in (caller['complexCaller'], caller['nativeHeavisidePrinter']['file'], caller['nativeHeavisidePrinter']['namespaceFile']):
        reader.retain(rec['logical'], rec['sha256'])
    require(caller['oldPrepareBasisRan'] and caller['nativeHeavisidePrinter']['actualTranslation'] == 'heaviside', 'completed old native basis and actual Heaviside convention')
    views = [metadata('rows/'+str(i)+'.json') for i in (51, 52, 53, 54)]
    selected = views[0]['route']; own = packet(selected['sourceInputs']['packet']); scalars = packet(selected['scalarInputs']); system = packet(selected['basis'])
    physical = reader.json(selected['ownPhysicalRoutes']['packet']['logical'], selected['ownPhysicalRoutes']['packet']['sha256'])
    common = packet(physical['context']); context = common['contextPair'][0]
    require(same(*common['contextPair']) and same(*common['basisPair']), 'whole physical context and field basis')
    for rec in physical.values():
        if isinstance(rec, dict) and 'logical' in rec: reader.retain(rec['logical'], rec['sha256'])
    rule_meta = views[0]['savedRulePacket']; prepared = packet(rule_meta['packet']); native_rule = prepared['original']
    nodes, weights = native_rule['source_nodes'], native_rule['source_weights']; size = native_rule['size']; settings = own['settings']
    node_route = {'packet': rule_meta['packet'], 'keys': ['original', 'source_nodes']}
    weight_route = {'packet': rule_meta['packet'], 'keys': ['original', 'source_weights']}
    require(same(settings, context['settings']) and same(settings, system['settings']) and
            same(rule_meta['sourceSettings'], json.loads(json.dumps(settings))) and size == len(system['nodes']) == 129 and
            nodes.shape == weights.shape == (1024,) and np.all(weights > 0), 'actual full source rule settings and trial size')
    require(complex(scalars['frequency']) == 1-.01j and settings['regulator'] == .1 and settings['sourceBound'] == 64 and settings['profileBound'] == 14, 'fixed actual practical setting')
    reader.retain(Path(p.__file__), 'c579b8881e890f0257fb7cebb5b50ca30dce356b3b2141f622f0b35deb7ad2ec')
    reader.retain(Path(saved.recovery.__file__), '7d9d446165f6e3301dabd5dc225f1585930f5879a1feb65fa58e96ba24f82a40')
    evaluate, evaluation_join = evaluator()
    rule_file = Path(inspect.getsourcefile(special.roots_legendre)); rule_module = ast.parse(rule_file.read_text())
    new_rule_bodies = {n.name: ast.unparse(n) for n in rule_module.body if isinstance(n, ast.FunctionDef) and n.name in ('roots_legendre', '_gen_roots_and_weights')}
    require(len(new_rule_bodies) == 2, 'whole new SciPy rule implementation')
    for path in (Path(__file__).resolve(), M/(NAME+'_plan.md'), Path(inputs.__file__), Path(saved.__file__),
                 Path(p.__file__), Path(saved.recovery.__file__), rule_file, M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):
        reader.retain(path)
    for path in (Path(__file__).resolve(), M/(NAME+'_plan.md')):
        target = base/'source'/path.name; target.parent.mkdir(exist_ok=True)
        with target.open('xb') as out: out.write(path.read_bytes())
        reader.retain(target, saved.digest(path))
    journal.json('native-and-new-numerical-callers.json', {'savedCaller': ready['artifacts']['native-caller-and-unsaved-internal-history.json'],
        'evaluator': evaluation_join, 'newSourceRecurrence': inspect.getsource(p.source_matrix),
        'selectedIndependentSource': inspect.getsource(p.independent_source), 'newRuleFile': reader.retain(rule_file),
        'newRuleBodies': new_rule_bodies, 'nativePrepareBasisCalled': False,
        'interpretation': 'New numerical recurrence and split SciPy384/768-node profile rules; no old native arrays, lambdas, derivatives or Gauss returns restored.'})
    scope = {'case': p.CASE, 'rows': [51, 52, 53, 54], 'sources': [5, 7, 9, 11], 'frequency': {'real': 1., 'imag': -.01},
        'sourceRule': {'nodes': node_route, 'weights': weight_route, 'size': size, 'bound': 64.},
        'sourceMethod': 'Numerical expression tree + simultaneous weighted Chebyshev recurrence, full own coefficients; selected source5/order1 and source11/order2 trigonometric comparisons.',
        'profileMethod': 'New scipy.special.roots_legendre384/768 on each of[-14,0],[0,14]; retain finite remainder including Heaviside.',
        'profileDifferences': {'first': -8., 'last': 8., 'count': 4097, 'spacing': 1/256},
        'profileComparisons': [-8., -95/256, 0., 95/256, 8.], 'profileAbsoluteTarget': 1e-8,
        'sourceScaledTarget': 2e-10, 'laterCallGate': '3x actual first-rule cost +60 inside remaining900s',
        'scope': 'Source/profile preparation only. No row/scattering result, Abel limit, rigorous error bound, physical pole or domain acceptance. Positive regulator term remains in full row input.'}
    journal.json('preparation-scope.json', scope)
    def forbidden(*args, **kwargs): raise RuntimeError('completed native science disabled in new numerical preparation')
    io.native.f.source_jets = io.native.f.polynomial_basis = io.native.f.BasisMomentum.prepare_basis = forbidden
    io.native.Pair.__init__ = io.native.maps = io.native.continue_pair = forbidden
    for name in ('diff', 'lambdify', 'cancel', 'expand', 'factor', 'solve', 'gcd', 'resultant', 'integrate'): setattr(sp, name, forbidden)
    # Previously saved new numerical source22 calls, not its unsaved internals.
    pilot_cp = reader.json(M/'S11c_d_remaining_case_frequency_row_pilot_checkpoint.json',
                           '0c2448600fbd33b82b5e9c1a0d6b16ddbf5dbb423decff540ceec7d221098e15')
    pc = reader.json(Path(pilot_cp['runDirectory'])/'checks.json', pilot_cp['checksSha256'])
    arts = pc['artifacts']; prior_input = packet(arts['source-action/input.pickle']); available = []
    for name, rec in arts.items():
        if name.startswith('source-action/coefficients/') and name.endswith('/input.pickle'):
            folder = name.rsplit('/', 1)[0]; arg = packet(rec); val = arts[folder+'/value.pickle']
            available.append({'argument': arg, 'input': rec, 'value': val, 'ownerInput': arts['source-action/input.pickle'],
                              'units': (prior_input['jet']['amplitudeUnit'], prior_input['jet']['integralUnit'])})
    coefficients_cache = []; outputs = []; new_coefficients = 0
    for view in views:
        route = view['route']; ri, si = view['rowIndex'], view['sourceIndex']; row = own['bound']['rows'][ri]
        source = own['bound']['sources'][0, si]; jet = address(route['jets'][str(si)])
        require(row['index'] == ri and len(row['factors']) == 1 and row['factors'][0]['sourceIndex'] == si and
                same(jet['originalBoundAmplitude'], scalars['actual']['source', si]) and
                same(jet['amplitudeUnit'], source['amplitudeUnit']) and same(jet['integralUnit'], source['integralUnit']) and
                jet['probe'].args == (context['zp'],) and all(v.free_symbols <= {context['zp']} for v in jet['coefficients']), 'full own source coefficients and variable/units')
        require(not view['savedBasisSearch']['savedCoefficientBasisCandidates'] and route['scalarInputs'] == selected['scalarInputs'] and
                route['sourceInputs']['packet'] == selected['sourceInputs']['packet'] and view['physicalRoutes'] == physical, 'whole actual saved inputs and existing source-array search')
        folder = 'sources/'+str(si)
        full_input = journal.write(folder+'/input.pickle', {'row': row, 'source': source, 'jet': jet, 'context': context,
            'coefficient': scalars['actual']['factor', ri, 0], 'settings': settings, 'fieldUnits': own['fieldUnits'],
            'equationUnits': own['equationUnits'], 'profileUnits': own['bound']['profileUnits'], 'abel': own['bound']['abel'],
            'pairs': own['bound']['pairs'], 'physicalRoutes': physical, 'nodeRoute': node_route, 'weightRoute': weight_route,
            'bound': settings['sourceBound'], 'size': size, 'sourceInputRoute': route})
        values, routes = [], []
        units = (jet['amplitudeUnit'], jet['integralUnit'])
        for order, expression in enumerate(jet['coefficients']):
            prefix = folder+'/coefficients/'+str(order)
            wanted = {'expression': expression, 'variable': context['zp'], 'nodes': node_route, 'units': units}
            arg = journal.write(prefix+'/input.pickle', {'requested': wanted, 'ownInput': full_input, 'order': order})
            matches = [v for v in coefficients_cache if same(v['requested'], wanted)]
            old_matches = [v for v in available if same(v['argument']['expression'], expression) and
                same(v['argument']['variable'], context['zp']) and v['argument']['nodes'] == node_route and
                same(v['units'], units) and same(prior_input['context'], context)]
            if matches:
                value, vr, disposition = matches[0]['value'], matches[0]['route'], 'COMPLETED_NEW_VALUE'
            elif old_matches:
                value, vr, disposition = packet(old_matches[0]['value']), old_matches[0]['value'], 'ACCEPTED_VALUE'
                require(all(same(value, packet(v['value'])) for v in old_matches), 'full matching saved coefficient returns agree')
            else:
                with np.errstate(over='raise', invalid='raise', divide='raise', under='ignore'):
                    value = np.broadcast_to(np.asarray(evaluate(expression, {context['zp']: nodes}), complex), nodes.shape)
                vr = journal.write(prefix+'/value.pickle', value); disposition = 'NEW'; new_coefficients += 1
            journal.json(prefix+'/completed.json', {'input': arg, 'value': vr, 'disposition': disposition,
                'oldMatches': [v['input'] for v in old_matches], 'precedingNewMatches': [v['input'] for v in matches]})
            require(value.shape == nodes.shape and np.isfinite(value).all(), 'finite full new coefficient array')
            coefficients_cache.append({'requested': wanted, 'input': arg, 'value': value, 'route': vr}); values.append(value); routes.append(vr)
        ai = journal.write(folder+'/action-input.pickle', {'coefficients': routes, 'nodes': node_route, 'weights': weight_route,
            'bound': settings['sourceBound'], 'size': size, 'units': units, 'ownInput': full_input,
            'method': 'new simultaneous value/derivative numerical recurrence'})
        matrix = p.source_matrix(values, nodes, weights, settings['sourceBound'], size)
        ar = journal.write(folder+'/action-value.pickle', matrix); journal.json(folder+'/action-completed.json', {'input': ai, 'value': ar})
        require(matrix.shape == (1024, 129) and matrix.dtype == np.dtype(complex) and np.isfinite(matrix).all(), 'full new numerical source action')
        comparison = None
        if si in (5, 11):
            ci = journal.write(folder+'/comparison-input.pickle', {'coefficientValues': routes, 'actionInput': ai,
                'actionValue': ar, 'method': 'independent closed trigonometric Chebyshev derivatives'})
            other = p.independent_source(values, nodes, weights, settings['sourceBound'], size)
            cv = journal.write(folder+'/comparison-value.pickle', other)
            spread = float(np.max(abs(matrix-other))/(1+np.max(abs(matrix))))
            comparison = {'input': ci, 'value': cv, 'scaledDifference': spread, 'target': 2e-10}
            journal.json(folder+'/comparison.json', comparison); require(spread < 2e-10, 'selected full source recurrence comparison')
        summary = {'rowIndex': ri, 'sourceIndex': si, 'ownInput': full_input, 'jet': route['jets'][str(si)],
            'coefficientValues': routes, 'actionInput': ai, 'actionValue': ar, 'comparison': comparison}
        journal.json(folder+'/summary.json', summary); outputs.append(summary)
    # Full original coefficient and its finite Integral stay saved, including Abel.
    coefficient = scalars['actual']['factor', 51, 0]; profiles = tuple(coefficient.atoms(sp.Integral))
    require(len(profiles) == 1, 'one actual finite profile remainder'); integral = profiles[0]
    k, q = (v[0] for v in own['bound']['rows'][51]['limits']); xi = context['xi']
    require(len(integral.limits) == 1 and same(integral.limits[0][0], xi) and
            tuple(map(float, integral.limits[0][1:])) == (-14., 14.) and integral.function.func is sp.Mul, 'literal finite profile and product')
    envelope = [v for v in integral.function.args if v.free_symbols <= {xi}]
    phase = [v for v in integral.function.args if v not in envelope]
    require(len(envelope) == len(phase) == 1 and phase[0].func is sp.exp, 'literal envelope and Fourier phase')
    choices = []
    for node in sp.preorder_traversal(phase[0]):
        if isinstance(node, sp.Add) and len(node.args) == 2 and k in node.args:
            rest = next(a for a in node.args if a != k)
            if isinstance(rest, sp.Mul) and set(rest.args) == {sp.S.NegativeOne, q}: choices.append(node)
    require(len(set(choices)) == 1, 'literal k-q phase subtree'); difference = choices[0]
    def covered(node):
        if node == difference: return
        require(node not in (k, q), 'no momentum hidden outside difference')
        for child in node.args: covered(child)
    covered(phase[0])
    for row in (52, 53, 54): require(same(coefficient, scalars['actual']['factor', row, 0]), 'actual shared coefficient expression retained; units remain own')
    pi = journal.write('profile/input.pickle', {'integral': integral, 'envelope': envelope[0], 'phase': phase[0],
        'differenceSubtree': difference, 'profileUnit': own['bound']['profileUnits'][integral],
        'context': context, 'settings': settings, 'ownRows': [v['ownInput'] for v in outputs], 'originalCoefficient': coefficient})
    difference_grid = np.arange(-2048, 2049, dtype=float)/256
    dg = journal.write('profile/difference-grid.pickle', difference_grid)
    selected_indices = np.array([0, 1953, 2048, 2143, 4096]); rule_results = {}; costs = []
    for order in (384, 768):
        if costs:
            remaining = 900-(time.monotonic()-started); reserve = 3*costs[0]+60
            journal.json('profile/next-rule-cost.json', {'firstSeconds': costs[0], 'remainingSeconds': remaining, 'requiredSeconds': reserve})
            require(reserve < remaining, 'bounded next profile rule fits')
        tick = time.monotonic(); folder = 'profile/rules/'+str(order)
        ri = journal.write(folder+'/rule-input.pickle', {'function': 'scipy.special.roots_legendre', 'args': (order,), 'kwargs': {},
            'source': reader.retain(rule_file), 'newAlgorithm': True})
        standard_nodes, standard_weights = special.roots_legendre(order)
        rv = journal.write(folder+'/rule-value.pickle', (standard_nodes, standard_weights)); journal.json(folder+'/rule-completed.json', {'input': ri, 'value': rv})
        si = journal.write(folder+'/split-input.pickle', {'standardRule': rv, 'intervals': ((-14., 0.), (0., 14.)), 'nativeProfile': pi})
        split_nodes = np.concatenate((7*standard_nodes-7, 7*standard_nodes+7)); split_weights = np.tile(7*standard_weights, 2)
        sv = journal.write(folder+'/split-value.pickle', {'nodes': split_nodes, 'weights': split_weights}); journal.json(folder+'/split-completed.json', {'input': si, 'value': sv})
        ei = journal.write(folder+'/envelope-input.pickle', {'expression': envelope[0], 'variable': xi, 'rule': sv, 'profile': pi})
        envelope_values = np.asarray(evaluate(envelope[0], {xi: split_nodes}), complex)
        ev = journal.write(folder+'/envelope-value.pickle', envelope_values); journal.json(folder+'/envelope-completed.json', {'input': ei, 'value': ev})
        wanted = difference_grid[selected_indices] if order == 384 else difference_grid
        values, routes = [], []
        for start in range(0, len(wanted), 64):
            coords = wanted[start:start+64]; prefix = folder+'/batches/'+str(start)
            ar = journal.write(prefix+'/input.pickle', {'profile': pi, 'differences': coords, 'rule': sv, 'envelope': ev, 'phase': phase[0]})
            phase_values = np.asarray(evaluate(phase[0], {xi: split_nodes[None, :], difference: coords[:, None]}), complex)
            pr = journal.write(prefix+'/phase-value.pickle', phase_values)
            integrand = phase_values*envelope_values[None, :]
            ir = journal.write(prefix+'/integrand-value.pickle', integrand)
            cr = journal.write(prefix+'/contraction-input.pickle', {'integrand': ir, 'rule': sv, 'profile': pi, 'operation': 'integrand @ actual split weights'})
            value = integrand@split_weights
            vr = journal.write(prefix+'/value.pickle', value)
            journal.json(prefix+'/completed.json', {'input': ar, 'phase': pr, 'integrand': ir, 'contraction': cr, 'value': vr})
            require(np.isfinite(value).all(), 'finite actual remainder profile value')
            values.extend(value); routes.extend({'packet': vr, 'keys': [i]} for i in range(len(coords)))
        rule_results[order] = np.asarray(values)
        cost = time.monotonic()-tick; costs.append(cost)
        journal.json(folder+'/summary.json', {'orderPerLeg': order, 'profile': pi, 'differenceGrid': dg,
            'selectedIndices': selected_indices.tolist() if order == 384 else None, 'values': routes, 'wallSeconds': cost})
    ci = journal.write('profile/comparison-input.pickle', {'first': journal.artifacts['profile/rules/384/summary.json'],
        'second': journal.artifacts['profile/rules/768/summary.json'], 'selectedIndices': selected_indices})
    delta = rule_results[768][selected_indices]-rule_results[384]
    cv = journal.write('profile/comparison-value.pickle', delta); absolute = float(np.max(abs(delta)))
    comparison = {'input': ci, 'value': cv, 'maximumAbsoluteDifference': absolute, 'target': 1e-8, 'targetMet': absolute <= 1e-8}
    journal.json('profile/comparison.json', comparison)
    reader.postcheck()
    for rec in journal.artifacts.values(): require(saved.digest(rec['path']) == rec['sha256'], 'all new bytes unchanged')
    journal.json('inputs.json', {'savedInspection': reader.retain(READY/'checks.json', READY_SHA), 'sourceNodes': node_route,
        'sourceWeights': weight_route, 'consumedRoutes': reader.routes})
    checks = {'status': 'COMPLETED_NEW_REMAINDER_NUMERICAL_PREPARATION', 'sourceActions': outputs,
        'newSourceActions': len(outputs), 'newCoefficientArrays': new_coefficients, 'profile': comparison,
        'profileValueCounts': {str(k): len(v) for k, v in rule_results.items()}, 'ruleSeconds': costs,
        'newNativeBasisOrCompilerCalls': 0, 'newRows': 0, 'allConsumedHashesUnchanged': True,
        'artifacts': dict(journal.artifacts), 'wallSeconds': time.monotonic()-started, 'scope': scope['scope']}
    journal.json('checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__': main()
