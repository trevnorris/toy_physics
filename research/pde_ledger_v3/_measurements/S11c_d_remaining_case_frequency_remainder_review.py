#!/usr/bin/env python3
"""Bounded saved source/profile review, with no numerical reconstruction."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time

import S11c_d_remaining_case_frequency_remainder_prepare as producer

saved, p, sp, np = producer.saved, producer.p, producer.sp, producer.np
M, F, require, same = producer.M, producer.F, producer.require, producer.same
ROOT = F/'remainder-preparation'
CHECKS_SHA = '352e03a883166253634761279353d7d32a64fb2a9991158663bf995aeac43950'
HELPER_SHA = '685b35715c212ea8cc5a050db276bd849de8c992c3285647d583fac5c201ee9e'
NAME = 'S11c_d_remaining_case_frequency_remainder_review'


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory', type=Path, required=True)
    base = ap.parse_args().run_directory.resolve(); base.relative_to(p.REPO/'_scratch/s11c'); base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900)
    started = time.monotonic(); saved.io.digest = saved.digest
    reader, journal = saved.Reader(), producer.inputs.MetadataJournal(base)
    checks = reader.json(ROOT/'complete/checks.json', CHECKS_SHA); arts = checks['artifacts']
    require(checks['status'] == 'COMPLETED_NEW_REMAINDER_NUMERICAL_PREPARATION' and checks['newRows'] == 0, 'completed preparation only')
    guard = reader.json(ROOT/'resource-guard/outcome.json'); child = reader.json(ROOT/'resource-guard/child-outcome.json')
    invocation = reader.json(ROOT/'frequency_remainder_prepare.invocation.json'); limits = reader.json(ROOT/'resource-guard/effective-limits.json')
    require(guard['exitCode'] == child['exitCode'] == invocation['exitCode'] == 0 and child['guardReason'] is None and guard['limitsVerified'], 'actual final producer exits')
    require(limits['memory.max'] == '2147483648' and limits['memory.swap.max'] == '0' and limits['pids.max'] == '32' and limits['nice'] == 15 and len(limits['affinity']) == 1 and set(limits['threads'].values()) == {'1'}, 'mandatory whole-job limits')
    for name in ('frequency_remainder_prepare.stderr', 'guard.stderr', 'resource-guard/stderr'):
        rec = reader.retain(ROOT/name); require(rec['bytes'] == 0, 'empty producer strict stderr')
    require(reader.retain(ROOT/'frequency_remainder_prepare.stdout')['sha256'] == CHECKS_SHA, 'producer checks/stdout byte identity')
    samples = ROOT/'resource-guard/resource-samples.jsonl'; reader.retain(samples)
    for line in samples.read_text().splitlines():
        s = json.loads(line); events = dict(x.split() for x in s['memory.events'].splitlines())
        require(int(s['memory.swap.current']) == 0 and int(s['memory.peak']) <= 2*1024**3 and all(int(events[k]) == 0 for k in ('max', 'oom', 'oom_kill')), 'actual producer cap/OOM/swap telemetry')
    reader.retain(Path(producer.__file__), HELPER_SHA)
    def meta(name): return reader.json(arts[name]['path'], arts[name]['sha256'])
    def packet(rec): return reader.packet(rec.get('path', rec.get('logical')), rec['sha256'])
    def value(name): return packet(arts[name])
    def address(route):
        v = packet(route['packet'])
        for key in route['keys']: v = v[key]
        return v
    def array(v, shape, dtype=complex):
        require(isinstance(v, np.ndarray) and v.shape == shape and v.dtype == np.dtype(dtype) and np.isfinite(v).all(), 'saved complete finite array/type/shape')
    inputs = meta('inputs.json')
    for logical, rec in inputs['consumedRoutes'].items(): reader.retain(logical, rec['sha256'])
    prepared = packet(inputs['sourceNodes']['packet'])
    nodes, weights = prepared, prepared
    for key in inputs['sourceNodes']['keys']: nodes = nodes[key]
    for key in inputs['sourceWeights']['keys']: weights = weights[key]
    array(nodes, (1024,), float); array(weights, (1024,), float)
    require(np.all(weights > 0), 'complete accepted source nodes and positive weights')
    del prepared
    caller = meta('native-and-new-numerical-callers.json'); rule_file = caller['newRuleFile']
    tree = ast.parse(Path(rule_file['canonical']).read_text())
    for name, body in caller['newRuleBodies'].items(): require(ast.dump(next(n for n in tree.body if getattr(n, 'name', None) == name)) == ast.dump(ast.parse(body).body[0]), 'whole new rule source body')
    source_module = ast.parse(Path(p.__file__).read_text())
    for name, key in (('source_matrix', 'newSourceRecurrence'), ('independent_source', 'selectedIndependentSource')):
        require(ast.dump(next(n for n in source_module.body if getattr(n, 'name', None) == name)) == ast.dump(ast.parse(caller[key]).body[0]), 'whole source numerical algorithm body')
    scope = meta('preparation-scope.json'); require(scope['frequency'] == {'real': 1., 'imag': -.01} and scope['sourceRule']['size'] == 129, 'scoped actual fixed input')
    def forbidden(*args, **kwargs): raise RuntimeError('saved review disables scientific calculation')
    producer.main = producer.evaluator = producer.special.roots_legendre = forbidden
    p.source_matrix = p.independent_source = forbidden
    saved.io.native.f.source_jets = saved.io.native.f.polynomial_basis = saved.io.native.f.BasisMomentum.prepare_basis = forbidden
    saved.io.native.Pair.__init__ = saved.io.native.maps = saved.io.native.continue_pair = forbidden
    journal.write = forbidden
    for name in ('diff', 'lambdify', 'cancel', 'expand', 'factor', 'solve', 'gcd', 'resultant', 'integrate'): setattr(sp, name, forbidden)
    sources = []; coefficient_count = 0; actual_context = None
    for summary in checks['sourceActions']:
        si, ri = summary['sourceIndex'], summary['rowIndex']; folder = 'sources/'+str(si); own = packet(summary['ownInput'])
        require(own['row']['index'] == ri and own['row']['factors'][0]['sourceIndex'] == si and
                same(own['jet'], address(summary['jet'])) and same(own['source']['integralUnit'], own['jet']['integralUnit']) and
                same(own['source']['amplitudeUnit'], own['jet']['amplitudeUnit']) and own['nodeRoute'] == inputs['sourceNodes'] and own['weightRoute'] == inputs['sourceWeights'], 'actual own source, units and saved-rule routes')
        selected = own['sourceInputRoute']; raw = packet(selected['sourceInputs']['packet']); scalars = packet(selected['scalarInputs'])
        require(same(raw['bound']['rows'][ri], own['row']) and same(raw['bound']['sources'][0, si], own['source']) and
                same(scalars['actual']['factor', ri, 0], own['coefficient']) and same(scalars['actual']['source', si], own['jet']['originalBoundAmplitude']) and
                same(raw['settings'], own['settings']) and same(raw['fieldUnits'], own['fieldUnits']) and same(raw['equationUnits'], own['equationUnits']) and
                same(raw['bound']['profileUnits'], own['profileUnits']) and same(raw['bound']['abel'], own['abel']) and same(raw['bound']['pairs'], own['pairs']), 'entire own physical input fields and scalar returns')
        common = packet(own['physicalRoutes']['context'])
        require(same(*common['contextPair']) and same(*common['basisPair']) and same(own['context'], common['contextPair'][0]), 'own context and field/equation basis')
        actual_context = own['context'] if actual_context is None else actual_context
        require(same(own['context'], actual_context), 'shared physical context while preserving own source')
        for order, expression in enumerate(own['jet']['coefficients']):
            prefix = folder+'/coefficients/'+str(order); arg = value(prefix+'/input.pickle'); receipt = meta(prefix+'/completed.json'); request = arg['requested']
            require(same(request['expression'], expression) and same(request['variable'], own['context']['zp']) and
                    request['nodes'] == own['nodeRoute'] and same(request['units'], (own['jet']['amplitudeUnit'], own['jet']['integralUnit'])) and
                    arg['ownInput'] == summary['ownInput'] and arg['order'] == order and receipt['input'] == arts[prefix+'/input.pickle'] and receipt['value'] == summary['coefficientValues'][order], 'whole actual numerical coefficient input/owner/unit/receipt')
            require(receipt['disposition'] == 'NEW' and not receipt['oldMatches'] and not receipt['precedingNewMatches'], 'actual captured new-call branch; no cache inferred from restored keys')
            array(packet(receipt['value']), (1024,)); coefficient_count += 1
        ai = packet(summary['actionInput']); ar = meta(folder+'/action-completed.json'); matrix = packet(summary['actionValue'])
        require(ai['coefficients'] == summary['coefficientValues'] and ai['nodes'] == own['nodeRoute'] and ai['weights'] == own['weightRoute'] and
                ai['bound'] == own['bound'] == 64 and ai['size'] == own['size'] == 129 and ai['ownInput'] == summary['ownInput'] and
                same(ai['units'], (own['jet']['amplitudeUnit'], own['jet']['integralUnit'])) and ar == {'input': summary['actionInput'], 'value': summary['actionValue']}, 'full source recurrence input/output receipt')
        array(matrix, (1024, 129))
        comparison = summary['comparison']
        if comparison is not None:
            require(si in (5, 11) and comparison == meta(folder+'/comparison.json') and 0 <= comparison['scaledDifference'] < 2e-10, 'saved selected source comparison')
            ci = packet(comparison['input']); require(ci['coefficientValues'] == summary['coefficientValues'] and ci['actionInput'] == summary['actionInput'] and ci['actionValue'] == summary['actionValue'], 'actual independent comparison arguments')
            array(packet(comparison['value']), (1024, 129))
        else: require(si in (7, 9), 'declared selected source comparisons')
        require(meta(folder+'/summary.json') == summary, 'complete source summary')
        journal.json('sources/'+str(si)+'.json', {'source': summary, 'fullTypedInputAndUnitsJoined': True,
            'sourceArrayFinite': True, 'newScienceInReview': 0, 'comparisonNotRecomputed': True})
        sources.append(summary)
    profile = value('profile/input.pickle'); grid = value('profile/difference-grid.pickle'); array(grid, (4097,), float)
    require(profile['ownRows'] == [v['ownInput'] for v in sources] and same(profile['context'], actual_context), 'own profile physical input routes')
    own = packet(sources[0]['ownInput']); require(same(profile['originalCoefficient'], own['coefficient']) and
        profile['integral'] in own['profileUnits'] and same(profile['profileUnit'], own['profileUnits'][profile['integral']]), 'actual finite profile in unmodified coefficient including Abel')
    require(len(profile['integral'].limits) == 1 and tuple(map(float, profile['integral'].limits[0][1:])) == (-14., 14.) and
            profile['envelope'] in profile['integral'].function.args and profile['phase'] in profile['integral'].function.args, 'literal finite Integral operand tree')
    phase_count = 0; profile_count = 0
    for order, count in ((384, 5), (768, 4097)):
        folder = 'profile/rules/'+str(order); arg = value(folder+'/rule-input.pickle'); receipt = meta(folder+'/rule-completed.json')
        require(arg['function'] == 'scipy.special.roots_legendre' and arg['args'] == (order,) and arg['kwargs'] == {} and arg['source'] == rule_file and
                receipt == {'input': arts[folder+'/rule-input.pickle'], 'value': arts[folder+'/rule-value.pickle']}, 'actual new alternative rule call and source')
        standard = value(folder+'/rule-value.pickle'); array(standard[0], (order,), float); array(standard[1], (order,), float)
        split_arg = value(folder+'/split-input.pickle'); split = value(folder+'/split-value.pickle')
        require(split_arg['standardRule'] == receipt['value'] and split_arg['intervals'] == ((-14., 0.), (0., 14.)) and split_arg['nativeProfile'] == arts['profile/input.pickle'] and
                meta(folder+'/split-completed.json') == {'input': arts[folder+'/split-input.pickle'], 'value': arts[folder+'/split-value.pickle']}, 'actual finite split rule input/receipt')
        array(split['nodes'], (2*order,), float); array(split['weights'], (2*order,), float)
        require(np.all(split['weights'] > 0), 'positive saved finite profile weights')
        envarg = value(folder+'/envelope-input.pickle'); envelope = value(folder+'/envelope-value.pickle')
        require(same(envarg['expression'], profile['envelope']) and same(envarg['variable'], profile['integral'].limits[0][0]) and envarg['rule'] == arts[folder+'/split-value.pickle'] and envarg['profile'] == arts['profile/input.pickle'], 'actual envelope inputs')
        array(envelope, (2*order,)); require(meta(folder+'/envelope-completed.json') == {'input': arts[folder+'/envelope-input.pickle'], 'value': arts[folder+'/envelope-value.pickle']}, 'envelope receipt')
        summary = meta(folder+'/summary.json'); wanted = grid[summary['selectedIndices']] if order == 384 else grid
        require(len(summary['values']) == count and summary['profile'] == arts['profile/input.pickle'] and summary['differenceGrid'] == arts['profile/difference-grid.pickle'], 'full profile return catalogue')
        for start in range(0, count, 64):
            prefix = folder+'/batches/'+str(start); arg = value(prefix+'/input.pickle'); receipt = meta(prefix+'/completed.json'); length = len(arg['differences'])
            require(same(arg['differences'], wanted[start:start+64]) and arg['rule'] == arts[folder+'/split-value.pickle'] and arg['envelope'] == arts[folder+'/envelope-value.pickle'] and arg['profile'] == arts['profile/input.pickle'] and same(arg['phase'], profile['phase']), 'actual profile phase arguments')
            array(packet(receipt['phase']), (length, 2*order)); array(packet(receipt['integrand']), (length, 2*order)); array(packet(receipt['value']), (length,))
            require(receipt == {'input': arts[prefix+'/input.pickle'], 'phase': arts[prefix+'/phase-value.pickle'],
                'integrand': arts[prefix+'/integrand-value.pickle'], 'contraction': arts[prefix+'/contraction-input.pickle'], 'value': arts[prefix+'/value.pickle']}, 'actual profile batch receipts')
            contraction = packet(receipt['contraction']); require(contraction['integrand'] == receipt['integrand'] and contraction['rule'] == arg['rule'] and contraction['profile'] == arg['profile'], 'saved full contraction arguments')
            require(summary['values'][start:start+length] == [{'packet': receipt['value'], 'keys': [i]} for i in range(length)], 'exact saved output element routes')
            journal.json('profile/'+str(order)+'/batch-'+str(start)+'.json', {'receipt': receipt, 'count': length, 'nativeInputSourceAndShapesJoined': True, 'scientificRecomputation': False})
            phase_count += 1; profile_count += length
        journal.json('profile/'+str(order)+'/summary.json', {'saved': arts[folder+'/summary.json'], 'values': count, 'newScienceInReview': 0})
    comparison = meta('profile/comparison.json'); ci = packet(comparison['input']); delta = packet(comparison['value']); array(delta, (5,))
    require(comparison == checks['profile'] and comparison['targetMet'] and 0 <= comparison['maximumAbsoluteDifference'] <= comparison['target'] == 1e-8 and
            ci['first'] == arts['profile/rules/384/summary.json'] and ci['second'] == arts['profile/rules/768/summary.json'] and
            ci['selectedIndices'].tolist() == [0, 1953, 2048, 2143, 4096], 'actual saved selected profile comparison result')
    require(coefficient_count == checks['newCoefficientArrays'] == 9 and len(sources) == checks['newSourceActions'] == 4 and profile_count == 4102, 'actual saved receipt counts')
    journal.json('validated-preparation.json', {'sources': sources, 'profileComparison': comparison,
        'profileInput': arts['profile/input.pickle'], 'profileFineCatalogue': arts['profile/rules/768/summary.json'],
        'sourceCalls': len(sources), 'coefficientCalls': coefficient_count, 'profileBatches': phase_count,
        'nativeBasisReconstructed': False, 'newScientificCalls': 0, 'scope': checks['scope']})
    for name, rec in arts.items(): reader.retain(rec['path'], rec['sha256'])
    for path in (Path(__file__).resolve(), M/(NAME+'_plan.md')):
        reader.retain(path); dest = base/'source'/path.name; dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out: out.write(path.read_bytes())
        reader.retain(dest, saved.digest(path))
    reader.postcheck(); journal.json('validated-paths.json', reader.routes)
    result = {'status': 'PASSED_BOUNDED_SAVED_REMAINDER_PREPARATION_REVIEW', 'sourceActions': 4, 'coefficientCalls': 9,
        'profileBatches': phase_count, 'profileValues': profile_count, 'newScientificCalls': 0,
        'consumedLogicalPaths': len(reader.routes), 'allConsumedHashesUnchanged': True,
        'artifacts': dict(journal.artifacts), 'wallSeconds': time.monotonic()-started}
    journal.json('checks.json', result); signal.alarm(0); print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
