#!/usr/bin/env python3
"""Read saved source and mixed-profile operands for four unfinished rows."""
import argparse
import ast
import inspect
import json
from pathlib import Path
import resource
import signal
import time

import S11c_d_remaining_case_frequency_row_2d_pilot as saved

p, sp, io = saved.p, saved.sp, saved.io
M, F, require, same = saved.M, saved.F, saved.require, saved.same
NAME = 'S11c_d_remaining_case_frequency_remainder_inputs'
CP = M/'S11c_d_remaining_case_frequency_row_47_checkpoint.json'
CP_SHA = '2d2c98959f64f2edf48c10decd7f25051f8ac660bcd96ff4855966987eb7ad18'


def tree(value):
    return {'type': type(value).__name__, 'expression': str(value),
            'arguments': [tree(child) for child in value.args]}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory', type=Path, required=True)
    base = ap.parse_args().run_directory.resolve(); base.relative_to(p.REPO/'_scratch/s11c')
    base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900)
    started = time.monotonic(); io.digest = saved.digest
    reader, journal, packets = saved.Reader(), saved.Journal(base), {}
    def packet(rec):
        path = rec.get('path', rec.get('logical'))
        if path not in packets: packets[path] = reader.packet(path, rec['sha256'])
        return packets[path]
    def address(route):
        result = packet(route['packet'])
        for key in route['keys']: result = result[key]
        return result
    cp = reader.json(CP, CP_SHA)
    require(cp['status'] == 'ACCEPTED_BOUNDED_CASE_FREQUENCY_ROW_47', 'latest accepted paired row')
    inventory = reader.json(p.READY/'completed-input-artifact-inventory.json', p.INVENTORY_SHA)
    def metadata(name): return reader.json(p.READY/'complete'/name, inventory[name]['sha256'])
    reader.retain(p.READY/'complete/checks.json', p.READY_SHA)
    native = metadata('native-callers.json')
    for item in native.values():
        record = reader.retain(item['file']['logical'], item['file']['sha256'])
        module = ast.parse(Path(record['canonical']).read_text())
        for name, body in item['bodies'].items():
            node = next(n for n in module.body if getattr(n, 'name', None) == name)
            require(ast.dump(node) == ast.dump(ast.parse(body).body[0]), 'whole native caller input')
    from sympy.printing.numpy import NumPyPrinter
    printer_file = Path(inspect.getsourcefile(NumPyPrinter))
    caller = M/'S11c_d_frequency_matrix.py'
    reader.retain(caller, '97797ab676a58f16785bf4234d41eb5cd01abc0a42f39edc4e20f926b2c2b336')
    module = ast.parse(caller.read_text())
    complex_body = next(n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == 'complex_case')
    require(any(isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == 'prepare_basis' for n in ast.walk(complex_body)), 'old complex prepare_basis actually called')
    journal.json('native-caller-and-unsaved-internal-history.json', {
        'native': native, 'complexCaller': reader.retain(caller), 'complexBody': ast.unparse(complex_body),
        'oldPrepareBasisRan': True, 'missingSerializedActionDoesNotLicenseReplay': True,
        'nativeHeavisidePrinter': {'file': reader.retain(printer_file), 'body': inspect.getsource(NumPyPrinter._print_Heaviside)},
        'scientificCalls': 0})
    def forbidden(*args, **kwargs): raise RuntimeError('saved-input reader cannot call scientific algorithms')
    journal.write = forbidden
    io.native.f.source_jets = io.native.f.polynomial_basis = io.native.f.BasisMomentum.prepare_basis = forbidden
    io.native.Pair.__init__ = io.native.maps = io.native.continue_pair = forbidden
    for name in ('diff', 'lambdify', 'cancel', 'expand', 'factor', 'solve', 'gcd', 'resultant', 'integrate'):
        setattr(sp, name, forbidden)
    basis = metadata('cases/'+p.CASE+'/source-basis-inputs.json')
    saved_basis = metadata('saved-prepared-basis/'+p.CASE+'.json')
    summaries = []
    for row_index, si in ((51, 5), (52, 7), (53, 9), (54, 11)):
        route = metadata('cases/'+p.CASE+'/row-'+str(row_index)+'.json')
        require(route['status'] == 'UNMATCHED_FULL_ROW_INPUT' and route['firstOwner'] == [p.CASE, row_index] and not route['savedCompleteMatches'], 'actual unfinished full row input')
        own, scalars = packet(route['sourceInputs']['packet']), packet(route['scalarInputs'])
        physical = reader.json(route['ownPhysicalRoutes']['packet']['logical'], route['ownPhysicalRoutes']['packet']['sha256'])
        common = packet(physical['context']); context = common['contextPair'][0]
        require(same(*common['contextPair']) and same(*common['basisPair']), 'own physical context and basis')
        row = own['bound']['rows'][row_index]; source = own['bound']['sources'][0, si]
        jet_route = route['jets'][str(si)]; jet = address(jet_route)
        require(len(row['factors']) == 1 and row['factors'][0]['sourceIndex'] == si and
                same(jet['originalBoundAmplitude'], scalars['actual']['source', si]) and
                same(jet['amplitudeUnit'], source['amplitudeUnit']) and same(jet['integralUnit'], source['integralUnit']), 'actual saved whole source coefficients and units')
        candidate = next(v for v in basis if v['sourceIndex'] == si)
        require(candidate['coefficientRoute'] == jet_route and not candidate['savedCoefficientBasisCandidates'], 'saved complete array search result retained')
        coefficient = scalars['actual']['factor', row_index, 0]
        profiles = list(coefficient.atoms(sp.Integral))
        info = {'rowIndex': row_index, 'sourceIndex': si, 'route': route,
                'sourceFrequency': str(source['frequency']), 'sourceAmplitude': str(jet['originalBoundAmplitude']),
                'probe': str(jet['probe']), 'coefficients': [tree(v) for v in jet['coefficients']],
                'amplitudeUnit': jet['amplitudeUnit'], 'integralUnit': jet['integralUnit'],
                'rowFactorUnit': row['factors'][0]['unit'], 'coefficientTree': tree(coefficient),
                'profiles': [{'tree': tree(v), 'unit': own['bound']['profileUnits'][v]} for v in profiles],
                'coefficientFreeSymbols': sorted(map(str, coefficient.free_symbols)),
                'sourceCoefficientFreeSymbols': [sorted(map(str, v.free_symbols)) for v in jet['coefficients']],
                'savedBasisSearch': candidate, 'savedRulePacket': saved_basis,
                'settings': own['settings'], 'physicalRoutes': physical, 'newScientificCalls': 0}
        journal.json('rows/'+str(row_index)+'.json', info)
        summaries.append({'rowIndex': row_index, 'sourceIndex': si, 'coefficientOrders': len(jet['coefficients']),
                          'sourceCoefficients': [str(v) for v in jet['coefficients']],
                          'coefficientTypes': [type(v).__name__ for v in jet['coefficients']],
                          'input': journal.artifacts['rows/'+str(row_index)+'.json']})
    for path in (Path(__file__).resolve(), M/(NAME+'_plan.md'), Path(saved.__file__),
                 M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):
        reader.retain(path)
    for path in (Path(__file__).resolve(), M/(NAME+'_plan.md')):
        dest = base/'source'/path.name; dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out: out.write(path.read_bytes())
        reader.retain(dest, saved.digest(path))
    journal.json('source-coefficient-summary.json', summaries)
    reader.postcheck()
    for rec in journal.artifacts.values(): require(saved.digest(rec['path']) == rec['sha256'], 'stable new metadata')
    journal.json('inputs.json', {'latestAcceptedRow': reader.retain(CP, CP_SHA), 'consumedRoutes': reader.routes})
    checks = {'status': 'COMPLETED_SAVED_REMAINDER_INPUT_INSPECTION', 'rows': summaries,
              'newScientificCalls': 0, 'oldBasisReconstructed': False, 'allConsumedHashesUnchanged': True,
              'artifacts': dict(journal.artifacts), 'wallSeconds': time.monotonic()-started}
    journal.json('checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__': main()
