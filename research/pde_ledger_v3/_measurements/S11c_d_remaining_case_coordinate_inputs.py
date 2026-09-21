#!/usr/bin/env python3
"""Saved case source joins, chart-reuse census and live density operands."""
import argparse
import json
from pathlib import Path
import resource
import shutil
import signal
import time

import sympy as sp
import S11c_d_coordinate_source as coordinate
import S11c_d_remaining_case_first_jet_sources as source

f, engine, native, modes = source.f, source.engine, source.native, source.modes
PLAN = f.M/'S11c_d_remaining_case_coordinate_inputs_plan.md'
FCP = f.M/'S11c_d_remaining_case_first_jet_sources_preflight.json'
CCP = f.M/'S11c_d_coordinate_source_checkpoint.json'
SCOPE = ('Saved four-case source/coordinate reuse census and actual live volume-density '
         'gradient operands. No new coordinate images, operator mutation, numerical '
         'binding, quadrature, current/mode construction or response solve.')


def load(base):
    old, checks, cp = source.provenance.accepted(FCP, 'ACCEPTED_CASE_FIRST_JET_SOURCE_INPUTS')
    cr, cc, publication = source.provenance.accepted(CCP, 'PUBLISHED_ANNEX_VERIFIED')
    path = f.ROOT/publication['publication']['path']
    f.require(path.is_symlink() and f.digest(path) == publication['publication']['sha256'],
              'actual accepted coordinate-source publication')
    pins = dict(checks['sourceFiles'])
    for name, digest in cc['sourceFiles'].items():
        f.require(name not in pins or pins[name] == digest, 'same consumed source versions')
        pins[name] = digest
    for path in (Path(__file__).resolve(), PLAN, FCP, CCP, Path(source.__file__).resolve(),
                 Path(coordinate.__file__).resolve()):
        pins[str(path.relative_to(f.ROOT))] = f.digest(path)
    manifest = {'runDirectory': str(base), 'sourceFiles': pins, 'inputPackets': {},
                'copiedInputs': {}, 'scope': SCOPE, 'input': checks['input'],
                'settings': checks['settings'], 'acceptedSourcePreflight': cp['checksSha256'],
                'acceptedCoordinateChecks': publication['checksSha256'],
                'nativeJoins': {name: source.body(fn) for name, fn in (
                    ('literalRecordJoin', source.require_record), ('literalEquality', native.same),
                    ('selectiveExports', coordinate.selected_exports),
                    ('coordinateChange', coordinate.coordinate_change),
                    ('sourceJets', coordinate.linear_source), ('fieldJets', coordinate.jets))}}
    for origin, result in ((old, checks), (cr, cc)):
        manifest['inputPackets'][str(origin/'checks.json')] = f.digest(origin/'checks.json')
        for name, digest in result['inputPackets'].items():
            f.require(name not in manifest['inputPackets'] or manifest['inputPackets'][name] == digest,
                      'same inherited input versions')
            manifest['inputPackets'][name] = digest
    labels = tuple(checks['cases'])
    for label in labels:
        for relative in (f'accepted-bindings/{label}/case-binding.pickle',
                         *(f'accepted-cases/{label}/{name}.pickle' for name in
                           ('reduced-action', 'actions', 'assembly', 'factorization')),
                         f'input-joins/{label}/source-record-pairs.pickle'):
            modes.retain(old/relative, base/relative, manifest)
    modes.retain(old/'preflight.json', base/'accepted-source-preflight.json', manifest)
    for name, item in cc['artifacts'].items():
        modes.retain(cr/name, base/'accepted-coordinate'/name, manifest, item['sha256'])
    for name, digest in pins.items():
        target = base/'source'/name; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/name, target)
        f.require(f.digest(target) == digest, 'frozen current source join')
    f.save(base/'inputs.json', manifest)
    return manifest, labels


def rational_proof(value):
    numerator, denominator = sp.together(value).as_numer_denom()
    return {'raw': value, 'numerator': numerator, 'denominator': denominator,
            'expandedNumerator': sp.expand(numerator)}


def density_source(density, rule):
    row = next(value for key, value in density if tuple(map(str, key)) == (rule,))
    equations = tuple(dict(row)[sp.core.symbol.Str('VALUE')][0])
    name = {'RHO4_CONSTANT': 'rho_4D_bg_rho4_constant',
            'RHOBR_CONSTANT': 'rho_4D_bg_rhobr_constant'}[rule]
    selected = next(eq for eq in equations if isinstance(eq, sp.Equality) and str(eq.lhs) == name)
    thickness = next(eq for eq in equations if isinstance(eq, sp.Equality) and str(eq.lhs) == 'W_bg')
    return row, selected, thickness


def density_operands(folder, r, geometry, density, rule, accepted):
    row, equation, thickness = density_source(density, rule)
    profile_symbol = next(s for s in thickness.rhs.free_symbols if str(s) == 'w1_profile')
    profile = r.profiles['w'](r.z/r.ell)
    thickness_profile = thickness.rhs.xreplace({profile_symbol: profile})
    rho = equation.rhs.xreplace({thickness.lhs: thickness_profile})
    displacement = geometry['displacement']; X = geometry['coordinates']
    pullback = {r.z: geometry['ansatz'][2]}
    if rule == 'RHO4_CONSTANT':
        f.require(native.same((row, equation, displacement),
                  (accepted['sourceCase'], accepted['densityEquation'], accepted['displacement'])),
                  'exact source and chart operands of accepted density absence')
        gradient = accepted['gradient']; factor = accepted['factor']
        chain = None; material_gradient = None; proofs = {}
        pulled_rho = rho; pulled_gradient = gradient
    else:
        gradient = sp.Matrix([sp.S.Zero, sp.S.Zero, sp.diff(rho, r.z)])
        # Independent density chain rule retains the actual W_bg operand.
        chain = sp.Matrix([0, 0, sp.diff(equation.rhs, thickness.lhs)
                           .xreplace({thickness.lhs: thickness_profile})*sp.diff(thickness_profile, r.z)])
        pulled_rho = rho.xreplace(pullback); pulled_gradient = gradient.xreplace(pullback)
        material_gradient = sp.Matrix([sp.diff(pulled_rho, x) for x in X])
        factor = (displacement.T*pulled_gradient)[0]/pulled_rho
        proofs = {'densityChain': tuple(rational_proof(v) for v in gradient-chain),
                  'covector': tuple(rational_proof(v) for v in
                                   material_gradient-geometry['F'].T*pulled_gradient)}
    dimensions = engine.PHYSICAL_METADATA.dimensions
    dimensions.known[thickness.lhs] = dimensions.measure(thickness_profile)
    dimensions.known[profile_symbol] = (0, 0, 0)
    units = {'density': dimensions.measure(rho), 'thickness': dimensions.measure(thickness_profile),
             'factor': dimensions.measure(factor) if factor != 0 else (0, 0, 0)}
    raw = {'rule': rule, 'sourceRow': row, 'densityEquation': equation,
           'thicknessEquation': thickness, 'profileSymbol': profile_symbol,
           'profile': profile, 'thicknessProfile': thickness_profile, 'physicalDensity': rho,
           'physicalNormalCoordinate': r.z, 'physicalGradient': gradient,
           'chainGradient': chain, 'materialDensity': pulled_rho,
           'pulledPhysicalGradient': pulled_gradient, 'materialGradient': material_gradient,
           'displacement': displacement, 'factor': factor, 'proofs': proofs, 'units': units,
           'originalDenominators': (sp.denom(equation.rhs), sp.denom(rho), sp.denom(factor)),
           'operatorInsertion': 'Uncomputed; this source factor does not isolate a closed-operator term.'}
    if rule == 'RHO4_CONSTANT':
        raw['acceptedAbsence'] = accepted
    else:
        raw['factorMutation'] = {'base': factor, 'reversed': -factor,
                                 'reversalDifference': -2*factor,
                                 'omitted': sp.S.Zero, 'omissionDifference': -factor}
    f.atomic_pickle(folder/'density-operands.pickle', raw)
    f.require(all(p['expandedNumerator'] == 0 for group in proofs.values() for p in group),
              'actual density chain and physical/material covector equations')
    f.require(units['thickness'] == (1, 0, 0) and units['density'] == (-4, 0, 1)
              and units['factor'] == (0, 0, 0), 'actual volume-density and dimensionless advection units')
    if rule == 'RHO4_CONSTANT':
        f.require(native.same(accepted['sourceCase'], row)
                  and native.same(accepted['densityEquation'], equation)
                  and native.same(accepted['gradient'], gradient)
                  and accepted['factor'] == factor == 0, 'source-identical computed structural absence')
    else:
        f.require(gradient[2] != 0 and factor != 0 and equation.rhs.has(thickness.lhs),
                  'actual nonconstant volume density kept live')
        f.require(rational_proof(factor-(-factor))['expandedNumerator'] != 0,
                  'actual reversed density factor responds')
    controls = {'changedDensityFamily': not native.same(equation, density_source(
        density, 'RHOBR_CONSTANT' if rule == 'RHO4_CONSTANT' else 'RHO4_CONSTANT')[1]),
        'changedDensityUnit': units['density'] != (-3, 0, 1),
        'ignoredThicknessChainResponds': gradient[2] != 0 if rule == 'RHOBR_CONSTANT' else None}
    f.require(controls['changedDensityFamily'] and controls['changedDensityUnit']
              and controls['ignoredThicknessChainResponds'] is not False, 'actual density controls')
    f.save(folder/'density-checks.json', {'rule': rule, 'structuralAbsence': factor == 0,
           'gradientNonzeroEntries': sum(v != 0 for v in gradient), 'controls': controls,
           'proofScalars': sum(len(g) for g in proofs.values())})
    return raw


def reuse_signature(kind, expression, unit):
    # The hash only narrows candidates; a full literal comparison is mandatory.
    return kind, hash(expression), tuple(unit)


def prepare(base, manifest, labels):
    accepted = f.unpickle(base/'accepted-coordinate/coordinate-source.pickle')
    geometry = accepted['chart']; atlas = {}; cases = {}; first_context = None; density_families = {}
    selected, census = coordinate.selected_exports(f.ROOT/'scripts/S11c_b_exports.py', {'background_density_map'})
    density_row = selected['background_density_map']; density = density_row['value']
    f.atomic_pickle(base/'density-export.pickle', density_row)
    f.save(base/'density-export-census.json', {'owner': 'scripts/S11c_b_exports.py', 'keys': census})
    for key, item in accepted['records'].items():
        atlas.setdefault(reuse_signature(item['address'][0], item['original'], item['unit']), []).append(
            ({'kind': 'accepted-coordinate', 'case': native.BASELINE, 'key': key,
              'address': item['address']}, item['original'], item['unit']))
    for label in labels:
        folder = base/'cases'/label; folder.mkdir(parents=True)
        case = f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle'); grade = case['grades']
        reduction = f.unpickle(base/'accepted-cases'/label/'reduced-action.pickle')
        r, dimensions = f.prior.domain.momentum.source.native.source.restore_context(reduction)
        dimensions.__dict__.update(grade['dimensionState'])
        for atom, unit in accepted['dimensionState']['known'].items():
            f.require(atom not in dimensions.known or dimensions.known[atom] == unit,
                      'common chart source/dimension identity')
            dimensions.known[atom] = unit
        context = {name: getattr(r, name) for name in ('z', 'zp', 'xi', 'ell', 'omega', 'tangents',
                    'profiles', 'normal_map', 'momentum_groups')}
        context.update(fieldUnits=grade['fieldUnits'], equationUnits=grade['equationUnits'],
                       generators=grade['generators'])
        if first_context is None: first_context = context
        f.atomic_pickle(folder/'chart-context-pair.pickle', (context, first_context))
        f.require(native.same(context, first_context) and native.same(grade['fieldUnits'], accepted['fieldUnits'])
                  and native.same(grade['equationUnits'], accepted['equationUnits']),
                  'actual common chart variables, profiles, fields, units and grade generators')
        pairs = f.unpickle(base/'input-joins'/label/'source-record-pairs.pickle')
        f.require(len(pairs) == len(grade['records']) and {p[0] for p in pairs} == set(grade['records']),
                  'complete accepted native source/address proof census')
        routes = {}; full_pairs = []; baseline_reuse = new_union = cross_case = 0
        for key, address, expression, unit, old in pairs:
            source.require_record(grade['records'][key], address, expression, unit)
            source.require_record(old, address, expression, unit)
            if label == native.BASELINE:
                native_record = accepted['records'][key]
                f.require(native.same((address, expression, unit),
                          (native_record['address'], native_record['original'], native_record['unit'])),
                          'complete original baseline coordinate record join')
            signature = reuse_signature(address[0], expression, unit)
            owner = next((v for v in atlas.get(signature, ()) if native.same(expression, v[1])
                          and tuple(unit) == tuple(v[2])), None)
            if owner is None:
                owner = ({'kind': 'unmapped-source', 'case': label, 'key': key, 'address': address}, expression, unit)
                atlas.setdefault(signature, []).append(owner); new_union += 1
            elif owner[0]['kind'] == 'accepted-coordinate': baseline_reuse += 1
            else: cross_case += 1
            routes[key] = {'address': address, 'sourceUnit': unit, 'owner': owner[0],
                           'coordinateImageExists': owner[0]['kind'] == 'accepted-coordinate'}
            full_pairs.append((key, address, expression, unit, owner))
        f.atomic_pickle(folder/'source-reuse-pairs.pickle', full_pairs)
        f.atomic_pickle(folder/'source-routes.pickle', routes)
        f.atomic_pickle(folder/'term-inputs.pickle', grade['termJoins'])
        changed = next(v for v in pairs if v[2] != 0)
        k, address, expression, unit, old = changed
        controls = {'changedAddress': source.rejects(lambda: source.require_record(old,
                    (address[0], 999, *address[2:]), expression, unit)),
                    'changedCoefficient': source.rejects(lambda: source.require_record(old, address, 2*expression, unit)),
                    'changedUnit': source.rejects(lambda: source.require_record(old, address, expression,
                                                                (unit[0]+1, *unit[1:])))}
        term = grade['termJoins'][0]; limit = term['sourceLimit']; changed_limit = (limit[0], limit[1], limit[2]+r.ell)
        f.atomic_pickle(folder/'mutation-operands.pickle', {'record': changed, 'sourceLimit': limit,
            'changedSourceLimit': changed_limit, 'changedCoefficient': 2*expression})
        controls['changedSourceLimit'] = not native.same(limit, changed_limit)
        f.require(all(controls.values()), 'actual record/address/coefficient/unit/ordered-limit mutation controls')
        rule = label.split('__')[1]; row, equation, thickness = density_source(density, rule)
        signature = (row, equation, thickness, r.profiles['w'], r.z, r.ell)
        if rule in density_families:
            owner_label, old_signature, value = density_families[rule]
            f.atomic_pickle(folder/'density-family-pair.pickle', (signature, old_signature))
            f.require(native.same(signature, old_signature), 'full same-density source/profile family join')
            modes.retain(base/'cases'/owner_label/'density-operands.pickle', folder/'density-operands.pickle', manifest)
            modes.retain(base/'cases'/owner_label/'density-checks.json', folder/'density-checks.json', manifest)
        else:
            value = density_operands(folder, r, geometry, density, rule, accepted['densityAdvection'])
            density_families[rule] = label, signature, value
        cases[label] = {'records': len(routes), 'rows': len(case['binding']['bound']['rows']),
                        'terms': len(grade['termJoins']),
                        'sources': sum(v['address'][0] == 'source' for v in grade['records'].values()),
                        'acceptedCoordinateRecordReuses': baseline_reuse, 'newUnionOperands': new_union,
                        'unmappedSharedOperands': cross_case, 'controls': controls,
                        'density': json.loads((folder/'density-checks.json').read_text())}
        f.save(base/'case-inventory.json', cases)
        f.atomic_pickle(folder/'dimension-state.pickle', dict(vars(dimensions)))
    f.save(base/'inputs.json', manifest)
    return cases


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args(); base = args.run_directory.resolve(); base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900); started = time.monotonic()
    manifest, labels = load(base); cases = prepare(base, manifest, labels)
    for name, digest in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == digest, 'source/frozen pre/post identity')
    for name, digest in manifest['inputPackets'].items():
        f.require(f.digest(Path(name)) == digest, 'original input pre/post identity')
    for name, digest in manifest['copiedInputs'].items():
        f.require(f.digest(base/name) == digest, 'complete copied operand identity')
    artifacts = {str(p.relative_to(base)): {'sha256': f.digest(p), 'bytes': p.stat().st_size}
                 for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts
                 and p.name not in ('inputs.json', 'checks.json')}
    checks = {**manifest, 'status': 'COMPLETED_CASE_COORDINATE_SOURCE_INPUTS', 'cases': cases,
              'newCoordinateImages': 0, 'newSourceJets': 0, 'newNumericalWork': 0,
              'artifacts': artifacts, 'wallSeconds': time.monotonic()-started,
              'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
