#!/usr/bin/env python3
"""Saved-input material end maps, constructed before any response solve."""
import argparse
import ast
import copy
import gc
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import textwrap
import time
from types import SimpleNamespace

import numpy as np
import scipy.linalg as la
import sympy as sp
import S11c_d_remaining_case_coordinate_boundary_inputs as inputs
import S11c_d_finite_scattering_domain as domain

f, h, c, engine, modes = inputs.f, inputs.h, inputs.c, inputs.engine, inputs.modes
matrices, ends_native, b = inputs.matrices, inputs.ends_native, inputs.c.boundary
CP = f.M/'S11c_d_remaining_case_coordinate_boundary_focused.json'
PLAN = f.M/'S11c_d_remaining_case_coordinate_boundary_production_plan.md'
WIRING = f.M/'S11c_d_remaining_case_coordinate_boundary_production_wiring.json'
SCOPE = ('Actual finite and independent-grade material end/current/phase maps in common Eulerian '
         'coordinates before boundary replacement. No quadrature, modal/current closure or response solve. '
         'Finite current is open-only; closed matching amplitudes are separate. Historical material response is continuum only.')


def tree(fun):
    return ast.parse(textwrap.dedent(inspect.getsource(fun))).body[0]


def ast_sha(node):
    return hashlib.sha256(ast.dump(node).encode()).hexdigest()


def compile_function(node, namespace):
    module = ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[]))
    env = dict(namespace)
    exec(compile(module, '<saved-material-boundary-adapter>', 'exec'), env)
    return env[node.name]


def native_adapters():
    original = tree(f.boundary_map)
    start = next(i for i, n in enumerate(original.body)
                 if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]) == 'right')
    tail = copy.deepcopy(original.body[start:])
    wrapper = ast.parse('def finite_tail(outgoing, incoming, excluded, packet, orientation):\n pass').body[0]
    wrapper.body = tail
    f.require(ast.dump(ast.Module(body=wrapper.body, type_ignores=[])) ==
              ast.dump(ast.Module(body=original.body[start:], type_ignores=[])), 'entire native finite boundary assembly tail')
    finite_tail = compile_function(wrapper, vars(f))
    original_continuum = tree(c.material_ends)
    transformed = copy.deepcopy(original_continuum)
    # Checkpoints only. Removing them recovers the entire native function AST.
    class Observe(ast.NodeTransformer):
        def visit_Assign(self, node):
            key = ast.unparse(node.targets[0])
            if key == 'all_symbols':
                before = ast.parse("observe(base, end+'-prepared', dict(material=material, common=common, md=md, ed=ed, mi=mi, ei=ei, X=X, U=U, T=T, phase_in=phase_in, phase_out=phase_out, transported_trace=transported_trace, transported_insert=transported_insert))").body[0]
                return [before, node]
            if key == 'pair' and isinstance(node.value, ast.Call) and ast.unparse(node.value.func) == 'boundary.current_pair':
                before = ast.parse("observe(base, end+'-'+name+'-pair-'+str(i)+'-'+str(j)+'-input', dict(table=transformed, left=left, right=right))").body[0]
                after = ast.parse("observe(base, end+'-'+name+'-pair-'+str(i)+'-'+str(j)+'-value', pair)").body[0]
                return [before, node, after]
            return self.generic_visit(node)
        def visit_Expr(self, node):
            if isinstance(node.value, ast.Call) and ast.unparse(node.value.func) == 'f.require' and 'complete common-frame material boundary/phase/current route' in ast.unparse(node):
                before = ast.parse("observe(base, end+'-complete-before-guard', dict(material=material, common=common, md=md, ed=ed, currents=currents, material_currents=material_currents, differences=differences, mi=mi, ei=ei, X=X, U=U, T=T, phase_in=phase_in, phase_out=phase_out))").body[0]
                return [before, node]
            return self.generic_visit(node)
    transformed = Observe().visit(transformed)
    class Remove(ast.NodeTransformer):
        def visit_Expr(self, node):
            if isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == 'observe':
                return None
            return self.generic_visit(node)
    reverse = Remove().visit(copy.deepcopy(transformed))
    f.require(ast.dump(reverse) == ast.dump(original_continuum), 'whole native material_ends reverse AST: observers only')
    def observe(base, name, value):
        (base/'native-checkpoints').mkdir(exist_ok=True)
        f.atomic_pickle(base/'native-checkpoints'/(name+'.pickle'), value)
    continuum = compile_function(transformed, dict(vars(c), observe=observe))
    # Literal native Piola/field/normal-covector table transformation, including
    # its Taylor derivative factors. A finite evaluated source uses order (0,0).
    assignments = [n for n in ast.walk(original_continuum) if isinstance(n, ast.Assign)
                   and ast.unparse(n.targets[0]) == 'transformed']
    f.require(len(assignments) == 1, 'one actual native material-current pullback')
    pull = ast.parse('def pullback(table, chart, covectors):\n pass').body[0]
    pull.body = [copy.deepcopy(assignments[0]), ast.Return(value=ast.Name(id='transformed', ctx=ast.Load()))]
    pullback = compile_function(pull, vars(c))
    joins = {'finiteBoundaryWholeSource': ast_sha(original), 'finiteBoundaryTail': ast_sha(wrapper),
             'finiteTailUnchanged': True, 'materialEndsWholeSource': ast_sha(original_continuum),
             'materialEndsObserved': ast_sha(transformed), 'materialEndsReverseAST': True,
             'nativeCurrentPullback': ast_sha(assignments[0]), 'currentPair': matrices.body(b.current_pair),
             'domainPhases': matrices.body(domain.phases), 'arrayUnits': matrices.body(ends_native.array_units),
             'restore': matrices.body(h.restore), 'restoreChart': matrices.body(h.restore_chart),
             'chartMaps': matrices.body(c.Chart.maps)}
    joins['acceptedNativeBodies'] = {
        'material_ends': ast_sha(ast.Module(body=[original_continuum], type_ignores=[])),
        'boundary_map': ast_sha(ast.Module(body=[original], type_ignores=[])),
        'current_pair': matrices.body(b.current_pair), 'array_units': matrices.body(ends_native.array_units),
        'restore': matrices.body(h.restore), 'restore_chart': matrices.body(h.restore_chart)}
    return finite_tail, continuum, pullback, joins


def load(base, resume):
    origin, checks, checkpoint = matrices.accepted(CP, 'ACCEPTED_CASE_MATERIAL_BOUNDARY_INPUTS')
    f.require(origin == resume and checks['result']['newBoundaryConstructions'] == 0,
              'accepted saved material boundary input directory')
    manifest = {'runDirectory': str(base), 'sourceFiles': dict(checks['sourceFiles']),
                'inputPackets': dict(checks['inputPackets']), 'copiedInputs': {},
                'input': checks['input'], 'settings': checks['settings'], 'scope': SCOPE}
    for n, sha in checks['sourceFiles'].items():
        f.require(f.digest(f.ROOT/n) == f.digest(origin/'source'/n) == sha, 'accepted current/frozen native implementation')
    for name, item in checks['artifacts'].items():
        modes.retain(origin/name, base/name, manifest, item['sha256'])
    for name, target in (('checks.json', 'accepted-focus-checks.json'), ('inputs.json', 'accepted-focus-inputs.json')):
        modes.retain(origin/name, base/target, manifest)
    for path in (Path(__file__).resolve(), PLAN, WIRING, CP, Path(domain.__file__).resolve()):
        manifest['sourceFiles'][str(path.relative_to(f.ROOT))] = f.digest(path)
    manifest['completedFocusReuse'] = {'directory': str(origin), 'checksSha256': checkpoint['checksSha256'],
                                      'artifacts': len(checks['artifacts']), 'validator': checkpoint['validator']}
    for name, sha in manifest['sourceFiles'].items():
        target = base/'source'/name; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/name, target); f.require(f.digest(target) == sha, 'frozen new boundary implementation')
    f.save(base/'inputs.json', manifest)
    matrices.hash_check(base, manifest)
    return manifest, tuple(checks['result']['cases']), checks['result']


def restrict_science():
    def forbidden(*a, **kw):
        raise RuntimeError('source/binding/quadrature/modal/current closure/response construction prohibited in boundary adapter')
    for cls in (engine.NumericalReducedAction, engine.ModalCurrentSubspaces, engine.TwoEndedMatchingChannels):
        cls.__init__ = forbidden
    engine.NumericalReducedAction.bind = forbidden
    c.Chart.__init__ = forbidden; c.Chart.maps = forbidden; c.Chart.image = forbidden; c.Chart.basis = forbidden
    c.source.coordinate_change = forbidden; c.bind = forbidden; c.MaterialMomentum.prepare_basis = forbidden
    matrices.interior.assemble = forbidden; matrices.matrices.direct_cells = forbidden
    f.boundary_map = forbidden; f.source_jets = forbidden; b.construct_end = forbidden
    c.response.systems = forbidden; c.response.solve = forbidden; c.response.channels = forbidden; c.response.open_flux = forbidden
    # The original boundary constructors need five-dimensional coordinate
    # inverses, rank/condition checks only. Never enable a response or root solve.
    allowed = {}
    for module in (np.linalg, la):
        for name in ('solve', 'inv', 'svd', 'matrix_rank', 'cond'):
            if not hasattr(module, name): continue
            original = getattr(module, name)
            def bounded(a, *args, _original=original, _name=name, **kw):
                f.require(np.shape(a) == (5, 5), ('only actual five-field boundary coordinate algebra', _name))
                return _original(a, *args, **kw)
            setattr(module, name, bounded); allowed[module.__name__+'.'+name] = '5x5 only'
        for name in ('pinv', 'eig', 'eigh', 'eigvals', 'eigvalsh', 'lu_factor', 'lu_solve'):
            if hasattr(module, name): setattr(module, name, forbidden)
    return allowed


def chart_maps(base, state, position):
    end = 'left' if position < 0 else 'right'
    saved = f.unpickle(base/'baseline-material'/(end+'-material-boundary.pickle'))['material']
    g = state['values']; X, U, T = saved['position'], saved['materialFromEulerian'], saved['eulerianFromMaterial']
    f.require(float(position) == (-1 if end == 'left' else 1)*64.0 and
              X == position/float(g['d']) and saved['tangentialMeasure'] == float(g['A']),
              'approved saved common chart endpoint')
    f.require(b.norm([T@U-np.eye(5), U@T-np.eye(5)]) < 1e-12, 'saved actual chart inverse maps')
    return X, U.copy(), T.copy()


def finite_current(folder, channel, state, X, U, T, pullback):
    chart = SimpleNamespace(**state['values']); A, d, kappa = map(float, (chart.A, chart.d, chart.kappa))
    candidates = {item['INFO']['INDEX']: item for item in channel['CANDIDATES']}
    labels = channel['CHANNELS']; size = len(labels)
    f.require(size == len(channel['OUTWARD_CURRENT']) and
              sorted(item['MATRIX_COLUMN'] for item in labels) == list(range(size)), 'complete actual open current columns')
    # All elements must be supplied by actual root pairs. No absent entry is zero.
    field = np.empty((size, size), complex); flux = np.empty_like(field)
    slab_form = np.empty_like(field); bulk_form = np.empty_like(field); covered = np.zeros((size, size), int)
    pair_inventory = []; differences = []; coefficient_controls = []
    for index, pair in enumerate(channel['PAIR_OPERANDS']):
        path = folder/'pairs'/str(index); path.mkdir(parents=True)
        left, right = (candidates[pair[n]] for n in ('LEFT_RECORD_INDEX', 'RIGHT_RECORD_INDEX'))
        li = sorted((v for v in labels if v['RECORD_INDEX'] == left['INFO']['INDEX']), key=lambda x:x['MATRIX_COLUMN'])
        ri = sorted((v for v in labels if v['RECORD_INDEX'] == right['INFO']['INDEX']), key=lambda x:x['MATRIX_COLUMN'])
        f.require(len(li) == left['INFO']['NULLITY'] and len(ri) == right['INFO']['NULLITY'], 'complete degenerate open bases')
        physical = (complex(left['INFO']['K']), complex(right['INFO']['K']))
        km = tuple(d*k+kappa for k in physical)
        q = tuple(complex(v['INFO']['PHYSICAL_Q']) for v in (left, right))
        f.atomic_pickle(path/'input.pickle', {'pair': pair, 'left': left, 'right': right, 'leftChannels': li,
            'rightChannels': ri, 'physicalK': physical, 'materialK': km, 'q': q, 'U': U, 'T': T,
            'chart': state, 'zeroTaylorDerivativeOrder': (0, 0)})
        f.require(all(v.imag > 0 for v in q) and b.norm([(k-kappa)/d-p for k,p in zip(km, physical)]) < 1e-12,
                  'actual inverse normal-covector evaluation and unchanged bulk decay domain')
        f.require(np.isfinite(pair['DEPTH_INTEGRAL']) and np.isfinite(pair['DEPTH_RATE']), 'saved finite depth operands')
        material = {}
        for name in ('CURRENT_SLAB', 'CURRENT_BULK'):
            value = pullback({(0, 0, 0, 0): sp.ImmutableMatrix(pair[name])}, chart, {})[(0, 0, 0, 0)]
            material[name] = np.asarray(value, complex)
        material['COMPOSED_CURRENT'] = material['CURRENT_SLAB']+pair['DEPTH_INTEGRAL']*material['CURRENT_BULK']
        lb, rb = left['BASES'], right['BASES']; lf, rf = U@lb['FLUX_RIGHT'], U@rb['FLUX_RIGHT']
        values = {'field': (U@lb['RIGHT']).conj().T@material['COMPOSED_CURRENT']@(U@rb['RIGHT']),
                  'flux': lf.conj().T@material['COMPOSED_CURRENT']@rf,
                  'slab': lf.conj().T@material['CURRENT_SLAB']@rf,
                  'bulk': lf.conj().T@(pair['DEPTH_INTEGRAL']*material['CURRENT_BULK'])@rf}
        f.atomic_pickle(path/'value.pickle', {'materialSources': material, 'materialContractions': values,
                                            'commonEulerianContractions': {n:A*v for n,v in values.items()}})
        if not coefficient_controls and b.norm(values['flux']) > 1e-10:
            altered_source = 2*pair['COMPOSED_CURRENT']
            altered_material = np.asarray(pullback({(0, 0, 0, 0): sp.ImmutableMatrix(altered_source)}, chart, {})[(0, 0, 0, 0)], complex)
            altered_contraction = A*(lf.conj().T@altered_material@rf)
            control = {'originalSource': pair['COMPOSED_CURRENT'], 'changedSource': altered_source,
                       'originalMaterialSource': material['COMPOSED_CURRENT'], 'changedMaterialSource': altered_material,
                       'originalContraction': A*values['flux'], 'changedContraction': altered_contraction,
                       'response': altered_contraction-A*values['flux']}
            f.atomic_pickle(path/'coefficient-control.pickle', control)
            coefficient_controls.append(control)
        differences.append(pair['COMPOSED_CURRENT']-pair['CURRENT_SLAB']-pair['DEPTH_INTEGRAL']*pair['CURRENT_BULK'])
        for a in li:
            for z in ri:
                i, j = a['MATRIX_COLUMN'], z['MATRIX_COLUMN']; ia, ja = a['BASIS_COLUMN'], z['BASIS_COLUMN']
                field[i,j] = values['field'][ia,ja]; flux[i,j] = values['flux'][ia,ja]
                slab_form[i,j] = values['slab'][ia,ja]; bulk_form[i,j] = values['bulk'][ia,ja]; covered[i,j] += 1
        pair_inventory.append({'index': index, 'records': (left['INFO']['INDEX'], right['INFO']['INDEX']),
                               'inputSha256': f.digest(path/'input.pickle'), 'valueSha256': f.digest(path/'value.pickle')})
    f.require(np.all(covered == 1), 'all actual finite open-current entries supplied exactly once')
    orientation = labels[0]['END']; sign = {'LEFT': -1, 'RIGHT': 1}[orientation]
    common = {'FIELD_CURRENT': A*field, 'FLUX_CURRENT': A*flux, 'OUTWARD_CURRENT': sign*A*flux,
              'slab': sign*A*slab_form, 'bulk': sign*A*bulk_form}
    difference = {'sourceComposition': differences,
                  'field': common['FIELD_CURRENT']-channel['FIELD_CURRENT'],
                  'flux': common['FLUX_CURRENT']-channel['FLUX_CURRENT'],
                  'outward': common['OUTWARD_CURRENT']-channel['OUTWARD_CURRENT'],
                  'basisChange': common['FLUX_CURRENT']-channel['FIELD_TO_FLUX_MAP'].conj().T@common['FIELD_CURRENT']@channel['FIELD_TO_FLUX_MAP'],
                  'hermitian': common['OUTWARD_CURRENT']-common['OUTWARD_CURRENT'].conj().T,
                  'decomposition': common['OUTWARD_CURRENT']-common['slab']-common['bulk']}
    omitted = sign*flux-common['OUTWARD_CURRENT']
    f.require(len(coefficient_controls) == 1, 'actual responding finite source coefficient operand')
    mutation = coefficient_controls[0]['response']
    result = {'materialField': field, 'materialFlux': flux, 'materialOutward': sign*flux,
              'commonEulerian': common, 'coverage': covered, 'pairs': pair_inventory,
              'differences': difference, 'omittedTangentialMeasure': omitted, 'coefficientMutation': mutation,
              'scope': 'Open-only native finite current form. No closed-current entries constructed or supplied.'}
    f.atomic_pickle(folder/'open-current-route.pickle', result)
    f.require(b.norm(difference) < 1e-8 and b.norm(omitted) > 1e-10 and b.norm(mutation) > 1e-10,
              'finite source current pullback, actual measure and coefficient controls')
    return result


def finite_construct(folder, original, channels, state, settings, maps, finite_tail, pullback):
    folder.mkdir(parents=True)
    X, U, T = maps; g = state['values']; d, kappa = float(g['d']), float(g['kappa'])
    x = original['orientation']*settings['sourceBound']
    candidates = {v['INFO']['INDEX']: v for v in channels['CANDIDATES']}
    f.atomic_pickle(folder/'input.pickle', {'finite': original, 'channels': channels, 'chart': state,
                                          'settings': settings, 'position': x, 'X': X, 'U': U, 'T': T})
    material, common = {}, {}; directions = []
    for name in ('outgoing', 'incoming'):
        material[name], common[name] = [], []
        for item in original[name]:
            native = candidates[item['RECORD_INDEX']]; basis_name = 'FLUX_RIGHT' if item['kind'] == 'open' else 'RIGHT'
            f.require(item['k'] == complex(native['INFO']['K']) and
                      np.array_equal(item['vector'], native['BASES'][basis_name][:,item['BASIS_COLUMN']]),
                      'actual finite root, degenerate basis and physical address')
            mv = U@item['vector']; mk = d*item['k']+kappa
            material[name].append(dict(item, vector=mv, k=mk))
            common[name].append(dict(item, vector=T@mv, k=(mk-kappa)/d))
            directions.append({'address': (name, item['RECORD_INDEX'], item['BASIS_COLUMN']),
                'original': item, 'material': material[name][-1], 'common': common[name][-1],
                'materialDerivative': 1j*mk*mv,
                'commonDerivative': T@(1j*mk*mv-1j*kappa*mv)/d})
    # Preserve every candidate and complete original dual basis. Only right
    # field components are transported here; no new adjoint normalization.
    full = [{'original': v, 'materialRightBases': {n:U*a for n,a in v['BASES'].items() if n in ('RIGHT','FLUX_RIGHT')},
             'materialK': d*complex(v['INFO']['K'])+kappa,
             'commonRightBases': {n:T@(U@a) for n,a in v['BASES'].items() if n in ('RIGHT','FLUX_RIGHT')}}
            for v in channels['CANDIDATES']]
    f.atomic_pickle(folder/'transported-candidates.pickle', full)
    f.atomic_pickle(folder/'selected-directions.pickle', directions)
    currents = finite_current(folder, channels, state, X, U, T, pullback)
    # Full actual first-write operands are checkpointed before either small
    # native five-field boundary-coordinate solve. These are not response solves.
    finite = {}
    for name, records, metric in (('material', material, currents['materialOutward']),
                                 ('common', common, currents['commonEulerian']['OUTWARD_CURRENT'])):
        request = {'outgoing': records['outgoing'], 'incoming': records['incoming'],
                   'excluded': original['excluded'], 'packet': {'OUTWARD_CURRENT': metric},
                   'orientation': original['orientation']}
        f.atomic_pickle(folder/(name+'-boundary-input.pickle'), request)
        finite[name] = finite_tail(**request)
        f.atomic_pickle(folder/(name+'-boundary.pickle'), finite[name])
    md, ed = finite['material'], finite['common']
    phases = {}
    for name, direction, sense in (('incoming', 'incoming', 1), ('outgoing', 'outgoing', -1)):
        selected = [v for v in material[direction] if v['kind'] == 'open']
        physical = [v for v in original[direction] if v['kind'] == 'open']
        f.require(all(v['k'].imag == 0 for v in physical), 'native common-origin phase is open-only')
        phases[name] = {'material': np.asarray([np.exp(sense*1j*v['k']*X) for v in selected]),
            'shearFactor': np.exp(-sense*1j*kappa*X),
            'nativePhysical': np.asarray([np.exp(sense*1j*v['k']*x) for v in physical])}
        phases[name]['commonEulerian'] = phases[name]['material']*phases[name]['shearFactor']
    differences = {n: ed[n]-original[n] for n in ('right','derivative','traceMap','incomingValues',
                    'incomingDerivative','incomingBoundaryData','current')}
    differences.update(traceRoute=T@md['traceMap']@U/d-1j*kappa/d*np.eye(5)-ed['traceMap'],
                       insertionRoute=T@md['incomingBoundaryData']/d-ed['incomingBoundaryData'],
                       derivativeRoute=T@(md['derivative']-1j*kappa*md['right'])/d-ed['derivative'],
                       incomingDerivativeRoute=T@(md['incomingDerivative']-1j*kappa*md['incomingValues'])/d-ed['incomingDerivative'],
                       phases={n:v['commonEulerian']-v['nativePhysical'] for n,v in phases.items()},
                       fullCandidateBases=[v['commonRightBases'][n]-v['original']['BASES'][n] for v in full for n in v['commonRightBases']])
    controls = {'omittedDerivativeCovector': T@md['derivative']/d-ed['derivative'],
                'omittedShearPhase': {n:v['material']-v['nativePhysical'] for n,v in phases.items()},
                'omittedMeasure': currents['omittedTangentialMeasure'],
                'currentCoefficient': currents['coefficientMutation']}
    result = {'material': md, 'commonEulerian': ed, 'phases': phases, 'differences': differences,
              'controls': controls, 'position': x, 'materialPosition': X,
              'openCurrentPath': str(folder/'open-current-route.pickle'),
              'closedMatchingDirections': [v for v in ed['outgoing'] if v['kind'] != 'open'],
              'scope': SCOPE}
    f.atomic_pickle(folder/'finite-material-boundary.pickle', result)
    f.require(b.norm(differences) < 1e-8 and all(b.norm(v) > 1e-10 for v in controls.values()),
              'complete material finite boundary/current/phase and physical omission controls')
    f.require(len(ed['incoming']) == 2 and len(ed['outgoing']) == 5 and len(result['closedMatchingDirections']) == 3,
              'actual saved finite open/closed matching census')
    return result


def construct(base, manifest, labels, accepted, finite_tail, continuum_constructor, pullback):
    state = f.unpickle(base/'bindings/chart-state.pickle')
    inventory = {}; families = {}; outputs = {}; new_continuum = 0; new_finite = 0
    # Use the saved input-family census, then prove complete input identity again
    # at each reuse. No symbolic family label supplies a computed boundary map.
    for label in labels:
        case = inputs.actual_address(base, label, 'LEFT', labels)
        outputs[label] = {}
        for end in ('LEFT', 'RIGHT'):
            address = label+'__'+end; route = accepted['ends'][address]; owner = route['owner']
            actual = f.unpickle(base/'preparation'/label/end.lower()/'end-inputs.pickle')
            expected = inputs.signature(case['ends'][end], case['finite'][end], state,
                                        actual['settings'], {n:case[n] for n in ('fieldUnits','rowUnits','currentUnit')})
            f.require(modes.same(actual, expected), 'full actual own-case boundary input before construction')
            directory = base/'material-families'/owner
            channel = f.unpickle(base/'boundary-cases'/label/end.lower()/'channels.pickle')
            if owner not in families:
                directory.mkdir(parents=True)
                f.atomic_pickle(directory/'input.pickle', actual)
                maps = chart_maps(base, state, actual['finite']['orientation']*actual['settings']['sourceBound'])
                f.atomic_pickle(directory/'saved-chart-map-reuse.pickle', {'maps': maps, 'chart': state,
                    'source': str(base/'baseline-material'/(end.lower()+'-material-boundary.pickle'))})
                finite = finite_construct(directory/'finite', actual['finite'], channel, state,
                                          actual['settings'], maps, finite_tail, pullback)
                new_finite += 1
                continuum_dir = directory/'continuum'; continuum_dir.mkdir()
                if owner in ('LEFT', 'RIGHT'):
                    for name in ('material-boundary', 'slab-material-current-tables', 'bulk-material-current-tables'):
                        src = base/'baseline-material'/(owner.lower()+'-'+name+'.pickle')
                        modes.retain(src, continuum_dir/src.name, manifest)
                    continuum = f.unpickle(continuum_dir/(owner.lower()+'-material-boundary.pickle'))
                    disposition = 'accepted baseline material end and full source tables reused'
                else:
                    data, bound, packets = h.restore(base/'bindings', label, manifest, saved=True)
                    f.require(modes.same(data['system']['settings'], actual['settings']) and
                              modes.same(data['coordinate']['chart'], state['values']['g']) and
                              modes.same(data['coordinate']['fieldJets'], state['values']['jets']),
                              'actual restored own-source settings/chart/field jets before native material end')
                    chart = h.restore_chart(data, state)
                    chart.maps = lambda position: chart_maps(base, state, position)
                    data['ends'] = {n:case[n] for n in ('fieldUnits','rowUnits','currentUnit')}
                    data['ends']['ends'] = {end: actual['continuum']}
                    f.atomic_pickle(continuum_dir/'native-input.pickle', {'ends': data['ends'], 'settings': data['system']['settings'],
                        'chart': state, 'case': label, 'end': end})
                    whole = continuum_constructor(continuum_dir, data, chart)
                    f.atomic_pickle(continuum_dir/'dimension-state.pickle', engine.PHYSICAL_METADATA.dimensions.__dict__)
                    continuum = f.unpickle(continuum_dir/(end.lower()+'-material-boundary.pickle'))
                    f.require(modes.same(whole['ends']['ends'][end], continuum['commonEulerian']), 'complete new native continuum result')
                    new_continuum += 1; disposition = 'new native material end from own complete current/source input'
                    del data, bound, packets, chart, whole
                unit = ends_native.array_units(finite['commonEulerian'], continuum['commonEulerian'],
                                               case['fieldUnits'], case['rowUnits'], case['currentUnit'])
                f.atomic_pickle(directory/'array-units.pickle', unit)
                field_units = case['fieldUnits']; current_unit = case['currentUnit']
                f.atomic_pickle(directory/'source-current-units.pickle', {
                    'integrated': tuple(tuple(tuple(q-a-z for q,a,z in zip(current_unit,u,v))
                                             for v in field_units) for u in field_units),
                    'bulkBeforeDepthIntegration': tuple(tuple(tuple(q-a-z-(1 if i == 0 else 0)
                        for i,(q,a,z) in enumerate(zip(current_unit,u,v))) for v in field_units) for u in field_units),
                    'normalMomentumUnit': (-1,0,0), 'materialFieldUnits': field_units,
                    'chartMapUnits': 'dimensionless', 'zeroTaylorDerivativeOrder': (0,0)})
                families[owner] = {'input': actual, 'finite': finite, 'continuum': continuum, 'units': unit,
                                   'channels': channel, 'sourceLabel': label, 'sourceEnd': end, 'directory': str(directory)}
                summary = {'owner': owner, 'sourceLabel': label, 'sourceEnd': end, 'directory': str(directory),
                    'finiteMaximum': b.norm(finite['differences']), 'finiteControls': {n:b.norm(v) for n,v in finite['controls'].items()},
                    'continuumMaximum': b.norm(continuum['proofs']['differences']), 'unitArrays': len(unit),
                    'continuumDisposition': disposition}
                f.save(directory/'checks.json', summary)
            else:
                f.require(modes.same(actual, families[owner]['input']), 'entire family input identity before completed map reuse')
                f.require(modes.same(channel, families[owner]['channels']), 'full candidate, dual basis and finite current packet before family reuse')
            result = families[owner]
            outputs[label][end] = {'finite': result['finite']['commonEulerian'], 'continuum': result['continuum']['commonEulerian'],
                                   'finitePhases': result['finite']['phases'], 'owner': owner}
            inventory[address] = {'owner': owner, 'sourceLabel': result['sourceLabel'], 'sourceEnd': result['sourceEnd'],
                'inputSha256': f.digest(base/'preparation'/label/end.lower()/'end-inputs.pickle'),
                'familyInputSha256': f.digest(Path(result['directory'])/'input.pickle'),
                'finitePacketSha256': f.digest(Path(result['directory'])/'finite/finite-material-boundary.pickle'),
                'continuumPacketSha256': f.digest(Path(result['directory'])/'continuum'/(result['sourceEnd'].lower()+'-material-boundary.pickle')),
                'finiteIncoming': len(result['finite']['commonEulerian']['incoming']),
                'finiteOutgoing': len(result['finite']['commonEulerian']['outgoing']),
                'finiteCurrentShape': list(result['finite']['commonEulerian']['current'].shape),
                'continuumDirections': int(result['continuum']['commonEulerian']['offsets'][-1])}
            f.save(base/'material-boundary-inventory.json', inventory)
        packet = {n:case[n] for n in ('fieldUnits','rowUnits','currentUnit')}
        packet.update(case=label, finite={e:v['finite'] for e,v in outputs[label].items()},
            ends={e:v['continuum'] for e,v in outputs[label].items()},
            finitePhases={e:v['finitePhases'] for e,v in outputs[label].items()},
            materialFamilyRoutes={e:v['owner'] for e,v in outputs[label].items()}, scope=SCOPE)
        (base/'material-cases'/label).mkdir(parents=True)
        f.atomic_pickle(base/'material-cases'/label/'case-material-boundary.pickle', packet)
        del case, actual, expected, packet
        gc.collect()
    result = {'cases': list(labels), 'addresses': inventory, 'newContinuumFamilies': new_continuum,
              'newFiniteFamilies': new_finite, 'originalContinuumReuseAddresses': sum(v['owner'] in ('LEFT','RIGHT') for v in inventory.values()),
              'newResponseSolves': 0, 'baselineMaterialResponseIsContinuumOnly': True,
              'finiteMaximum': max(b.norm(v['finite']['differences']) for v in families.values()),
              'continuumMaximum': max(b.norm(v['continuum']['proofs']['differences']) for v in families.values())}
    f.atomic_pickle(base/'remaining-case-material-boundary.pickle', {'result': result, 'caseRoutes': outputs, 'scope': SCOPE})
    f.require(new_continuum == len(accepted['newMaterialFamilies']) == 1 and new_finite == len(families) == 3
              and len(inventory) == 8, 'all complete actual families and end addresses')
    return result


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory', type=Path, required=True)
    ap.add_argument('--resume-focused', type=Path, required=True); args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900); start = time.monotonic()
    base = args.run_directory.resolve(); base.relative_to(f.STORE); base.mkdir(parents=True, exist_ok=False)
    finite_tail, continuum, pullback, joins = native_adapters()
    restrictions = restrict_science()
    manifest, labels, accepted = load(base, args.resume_focused.resolve())
    inherited = json.loads((base/'native-boundary-joins.json').read_text())
    f.require(all(inherited[n] == value for n,value in joins['acceptedNativeBodies'].items()),
              'consumed native definitions joined to accepted input preparation')
    f.save(base/'material-constructor-joins.json', {'native': joins, 'allowedCoordinateAlgebra': restrictions})
    result = construct(base, manifest, labels, accepted, finite_tail, continuum, pullback)
    f.save(base/'inputs.json', manifest); matrices.hash_check(base, manifest)
    artifacts = {str(p.relative_to(base)): {'sha256': f.digest(p), 'bytes': p.stat().st_size}
        for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts
        and p not in (base/'inputs.json', base/'checks.json')}
    checks = {**manifest, 'status': 'COMPLETED_CASE_MATERIAL_BOUNDARY_MAPS', 'result': result,
              'nativeJoins': joins, 'artifacts': artifacts, 'wallSeconds': time.monotonic()-start,
              'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__': main()
