#!/usr/bin/env python3
"""Source-bound one-case controls for the native closed current builders."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import pickle
import resource
import time

import sympy as sp

from S11c_d_modal_subspace_check import load_current, engine
from S11c_d_modal_current_check import digest, source_node
from S11c_d_joint_sheet_check import decoded_lines, _restore


def reuse_control(args, label, provenance, replacements, extra_bindings):
    """Reuse only completed source packets with identical non-modal constructors."""
    if args.reuse_checkpoint is None:
        return None, None
    checkpoint = json.loads(args.reuse_checkpoint.read_text())
    if label not in checkpoint['controls']:
        return None, None
    base = Path(checkpoint['runDirectory'])
    for name,sha in checkpoint['sourceFiles'].items():
        if digest(base/'source'/name)!=sha:
            raise ValueError(('source-control snapshot mismatch',name))
    for key,value in checkpoint['provenance'].items():
        if key not in ('engineSha256','instrumentSha256') and json.loads(json.dumps(provenance.get(key)))!=value:
            raise ValueError(('source-control cache dependency mismatch',key))
    def nonmodal_source(path):
        module = ast.parse(path.read_text())
        module.body = [node for node in module.body if not isinstance(node,ast.ClassDef) or
                       node.name not in ('ModalCurrentSubspaces','NormalRealityCoverage')]
        return ast.dump(module,include_attributes=False)
    if nonmodal_source(base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py')!=nonmodal_source(Path(engine.__file__)):
        raise ValueError('source-control constructors changed outside modal certificate repair')
    report = checkpoint['controls'][label]
    if report['sourceReplacements']!={str(s):str(v) for s,v in replacements.items()} or \
       report['controlBindings']!={str(s):str(v) for s,v in extra_bindings.items()}:
        raise ValueError(('source-control cache specification mismatch',label))
    payload = base/label/'objects.pickle'
    if digest(payload)!=report['objectsSha256']:
        raise ValueError(('source-control cache payload mismatch',label))
    objects = pickle.loads(payload.read_bytes())
    frozen = base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    reuse_modes = ('NORMAL_REALITY_COVERAGE' in objects.get('modes',{}) and all(
        source_node(frozen.read_text(),name)==source_node(Path(engine.__file__).read_text(),name)
        for name in ('ModalCurrentSubspaces','NormalRealityCoverage')))
    if reuse_modes:
        transcript = base/'full.out'
        if digest(transcript)!=checkpoint['artifacts']['full.out']['sha256']:
            raise ValueError('cached modal extraction transcript pin')
        prefix = 'PY_S11CD_CURRENT_SOURCE_CONTROL_'+label+'_MODAL_SUBSPACE_QUADRATIC_EXTRACTION_'
        extraction = {}
        for line in decoded_lines(transcript):
            tag,_,body = line.rstrip('\n').partition(': ')
            if tag.startswith(prefix):extraction[tag.removeprefix(prefix)] = _restore(body)
        cls = next(node for node in ast.parse(frozen.read_text()).body if isinstance(node,ast.ClassDef) and node.name=='ModalCurrentSubspaces')
        prepare = next(node for node in cls.body if isinstance(node,ast.FunctionDef) and node.name=='prepare')
        expected = next(ast.literal_eval(node.value) for node in prepare.body if isinstance(node,ast.Assign) and
                        any(isinstance(target,ast.Name) and target.id=='sources' for target in node.targets))
        if set(extraction)!=set(expected.values()):raise ValueError('cached quadratic extraction coverage')
        objects['coefficientResiduals'] = extraction
    engine.PHYSICAL_METADATA.dimensions.known.update(objects['knownDimensions'])
    return objects, {'checkpointSha256':digest(args.reuse_checkpoint),'control':label,
        'objectsSha256':digest(payload),'frozenEngineSha256':checkpoint['provenance']['engineSha256'],
        'reusedChecks':True,'reusedModalPacket':reuse_modes,
        'reused':'source current, exact balance checks and isolated spectrum; modal reuse requires identical modal/certificate source and pinned extraction transcript'}


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--current-manifest', type=Path, required=True)
    parser.add_argument('--current-checkpoint', type=Path, required=True)
    parser.add_argument('--run-directory', type=Path, required=True)
    parser.add_argument('--reuse-checkpoint', type=Path)
    parser.add_argument('--control', choices=('massTransferOff','bulkDensityChange','densityGradientOff','all'), default='all')
    parser.add_argument('--emit', action='store_true')
    parser.add_argument('--modes', action='store_true')
    args = parser.parse_args()
    args.end = 'REFERENCE'
    args.run_directory.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    def progress(record):
        record['elapsedSeconds'] = time.monotonic()-started
        with (args.run_directory/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps(record)+'\n')
        if not args.emit:
            print(json.dumps(record), flush=True)
    baseline, base_result, inputs, provenance = load_current(args)
    provenance['instrumentSha256'] = digest(Path(__file__))
    source_manifest = json.loads(args.manifest.read_text())
    source_energy = pickle.loads((Path(source_manifest['run_directory'])/'symbols/REFERENCE_LAB_HELD_RHO4_CONSTANT.pickle').read_bytes())[5]
    factory = engine.CurrentSourceControls(baseline, source_energy)
    d, r = engine.PHYSICAL_METADATA.dimensions, baseline.r
    density = sp.Symbol('s11cdControlBulkDensity', positive=True)
    d.known[density] = d.measure(r.symbols['rho_m'])
    specifications = {
        'massTransferOff': ({r.symbols[name]:sp.S.Zero for name in ('Lambda_A_0','Lambda_V_0')}, {}),
        'bulkDensityChange': ({r.symbols['rho_m']:density}, {density:2*inputs.parameters['rho_m']}),
        'densityGradientOff': ({r.symbols[name]:sp.S.Zero for name in ('kappa_theta','kappa_theta_W')}, {}),
    }
    selected = specifications if args.control == 'all' else {args.control:specifications[args.control]}
    reports = {}
    for label,(replacements, extra_bindings) in selected.items():
        progress({'control':label,'stage':'source_construction'})
        builder, controlled_energy = factory.build(replacements)
        reused, reuse_provenance = reuse_control(args,label,provenance,replacements,extra_bindings)
        if reused is not None:
            builder.construct = lambda anchoring,end:reused['pairing']
            builder.balance.construct = lambda anchoring,end:reused['slab']
            builder.c.construct = lambda anchoring,end:reused['conservative']
            builder.acoustic.construct = lambda anchoring,end:reused['acoustic']
        m = builder.modes
        algebraic, relation, _ = m.analytic(builder.acoustic.strong)
        parameters = {**inputs.parameters, **{symbol.name:value for symbol,value in extra_bindings.items()}}
        binding = {symbol:parameters[symbol.name] for symbol in algebraic.free_symbols | relation.free_symbols
                   if symbol not in (m.k,m.q,m.eta,m.sigma,r.omega)}
        binding.update({m.eta:sp.S.Zero,m.sigma:sp.S.Zero})
        destination = args.run_directory/label
        destination.mkdir(exist_ok=True)
        result = builder.construct('LAB_HELD', None)
        slab = builder.balance.construct('LAB_HELD',None)
        conservative = builder.c.construct('LAB_HELD',None)
        progress({'control':label,'stage':'current_constructed'})
        checks = dict(reused['checks']) if reused is not None else builder.split_balance_checks(
            result, binding, lambda record:progress({'control':label,**record}))
        same_frequency = dict.fromkeys(builder.frequencies, r.omega)
        checks['SLAB_CURRENT_DIAGONAL_JOIN_RESIDUAL'] = (result['SLAB_CURRENT_MATRIX'].xreplace(same_frequency)-
            slab['SLAB_CURRENT_MATRIX']).applyfunc(sp.expand)
        checks['REDUCED_MASS_FACE_JOIN_RESIDUAL'] = builder.acoustic.construct('LAB_HELD',None)['REDUCED_MASS_FACE_JOIN_RESIDUAL']
        checks['REDUCED_MECHANICAL_FACE_JOIN_RESIDUAL'] = builder.acoustic.construct('LAB_HELD',None)['REDUCED_MECHANICAL_FACE_JOIN_RESIDUAL']
        progress({'control':label,'stage':'checks_constructed'})
        prefix = 'CURRENT_SOURCE_CONTROL_'+label
        def output(name,value,units=None,heavy=False):
            if not args.emit:
                return
            body = engine.cas(value)
            if units is None:
                units = {path:d.measure(value) for path,value in engine.leaves(body)
                         if not isinstance(value,sp.core.symbol.Str)}
            payload = engine.carrier_fingerprint(body) if heavy else body
            engine.emit(prefix+'_'+name,payload)
            engine.emit('METADATA_'+prefix+'_'+name,m.numeric_metadata(body,lambda path:units[path]))
        def uniform_units(value,unit):
            return {path:unit for path,v in engine.leaves(engine.cas(value)) if not isinstance(v,sp.core.symbol.Str)}
        output('PROVENANCE',provenance,uniform_units(provenance,d.zero))
        if reuse_provenance is not None:
            output('CACHE_PROVENANCE',reuse_provenance,uniform_units(reuse_provenance,d.zero))
        for name,value in (('SOURCE_REPLACEMENTS',{str(s):v for s,v in replacements.items()}),
                           ('CONTROL_PARAMETER_BINDINGS',{str(s):v for s,v in extra_bindings.items()})):
            symbols = replacements if name=='SOURCE_REPLACEMENTS' else extra_bindings
            output(name,value,{(str(s),):d.measure(s) for s in symbols})
        output('MATERIAL_BINDINGS',{str(s):v for s,v in binding.items()},
               {(str(s),):d.measure(s) for s in binding})
        source_unit = d.measure(baseline.c.construct('LAB_HELD',None)['UNIFORM_SOURCE_ENERGY'])
        for name,value in (('BASE_SOURCE_ENERGY',baseline.c.construct('LAB_HELD',None)['UNIFORM_SOURCE_ENERGY']),
                           ('CONTROL_SOURCE_ENERGY',conservative['UNIFORM_SOURCE_ENERGY'])):
            output(name,value,{():source_unit},heavy=True)
        field_units = baseline.c.field_units
        row_units = [d.measure(value) for value in baseline.residual_amplitudes[0]]
        pencil_units = {(5*i+j,):tuple(a-b for a,b in zip(row_units[i],field_units[j])) for i in range(5) for j in range(5)}
        output('BASE_REDUCED_PENCIL',baseline.acoustic.strong,pencil_units,heavy=True)
        output('CONTROL_REDUCED_PENCIL',builder.acoustic.strong,pencil_units,heavy=True)
        selected_keys = ('CLOSED_PENCIL_LEGS','STORED_ENERGY_MATRIX','KINETIC_ENERGY_MATRIX','SLAB_ENERGY_MATRIX',
            'SLAB_CURRENT_MATRIX','BULK_ENERGY_DENSITY_MATRIX','BULK_NORMAL_CURRENT_DENSITY_MATRIX',
            'BULK_DEPTH_CURRENT_MATRIX','PORT_POWER_MATRIX','INTERFACE_POWER_MATRIX','SOURCE_POWER_MATRIX',
            'PLUS_ROW_POWER_MAP','MINUS_ROW_POWER_MAP','FINITE_TOTAL_ENERGY_MATRIX','FINITE_TOTAL_CURRENT_MATRIX')
        comparisons = {}
        for key in selected_keys:
            base_value, value = base_result[key], result[key]
            units = baseline.output_units(key,base_value)
            output('BASE_'+key,base_value,units,heavy=True)
            output('CONTROL_'+key,value,units,heavy=True)
            # Rebuilt-minus-baseline is a sensitivity operand, not a zero-target identity.
            difference = tuple(v-b for v,b in zip(value,base_value)) if isinstance(value,tuple) else value-base_value
            output('CHANGE_'+key,difference,units,heavy=True)
            comparisons[key] = {'base':engine.carrier_fingerprint(engine.cas(base_value)),
                'control':engine.carrier_fingerprint(engine.cas(value)),
                'change':engine.carrier_fingerprint(engine.cas(difference))}
        current_units = baseline.output_units('SLAB_CURRENT_MATRIX',base_result['SLAB_CURRENT_MATRIX'])
        output('MASS_RATE_BOUNDARY_CORRECTION_MATRIX',slab['MASS_RATE_CORRECTION_MATRIX'],current_units,heavy=True)
        material = {r.symbols[name]:value for name,value in inputs.parameters.items() if name in r.symbols}
        material.update(extra_bindings)
        material.update({m.eta:sp.S.Zero,m.sigma:sp.S.Zero})
        sample = {**dict.fromkeys(builder.frequencies,inputs.parameters['omega']),
                  **dict.fromkeys(builder.c.leg_momenta,inputs.parameters[r.tangents[0].name]),
                  builder.acoustic.height:inputs.parameters['W_0']}
        qleft,qright = builder.acoustic.qlegs
        wave = base_result['ACOUSTIC_WAVE_ROWS'][0].xreplace({**material,**sample})
        roots = sp.solve(wave,qright)
        outgoing = [q for q in roots if sp.im(q).is_positive]
        if not outgoing:
            outgoing = [q for q in roots if sp.im(q)==0 and sp.re(q).is_positive]
        if len(outgoing)!=1:
            raise ValueError('source sensitivity point needs one regular outgoing bulk lift')
        sample.update({qright:outgoing[0],qleft:sp.conjugate(outgoing[0])})
        output('SENSITIVITY_POINT',{str(s):v for s,v in sample.items()},
               {(str(s),):d.measure(s) for s in sample})
        output('SENSITIVITY_WAVE_RESIDUAL',sp.simplify(wave.subs(qright,outgoing[0])),
               {():tuple(2*v for v in d.measure(r.omega))})
        bound_comparisons, sensitivity = {}, {}
        epsilon = r.symbols['epsilon_shape']
        for key in ('SLAB_ENERGY_MATRIX','SLAB_CURRENT_MATRIX','FINITE_TOTAL_CURRENT_MATRIX',
                    'INTERFACE_POWER_MATRIX','SOURCE_POWER_MATRIX'):
            left = base_result[key].xreplace({**material,**sample}).evalf(17)
            right = result[key].xreplace({**material,**sample}).evalf(17)
            change = (right-left).applyfunc(sp.expand)
            bound_comparisons[key] = {'BASE':left,'CONTROL':right,'CHANGE':change}
            units = baseline.output_units(key,base_result[key])
            for side,value in bound_comparisons[key].items():
                output('BOUND_'+side+'_'+key,value,units,heavy=side!='CHANGE')
            coefficients = [complex(v.subs(epsilon,1)) for v in change]
            sensitivity[key] = {'nonzeroScalars':sum(v!=0 for v in coefficients),
                                'maximumAbsoluteCoefficient':max(map(abs,coefficients))}
        residual_inventory = {}
        for key,value in checks.items():
            if key in ('REDUCED_MASS_FACE_JOIN_RESIDUAL','REDUCED_MECHANICAL_FACE_JOIN_RESIDUAL'):
                index = 3 if key.startswith('REDUCED_MASS') else 4
                units = {(j,):pencil_units[(5*index+j,)] for j in range(5)}
            elif key=='SLAB_CURRENT_DIAGONAL_JOIN_RESIDUAL':
                units = current_units
            elif key in ('EQUAL_DEPTH_WAVE_ROWS','EQUAL_DEPTH_WAVE_ELIMINANT'):
                units = uniform_units(value,tuple(2*v for v in d.measure(r.omega)))
            else:
                units = builder.check_output_units(key,value,checks)
            output('CHECK_'+key,value,units,heavy=not key.endswith('_RESIDUAL'))
            if key.endswith('_RESIDUAL'):
                entries = [v for _,v in engine.leaves(engine.cas(value))]
                residual_inventory[key] = {'scalars':len(entries),'nonzeroScalars':sum(v!=0 for v in entries)}
        record = {'provenance':provenance,'control':label,'residualInventory':residual_inventory,
                  'sourceReplacements':{str(s):str(v) for s,v in replacements.items()},
                  'controlBindings':{str(s):str(v) for s,v in extra_bindings.items()},
                  'sensitivity':sensitivity,'cacheReuse':reuse_provenance if reuse_provenance is not None else 'FRESH_SOURCE_CONSTRUCTION'}
        objects = {'pairing':result,'checks':checks,'slab':slab,'conservative':conservative,
                   'acoustic':builder.acoustic.construct('LAB_HELD',None),'comparisons':comparisons,
                   'boundComparisons':bound_comparisons,'sensitivityPoint':sample}
        if args.modes:
            physical_binding = {**binding,r.omega:inputs.parameters['omega']}
            spectrum = reused['spectrum'] if reused is not None else factory.spectrum(builder,physical_binding)
            if reused is not None:
                actual = algebraic.xreplace(physical_binding).applyfunc(sp.cancel)
                if actual!=spectrum['PHYSICAL_PENCIL'] or relation.xreplace(physical_binding)!=spectrum['RADICAL_RELATION']:
                    raise ValueError(('source-control cached spectrum pencil mismatch',label))
            objects['spectrum'] = spectrum
            output('SPECTRUM_DEFINED',spectrum['DEFINED'],{():d.zero})
            if spectrum['DEFINED']:
                for key in ('ELIMINATION_OPERANDS','EXCEPTION_POLYNOMIALS','EXCEPTION_DEGREES','COVERAGE','RECORDS','SHEET_PATHS'):
                    value = spectrum[key]
                    # Polynomial coefficients, row-denominator and SVD diagnostics
                    # use the reference-unit coordinate frame. Root/path units are restored.
                    units = uniform_units(value,d.zero)
                    for path in units:
                        if path[-1] in ('K','START_K','END_K','MINIMUM_BRANCH_POINT_DISTANCE','GEOMETRIC_RESOLUTION'):
                            units[path] = d.measure(m.k)
                        elif path[-1] in ('Q','CENTER','RADIUS','PRECISION_REFINEMENT_DIFFERENCE','SEED_Q','END_Q',
                                        'REFINEMENT_DIFFERENCE','SHEET_DIFFERENCE','OPPOSITE_SHEET_DIFFERENCE'):
                            units[path] = d.measure(m.q)
                        elif path[-1] in ('RADICAL_RESIDUAL','MAXIMUM_RADICAL_RESIDUAL'):
                            units[path] = tuple(2*v for v in d.measure(m.q))
                        elif 'BRANCH_POINTS' in path:
                            units[path] = d.measure(m.k)
                    output('SPECTRUM_'+key,value,units,heavy=key in ('ELIMINATION_OPERANDS','EXCEPTION_POLYNOMIALS'))
                regular = (spectrum['COVERAGE']['FINITE_POLYNOMIAL_ROOT_COVERAGE'] and
                    all(degree==0 for degree in spectrum['EXCEPTION_DEGREES'].values()) and
                    all(v['FINITE_PENCIL'] and v.get('NULLITY',0)>0 for v in spectrum['RECORDS']))
                output('REGULAR_MODE_CONTRACTION_DOMAIN',regular,{():d.zero})
                record['spectrum'] = {'candidateCount':len(spectrum['RECORDS']),
                    'rootCoverage':spectrum['COVERAGE']['FINITE_POLYNOMIAL_ROOT_COVERAGE'],
                    'exceptionDegrees':spectrum['EXCEPTION_DEGREES'],'regularContractionDomain':regular}
                if regular:
                    modal = engine.ModalCurrentSubspaces(builder,result,physical_binding)
                    coverage = {str(k):v for k,v in engine.cas(spectrum['COVERAGE'])}
                    if reused is not None and 'coefficientResiduals' in reused:
                        modes = reused['modes']
                        modal.coefficient_residuals = reused['coefficientResiduals']
                    else:
                        modes = modal.construct(spectrum['RECORDS'],coverage,float(inputs.parameters['omega']),
                            float(inputs.parameters['W_0']),lambda item:progress({'control':label,**item}))
                    objects['modes'] = modes
                    if args.emit:
                        modal.emit(modes,record,prefix=prefix+'_MODAL_SUBSPACE')
                    record['spectrum']['normalizationCount'] = sum(v.get('PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED',False) for v in modes['RECORDS'])
                    record['spectrum']['nullityCounts'] = dict(engine.Counter(str(v['NULLITY']) for v in modes['RECORDS']))
                    record['spectrum']['residualNormMaxima'] = {name:max(v['RESIDUAL_NORMS'].get(name,0) for v in modes['RECORDS'])
                        for name in sorted(set().union(*(v['RESIDUAL_NORMS'] for v in modes['RECORDS'])))}
                    import numpy as np
                    correction = slab['MASS_RATE_CORRECTION_MATRIX'].applyfunc(
                        lambda v:dict(engine.polynomial_terms(v,(epsilon,))).get((2,),sp.S.Zero))
                    correction = correction.xreplace(material)
                    corrections = []
                    current_unit = tuple(a+b+c for a,b,c in zip(
                        d.measure(conservative['TANGENTIAL_ENERGY_REDUCTION']),d.measure(r.omega),d.measure(r.z)))
                    for mode in modes['RECORDS']:
                        k = complex(mode['K'])
                        matrix = np.asarray(correction.subs(dict(zip(builder.c.leg_momenta,(k.conjugate(),k)))),dtype=complex)
                        right = mode['FORMS']['RIGHT']
                        contracted = right.conj().T@matrix@right
                        body = epsilon**2*sp.ImmutableMatrix(contracted.shape[0],contracted.shape[1],
                            [m.number(v) for v in contracted.ravel()])
                        output('MODE_'+str(mode['INDEX'])+'_MASS_RATE_BOUNDARY_CORRECTION',body,
                               uniform_units(body,current_unit),heavy=True)
                        corrections.append({'index':mode['INDEX'],'matrix':contracted,'norm':float(np.linalg.norm(contracted))})
                    objects['modalBoundaryCorrections'] = corrections
                    record['spectrum']['maximumBoundaryCorrectionNorm'] = max(v['norm'] for v in corrections)
        objects['knownDimensions'] = dict(d.known)
        payload = destination/'objects.pickle'
        payload.write_bytes(pickle.dumps(objects,protocol=5))
        record.update({'objectsSha256':digest(payload),'elapsedSeconds':time.monotonic()-started,
                       'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
        (destination/'checks.json').write_text(json.dumps(record,indent=2)+'\n')
        reports[label] = record
        progress({'control':label,'stage':'completed','record':record})
    output('DIMENSION_CONSTRAINTS',tuple(d.constraints),uniform_units(tuple(d.constraints),d.zero))
    summary = {'provenance':provenance,'controls':reports,'wallSeconds':time.monotonic()-started,
               'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               'dimensionConstraints':[sp.srepr(v) for v in d.constraints]}
    (args.run_directory/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    if args.emit:
        index = engine.emission_index(engine.EMISSION_LINES)
        engine.emit('CURRENT_SOURCE_CONTROL_EMISSION_LINES',index)
        engine.emit('METADATA_CURRENT_SOURCE_CONTROL_EMISSION_LINES',m.numeric_metadata(engine.cas(index),lambda path:d.zero))
    progress({'stage':'completed'})
    if any(v['nonzeroScalars'] for report in reports.values() for v in report['residualInventory'].values()):
        raise ValueError('source-control reconstruction residual; see emitted operands')


if __name__ == '__main__':
    run()
