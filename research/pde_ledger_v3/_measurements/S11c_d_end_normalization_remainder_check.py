#!/usr/bin/env python3
"""Trace saved pairing truncation remainders through every native mode space.

This computes a second route from the independently validated two-frequency
pairing certificate. It never changes the current, pencil, or cached residual.
"""
import argparse
from contextlib import redirect_stdout
import json
from pathlib import Path
import pickle
from types import SimpleNamespace

import numpy as np
import sympy as sp

import S11c_d_end_normalization_check as runner


def compute(base, pairing, builder, modal, adjoint, binding):
    engine = runner.engine
    args = json.loads((base/'arguments.json').read_text())
    checkpoint_path = Path(args['pairing_checkpoint'])
    checkpoint = json.loads(checkpoint_path.read_text())
    source = runner.checked_artifacts(checkpoint)
    packet, _ = pickle.loads((source/'complete.pickle').read_bytes())
    if runner.digest(source/'complete.pickle') != modal_provenance(base)['pairingCacheSha256']:
        raise ValueError('normalization/remainder pairing cache join')
    epsilon = builder.epsilon
    eta, sigma = (pairing.r.symbols[n] for n in ('eta_bg','sigma_W'))
    qleft, qright = pairing.acoustic.qlegs
    expressions, algebraic_checks = {}, {}
    for label, name in (('FINITE','FINITE_COMPOSED_BALANCE_RESIDUAL'),
                        ('EQUAL','EQUAL_DEPTH_COMPOSED_BALANCE_RESIDUAL')):
        raw, kept, rest = (packet[group][name] for group in ('checks','retained','remainders'))
        algebraic_checks[label+'_DECOMPOSITION'] = (raw-kept-rest).applyfunc(sp.expand)
        expressions[label] = rest.applyfunc(lambda v:dict(engine.polynomial_terms(v,(epsilon,))).get((2,),sp.S.Zero))
        algebraic_checks[label+'_EPSILON_COEFFICIENT'] = (rest-epsilon**2*expressions[label]).applyfunc(sp.expand)
        algebraic_checks[label+'_RETAINED'] = kept
        algebraic_checks[label+'_REPROJECTED'] = expressions[label].applyfunc(pairing.c.retained)
    # Differentiate the actual discarded operand, before the equal-depth limit.
    # The pairing packet separately certifies the regular-sheet wave quotient
    # and its tangent, so this is not differentiation of a substituted zero.
    for label, variable in (('NORMAL',pairing.c.leg_momenta[1]),('FREQUENCY',pairing.frequencies[1])):
        transport = packet['checks'][label+'_RADICAL_TRANSPORT'].xreplace(packet['bindings'])
        expressions[label] = (expressions['FINITE'].diff(variable)+transport*expressions['FINITE'].diff(qright)).applyfunc(sp.cancel)
        algebraic_checks[label+'_REPROJECTED'] = expressions[label].applyfunc(pairing.c.retained)
        expressions[label+'_EQUAL'] = expressions[label].applyfunc(lambda v:sp.limit(v,qleft,qright))
    algebraic_checks['EQUAL_BRANCH_JOIN'] = (expressions['FINITE'].applyfunc(lambda v:sp.limit(v,qleft,qright))-
                                            expressions['EQUAL']).applyfunc(sp.cancel)
    coefficients = {}
    for name, matrix in expressions.items():
        coefficients[name] = sorted({grade for value in matrix for grade,term in engine.polynomial_terms(value,(eta,sigma)) if term!=0})
    evaluate = {name:sp.lambdify(builder.variables,matrix.xreplace(binding),'numpy',cse=True)
                for name,matrix in expressions.items()}
    scalars = modal['SCALAR_OPERANDS']
    weights = {label:sp.lambdify(builder.variables,scalars[label+suffix].xreplace(binding),'numpy',cse=True)
               for label,suffix in (('NORMAL','_NORMAL_WEIGHT'),('FREQUENCY','_TIME_WEIGHT'))}
    records, tensors = [], []
    maxima = {}
    for mode, dual in zip(modal['RECORDS'],adjoint['RECORDS']):
        if mode['INDEX']!=dual['INDEX'] or mode['NULLITY']!=dual['NULLITY']:
            raise ValueError('remainder modal/adjoint subspace join')
        right = mode['FORMS']['RIGHT']
        k,q,w,h = (complex(mode[key]) for key in ('K','PHYSICAL_Q','OMEGA','DEPTH_CUTOFF'))
        point = (w,w,k.conjugate(),k,q.conjugate(),q,h)
        equal = mode['EQUAL_DEPTH_BRANCH']
        values = {name:np.asarray(evaluate[name+('_EQUAL' if equal and name!='FINITE' else '')](*point),dtype=complex)
                  for name in ('NORMAL','FREQUENCY')}
        values['FINITE'] = np.asarray(evaluate['EQUAL' if equal else 'FINITE'](*point),dtype=complex)
        expected = {'MODAL_FINITE_BALANCE':right.conj().T@values['FINITE']@right}
        actual = {'MODAL_'+name:value for name,value in mode['RESIDUALS'].items()}
        actual.update({'ADJOINT_'+item['NAME']:item['VALUE'] for item in dual['ITEMS'] if item['GROUP']=='RESIDUALS'})
        field = next((item['VALUE'] for item in dual['ITEMS'] if item['NAME']=='ADJOINT_FIELD'),None)
        for label in ('NORMAL','FREQUENCY'):
            weight = complex(weights[label](*point))
            derivative = right.conj().T@values[label]@right
            expected['MODAL_'+label+'_BALANCE_DERIVATIVE'] = derivative
            expected['MODAL_'+label+'_CURRENT_ENERGY_RECONSTRUCTION'] = -derivative/weight
            if field is not None:
                mixed = field.conj().T@values[label]@right
                expected['ADJOINT_'+label+'_MIXED_BALANCE_DERIVATIVE'] = mixed
                expected['ADJOINT_'+label+'_MIXED_CURRENT_ENERGY_RECONSTRUCTION'] = -mixed/weight
        norms = {}
        for name,value in expected.items():
            difference = actual[name]-value
            if not all(np.isfinite(v).all() for v in (actual[name],value,difference)):
                raise ValueError(('nonfinite remainder contraction',mode['INDEX'],name))
            norm = float(np.linalg.norm(difference));norms[name]=norm
            maxima[name]=max(maxima.get(name,0.),norm)
            tensors.append({'index':mode['INDEX'],'name':name,'actual':actual[name],
                            'remainder':value,'difference':difference})
        records.append({'index':mode['INDEX'],'rootDiskIndex':mode['ROOT_DISK_INDEX'],
            'normalLiftSign':mode['NORMAL_LIFT_SIGN'],'nullity':mode['NULLITY'],
            'rightRank':int(np.linalg.matrix_rank(right,tol=1e-9)),
            'adjointRank':None if field is None else int(np.linalg.matrix_rank(field,tol=1e-9)),
            'residualMinusRemainderNorms':norms})
    summary = {'pairingCheckpointSha256':runner.digest(checkpoint_path),
        'pairingCacheSha256':runner.digest(source/'complete.pickle'),
        'sourceSha256':runner.digest(Path(__file__)),
        'coefficientGradesEtaSigma':coefficients,
        'algebraicChecks':{name:{'scalars':len(value),'nonzeroScalars':sum(v!=0 for v in value)}
                           for name,value in algebraic_checks.items()},
        'records':records,'residualMinusRemainderNormMaxima':maxima}
    return {'summary':summary,'expressions':expressions,'algebraicChecks':algebraic_checks,'tensors':tensors}


def modal_provenance(base):
    return json.loads((base/'checks.json').read_text())['provenance']


def emit_evidence(base, builder, computed):
    """Emit both contraction operands and their difference with restored units."""
    engine=runner.engine
    summary=json.loads((base/'checks.json').read_text())
    prefix='END_REMAINDER_'+summary['end']+'_'+runner.CASE
    packets=[]
    def put(name, value, unit=lambda path:(0,0,0), carrier=False):
        body=engine.cas(value)
        payload=engine.carrier_fingerprint(body) if carrier else body
        metadata=builder.modes.numeric_metadata(body,unit)
        packets.extend(((prefix+'_'+name,payload),('METADATA_'+prefix+'_'+name,metadata)))
    put('PROVENANCE',{'END':summary['end'],'CASE':runner.CASE,
        'NORMALIZATION_MODAL_SHA256':summary['objects']['modal.pickle'],
        'NORMALIZATION_ADJOINT_SHA256':summary['objects']['adjoint.pickle'],
        'PAIRING_CACHE_SHA256':computed['summary']['pairingCacheSha256'],
        'SOURCE_SHA256':computed['summary']['sourceSha256'],
        'FINITE_CONTRAST_EVALUATION':True,'CONTINUUM_REEXPANSION_PERFORMED':False})
    def matrix_unit(quantity):
        return lambda path:tuple(a-b-c for a,b,c in zip(quantity,
            builder.field_units[path[0]//5],builder.field_units[path[0]%5]))
    for name,matrix in computed['expressions'].items():
        unit=builder.current_unit if name.startswith('NORMAL') else builder.energy_unit if name.startswith('FREQUENCY') else builder.power_unit
        put(name+'_DISCARDED_MATRIX',builder.epsilon**2*matrix,matrix_unit(unit),True)
    for name,matrix in computed['algebraicChecks'].items():
        unit=builder.current_unit if name.startswith('NORMAL') else builder.energy_unit if name.startswith('FREQUENCY') else builder.power_unit
        coefficient=name.endswith('_REPROJECTED') or name=='EQUAL_BRANCH_JOIN'
        put(name+'_RESIDUAL',builder.epsilon**2*matrix if coefficient else matrix,matrix_unit(unit))
    for record in computed['summary']['records']:
        put('MODE_'+str(record['index'])+'_SUBSPACE',{key:value for key,value in record.items()
            if key!='residualMinusRemainderNorms' and value is not None})
    for tensor in computed['tensors']:
        name=tensor['name']
        unit=(tuple(a+b for a,b in zip(builder.length_unit,builder.frequency_unit)) if 'NORMAL' in name else (0,0,0)) if name.startswith('ADJOINT_') else (
            builder.current_unit if 'NORMAL' in name else builder.energy_unit if 'FREQUENCY' in name else builder.power_unit)
        tag='MODE_'+str(tensor['index'])+'_'+name
        for label in ('actual','remainder','difference'):
            array=tensor[label]
            body=builder.epsilon**2*sp.ImmutableMatrix(*array.shape,[builder.modes.number(v) for v in array.ravel()])
            put(tag+'_'+label.upper(),body,lambda path,unit=unit:unit)
    names=[name for name,_ in packets if not name.startswith('METADATA_')]
    keys={name:'s11cd'+''.join(word.title() for word in name.split('_')) for name in names}
    put('WRITE_KEYS',keys)
    if len(set(keys.values()))!=len(keys) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('remainder evidence key collision')
    path=base/'remainders.out'
    if path.exists():raise ValueError('remainder evidence target exists')
    with path.open('x') as stream,redirect_stdout(stream):
        for name,value in packets:engine.emit(name,value)
    decoded={tag:runner._restore(body) for line in runner.decoded_lines(path)
             for tag,sep,body in (line.rstrip('\n').partition(': '),) if sep}
    expected={'PY_S11CD_'+name:value for name,value in packets}
    if decoded!=expected or len(decoded)!=len(packets):raise ValueError('remainder evidence serialization/census')
    paths=sum(len(engine.named(group,'PATHS')) for name,body in packets if name.startswith('METADATA_') for group in body)
    return {'path':str(path),'bytes':path.stat().st_size,'sha256':runner.digest(path),
            'tagCount':len(packets),'metadataPaths':paths,'writeKeyCount':len(keys)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    args=parser.parse_args();base=args.run_directory.resolve()
    args.destination.resolve().relative_to(runner.ROOT.parents[1]/'_scratch/s11c')
    args.destination.mkdir(parents=True,exist_ok=False)
    values=json.loads((base/'arguments.json').read_text())
    constructor_args=SimpleNamespace(**{k:Path(v) if v is not None and k!='end' else v for k,v in values.items()})
    pairing,current,inputs,provenance=runner.load_pairing(constructor_args)
    native,coverage,binding,physical=runner.load_native(constructor_args,pairing,inputs,provenance)
    modal,known=pickle.loads((base/'modal.pickle').read_bytes())
    adjoint,dual_known=pickle.loads((base/'adjoint.pickle').read_bytes())
    runner.engine.PHYSICAL_METADATA.dimensions.known.update(known)
    runner.engine.PHYSICAL_METADATA.dimensions.known.update(dual_known)
    builder=runner.engine.ModalCurrentSubspaces(pairing,current,binding)
    result=compute(base,pairing,builder,modal,adjoint,binding)
    (args.destination/'objects.pickle').write_bytes(pickle.dumps(result,protocol=5))
    (args.destination/'checks.json').write_text(json.dumps(result['summary'],indent=2)+'\n')
    print(json.dumps(result['summary'],indent=2),flush=True)


if __name__=='__main__':main()
