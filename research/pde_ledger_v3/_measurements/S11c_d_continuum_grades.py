#!/usr/bin/env python3
"""Extract independent continuum grades from saved native coefficient operands."""
import argparse
import ast
import contextlib
from functools import lru_cache
import inspect
import json
from pathlib import Path
import resource
import signal
import time

import sympy as sp
import S11c_d_finite_scattering as f
from S11c_d_output_codec import PayloadEncoder, decoded_lines, restore_emission_index
from ledger_fold import _restore

engine=f.engine
PLAN=f.M/'S11c_d_continuum_grade_plan.md'
PREFIX='CONTINUUM_OPERATOR_GRADES_LAB_HELD_RHO4_CONSTANT'


def packet(checkpoint_name,filename):
    checkpoint=f.M/checkpoint_name
    record=json.loads(checkpoint.read_text());base=Path(record['runDirectory'])
    path=base/filename
    f.require(f.digest(path)==record['artifacts'][filename]['sha256'],('accepted packet',filename))
    publications=record.get('publications',{'main':record.get('publication')})
    for publication in publications.values():
        f.require(publication is not None and (f.ROOT/publication['path']).is_symlink() and
                  f.digest(f.ROOT/publication['path'])==publication['sha256'],'accepted annex publication')
    for name,sha in record['sourceFiles'].items():
        f.require(f.digest(base/'source'/name)==sha,('frozen source',name))
    return f.unpickle(path),record,path,checkpoint


def method(source,cls,name):
    node=next(n for n in ast.parse(source).body if isinstance(n,ast.ClassDef) and n.name==cls)
    return ast.dump(next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name==name))


def load(base):
    assembly,a,ap,ac=packet('S11c_d_reduced_action_assembly_checkpoint.json','assembly.pickle')
    fourier,b,bp,bc=packet('S11c_d_source_fourier_factorization_checkpoint.json','source-factorization.pickle')
    source,c,cp,cc=packet('S11c_d_reduced_action_source_checkpoint.json','reduced-action.pickle')
    actions_path=Path(c['runDirectory'])/'actions.pickle'
    f.require(f.digest(actions_path)==c['artifacts']['actions.pickle']['sha256'],'accepted action packet')
    actions=f.unpickle(actions_path)
    f.require(assembly['provenance']==a['provenance'] and assembly['provenance']['ACTIONS_SHA256']==f.digest(actions_path)
        and actions['reducedActionSha256']==f.digest(cp),'assembly/action/reduction provenance')
    f.require(f.digest(bp)==b['sourcePacketSha256'],'accepted factorization source identity')
    f.require(len(fourier['result']['ROWS'])==len(assembly['result']['NONLOCAL_INTEGRALS'])==80,'all native integrals')
    for row,original in zip(fourier['result']['ROWS'],assembly['result']['NONLOCAL_INTEGRALS']):
        f.require(row['ORIGINAL']==original,'native integral/factor row join')
    engine_name=str(engine.HERE.relative_to(f.ROOT))
    frozen=(Path(b['sourceRunDirectory'])/'source'/engine_name).read_text()
    helper_joins={name:method(frozen,'PhysicalMetadata',name)==method(engine.HERE.read_text(),'PhysicalMetadata',name)
                  for name in ('coefficients','record')}
    f.require(all(helper_joins.values()),'unchanged native grade/metadata helpers')
    stable=('scripts/S11c_b_exports.py','scripts/S11c_c1_exports.py','scripts/S11c_c2_exports.py',
            'directives/S11c_d_SHARED_PHYSICS.md','directives/S11b_SHARED_PHYSICS.md','scripts/ledger_fold.py')
    for name in stable:
        f.require(a['sourceFiles'][name]==b['sourceFiles'][name]==c['sourceFiles'][name]==f.digest(f.ROOT/name),('consumed physics identity',name))
    restore=f.prior.domain.momentum.source.native.source.restore_context
    r,dimensions=restore(source);dimensions.__dict__.update(fourier['dimensionState'])
    field_units=[dimensions.known[v] for v in actions['fields']]
    equation_units=[actions['columnUnits'][(0,i)] for i in range(5)]
    paths=(Path(__file__),PLAN,ac,bc,cc,engine.HERE,Path(inspect.getfile(restore)),f.ACCEPTANCE,
           f.M/'S11c_d_finite_scattering_domain_checkpoint.json',f.M/'S11c_d_variable_profile_development_input.json',
           f.ROOT/'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md',*(f.ROOT/n for n in stable))
    pins={str(p.resolve().relative_to(f.ROOT)):f.digest(p) for p in paths}
    for name in pins:
        destination=base/'source'/name;destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_bytes((f.ROOT/name).read_bytes())
    inputs={str(p):f.digest(p) for p in (ap,bp,cp,actions_path)}
    f.save(base/'inputs.json',{'sourceFiles':pins,'inputPackets':inputs,'nativeHelperJoins':helper_joins,
        'fieldUnits':[list(map(str,u)) for u in field_units],'equationUnits':[list(map(str,u)) for u in equation_units],
        'gradeGenerators':[str(v) for v in engine.PHYSICAL_METADATA.generators],
        'scope':'Coefficient operators before independent background-grade response construction; no new numerical integration.'})
    return r,assembly['result'],fourier['result'],field_units,equation_units,pins,inputs


def split(expression,generators,unit):
    """Reversible grade-free carriers keep polynomial collection small."""
    carriers={}
    @lru_cache(maxsize=None)
    def encode(node):
        if node in generators:return node
        if not node.has(*generators):
            if node.is_number:return node
            if node not in carriers:carriers[node]=sp.Dummy('s11cdContinuumGradeCarrier'+str(len(carriers)))
            return carriers[node]
        return node.func(*(encode(v) for v in node.args))
    encoded=encode(expression);restore={symbol:value for value,symbol in carriers.items()}
    replay=engine.memo_xreplace(encoded,restore)
    polynomial=sp.Poly(encoded,*generators,domain='EX')
    encoded_coefficients={g:c for g,c in polynomial.terms() if c!=0}
    coefficients={g:engine.memo_xreplace(c,restore) for g,c in encoded_coefficients.items()}
    monomial=lambda g:sp.prod(v**p for v,p in zip(generators,g))
    reconstructed=sp.Add(*(c*monomial(g) for g,c in encoded_coefficients.items()))
    residual=sp.expand(encoded-reconstructed)
    zero=dict.fromkeys(generators,sp.S.Zero)
    derivatives={g:sp.expand(sp.diff(encoded,*(item for v,p in zip(generators,g) for item in (v,p))).xreplace(zero)
                   /sp.prod(sp.factorial(p) for p in g)-c) for g,c in encoded_coefficients.items()}
    native={g:sp.expand(c-encoded_coefficients.get(g,sp.S.Zero)) for g,c in engine.PHYSICAL_METADATA.coefficients(encoded).items()}
    components={g:c*monomial(g) for g,c in coefficients.items()}
    retained={g:v for g,v in components.items() if all(p<=1 for p in g)}
    omitted={g:v for g,v in components.items() if g not in retained}
    return {'ORIGINAL':expression,'ENCODED':encoded,'CARRIERS':restore,'REPLAY':replay,
        'COEFFICIENTS':coefficients,'COMPONENTS':components,'RETAINED':retained,'OUTSIDE_RECTANGLE':omitted,
        'UNIT':tuple(unit),'RECONSTRUCTION_RESIDUAL':residual,'ROUND_TRIP_RESIDUAL':replay-expression,
        'DERIVATIVE_RESIDUALS':derivatives,'NATIVE_COLLECTOR_RESIDUALS':native}


def check(record):
    f.require(record['ROUND_TRIP_RESIDUAL']==0 and record['RECONSTRUCTION_RESIDUAL']==0,'exact carrier round-trip/reconstruction')
    f.require(all(v==0 for v in (*record['DERIVATIVE_RESIDUALS'].values(),*record['NATIVE_COLLECTOR_RESIDUALS'].values())),
              'independent derivative/native grade residuals')
    f.require(all(not v.has(*engine.PHYSICAL_METADATA.generators) for v in record['COEFFICIENTS'].values()),'grade-free coefficients')


def specs(assembly,fourier,field_units,equation_units,r):
    dimensions=engine.PHYSICAL_METADATA.dimensions
    for order,matrix in sorted(assembly['LOCAL_MATRICES'].items()):
        for i in range(5):
            for j in range(5):
                unit=tuple(a-b+order*c for a,b,c in zip(equation_units[i],field_units[j],dimensions.measure(r.z)))
                yield f'local{order}Row{i}Column{j}',matrix[i,j],unit,('local',int(order),i,j)
    for cell in assembly['ROWS']:
        for index,(integral,coefficient) in enumerate(cell['NONLOCAL']):
            unit=tuple(a-b for a,b in zip(equation_units[cell['ROW']],dimensions.measure(integral)))
            yield f"cellRow{cell['ROW']}Column{cell['COLUMN']}Term{index}",coefficient,unit,('cell',cell['ROW'],cell['COLUMN'],index)
    for row in fourier['ROWS']:
        for index,factor in enumerate(row['FACTORS']):
            yield f"integral{row['INDEX']}Factor{index}Coefficient",factor['COEFFICIENT'],dimensions.measure(factor['COEFFICIENT']),('factor',row['INDEX'],index)
    for index,integral in enumerate(fourier['SOURCE_INTEGRALS']):
        matches=[f for row in fourier['ROWS'] for f in row['FACTORS'] if f['SOURCE_INTEGRAL']==integral]
        f.require(matches and all(v['AMPLITUDE']==matches[0]['AMPLITUDE'] and v['CHARACTER']==matches[0]['CHARACTER'] for v in matches),'distinct source amplitude/character join')
        f.require(not matches[0]['CHARACTER'].has(*engine.PHYSICAL_METADATA.generators),'grade-independent source character')
        yield f'sourceAmplitude{index}',matches[0]['AMPLITUDE'],dimensions.measure(matches[0]['AMPLITUDE']),('source',index)


def term_joins(assembly,fourier,records):
    """All native terms, with the full grade convolution before truncation."""
    by_address={tuple(v['address']):v for v in records.values()}
    integrals={v['ORIGINAL']:v for v in fourier['ROWS']}
    sources={v:i for i,v in enumerate(fourier['SOURCE_INTEGRALS'])}
    result=[]
    for cell in assembly['ROWS']:
        for index,(integral,coefficient) in enumerate(cell['NONLOCAL']):
            row=integrals[integral];outer=by_address[('cell',cell['ROW'],cell['COLUMN'],index)]
            factors=[]
            for fi,factor in enumerate(row['FACTORS']):
                source_index=sources[factor['SOURCE_INTEGRAL']]
                middle=by_address[('factor',row['INDEX'],fi)];source=by_address[('source',source_index)]
                combinations=[]
                for ga,ca in outer['record']['COEFFICIENTS'].items():
                    for gb,cb in middle['record']['COEFFICIENTS'].items():
                        for gc,cc in source['record']['COEFFICIENTS'].items():
                            grade=tuple(a+b+c for a,b,c in zip(ga,gb,gc))
                            combinations.append({'grade':grade,'factorGrades':(ga,gb,gc),
                                'coefficients':(ca,cb,cc),'retained':all(v<=1 for v in grade)})
                factors.append({'factor':fi,'sourceIndex':source_index,'character':factor['CHARACTER'],
                    'frequency':factor['FREQUENCY'],'combinations':combinations})
            result.append({'row':cell['ROW'],'column':cell['COLUMN'],'term':index,'integralIndex':row['INDEX'],
                'originalIntegral':integral,'cellCoefficient':coefficient,'sourceLimit':row['SOURCE_LIMIT'],
                'remainingLimits':row['REMAINING_LIMITS'],'factors':factors})
    return result


def emit_records(records,joins,r):
    for key,item in records.items():
        record=item['record'];unit=record['UNIT'];name=PREFIX+'_'+key
        values=sp.Tuple(record['ORIGINAL'],*record['COMPONENTS'].values())
        engine.physical(name+'_OPERANDS',values,zero_dimensions={p:unit for p,_ in engine.leaves(values)})
        residuals=sp.Tuple(record['ROUND_TRIP_RESIDUAL'],record['RECONSTRUCTION_RESIDUAL'],
            *record['DERIVATIVE_RESIDUALS'].values(),*record['NATIVE_COLLECTOR_RESIDUALS'].values())
        engine.physical(name+'_RESIDUALS',residuals,zero_dimensions={p:unit for p,_ in engine.leaves(residuals)})
        structural(name+'_GRADE_INDEX',tuple(record['COMPONENTS']))
    # Structural addresses remain lossless; the corresponding physical factors
    # and their independent grade metadata were emitted above.
    census=[{'row':v['row'],'column':v['column'],'term':v['term'],'integralIndex':v['integralIndex'],
        'factorGrades':[[c['factorGrades'] for c in factor['combinations']] for factor in v['factors']]} for v in joins]
    structural(PREFIX+'_TERM_GRADE_JOINS',census)


def structural(name,value):
    body=engine.cas(value);engine.emit(name,body)
    zeros={p:(0,0,0) for p,_ in engine.leaves(body)}
    engine.emit('METADATA_'+name,engine.PHYSICAL_METADATA.record(body,zeros))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--focused',action='store_true');args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def budget(*_):raise TimeoutError('continuum grade initial budget; preserve completed records')
    signal.signal(signal.SIGALRM,budget);signal.alarm(900);started=time.monotonic()
    r,assembly,fourier,field_units,equation_units,pins,inputs=load(base)
    generators=engine.PHYSICAL_METADATA.generators;records={};inventory=[]
    (base/'records').mkdir()
    chosen=list(specs(assembly,fourier,field_units,equation_units,r))
    if args.focused:
        selected=[]
        for kind in ('local','cell','factor','source'):
            available=[v for v in chosen if v[3][0]==kind and v[1]!=0]
            candidates=[v for v in available if v[1].has(*generators)]
            f.require(available,('actual focused operand',kind))
            selected.append((candidates or available)[0])
        chosen=selected
    for key,value,unit,address in chosen:
        f.require(unit is not None and len(unit)==3,('restored coefficient unit',key))
        record=split(value,generators,unit);item={'address':address,'record':record}
        path=base/'records'/(key+'.pickle');f.atomic_pickle(path,item)
        inventory.append({'key':key,'address':address,'path':str(path.relative_to(base)),'sha256':f.digest(path),
                          'grades':[list(g) for g in record['COMPONENTS']],'carriers':len(record['CARRIERS'])})
        f.save(base/'record-inventory.json',inventory)
        check(record);records[key]=item
        with (base/'progress.jsonl').open('a') as stream:stream.write(json.dumps({'records':len(records),'key':key,'wallSeconds':time.monotonic()-started})+'\n')
    controls=[]
    for key,item in records.items():
        record=item['record']
        if not record['COEFFICIENTS']:continue
        polynomial=sp.Poly(record['ENCODED'],*generators,domain='EX')
        grade,coefficient=next((g,c) for g,c in polynomial.terms() if c!=0)
        mutation=sp.expand(record['ENCODED']-(record['ENCODED']+coefficient*sp.prod(v**p for v,p in zip(generators,grade))/100))
        controls.append({'key':key,'grade':grade,'residual':mutation})
    f.atomic_pickle(base/'coefficient-controls.pickle',controls)
    f.require(controls and all(v['residual']!=0 for v in controls),'one-sided coefficient controls')
    if args.focused:
        # Test a mixed slot with symbols from the actual instance in the same
        # algorithm; this synthetic diagnostic is not a physical operator.
        eta,sigma=generators[1:];mixed=split((1+eta)*(1+sigma),generators,(0,0,0));check(mixed)
        mixed_residual=sp.expand(mixed['ENCODED']-(mixed['ENCODED']-eta*sigma))
        f.atomic_pickle(base/'mixed-control.pickle',{'record':mixed,'omissionResidual':mixed_residual})
        f.require(mixed_residual!=0 and (0,1,1) in mixed['COMPONENTS'],'missing mixed-grade control')
        joins=[]
    else:joins=term_joins(assembly,fourier,records)
    result={'records':records,'termJoins':joins,'generators':generators,'sourceFiles':pins,'inputPackets':inputs,
        'fieldUnits':field_units,'equationUnits':equation_units,'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions)),
        'finiteCutoffs':fourier['CUTOFFS'],'scope':'Unintegrated native coefficient operators; continuum solution and channel-map expansion remain.'}
    f.atomic_pickle(base/'continuum-grades.pickle',result);before=f.digest(base/'continuum-grades.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_records(records,joins,r)
        keys={tag:'s11cdContinuumGrade'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')}
        structural(PREFIX+'_WRITE_KEYS',keys)
        index=engine.emission_index(engine.EMISSION_LINES)
        structural(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ')
        f.require(tag not in entries,'unique emitted tag');entries[tag]=_restore(body)
    seen=set();original=engine.emit
    def compare(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('emission replay',tag));seen.add(tag)
    engine.emit=compare
    try:
        emit_records(records,joins,r);structural(PREFIX+'_WRITE_KEYS',keys)
        structural(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=original
    f.require(seen==set(entries) and len(set(keys.values()))==len(keys) and not set(keys.values())&set(engine.IMPORT_KEYS),'complete emission/key replay')
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    metadata_paths=0
    for tag,value in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        for path,descriptor in value:
            fields={str(k):v for k,v in descriptor};unit=fields['DIMENSION_L_T_M']
            f.require(isinstance(unit,sp.Tuple) and len(unit)==3 and all(not v.free_symbols for v in unit),'resolved emitted units')
            f.require('MULTIGRADE' in fields and 'EPSILON_LAMBDA_SUPPORT' in fields,'emitted independent/path grades');metadata_paths+=1
    f.require(before==f.digest(base/'continuum-grades.pickle'),'unchanged pre/post packet')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()) and all(f.digest(Path(n))==h for n,h in inputs.items()),'unchanged consumed sources/packets')
    f.require(all(f.digest(base/v['path'])==v['sha256'] for v in inventory),'saved record hashes')
    summary={'runDirectory':str(base),'focused':args.focused,'sourceFiles':pins,'inputPackets':inputs,
        'records':len(records),'nativeTermJoins':len(joins),'sourceIntegrals':len(fourier['SOURCE_INTEGRALS']),
        'recordKinds':{kind:sum(v['address'][0]==kind for v in inventory) for kind in ('local','cell','factor','source')},
        'coefficientControls':len(controls),'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':metadata_paths,
        'allGrades':sorted({tuple(g) for v in inventory for g in v['grades']}),
        'combinedTermGrades':sorted({c['grade'] for v in joins for factor in v['factors'] for c in factor['combinations']}),
        'outsideRectangleFactorGrades':sum(len(v['record']['OUTSIDE_RECTANGLE']) for v in records.values()),
        'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'continuum-grades.pickle'),
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.iterdir() if p.suffix in ('.out','.pickle')},
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope':result['scope']}
    f.save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
