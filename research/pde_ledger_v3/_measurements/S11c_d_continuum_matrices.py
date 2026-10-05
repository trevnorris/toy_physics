#!/usr/bin/env python3
"""Assemble native continuum coefficient matrices from accepted quadratures."""
import argparse
import contextlib
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from types import SimpleNamespace

import numpy as np
import sympy as sp
import S11c_d_continuum_grades as grades
import S11c_d_finite_scattering_domain as domain

f=grades.f
engine=f.engine
PLAN=f.M/'S11c_d_continuum_matrix_plan.md'
GRADE_CHECKPOINT=f.M/'S11c_d_continuum_grade_checkpoint.json'
FINITE_CHECKPOINT=f.M/'S11c_d_finite_scattering_domain_checkpoint.json'
PREFIX='CONTINUUM_INTERIOR_MATRICES_LAB_HELD_RHO4_CONSTANT'


def load(base):
    packet,accepted,path=f.accepted_packet(GRADE_CHECKPOINT,'continuum-grades.pickle')
    finite=json.loads(FINITE_CHECKPOINT.read_text())
    f.require(finite['status']=='VALIDATED_FINITE_BOUNDARY_REGULATOR_RESPONSE','accepted finite response')
    directory=Path(finite['checks']['records'][-1]['checks']['runDirectory'])
    f.require(f.digest(Path(finite['runDirectory'])/'checks.json')==finite['checksSha256'],'finite acceptance join')
    _,_,inspect_case,_=domain.adapters();observed=domain.observable(directory,inspect_case)
    f.require(observed['checks']==finite['checks']['records'][-1]['checks'],'finite case checkpoint join')
    bound=f.unpickle(directory/'domain-binding.pickle')['bound']
    source=f.unpickle(directory/'source-binding.pickle')
    system=f.unpickle(directory/'finite-system.pickle');solution=f.unpickle(directory/'finite-solution.pickle')
    original_inputs=json.loads((directory/'inputs.json').read_text())
    reduction_path=next(Path(n) for n in packet['inputPackets'] if Path(n).name=='reduced-action.pickle')
    f.require(f.digest(reduction_path)==packet['inputPackets'][str(reduction_path)],'reduction context hash')
    r,dimensions=f.prior.domain.momentum.source.native.source.restore_context(f.unpickle(reduction_path))
    dimensions.__dict__.update(packet['dimensionState'])
    specification=json.loads((f.M/'S11c_d_variable_profile_development_input.json').read_text())
    f.require(specification==original_inputs['input'],'approved input identity')
    adapter=engine.NumericalReducedAction(SimpleNamespace(r=r),{},specification)
    by_address={tuple(v['address']):v['record'] for v in packet['records'].values()}
    rows={row['index']:row for row in bound['rows']}
    for (kind,*address),record in by_address.items():
        if kind in ('factor','source'):
            f.require(set(record['COEFFICIENTS'])=={(0,0,0)},('actual grade-independent saved operand',kind,address))
        if kind=='factor':
            factor=rows[address[0]]['factors'][address[1]]
            f.require(factor['symbolicCoefficient']==record['ORIGINAL'],'native momentum coefficient source')
        elif kind=='source':
            f.require(bound['sources'][(0,address[0])]['symbolicAmplitude']==record['ORIGINAL'],'native source amplitude')
    for term in packet['termJoins']:
        row=rows[term['integralIndex']]
        f.require(row['original']==term['originalIntegral'],'native original integral identity')
        expected_limits=tuple(tuple(v.xreplace(bound['cutoffBindings']) for v in limit) for limit in term['remainingLimits'])
        f.require(row['limits']==expected_limits and row['sourceLimit'][1:]==(-64,64),'actual ordered integration limits')
        f.require(all(c['factorGrades'][1:]==((0,0,0),(0,0,0)) for factor in term['factors'] for c in factor['combinations']),
                  'reuse only grade-independent momentum/source factors')
    row_matrices={}
    for index,group in enumerate(solution['groups'],1):
        f.require(f.digest(directory/f'layout-{index}.pickle')==observed['checks']['artifacts'][f'layout-{index}.pickle']['sha256'],'accepted layout hash')
        f.require(group['setting']==system['settings'] and np.max(abs(group['actionResidual']))<1e-10,'same quadrature and direct action')
        row_matrices.update(zip(group['rowIndices'],group['matrices']))
    f.require(set(row_matrices)==set(range(80)) and len(source['jets'])==35,'complete reused rows and sources')
    pins=dict(accepted['sourceFiles'])
    for name,sha in pins.items():f.require(f.digest(f.ROOT/name)==sha,'current grade source')
    for p in (Path(__file__),PLAN,GRADE_CHECKPOINT,FINITE_CHECKPOINT,Path(domain.__file__),Path(f.__file__),
              Path(domain.study.__file__),f.ACCEPTANCE):pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    operands={str(path):f.digest(path),str(directory/'checks.json'):f.digest(directory/'checks.json'),
              str(directory/'inputs.json'):f.digest(directory/'inputs.json'),str(reduction_path):f.digest(reduction_path)}
    operands.update({str(directory/n):v['sha256'] for n,v in observed['checks']['artifacts'].items()})
    for name in pins:
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((f.ROOT/name).read_bytes())
    f.save(base/'inputs.json',{'sourceFiles':pins,'inputPackets':operands,'finiteDirectory':str(directory),
        'settings':system['settings'],'gradeGenerators':list(map(str,packet['generators'])),
        'gradeOrigin':{str(k):str(v) for k,v in adapter.input.origin.items()},'fieldUnits':original_inputs['fieldUnits'],
        'equationUnits':original_inputs['equationUnits'],'operatorBlockUnits':original_inputs['operatorBlockUnits'],
        'scope':'Interior coefficient matrices before boundary replacement; existing momentum matrices reused.'})
    return r,packet,bound,system,row_matrices,adapter,pins,operands


def values(expression,r,system):
    nodes=system['nodes'];regulator=system['settings']['regulator']
    f.require(not (expression.free_symbols-{r.z,r.regulator}),'fully bound coefficient arguments')
    result=np.broadcast_to(np.asarray(sp.lambdify((r.z,r.regulator),expression,'numpy',cse=True)(nodes,regulator),complex),nodes.shape)
    f.require(np.isfinite(result).all(),'finite coefficient values');return result


def assemble(records,terms,adapter,r,system,row_matrices,graded):
    size=len(system['nodes']);shape=system['unreplacedOperator'].shape
    support=sorted({g for item in records.values() for g in item['record']['COEFFICIENTS']}) if graded else [(0,0,0)]
    local={g:np.zeros(shape,complex) for g in support};nonlocal_={g:np.zeros(shape,complex) for g in support}
    addresses={tuple(v['address']):v['record'] for v in records.values()}
    bound_records={}
    for address,record in addresses.items():
        if address[0] not in ('local','cell'):continue
        coefficients=record['COEFFICIENTS'] if graded else {(0,0,0):record['ORIGINAL']}
        bound_records[address]={g:adapter.bind(c) for g,c in coefficients.items()}
    for address,coefficients in bound_records.items():
        if address[0]!='local':continue
        _,order,i,j=address
        for grade,coefficient in coefficients.items():
            local[grade][i*size:(i+1)*size,j*size:(j+1)*size]+=values(coefficient,r,system)[:,None]*system['derivativeMatrices'][order]
    for term in terms:
        i,j,index=term['row'],term['column'],term['integralIndex']
        for grade,coefficient in bound_records[('cell',i,j,term['term'])].items():
            nonlocal_[grade][i*size:(i+1)*size,j*size:(j+1)*size]+=values(coefficient,r,system)[:,None]*row_matrices[index]
    return {'local':local,'nonlocal':nonlocal_,'total':{g:local[g]+nonlocal_[g] for g in support},'boundCoefficients':bound_records}


def recombine(matrices,origin,generators):
    weights={g:float(sp.prod(origin.get(v,sp.S.One)**p for v,p in zip(generators,g))) for g in matrices}
    return sum((weights[g]*a for g,a in matrices.items()),np.zeros_like(next(iter(matrices.values()))))


def differences(actual,expected,size):
    residual=actual-expected
    blocks=[]
    for i in range(5):
        for j in range(5):
            rows=slice(i*size,(i+1)*size);columns=slice(j*size,(j+1)*size)
            delta=residual[rows,columns];reference=expected[rows,columns]
            absolute=float(np.max(abs(delta)));scale=1+float(np.max(abs(reference)))
            blocks.append({'row':i,'column':j,'absolute':absolute,'scaledReferenceFrame':absolute/scale})
    return {'residual':residual,'blocks':blocks,'maximumScaledReferenceFrame':max(v['scaledReferenceFrame'] for v in blocks)}


def fingerprint(array):
    array=np.ascontiguousarray(array,dtype=np.complex128)
    rows,columns=np.indices(array.shape,dtype=np.int64);projections=[]
    for prime in (101,103,107):
        weights=(((rows+1)*(columns+prime)+(rows+prime)**2)%prime+1)/prime
        projections.append(complex(np.sum(array*weights)))
    return {'SHAPE':array.shape,'DTYPE':array.dtype.str,'C_ORDER_SHA256':hashlib.sha256(array.tobytes()).hexdigest(),
            'NUMERIC_TENSOR_PROJECTIONS':sp.Tuple(*(engine.FullPencilModes.number(v) for v in projections))}


def emit_result(result,r):
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r
    size=result['size'];generators=result['generators']
    zero=(0,0,0)
    def numeric(name,value,unit):
        body=engine.cas(value);engine.emit(PREFIX+'_'+name,body)
        engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(body,unit))
    for kind,matrices in result['matrices'].items():
        for grade,matrix in matrices.items():
            weight=sp.prod(v**p for v,p in zip(generators,grade))
            for i in range(5):
                for j in range(5):
                    block=matrix[i*size:(i+1)*size,j*size:(j+1)*size]
                    unit=result['blockUnits'][i][j];body=fingerprint(block)
                    body['NUMERIC_TENSOR_PROJECTIONS']=sp.Tuple(*(v*weight for v in body['NUMERIC_TENSOR_PROJECTIONS']))
                    body['COEFFICIENT_ENTRY_UNIT']=unit;body['COEFFICIENT_GRADE']=grade
                    body['COMPONENT_WEIGHT']=weight;body['LAMBDA_ORDER']=grade[1]+grade[2]
                    numeric(kind+'_'+''.join(map(str,grade))+f'_ROW_{i}_COLUMN_{j}',body,
                        lambda path,u=unit:u if path and path[0]=='NUMERIC_TENSOR_PROJECTIONS' else zero)
    for name,test in result['comparisons'].items():
        for block in test['blocks']:
            i,j=block['row'],block['column'];unit=result['blockUnits'][i][j]
            body={'MAXIMUM_ABSOLUTE_RESIDUAL':engine.FullPencilModes.number(block['absolute']),
                  'SCALED_REFERENCE_FRAME_RESIDUAL':engine.FullPencilModes.number(block['scaledReferenceFrame'])}
            numeric('COMPARISON_'+name+f'_ROW_{i}_COLUMN_{j}',body,
                    lambda path,u=unit:u if path and path[0]=='MAXIMUM_ABSOLUTE_RESIDUAL' else zero)
    length=engine.PHYSICAL_METADATA.dimensions.measure(r.z)
    setting_units={'sourceBound':length,'momentumBound':tuple(-v for v in length)}
    numeric('SETTINGS',result['settings'],lambda path:setting_units.get(path[0],zero))
    grades.structural(PREFIX+'_INPUT_MANIFEST',{'reusedLayoutHashes':result['reusedLayoutHashes'],
        'finiteCase':result['finiteCase'],'coefficientMatrixShape':next(iter(result['matrices']['total'].values())).shape})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('continuum matrix budget; preserve completed operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);started=time.monotonic()
    r,packet,bound,system,row_matrices,adapter,pins,operands=load(base)
    assembled=assemble(packet['records'],packet['termJoins'],adapter,r,system,row_matrices,True)
    f.atomic_pickle(base/'bound-coefficients.pickle',assembled['boundCoefficients'])
    matrices={kind:assembled[kind] for kind in ('local','nonlocal','total')}
    for grade in matrices['total']:
        f.atomic_pickle(base/('component-'+''.join(map(str,grade))+'.pickle'),{kind:data[grade] for kind,data in matrices.items()})
    size=len(system['nodes']);origin=adapter.input.origin;generators=packet['generators']
    comparisons={'approvedTotal':differences(recombine(matrices['total'],origin,generators),system['unreplacedOperator'],size),
                 'approvedLocal':differences(recombine(matrices['local'],origin,generators),system['localMatrix'],size)}
    for name,eta,sigma in [('zero',sp.S.Zero,sp.S.Zero),('independent',sp.Rational(1,137),sp.Rational(1,911))]:
        other=engine.NumericalReducedAction(SimpleNamespace(r=r),{},adapter.input.specification)
        other.input.origin.update({generators[1]:eta,generators[2]:sigma})
        native=assemble(packet['records'],packet['termJoins'],other,r,system,row_matrices,False)
        f.atomic_pickle(base/('direct-'+name+'.pickle'),native)
        comparison=differences(recombine(matrices['total'],other.input.origin,generators),native['total'][(0,0,0)],size)
        comparison['formalGradePoint']=(eta,sigma);comparisons[name]=comparison
    mixed=(0,1,1);f.require(mixed in matrices['total'],'computed mixed coefficient')
    full=recombine(matrices['total'],origin,generators)
    mutated=recombine({g:a for g,a in matrices['total'].items() if g!=mixed},origin,generators)
    comparisons['mixedOmission']=differences(mutated,full,size)
    f.atomic_pickle(base/'comparisons.pickle',comparisons)
    f.require(all(v['maximumScaledReferenceFrame']<1e-11 for k,v in comparisons.items() if k!='mixedOmission'),'literal original/coefficient recombination')
    f.require(comparisons['mixedOmission']['maximumScaledReferenceFrame']>0,'actual mixed-term omission response')
    inputs=json.loads((base/'inputs.json').read_text())
    block_units=[[tuple(map(sp.Rational,u)) for u in row] for row in inputs['operatorBlockUnits']]
    result={'matrices':matrices,'comparisons':comparisons,'blockUnits':block_units,'size':size,'generators':generators,
        'gradeOrigin':origin,'settings':system['settings'],'finiteCase':inputs['finiteDirectory'],
        'reusedLayoutHashes':{n:h for n,h in operands.items() if Path(n).name.startswith('layout-')},
        'sourceFiles':pins,'inputPackets':operands,'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions)),
        'scope':'Finite interior coefficient operators, before actual boundary/channel response expansion.'}
    f.atomic_pickle(base/'continuum-matrices.pickle',result);before=f.digest(base/'continuum-matrices.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result,r)
        keys={tag:'s11cdContinuumMatrix'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')}
        grades.structural(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES)
        grades.structural(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique matrix tags');entries[tag]=grades._restore(body)
    original=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('matrix emission replay',tag));seen.add(tag)
    engine.emit=replay
    try:
        emit_result(result,r);grades.structural(PREFIX+'_WRITE_KEYS',keys);grades.structural(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=original
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'full matrix emission/key replay')
    metadata_paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        structural=tag.endswith(('_INPUT_MANIFEST','_WRITE_KEYS','_EMISSION_LINES'))
        for item in body:
            descriptor=item[1] if structural else item
            fields={str(k):v for k,v in descriptor};unit=fields['DIMENSION_L_T_M']
            f.require(isinstance(unit,sp.Tuple) and len(unit)==3 and all(not v.free_symbols for v in unit),'resolved coefficient tensor units')
            f.require('MULTIGRADE' in fields and 'EPSILON_LAMBDA_SUPPORT' in fields,'independent coefficient/path grades')
            metadata_paths+=1 if structural else len(fields['PATHS'])
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    f.require(before==f.digest(base/'continuum-matrices.pickle'),'matrix pre/post packet hash')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()) and all(f.digest(Path(n))==h for n,h in operands.items()),'source/input post hashes')
    summary={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':operands,'size':size,'unknowns':5*size,
        'grades':list(matrices['total']),'nativeTerms':len(packet['termJoins']),'reusedMomentumRows':len(row_matrices),
        'newIntegrationNodes':0,'settings':system['settings'],'comparisons':{k:v['maximumScaledReferenceFrame'] for k,v in comparisons.items()},
        'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':metadata_paths,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'continuum-matrices.pickle'),
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
