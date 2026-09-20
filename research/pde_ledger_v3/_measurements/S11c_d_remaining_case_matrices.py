#!/usr/bin/env python3
"""Evaluate only new case rows, then assemble actual interior coefficient matrices."""
import argparse
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import time
from types import SimpleNamespace
import numpy as np
import sympy as sp
import S11c_d_remaining_case_bindings as binding

f, engine, interior = binding.f, binding.engine, binding.interior
PLAN=f.M/'S11c_d_remaining_case_matrices_plan.md'
CHECKPOINT=f.M/'S11c_d_remaining_case_bindings_checkpoint.json'
BASELINE=binding.BASELINE


def copy_packet(origin, destination, expected, inputs):
    f.require(f.digest(origin)==expected, ('accepted packet',str(origin)))
    destination.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(origin,destination)
    f.require(f.digest(destination)==expected,'byte-identical input copy')
    inputs[str(origin)]=expected
    return f.unpickle(destination)


def context(base,label,case,specification):
    packets={kind:f.unpickle(base/'accepted-cases'/label/(kind+'.pickle'))
             for kind in ('reduced-action','actions','assembly','factorization')}
    r,dimensions,pencil=binding.factors.context(packets['reduced-action'],packets['actions'],packets['assembly'])
    dimensions.__dict__.update(case['grades']['dimensionState'])
    adapter=engine.NumericalReducedAction(pencil,packets['assembly']['result'],specification)
    return r,adapter,packets


def load(base,resume=None):
    cp=json.loads(CHECKPOINT.read_text());origin=Path(cp['runDirectory'])
    f.require(cp['status']=='ACCEPTED_BINDINGS_AND_GRADES','accepted four-case bindings')
    f.require(f.digest(origin/'checks.json')==cp['checksSha256'],'binding final checks')
    checks=json.loads((origin/'checks.json').read_text())
    for name,h in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==f.digest(origin/'source'/name)==h,('binding source',name))
    if resume is not None:
        old=json.loads((resume/'checks.json').read_text())
        f.require(old['mode']=='focused' and old['status']=='VALIDATED_ROW_REUSE_AND_ASSEMBLY', 'accepted focused operands')
        for name,h in old['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(resume/'source'/name)==h,('focused source',name))
        for name,item in old['artifacts'].items():
            f.require(f.digest(resume/name)==item['sha256'],'completed focused artifact')
            target=base/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(resume/name,target)
            f.require(f.digest(target)==item['sha256'],'exact focused artifact reuse')
        inputs={str(resume/'checks.json'):f.digest(resume/'checks.json')}
        inputs.update({str(resume/name):item['sha256'] for name,item in old['artifacts'].items()})
        cases={label:f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle') for label in cp['cases']}
        system=f.unpickle(base/'accepted-finite-system.pickle')
        matrices=f.unpickle(base/'accepted-continuum-matrices.pickle')
        baseline_rows=f.unpickle(base/'accepted-row-matrices.pickle')
    else:
        inputs={str(origin/'checks.json'):f.digest(origin/'checks.json')};cases={}
        for label in cp['cases']:
            name='cases/'+label+'/case-binding.pickle'
            cases[label]=copy_packet(origin/name,base/'accepted-bindings'/label/'case-binding.pickle',cp['artifacts'][name]['sha256'],inputs)
            for kind in ('reduced-action','actions','assembly'):
                name='accepted-cases/'+label+'/'+kind+'.pickle'
                copy_packet(origin/name,base/name,cp['artifacts'][name]['sha256'],inputs)
            name='cases/'+label+'/factorization.pickle'
            copy_packet(origin/name,base/'accepted-cases'/label/'factorization.pickle',cp['artifacts'][name]['sha256'],inputs)
        system=copy_packet(origin/'accepted-finite-system.pickle',base/'accepted-finite-system.pickle',cp['artifacts']['accepted-finite-system.pickle']['sha256'],inputs)
        mc=json.loads((f.M/'S11c_d_continuum_matrix_checkpoint.json').read_text())
        f.require(mc['status']=='PUBLISHED_ANNEX_VERIFIED','accepted baseline interior matrices')
        md=Path(mc['checks']['runDirectory']);name='continuum-matrices.pickle'
        matrices=copy_packet(md/name,base/('accepted-'+name),mc['checks']['artifacts'][name]['sha256'],inputs)
        directory=Path(checks['acceptedFiniteDirectory']);baseline_rows={}
        for count in (1,2,3):
            path=directory/f'layout-{count}.pickle'
            group=copy_packet(path,base/'accepted-layouts'/path.name,mc['checks']['inputPackets'][str(path)],inputs)
            f.require(group['setting']==system['settings'],'same complete baseline quadrature rule')
            f.require(np.isfinite(group['matrices']).all() and abs(group['massResidual'])<1e-9*(1+abs(group['mass'])),'accepted full layout measures')
            f.require(np.max(abs(group['actionResidual']))/(1+np.max(abs(group['direct'])))<1e-10,'accepted full layout direct actions')
            baseline_rows.update(zip(group['rowIndices'],group['matrices']))
        f.require(set(baseline_rows)==set(range(len(cases[BASELINE]['binding']['bound']['rows']))),'complete accepted baseline rows')
        f.atomic_pickle(base/'accepted-row-matrices.pickle',baseline_rows)
    pins=dict(cp['sourceFiles'])
    for path in (Path(__file__).resolve(),PLAN,CHECKPOINT,f.M/'S11c_d_continuum_matrix_checkpoint.json'):
        pins[str(path.relative_to(f.ROOT))]=f.digest(path)
    for name,h in pins.items():
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,target)
        f.require(f.digest(target)==h,'new frozen source')
    manifest={'sourceFiles':pins,'inputPackets':inputs,'settings':system['settings'],'input':checks['input'],
        'bindingChecksSha256':cp['checksSha256'],'fieldCoefficients':len(system['nodes']),
        'scope':'All four finite interior operators and independent grades, before fresh case-specific boundary replacement. No new modes, currents, scattering or pole search.',
        'nativeBodies':{name:f.digest(Path(module.__file__)) for name,module in [('finite',f),('interior',interior)]},
        'resumeFrom':str(resume) if resume else None}
    manifest['copiedInputs']={str(p.relative_to(base)):f.digest(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
    f.save(base/'inputs.json',manifest)
    return cases,system,matrices,baseline_rows,manifest


def preflight(base,cases,system):
    summaries={};new_total=0;reuse_total=0
    for label,case in cases.items():
        data=case['binding'];bound=data['bound'];grades=case['grades']
        f.require(data['settings']==system['settings'],'full settings identity before reuse')
        f.require(len(data['fieldUnits'])==len(data['equationUnits'])==5,'actual five field/equation slots')
        for item in grades['records'].values():
            kind,*address=item['address'];record=item['record']
            if kind in ('factor','source'):
                f.require(set(record['COEFFICIENTS'])<={(0,0,0)},('actual grade-free quadrature operand',label,kind,address))
            if kind=='factor':
                actual=bound['rows'][address[0]]['factors'][address[1]]['symbolicCoefficient']
                f.require(binding.same(record['ORIGINAL'],actual),'factor grade/binding identity')
            elif kind=='source':
                f.require(binding.same(record['ORIGINAL'],bound['sources'][0,address[0]]['symbolicAmplitude']),'source grade/binding identity')
        for term in grades['termJoins']:
            row=bound['rows'][term['integralIndex']]
            f.require(binding.same(term['originalIntegral'],row['original']),'full native integral in term convolution')
            f.require(all(c['factorGrades'][1:]==((0,0,0),(0,0,0)) for factor in term['factors'] for c in factor['combinations']), 'no omitted quadrature grade')
        reused={v['row'] for v in case['reusedRows']};fresh=set(case['newRows'])
        f.require(not reused&fresh and reused|fresh==set(range(len(bound['rows']))),'complete disjoint row partition')
        for item in case['reusedRows']:
            old=cases[item['fromCase']]['binding'];row=bound['rows'][item['row']];prior=old['bound']['rows'][item['fromRow']]
            f.require(binding.same(binding.signature(row,data),binding.signature(prior,old)),'literal numerical reuse signature')
            f.require(data['settings']==old['settings'] and binding.same(bound['profiles'],old['bound']['profiles'])
                      and binding.same(data['fieldUnits'],old['fieldUnits']) and binding.same(bound['abel'],old['bound']['abel']), 'full cache, field and measure reuse')
            changed=dict(row);limit=row['limits'][0]
            changed['limits']=(sp.Tuple(limit[0],limit[1],limit[2]+1),*row['limits'][1:])
            f.require(not binding.same(binding.signature(changed,data),binding.signature(prior,old)),'wrong-limit numerical reuse rejects')
        for index in fresh:f.require(len(bound['rows'][index]['limits'])<=2,'measured new-layout scope excludes triple integration')
        summaries[label]={'rows':len(bound['rows']),'sources':len(data['jets']),'terms':len(grades['termJoins']),
            'newRows':case['newRows'],'reusedRows':case['reusedRows'],'fieldUnits':list(map(str,data['fieldUnits'])),
            'equationUnits':list(map(str,data['equationUnits']))}
        new_total+=len(fresh);reuse_total+=len(reused)
    result={'cases':summaries,'newUnionRows':new_total,'reusedCaseRows':reuse_total,
        'nativeAccumulatorSha256':__import__('hashlib').sha256(inspect.getsource(f.BasisMomentum).encode()).hexdigest(),
        'nativeAssemblySha256':__import__('hashlib').sha256(inspect.getsource(interior.assemble).encode()).hexdigest()}
    f.save(base/'preflight.json',result);return result


def direct_cells(case,packets,adapter,r,system,rows):
    size=len(system['nodes']);local=np.zeros((5*size,5*size),complex);nonlocal_=np.zeros_like(local)
    for order,matrix in case['binding']['local'].items():
        for i in range(5):
            for j in range(5):local[i*size:(i+1)*size,j*size:(j+1)*size]+=interior.values(matrix[i,j],r,system)[:,None]*system['derivativeMatrices'][order]
    assembly=packets['assembly']['result'];fourier=packets['factorization']['result']
    lookup=binding.integral_addresses(assembly,fourier);terms=0
    for cell in assembly['ROWS']:
        i,j=cell['ROW'],cell['COLUMN']
        for integral,coefficient in cell['NONLOCAL']:
            index=lookup[id(integral)]
            for factor in case['binding']['bound']['rows'][index]['factors']:
                f.require(case['binding']['jets'][factor['sourceIndex']]['column']==j,'actual generic source belongs to native cell column')
            nonlocal_[i*size:(i+1)*size,j*size:(j+1)*size]+=interior.values(adapter.bind(coefficient),r,system)[:,None]*rows[index]
            terms+=1
    f.require(terms==len(case['grades']['termJoins']),'complete independent native cell contraction')
    return {'local':local,'nonlocal':nonlocal_,'total':local+nonlocal_,'terms':terms}


def assemble_case(target,case,packets,r,adapter,system,rows,baseline=None):
    grade=case['grades'];size=len(system['nodes'])
    matrices=interior.assemble(grade['records'],grade['termJoins'],adapter,r,system,rows,True)
    f.atomic_pickle(target/'coefficient-matrices.pickle',matrices)
    direct=direct_cells(case,packets,adapter,r,system,rows);f.atomic_pickle(target/'direct-native-cells.pickle',direct)
    unsplit=interior.assemble(grade['records'],grade['termJoins'],adapter,r,system,rows,False)
    f.atomic_pickle(target/'direct-unsplit.pickle',unsplit)
    comparisons={}
    for kind in ('local','nonlocal','total'):
        comparisons['native_'+kind]=interior.differences(unsplit[kind][(0,0,0)],direct[kind],size)
        comparisons['approved_'+kind]=interior.differences(interior.recombine(matrices[kind],adapter.input.origin,grade['generators']),direct[kind],size)
        if baseline is not None:
            for g,array in matrices[kind].items():comparisons['baseline_'+kind+str(g)]=interior.differences(array,baseline['matrices'][kind][g],size)
    # A formal arithmetic control of the assembled polynomial; no new physical case.
    other=engine.NumericalReducedAction(SimpleNamespace(r=r),{},adapter.input.specification)
    other.input.origin.update({grade['generators'][1]:sp.Rational(1,137),grade['generators'][2]:sp.Rational(1,911)})
    control=interior.assemble(grade['records'],grade['termJoins'],other,r,system,rows,False)
    f.atomic_pickle(target/'direct-independent.pickle',control)
    comparisons['independent']=interior.differences(interior.recombine(matrices['total'],other.input.origin,grade['generators']),control['total'][(0,0,0)],size)
    mixed=(0,1,1);full=interior.recombine(matrices['total'],adapter.input.origin,grade['generators'])
    omitted=interior.recombine({g:a for g,a in matrices['total'].items() if g!=mixed},adapter.input.origin,grade['generators'])
    mutation=interior.differences(omitted,full,size)
    f.atomic_pickle(target/'comparisons.pickle',{'comparisons':comparisons,'mixedOmission':mutation})
    f.require(all(v['maximumScaledReferenceFrame']<1e-10 for v in comparisons.values()),'independent original-cell and grade matrix comparisons')
    if mixed in matrices['total'] and np.any(matrices['total'][mixed]):
        f.require(mutation['maximumScaledReferenceFrame']>0,'actual mixed coefficient omission responds')
    block_units=[[tuple(a-b) for b in map(sp.ImmutableMatrix,case['binding']['fieldUnits'])]
                 for a in map(sp.ImmutableMatrix,case['binding']['equationUnits'])]
    result={'matrices':matrices,'direct':direct,'comparisons':comparisons,'mixedOmission':mutation,
        'size':size,'generators':grade['generators'],'gradeOrigin':adapter.input.origin,'blockUnits':block_units,
        'fieldUnits':case['binding']['fieldUnits'],'equationUnits':case['binding']['equationUnits'],
        'dimensionState':grade['dimensionState'],'settings':system['settings'],
        'scope':'Actual finite interior matrices before new case-specific boundary/current construction.'}
    f.atomic_pickle(target/'interior-matrices.pickle',result)
    return {'nativeTerms':direct['terms'],'grades':list(matrices['total']),
        'maximumComparison':max(v['maximumScaledReferenceFrame'] for v in comparisons.values()),
        'mixedOmission':mutation['maximumScaledReferenceFrame'],'unknowns':5*size}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',required=True,type=Path)
    parser.add_argument('--focused',action='store_true');parser.add_argument('--resume-from',type=Path)
    args=parser.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));start=time.monotonic()
    def timeout(*_):raise TimeoutError('case matrix budget; retain every completed layout and matrix')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    f.require(args.focused!=(args.resume_from is not None),'focused preparation or accepted focused resume')
    cases,system,baseline,baseline_rows,manifest=load(base,args.resume_from)
    if args.focused:planned=preflight(base,cases,system)
    else:
        planned=json.loads((args.resume_from/'preflight.json').read_text());shutil.copyfile(args.resume_from/'preflight.json',base/'preflight.json')
        f.require(f.digest(base/'preflight.json')==f.digest(args.resume_from/'preflight.json'),'saved preflight identities')
    all_rows={BASELINE:baseline_rows};summaries={};new_nodes=0
    for label,case in cases.items():
        if args.focused and label!=BASELINE:continue
        target=base/'cases'/label;target.mkdir(parents=True,exist_ok=True)
        if label==BASELINE and not args.focused:
            prior=json.loads((args.resume_from/'checks.json').read_text());summaries[label]=prior['cases'][label];continue
        r,adapter,packets=context(base,label,case,manifest['input']);data=case['binding'];bound=data['bound']
        rows={v['row']:all_rows[v['fromCase']][v['fromRow']] for v in case['reusedRows']}
        selected=[bound['rows'][i] for i in case['newRows']];groups=[]
        if selected:
            worker=f.BasisMomentum(selected,bound['sources'],r)
            worker.prepare_basis(data['jets'],system['nodes'],system['settings']['sourceBound'],len(system['nodes']),system['settings']['sourceNodes'])
            width=float(bound['abel']['width'].subs(r.regulator,system['settings']['regulator']))
            for variables in sorted({tuple(l[0] for l in row['limits']) for row in selected},key=lambda v:(len(v),str(v))):
                work=target/('layout-'+str(len(groups)));work.mkdir()
                group=worker.matrix_group(variables,system['settings'],bound['pairs'],width,system['nodes'],work)
                groups.append(group);rows.update(zip(group['rowIndices'],group['matrices']));new_nodes+=group['nodes']
                f.save(target/'layout-inventory.json',{str(i):{'rowIndices':v['rowIndices'],'nodes':v['nodes'],
                    'sha256':f.digest(target/('layout-'+str(i))/f'layout-{len(v["variables"])}.pickle')} for i,v in enumerate(groups)})
        f.require(set(rows)==set(range(len(bound['rows']))),'complete actual case row arrays')
        for v in case['reusedRows']:f.require(np.array_equal(rows[v['row']],all_rows[v['fromCase']][v['fromRow']]),'unchanged reused matrix arrays')
        all_rows[label]=rows
        f.atomic_pickle(target/'row-matrices.pickle',{'rows':rows,'reusedRows':case['reusedRows'],'newGroups':groups,'settings':system['settings']})
        summary=assemble_case(target,case,packets,r,adapter,system,rows,baseline if label==BASELINE else None)
        summary.update(rows=len(rows),newRows=case['newRows'],newNodes=sum(v['nodes'] for v in groups),
            layoutCount=len(groups),reusedRows=len(case['reusedRows']))
        summaries[label]=summary;f.save(base/'case-inventory.json',summaries)
    for name,h in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==h,'unchanged current/frozen source')
    for name,h in manifest['inputPackets'].items():f.require(f.digest(Path(name))==h,'unchanged original input')
    for name,h in manifest['copiedInputs'].items():f.require(f.digest(base/name)==h,'unchanged copied input')
    artifacts={str(p.relative_to(base)):binding.factors.artifact(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
    checks={**manifest,'runDirectory':str(base),'mode':'focused' if args.focused else 'production',
        'status':'VALIDATED_ROW_REUSE_AND_ASSEMBLY' if args.focused else 'COMPLETED_FOUR_CASE_INTERIOR_MATRICES',
        'preflight':planned,'cases':summaries,'newNodes':new_nodes,'artifacts':artifacts,
        'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
