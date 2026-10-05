#!/usr/bin/env python3
"""Complete finite-action momentum quadrature from accepted source transforms."""
import argparse
import ast
import contextlib
import json
from pathlib import Path
import resource
import shutil
import time

import numpy as np
import sympy as sp
from sympy.core.function import AppliedUndef
import S11c_d_source_fourier_quadrature_check as source
from S11c_d_source_fourier_quadrature_check import ROOT,STORE,engine,digest,save,atomic_pickle,unpickle
from S11c_d_numerical_action_check import number,progress
from S11c_d_output_codec import decoded_lines,restore_emission_index
from ledger_fold import _restore

M=ROOT/'_measurements'
SOURCE_CHECKPOINT=M/'S11c_d_source_fourier_quadrature_checkpoint.json'
DOMAIN_CHECKPOINT=M/'S11c_d_quadrature_domain_checkpoint.json'
PLAN=M/'S11c_d_momentum_action_plan.md'
PREFIX='MOMENTUM_ACTION_LAB_HELD_RHO4_CONSTANT'


def load():
    accepted,base=source.accepted(SOURCE_CHECKPOINT)
    domains,domain_base=source.accepted(DOMAIN_CHECKPOINT)
    engine_name=str(engine.HERE.relative_to(ROOT))
    for name,sha in accepted['sourceFiles'].items():
        if digest(base/'source'/name)!=sha:raise ValueError(('source quadrature snapshot changed',name))
        if name!=engine_name and digest(ROOT/name)!=sha:raise ValueError(('consumed quadrature helper changed',name))
    old=ast.parse((base/'source'/engine_name).read_text());new=ast.parse(engine.HERE.read_text())
    cls=next(n for n in new.body if getattr(n,'name',None)=='BoundedSourceFourierQuadrature')
    additions=[n for n in cls.body if getattr(n,'name',None)=='FiniteMomentum']
    if len(additions)!=1:raise ValueError('expected one finite momentum helper')
    cls.body.remove(additions[0])
    if ast.dump(old)!=ast.dump(new):raise ValueError('engine changed beyond finite momentum helper')
    pencil,assembly,factors,tests,numerical,provenance,joins=source.load()
    bound=unpickle(base/'bound-sources.pickle')
    if bound['provenance']!=accepted['provenance'] or provenance!=accepted['provenance']:
        raise ValueError('accepted source quadrature provenance differs')
    engine.PHYSICAL_METADATA.dimensions.__dict__.update(bound['dimensionState'])
    domain=unpickle(domain_base/'domains.pickle')
    if domain['provenance']!=domains['provenance']:raise ValueError('domain packet provenance')
    denominator=unpickle(domain_base/'denominators.pickle')
    pairs=denominator['pairs'];abel=domain['result']['abel']
    if (abel['primitiveResidual']!=0 or abel['widthResidual']!=0 or
            domain['provenance']['NUMERICAL_ACTION_CHECKPOINT_SHA256']!=digest(source.NUMERICAL)):
        raise ValueError('source-derived width/source join')
    provenance=dict(provenance,SOURCE_QUADRATURE_CHECKPOINT_SHA256=digest(SOURCE_CHECKPOINT),
        DOMAIN_CHECKPOINT_SHA256=digest(DOMAIN_CHECKPOINT),BOUND_SOURCE_PACKET_SHA256=digest(base/'bound-sources.pickle'))
    return pencil,assembly,factors,tests,numerical,bound['result'],domain['result'],pairs,provenance,joins,accepted


def bind(base,pencil,assembly,factors,tests,numerical,bound,domains,pairs,provenance):
    adapter=engine.NumericalReducedAction(pencil,assembly,json.loads(source.native.INPUT.read_text()))
    r=pencil.r;settings=numerical['results'][0]['settings']
    cutoffs={symbol:sp.Rational(str(settings['profileBound'] if variable==r.xi else
              settings['sourceBound'] if variable==r.zp else settings['momentumBound'])) for variable,symbol in factors['CUTOFFS'].items()}
    rows=[];inventory=[]
    for row in factors['ROWS']:
        record={'index':row['INDEX'],'original':row['ORIGINAL'],'symbolicLimits':row['REMAINING_LIMITS'],
            'limits':tuple(engine.memo_xreplace(limit,cutoffs) for limit in row['REMAINING_LIMITS']),
            'sourceLimit':engine.memo_xreplace(row['SOURCE_LIMIT'],cutoffs),
            'unit':engine.PHYSICAL_METADATA.dimensions.measure(row['ORIGINAL']),'factors':[]}
        for f in row['FACTORS']:
            si=factors['SOURCE_INTEGRALS'].index(f['SOURCE_INTEGRAL'])
            coefficient=engine.memo_xreplace(adapter.bind(f['COEFFICIENT']),cutoffs)
            allowed={r.z,r.regulator,*(lim[0] for lim in record['limits'])}
            if engine.dag_free_symbols(coefficient)-allowed or coefficient.has(AppliedUndef,sp.Derivative,sp.Subs):
                raise ValueError(('unbound momentum coefficient',row['INDEX']))
            record['factors'].append({'sourceIndex':si,'symbolicCoefficient':f['COEFFICIENT'],
                'coefficient':coefficient,'unit':engine.PHYSICAL_METADATA.dimensions.measure(f['COEFFICIENT'])})
        for limit in record['limits']:
            if len(limit)!=3 or tuple(limit[1:])!=(-sp.Rational(str(settings['momentumBound'])),sp.Rational(str(settings['momentumBound']))):
                raise ValueError(('finite momentum limit differs from quadrature domain',row['INDEX']))
        rows.append(record)
        path=base/'rows'/f'{row["INDEX"]:03}.pickle';path.parent.mkdir(exist_ok=True)
        atomic_pickle(path,record);inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)})
        save(base/'row-inventory.json',inventory)
    sources={}
    for record in bound['records']:
        ti,si=record['test'],record['sourceIndex']
        if record['originalSourceIntegral']!=factors['SOURCE_INTEGRALS'][si]:raise ValueError('source/factor join')
        sources[(ti,si)]=dict(record,testWidth=tests[ti]['width'],profileWidth=bound['profileWidth'])
    profile_units={sp.Integral(p['bound'].function,(r.xi,-sp.Rational(str(settings['profileBound'])),
                       sp.Rational(str(settings['profileBound'])))):p['unit'] for p in domains['profiles']}
    actual_profiles=set().union(*(f['coefficient'].atoms(sp.Integral) for row in rows for f in row['factors']))
    atomic_pickle(base/'profile-join.pickle',{'actual':actual_profiles,'expected':profile_units})
    if actual_profiles!=set(profile_units):raise ValueError('complete nested profile operand join')
    row_ids={row['ORIGINAL']:row['INDEX'] for row in factors['ROWS']}
    cells=[]
    for ti,test in enumerate(tests):
        for original in assembly['ROWS']:
            j,i=original['COLUMN'],original['ROW'];cached=next(v for v in test['assembled'] if (v['COLUMN'],v['ROW'])==(j,i))
            if len(original['NONLOCAL'])!=len(cached['NONLOCAL']):raise ValueError('action term census')
            cells.append({'test':ti,'column':j,'row':i,'local':cached['LOCAL'],
                'terms':[(row_ids[g],cached['NONLOCAL'][n][0]) for n,(g,_) in enumerate(original['NONLOCAL'])]})
    source_checkpoint,sb=source.accepted(source.SOURCE);actions=unpickle(sb/'actions.pickle')
    if set(sources)!={(ti,si) for ti in range(2) for si in range(len(factors['SOURCE_INTEGRALS']))}:
        raise ValueError('incomplete bound source census')
    if [row['index'] for row in rows]!=list(range(len(factors['ROWS']))):raise ValueError('factor row index census')
    result={'rows':rows,'sources':sources,'cells':cells,'positions':numerical['positions'],'tests':numerical['tests'],
        'pairs':pairs,'abel':domains['abel'],'profiles':domains['profiles'],'profileUnits':profile_units,
        'rowInventory':inventory,'referenceResults':numerical['results'],
        'equationUnits':[actions['columnUnits'][(0,i)] for i in range(5)],
        'momentumUnit':engine.PHYSICAL_METADATA.dimensions.measure(bound['momenta'][0]),
        'cutoffBindings':cutoffs,'provenance':provenance}
    atomic_pickle(base/'bound-momentum.pickle',{'result':result,'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    return result


def settings(result):
    output=[]
    for value in result['referenceResults']:
        if value['test']==0:output.append(dict(value['settings'],kind='legacy',referenceGrid=value['grid']))
    common=result['referenceResults'][0]['settings']
    for panel,outer in ((8,24),(12,40),(16,64)):
        output.append(dict(common,kind='concentrated',panelOrder=panel,outerOrder=outer,sourceNodes=128,profileNodes=128))
    return output


def evaluate(base,pencil,result):
    contractor=engine.BoundedSourceFourierQuadrature.FiniteMomentum(result['rows'],result['sources'],pencil.r)
    variables=tuple(dict.fromkeys(tuple(l[0] for l in row['limits']) for row in result['rows']))
    records=[];inventory=[];partial_inventory=[]
    for gi,setting in enumerate(settings(result)):
        width=float(result['abel']['width'].subs(pencil.r.regulator,setting['regulator']))
        if not np.isfinite(width) or width<=0:raise ValueError('invalid computed Abel width')
        for ti in range(2):
            integral_values={};groups=[]
            for li,layout in enumerate(variables):
                def batch_progress(state):
                    if state['batchCount']%64==0:
                        path=base/'partials'/f'{gi}-{ti}-{li}-{state["batchCount"]:08}.pickle'
                        path.parent.mkdir(exist_ok=True)
                        atomic_pickle(path,dict(state,boundPacketSha256=digest(base/'bound-momentum.pickle')))
                        partial_inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)})
                        save(base/'partial-inventory.json',partial_inventory)
                        progress(base,'momentum_batches',test=ti,grid=gi,layout=li,batch=state['batchCount'],nodes=state['nodeCount'])
                group=contractor.group(ti,layout,setting,result['pairs'],width,result['positions'],batch_progress)
                path=base/'groups'/f'{gi}-{ti}-{li}.pickle';path.parent.mkdir(exist_ok=True)
                atomic_pickle(path,group);inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)})
                save(base/'group-inventory.json',inventory)
                if abs(group['volumeResidual'])>1e-10*(1+abs(group['boxVolume'])):raise ValueError('momentum quadrature box mass')
                if np.any((np.abs(group['values'])>1e-9)&(np.abs(group['measureMutationResidual'])<1e-12)):
                    raise ValueError('momentum measure mutation lost a resolved nonzero contribution')
                for i,index in enumerate(group['rowIndices']):integral_values[index]=group['values'][i]
                groups.append(group);progress(base,'momentum_group_saved',test=ti,grid=gi,layout=li,nodes=group['nodeCount'])
            if set(integral_values)!=set(range(len(result['rows']))):raise ValueError('incomplete momentum row coverage')
            local=np.zeros((len(result['positions']),5,5),dtype=complex);nonlocal_values=np.zeros_like(local);terms=[]
            compiler=engine.BoundedActionQuadrature({})
            for cell in (c for c in result['cells'] if c['test']==ti):
                j,i=cell['column'],cell['row']
                for pi,z in enumerate(result['positions']):
                    env={pencil.r.z:z,pencil.r.regulator:setting['regulator']}
                    local[pi,j,i]=complex(compiler.evaluate(cell['local'],env))
                    values=[complex(compiler.evaluate(c,env))*integral_values[index][pi] for index,c in cell['terms']]
                    nonlocal_values[pi,j,i]=sum(values)
                    terms.append({'index':(pi,j,i),'terms':values})
            item={'test':ti,'grid':gi,'setting':setting,'width':width,'groups':groups,
                'local':local,'nonlocal':nonlocal_values,'action':local+nonlocal_values,'contributions':terms}
            if setting['kind']=='legacy':
                old=next(r for r in result['referenceResults'] if (r['test'],r['grid'])==(ti,setting['referenceGrid']))
                item['nativeDirect']=old['arrays']['direct'];item['nativeAssembled']=old['arrays']['local']+old['arrays']['nonlocal']
                item['nativeDirectResidual']=item['action']-item['nativeDirect']
                item['nativeAssembledResidual']=item['action']-item['nativeAssembled']
                joined={tuple(r['index']):r['terms'] for r in old['contributions']}
                item['termResiduals']=[{'index':r['index'],'residual':np.asarray(r['terms'])-np.asarray(joined[tuple(r['index'])])} for r in terms]
            path=base/'grids'/f'{gi}-{ti}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,item)
            inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)});save(base/'group-inventory.json',inventory)
            records.append(item);progress(base,'momentum_grid_saved',test=ti,grid=gi)
            if setting['kind']=='legacy':
                scale=1+np.abs(item['nativeDirect'])
                if (np.max(np.abs(item['nativeDirectResidual'])/scale)>1e-9 or
                        max(float(np.max(np.abs(r['residual']))) if len(r['residual']) else 0. for r in item['termResiduals'])>1e-9):
                    raise ValueError('factorized/native whole-action common-grid residual')
    if set(contractor.source_frequency_census)!=set(result['sources']) or contractor.profile_integrals!=set(result['profileUnits']):
        raise ValueError('incomplete source/profile evaluation census')
    changes=[]
    for ti in range(2):
        refined=[r for r in records if r['test']==ti and r['setting']['kind']=='concentrated']
        changes.extend({'test':ti,'fromGrid':a['grid'],'toGrid':b['grid'],'difference':b['action']-a['action']}
                       for a,b in zip(refined,refined[1:]))
    output={'records':records,'refinements':changes,'groupArtifacts':inventory,'partialArtifacts':partial_inventory,
        'sourceFrequencyCensus':contractor.source_frequency_census,
        'evaluatedProfileIntegrals':tuple(sorted(contractor.profile_integrals,key=sp.default_sort_key)),
        'phaseWorkspaceBudgetBytes':contractor.workspace_bytes,'batchCacheBudgetBytes':contractor.workspace_bytes,
        'peakPhaseWorkspaceEstimateBytes':contractor.peak_phase_workspace_estimate,
        'peakBatchCacheEstimateBytes':contractor.peak_batch_cache_estimate,
        'peakWorkspaceEstimateBytes':contractor.peak_workspace_estimate}
    atomic_pickle(base/'momentum-action.pickle',{'result':output,'boundPacketSha256':digest(base/'bound-momentum.pickle'),
        'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    return output


def emit_bound(bound,pencil,provenance):
    metadata=engine.FullPencilModes.__new__(engine.FullPencilModes);metadata.r=pencil.r
    zero=engine.PHYSICAL_METADATA.dimensions.zero
    def numeric(name,value,unit,literal=False):
        body=number(value)
        engine.emit(PREFIX+'_'+name,body if literal else metadata.compact_fingerprint(body))
        engine.emit('METADATA_'+PREFIX+'_'+name,metadata.numeric_metadata(body,unit))
    engine.physical(PREFIX+'_PROVENANCE',provenance)
    engine.physical(PREFIX+'_ABEL_TRANSFER_PAIRS',bound['pairs'])
    engine.physical(PREFIX+'_SOURCE_DERIVED_ABEL_WIDTH',bound['abel']['width'])
    numeric('TEST_POSITIONS',bound['positions'],lambda p:engine.PHYSICAL_METADATA.dimensions.measure(pencil.r.z),True)
    numeric('TEST_WIDTHS_MOMENTA',bound['tests'],lambda p:tuple((1 if p[1]==0 else -1)*v for v in engine.PHYSICAL_METADATA.dimensions.measure(pencil.r.z)),True)
    for row in bound['rows']:
        index=str(row['index']);engine.fingerprinted(PREFIX+'_ORIGINAL_INTEGRAL_'+index,row['original'])
        engine.physical(PREFIX+'_SYMBOLIC_REMAINING_LIMITS_'+index,row['symbolicLimits'])
        numeric('FINITE_REMAINING_LIMITS_'+index,row['limits'],lambda p:bound['momentumUnit'],True)
        for fi,f in enumerate(row['factors']):
            name=index+'_'+str(fi)
            engine.fingerprinted(PREFIX+'_SYMBOLIC_COEFFICIENT_'+name,f['symbolicCoefficient'])
            engine.emit(PREFIX+'_BOUND_COEFFICIENT_'+name,engine.carrier_fingerprint(f['coefficient']))
            engine.emit('METADATA_'+PREFIX+'_BOUND_COEFFICIENT_'+name,metadata.numeric_metadata(f['coefficient'],lambda p:f['unit']))
            numeric('SOURCE_INDEX_'+name,f['sourceIndex'],lambda p:zero,True)


def emit(result,bound,pencil,provenance):
    emit_bound(bound,pencil,provenance)
    metadata=engine.FullPencilModes.__new__(engine.FullPencilModes);metadata.r=pencil.r
    zero=engine.PHYSICAL_METADATA.dimensions.zero
    def numeric(name,value,unit,literal=False):
        body=number(value)
        engine.emit(PREFIX+'_'+name,body if literal else metadata.compact_fingerprint(body))
        engine.emit('METADATA_'+PREFIX+'_'+name,metadata.numeric_metadata(body,unit))
    for item in result['records']:
        key=str(item['grid'])+'_'+str(item['test'])
        numeric('COMPUTED_WIDTH_'+key,item['width'],lambda p:bound['momentumUnit'],True)
        numeric('SOURCE_BOUNDS_'+key,[-item['setting']['sourceBound'],item['setting']['sourceBound']],lambda p:engine.PHYSICAL_METADATA.dimensions.measure(pencil.r.zp),True)
        numeric('PROFILE_BOUNDS_'+key,[-item['setting']['profileBound'],item['setting']['profileBound']],lambda p:engine.PHYSICAL_METADATA.dimensions.measure(pencil.r.xi),True)
        numeric('MOMENTUM_BOUNDS_'+key,[-item['setting']['momentumBound'],item['setting']['momentumBound']],lambda p:bound['momentumUnit'],True)
        numeric('FINITE_ABEL_REGULATOR_'+key,item['setting']['regulator'],lambda p:engine.PHYSICAL_METADATA.dimensions.measure(pencil.r.regulator),True)
        orders={k:v for k,v in item['setting'].items() if k.endswith('Nodes') or k.endswith('Order')}
        if item['setting']['kind']!='legacy':orders.pop('momentumNodes')
        numeric('GAUSS_ORDERS_'+key,orders,lambda p:zero,True)
        engine.physical(PREFIX+'_RULE_KIND_'+key,item['setting']['kind'])
        for name in ('local','nonlocal','action','nativeDirect','nativeAssembled','nativeDirectResidual','nativeAssembledResidual'):
            if name in item:numeric(name.upper()+'_'+key,item[name],lambda p:bound['equationUnits'][p[2]],name.endswith('Residual'))
        for li,group in enumerate(item['groups']):
            suffix=key+'_'+str(li)
            numeric('INTEGRAL_ROW_INDICES_'+suffix,group['rowIndices'],lambda p:zero,True)
            engine.physical(PREFIX+'_ORDERED_VARIABLES_'+suffix,group['variables'])
            for name in ('values','measureMutationValues','measureMutationResidual'):
                numeric(name.upper()+'_'+suffix,group[name],lambda p:bound['rows'][group['rowIndices'][p[0]]]['unit'],name.endswith('Residual'))
            volume_unit=tuple(len(group['variables'])*v for v in bound['momentumUnit'])
            numeric('VOLUME_OPERANDS_'+suffix,(group['quadratureMass'],group['boxVolume'],group['volumeResidual']),lambda p:volume_unit,True)
            numeric('RULE_COUNTS_'+suffix,(group['nodeCount'],group['batchCount']),lambda p:zero,True)
        for ni,record in enumerate(item['contributions']):
            numeric('CELL_INDEX_'+key+'_'+str(ni),record['index'],lambda p:zero,True)
            numeric('NONLOCAL_TERMS_'+key+'_'+str(ni),record['terms'],lambda p:bound['equationUnits'][record['index'][2]])
        for ni,record in enumerate(item.get('termResiduals',[])):
            numeric('NATIVE_TERM_RESIDUAL_'+key+'_'+str(ni),record['residual'],lambda p:bound['equationUnits'][record['index'][2]],True)
    for i,change in enumerate(result['refinements']):
        numeric('REFINEMENT_TEST_GRIDS_'+str(i),(change['test'],change['fromGrid'],change['toGrid']),lambda p:zero,True)
        numeric('MOMENTUM_REFINEMENT_'+str(i),change['difference'],lambda p:bound['equationUnits'][p[2]],True)
    numeric('WORKSPACE_BUDGETS_AND_ESTIMATE',(result['phaseWorkspaceBudgetBytes'],result['batchCacheBudgetBytes'],result['peakWorkspaceEstimateBytes']),lambda p:zero,True)
    numeric('PHASE_CACHE_ESTIMATES',(result['peakPhaseWorkspaceEstimateBytes'],result['peakBatchCacheEstimateBytes']),lambda p:zero,True)
    for (ti,si),entry in sorted(result['sourceFrequencyCensus'].items()):
        numeric('SOURCE_FREQUENCY_CENSUS_'+str(ti)+'_'+str(si),entry,lambda p:zero if p[0]==2 else bound['momentumUnit'],True)
    for i,integral in enumerate(result['evaluatedProfileIntegrals']):
        engine.emit(PREFIX+'_EVALUATED_PROFILE_'+str(i),engine.carrier_fingerprint(integral))
        engine.emit('METADATA_'+PREFIX+'_EVALUATED_PROFILE_'+str(i),metadata.numeric_metadata(integral,
            lambda p:bound['profileUnits'][integral]))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    pencil,assembly,factors,tests,numerical,sources,domains,pairs,provenance,joins,accepted=load()
    paths=tuple(dict.fromkeys((*source.SOURCES,*(ROOT/name for name in accepted['sourceFiles']),
        SOURCE_CHECKPOINT,DOMAIN_CHECKPOINT,PLAN,Path(__file__).resolve(),M/'S11c_d_momentum_action_preflight.py')))
    pins={str(p.relative_to(ROOT)):digest(p) for p in paths}
    for p in paths:
        target=base/'source'/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    save(base/'preflight.json',{'sourceFiles':pins,'provenance':provenance,'nativeEngineAstJoin':True,'originalLimitJoins':joins})
    bound=bind(base,pencil,assembly,factors,tests,numerical,sources,domains,pairs,provenance)
    before_bound=digest(base/'bound-momentum.pickle');progress(base,'bound_momentum_saved',rows=len(bound['rows']))
    # Exercise every actual operand's emitter before expensive integration.
    with (base/'bound-operands.out').open('x') as out,contextlib.redirect_stdout(out):
        emit_bound(bound,pencil,provenance)
    if engine.PHYSICAL_METADATA.dimensions.constraints:raise ValueError('bound operand dimensions')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=engine.PayloadEncoder()
    progress(base,'bound_metadata_saved')
    result=evaluate(base,pencil,bound);before=digest(base/'momentum-action.pickle')
    with (base/'full.out').open('x') as out,contextlib.redirect_stdout(out):
        emit(result,bound,pencil,provenance)
        keys={tag:'s11cd'+''.join(w.title() for w in tag.removeprefix('PY_S11CD_').split('_')) for tag in engine.EMISSION_LINES if not tag.startswith('PY_S11CD_METADATA_')}
        source.emit_manifest(PREFIX+'_WRITE_KEYS',keys)
        index=engine.emission_index(engine.EMISSION_LINES);units={p:(0,0,0) for p,_ in engine.leaves(engine.cas(index))}
        source.emit_manifest(PREFIX+'_EMISSION_LINES',index,zero_dimensions=units)
    entries={}
    for line in decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ')
        if tag in entries:raise ValueError('duplicate momentum emission')
        entries[tag]=_restore(body)
    seen=set();old_emit=engine.emit
    def compare(name,value):
        tag='PY_S11CD_'+name
        if tag in seen or entries.get(tag)!=engine.cas(value):raise ValueError(('momentum replay',tag))
        seen.add(tag)
    engine.emit=compare
    try:
        emit(result,bound,pencil,provenance);source.emit_manifest(PREFIX+'_WRITE_KEYS',keys)
        source.emit_manifest(PREFIX+'_EMISSION_LINES',index,zero_dimensions=units)
    finally:engine.emit=old_emit
    if seen!=set(entries) or len(keys)!=len(set(keys.values())) or set(keys.values())&set(engine.IMPORT_KEYS):raise ValueError('momentum emission/key census')
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES';restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    paths_count=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        for record in body:
            if len(record)==2 and isinstance(record[0],sp.Tuple):fields={str(k):v for k,v in record[1]};count=1
            else:fields={str(k):v for k,v in record};count=len(fields['PATHS'])
            unit=fields['DIMENSION_L_T_M']
            if len(unit)!=3 or any(v.free_symbols for v in unit) or any(k not in fields for k in ('MULTIGRADE','EPSILON_LAMBDA_SUPPORT')):raise ValueError('momentum metadata')
            paths_count+=count
    summary={'runDirectory':str(base),'sourceFiles':pins,'provenance':provenance,'nativeEngineAstJoin':True,'originalLimitJoins':joins,
        'rows':len(bound['rows']),'sources':len(bound['sources']),'transferPairs':[[str(k) for k in pair] for pair in pairs],
        'gridRecords':len(result['records']),'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':paths_count,
        'legacyMaxScaledResiduals':[float(np.max(np.abs(v['nativeDirectResidual'])/(1+np.abs(v['nativeDirect'])))) for v in result['records'] if 'nativeDirect' in v],
        'refinementNorms':[float(np.max(np.abs(v['difference']))) for v in result['refinements']],
        'boundPacketSha256BeforeEmission':before_bound,'boundPacketSha256AfterEmission':digest(base/'bound-momentum.pickle'),
        'packetSha256BeforeEmission':before,'packetSha256AfterEmission':digest(base/'momentum-action.pickle'),
        'rowArtifacts':bound['rowInventory'],'groupArtifacts':result['groupArtifacts'],
        'partialArtifacts':result['partialArtifacts'],'evaluatedProfileCount':len(result['evaluatedProfileIntegrals']),
        'evaluatedSourceCount':len(result['sourceFrequencyCensus']),
        'legacyComponentCount':sum(v['action'].size for v in result['records'] if 'nativeDirect' in v),
        'legacyTermCount':sum(len(t['terms']) for v in result['records'] if 'nativeDirect' in v for t in v['contributions']),
        'legacyMaxTermResiduals':[max(float(np.max(np.abs(t['residual']))) if len(t['residual']) else 0. for t in v['termResiduals']) for v in result['records'] if 'nativeDirect' in v],
        'measureMutationMaxima':[[float(np.max(np.abs(g['measureMutationResidual']))) for g in v['groups']] for v in result['records']],
        'phaseWorkspaceBudgetBytes':result['phaseWorkspaceBudgetBytes'],'batchCacheBudgetBytes':result['batchCacheBudgetBytes'],
        'peakPhaseWorkspaceEstimateBytes':result['peakPhaseWorkspaceEstimateBytes'],'peakBatchCacheEstimateBytes':result['peakBatchCacheEstimateBytes'],
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'scope':'Complete finite-source/profile/momentum action at positive regulator, common-grid native comparisons and recorded momentum refinements. No infinite-domain interchange, full-action tail or Abel weak limit, matching, scattering or pole solve.'}
    save(base/'checks.json',summary)
    for item in (*bound['rowInventory'],*result['groupArtifacts'],*result['partialArtifacts']):
        if digest(base/item['path'])!=item['sha256']:raise ValueError('momentum saved operand changed')
    if pins!={str(p.relative_to(ROOT)):digest(p) for p in paths} or before_bound!=digest(base/'bound-momentum.pickle') or before!=digest(base/'momentum-action.pickle') or engine.PHYSICAL_METADATA.dimensions.constraints:raise ValueError('momentum source/packet/dimension final guard')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
