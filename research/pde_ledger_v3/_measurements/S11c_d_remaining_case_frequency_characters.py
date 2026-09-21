#!/usr/bin/env python3
"""Saved source-character routes and native row-frequency input census."""
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
import time

import S11c_d_remaining_case_frequency_sources as source

f,native,engine,sp,q=source.f,source.native,source.engine,source.sp,source.q
PLAN=f.M/'S11c_d_remaining_case_frequency_characters_plan.md'
CP=f.M/'S11c_d_remaining_case_frequency_sources_checkpoint.json'
SCOPE=('Saved Fourier-character binding reuse and full own-source/row input joins only. '
       'No new binding, differentiation, analytic lift, chart/end/mode/current, quadrature, '
       'matrix, solve, contour point or physical output. Numerical row reuse remains unclassified.')


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory']);v=Path(cp['validation']['runDirectory'])
    f.require(cp['status']=='ACCEPTED_CASE_FREQUENCY_SOURCE_RECORDS' and f.digest(origin/'checks.json')==cp['checksSha256'],
              'accepted complete native source records')
    source.receipts.inspect_guard(v,'validate')
    f.require(f.digest(v/'checks.json')==cp['validation']['checksSha256'] and
              (v/'checks.json').read_bytes()==(v/'validate.stdout').read_bytes(),'accepted saved record validator')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':dict(cp['inputPackets']),
              'referencedInputs':{},'input':cp['input'],'settings':cp['settings'],'scope':SCOPE,
              'acceptedSources':{'checkpoint':str(CP),'runDirectory':str(origin),'checksSha256':cp['checksSha256'],
                                 'validatorChecksSha256':cp['validation']['checksSha256']}}
    for name,item in cp['artifacts'].items():source.reference(base,manifest,origin/name,'accepted/'+name,item['sha256'])
    for path,name in ((CP,'accepted-source-checkpoint.json'),(origin/'checks.json','accepted-source-checks.json'),
                      (origin/'inputs.json','accepted-source-inputs.json'),(v/'checks.json','accepted-source-validation.json')):
        source.reference(base,manifest,path,name,f.digest(path))
    for path in (Path(__file__).resolve(),PLAN,CP):
        name=str(path.relative_to(f.ROOT));value=f.digest(path)
        f.require(name not in manifest['sourceFiles'] or manifest['sourceFiles'][name]==value,'exact additional helper pin')
        manifest['sourceFiles'][name]=value
    for name,value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==value,('current source',name))
        if name in cp['sourceFiles']:f.require(f.digest(origin/'source'/name)==value,'accepted frozen source')
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,dest)
        f.require(f.digest(dest)==value,'current frozen helper')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'input prehash')
    f.save(base/'inputs.json',manifest)
    return manifest,tuple(cp['cases'])


def row_constructor():
    tree=ast.parse(inspect.getsource(q.sources));fn=tree.body[0]
    loop=next(n for n in fn.body if isinstance(n,ast.For) and isinstance(n.target,ast.Name) and n.target.id=='row')
    changed=copy.deepcopy(loop)
    f.require(isinstance(changed.body[-1],ast.Expr) and ast.unparse(changed.body[-1].value.func)=='f.require',
              'native final source-column row guard')
    observer=ast.parse('save_row(base,row,row_census[-1],factors,used,source_frequencies,addresses,native)').body[0]
    changed.body.insert(-1,observer)
    reverse=copy.deepcopy(changed);del reverse.body[-2]
    f.require(ast.dump(reverse)==ast.dump(loop),'whole native row census reverse AST')
    args=ast.arguments(posonlyargs=[],args=[ast.arg(arg=n) for n in ('base','bound','addresses','source_frequencies','native')],
                       kwonlyargs=[],kw_defaults=[],defaults=[])
    wrapper=ast.FunctionDef(name='row_census',args=args,body=ast.parse('row_census=[]').body+[changed]+ast.parse('return row_census').body,decorator_list=[])
    module=ast.fix_missing_locations(ast.Module(body=[wrapper],type_ignores=[]));ns=dict(vars(q),save_row=save_row)
    exec(compile(module,'<native-saved-frequency-row-census>','exec'),ns)
    character=next(n for n in fn.body if isinstance(n,ast.For) and isinstance(n.target,ast.Name) and n.target.id=='si')
    return ns['row_census'],{'wholeNativeRowLoopReverseAST':True,'persistenceObservers':1,
        'nativeRowLoopSha256':hashlib.sha256(ast.dump(loop).encode()).hexdigest(),
        'nativeCharacterLoopSha256':hashlib.sha256(ast.dump(character).encode()).hexdigest(),
        'wholeNativeSourcesSha256':source.inputs.body(q.sources),
        'characterBindingInputs':'Only symbolicFrequency, complete saved live binder context and cutoff substitution; full physical aliases retained separately.',
        'newCharacterBindingCalls':0,'newDerivativeCalls':0}


def save_row(base,row,census,factors,used,characters,addresses,binding):
    target=base/'rows'/f'{row["index"]:03}.pickle';target.parent.mkdir(exist_ok=True)
    f.atomic_pickle(target,{'row':row,'census':census,'factorRecords':factors,
        'sourceRecords':{i:addresses['source',i] for i in used},'characters':{i:characters[i] for i in used},
        'sourceJets':{i:binding['jets'][i] for i in used}})


def prohibit():
    source.prohibit();source.inputs.prohibit()
    def forbidden(*args,**kwargs):raise RuntimeError('saved frequency-character input stage cannot construct scientific operands')
    source.load=source.restore=source.adapter=source.construct=forbidden
    sp.diff=sp.integrate=sp.lambdify=forbidden;sp.Basic.diff=forbidden
    native.factors.context=native.factors.cases.source.restore_context=forbidden


def signature(source_record,context,unit):
    return {'original':source_record['symbolicFrequency'],'accepted':source_record['frequency'],
            'unit':unit,'context':context}


def require_route(actual,expected):
    f.require(native.same(actual,expected),'literal character expression/unit/context/source address route')


def prepare(base,manifest,labels,row_run):
    accepted=base/'accepted';old=f.unpickle(accepted/'accepted-frequency/frequency-sources.pickle')
    baseline=f.unpickle(accepted/'accepted-frequency/baseline-binding.pickle')
    common=f.unpickle(accepted/'contexts'/native.BASELINE/'native-input.pickle')
    atlas=[];cases={};pending=[];joined=0;baseline_source_count=len(old['sourceFrequencies'])
    # The saved spatial source coordinate and native character equation fix
    # inverse-length units. Original context/unit/source packets remain raw.
    factor_base=f.unpickle(accepted/'accepted-cases'/native.BASELINE/'factorization.pickle')
    base_unit=tuple(factor_base['dimensionState']['known'][common['zp']])
    f.require(base_unit==(1,0,0),'actual accepted source-coordinate length declaration')
    character_unit=tuple(-v for v in base_unit)
    for si,value in old['sourceFrequencies'].items():
        raw=baseline['bound']['sources'][0,si]
        f.require(native.same((value['original'],value['accepted']),(raw['symbolicFrequency'],raw['frequency']))
                  and value['residual']==0,'baseline saved character producer/input/seed join')
        atlas.append({'owner':{'case':native.BASELINE,'sourceIndex':si,'kind':'accepted-baseline-character'},
                      'signature':signature(raw,common,character_unit),'value':value})
    f.atomic_pickle(base/'baseline-character-owners.pickle',atlas)
    for label in labels:
        target=base/'character-cases'/label;target.mkdir(parents=True)
        case=f.unpickle(accepted/'accepted-bindings'/label/'case-binding.pickle');b,g=case['binding'],case['grades'];bound=b['bound']
        context=f.unpickle(accepted/'contexts'/label/'native-input.pickle')
        packet=f.unpickle(accepted/'frequency-cases'/label/'frequency-source.pickle');records=packet['records']
        factors_packet=f.unpickle(accepted/'accepted-cases'/label/'factorization.pickle');fourier=factors_packet['result']
        f.atomic_pickle(target/'context-pair.pickle',(context,common))
        f.require(native.same(context,common),'complete original/live grade/frequency/field/profile/settings context')
        raw_units={'sourceCoordinate':context['zp'],'factorDimensionState':factors_packet['dimensionState'],
                   'recordDimensionState':packet['dimensionState'],'sourceCoordinateUnit':base_unit,'characterUnit':character_unit}
        f.atomic_pickle(target/'source-unit-inputs.pickle',raw_units)
        f.require(tuple(factors_packet['dimensionState']['known'][context['zp']])==base_unit,'own source-coordinate unit producer')
        row_inputs={'rows':bound['rows'],'sources':bound['sources'],'profiles':bound['profiles'],'profileUnits':bound['profileUnits'],
                    'jets':b['jets'],'abelSourceJoin':b['abelSourceJoin'],'settings':b['settings']}
        accepted_rows=f.unpickle(accepted/'cases'/label/'row-source-profile-inputs.pickle')
        f.atomic_pickle(target/'full-row-input-pair.pickle',(row_inputs,accepted_rows))
        f.require(native.same(row_inputs,accepted_rows),'whole actual original row/source/profile/Abel/settings packets')
        f.atomic_pickle(target/'term-inputs.pickle',g['termJoins'])
        f.require(native.same(g['termJoins'],f.unpickle(accepted/'frequency-cases'/label/'term-inputs.pickle')),'all actual independent-grade term convolution inputs')
        source_values={};routes={};raws={};uncomputed=[]
        for si,jet in b['jets'].items():
            raw=bound['sources'][0,si]
            uses=[(r['INDEX'],fi,v) for r in fourier['ROWS'] for fi,v in enumerate(r['FACTORS'])
                  if native.same(v['SOURCE_INTEGRAL'],raw['originalSourceIntegral'])]
            physical={'case':label,'sourceIndex':si,'source':raw,'jet':jet,'uses':uses,
                'allTests':{ti:s for (ti,j),s in bound['sources'].items() if j==si},
                'record':next(v for v in records.values() if v['address']==('source',si)),
                'rows':{r['index']:r for r in bound['rows'] if any(t['sourceIndex']==si for t in r['factors'])},
                'context':context,'unit':character_unit}
            raw_path=target/'source-inputs'/f'{si:03}.pickle';raw_path.parent.mkdir(exist_ok=True);f.atomic_pickle(raw_path,physical)
            f.require(uses and [(ri,fi) for ri,fi,_ in uses]==raw['uses'],'actual physical source-to-factor addresses')
            for _,_,v in uses:
                f.require(native.same((v['AMPLITUDE'],v['FREQUENCY']),(raw['symbolicAmplitude'],raw['symbolicFrequency']))
                    and v['CHARACTER_EQUATION_RESIDUAL']==v['CHARACTER_NORMALIZATION_RESIDUAL']==v['AMPLITUDE_RECONSTRUCTION_RESIDUAL']==0,
                    'complete saved actual Fourier character/amplitude equations')
            for test in physical['allTests'].values():
                f.require(native.same((test['symbolicFrequency'],test['frequency'],test['symbolicAmplitude']),
                           (raw['symbolicFrequency'],raw['frequency'],raw['symbolicAmplitude'])),'all Gaussian-test character/source addresses')
            rec=physical['record']
            f.require(native.same((rec['original'],tuple(rec['unit'])),(raw['symbolicAmplitude'],tuple(raw['amplitudeUnit'])))
                and native.same(rec['acceptedBinding'],jet['originalBoundAmplitude']),'actual source-record/jet/column/amplitude join')
            sig=signature(raw,context,character_unit)
            # A scalar character result may be shared across different amplitudes;
            # every source/row/domain alias above stays its own physical input.
            owner=next((v for v in atlas if native.same(v['signature'],sig)),None)
            route={'case':label,'sourceIndex':si,'signature':sig,'owner':None if owner is None else owner['owner'],
                   'physicalSourceInput':str(raw_path),'sourceColumn':jet['column'],'numericalRowReuse':False}
            raws[si]=physical;routes[si]=route
            if owner is None:
                uncomputed.append(si);pending.append(route)
            else:
                source_values[si]=owner['value'];joined+=1
        f.atomic_pickle(target/'character-routes.pickle',routes)
        f.atomic_pickle(target/'source-frequencies.pickle',source_values)
        f.atomic_pickle(target/'uncomputed-character-inputs.pickle',[routes[i] for i in uncomputed])
        # Saved actual operands test the same equality used to accept a route.
        si=next(iter(routes));route=routes[si];changed=[]
        for key,value in (('sourceIndex',si+10000),('sourceColumn',route['sourceColumn']+1)):
            v=copy.deepcopy(route);v[key]=value;changed.append((key,v))
        for key,value in (('original',route['signature']['original']+1),('accepted',route['signature']['accepted']+1),
                          ('unit',(character_unit[0]+1,*character_unit[1:]))):
            v=copy.deepcopy(route);v['signature'][key]=value;changed.append((key,v))
        v=copy.deepcopy(route);v['signature']['context']['parameters']['omega']+=1;changed.append(('frequency',v))
        v=copy.deepcopy(route);v['signature']['context']['settings']['sourceBound']+=1;changed.append(('sourceBound',v))
        row=bound['rows'][0];limit=row['limits'][0];wrong_limit=sp.Tuple(limit[0],limit[1],limit[2]+1)
        f.atomic_pickle(target/'mutation-operands.pickle',{'original':route,'changed':changed,'limit':limit,'changedLimit':wrong_limit})
        controls={name:source.inputs.source.rejects(lambda bad=v:require_route(bad,route)) for name,v in changed}
        controls['orderedLimit']=not native.same(limit,wrong_limit)
        f.save(target/'mutation-controls.json',controls);f.require(all(controls.values()),'actual character/address/unit/context/limit mutations reject')
        addresses={tuple(v['address']):v for v in records.values()}
        if not uncomputed:
            if label==native.BASELINE:
                f.require(native.same(bound['rows'],baseline['bound']['rows']) and native.same(b['jets'],baseline['jets'])
                          and native.same(source_values,old['sourceFrequencies']),'full baseline completed row/jet/character reuse')
                census=old['rowCensus'];f.atomic_pickle(target/'baseline-row-census-reuse.pickle',census)
            else:census=row_run(target,bound,addresses,source_values,b)
            f.atomic_pickle(target/'row-census.pickle',census)
            complete=True
        else:census=None;complete=False
        summary={'records':len(records),'rows':len(bound['rows']),'terms':len(g['termJoins']),'sources':len(b['jets']),
            'savedCharacterUses':len(source_values),'uncomputedCharacterBindings':len(uncomputed),'rowCensusComplete':complete,
            'frequencyDependentRows':None if census is None else [v['index'] for v in census if v['frequencyDependent']],
            'characterFrequencyDependence':{str(i):v['dependsOnFrequency'] for i,v in source_values.items()},
            'controls':controls,'numericalRowReuseAccepted':False}
        result={'sourceFrequencies':source_values,'routes':routes,'rowCensus':census,'summary':summary,
                'frequency':old['frequency'],'referenceFrequency':old['referenceFrequency'],'origin':old['origin'],
                'scope':SCOPE,'analyticDomainAndEndMapsPending':True}
        f.atomic_pickle(target/'frequency-characters.pickle',result);cases[label]=summary
        f.save(base/'character-case-inventory.json',cases)
        del case,b,g,bound,records,raws,physical,source_values,routes,result,factors_packet,fourier
        gc.collect()
    f.atomic_pickle(base/'uncomputed-character-operands.pickle',pending)
    f.require(sum(v['records'] for v in cases.values())==1467 and sum(v['rows'] for v in cases.values())==300
              and sum(v['terms'] for v in cases.values())==647 and sum(v['sources'] for v in cases.values())==120,'complete actual four-case inputs')
    return {'cases':cases,'baselineCharacters':baseline_source_count,'savedCharacterUses':joined,'uncomputedCharacterBindings':len(pending),
            'newBindings':0,'newDerivatives':0,'newNumericalWork':0,'sourceRecordAcceptance':str(CP)}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',required=True,type=Path);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    manifest,labels=load(base);row_run,join=row_constructor();f.save(base/'native-character-row-joins.json',join)
    prohibit();result=prepare(base,manifest,labels,row_run)
    for name,value in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(base/'source'/name)==value,'current/frozen posthash')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'input pre/posthash')
    for name,item in manifest['referencedInputs'].items():
        p=base/name;f.require(p.is_symlink() and str(p.readlink())==item['original'] and str(p.resolve())==item['resolvedOriginal']
                           and p.stat().st_size==item['bytes'] and f.digest(p)==item['sha256'],'literal original/reference identity')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*')
               if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,**result,'status':'COMPLETED_CASE_FREQUENCY_CHARACTER_INPUTS','nativeJoins':join,'artifacts':artifacts,
            'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
