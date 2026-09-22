#!/usr/bin/env python3
"""Read full actual frequency-row inputs and saved source-basis/row returns."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time

import numpy as np
import sympy as sp
import S11c_d_remaining_case_frequency_end_continuation_inputs as io

M,REPO=io.M,io.REPO
CP=M/'S11c_d_remaining_case_frequency_source_coefficients_checkpoint.json'
SHA='4330b4ed945b3cf6f64f74e08e850ee33676fdce163945a899e832c441370782'
require,same=io.require,io.same
BASELINE='LAB_HELD__RHO4_CONSTANT'


def signature(row,binding):
    """Literal native consumed fields plus full settings/field/measure context."""
    bound=binding['bound'];jets=binding['jets'];factors=[];used_profiles=set()
    for factor in row['factors']:
        source=bound['sources'][0,factor['sourceIndex']];jet=jets[factor['sourceIndex']]
        factors.append((factor['symbolicCoefficient'],factor['coefficient'],factor['unit'],
            source['originalSourceIntegral'],source['symbolicAmplitude'],source['frequency'],
            jet['column'],jet['probe'],tuple(jet['coefficients']),jet['amplitudeUnit'],jet['integralUnit'],
            source['boundCharacter']))
        used_profiles.update(factor['coefficient'].atoms(sp.Integral))
    require(used_profiles<=set(bound['profileUnits']),'complete consumed profile unit operands')
    return {'row':(row['original'],row['symbolicLimits'],row['limits'],row['sourceLimit'],row['unit'],tuple(factors)),
            'fieldUnits':binding['fieldUnits'],'equationUnits':binding.get('equationUnits',binding.get('rowUnits')),
            'profileUnits':{p:bound['profileUnits'][p] for p in used_profiles},'abel':bound['abel'],
            'settings':binding['settings'],'pairs':bound['pairs']}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path)
    base=ap.parse_args().run_directory.resolve();base.relative_to(REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic();reader=io.Reader();cache={}
    cp=reader.json(CP,SHA);require(cp['status']=='ACCEPTED_CASE_FREQUENCY_SOURCE_COEFFICIENTS','accepted complete source coefficients');origin=Path(cp['runDirectory'])
    reader.retain(origin/'checks.json',cp['checksSha256']);scalar_ref=cp['upstreamCheckpoint'];scalar_cp=reader.json(scalar_ref['logical'],scalar_ref['sha256']);sr=Path(scalar_cp['runDirectory'])
    reader.retain(sr/'checks.json',scalar_cp['checksSha256'])
    def load(record):
        p=record.get('logical',record.get('path'));r=reader.retain(p,record['sha256'])
        if r['canonical'] not in cache:cache[r['canonical']]=reader.packet(p,record['sha256'])
        return cache[r['canonical']]
    def accepted(c,name):
        p=Path(c['runDirectory'])/name;record=reader.retain(p,c['artifacts'][name]['sha256']);return load(record),record
    mref=scalar_cp['upstreamCheckpoints']['baselineFrequencyMatrix'];mc=reader.json(mref['logical'],mref['sha256']);reader.retain(Path(mc['runDirectory'])/'checks.json',mc['checksSha256'])
    ecp=reader.json(M/'S11c_d_remaining_case_matrices_checkpoint.json');require(ecp['status']=='ACCEPTED_FOUR_CASE_INTERIOR_MATRICES','accepted own Eulerian row outputs')
    reader.retain(Path(ecp['runDirectory'])/'checks.json',ecp['checksSha256'])
    ccp=reader.json(M/'S11c_d_remaining_case_coordinate_matrices_checkpoint.json');require(ccp['status']=='ACCEPTED_CASE_MATERIAL_INTERIOR_MATRICES','accepted saved source preparation')
    reader.retain(Path(ccp['runDirectory'])/'checks.json',ccp['checksSha256'])
    sources={}
    for name,names,c in (
        ('S11c_d_finite_scattering.py',('BasisMomentum','polynomial_basis','source_jets'),mc),
        ('S11c_d_frequency_matrix.py',('bind','complex_case','assemble','system_and_solve'),mc),
        ('S11c_d_remaining_case_bindings.py',('signature',),ecp),
        ('S11c_d_remaining_case_matrices.py',('main','preflight'),ecp),
        ('S11c_d_remaining_case_coordinate_matrices.py',('construct',),ccp)):
        path=M/name;h=c['sourceFiles']['_measurements/'+name];reader.retain(path,h)
        frozen=Path(c['runDirectory'])/'source/_measurements'/name;reader.retain(frozen,h)
        module=ast.parse(path.read_text());sources[name]={'file':reader.retain(path),'bodies':{n.name:ast.unparse(n) for n in module.body if getattr(n,'name',None) in names}}
    for p in (Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_row_inputs_plan.md',Path(io.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(p)
    io.save(base,'native-callers.json',sources)
    # Baseline complex full row outputs and original own real-frequency rows.
    baseline,baseline_route=accepted(mc,'complex/frequency-binding.pickle');br,br_route=accepted(mc,'complex/frequency-rows.pickle')
    system,system_route=accepted(ecp,'accepted-finite-system.pickle')
    require(len(system['nodes'])==129 and baseline['settings']==system['settings'],'full approved quadrature settings and129 field coefficients')
    atlas=[];oldcases={};prepared=[]
    for row in baseline['bound']['rows']:
        index=row['index'];array=br['rows'][index]
        require(array.shape==(129,129) and np.isfinite(array).all(),'actual complete baseline frequency row')
        atlas.append({'signature':signature(row,baseline),'binding':baseline_route,'owner':['baseline-complex',index],
                      'result':{'packet':br_route,'keys':['rows',index]},'array':array})
    for label in ecp['cases']:
        old,old_route=accepted(ecp,'accepted-bindings/'+label+'/case-binding.pickle');data=old['binding'];oldcases[label]=(data,old_route)
        rows,rows_route=accepted(ecp,'cases/'+label+'/row-matrices.pickle')
        require(data['settings']==system['settings'],'own saved original row settings')
        for row in data['bound']['rows']:
            index=row['index'];array=rows['rows'][index];require(array.shape==(129,129) and np.isfinite(array).all(),'actual complete original case row')
            atlas.append({'signature':signature(row,data),'binding':old_route,'owner':[label,index],
                          'result':{'packet':rows_route,'keys':['rows',index]},'array':array})
        if label!=BASELINE:
            name='cases/'+label+'/prepared-source-basis.pickle';saved,saved_route=accepted(ccp,name)
            original=saved['original'];require(original['size']==129 and original['source_nodes'].shape==original['source_weights'].shape,'saved complete original source basis')
            for si,jet in original['jet_data'].items():
                a=original['amplitudes'][si];require(a.shape==(len(original['source_nodes']),129) and np.isfinite(a).all(),'complete actual prepared source array')
                prepared.append({'jet':jet,'nodes':original['source_nodes'],'weights':original['source_weights'],
                    'size':original['size'],'result':{'packet':saved_route,'keys':['original','amplitudes',si]},
                    'input':{'packet':saved_route,'keys':['original','jet_data',si]},'owner':[label,si],
                    'nodeRoute':{'packet':saved_route,'keys':['original','source_nodes']},
                    'weightRoute':{'packet':saved_route,'keys':['original','source_weights']}})
            io.save(base,'saved-prepared-basis/'+label+'.json',{'packet':saved_route,'nodeCount':len(original['source_nodes']),
                'sourceIndices':list(original['jet_data']),'basisShape':[len(original['source_nodes']),129],
                'retainedFields':list(original),'sourceSettings':system['settings'],'nativeCaller':sources['S11c_d_remaining_case_coordinate_matrices.py']})
    io.save(base,'saved-row-atlas.json',[{k:v for k,v in a.items() if k not in ('signature','array')} for a in atlas])
    pending=[];case_counts={};basis_counts={}
    for name in cp['artifacts']:
        if not name.endswith('/source-coefficient-routes.json'):continue
        routes=reader.json(origin/name,cp['artifacts'][name]['sha256']);label=routes[0]['consumer']['case'];jets={};jet_routes={}
        for route in routes:
            target=route['result'];jet=load(target['packet'])
            for k in target['keys']:jet=jet[k]
            si=route['consumer']['sourceIndex'];jets[si]=jet;jet_routes[si]=target
        request=routes[0]['consumer']['fullInputRoutes'];raw=load(request['nativeSource']['packet'])
        scalars,sr_route=accepted(scalar_cp,'cases/'+label+'/scalar-bindings.pickle')
        bound=raw['bound'];rows=[]
        for oldrow in bound['rows']:
            factors=[dict(v,coefficient=scalars['actual']['factor',oldrow['index'],i]) for i,v in enumerate(oldrow['factors'])]
            rows.append(dict(oldrow,factors=factors))
        data={'bound':dict(bound,rows=rows),'jets':jets,'settings':raw['settings'],'fieldUnits':raw['fieldUnits'],'equationUnits':raw['equationUnits']}
        require(data['settings']==system['settings'],'desired full row settings')
        own_routes={'packet':reader.retain(sr/'cases'/label/'input-routes.json',scalar_cp['artifacts']['cases/'+label+'/input-routes.json']['sha256']),'keys':[]}
        row_results=[];new_rows=0;saved_rows=0
        for row in rows:
            index=row['index'];sig=signature(row,data)
            # Preserve separate dimensions/type diagnostics. No coercion and no
            # declaration of numerical reuse from partial matches.
            candidates=[v for v in atlas if same(v['signature'],sig)]
            result={'case':label,'rowIndex':index,'sourceInputs':request['nativeSource'],'scalarInputs':sr_route,
                    'factorAddresses':[['factor',index,i] for i in range(len(row['factors']))],
                    'jets':{str(v['sourceIndex']):jet_routes[v['sourceIndex']] for v in row['factors']},
                    'ownPhysicalRoutes':own_routes,'settings':data['settings'],'basis':system_route,
                    'savedCompleteMatches':[{'owner':v['owner'],'binding':v['binding'],'result':v['result']} for v in candidates],
                    'layoutDimension':len(row['limits']),'integrationVariables':[str(v[0]) for v in row['limits']]}
            if candidates:
                require(all(np.array_equal(candidates[0]['array'],v['array']) for v in candidates),'all full saved row results agree')
                result['status']='SAVED_COMPLETE_NUMERICAL_ROW';saved_rows+=1
            else:
                siblings=[i for i,v in enumerate(pending) if same(v['signature'],sig)]
                family=siblings[0] if siblings else len(pending)
                if not siblings:pending.append({'signature':sig,'owner':[label,index]});new_rows+=1
                result.update(status='UNMATCHED_FULL_ROW_INPUT',pendingFamily=family,firstOwner=pending[family]['owner'])
                # Diagnose only full component equality, never infer science reuse.
                result['partialComponentCandidates']=[{'owner':v['owner'],'equalComponents':[k for k in sig if same(v['signature'][k],sig[k])]} for v in atlas if same(v['signature']['row'],sig['row'])]
            io.save(base,'cases/'+label+'/row-'+str(index)+'.json',result);row_results.append({'rowIndex':index,'status':result['status'],'savedMatches':len(candidates),'pendingFamily':result.get('pendingFamily')})
        source_results=[]
        for si,jet in jets.items():
            key=(jet['probe'],tuple(jet['coefficients']),jet['amplitudeUnit'],jet['integralUnit'])
            matches=[v for v in prepared if same(key,(v['jet']['probe'],tuple(v['jet']['coefficients']),v['jet']['amplitudeUnit'],v['jet']['integralUnit']))]
            source_results.append({'sourceIndex':si,'coefficientRoute':jet_routes[si],
                'savedCoefficientBasisCandidates':[{k:v[k] for k in ('owner','result','input','nodeRoute','weightRoute')} for v in matches],
                'scope':'Candidate full coefficients on stored nodes/weights; actual future source-rule and full basis caller join remains required.'})
        io.save(base,'cases/'+label+'/source-basis-inputs.json',source_results)
        counts={'rows':len(rows),'savedCompleteRows':saved_rows,'unmatchedUses':len(rows)-saved_rows,'newInputFamilies':new_rows,
                'sources':len(jets),'savedBasisCoefficientCandidates':sum(bool(v['savedCoefficientBasisCandidates']) for v in source_results)}
        io.save(base,'cases/'+label+'/summary.json',{'counts':counts,'rows':row_results})
        case_counts[label]=counts
    io.save(base,'pending-row-families.json',[{'family':i,'firstOwner':v['owner'],'layoutDimension':len(v['signature']['row'][2]),'variables':[str(x[0]) for x in v['signature']['row'][2]]} for i,v in enumerate(pending)])
    reader.postcheck();io.save(base,'inputs.json',{'acceptedCoefficientCheckpoint':reader.retain(CP,SHA),'savedSystem':system_route,'consumedRoutes':reader.routes})
    checks={'status':'COMPLETED_SAVED_FREQUENCY_ROW_INPUT_INSPECTION','cases':case_counts,'savedRowInputs':len(atlas),
        'savedPreparedSourceArrays':len(prepared),'pendingLiteralRowFamilies':len(pending),'newScientificCalls':0,'allConsumedHashesUnchanged':True,
        'wallSeconds':time.monotonic()-start,'scope':'Actual full row and saved source basis inputs only. No new quadrature, source evaluator, basis, matrix, solve or numerical reuse acceptance. Use actual results to scope the next small integration pilot.'}
    io.save(base,'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
