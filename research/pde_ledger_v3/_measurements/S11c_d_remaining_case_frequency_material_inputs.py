#!/usr/bin/env python3
"""Inspect only the remaining material row operands and saved source-call routes."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time
import S11c_d_remaining_case_frequency_remainder_inputs_finish as metadata

saved=metadata.original.saved
p,io,sp=saved.p,saved.io,saved.sp
M,F,require,same=saved.M,saved.F,saved.require,saved.same
NAME='S11c_d_remaining_case_frequency_material_inputs'
SELECTION={'MATERIAL_ADVECTED__RHO4_CONSTANT':(2,3,4,5,6,7,8,9,10,11,21,22,23,24,25,60,61,62,63),
           'MATERIAL_ADVECTED__RHOBR_CONSTANT':(50,51,52,53)}


def expression_view(value):
    return {'type':type(value).__name__,'expression':str(value),'arguments':[expression_view(v) for v in value.args]}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();io.digest=saved.digest
    reader,journal,cache=saved.Reader(),metadata.MetadataJournal(base),{}
    def packet(rec):
        r=reader.retain(rec.get('logical',rec.get('path')),rec['sha256'])
        if r['canonical'] not in cache:cache[r['canonical']]=reader.packet(r['logical'],r['sha256'])
        return cache[r['canonical']]
    def meta(rec):return reader.json(rec.get('logical',rec.get('path')),rec['sha256'])
    def address(route):
        v=packet(route['packet'])
        for key in route['keys']:v=v[tuple(key) if isinstance(key,list) else key]
        return v
    latest=reader.json(M/'S11c_d_remaining_case_frequency_lab_row_sensitivity_checkpoint.json','cf12257fe054d83e35fa48db566db8a98698596d29cac017df5552152df5d580')
    require(latest['status']=='ACCEPTED_BOUNDED_LAB_ROW46_RESPONSE_SENSITIVITY','accepted LAB response sensitivity before material work')
    inventory=reader.json(p.READY/'completed-input-artifact-inventory.json',p.INVENTORY_SHA)
    reader.retain(p.READY/'complete/checks.json',p.READY_SHA)
    def prior(name):return meta(inventory[name])
    prep=reader.json(M/'S11c_d_remaining_case_frequency_remainder_preparation_checkpoint.json','3df3fb8b3bf9796eb32f22b8d86ab68f3708811bc24119e2757b93866a3a3a25')
    pc=meta(prep['manifestReferences']['artifacts']['file']);nc=meta(pc['artifacts']['native-and-new-numerical-callers.json']);caller=meta(nc['savedCaller'])
    for entry in caller['native'].values():
        r=reader.retain(entry['file']['logical'],entry['file']['sha256']);parsed=ast.parse(Path(r['canonical']).read_text())
        for name,body in entry['bodies'].items():
            require(ast.dump(next(n for n in parsed.body if getattr(n,'name',None)==name))==ast.dump(ast.parse(body).body[0]),'whole native source/basis/row caller body')
    for rec in [caller['complexCaller'],caller['nativeHeavisidePrinter']['file'],caller['nativeHeavisidePrinter']['namespaceFile']]:reader.retain(rec['logical'],rec['sha256'])
    require(caller['oldPrepareBasisRan'],'old accepted native preparation actually ran; missing internal packet is not absence of a call')
    journal.json('native-callers-and-prior-preparation.json',{'savedCaller':nc['savedCaller'],'oldPrepareBasisRan':True,'oldInternalReturnAbsenceIsNotReplayAuthorization':True,
        'previousNewNumericalMethods':pc['artifacts']['native-and-new-numerical-callers.json'],'currentWholeBodyJoins':True,'nativeMethodsCalled':0})
    def forbidden(*args,**kwargs):raise RuntimeError('saved material input inspection cannot execute science')
    journal.write=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    # Actual new source actions are extra candidates beyond the original85 arrays.
    available=[]
    for v in prep['sourceActions']:
        old=packet(v['ownInput']);ai=packet(v['actionInput'])
        rec=pc['artifacts']['sources/'+str(v['sourceIndex'])+'/action-completed.json'];receipt=meta(rec)
        require(receipt=={'input':v['actionInput'],'value':v['actionValue']},'actual new source action receipt')
        available.append({'owner':['LAB_HELD__RHOBR_CONSTANT',v['sourceIndex']],'input':v['ownInput'],'actionInput':v['actionInput'],'return':v['actionValue'],
            'receipt':rec,'jet':old['jet'],'context':old['context'],'nodes':address(ai['nodes']),'weights':address(ai['weights']),'bound':ai['bound'],'size':ai['size']})
    first=reader.json(M/'S11c_d_remaining_case_frequency_row_pilot_checkpoint.json','0c2448600fbd33b82b5e9c1a0d6b16ddbf5dbb423decff540ceec7d221098e15')
    fc=reader.json(Path(first['runDirectory'])/'checks.json',first['checksSha256']);fa=fc['artifacts'];old=packet(fa['source-action/input.pickle']);fr=meta(fa['source-action/completed.json'])
    require(fr['value']==fa['source-action/value.pickle'],'actual completed source22 recurrence')
    available.append({'owner':['LAB_HELD__RHOBR_CONSTANT',22],'input':fa['source-action/input.pickle'],'return':fa['source-action/value.pickle'],
        'receipt':fa['source-action/completed.json'],'jet':old['jet'],'context':old['context'],'nodes':old['nodes'],'weights':old['weights'],'bound':old['bound'],'size':old['size']})
    source_keys=('probe','column','coefficients','amplitudeUnit','integralUnit','originalBoundAmplitude')
    summaries=[];sources=[];source_inputs={};new_source_families=[]
    for case,rows in SELECTION.items():
        basis=prior('cases/'+case+'/source-basis-inputs.json');rule_meta=prior('saved-prepared-basis/'+case+'.json');prepared=packet(rule_meta['packet'])['original']
        nodes,weights=prepared['source_nodes'],prepared['source_weights'];node_route={'packet':rule_meta['packet'],'keys':['original','source_nodes']};weight_route={'packet':rule_meta['packet'],'keys':['original','source_weights']}
        for ri in rows:
            route=prior('cases/'+case+'/row-'+str(ri)+'.json');require(route['firstOwner']==[case,ri] and not route['savedCompleteMatches'] and route['status']=='UNMATCHED_FULL_ROW_INPUT','actual unfinished material first-owner input')
            own,scalars,system=packet(route['sourceInputs']['packet']),packet(route['scalarInputs']),packet(route['basis'])
            physical=meta(route['ownPhysicalRoutes']['packet']);common=packet(physical['context']);context=common['contextPair'][0]
            require(same(*common['contextPair']) and same(*common['basisPair']),'full actual own source/context/field basis')
            for rec in physical.values():
                if isinstance(rec,dict) and 'logical' in rec:reader.retain(rec['logical'],rec['sha256'])
            require(same(own['settings'],system['settings']) and same(own['settings'],context['settings']) and
                    same(rule_meta['sourceSettings'],json.loads(json.dumps(own['settings']))) and complex(scalars['frequency'])==1-.01j,'full case source-rule and fixed frequency settings')
            require(prepared['size']==len(system['nodes'])==129 and nodes.shape==weights.shape==(1024,),'actual source rule/trial dimensions')
            row=own['bound']['rows'][ri];require(row['index']==ri and len(row['factors'])==1,'entire actual one-factor row')
            factor=row['factors'][0];si=factor['sourceIndex'];source=own['bound']['sources'][0,si];jet_route=route['jets'][str(si)];jet=address(jet_route);coefficient=scalars['actual']['factor',ri,0]
            require(same(jet['originalBoundAmplitude'],scalars['actual']['source',si]) and same(jet['amplitudeUnit'],source['amplitudeUnit']) and
                    same(jet['integralUnit'],source['integralUnit']) and jet['probe'].args==(context['zp'],),'actual full source coefficient/field/unit input')
            source_key=(case,si)
            if source_key not in source_inputs:
                candidate=next(v for v in basis if v['sourceIndex']==si);require(candidate['coefficientRoute']==jet_route,'saved source-basis candidate input identity')
                comparisons=[]
                for v in available:
                    flags={'wholeJet':same({k:jet[k] for k in source_keys},{k:v['jet'][k] for k in source_keys}),'context':same(context,v['context']),
                        'nodes':same(nodes,v['nodes']),'weights':same(weights,v['weights']),'bound':same(own['settings']['sourceBound'],v['bound']),'size':prepared['size']==v['size']}
                    comparisons.append({'owner':v['owner'],'input':v['input'],'return':v['return'],'receipt':v['receipt'],'components':flags,'wholeSourceCallMatch':all(flags.values())})
                descriptor={'jet':{k:jet[k] for k in source_keys},'context':context,'nodes':nodes,'weights':weights,'bound':own['settings']['sourceBound'],'size':prepared['size']}
                family_matches=[j for j,v in enumerate(new_source_families) if same(v['descriptor'],descriptor)]
                fi=family_matches[0] if family_matches else len(new_source_families)
                if not family_matches:new_source_families.append({'descriptor':descriptor,'owner':[case,si]})
                info={'case':case,'sourceIndex':si,'sourceJet':jet_route,'nativeSource':{'packet':route['sourceInputs']['packet'],'keys':['bound','sources',[0,int(si)]],'keyEncoding':'last key is exact tuple'},
                    'coefficients':[expression_view(v) for v in jet['coefficients']],'order':len(jet['coefficients'])-1,'probe':str(jet['probe']),'column':jet['column'],
                    'sourceFrequency':expression_view(source['frequency']),'amplitudeUnit':jet['amplitudeUnit'],'integralUnit':jet['integralUnit'],
                    'physicalRoutes':physical,'settings':own['settings'],'nodeRoute':node_route,'weightRoute':weight_route,'size':prepared['size'],
                    'originalPreparedArrayCandidates':candidate,'newSavedActionComparisons':comparisons,'sourceInputFamily':fi,'sourceInputFamilyOwner':new_source_families[fi]['owner'],
                    'scope':'Saved input comparison only. A partial match is not a source action or row return. Any missing action uses an explicit new numerical method, never recreation of old native internal derivatives/lambdas.'}
                path='cases/'+case+'/sources/'+str(si)+'.json';journal.json(path,info);source_inputs[source_key]=journal.artifacts[path]
                sources.append({'case':case,'sourceIndex':si,'order':len(jet['coefficients'])-1,'input':source_inputs[source_key],
                    'originalArrayCandidates':len(candidate['savedCoefficientBasisCandidates']),'newFullActionMatches':sum(v['wholeSourceCallMatch'] for v in comparisons),'inputFamily':fi})
            profiles=list(coefficient.atoms(sp.Integral))
            info={'case':case,'rowIndex':ri,'sourceIndex':si,'route':route,'sourceView':source_inputs[source_key],'coefficientTree':expression_view(coefficient),
                'coefficientUnit':factor['unit'],'coefficientFreeSymbols':sorted(map(str,coefficient.free_symbols)),
                'profiles':[{'tree':expression_view(v),'unit':own['bound']['profileUnits'][v]} for v in profiles],
                'sourceFrequency':expression_view(source['frequency']),'rowLimits':[str(v) for v in row['limits']],'layoutDimension':route['layoutDimension'],
                'fieldUnits':own['fieldUnits'],'equationUnits':own['equationUnits'],'settings':own['settings'],'physicalRoutes':physical,
                'abelSource':{'packet':route['sourceInputs']['packet'],'keys':['bound','abel']},'orderedLimitsPreserved':True,'newScientificCalls':0}
            path='cases/'+case+'/rows/'+str(ri)+'.json';journal.json(path,info)
            summaries.append({'case':case,'rowIndex':ri,'sourceIndex':si,'layoutDimension':route['layoutDimension'],'profiles':len(profiles),'sourceOrder':len(jet['coefficients'])-1,'input':journal.artifacts[path]})
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md'),Path(saved.__file__),Path(metadata.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md')):
        dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out:out.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    journal.json('material-source-summary.json',sources);journal.json('material-row-summary.json',summaries)
    journal.json('inputs.json',{'consumedRoutes':dict(reader.routes),'latestLABSensitivity':reader.retain(M/'S11c_d_remaining_case_frequency_lab_row_sensitivity_checkpoint.json'),
        'originalRowInputInventory':reader.retain(p.READY/'completed-input-artifact-inventory.json'),'sourceActionCheckpoint':reader.retain(M/'S11c_d_remaining_case_frequency_remainder_preparation_checkpoint.json')})
    reader.postcheck()
    checks={'status':'COMPLETED_SAVED_MATERIAL_ROW_INPUT_INSPECTION','rows':summaries,'sources':sources,'sourceInputFamilies':len(new_source_families),
        'newScientificCalls':0,'newPhysicalPickleWrites':0,'consumedPaths':len(reader.routes),'allConsumedHashesUnchanged':True,
        'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-started,'scope':'Exact remaining material inputs and actual existing source-call matches only; no numerical reuse acceptance from partial input families.'}
    journal.json('checks.json',checks);signal.alarm(0);print((base/'checks.json').read_text(),end='')


if __name__=='__main__':main()
