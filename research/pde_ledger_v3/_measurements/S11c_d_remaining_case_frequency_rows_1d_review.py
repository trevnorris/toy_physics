#!/usr/bin/env python3
"""Bounded streaming review of saved two-row results, without new science."""
import argparse
import ast
import gc
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import time

import S11c_d_remaining_case_frequency_rows_1d as producer
p=producer.p
io,np,sp=p.io,p.np,p.sp
M,F,require,same=p.M,p.F,p.require,p.same
ROOT=F/'rows-1d'
CHECKS_SHA='d0c508a88fa9196e196072be1cb8b354ab82eedc1d1e5e96eeafecee2e86a194'
HELPER_SHA='45085deed3cb1bbec5bd0b9e5a573f2ce855b3a7cd4f8c54f9673f6e20886c90'


def release(stream):
    # This hints that clean pages just read need not remain cached. It changes
    # neither file bytes nor stored scientific values and raises no resource cap.
    os.posix_fadvise(stream.fileno(),0,0,os.POSIX_FADV_DONTNEED)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
        release(stream)
    return h.hexdigest()


class StreamingReader(io.Reader):
    def packet(self,path,expected):
        route=self.retain(path,expected)
        with Path(route['canonical']).open('rb') as stream:
            value=pickle.load(stream);release(stream)
        return value


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic()
    io.digest=digest;reader=StreamingReader();large={}
    def packet(record):return reader.packet(record.get('path',record.get('logical')),record['sha256'])
    def address(route):
        r=route['packet'];path=r.get('path',r.get('logical'))
        if path not in large:large[path]=packet(r)
        value=large[path]
        for key in route['keys']:value=value[key]
        return value
    inspection=reader.json(ROOT/'completion-inspection.json');checks=reader.json(ROOT/'complete/checks.json',CHECKS_SHA)
    require(inspection['actualExits']==[0,0,0] and inspection['guardInterruption'] is None and inspection['checksStdoutIdentity'] and inspection['emptyStrictStderr'] and inspection['zeroOOMSwap'],'final producer outcome and contained memory qualification')
    require(inspection['peakBytes']==2147483648 and inspection['maximumMemoryMaxEvents']==176 and not inspection['zeroCapEvents'],'preserve actual nonzero memory-cap contacts')
    for r in inspection['evidence'].values():reader.retain(r['path'],r['sha256'])
    outcomes=[reader.json(ROOT/q) for q in ['resource-guard/outcome.json','resource-guard/child-outcome.json','frequency_rows_1d.invocation.json']]
    require(all(q['exitCode']==0 for q in outcomes),'actual final producer exits')
    reader.retain(Path(producer.__file__),HELPER_SHA)
    for path in (Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_rows_1d_review_plan.md',producer.PLAN,Path(p.__file__),Path(producer.recovery.__file__),Path(io.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):reader.retain(path)
    artifacts=checks['artifacts']
    def saved(name):
        r=artifacts[name]
        return reader.json(r['path'],r['sha256']) if name.endswith('.json') else packet(r)
    for r in artifacts.values():require(reader.retain(r['path'],r['sha256'])['bytes']==r['bytes'],'all actual produced bytes')
    inputs=saved('inputs.json')
    for logical,r in inputs['consumedRoutes'].items():require(reader.retain(logical,r['sha256'])==r,'complete consumed logical/canonical routes')
    scope=saved('batch-scope.json');callers=saved('native-and-numerical-callers.json')
    for entry in callers['native'].values():
        r=reader.retain(entry['file']['logical'],entry['file']['sha256']);tree=ast.parse(Path(r['canonical']).read_text())
        for name,body in entry['bodies'].items():require(ast.dump(next(n for n in tree.body if getattr(n,'name',None)==name))==ast.dump(ast.parse(body).body[0]),'whole actual native caller body')
    module,join=producer.recovery.adapted(Path(p.__file__).read_text());compile(module,'<unexecuted-row-adapter>','exec')
    require(join==callers['evaluatorJoin'] and ast.unparse(ast.Module(body=[module.body[0]],type_ignores=[]))==callers['evaluator'],'saved unchanged numerical evaluator source')
    def forbidden(*args,**kwargs):raise RuntimeError('scientific operation forbidden in saved row review')
    producer.main=p.evaluate=p.source_matrix=p.independent_source=p.integrate.quad_vec=forbidden
    p.storage.Journal.write=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    cp=reader.json(producer.CP,producer.CP_SHA);pc=reader.json(Path(cp['runDirectory'])/'checks.json',cp['checksSha256']);pilot=packet(cp['rowInput'])
    catalogue=reader.json(F/'row-pilot-review/complete/validated-point-catalogue.json',producer.CATALOGUE_SHA);old_points={v['momentumHex']:v for v in catalogue}
    inventory=reader.json(p.READY/'completed-input-artifact-inventory.json',p.INVENTORY_SHA)
    def metadata(name):return reader.json(p.READY/'complete'/name,inventory[name]['sha256'])
    source_meta=metadata('saved-prepared-basis/'+producer.CASE+'.json');source_saved=address({'packet':source_meta['packet'],'keys':['original']})
    nodes,weights=source_saved['source_nodes'],source_saved['source_weights']
    require(source_saved['size']==129 and nodes.shape==weights.shape==(1024,),'actual full saved source-rule dimensions')
    basis_meta=metadata('cases/'+producer.CASE+'/source-basis-inputs.json')
    selected=metadata('cases/'+producer.CASE+'/row-29.json');raw=packet(selected['sourceInputs']['packet']);scalars=packet(selected['scalarInputs']);system=packet(selected['basis'])
    physical=reader.json(selected['ownPhysicalRoutes']['packet']['logical'],selected['ownPhysicalRoutes']['packet']['sha256']);common=packet(physical['context']);context=common['contextPair'][0]
    known_scalars={};known_fourier={};summaries=[];point_catalogue=[]
    for expected in checks['rows']:
        ri=expected['rowIndex'];si=expected['sourceIndex'];prefix='rows/'+str(ri);row=saved(prefix+'/input.pickle');routes=saved(prefix+'/source-basis-route.json');pairs=saved(prefix+'/saved-scalar-input-pairs.pickle');summary=saved(prefix+'/summary.json')
        require(summary==expected and (ri,si) in ((29,21),(30,24)),'actual two saved row owners')
        selected=metadata('cases/'+producer.CASE+'/row-'+str(ri)+'.json');jet=address(selected['jets'][str(si)])
        require(same(row['row'],raw['bound']['rows'][ri]) and same(row['coefficient'],scalars['actual']['factor',ri,0]) and same(row['source'],raw['bound']['sources'][0,si]),'actual whole row/scalar/source operands')
        require(same(row['jet'],jet) and same(jet['originalBoundAmplitude'],scalars['actual']['source',si]) and same(row['context'],context) and same(*common['contextPair']) and same(*common['basisPair']),'complete source jet and own context')
        for k in ('fieldUnits','equationUnits','settings'):require(same(row[k],raw[k]),'actual field-equation/settings input')
        require(same(row['settings'],context['settings']) and same(row['settings'],system['settings']) and same(json.loads(json.dumps(row['settings'])),source_meta['sourceSettings']),'whole saved source-rule and trial caller settings')
        for k in ('abel','pairs','profileUnits'):require(same(row[k],raw['bound'][k]),'actual measure/profile/ordered-limit input')
        require(same(row['physicalRoutes'],physical) and same(row['frequency'],scalars['frequency']) and same(row['positions'],system['nodes']),'complete frequency/position/physical routes')
        candidates=next(v for v in basis_meta if v['sourceIndex']==si)['savedCoefficientBasisCandidates'];own=[v for v in candidates if v['owner'][0]==producer.CASE]
        require(routes['matches']==own and routes['chosen']==row['sourceAction']==own[0]['result'] and routes['newSourceActionCalls']==0,'actual saved source action route')
        key=(jet['probe'],tuple(jet['coefficients']),jet['amplitudeUnit'],jet['integralUnit']);matrices=[]
        for candidate in own:
            j=address(candidate['input']);require(same(key,(j['probe'],tuple(j['coefficients']),j['amplitudeUnit'],j['integralUnit'])),'full saved source-basis typed inputs')
            require(same(address(candidate['nodeRoute']),nodes) and same(address(candidate['weightRoute']),weights),'exact saved source-rule nodes/measure')
            matrices.append(address(candidate['result']))
        require(all(v.shape==(1024,129) and v.dtype==np.dtype(complex) and np.isfinite(v).all() and same(v,matrices[0]) for v in matrices),'complete saved source arrays')
        variable=row['row']['limits'][0][0]
        common_match=same(context,pilot['context']) and same(row['positions'],pilot['positions']) and same(row['fieldUnits'],pilot['fieldUnits']) and same(row['equationUnits'],pilot['equationUnits']) and same(row['settings'],pilot['settings']) and same(variable,pilot['row']['limits'][0][0])
        flags={'coefficient':common_match and same(row['coefficient'],pilot['coefficient']) and same(row['row']['factors'][0]['unit'],pilot['row']['factors'][0]['unit']),'frequency':common_match and same(row['source']['frequency'],pilot['source']['frequency'])}
        require(same(pairs['requested'],row) and same(pairs['accepted'],pilot) and pairs['acceptedRowInput']==cp['rowInput'] and pairs['matches']==flags,'actual full saved scalar candidate arguments')
        io.save(base,prefix+'/validated-input-and-source.json',{'rowIndex':ri,'sourceIndex':si,'sourceAction':row['sourceAction'],'fullNativeSourceUnitsContextJoined':True,'savedScalarMatches':flags})
        count=summary['counts']['newPoints'];observed={'newPoints':0,'newCoefficients':0,'newFrequencies':0,'savedScalars':0,'newFourier':0,'savedFourier':0};seen=set()
        expressions={'coefficient':row['coefficient'],'frequency':row['source']['frequency']};units={'coefficient':row['row']['factors'][0]['unit'],'frequency':('source-frequency',variable)}
        for index in range(count):
            folder=prefix+'/points/'+str(index);inp=saved(folder+'/input.pickle');receipt=saved(folder+'/completed.json');operands=saved(folder+'/operand-routes.json')
            require(receipt['input']==artifacts[folder+'/input.pickle'] and receipt['value']==artifacts[folder+'/value.pickle'],'actual full point immediate receipt')
            require(type(inp['momentum']) is float and -4<=inp['momentum']<=4 and same(inp,{'momentum':inp['momentum'],'variable':variable,'positions':selected['basis'],'rowInput':summary['rowInput'],'sourceAction':row['sourceAction']}),'complete literal point caller')
            h=inp['momentum'].hex();require(h not in seen,'unique completed point owner');seen.add(h);values={}
            for kind in ('coefficient','frequency'):
                ref=operands[kind];environment={variable:inp['momentum']}
                if kind=='coefficient':environment.update({context['z']:row['positions'],context['regulator']:row['settings']['regulator']})
                request={'expression':expressions[kind],'environment':environment,'unit':units[kind]}
                if folder+'/'+kind+'-input.pickle' in artifacts:
                    ci=saved(folder+'/'+kind+'-input.pickle');cr=saved(folder+'/'+kind+'-completed.json')
                    require(same(ci,dict(request,rowInput=summary['rowInput'])) and cr=={'input':artifacts[folder+'/'+kind+'-input.pickle'],'value':ref},'actual new scalar before/after receipt')
                    require(ref==artifacts[folder+'/'+kind+'-value.pickle'] and ref['path'] not in known_scalars,'first actual new scalar owner')
                    known_scalars[ref['path']]=(request,ref);observed['newCoefficients' if kind=='coefficient' else 'newFrequencies']+=1
                elif ref['path'] in known_scalars:
                    old_request,old_ref=known_scalars[ref['path']];require(same(request,old_request) and ref==old_ref,'exact preceding scalar input and actual saved return');observed['savedScalars']+=1
                else:
                    require(flags[kind] and h in old_points,'actual accepted scalar input match')
                    owner=old_points[h];old_input=packet(owner['input']);require(old_input['momentum'].hex()==h and same(old_input['variable'],variable) and old_input['rowInput']==cp['rowInput'],'actual full original scalar point caller')
                    old_ref=pc['artifacts']['integrand/'+str(owner['index'])+'/'+kind+'-value.pickle'];require(ref==old_ref,'actual accepted scalar value address')
                    known_scalars[ref['path']]=(request,ref);observed['savedScalars']+=1
                values[kind]=packet(ref)
            coefficient,frequency=values['coefficient'],values['frequency']
            require(coefficient.shape==(129,) and coefficient.dtype==np.dtype(complex) and np.isfinite(coefficient).all() and type(frequency) is complex and abs(frequency.imag)<1e-14,'complete saved scalar return types and sizes')
            ref=operands['sourceFourier']
            if folder+'/fourier-input.pickle' in artifacts:
                fi=saved(folder+'/fourier-input.pickle');fr=saved(folder+'/fourier-completed.json')
                require(same(fi,{'frequency':frequency,'sourceAction':row['sourceAction'],'sourceNodes':routes['nodeRoute'],'unit':jet['integralUnit'],'method':'literal minus-phase vector times accepted weighted source action'}),'actual full Fourier input')
                require(fr=={'input':artifacts[folder+'/fourier-input.pickle'],'value':ref} and ref==artifacts[folder+'/fourier-value.pickle'],'actual Fourier immediate completion')
                known_fourier[ref['path']]=(fi,ref);observed['newFourier']+=1
            else:
                require(ref['path'] in known_fourier and known_fourier[ref['path']][1]==ref,'preceding actual full Fourier owner')
                require(same(known_fourier[ref['path']][0],{'frequency':frequency,'sourceAction':row['sourceAction'],'sourceNodes':routes['nodeRoute'],'unit':jet['integralUnit'],'method':'literal minus-phase vector times accepted weighted source action'}),'whole preceding Fourier inputs')
                observed['savedFourier']+=1
            source_value=packet(ref);value=saved(folder+'/value.pickle')
            require(source_value.shape==(129,) and value.shape==(129,129) and source_value.dtype==value.dtype==np.dtype(complex) and np.isfinite(source_value).all() and np.isfinite(value).all(),'saved complete Fourier and full row values; no product reconstruction')
            observed['newPoints']+=1;point_catalogue.append({'rowIndex':ri,'pointIndex':index,'momentumHex':h,'input':receipt['input'],'operands':operands,'value':receipt['value']})
            if (index+1)%512==0:io.save(base,prefix+'/validated-points/'+str(index+1)+'.json',{'firstIndex':index-511,'lastIndex':index,'fullInputsReturnsAndOwnerReceiptsJoined':True})
        require(all(summary['counts'][k]==v for k,v in observed.items()),'counts from actual operation presence and reference uses, not rebuilt live cache membership')
        require(sum(v['evaluations'] for v in summary['rules'])-count==summary['counts']['savedPoints'],'actual native evaluation census and full point return uses')
        rules=[]
        for result in summary['rules']:
            rule=result['rule'];folder=prefix+'/quadrature/'+rule;inp=saved(folder+'/input.pickle');value=saved(folder+'/value.pickle');receipt=saved(folder+'/completed.json');cost=saved(folder+'/cost-decision.json')
            require(receipt['input']==artifacts[folder+'/input.pickle'] and receipt['value']==result['value']==artifacts[folder+'/value.pickle'],'actual quadrature receipt')
            require(inp=={'rowInput':summary['rowInput'],'sourceAction':row['sourceAction'],'interval':(-4.,4.),'rule':rule,'epsabs':1e-8 if rule=='gk21' else 5e-9,'epsrel':1e-6 if rule=='gk21' else 5e-7,'workers':1,'norm':'max','limit':256,'cacheSize':32*1024**2},'full actual quadrature input')
            require(cost['requiredReserveSeconds']<cost['remainingSeconds'] and value['success'] and value['status']==0 and np.isfinite(value['row']).all() and value['row'].shape==(129,129) and value['row'].dtype==np.dtype(complex),'actual budget and complete successful row')
            intervals=value['intervals'];errors=value['intervalErrors'];cache=value['intervalValues'];finite=np.isfinite(cache).all(axis=(1,2));evicted=np.isnan(cache).all(axis=(1,2))
            require(intervals.shape==(len(errors),2) and cache.shape==(len(errors),129,129) and np.isfinite(intervals).all() and np.isfinite(errors).all() and np.all(finite|evicted),'complete or explicit evicted optional interval diagnostics')
            require(value['evaluations']==result['evaluations'] and value['wallSeconds']==result['wallSeconds'] and value['errorEstimate']==result['errorEstimate'],'saved native result statistics')
            rules.append({'rule':rule,'intervalCount':len(intervals),'evictedDiagnosticIntervals':int(np.count_nonzero(evicted)),**result})
            io.save(base,folder+'/validated-result.json',rules[-1])
            final_row=value['row'];del value,cache;gc.collect()
        require(summary['rowMaxNorm']==float(np.max(abs(final_row))) and summary['rowReturn']=={'packet':summary['rules'][-1]['value'],'keys':['row']},'actual final consumer row and norm')
        if ri==30:
            route=saved(prefix+'/comparison-input.json');value=saved(prefix+'/comparison-value.pickle')
            require(route=={'first':summary['rules'][0]['value'],'second':summary['rules'][1]['value']} and value['difference'].shape==(129,129) and np.isfinite(value['difference']).all(),'complete saved rule comparison')
            require(all(value[k]==v for k,v in summary['comparison'].items()) and value['scaledDifference']<1e-6,'saved comparison statistics without arithmetic replay')
        else:require(summary['comparison'] is None,'no second-rule comparison invented for row29')
        validated={'rowIndex':ri,'sourceIndex':si,'rules':rules,'counts':summary['counts'],'comparison':summary['comparison'],'rowMaxNorm':summary['rowMaxNorm'],'rowReturn':summary['rowReturn']};summaries.append(validated);io.save(base,prefix+'/validated-summary.json',validated)
    # Small next-input view: read one untouched 2D row already present in the same
    # packets, without generating any profile transform, rule or evaluator.
    next_row=metadata('cases/'+producer.CASE+'/row-46.json');native_row=raw['bound']['rows'][46]
    next_view={'route':next_row,'layoutDimension':len(native_row['limits']),'coefficientOperands':[],'newScientificCalls':0,'scope':'Saved input view for choosing the next bounded pilot, not a row or profile value.'}
    for index,factor in enumerate(native_row['factors']):
        expression=scalars['actual']['factor',46,index];integrals=expression.atoms(sp.Integral)
        next_view['coefficientOperands'].append({'factorIndex':index,'sourceIndex':factor['sourceIndex'],'sourceFrequency':str(raw['bound']['sources'][0,factor['sourceIndex']]['frequency']),'expression':str(expression),'freeSymbols':[str(v) for v in sorted(expression.free_symbols,key=str)],'integrals':[{'expression':str(v),'unit':repr(raw['bound']['profileUnits'][v]) if v in raw['bound']['profileUnits'] else None} for v in sorted(integrals,key=str)]})
    io.save(base,'next-pending-2d-input.json',next_view)
    io.save(base,'validated-point-catalogue.json',point_catalogue)
    reader.postcheck();io.save(base,'validated-paths.json',reader.routes)
    result={'status':'PASSED_BOUNDED_SAVED_TWO_ROW_REVIEW','producerChecksSha256':CHECKS_SHA,'case':producer.CASE,'rows':summaries,'fullPointReceipts':len(point_catalogue),'producerMemoryQualification':{k:inspection[k] for k in ('peakBytes','maximumMemoryMaxEvents','zeroCapEvents','zeroOOMSwap')},'newScientificCalls':0,'allConsumedHashesUnchanged':True,'consumedLogicalPaths':len(reader.routes),'wallSeconds':time.monotonic()-started,'scope':'Saved full arguments/results/units/source/receipt review of two129x129 rows at1-0.01i. Producer cap contacts preserved. No numerical evaluator, source, Fourier, products, quadrature, controls or end maps repeated; no scattering error bound or pole/domain acceptance.'}
    io.save(base,'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
