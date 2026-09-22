#!/usr/bin/env python3
"""Bounded saved review of the material1D batch and its actual reuse receipts."""
import argparse
import ast
from pathlib import Path
import resource
import signal
import time
import S11c_d_remaining_case_frequency_material_rows_1d as batch

saved,p,io,sp,np=batch.saved,batch.p,batch.io,batch.sp,batch.np
M,F,require,same=batch.M,batch.F,batch.require,batch.same
NAME='S11c_d_remaining_case_frequency_material_rows_1d_review'
ROOT=F/'material-rows-1d'
HELPER_SHA='c252fd849c870c8518f5e6868b5f405815e2f2aa16ee93ca7baac7e410e86dde'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);base=ap.parse_args().run_directory.resolve();base.relative_to(p.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();io.digest=saved.digest
    reader,journal,cache=saved.Reader(),batch.pilot.prior.inputs.MetadataJournal(base),{}
    def packet(rec):
        r=reader.retain(rec.get('logical',rec.get('path')),rec['sha256'])
        if r['canonical'] not in cache:cache[r['canonical']]=reader.packet(r['logical'],r['sha256'])
        return cache[r['canonical']]
    def meta(rec):return reader.json(rec.get('logical',rec.get('path')),rec['sha256'])
    inspection=reader.json(base.parent/'prelaunch-production-inspection.json');require(inspection['outcome']['actualExits']==[0,0,0] and inspection['outcome']['emptyStrictStderr'] and inspection['outcome']['checksStdoutIdentity'],'actual clean final producer before saved review')
    c=reader.json(ROOT/'complete/checks.json',inspection['outcome']['checksSha256']);require(c['status']=='COMPLETED_BOUNDED_MATERIAL_1D_ROWS','completed batch')
    reader.retain(batch.__file__,HELPER_SHA);reader.retain(ROOT/'complete/source'/Path(batch.__file__).name,HELPER_SHA);reader.retain(batch.pilot.__file__,batch.PILOT_SHA)
    # Compile/reverse the exact source loader only. Do not call the returned load.
    compile(batch.loader_adapter(Path(batch.pilot.__file__).read_text()),'<saved loader source review only>','exec')
    cp=reader.json(M/'S11c_d_remaining_case_frequency_material_row_pilot_checkpoint.json',batch.CP_SHA);pc=meta(cp['artifacts']['manifest']);pa=pc['artifacts'];old=packet(cp['rowInput']);old_source=packet(old['sourceAction']);merge=meta(pa['grids/1024/exact-coarse-reuse-completed.json']);route=meta(pa['saved-source-and-rule-routes.json']);source_nodes=packet(route['sourceNodes']['packet'])
    for key in route['sourceNodes']['keys']:source_nodes=source_nodes[key]
    sources_cp=meta(cp['sourceCheckpoint']);rules={n:packet(route['rules'][str(n)]) for n in (512,1024)};momentum=rules[1024]['nodes']
    require(same(momentum[::2],rules[512]['nodes']),'actual fine/coarse address identity')
    def physical(inp):return {k:inp[k] for k in ('context','settings','fieldUnits','equationUnits','abel','pairs','profileUnits')}
    def forbidden(*args,**kwargs):raise RuntimeError('saved material row review cannot execute numerical or symbolic science')
    journal.write=forbidden;batch.main=batch.pilot.main=batch.pilot.prior.main=saved.main=p.main=forbidden
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden;io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    completed={};counts={'newCoefficientBatches':0,'savedCoefficientGridUses':0,'newFourierBatches':0,'savedFourierGridUses':0,'newWeightedCoefficientCalls':0,'savedWeightedCoefficientUses':0,'newRowContractions':0,'newSourceProfileRuleEndCalls':0};results=[]
    def validate_cache(a,folder,kind,key_check):
        ar=a[folder+'/input.pickle'];arg=packet(ar);receipt=meta(a[folder+'/completed.json']);value=packet(receipt['value']);require(receipt['input']==ar,'actual cache input receipt');key=arg['requested'];key_check(key)
        matches=arg['actualMatches'];require(value.dtype==np.dtype(complex) and np.isfinite(value).all(),'full finite saved cache value')
        if receipt['disposition']=='NEW':
            require(not matches and receipt['value']==a[folder+'/value.pickle'],'actual new branch and immediate return');counts[{'CoefficientGrid':'newCoefficientBatches','FourierGrid':'newFourierBatches','WeightedCoefficient':'newWeightedCoefficientCalls'}[kind]]+=1
        else:
            require(matches and receipt['disposition']=='REUSED_'+matches[0]['origin'] and receipt['value']==matches[0]['value'],'actual captured reuse branch')
            counts['saved'+kind+'Uses']+=1
            for match in matches:
                prior_receipt=meta(match['receipt']);owner_value=packet(match['value']);require(same(value,owner_value),'full matching saved return')
                if match['origin']=='COMPLETED_NEW':
                    require(match['receipt']['path'] in completed,'reuse refers to preceding completed actual operation');prior_arg=packet(match['input']);require(same(prior_arg['requested'],key) and prior_receipt['input']==match['input'] and prior_receipt['value']==match['value'],'whole typed prior operation key/value/receipt')
                elif match['origin']=='ACCEPTED_MERGED_GRID':
                    require(prior_receipt==merge and match['input']==merge['input'],'actual accepted full fine merge owner')
                    if kind=='CoefficientGrid':
                        require(match['value']==merge['coefficient'] and same(key['expression'],old['coefficient']) and same(key['unit'],old['row']['factors'][0]['unit']) and same(key['physical'],physical(old)),'accepted whole coefficient physical input')
                        require(same(key['environment'][old['row']['limits'][0][0]],momentum[:,None]) and same(key['environment'][old['context']['z']],old['positions'][None,:]),'accepted coefficient complete argument coordinates')
                    else:
                        require(kind=='FourierGrid' and match['value']==merge['fourier'] and same(key['source'],old['source']) and same(key['jet'],old['jet']) and same(key['physical'],physical(old)) and same(key['frequencies'],momentum) and same(key['sourceNodes'],source_nodes),'accepted full Fourier caller input')
                        identity=key['weightedSourceActionIdentity'];require(identity['canonical']==str(Path(old['sourceAction']['path']).resolve()) and identity['sha256']==old['sourceAction']['sha256'] and identity['bytes']==old['sourceAction']['bytes'],'accepted complete source action byte route')
                elif match['origin']=='ACCEPTED':
                    require(kind=='WeightedCoefficient' and prior_receipt=={'input':match['input'],'value':match['value']},'actual accepted weighted operation');oi=packet(match['input']);ov=packet(oi['coefficient']);ow=packet(oi['rule'])['weights'];require(same(key['coefficient'],ov) and same(key['weights'],ow) and same(key['unit'],old['row']['factors'][0]['unit']) and same(key['physical'],physical(old)),'full accepted weighted arguments and unit/context')
                else:raise AssertionError('unknown actual cache owner origin')
        completed[a[folder+'/completed.json']['path']]=receipt
        return value,receipt['value'],key
    require([v['rowIndex'] for v in c['rows']]==list(batch.ROWS),'actual requested14row order')
    for row_result in c['rows']:
        ri=row_result['rowIndex'];rc=reader.json(ROOT/f'complete/row-{ri}/checks.json');a=rc['artifacts'];inp=packet(row_result['rowInput']);require(rc['result']==row_result and row_result['savedReadback'] and inp['row']['index']==ri and inp['row']['factors'][0]['sourceIndex']==row_result['sourceIndex'],'full own result and source address')
        view=meta(inp['rowView']);raw=packet(view['route']['sourceInputs']['packet']);scalars=packet(view['route']['scalarInputs']);si=row_result['sourceIndex']
        require(same(inp['row'],raw['bound']['rows'][ri]) and same(inp['coefficient'],scalars['actual']['factor',ri,0]) and same(inp['source'],raw['bound']['sources'][0,si]) and same(inp['settings'],raw['settings']) and same(inp['fieldUnits'],raw['fieldUnits']) and same(inp['equationUnits'],raw['equationUnits']) and same(inp['abel'],raw['bound']['abel']) and same(inp['pairs'],raw['bound']['pairs']) and same(inp['profileUnits'],raw['bound']['profileUnits']),'whole original source/row/unit/settings/Abel inputs')
        accepted=next(v for v in sources_cp['sourceActions'] if v['case']==row_result['case'] and v['sourceIndex']==si);accepted_input=packet(accepted['ownInput'])
        require(inp['sourceAction']==accepted['actionValue'] and same(inp['jet'],accepted_input['jet']) and same(inp['source'],accepted_input['source']) and same(inp['context'],accepted_input['context']),'actual accepted own source action and complete typed input')
        source=packet(inp['sourceAction']);require(source.shape==(1024,129) and source.dtype==np.dtype(complex) and np.isfinite(source).all(),'full actual source action saved value')
        def ck(key):
            require(same(key['expression'],inp['coefficient']) and same(key['unit'],inp['row']['factors'][0]['unit']) and same(key['physical'],physical(inp)),'own coefficient expression/unit/context')
            expected={inp['row']['limits'][0][0]:momentum[:,None],inp['context']['z']:inp['positions'][None,:],inp['context']['regulator']:inp['settings']['regulator']};require(same(key['environment'],expected),'actual full batched coefficient environment')
        cv,cr,ckey=validate_cache(a,'coefficient-grid','CoefficientGrid',ck)
        def fk(key):
            require(same(key['frequencies'],momentum) and same(key['sourceNodes'],source_nodes) and same(key['source'],inp['source']) and same(key['frequencyExpression'],inp['source']['frequency']) and same(key['jet'],inp['jet']) and same(key['unit'],inp['jet']['integralUnit']) and same(key['physical'],physical(inp)),'whole actual own Fourier source/unit/argument tuple')
            ident=key['weightedSourceActionIdentity'];require(ident=={'canonical':str(Path(inp['sourceAction']['path']).resolve()),'sha256':inp['sourceAction']['sha256'],'bytes':inp['sourceAction']['bytes']},'actual canonical source array reference')
        fv,fr,fkey=validate_cache(a,'fourier-grid','FourierGrid',fk);require(cv.shape==fv.shape==(1025,129),'full saved fine grid dimensions')
        if 'fourier-grid/contraction-input.pickle' in a:
            fi=packet(a['fourier-grid/contraction-input.pickle']);phase=packet(fi['phase']);receipt=meta(a['fourier-grid/contraction-completed.json']);require(phase.shape==(1025,1024) and np.isfinite(phase).all() and fi['sourceAction']==inp['sourceAction'] and receipt=={'input':a['fourier-grid/contraction-input.pickle'],'value':fr},'actual phase/intermediate/full Fourier receipt')
        for grid in row_result['grids']:
            panels=grid['panels'];folder='grids/'+str(panels);selection=packet(a[folder+'/operand-routes.pickle']);indices=selection['indices'];require(selection['coefficient']==cr and selection['sourceFourier']==fr and selection['ownRowInput']==row_result['rowInput'] and selection['rule']==route['rules'][str(panels)] and list(indices)==list(range(0,1025,2) if panels==512 else range(1025)),'saved fine/coarse operand addresses')
            def wk(key):require(same(key['coefficient'],cv[indices]) and same(key['weights'],rules[panels]['weights']) and same(key['unit'],inp['row']['factors'][0]['unit']) and same(key['physical'],physical(inp)),'full saved weighted operation arguments and units')
            wv,wr,wkey=validate_cache(a,folder+'/weighted-coefficient','WeightedCoefficient',wk)
            ar=packet(grid['input']);value=packet(grid['value']);receipt=meta(a[folder+'/row-completed.json']);require(ar['weightedCoefficient']==wr and ar['sourceFourier']==fr and same(ar['sourceIndices'],indices) and ar['ownRowInput']==row_result['rowInput'] and ar['operation']=='weightedCoefficient.T @ sourceFourier' and ar['transposeConjugates'] is False and receipt=={'input':grid['input'],'value':grid['value']},'actual full native row-contraction operand/receipt joins')
            require(value.shape==(129,129) and value.dtype==np.dtype(complex) and np.isfinite(value).all(),'finite actual full row return');counts['newRowContractions']+=1
        comparison=row_result['comparison']
        if comparison:
            require(ri in batch.COMPARE and comparison==meta(a['comparison.json']),'selected saved comparison');ci=packet(comparison['input']);v=packet(comparison['value']);require(ci['coarse']==row_result['grids'][0]['value'] and ci['fine']==row_result['grids'][1]['value'] and all(same(v[k],comparison[k]) for k in ('absoluteSpread','fineMaxNorm','relativeSpread','target','withinTarget')) and v['difference'].shape==(129,129) and np.isfinite(v['difference']).all(),'full actual saved comparison packet; no subtraction repeated')
        journal.json('validated-rows/'+str(ri)+'.json',{'result':row_result,'fullNativeCacheArgumentsAndReceiptsJoined':True,'finiteFullRows':True,'scienceRecomputed':0});results.append(row_result)
        for path,rec in rc['consumedRoutes'].items():reader.retain(path,rec['sha256'])
        for rec in a.values():reader.retain(rec['path'],rec['sha256'])
    require(counts==c['counts'],'counts from actual saved branch receipts')
    for path,rec in c['consumedRoutes'].items():reader.retain(path,rec['sha256'])
    for rec in c['artifacts'].values():reader.retain(rec['path'],rec['sha256'])
    for path in (Path(__file__).resolve(),M/(NAME+'_plan.md')):
        reader.retain(path);dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out:out.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    reader.postcheck();journal.json('validated-paths.json',dict(reader.routes));checks={'status':'PASSED_BOUNDED_SAVED_MATERIAL_1D_ROWS_REVIEW','producerChecksSha256':inspection['outcome']['checksSha256'],'rows':results,'counts':counts,'consumedPaths':len(reader.routes),'allConsumedHashAndLinkIdentities':True,'scientificRecomputation':0,'wallSeconds':time.monotonic()-started,'scope':'Bounded actual saved cache/operand/unit/row/result review only. No evaluator, Fourier, product, contraction, source, rule, comparison, map or ancestral science replay.'}
    journal.json('checks.json',checks);signal.alarm(0);print((base/'checks.json').read_text(),end='')


if __name__=='__main__':main()
