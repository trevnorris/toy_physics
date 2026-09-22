#!/usr/bin/env python3
"""Bounded saved LAB response review; no numerical response is recomputed."""
import argparse
import ast
import json
from pathlib import Path
import resource
import signal
import time
import S11c_d_remaining_case_frequency_lab_operator as producer
import S11c_d_remaining_case_frequency_lab_operator_recover as recovery

saved,io,np,sp=producer.saved,producer.io,producer.np,producer.sp
M,F,require,same=producer.M,producer.F,producer.require,producer.same
NAME='S11c_d_remaining_case_frequency_lab_operator_review'
ROOT=F/'lab-operator-recovery-01'
CHECKS_SHA='6b599efce28034bb45a7c5542cb509c0c5d97955d9863e2f04af8c612215f003'
PRODUCER_SHA='06aa307c7b43c97db370cbebbb7a179042d41a374203ca16ee76bc31bf5869db'
RECOVERY_SHA='c1b8fa71c4b5e9c666bf85c81707aef023aa1eae0068381af9e6aedf0f38ed27'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True)
    base=ap.parse_args().run_directory.resolve();base.relative_to(io.REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();io.digest=saved.digest
    reader,journal=saved.Reader(),producer.inputs.metadata.MetadataJournal(base)
    read_packets=set()
    reader.retain(producer.__file__,PRODUCER_SHA);reader.retain(recovery.__file__,RECOVERY_SHA)
    def retain(rec):return reader.retain(rec.get('logical',rec.get('path')),rec['sha256'])
    def packet(rec):
        route=retain(rec);read_packets.add(route['canonical'])
        return reader.packet(rec.get('logical',rec.get('path')),rec['sha256'])
    def meta(rec):return reader.json(rec.get('logical',rec.get('path')),rec['sha256'])
    checks=reader.json(ROOT/'complete/checks.json',CHECKS_SHA);arts=checks['artifacts']
    require(checks['status']=='COMPLETED_BOUNDED_LAB_FREQUENCY_RESPONSE' and checks['completedRowOrEndCalls']==0,'completed single response only')
    outcome=reader.json(ROOT/'completion-inspection.json')
    require(outcome['actualExits']==[0,0,0] and outcome['checksStdoutIdentity'] and outcome['checksSha256']==CHECKS_SHA and outcome['zeroCapOOMSwap'],'actual final clean producer')
    for rec in outcome['evidence'].values():retain(rec)
    for rec in checks['consumedRoutes'].values():
        actual=retain(rec);require(actual==rec,'exact consumed logical/canonical/hash/size/link identity')
    def read(name):return packet(arts[name])
    def info(name):return meta(arts[name])
    def array(value,shape):
        require(isinstance(value,np.ndarray) and value.dtype==np.dtype(complex) and value.shape==shape and np.isfinite(value).all(),'actual finite full saved numerical array')
    def no_science(*args,**kwargs):raise RuntimeError('saved response review disables new science')
    producer.main=recovery.main=producer.prep.main=producer.prep.evaluator=no_science
    journal.write=no_science
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=no_science
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=no_science
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,no_science)
    # All scientific numpy/scipy calls used by the constructor are disabled; reading classes/codec stays native.
    for name in ('svd','solve','lstsq'):setattr(np.linalg,name,no_science)
    producer.la.lu_factor=producer.la.lu_solve=no_science
    caller=info('native-and-numerical-callers.json')
    for entry in caller['native'].values():
        for name in ('current','frozen'):retain(entry[name])
        tree=ast.parse(Path(entry['current']['canonical']).read_text())
        for name,body in entry['bodies'].items():require(ast.dump(next(n for n in tree.body if getattr(n,'name',None)==name))==ast.dump(ast.parse(body).body[0]),'complete native source body')
    for rec in caller['linearAlgebraSources'].values():retain(rec)
    retain(caller['nativeNorm']['file'])
    actual_text=Path(caller['native']['S11c_d_frequency_matrix.py']['current']['canonical']).read_text()
    module,join=producer.solve_adapter(actual_text);compile(module,'<review compile only>','exec')
    require(same(join,caller['solveAdapter']),'whole actual adapted native solve')
    compile(recovery.adapted(Path(producer.__file__).read_text()),'<review continuation compile only>','exec')
    journal.json('validated-native-and-continuation.json',{'nativeCaller':arts['native-and-numerical-callers.json'],'wholeSolveReverseAST':True,'wholeRecoveryReverseAST':True,'newScientificCalls':0})
    baseline='LAB_HELD__RHO4_CONSTANT';own='LAB_HELD__RHOBR_CONSTANT';case_data={}
    for label in (baseline,own):
        prefix='cases/'+label;arg=read(prefix+'/operator-input.pickle');sc=packet(arg['scalar']);ctx=packet(arg['physicalRoutes']['context']);graded=packet(arg['independentGrades']);ends=packet(arg['boundaryUnitFrame']);system=packet(arg['trialSystem']);seed=packet(arg['fixedScales'])
        for rec in arg['physicalRoutes'].values():
            if isinstance(rec,dict) and 'sha256' in rec:retain(rec)
        require(arg['case']==label and same(arg['frequency'],sc['frequency']) and complex(arg['frequency'])==1-.01j and same(*ctx['contextPair']) and same(arg['settings'],system['settings']),'whole actual case setting/frequency/context')
        require(tuple(map(tuple,ends['fieldUnits']))==tuple(map(tuple,sc['fieldUnits'])) and tuple(map(tuple,ends['rowUnits']))==tuple(map(tuple,sc['equationUnits'])),'own physical coefficient/end unit frame')
        require(set(v['rowIndex'] for v in arg['rows'])==set(range(checks['rows'][label])),'all complete own row routes')
        for row in arg['rows']:retain(row['input']);retain(row['value']['packet']);meta(row['savedFullRowInputJoin'])
        routes=info(prefix+'/end-map-consumer-routes.json');maps={}
        for side,route in routes.items():
            value=packet(route['actualMap']['file'])
            for key in route['actualMap']['keys']:value=value[key]
            selection=packet(route['selection']);full=packet(route['physicalInput']);meta(route['existingWholeCallerJoin'])
            require(full['address']==(label,side) and route['openIndices']==[i for i,x in enumerate(selection['channel']['outgoing']) if x['kind']=='open'],'own actual physical open end addresses')
            for name,shape in [('trace',(5,5)),('forcingAtCommonOrigin',(5,2)),('observationAtCommonOrigin',(5,5)),('directIncomingSubtraction',(5,2))]:array(value[name],shape)
            maps[side]=dict(value,openIndices=route['openIndices'])
        case_data[label]={'input':arg,'scalar':sc,'context':ctx['contextPair'][0],'ends':ends,'system':system,'seed':seed,'maps':maps,
                          'joins':{tuple(v['address']):v for v in sc['joins']}}
        journal.json('cases/'+label+'/validated-physical-input.json',{'input':arts[prefix+'/operator-input.pickle'],'rows':len(arg['rows']),'ownEndRoutes':arts[prefix+'/end-map-consumer-routes.json'],'gradeSource':arg['independentGrades']})
    coef_values={};coef_inputs={}
    names=sorted((n for n in arts if n.startswith('coefficient-values/') and n.endswith('/completed.json')),key=lambda n:int(n.split('/')[1]))
    for name in names:
        prefix=name.rsplit('/',1)[0];receipt=info(name);require(receipt=={'input':arts[prefix+'/input.pickle'],'value':arts[prefix+'/value.pickle']},'actual coefficient completion hashes')
        ar=read(prefix+'/input.pickle');value=read(prefix+'/value.pickle');request=ar['call'];case=case_data[own];adr=ar['address'];ctx=case['context']
        require(same(request['expression'],case['scalar']['actual'][adr]) and same(request['unit'],case['joins'][adr]['unit']) and
                same(request['environment'][ctx['z']],case['system']['nodes']) and request['environment'][ctx['regulator']]==case['system']['settings']['regulator'], 'actual coefficient expression/unit/environment before call')
        array(value,(129,));coef_values[receipt['value']['path']]=value;coef_inputs[receipt['value']['path']]=request
        journal.json('validated-'+prefix+'.json',{'input':receipt['input'],'value':receipt['value'],'address':list(adr),'finiteFullSavedReturn':True,'recomputed':False})
    require(len(names)==checks['counts']['newCoefficientCalls']==66,'actual new coefficient receipts')
    new_blocks=saved_blocks=products=0;blocks={};comparison_count=0
    for label in (baseline,own):
        matches=[]
        for kind in ('local','nonlocal'):
            for i in range(5):
                for j in range(5):
                    prefix='cases/'+label+'/blocks/'+kind+f'/{i}-{j}';ar=read(prefix+'/input.pickle');a,b=ar['items'],ar['baselineItems']
                    match=len(a)==len(b) and all(same({k:v for k,v in x.items() if k!='operandAddress'},
                             {k:v for k,v in y.items() if k!='operandAddress'}) for x,y in zip(a,b));matches.append([kind,i,j,match])
                    require(ar['block']==(kind,i,j),'whole actual block index')
                    if label==baseline:require(match,'each actual baseline whole block input matches')
                    else:
                        receipt=info(prefix+'/completed.json');require(receipt['input']==arts[prefix+'/input.pickle'] and receipt['block']==[kind,i,j],'block completed source address')
                        if match:
                            require(receipt['disposition']=='ACCEPTED_COMPLETE_BLOCK','actual full saved block branch')
                            old=packet(receipt['value']['packet']);rs=receipt['value']['rowSlice'];cs=receipt['value']['columnSlice'];require(rs==[i*129,(i+1)*129] and cs==[j*129,(j+1)*129] and receipt['value']['keys']==[kind],'original full block return address')
                            value=old[kind][slice(*rs),slice(*cs)];saved_blocks+=1
                        else:
                            require(receipt['disposition']=='NEW_ORDERED_NATIVE_BLOCK' and receipt['value']==arts[prefix+'/value.pickle'],'actual missing full block branch')
                            value=read(prefix+'/value.pickle');new_blocks+=1
                            for ti,item in enumerate(a):
                                tp=prefix+'/terms/'+str(ti);tr=info(tp+'/completed.json');ar2=read(tp+'/input.pickle');pv=read(tp+'/value.pickle');acc=read(tp+'/accumulated.pickle')
                                require(tr=={'input':arts[tp+'/input.pickle'],'product':arts[tp+'/value.pickle'],'accumulated':arts[tp+'/accumulated.pickle']} and
                                        ar2['blockInput']==arts[prefix+'/input.pickle'] and ar2['item']==ti and same(ar2['operand'],item['operandAddress']),'actual ordered native product/accumulation routes')
                                cr=ar2['coefficient'];request=coef_inputs[cr['path']];require(same(request['expression'],item['coefficient']) and same(request['unit'],item['unit']),'exact completed new scalar input reuse')
                                array(pv,(129,129));array(acc,(129,129));products+=1
                            if a:require(same(acc,value),'actual last saved ordered accumulation is full block value')
                            comparisons=read(prefix+'/selected-direct-comparison.pickle');require(len(comparisons)==3,'selected direct scalar comparisons')
                            for v in comparisons:
                                require(v['assembled']==value[v['row'],v['column']] and v['absolute']<1e-11*(1+abs(v['direct'])),'actual saved scalar comparison result')
                            comparison_count+=1
                        array(value,(129,129));blocks[kind,i,j]=value
                    journal.json('validated-blocks/'+label+'/'+kind+f'/{i}-{j}.json',{'input':arts[prefix+'/input.pickle'],'fullInputMatchesSavedBaseline':match,'scienceRecomputed':False})
        wanted=info('cases/'+label+'/whole-baseline-operator-match.json')
        require(wanted['blockMatches']==matches and wanted['all50Blocks']==all(x[-1] for x in matches),'actual whole caller match aggregate')
    require(new_blocks==27 and saved_blocks==23,'actual27new/23whole saved blocks')
    prefix='cases/'+own;interior=read(prefix+'/frequency-interior.pickle')
    for kind in ('local','nonlocal'):
        array(interior[kind],(645,645))
        for i in range(5):
            for j in range(5):require(same(interior[kind][i*129:(i+1)*129,j*129:(j+1)*129],blocks[kind,i,j]),'all saved full blocks occur at actual whole matrix addresses')
    array(interior['total'],(645,645));journal.json('validated-interior.json',{'input':arts[prefix+'/interior-input.pickle'],'value':arts[prefix+'/frequency-interior.pickle'],'newBlocks':new_blocks,'savedBlocks':saved_blocks,'products':products,'selectedComparisons':comparison_count,'sumNotRecomputed':True})
    solve=prefix+'/solve';sy=read(solve+'/frequency-system.pickle');result=read(solve+'/frequency-solution.pickle');summary=info(prefix+'/response-summary.json')
    require(sy['rank']==result['rank']==summary['rank']==645 and result['fixedFrameCondition']==summary['fixedFrameCondition'],'actual native regular finite rank/condition')
    for name,shape in [('matrix',(645,645)),('balanced',(645,645)),('rhs',(645,4))]:array(sy[name],shape)
    for name,shape in [('coefficients',(645,4)),('independentDifference',(645,4)),('scaledEquationResidual',(645,4)),('fields',(5,129,4)),('openOriginScattering',(4,4))]:array(result[name],shape)
    require(same(sy['fixedRowScale'],case_data[own]['seed']['fixedRowScale']) and same(sy['fixedColumnScale'],case_data[own]['seed']['fixedColumnScale']),'actual existing fixed seed scales')
    require(summary['maximumScaledEquationResidual']<1e-9 and all(v<1e-9 for v in summary['maximumBoundaryResidual'].values()),'saved equation/boundary checks within declared native scale')
    local={}
    for obs in join['localObservations']:
        rec=arts[solve+'/locals/'+str(obs['site'])+'.pickle'];v=packet(rec);local[str(obs['site'])]=v
        require(set(v)<=set(obs['names']),'actual native assignment/whole-loop local names')
        journal.json('validated-locals/'+str(obs['site'])+'.json',{'actualLocal':rec,'names':list(v),'arithmeticRecomputed':False})
    require(same(local['3']['matrix'],interior['total']) and same(local['6']['matrix'],sy['matrix']) and same(local['6']['rhs'],sy['rhs']), 'actual matrix copy and boundary loop outputs')
    require(same(local['13']['system'],sy) and same(local['26']['result'],result),'whole native saved final packet joins')
    for key,site in [('coefficients','16'),('residual','18'),('difference','19'),('fields','21')]:
        output={'residual':'scaledEquationResidual','difference':'independentDifference'}.get(key,key)
        require(same(local[site][key],result[output]),'saved assignment -> final native result')
    for i,name in enumerate(('np.linalg.svd','la.lu_factor','la.lu_solve')):
        op=solve+'/operations/'+str(i);ar=read(op+'/input.pickle');val=read(op+'/value.pickle');receipt=info(op+'/completed.json')
        require(ar['function']==name and ar['site']==str(i) and receipt=={'input':arts[op+'/input.pickle'],'value':arts[op+'/value.pickle']},'full native call/input/value receipt')
        if i==0:
            require(same(ar['args'],(sy['balanced'],)) and ar['kwargs']=={'full_matrices':False} and same(val,(local['10']['u'],local['10']['s'],local['10']['vh'])),'actual whole native SVD arguments/return')
        elif i==1:
            require(same(ar['args'],(sy['balanced'],)) and ar['kwargs']=={},'actual native LU factor arguments');lu=val
        else:require(same(ar['args'][0],lu) and ar['kwargs']=={},'actual native LU factor -> solve argument')
        journal.json('validated-operations/'+str(i)+'.json',{'input':receipt['input'],'value':receipt['value'],'function':name,'recomputed':False})
    for side in ('LEFT','RIGHT'):
        array(result['modalAmplitudes'][side],(5,4));array(result['boundaryResiduals'][side],(5,4))
        require(same(result['modalAmplitudes'][side],local['25']['modal'][side]) and same(result['boundaryResiduals'][side],local['25']['boundary_residuals'][side]),'whole native boundary/observation outputs')
    reuse=info('cases/'+baseline+'/completed-response-reuse.json');oldresult=packet(reuse['solution']);oldsystem=packet(reuse['system']);oldinterior=packet(reuse['interior'])
    require(reuse['newScientificCalls']==0 and oldsystem['rank']==645 and oldresult['rank']==645,'whole completed baseline response preserved')
    contrast=read('LAB-response-contrast.pickle');require(same(contrast['own'],result['openOriginScattering']) and same(contrast['baseline'],oldresult['openOriginScattering']),'actual own/baseline amplitude contrast inputs');array(contrast['difference'],(4,4))
    journal.json('validated-system-and-response.json',{'system':arts[solve+'/frequency-system.pickle'],'response':arts[solve+'/frequency-solution.pickle'],
        'nativeInputs':arts[solve+'/input.pickle'],'actualSummary':summary,'baselineReuse':arts['cases/'+baseline+'/completed-response-reuse.json'],
        'newScientificCalls':0,'interpretation':'Saved native before/after packets and unchanged full caller; no independent response recomputation or observable error bound.'})
    # Read any remaining output packets once, retain all exact paths and bytes.
    for name,rec in arts.items():
        r=retain(rec)
        if name.endswith('.pickle') and r['canonical'] not in read_packets:
            value=packet(rec);del value
    for path in (Path(__file__),M/(NAME+'_plan.md')):
        reader.retain(path);dest=base/'source'/path.name;dest.parent.mkdir(exist_ok=True)
        with dest.open('xb') as out:out.write(path.read_bytes())
        reader.retain(dest,saved.digest(path))
    journal.json('validated-paths.json',dict(reader.routes));reader.postcheck()
    result_checks={'status':'PASSED_BOUNDED_SAVED_LAB_FREQUENCY_RESPONSE_REVIEW','producerChecksSha256':CHECKS_SHA,'newScientificCalls':0,
        'coefficientCalls':66,'newBlocks':new_blocks,'savedBlocks':saved_blocks,'productReceipts':products,'nativeLinearAlgebraCalls':3,'rank':645,
        'consumedPaths':len(reader.routes),'wallSeconds':time.monotonic()-started,'artifacts':dict(journal.artifacts),
        'scope':'One saved finite complex LAB response. No response/array arithmetic/science replay; observable precision remains to be assessed by selected comparisons.'}
    journal.json('checks.json',result_checks);signal.alarm(0);print((base/'checks.json').read_text(),end='')


if __name__=='__main__':main()
