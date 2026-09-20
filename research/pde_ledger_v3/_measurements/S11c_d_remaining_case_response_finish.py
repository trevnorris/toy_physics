#!/usr/bin/env python3
"""Validate saved case responses and join independently encoded transcripts."""
import argparse, ast, gc, hashlib, itertools, json, resource, shutil, signal, time
from pathlib import Path
import numpy as np
import S11c_d_remaining_case_response as h
from S11c_d_output_codec import PayloadEncoder, PayloadDecoder, decoded_lines, restore_emission_index

f,b,m,r=h.f,h.b,h.m,h.response
PLAN=f.M/'S11c_d_remaining_case_response_finish_plan.md'


def inventory(root):
    return {str(p.relative_to(root)):{'sha256':f.digest(p),'bytes':p.stat().st_size}
            for p in sorted(root.rglob('*')) if p.is_file()}


def equal(actual,expected,label):
    f.require(m.same(actual,expected),label)


def load(base,origin):
    inv=json.loads((origin.parent/'response_construct.invocation.json').read_text())
    equal(inv,json.loads((origin.parent/'active.json').read_text()),'actual final original outcome')
    f.require(inv['exitCode']==1 and inv['status']=='failed' and not (origin/'checks.json').exists(),'known incomplete aggregation')
    error=(origin.parent/'response_construct.stderr').read_text()
    f.require(error.endswith("ValueError: invalid shared-payload definition\n") and 'decoded_lines' in error,'actual recorded aggregation failure')
    old=json.loads((origin/'inputs.json').read_text());saved=inventory(origin)
    manifest={'runDirectory':str(base),'originalDirectory':str(origin),'sourceFiles':dict(old['sourceFiles']),
              'inputPackets':dict(old['inputPackets']),'copiedInputs':{},'originalArtifacts':saved,
              'settings':old['settings'],'input':old['input'],'originalOutcome':inv,
              'scope':old['scope'],'originalCopiedInputs':old['copiedInputs']}
    for n,v in old['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(origin/'source'/n)==v,('original current/frozen source',n))
    for n,v in old['inputPackets'].items():f.require(f.digest(Path(n))==v,('original input',n))
    for n,v in old['copiedInputs'].items():f.require(f.digest(origin/n)==v,('focused copy',n))
    for n,v in saved.items():
        target='original-combined.out' if n=='full.out' else 'production-inputs.json' if n=='inputs.json' else n
        p=base/target;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(origin/n,p)
        f.require(f.digest(p)==v['sha256'],('unchanged completed operand',n));manifest['copiedInputs'][target]=v['sha256']
    for name in ('active.json','response_construct.invocation.json','response_construct.stderr','response_construct.stdout'):
        p=origin.parent/name;manifest['inputPackets'][str(p)]=f.digest(p)
    for p in (Path(__file__).resolve(),PLAN):
        n=str(p.relative_to(f.ROOT));manifest['sourceFiles'][n]=f.digest(p)
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    f.save(base/'inputs.json',manifest)
    return old,manifest


def verify_case(base,label,manifest):
    target=base/'cases'/label;source=base/'interiors/cases'/label
    system=f.unpickle(target/'finite/finite-system.pickle');sol=f.unpickle(target/'finite/finite-solution.pickle');view=f.unpickle(target/'finite/observable.pickle')
    ends=f.unpickle(base/'boundary-cases'/label/'case-boundary.pickle');basis=f.unpickle(base/'interiors/accepted-finite-system.pickle')
    coefficient=f.unpickle(source/'interior-matrices.pickle');direct=h.unsplit_arrays(f.unpickle(source/'direct-unsplit.pickle'))
    rows=f.unpickle(source/'row-matrices.pickle');binding=f.unpickle(base/'interiors/accepted-bindings'/label/'case-binding.pickle');native=f.unpickle(source/'direct-native-cells.pickle')
    equal(system['channels'],ends['finite'],'complete actual end coordinates')
    equal(system['nodes'],basis['nodes'],'same field nodes');equal(system['derivativeMatrices'],basis['derivativeMatrices'],'all actual derivatives')
    equal(system['settings'],binding['binding']['settings'],'actual setting identity');equal(rows['settings'],system['settings'],'row settings')
    equal(coefficient['fieldUnits'],ends['fieldUnits'],'field units');equal(coefficient['equationUnits'],ends['rowUnits'],'equation units')
    f.require(set(rows['rows'])==set(range(len(binding['binding']['bound']['rows']))) and len(binding['grades']['termJoins'])==native['terms'],'full native row/source/term census')
    f.require(b.norm(direct['total']-native['total'])/(1+b.norm(native['total']))<1e-10,'original native cell assembly')
    _,prepare,join=h.finite_tail()
    a,rhs,original=prepare(target/'finite',direct['total'].copy(),ends['finite'],len(basis['nodes']),basis['derivativeMatrices'],basis['nodes'],direct['local'],basis['settings'],{}, {},sol['polynomialDerivativeResiduals'],sol['groups'])
    matrix_delta=b.norm((a-system['matrix'])/sol['rowScale'][:,None]);f.require(matrix_delta<1e-10,'complete boundary-replaced matrix')
    equal(rhs,system['rhs'],'actual incoming forcing')
    if label!=h.BASELINE:
        equal(original,system['unreplacedOperator'],'retained interior operator');equal(sol['groups'],[rows],'all original case rows retained')
        equal(system['sourceFiles'],manifest['sourceFiles'],'finite construction source pins');equal(system['operandHashes'],manifest['inputPackets'],'finite construction operand pins')
    actual=h.finite_view(system,sol,False)
    for k in actual:
        if k not in ('independentCoefficients','independentDifference','independentStatus'):equal(actual[k],view[k],('saved finite observable',k))
    if label!=h.BASELINE:
        equal(view['independentCoefficients']-sol['coefficients'],view['independentDifference'],'saved independent finite difference')
        f.require(b.norm((a@view['independentCoefficients']-rhs)/sol['rowScale'][:,None])<1e-9,'independent saved finite equation')
    f.require(sol['rank']==645 and len(sol['singularValues'])==645 and np.all(sol['singularValues']>0),'full finite rank')
    equal(sol['balancedCondition'],float(sol['singularValues'][0]/sol['singularValues'][-1]),'finite condition')
    del a,original,native,rows,binding,direct
    systems=f.unpickle(target/'continuum/coefficient-systems.pickle');solution=f.unpickle(target/'continuum/coefficient-solutions.pickle')
    channels=f.unpickle(target/'continuum/channel-response.pickle');packet=f.unpickle(target/'continuum/continuum-response.pickle')
    matrices,rhs=r.systems(coefficient,ends,basis);equal(matrices,systems['matrices'],'all original interior and endpoint grades');equal(rhs,systems['rhs'],'all actual forcing grades')
    for key,value in [('solve',solution),('response',channels)]:equal(packet[key],value,('complete packet component',key))
    for key in ('fieldUnits','rowUnits','currentUnit'):equal(packet[key],ends[key],('actual response unit',key))
    equal(packet['settings'],basis['settings'],'continuum settings')
    if label!=h.BASELINE:
        equal(packet['sourceFiles'],manifest['sourceFiles'],'continuum source pins');equal(packet['inputPackets'],manifest['inputPackets'],'continuum operand pins')
    residual=b.subtract(b.J.multiply(matrices,solution['coefficients']),rhs);equal(residual,solution['residual'],'literal coefficient equations')
    equal({g:v/solution['rowScale'][:,None] for g,v in residual.items()},solution['scaledResidual'],'scaled coefficient equations')
    independent=b.subtract(b.J.multiply(matrices,solution['independentCoefficients']),rhs)
    f.require(b.norm({g:v/solution['rowScale'][:,None] for g,v in independent.items()})<1e-9,'saved independent coefficient equations')
    equal(b.subtract(solution['coefficients'],solution['independentCoefficients']),solution['independentDifference'],'literal independent coefficient differences')
    previous={}
    for g in b.G:
        known=b.J.multiply(matrices,previous).get(g,np.zeros_like(rhs[g]));equal(rhs[g]-known,solution['forcing'][g],('baseline and mixed forcing',g));previous[g]=solution['coefficients'][g]
    mutation=matrices[(0,0)]@solution['mixedForcingMutation']-matrices[(1,1)]@solution['coefficients'][(0,0)]
    f.require(b.norm(mutation/solution['rowScale'][:,None])<1e-9 and b.norm(solution['mixedForcingMutation'])>0,'actual saved mixed omission solution')
    f.require(solution['rank']==645 and np.all(solution['singularValues']>0),'full continuum rank')
    equal(solution['condition'],float(solution['singularValues'][0]/solution['singularValues'][-1]),'continuum condition')
    # Recontract saved maps and current roots; no inverse, root or response solve.
    fields={g:np.stack([basis['derivativeMatrices'][0]@x[i*129:(i+1)*129] for i in range(5)]) for g,x in solution['coefficients'].items()}
    equal(fields,channels['fields'],'all field coefficients and derivatives')
    open_rows=[];incoming_current=[];outgoing_current=[];dip=[];dop=[]
    for offset,(end,index) in enumerate((('LEFT',0),('RIGHT',128))):
        e=ends['ends'][end];trace={g:x[:,index].copy() for g,x in fields.items()}
        for g in b.G:trace[g][:,2*offset:2*offset+2]-=e['incoming'][g]
        equal(b.J.multiply(e['outgoingInverse'],trace),channels['modalBoundary'][end],'full open/closed extraction')
        f.require(b.norm(b.subtract(b.J.multiply(e['outgoing'],channels['modalBoundary'][end]),trace))<1e-8,'all outgoing trace reconstruction')
        indices=[j for item in channels['labels'][end] if item['kind']=='open' for j in item['columns']]
        open_rows.append({g:x[indices] for g,x in channels['modalBoundary'][end].items()})
        incoming_current.append({g:-e['orientation']*x[5:,5:] for g,x in e['currents']['total'].items()})
        outgoing_current.append({g:e['orientation']*x[np.ix_(indices,indices)] for g,x in e['currents']['total'].items()})
        dip.append(e['incomingOriginPhase']);dop.append(e['outgoingOriginPhase'])
    equal({g:np.vstack([x[g] for x in open_rows]) for g in b.G},channels['openBoundaryScattering'],'complete open amplitude selection')
    equal(r.diagonal_series(dip),channels['incomingOriginPhase'],'actual incoming phases');equal(r.diagonal_series(dop),channels['outgoingOriginPhase'],'actual outgoing phases')
    equal(b.J.multiply(b.J.multiply(channels['outgoingOriginPhase'],channels['openBoundaryScattering']),channels['incomingOriginPhase']),channels['openOriginScattering'],'ordered phase application')
    equal(r.gram(channels['incomingOriginPhase'],r.diagonal_series(incoming_current)),channels['incomingCurrentOrigin'],'actual incoming current contraction')
    f.require(b.norm(b.subtract(r.gram(channels['outgoingOriginPhase'],channels['outgoingCurrentOrigin']),r.diagonal_series(outgoing_current)))<1e-8,'actual outgoing current contraction')
    for side in ('incoming','outgoing'):
        root=channels[side+'CurrentRoot'];f.require(b.norm(b.subtract(b.J.multiply(root,root),channels[side+'CurrentOrigin']))<1e-9,'saved current root equation')
        f.require(b.norm(b.subtract(root,r.adjoint(root)))<1e-9,'current Hermitian map')
    f.require(b.norm(b.subtract(b.J.multiply(channels['incomingCurrentRoot'],channels['incomingCurrentRootInverse']),{(0,0):np.eye(4)}))<1e-9,'saved current inverse equation')
    equal(b.J.multiply(b.J.multiply(channels['outgoingCurrentRoot'],channels['openOriginScattering']),channels['incomingCurrentRootInverse']),channels['fluxOriginScattering'],'full noncommuting current normalization')
    equal(r.open_flux(channels,packet['ratio']),packet['flux'],'all homotopy current denominators and coefficients')
    remainders=f.unpickle(target/'continuum/formal-remainders.pickle');equal(remainders,packet['remainders'],'complete formal remainder packet')
    for (eta,sigma),v in remainders.items():
        equal(b.evaluate(solution['coefficients'],eta,sigma),v['retained'],'actual retained formal polynomial')
        equal(v['direct']-v['retained'],v['difference'],'literal formal difference')
        q=b.evaluate(matrices,eta,sigma)@v['direct']-b.evaluate(rhs,eta,sigma)
        equal(q,v['directEquationResidual'],'saved formal direct equation');f.require(b.norm(q/solution['rowScale'][:,None])<1e-8,'formal direct residual')
        equal(b.norm(v['difference']),v['maximumReferenceFrame'],'formal diagnostic norm')
    return {'case':label,'finiteMatrixScaledDifference':matrix_delta,'finiteRank':sol['rank'],'continuumRank':solution['rank'],
            'finiteScaledResidual':b.norm(sol['scaledEquationResidual']),'continuumScaledResidual':b.norm(solution['scaledResidual']),
            'finiteCurrentRatios':view['totalCurrentRatio'].tolist(),'continuumMapResidual':b.norm(channels['residuals']),
            'newSolves':0,'newQuadratureNodes':0,'wholeNativeFinitePrefix':join}


def anchoring(base):
    saved=f.unpickle(base/'anchoring-comparisons.pickle');summary={}
    for density,record in saved.items():
        lab='LAB_HELD__'+density;mat='MATERIAL_ADVECTED__'+density
        le=f.unpickle(base/'boundary-cases'/lab/'case-boundary.pickle');me=f.unpickle(base/'boundary-cases'/mat/'case-boundary.pickle')
        equal(le['ends'],me['ends'],'same continuum coordinates');equal(le['finite'],me['finite'],'same finite coordinates')
        lv=f.unpickle(base/'cases'/lab/'finite/observable.pickle');mv=f.unpickle(base/'cases'/mat/'finite/observable.pickle')
        lc=f.unpickle(base/'cases'/lab/'continuum/channel-response.pickle');mc=f.unpickle(base/'cases'/mat/'continuum/channel-response.pickle')
        for name,key in [('finiteAmplitudeDifference','originScattering'),('finiteCurrentDifference','totalCurrentRatio'),('finiteFieldDifference','originFields')]:equal(mv[key]-lv[key],record[name],('actual anchoring difference',name))
        equal(b.subtract(mc['fluxOriginScattering'],lc['fluxOriginScattering']),record['continuumFluxAmplitudeDifferences'],'same-coordinate grade differences')
        summary[density]={k:b.norm(v) for k,v in record.items() if k!='scope'}
        summary[density]['evaluatedContinuumAmplitudeDifference']=b.norm(b.evaluate(record['continuumFluxAmplitudeDifferences'],.01,.001))
        summary[density]['continuumGradeDifferences']={str(g):b.norm(v) for g,v in record['continuumFluxAmplitudeDifferences'].items()}
    return summary


def aggregate(base,labels):
    encoder=PayloadEncoder();parts={};expected_hash=hashlib.sha256();raw_hash=hashlib.sha256();global_tags=set();global_keys=set();count=0
    with (base/'full.out').open('x') as output:
        for label in labels:
            path=base/'cases'/label/'continuum/full.out';raw_hash.update(path.read_bytes());tags=[];keys={};indices=0;digest=hashlib.sha256()
            for line in decoded_lines(path):
                tag,sep,body=line.rstrip('\n').partition(': ')
                f.require(tag not in global_tags,'global unique tag');global_tags.add(tag)
                if tag.endswith('_WRITE_KEYS') and not tag.startswith('PY_S11CD_METADATA_'):
                    keys={str(k):str(v) for k,v in r.grades._restore(body)}
                    f.require(len(keys)==len(set(keys.values())) and not global_keys&set(keys.values()),'actual disjoint export keys');global_keys.update(keys.values())
                if tag.endswith('_EMISSION_LINES') and not tag.startswith('PY_S11CD_METADATA_'):
                    index={str(k):v for k,v in r.grades._restore(body)};restored=restore_emission_index(index,tags)
                    f.require(all(v>0 for v in restored.values()),'actual local source lines');indices+=1
                encoded=encoder.encode(body);output.write(tag+sep+encoded+'\n')
                digest.update(line.encode());expected_hash.update(line.encode());tags.append(tag);count+=1
            f.require(indices==1,'one preserved local emission index per case')
            checks_path=path.parent/('checks.json' if label==h.BASELINE else 'emission-checks.json');checks=json.loads(checks_path.read_text())
            f.require(len(tags)==checks.get('tagCount',checks.get('tags')),'completed individual replay tag count')
            f.require(len(keys)==checks.get('writeKeys',len(checks.get('keys',{}))),'completed individual replay key count')
            if label!=h.BASELINE:equal(keys,checks['keys'],'exact original case export namespace')
            parts[label]={'sha256':f.digest(path),'decodedSha256':digest.hexdigest(),'tags':len(tags),'keys':len(keys),'metadataPaths':checks['metadataPaths'],
                          'emissionChecksSha256':f.digest(checks_path),'localEmissionIndexUnchanged':True}
    f.require(raw_hash.hexdigest()==f.digest(base/'original-combined.out'),'preserved original raw concatenation')
    originals=itertools.chain.from_iterable(decoded_lines(base/'cases'/label/'continuum/full.out') for label in labels)
    sentinel=object();actual_hash=hashlib.sha256();actual_count=0
    for actual,expected in itertools.zip_longest(decoded_lines(base/'full.out'),originals,fillvalue=sentinel):
        f.require(actual==expected,'literal original decoded payload identity');actual_hash.update(actual.encode());actual_count+=1
    f.require(actual_count==count and actual_hash.digest()==expected_hash.digest(),'complete global decoded stream')
    # Actual stream controls: bad first reference rejects; changing a decoded
    # physical payload changes its exact digest. Neither touches saved files.
    failed=False
    try:PayloadDecoder().decode("Tuple(Str('s11cdSharedPayloadReference'), Integer(0))")
    except ValueError:failed=True
    f.require(failed,'unresolved reference rejection')
    first=next(decoded_lines(base/'full.out'));mutated=first.replace(': ',': Integer(1) + ',1)
    f.require(hashlib.sha256(mutated.encode()).digest()!=hashlib.sha256(first.encode()).digest(),'actual changed payload detection')
    return {'tags':count,'keys':len(global_keys),'parts':parts,'decodedSha256':actual_hash.hexdigest(),
            'identicalDecodedPayloads':count,'metadataPaths':sum(p['metadataPaths'] for p in parts.values()),
            'referenceMutationRejected':failed,'payloadMutationDetected':True,'nativeCodecUnchanged':True,
            'scope':'Globally re-encoded storage references; every decoded case payload and local emission index unchanged.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);ap.add_argument('--resume-from',type=Path,required=True);args=ap.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False);origin=args.resume_from.resolve()
    old,manifest=load(base,origin);counts=json.loads((base/'case-inventory.json').read_text());labels=(h.BASELINE,*counts)
    results={}
    for label in labels:
        results[label]=verify_case(base,label,old);f.save(base/('validation-'+label+'.json'),results[label]);gc.collect()
    comparisons=anchoring(base)
    validation={'cases':results,'anchoring':comparisons,'sourceFiles':manifest['sourceFiles'],'originalArtifacts':manifest['originalArtifacts'],'newSolves':0,'individualEmissionsRepeated':0}
    f.save(base/'operand-validation.json',validation)
    combined=aggregate(base,labels);f.save(base/'aggregation-checks.json',combined)
    f.atomic_pickle(base/'remaining-case-response.pickle',{'cases':counts,'baseline':str(base/'cases'/h.BASELINE),'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':manifest['scope']})
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,('final current/frozen source',n))
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,('final original input',n))
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,('final copied artifact',n))
    equal(inventory(origin),manifest['originalArtifacts'],'every original artifact pre/post identity')
    artifacts={n:v for n,v in inventory(base).items() if not n.startswith('source/') and n not in ('inputs.json','checks.json')}
    checks={**manifest,'status':'COMPLETED_FOUR_CASE_RESPONSES','cases':counts,'validation':results,'anchoring':comparisons,
            'combinedOutput':combined,'artifacts':artifacts,'newFiniteCases':0,'newContinuumCases':0,'newQuadratureNodes':0,
            'completedNewFiniteCasesReused':3,'completedNewContinuumCasesReused':3,'individualEmissionsRepeated':0,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
