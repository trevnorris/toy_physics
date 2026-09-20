#!/usr/bin/env python3
"""Complete actual case mode spaces from saved current and isolated-root inputs."""
import argparse,ast,copy,json,re,resource,shutil,signal,time
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_remaining_case_currents as currents
import S11c_d_end_normalization_check as native
import S11c_d_end_normalization_remainder_check as remainder
f,engine=currents.f,currents.engine
CP=f.M/'S11c_d_remaining_case_currents_checkpoint.json'
PLAN=f.M/'S11c_d_remaining_case_modes_plan.md'
BASELINE='LAB_HELD__RHO4_CONSTANT'
ENDS=('REFERENCE','LEFT','RIGHT')


def same(a,b):
    if isinstance(a,np.ndarray):return isinstance(b,np.ndarray) and a.dtype==b.dtype and a.shape==b.shape and np.array_equal(a,b)
    if isinstance(a,dict):return isinstance(b,dict) and set(a)==set(b) and all(same(a[k],b[k]) for k in a)
    if isinstance(a,(tuple,list)):return type(a) is type(b) and len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return currents.same(a,b)


def node(path,name):return ast.dump(next(n for n in ast.parse(path.read_text()).body if getattr(n,'name',None)==name))


def retain(src,dst,manifest,expected=None):
    value=f.digest(src);f.require(expected is None or value==expected,('accepted packet hash',str(src)))
    dst.parent.mkdir(parents=True,exist_ok=True)
    f.require(not dst.exists(),'new immutable operand copy');shutil.copyfile(src,dst)
    f.require(f.digest(dst)==value,'copied packet byte identity')
    manifest['inputPackets'][str(src)]=value;manifest['copiedInputs'][str(dst.relative_to(Path(manifest['runDirectory'])))]=value


def load(base):
    cp=json.loads(CP.read_text());origin=Path(cp['runDirectory'])
    f.require(cp['status']=='ACCEPTED_FOUR_CASE_CURRENT_SOURCES' and f.digest(origin/'checks.json')==cp['checksSha256'],'accepted completed case currents')
    for n,v in cp['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(origin/'source'/n)==v,'accepted current source identity')
    for n,v in cp['inputPackets'].items():f.require(f.digest(Path(n))==v,'accepted current input identity')
    for n,v in cp['artifacts'].items():f.require(f.digest(origin/n)==v['sha256'],'accepted current artifact identity')
    manifest={'runDirectory':str(base),'sourceFiles':dict(cp['sourceFiles']),'inputPackets':{str(origin/'checks.json'):cp['checksSha256']},'copiedInputs':{},'input':json.loads((origin/'inputs.json').read_text())['input'],'nativeJoins':{},'scope':'Complete case modal/current inputs and subspaces; no scattering or continuum response yet.'}
    names=['remaining-case-end-sources.pickle','remaining-case-currents.pickle']
    names += [n for n in cp['artifacts'] if n.startswith('accepted-cases/')]
    names += ['cases/LAB_HELD__RHOBR_CONSTANT/right/pairing-raw-checks.pickle','cases/LAB_HELD__RHOBR_CONSTANT/right/pairing-retained-checks.pickle']
    for n in names:retain(origin/n,base/n,manifest,cp['artifacts'][n]['sha256'])
    sources=f.unpickle(base/'remaining-case-end-sources.pickle');values=f.unpickle(base/'remaining-case-currents.pickle')
    old={};producer=None
    for end in ENDS:
        path=f.M/('S11c_d_end_normalization_'+end.lower()+'_thickness_repair_checkpoint.json');oldcp=json.loads(path.read_text());root=Path(oldcp['runDirectory'])
        f.require(not oldcp['unaccountedResidualNormsAboveDiagnosticThreshold'],'accepted fully accounted original normalization')
        manifest['sourceFiles'][str(path.relative_to(f.ROOT))]=f.digest(path)
        frozen=root/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
        joins={n:node(frozen,n)==node(engine.HERE,n) for n in ('ModalCurrentSubspaces','AdjointCurrentMap','polynomial_terms')}
        f.require(all(joins.values()),'unchanged whole native modal constructors');manifest['nativeJoins'][end]=joins
        manifest['inputPackets'][str(frozen)]=f.digest(frozen)
        args=json.loads((root/'arguments.json').read_text());pcp_path=Path(args['pairing_checkpoint']);pcp_path=pcp_path if pcp_path.is_absolute() else f.ROOT/pcp_path
        pcp=json.loads(pcp_path.read_text());pair=Path(pcp['runDirectory'])/'complete.pickle'
        f.require(f.digest(pair)==oldcp['provenance']['pairingCacheSha256']==pcp['artifacts']['complete.pickle']['sha256'],'modal original complete pairing source join')
        manifest['sourceFiles'][str(pcp_path.relative_to(f.ROOT))]=f.digest(pcp_path)
        target=base/'accepted-modes'/end.lower()
        for kind in ('modal','adjoint'):retain(root/(kind+'.pickle'),target/(kind+'.pickle'),manifest,oldcp['artifacts'][kind+'.pickle']['sha256'])
        retain(pair,target/'pairing.pickle',manifest)
        modal,known=f.unpickle(target/'modal.pickle');adjoint,aknown=f.unpickle(target/'adjoint.pickle');pairing,pknown=f.unpickle(target/'pairing.pickle')
        old[end]={'modal':modal,'adjoint':adjoint,'known':known,'adjointKnown':aknown,'pairing':pairing,'pairingKnown':pknown,'checkpoint':oldcp}
        mpath=Path(args['manifest']);pmanifest=json.loads(mpath.read_text());transcript=Path(pmanifest['run_directory'])/'full.out'
        f.require(f.digest(mpath)==oldcp['provenance']['producerManifestSha256'] and f.digest(transcript)==oldcp['provenance']['producerTranscriptSha256'],'original all-case isolated roots')
        if producer is not None:f.require(producer==transcript,'one shared all-case root producer')
        producer=transcript;manifest['inputPackets'].update({str(mpath):f.digest(mpath),str(transcript):f.digest(transcript)})
    roots=read_roots(producer,tuple(values['cases']))
    (base/'native-roots').mkdir()
    for address,value in roots.items():f.atomic_pickle(base/'native-roots'/(address+'.pickle'),value)
    for path in (Path(__file__).resolve(),PLAN,CP,Path(native.__file__),Path(remainder.__file__),f.M/'S11c_d_variable_profile_development_input.json'):
        manifest['sourceFiles'][str(path.relative_to(f.ROOT))]=f.digest(path)
    for n,v in manifest['sourceFiles'].items():
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target);f.require(f.digest(target)==v,'frozen mode source')
    f.save(base/'inputs.json',manifest)
    return sources,values,old,roots,manifest


def read_roots(path,labels):
    prefixes={'PY_S11CD_END_SPECTRUM_INPUT_'+end+'_'+label.replace('__','_')+'_0':label+'__'+end for label in labels for end in ENDS}
    results={address:{'packet':{},'records':[]} for address in prefixes.values()}
    for line in native.decoded_lines(path):
        tag,sep,payload=line.partition(': ')
        if not tag.startswith('PY_S11CD_END_SPECTRUM_INPUT_'):continue
        matched=[prefix for prefix in prefixes if tag.startswith(prefix+'_')]
        if not matched:continue
        f.require(len(matched)==1,'unique actual case/end root address');prefix=matched[0];name=tag[len(prefix)+1:];value=results[prefixes[prefix]]
        if name in ('ROOT_COVERAGE','BOUND_CARRIERS','GRADE_ORIGIN','PHYSICAL_PENCIL'):
            f.require(name not in value['packet'],'unique root source key');value['packet'][name]=native._restore(payload)
        elif re.fullmatch(r'MODE_\d+_RECORD',name):
            f.require(int(name.split('_')[1])==len(value['records']),'complete native root record sequence')
            value['records'].append({str(k):v for k,v in native._restore(payload)})
    for value in results.values():
        f.require(set(value['packet'])=={'ROOT_COVERAGE','BOUND_CARRIERS','GRADE_ORIGIN','PHYSICAL_PENCIL'},'full original root operands')
        coverage={str(k):v for k,v in value['packet']['ROOT_COVERAGE']};records=value['records']
        expected={(i,s) for i in range(int(coverage['DISTINCT_ROOT_COUNT'])) for s in (-1,1)}
        f.require({(int(v['ROOT_DISK_INDEX']),int(v['NORMAL_LIFT_SIGN'])) for v in records}==expected and len(records)==len(expected),'all root/lift candidates retained')
        f.require(coverage['FINITE_POLYNOMIAL_ROOT_COVERAGE']==sp.true,'completed original root isolation');value['coverage']=coverage
    return results


def context(base,label,end,sources,values,manifest):
    c,balance,acoustic,inp,data=currents.context(base,label,end,sources,manifest);value=values['cases'][label][end]
    f.require(same(value['strong'],data['strong']) and same(value['profileBindings'],inp.limits),'actual current/pencil/limit address')
    a=label.split('__')[0];e={'REFERENCE':None,'LEFT':-sp.oo,'RIGHT':sp.oo}[end]
    for obj,key in ((c,'conservative'),(balance,'slab'),(acoustic,'acoustic')):obj.construct=currents.unchanged_cache(value[key],a,e)
    dims=engine.PHYSICAL_METADATA.dimensions
    for atom,unit in value['knownDimensions'].items():
        if atom in dims.known:f.require(dims.known[atom]==unit,'actual mode/current field units')
        else:dims.known[atom]=unit
    pair=engine.ClosedCurrentPairing(acoustic,a,e);pair.construct=currents.unchanged_cache(value['pairing'],a,e)
    return pair,inp,value


def joined_roots(pair,inp,root):
    m=pair.modes;algebraic,relation,_=m.analytic(pair.acoustic.strong)
    binding=inp.mapping(algebraic,relation,(m.k,m.q,m.eta,m.sigma));origin={m.eta:sp.S.Zero,m.sigma:sp.S.Zero} if pair.end is None else inp.origin
    f.require(dict(root['packet']['GRADE_ORIGIN'])==origin,'actual native independent grade origin')
    bound={str(k):v for k,v in root['packet']['BOUND_CARRIERS']}
    f.require(all(bound.get(str(k))==v for k,v in binding.items()),'actual original material/profile coefficients')
    binding.update(origin);physical=algebraic.xreplace(binding).applyfunc(sp.cancel)
    f.require(engine.carrier_fingerprint(engine.cas(physical))==root['packet']['PHYSICAL_PENCIL'],'actual full native physical pencil')
    return binding,physical,relation


def projection(pair,value):
    if isinstance(value,sp.MatrixBase):return value.applyfunc(pair.c.retained)
    if isinstance(value,(tuple,list)):return tuple(projection(pair,v) for v in value)
    return pair.c.retained(value)


def difference(a,b):
    if isinstance(a,(tuple,list)):return tuple(difference(x,y) for x,y in zip(a,b))
    return a-b


def derivative_packet(base,pair,value,inp,raw,manifest,progress):
    m=pair.modes;algebraic,relation,_=m.analytic(pair.acoustic.strong)
    mapping=inp.mapping(algebraic,relation,(m.k,m.q,m.eta,m.sigma,pair.r.omega))
    derived=pair.split_derivative_checks(value['pairing'],raw,mapping,progress)
    f.atomic_pickle(base/'pairing-derivatives.pickle',derived)
    combined={**raw,**derived}
    same_frequency=dict.fromkeys(pair.frequencies,pair.r.omega)
    combined['SLAB_CURRENT_DIAGONAL_JOIN_RESIDUAL']=(value['pairing']['SLAB_CURRENT_MATRIX'].xreplace(same_frequency)-value['slab']['SLAB_CURRENT_MATRIX']).applyfunc(sp.expand)
    kept={k:projection(pair,v) for k,v in combined.items() if k.endswith('_RESIDUAL')}
    rest={k:difference(combined[k],v) for k,v in kept.items()}
    packet={'result':value['pairing'],'checks':combined,'retained':kept,'remainders':rest,'bindings':mapping}
    target=base/'remainder-context';target.mkdir();f.atomic_pickle(target/'complete.pickle',(packet,dict(engine.PHYSICAL_METADATA.dimensions.known)))
    f.require(all(v==0 for v in currents.ends_source.scalars(kept)),'complete retained normal/frequency pairing and balance checks')
    cp={'runDirectory':str(target),'sourceFiles':{},'artifacts':{'complete.pickle':{'sha256':f.digest(target/'complete.pickle')}}}
    # The original remainder checker reads this local completed operand index;
    # it is not a production acceptance checkpoint.
    f.save(target/'pairing-packet-index.json',cp);f.save(target/'arguments.json',{'pairing_checkpoint':str(target/'pairing-packet-index.json')})
    f.save(target/'checks.json',{'provenance':{'pairingCacheSha256':f.digest(target/'complete.pickle')}})
    return packet,target


def construct_new(base,pair,value,inp,root,binding,manifest,progress):
    raw=f.unpickle(Path(manifest['runDirectory'])/'cases/LAB_HELD__RHOBR_CONSTANT/right/pairing-raw-checks.pickle')
    packet,rembase=derivative_packet(base,pair,value,inp,raw,manifest,progress)
    builder=engine.ModalCurrentSubspaces(pair,value['pairing'],binding)
    native_prepare=builder.prepare;prepared=[]
    def prepare_once():
        f.require(same(pair.construct(pair.anchoring,pair.end),value['pairing']),'unchanged complete preparation source')
        if prepared:return builder.coefficient_residuals
        result=native_prepare();prepared.append(True)
        f.atomic_pickle(base/'prepared-modal.pickle',{'symbolic':builder.symbolic_operands,'scalars':builder.scalar_operands,'residuals':result,'bindings':binding,'knownDimensions':dict(engine.PHYSICAL_METADATA.dimensions.known)})
        return result
    builder.prepare=prepare_once
    modal=builder.construct(root['records'],root['coverage'],float(inp.parameters['omega']),float(inp.parameters['W_0']),progress)
    modal.update(COEFFICIENT_RESIDUALS=builder.coefficient_residuals,SYMBOLIC_OPERANDS=builder.symbolic_operands,SCALAR_OPERANDS=builder.scalar_operands,NATIVE_RECORDS=root['records'],NATIVE_COVERAGE=root['coverage'])
    f.atomic_pickle(base/'modal.pickle',(modal,dict(engine.PHYSICAL_METADATA.dimensions.known)))
    adjoint=engine.AdjointCurrentMap(builder).construct(modal,progress)
    f.atomic_pickle(base/'adjoint.pickle',(adjoint,dict(engine.PHYSICAL_METADATA.dimensions.known)))
    f.require(len(prepared)==1,'one actual symbolic preparation')
    accounting=remainder.compute(rembase,pair,builder,modal,adjoint,binding)
    f.atomic_pickle(base/'remainders.pickle',accounting)
    summary=accounting['summary']
    f.require(all(v['nonzeroScalars']==0 for v in summary['algebraicChecks'].values()),'exact remainder decomposition and retained grades')
    f.require(not any(a<=1 and b<=1 for values in summary['coefficientGradesEtaSigma'].values() for a,b in values),'no retained remainder coefficient')
    f.require(all(v<=1e-8 for v in summary['residualMinusRemainderNormMaxima'].values()),'actual modal/adjoint remainder contractions')
    validate_modes(modal,adjoint,root,summary['residualMinusRemainderNormMaxima'])
    units={}
    for mode in modal['RECORDS']:
        for group in ('FORMS','RESIDUALS','OPERANDS'):
            for key,v in mode[group].items():
                array=np.asarray(v);units[(mode['INDEX'],group,key)]=tuple(builder.tensor_unit(group,key,(i,),mode['NULLITY']) for i in range(array.size))
    f.atomic_pickle(base/'modal-array-units.pickle',units)
    return modal,adjoint,{'pairingRetainedScalars':len(tuple(currents.ends_source.scalars(packet['retained']))),'remainderSummary':summary,'onePreparation':True}


def validate_modes(modal,adjoint,root,accounted):
    f.require(len(modal['RECORDS'])==len(adjoint['RECORDS'])==len(root['records']),'complete candidate subspaces')
    f.require(all(v==0 for v in currents.ends_source.scalars(modal['COEFFICIENT_RESIDUALS'])),'actual quadratic current extraction')
    f.require(all(v==0 for v in currents.ends_source.scalars(adjoint['SYMBOLIC_RESIDUALS'])),'actual adjoint product rules')
    for mode,dual,source in zip(modal['RECORDS'],adjoint['RECORDS'],root['records']):
        for key in ('ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','NULLITY'):f.require(mode[key]==dual[key]==int(source[key]),'actual candidate/rank identity')
        n=mode['NULLITY']
        for side in ('RIGHT','LEFT'):f.require(np.linalg.matrix_rank(mode['FORMS'][side],tol=1e-9)==n,'full degenerate left/right basis')
        for group in ('FORMS','OPERANDS','RESIDUALS'):
            for key,value in mode[group].items():f.require(np.isfinite(value).all(),('finite full mode array',key))
        for key,value in mode['RESIDUALS'].items():
            norm=float(np.linalg.norm(value));f.require(norm==mode['RESIDUAL_NORMS'][key],'literal modal norm replay')
            if norm>1e-8:f.require('MODAL_'+key in accounted,('unexplained modal residual',key,norm))
        for item in dual['ITEMS']:
            f.require(np.isfinite(item['VALUE']).all(),'finite full adjoint array')
            if item['GROUP']=='RESIDUALS':
                norm=float(np.linalg.norm(item['VALUE']));f.require(norm==dual['RESIDUAL_NORMS'][item['NAME']],'literal adjoint norm replay')
                if norm>1e-8:f.require('ADJOINT_'+item['NAME'] in accounted,('unexplained adjoint residual',item['NAME'],norm))
    native.check_certificate(modal['NORMAL_REALITY_COVERAGE'],modal['RECORDS'])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    def progress(info):
        with (base/'progress.jsonl').open('a') as stream:stream.write(json.dumps({'wallSeconds':time.monotonic()-start,**info})+'\n')
    sources,values,old,roots,manifest=load(base);inventory={};new=[];prepared={};controls=[]
    for label,ends in values['cases'].items():
        for end,value in ends.items():
            address=label+'__'+end;target=base/'cases'/label/end.lower();target.mkdir(parents=True,exist_ok=True)
            progress({'stage':'input_joins','address':address})
            pair,inp,value=context(base,label,end,sources,values,manifest);root=roots[address]
            binding,physical,relation=joined_roots(pair,inp,root)
            f.atomic_pickle(target/'mode-inputs.pickle',{'case':label,'end':end,'binding':binding,'physical':physical,'relation':relation,'rootPacket':root,'pairing':value['pairing'],'acoustic':value['acoustic'],'fieldUnits':pair.c.field_units,'frequency':inp.parameters['omega'],'depth':inp.parameters['W_0']})
            candidate=old[end];reusable=same(value['pairing'],candidate['pairing']['result'])
            if reusable:
                coordinates=(*pair.frequencies,*pair.c.leg_momenta,*pair.acoustic.qlegs,*pair.c.fields[0],*pair.c.fields[1])
                dims=engine.PHYSICAL_METADATA.dimensions
                f.require(all(atom in candidate['pairingKnown'] and dims.known[atom]==candidate['pairingKnown'][atom] for atom in coordinates),'full original mode/current coordinate units')
                f.require(same(value['acoustic'],values['cases'][BASELINE][end]['acoustic']),'complete original acoustic source before mode reuse')
                f.require(same(root['records'],candidate['modal']['NATIVE_RECORDS']) and same(root['coverage'],candidate['modal']['NATIVE_COVERAGE']),'exact original full candidate and coverage inputs')
                for key in ('PENCIL_PLUS','PENCIL_MINUS'):f.require(same(candidate['modal']['SYMBOLIC_OPERANDS'][key],value['pairing']['CLOSED_PENCIL_LEGS'][0 if key.endswith('PLUS') else 1]),'actual complete modal pencil source')
                for rec in candidate['modal']['RECORDS']:
                    f.require(complex(rec['OMEGA'])==complex(inp.parameters['omega']) and complex(rec['DEPTH_CUTOFF'])==complex(inp.parameters['W_0']),'actual modal numerical frequency/depth')
                # Exact prepared source, candidate and binding joins establish reuse;
                # all original remainder evidence remains in its accepted checkpoint.
                for kind in ('modal','adjoint'):retain(base/'accepted-modes'/end.lower()/(kind+'.pickle'),target/(kind+'.pickle'),manifest)
                result=candidate['modal'];mode='accepted-original-mode-and-adjoint'
            else:
                new.append(address);mode='new-current-family-pending';result=None
            inventory[address]={'mode':mode,'candidateCount':len(root['records']),'currentFamily':value['family'],'rootInputSha256':f.digest(target/'mode-inputs.pickle')}
            if address=='LAB_HELD__RHOBR_CONSTANT__RIGHT':
                for key in ('GRADE_ORIGIN','PHYSICAL_PENCIL'):
                    changed=copy.deepcopy(root);changed['packet'][key]=() if key=='GRADE_ORIGIN' else 'wrong-pencil-address'
                    try:joined_roots(pair,inp,changed)
                    except ValueError:controls.append(key)
                    else:raise ValueError('wrong root input accepted')
            f.save(base/'mode-inventory.json',inventory)
    f.require(len(inventory)==12 and len(controls)==2,'complete case/end joins and actual wrong-input controls')
    f.require(set(new)=={'LAB_HELD__RHOBR_CONSTANT__RIGHT','MATERIAL_ADVECTED__RHOBR_CONSTANT__RIGHT'},'actual two unmatched current addresses')
    lab=base/'cases/LAB_HELD__RHOBR_CONSTANT/right';mat=base/'cases/MATERIAL_ADVECTED__RHOBR_CONSTANT/right'
    li,mi=(f.unpickle(q/'mode-inputs.pickle') for q in (lab,mat))
    for key in ('binding','physical','relation','rootPacket','pairing','acoustic','fieldUnits','frequency','depth'):f.require(same(li[key],mi[key]),('actual shared new modal inputs',key))
    f.save(base/'inputs.json',manifest)
    pair,inp,value=context(base,'LAB_HELD__RHOBR_CONSTANT','RIGHT',sources,values,manifest)
    progress({'stage':'new_family_construction'})
    owner='LAB_HELD__RHOBR_CONSTANT__RIGHT'
    modal,adjoint,summary=construct_new(lab,pair,value,inp,roots[owner],li['binding'],manifest,progress)
    for kind in ('modal','adjoint'):retain(lab/(kind+'.pickle'),mat/(kind+'.pickle'),manifest)
    for address in new:inventory[address].update(mode='new-full-current-family' if address==owner else 'exact-shared-full-current-family',basisDirections=sum(v['NULLITY'] for v in modal['RECORDS']),physicalCurrentDirections=sum(v['NULLITY'] for v in modal['RECORDS'] if v.get('PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED',False)))
    f.atomic_pickle(base/'remaining-case-modes.pickle',{'inventory':inventory,'newFamily':summary,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets']})
    f.save(base/'mode-inventory.json',inventory);f.save(base/'inputs.json',manifest)
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'current/frozen mode source identity')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'mode input pre/post identity')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'mode copy pre/post identity')
    artifacts={str(z.relative_to(base)):currents.ends_source.binding.factors.artifact(z) for z in base.rglob('*.pickle') if 'source' not in z.relative_to(base).parts}
    checks={**manifest,'status':'COMPLETED_FOUR_CASE_MODE_SUBSPACES','inventory':inventory,'newFamily':summary,'newModeConstructions':1,'mutationControls':controls,'artifacts':artifacts,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
