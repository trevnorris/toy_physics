#!/usr/bin/env python3
"""Saved-input routing for the unchanged native case-current coordinator."""
import ast,copy,json,shutil,sys,time,resource,signal
from pathlib import Path
import S11c_d_remaining_case_currents as h
f,engine,same=h.f,h.engine,h.same
CHECKPOINT=f.M/'S11c_d_remaining_case_currents_preflight.json'
PLAN=f.M/'S11c_d_remaining_case_currents_production_plan.md'
STATE={}


def signature(data):
    return {'strong':data['strong'],'massSymbol':data['strong'][3,:],
        'massRow':data['ZERO_TRANSFER_MASS_ROW'],'energy':data['TANGENTIAL_ENERGY_REDUCTION'],
        'constraint':data['MATERIAL_CONSTRAINT_COEFFICIENTS']}


def check_pair_inputs(proofs,left,right):
    f.require(set(proofs)==set(left)==set(right),'complete acoustic input signature')
    for key,proof in proofs.items():
        f.require(same(left[key],proof['left']) and same(right[key],proof['right']),('saved actual acoustic proof operands',key))
        zero=all(v==0 for v in h.ends_source.scalars(proof['expandedDifference']))
        f.require(zero==proof['allZero'],'saved normalized proof status')
    return all(proof['allZero'] for proof in proofs.values())


def accepted_signature(old):
    return {'strong':old['strong'],'massSymbol':old['strong'][3,:],
        'massRow':old['conservative']['ZERO_TRANSFER_MASS_ROW'],
        'energy':old['conservative']['TANGENTIAL_ENERGY_REDUCTION'],
        'constraint':old['conservative']['MATERIAL_CONSTRAINT_COEFFICIENTS']}


def prepare(base):
    cp=json.loads(CHECKPOINT.read_text());origin=Path(cp['runDirectory']);proofroot=Path(cp['acceptanceDirectory'])
    f.require(cp['status']=='ACCEPTED_FOCUSED_CURRENT_INPUT_REUSE','accepted complete current reuse check')
    f.require(f.digest(origin/'checks.json')==cp['checksSha256'],'focused checks identity')
    f.require(f.digest(Path(cp['validationScript']))==cp['validationScriptSha256'],'saved operand proof method')
    sources,accepted,paired,manifest=h.load(base)
    f.require(cp['sourceFiles']==manifest['sourceFiles'],'unchanged focused/native helper sources')
    f.require(cp['inputPackets']==manifest['inputPackets'] and cp['copiedInputs']==manifest['copiedInputs'],'same original consumed current inputs')
    copies={};inputs={}
    def retain(source,dest):
        value=f.digest(source);dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists():f.require(f.digest(dest)==value,'existing focused copy identity')
        else:shutil.copyfile(source,dest)
        f.require(f.digest(dest)==value,'byte-identical focused operand reuse')
        inputs[str(source)]=value;copies[str(dest.relative_to(base))]=value
    for name,value in cp['sourceFiles'].items():
        f.require(f.digest(origin/'source'/name)==f.digest(f.ROOT/name)==value,'focused/current source identity')
    for name,record in cp['artifacts'].items():
        f.require(f.digest(origin/name)==record['sha256'],'focused artifact identity')
        dest=base/('accepted-focused/'+name if name=='remaining-case-currents.pickle' else name)
        retain(origin/name,dest)
    for name,value in cp['acceptanceArtifacts'].items():
        f.require(f.digest(proofroot/name)==value,'accepted semantic proof artifact')
        retain(proofroot/name,base/'accepted-reuse-proofs'/name)
    groups={};proofs={}
    for label,ends in sources['currentInputs'].items():
        for end,data in ends.items():
            address=label+'__'+end;path=base/'accepted-reuse-proofs'/(address+'.pickle')
            packet=f.unpickle(path);old=accepted['RIGHT' if end=='RIGHT' else 'LEFT']
            f.require(same(packet['profileBindings'],data['profileBindings']) and same(packet['fieldUnits'],data['fieldUnits']),'actual profile and field-unit proof inputs')
            reuse=check_pair_inputs(packet['pairs'],signature(data),accepted_signature(old))
            f.require(reuse==cp['caseInputComparisons'][address]['originalFullInputReusable'],'actual original-input disposition')
            groups[address]=('accepted-original','RIGHT' if end=='RIGHT' else 'LEFT') if reuse else ('new-input',address)
            proofs[address]=packet
            retain(origin/'cases'/label/end.lower()/'balance-reuse.json',base/'cases'/label/end.lower()/'balance-reuse.json')
    lab=sources['currentInputs']['LAB_HELD__RHOBR_CONSTANT']['RIGHT']
    mat=sources['currentInputs']['MATERIAL_ADVECTED__RHOBR_CONSTANT']['RIGHT']
    shared=f.unpickle(base/'accepted-reuse-proofs/right-rhobr-cross-anchoring.pickle')
    f.require(check_pair_inputs(shared,signature(mat),signature(lab)),'complete new-family cross-anchoring input join')
    f.require(same(lab['profileBindings'],mat['profileBindings']) and same(lab['fieldUnits'],mat['fieldUnits']),'new family actual profiles and field units')
    a='LAB_HELD__RHOBR_CONSTANT__RIGHT';b='MATERIAL_ADVECTED__RHOBR_CONSTANT__RIGHT'
    f.require(groups[a][0]==groups[b][0]=='new-input','actual new acoustic inputs retained')
    groups[b]=groups[a]
    f.require(sum(group[0]=='accepted-original' for group in groups.values())==10,'ten complete original current input matches')
    f.require(len(set(groups.values()))==3,'actual two original and one new acoustic source families')
    for p in (Path(__file__).resolve(),PLAN,CHECKPOINT):
        name=str(p.relative_to(f.ROOT));manifest['sourceFiles'][name]=f.digest(p)
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    inputs[str(origin/'checks.json')]=cp['checksSha256']
    inputs[cp['validationScript']]=cp['validationScriptSha256']
    manifest['inputPackets'].update(inputs);manifest['copiedInputs'].update(copies)
    manifest['focusedReuse']={'runDirectory':str(origin),'checksSha256':cp['checksSha256'],'copiedArtifacts':copies,
        'acousticInputGroups':groups,'wholeCoordinatorReverseAstJoin':coordinator()[1],
        'scope':'Copy completed tails and exact actual input certificates; compute only the one new acoustic/pairing family.'}
    f.save(base/'inputs.json',manifest)
    STATE.update(base=base,sources=sources,accepted=accepted,proofs=proofs,groups=groups,manifest=manifest)
    return sources,accepted,paired,manifest


def joined_address(data):
    address=data['case']+'__'+data['end']
    f.require(address in STATE['groups'],'known actual case/end address')
    original=STATE['sources']['currentInputs'][data['case']][data['end']]
    f.require(same(signature(data),signature(original)) and same(data['profileBindings'],original['profileBindings']) and same(data['fieldUnits'],original['fieldUnits']),'live full acoustic signature and address')
    return address


def original_current(data,old):
    address=joined_address(data)
    reuse=check_pair_inputs(STATE['proofs'][address]['pairs'],signature(data),accepted_signature(old))
    f.require(reuse==(STATE['groups'][address][0]=='accepted-original'),'validated original current routing')
    return reuse


def full_signature(data):
    return {'address':joined_address(data),'inputs':signature(data),'profileBindings':data['profileBindings'],'fieldUnits':data['fieldUnits']}


def matched_families(data,families):
    address=joined_address(data);matches=[]
    for index,(saved,_) in enumerate(families):
        owner=saved['address'];label,end=owner.rsplit('__',1)
        f.require(same(saved,full_signature(STATE['sources']['currentInputs'][label][end])),'previous complete family input identity')
        if STATE['groups'][owner]==STATE['groups'][address]:matches.append(index)
    f.require(len(matches)<=1,'one completed owner for actual current family')
    return matches


def saved_tail(base,label,end,current,data,old):
    f.require(joined_address(data)==label+'__'+end,'saved conservative tail context')
    result=f.unpickle(base/'conservative.pickle');proof=f.unpickle(base/'tail-reuse.pickle')
    f.require(proof['case']==label and proof['end']==end,'actual tail address')
    for key,pair in proof['pairs'].items():
        f.require(same(data[key],pair[0]) and same(old['conservative'][key],pair[1]),'actual saved energy/constraint pair')
    f.require(all(v==0 for v in h.ends_source.scalars(proof['residuals'])),'completed conservative tail proof')
    for key in proof['inputPayloadChanges']:f.require(same(result[key],data[key]),'preserved actual raw current input')
    f.require(same(result['MATERIAL_CONSTRAINT_RETAINED_RESIDUAL'],data['retainedConstraintResidual']),'actual retained constraint operand')
    for key in proof['unchangedTailKeys']:f.require(same(result[key],old['conservative'][key]),'preserved complete native tail')
    d=engine.PHYSICAL_METADATA.dimensions
    atoms=(*current.fields[0],*current.fields[1],*current.variations[0],*current.variations[1],*current.amplitudes[0],*current.amplitudes[1],*current.leg_momenta,current.phase_coordinate,current.virtual_parameter)
    f.require(all(atom in old['knownDimensions'] and d.known[atom]==old['knownDimensions'][atom] for atom in atoms),'actual current field/coordinate units')
    return result,proof


def retained_slab(path,slab):
    f.require(same(f.unpickle(path),slab),'complete original slab balance copy')


def coordinator():
    source=Path(h.__file__).read_text();node=next(n for n in ast.parse(source).body if getattr(n,'name',None)=='main')
    replacements=[
        ("[i for i,(signature,_) in enumerate(families) if same(signature['strong'],data['strong']) and same(signature['energy'],data['TANGENTIAL_ENERGY_REDUCTION']) and same(signature['constraint'],data['MATERIAL_CONSTRAINT_COEFFICIENTS'])]",'matched_families(data,families)'),
        ("same(data['strong'],old['strong'])",'original_current(data,old)'),
        ("{'strong':data['strong'],'energy':data['TANGENTIAL_ENERGY_REDUCTION'],'constraint':data['MATERIAL_CONSTRAINT_COEFFICIENTS']}",'full_signature(data)'),
        ("f.atomic_pickle(target/'slab.pickle',slab)","retained_slab(target/'slab.pickle',slab)")]
    pairs=[(ast.parse(a,mode='eval').body,ast.parse(b,mode='eval').body) for a,b in replacements]
    def convert(tree,pairs):
        hits=[0]*len(pairs)
        class Rewrite(ast.NodeTransformer):
            def visit(self,n):
                for i,(a,b) in enumerate(pairs):
                    if ast.dump(n)==ast.dump(a):hits[i]+=1;return copy.deepcopy(b)
                return super().visit(n)
        result=Rewrite().visit(copy.deepcopy(tree));f.require(hits==[1]*len(pairs),'four exact coordinator routing edits');return result
    generated=convert(node,pairs);restored=convert(generated,[(b,a) for a,b in pairs])
    f.require(ast.dump(restored)==ast.dump(node),'whole coordinator reverse AST join')
    namespace=dict(vars(h),load=prepare,join_tail=saved_tail,matched_families=matched_families,
        original_current=original_current,full_signature=full_signature,retained_slab=retained_slab)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[generated],type_ignores=[])),__file__,'exec'),namespace)
    fn=namespace['main'];f.require(fn.__globals__['load'] is prepare and fn.__globals__['join_tail'] is saved_tail,'actual completed-operand dispatch')
    return fn,{'wholeFunctionJoined':True,'expressionEdits':len(pairs),'namespaceOverrides':['load','join_tail'],'nativeHelperSha256':f.digest(Path(h.__file__))}


def focused(directory):
    directory=directory.resolve();directory.relative_to(f.STORE);directory.mkdir(exist_ok=False,parents=True)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(180);start=time.monotonic()
    sources,accepted,paired,manifest=prepare(directory);families=[];assignments={};controls=[]
    for label,ends in sources['currentInputs'].items():
        for end,data in ends.items():
            matches=matched_families(data,families);old=accepted['RIGHT' if end=='RIGHT' else 'LEFT']
            original=original_current(data,old)
            if not matches:families.append((full_signature(data),{}))
            assignments[label+'__'+end]={'family':matches[0] if matches else len(families)-1,'originalReusable':original}
    data=sources['currentInputs']['MATERIAL_ADVECTED__RHO4_CONSTANT']['RIGHT']
    for key in ('strong','ZERO_TRANSFER_MASS_ROW','TANGENTIAL_ENERGY_REDUCTION','MATERIAL_CONSTRAINT_COEFFICIENTS','profileBindings','case'):
        changed=dict(data)
        if key=='case':changed[key]='WRONG'
        elif key=='profileBindings':changed[key]={**data[key],h.sp.Symbol('wrongProfileAddress'):h.sp.Integer(1)}
        elif isinstance(data[key],tuple):changed[key]=tuple(2*v for v in data[key])
        else:changed[key]=2*data[key]
        try:joined_address(changed)
        except ValueError:controls.append(key)
        else:raise ValueError(('changed source/address accepted',key))
    f.require(len(families)==3 and len(controls)==6,'actual family/address mutation checks')
    for name,value in manifest['inputPackets'].items():f.require(f.digest(Path(name))==value,'unchanged focused input')
    for name,value in manifest['copiedInputs'].items():f.require(f.digest(directory/name)==value,'unchanged focused copy')
    checks={'status':'ACCEPTED_SAVED_CURRENT_REUSE_ROUTING','assignments':assignments,'mutationControls':controls,'coordinatorJoin':coordinator()[1],
        'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'copiedInputs':manifest['copiedInputs'],'wallSeconds':time.monotonic()-start,
        'scope':'Saved-input and dispatch checks only; no new current, pairing, mode, quadrature or solve.'}
    f.save(directory/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':
    if '--check-reuse' in sys.argv:
        import argparse
        p=argparse.ArgumentParser();p.add_argument('--check-reuse',action='store_true');p.add_argument('--run-directory',type=Path,required=True)
        args=p.parse_args();focused(args.run_directory)
    else:
        f.require('--focused' not in sys.argv,'production requires completed focused reuse evidence')
        coordinator()[0]()
        checks=json.loads((STATE['base']/'checks.json').read_text())
        f.require(checks['newAcousticConstructions']==checks['newTwoFrequencyConstructions']==1,'only the one genuinely new current/pairing family computed')
