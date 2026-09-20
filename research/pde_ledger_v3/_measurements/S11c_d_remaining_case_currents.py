#!/usr/bin/env python3
"""Actual case current closures, reusing only joined native generator inputs."""
import argparse,ast,copy,json,resource,shutil,signal,time
from pathlib import Path
import sympy as sp
import S11c_d_remaining_case_end_sources as ends_source
import S11c_d_end_current_resumable as rational
f,engine=ends_source.f,ends_source.engine
same=ends_source.binding.same
CHECKPOINT=f.M/'S11c_d_remaining_case_end_sources_checkpoint.json'
PLAN=f.M/'S11c_d_remaining_case_currents_plan.md'


def source_ast(path,name):
    return ast.dump(next(n for n in ast.parse(path.read_text()).body if getattr(n,'name',None)==name))


def unchanged_cache(value,anchoring,end):
    def cached(a,e):
        f.require((a,e)==(anchoring,end),'actual cached current context')
        return value
    return cached


def zero_difference(a,b):
    value=ends_source.expanded_difference(a,b)
    f.require(all(v==0 for v in ends_source.scalars(value)),'exact consumed current operand')
    return value


def load(base):
    cp=json.loads(CHECKPOINT.read_text());origin=Path(cp['runDirectory'])
    f.require(cp['status']=='ACCEPTED_CASE_END_AND_CURRENT_INPUT_SOURCES','accepted case end/current inputs')
    f.require(f.digest(origin/'checks.json')==cp['checksSha256'],'completed end-source checkpoint')
    for n,v in cp['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(origin/'source'/n)==v,('unchanged accepted source',n))
    inputs={str(origin/'checks.json'):cp['checksSha256']};copied={}
    names=['remaining-case-end-sources.pickle','accepted-left-current-source.pickle','accepted-right-current-source.pickle']
    for label in {k.rsplit('__',1)[0] for k in cp['familyAssignment']}:
        names.extend('accepted-cases/'+label+'/'+n+'.pickle' for n in ('reduced-action','actions','assembly'))
    for n in sorted(names):
        p=origin/n;f.require(f.digest(p)==cp['artifacts'][n]['sha256'],'accepted current input packet')
        dest=base/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest);inputs[str(p)]=f.digest(p);copied[n]=f.digest(dest)
    accepted={end:f.unpickle(base/('accepted-'+end.lower()+'-current-source.pickle')) for end in ('LEFT','RIGHT')}
    joins={};checkpoints=[];paired={}
    for end in ('left','right'):
        path=f.M/('S11c_d_end_current_source_'+end+'_thickness_repair_checkpoint.json');old=json.loads(path.read_text());directory=Path(old['runDirectory'])
        frozen=directory/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
        f.require(f.digest(directory/'objects.pickle')==copied['accepted-'+end+'-current-source.pickle'],'original full current source identity')
        joins[end]={n:source_ast(frozen,n)==source_ast(engine.HERE,n) for n in ('UniformSlabCurrent','SlabEnergyBalance','ClosedAcousticEnergy','ConstantEndPencil','polynomial_terms')}
        f.require(all(joins[end].values()),'unchanged whole native current constructors')
        inputs[str(frozen)]=f.digest(frozen);checkpoints.append(path)
        pairing_cp=f.M/('S11c_d_end_pairing_'+end+'_thickness_repair_checkpoint.json');pairing=json.loads(pairing_cp.read_text());directory=Path(pairing['runDirectory']);frozen=directory/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
        joins[end]['ClosedCurrentPairing']=source_ast(frozen,'ClosedCurrentPairing')==source_ast(engine.HERE,'ClosedCurrentPairing')
        f.require(joins[end]['ClosedCurrentPairing'] and pairing['retainedNonzeroScalars']==0,'unchanged accepted two-frequency current constructor')
        f.require(pairing['provenance']['sourcePacketSha256']==copied['accepted-'+end+'-current-source.pickle'],'actual pairing/current source join')
        src=directory/'complete.pickle';f.require(f.digest(src)==pairing['artifacts']['complete.pickle']['sha256'],'accepted pairing packet')
        dst=base/('accepted-'+end+'-pairing.pickle');shutil.copyfile(src,dst);inputs[str(src)]=f.digest(src);copied[dst.name]=f.digest(dst)
        value,known=f.unpickle(dst);f.require(all(v==0 for v in ends_source.scalars(value['retained'])),'accepted full pairing residuals')
        paired[end.upper()]={'packet':value,'knownDimensions':known}
        inputs[str(frozen)]=f.digest(frozen);checkpoints.append(pairing_cp)
    pins=dict(cp['sourceFiles'])
    for p in (Path(__file__).resolve(),PLAN,CHECKPOINT,Path(rational.__file__),*checkpoints):pins[str(p.relative_to(f.ROOT))]=f.digest(p)
    for n,v in pins.items():
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target);f.require(f.digest(target)==v,'frozen current helper')
    manifest={'sourceFiles':pins,'inputPackets':inputs,'copiedInputs':copied,'nativeConstructorJoins':joins,
        'input':json.loads((origin/'inputs.json').read_text())['input'],'sourceCheckpoint':str(CHECKPOINT),
        'scope':'Full native end-current closure and two-frequency matrices; no quadrature, root isolation, mode normalization or scattering solve.'}
    f.save(base/'inputs.json',manifest)
    return f.unpickle(base/'remaining-case-end-sources.pickle'),accepted,paired,manifest


def context(base,label,end,sources,manifest):
    packets={n:f.unpickle(base/'accepted-cases'/label/(n+'.pickle')) for n in ('reduced-action','actions','assembly')}
    r,dims,pencil=ends_source.binding.factors.context(packets['reduced-action'],packets['actions'],packets['assembly'])
    data=sources['currentInputs'][label][end];inp=engine.ChannelInput(r,manifest['input'])
    state=copy.copy(r);state.end_values={key:inp.limits[value] for key,value in r.end_values.items()}
    e=engine.ConstantEndPencil.__new__(engine.ConstantEndPencil);e.r=state;e.kn=sp.Symbol('s11cdSpectralNormalMomentum',real=True)
    dims.known[e.kn]=dims.measure(state.normal_map[state.momentum_groups[0][2]])
    result=sources['sources'][label];m=engine.FullPencilModes(e,result['curl'],result['units']['weak'])
    c=engine.UniformSlabCurrent(state,{'value':data['sourceEnergy']},e,data['strong'][3,:]);balance=engine.SlabEnergyBalance(c);acoustic=engine.ClosedAcousticEnergy(balance,m,data['strong'])
    f.require(same(data['profileBindings'],inp.limits) and tuple(c.field_units)==tuple(data['fieldUnits']),'actual profile and field-unit inputs')
    return c,balance,acoustic,inp,data


def join_tail(base,label,end,current,data,accepted):
    # The entire original suffix from virtual variation through polarization
    # consumes the reduced energy and normalized constraint, with native field,
    # phase, derivative and independent-grade coordinates held literally.
    old=accepted['conservative'];pairs={n:(data[n],old[n]) for n in ('TANGENTIAL_ENERGY_REDUCTION','MATERIAL_CONSTRAINT_COEFFICIENTS')}
    f.atomic_pickle(base/'tail-input-pairs.pickle',pairs)
    residuals={n:zero_difference(a,b) for n,(a,b) in pairs.items()}
    known=accepted['knownDimensions'];d=engine.PHYSICAL_METADATA.dimensions
    for atom in (*current.fields[0],*current.fields[1],*current.variations[0],*current.variations[1],*current.amplitudes[0],*current.amplitudes[1],*current.leg_momenta,current.phase_coordinate,current.virtual_parameter):
        f.require(atom in known and d.known[atom]==known[atom],'actual native current coordinate and unit')
    old_omega=next(s for s in accepted['strong'].free_symbols if s.name==current.r.omega.name)
    f.require(old_omega==current.r.omega,'literal harmonic frequency')
    names={s.name:s for s in old['TANGENTIAL_ENERGY_REDUCTION'].free_symbols}
    for s in (*current.r.tangents,current.r.z,current.r.t,*current.r.symbols.values()):
        if s.name in names:f.require(s==names[s.name],'literal current parameter assumptions')
    result=dict(old)
    input_keys=('SOURCE_ENERGY','UNIFORM_SOURCE_ENERGY','SOURCE_PARAMETER_ALIGNMENT','ZERO_TRANSFER_MASS_ROW','MATERIAL_CONSTRAINT_COEFFICIENTS','MATERIAL_CONSTRAINT_TRUNCATION_REMAINDER','TANGENTIAL_ENERGY_REDUCTION')
    for key in input_keys:result[key]=data[key]
    result['MATERIAL_CONSTRAINT_RETAINED_RESIDUAL']=data['retainedConstraintResidual']
    changed={key:not same(result[key],old[key]) for key in input_keys}
    tail_keys=set(old)-set(input_keys)-{'MATERIAL_CONSTRAINT_RETAINED_RESIDUAL'}
    f.require(all(same(result[k],old[k]) for k in tail_keys),'byte-preserved symbolic conservative tail operands')
    proof={'pairs':pairs,'residuals':residuals,'inputPayloadChanges':changed,'unchangedTailKeys':sorted(tail_keys),'case':label,'end':end}
    f.atomic_pickle(base/'tail-reuse.pickle',proof);f.atomic_pickle(base/'conservative.pickle',result)
    return result,proof


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--focused',action='store_true');args=p.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));start=time.monotonic()
    def timeout(*_):raise TimeoutError('case-current budget; preserve all source and rational/current packets')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    def progress(name):
        with (base/'progress.jsonl').open('a') as stream:stream.write(json.dumps({'stage':name,'wallSeconds':time.monotonic()-start})+'\n')
    sources,accepted,paired,manifest=load(base);outputs={};families=[];inventory={};rational_records={};new_acoustic=0;new_pairing=0;controls=[]
    labels=list(sources['sources'])
    for label in labels:
        outputs[label]={}
        for end in ('REFERENCE','LEFT','RIGHT'):
            target=base/'cases'/label/end.lower();target.mkdir(parents=True,exist_ok=True);progress(label+'_'+end+'_started')
            c,balance,acoustic,inp,data=context(base,label,end,sources,manifest);a=label.split('__')[0];e={'REFERENCE':None,'LEFT':-sp.oo,'RIGHT':sp.oo}[end]
            old=accepted['RIGHT' if end=='RIGHT' else 'LEFT']
            old_pair=paired['RIGHT' if end=='RIGHT' else 'LEFT']
            conservative,reuse=join_tail(target,label,end,c,data,old)
            c.construct=unchanged_cache(conservative,a,e)
            slab=old['slab'];balance.construct=unchanged_cache(slab,a,e);f.atomic_pickle(target/'slab.pickle',slab)
            # The native balance body uses only these actual conservative keys.
            used=('TANGENTIAL_ENERGY_REDUCTION','MATERIAL_CONSTRAINT_COEFFICIENTS','VIRTUAL_VARIATION','NORMAL_BOUNDARY_WORK','SLAB_CURRENT')
            for key in used:zero_difference(conservative[key],old['conservative'][key])
            f.save(target/'balance-reuse.json',{'consumedKeys':used,'wholeNativeConstructorJoined':True})
            matches=[i for i,(signature,_) in enumerate(families) if same(signature['strong'],data['strong']) and same(signature['energy'],data['TANGENTIAL_ENERGY_REDUCTION']) and same(signature['constraint'],data['MATERIAL_CONSTRAINT_COEFFICIENTS'])]
            original_same=same(data['strong'],old['strong'])
            if args.focused:
                if end=='RIGHT' and label=='LAB_HELD__RHOBR_CONSTANT':
                    energy=data['TANGENTIAL_ENERGY_REDUCTION'];f.atomic_pickle(target/'energy-mutation.pickle',(energy,2*energy,old['conservative']['TANGENTIAL_ENERGY_REDUCTION']))
                    try:zero_difference(2*energy,old['conservative']['TANGENTIAL_ENERGY_REDUCTION'])
                    except ValueError:controls.append('actual_same_unit_energy_coefficient_rejected')
                    else:raise ValueError('energy mutation did not reject')
                    try:c.construct(a+'_WRONG',e)
                    except ValueError:controls.append('wrong_current_cache_address_rejected')
                    else:raise ValueError('wrong current cache address did not reject')
                inventory[label+'__'+end]={'originalAcousticReusable':original_same,'conservativeTailReused':True,'balanceReused':True,'prefix':reuse['inputPayloadChanges']};continue
            if matches:
                family=matches[0];bulk= families[family][1]['acoustic'];pair=families[family][1]['pairing'];pair_checks=families[family][1]['checks'];origin=families[family][1]['origin'];mode='same-full-input-family'
            else:
                if original_same:
                    bulk=old['acoustic'];pair=old_pair['packet']['result'];pair_checks=old_pair['packet']['retained'];mode='accepted-original-current-source-and-pairing'
                else:
                    def checkpoint(name,operand,calculate):
                        path=target/'rational'/name;path.parent.mkdir(parents=True,exist_ok=True)
                        f.atomic_pickle(path.with_suffix('.input.pickle'),operand)
                        previous=old['cancellationProof'].get(name)
                        inherited=previous is not None and same(operand,previous['BEFORE'])
                        value=previous['AFTER'] if inherited else calculate()
                        f.atomic_pickle(path.with_suffix('.value.pickle'),value)
                        residual,fractions=rational.polynomial_identity(operand,value)
                        proof={'before':operand,'after':value,'residual':residual,'fractions':fractions,'originalSourceDenominator':operand.as_numer_denom()[1],'originalResultDenominator':value.as_numer_denom()[1],'acceptedOperationReused':inherited}
                        f.atomic_pickle(path.with_suffix('.proof.pickle'),proof);f.require(residual==0,'actual native rational operation proof');rational_records[str(path.relative_to(base))]={'acceptedOperationReused':inherited,'proofSha256':f.digest(path.with_suffix('.proof.pickle'))};f.save(base/'rational-inventory.json',rational_records)
                        return value
                    acoustic.rational_checkpoint=checkpoint;bulk=acoustic.construct(a,e);new_acoustic+=1;mode='new-actual-current-source'
                f.atomic_pickle(target/'acoustic.pickle',bulk);progress(label+'_'+end+'_acoustic_saved')
                acoustic.construct=unchanged_cache(bulk,a,e)
                if not original_same:
                    pairing=engine.ClosedCurrentPairing(acoustic,a,e)
                    pair=pairing.construct(a,e);new_pairing+=1;f.atomic_pickle(target/'pairing.pickle',pair);progress(label+'_'+end+'_pairing_saved')
                    algebraic,relation,_=acoustic.modes.analytic(data['strong'])
                    live=(acoustic.modes.k,acoustic.modes.q,acoustic.modes.eta,acoustic.modes.sigma,c.r.omega)
                    mapping=inp.mapping(algebraic,relation,live)
                    pair_checks=pairing.split_balance_checks(pair,mapping,lambda v:progress(label+'_'+end+'_'+v['stage']))
                    f.atomic_pickle(target/'pairing-raw-checks.pickle',pair_checks)
                    pair_checks={key:tuple(c.retained(v) for v in ends_source.scalars(value)) for key,value in pair_checks.items() if key.endswith('_RESIDUAL')}
                    f.atomic_pickle(target/'pairing-retained-checks.pickle',pair_checks)
                    f.require(all(v==0 for v in ends_source.scalars(pair_checks)),'actual retained two-frequency balance checks')
                family=len(families);origin=(label,end)
                families.append(({'strong':data['strong'],'energy':data['TANGENTIAL_ENERGY_REDUCTION'],'constraint':data['MATERIAL_CONSTRAINT_COEFFICIENTS']},{'acoustic':bulk,'pairing':pair,'checks':pair_checks,'origin':origin}))
            f.atomic_pickle(target/'acoustic.pickle',bulk);f.atomic_pickle(target/'pairing.pickle',pair)
            residuals={k:v for k,v in bulk.items() if k.endswith('_RESIDUAL')}
            retained={k:tuple(c.retained(v) for v in ends_source.scalars(value)) for k,value in residuals.items()}
            f.atomic_pickle(target/'acoustic-residuals.pickle',{'raw':residuals,'retained':retained})
            f.require(all(v==0 for v in ends_source.scalars(retained)),'actual retained acoustic source residuals')
            f.require(all(v==0 for v in pair['SOURCE_BRANCH_JOINS']),'actual pairing branches')
            f.require(pair['SLAB_CURRENT_MATRIX'].shape==(5,5) and pair['BULK_NORMAL_CURRENT_DENSITY_MATRIX'].shape==(5,5),'complete physical current arrays')
            for atom,unit in old_pair['knownDimensions'].items():
                d=engine.PHYSICAL_METADATA.dimensions
                if atom in d.known:f.require(d.known[atom]==unit,'inherited current unit identity')
                else:d.known[atom]=unit
            result={'conservative':conservative,'slab':slab,'acoustic':bulk,'pairing':pair,'pairingChecks':pair_checks,'strong':data['strong'],'profileBindings':inp.limits,'knownDimensions':dict(engine.PHYSICAL_METADATA.dimensions.known),'case':label,'end':end,'family':family,'origin':origin,'reuseMode':mode}
            f.atomic_pickle(target/'current-sources.pickle',result);outputs[label][end]=result
            inventory[label+'__'+end]={'family':family,'origin':origin,'reuseMode':mode,'retainedResidualScalars':sum(1 for _ in ends_source.scalars(retained))};f.save(base/'case-inventory.json',inventory)
    f.save(base/'case-inventory.json',inventory)
    if args.focused:f.require(len(controls)==2,'actual current reuse mutation controls')
    f.atomic_pickle(base/'remaining-case-currents.pickle',{'cases':outputs,'inventory':inventory,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets']})
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'current/frozen source identity')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original current input identity')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'copied current input identity')
    artifacts={str(path.relative_to(base)):ends_source.binding.factors.artifact(path) for path in base.rglob('*.pickle') if 'source' not in path.relative_to(base).parts}
    checks={**manifest,'status':'FOCUSED_CURRENT_INPUT_REUSE_CHECKED' if args.focused else 'COMPLETED_CASE_CURRENT_SOURCES','runDirectory':str(base),'caseInventory':inventory,'mutationControls':controls,'newAcousticConstructions':new_acoustic,'newTwoFrequencyConstructions':new_pairing,'rationalOperations':len(rational_records),'artifacts':artifacts,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
