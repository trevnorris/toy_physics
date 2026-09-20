#!/usr/bin/env python3
"""Case-specific profile bindings, with exact reuse of completed baseline FORM."""
import argparse,copy,json,resource,shutil,signal,time
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_remaining_case_bindings as native
import S11c_d_remaining_case_matrices as matrices
import S11c_d_remaining_case_modes as modes
import S11c_d_profile_form as profile
f,engine=native.f,native.engine
BASELINE=native.BASELINE
PLAN=f.M/'S11c_d_remaining_case_profile_bindings_plan.md'
BCP=f.M/'S11c_d_remaining_case_bindings_checkpoint.json'
PCP=f.M/'S11c_d_profile_form_checkpoint.json'


def accepted(path,status):
    cp=json.loads(path.read_text());root=Path(cp['runDirectory']);f.require(cp['status']==status and f.digest(root/'checks.json')==cp['checksSha256'],'actual accepted producer checkpoint')
    checks=json.loads((root/'checks.json').read_text())
    for n,v in checks['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(root/'source'/n)==v,'accepted source/frozen identity')
    for n,v in checks['inputPackets'].items():f.require(f.digest(Path(n))==v,'accepted original input')
    for n,v in checks['artifacts'].items():f.require(f.digest(root/n)==v['sha256'],'accepted completed artifact')
    return root,checks,cp


def load(base):
    origin,checks,bc=accepted(BCP,'ACCEPTED_BINDINGS_AND_GRADES');pr,pc,publication=accepted(PCP,'PUBLISHED_ANNEX_VERIFIED')
    pub=f.ROOT/publication['publication']['path'];f.require(pub.is_symlink() and f.digest(pub)==publication['publication']['sha256'],'actual accepted baseline FORM publication')
    pins=dict(checks['sourceFiles'])
    for n,v in pc['sourceFiles'].items():f.require(n not in pins or pins[n]==v,'same actual FORM source versions');pins[n]=v
    for p in (Path(__file__).resolve(),PLAN,BCP,PCP,Path(profile.__file__).resolve(),Path(matrices.__file__).resolve(),Path(modes.__file__).resolve()):pins[str(p.relative_to(f.ROOT))]=f.digest(p)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':{},'copiedInputs':{},'input':profile.shape_input(checks['input']),'baselineInput':checks['input'],'settings':checks['settings'],
      'scope':'Three new case profile-FORM source bindings; baseline FORM and original grades reused. No quadrature, modes or responses in this preparation.'}
    for root,source in ((origin,checks),(pr,pc)):
        manifest['inputPackets'][str(root/'checks.json')]=f.digest(root/'checks.json')
        for n,v in source['inputPackets'].items():f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n]==v,'shared original FORM inputs');manifest['inputPackets'][n]=v
    for label in checks['cases']:
        modes.retain(origin/'cases'/label/'case-binding.pickle',base/'accepted-bindings'/label/'case-binding.pickle',manifest)
        for name in ('reduced-action','actions','assembly'):
            modes.retain(origin/'accepted-cases'/label/(name+'.pickle'),base/'accepted-cases'/label/(name+'.pickle'),manifest)
        modes.retain(origin/'cases'/label/'factorization.pickle',base/'accepted-cases'/label/'factorization.pickle',manifest)
    for name in ('accepted-domain-binding.pickle','accepted-source-binding.pickle','accepted-finite-system.pickle'):
        modes.retain(origin/name,base/name,manifest)
    # Retain every accepted baseline FORM packet, including full layouts/partials
    # and original response operands. No producer is called again.
    for n,v in pc['artifacts'].items():modes.retain(pr/n,base/'accepted-profile'/n,manifest,v['sha256'])
    modes.retain(pr/'inputs.json',base/'accepted-profile/inputs.json',manifest)
    original_inputs=json.loads((base/'accepted-profile/inputs.json').read_text())
    f.require(manifest['input']==original_inputs['input'] and manifest['baselineInput']==original_inputs['baselineInput'],'same actual baseline and altered shapes')
    for n,v in pins.items():
        dst=base/'source'/n;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dst);f.require(f.digest(dst)==v,'current profile source snapshot')
    f.save(base/'inputs.json',manifest);return manifest,tuple(checks['cases'])


def full_signature(row,data):
    bound=data['bound']
    characters=tuple(bound['sources'][0,a['sourceIndex']]['boundCharacter'] for a in row['factors'])
    return native.signature(row,data),characters,data['fieldUnits'],data['equationUnits'],bound['profileUnits'],bound['abel'],data['settings']


def baseline_view(base,manifest):
    case=f.unpickle(base/'accepted-bindings'/BASELINE/'case-binding.pickle');altered=f.unpickle(base/'accepted-profile/altered-binding.pickle');result=f.unpickle(base/'accepted-profile/profile-form.pickle');system=f.unpickle(base/'accepted-profile/ablation/finite-system.pickle')
    f.require(system['settings']==manifest['settings'],'complete baseline FORM quadrature settings')
    f.require(native.same(result['fieldUnits'],case['binding']['fieldUnits']) and native.same(altered['bound']['equationUnits'],case['binding']['equationUnits']),'actual inherited FORM row and field units')
    view={k:altered[k] for k in ('bound','jets','local')};view.update(fieldUnits=case['binding']['fieldUnits'],equationUnits=case['binding']['equationUnits'],settings=system['settings'],dimensionState=case['grades']['dimensionState'],acceptedProfileBinding=str(base/'accepted-profile/altered-binding.pickle'))
    for old,row in zip(case['binding']['bound']['rows'],view['bound']['rows']):
        f.require(native.same((old['original'],old['symbolicLimits'],old['unit']),(row['original'],row['symbolicLimits'],row['unit'])),'original baseline FORM row/source address')
    rows={}
    for count in (1,2,3):
        group=f.unpickle(base/'accepted-profile/ablation'/f'layout-{count}.pickle')
        f.require(group['setting']==system['settings'] and np.isfinite(group['matrices']).all(),'actual complete baseline FORM layout')
        f.require(abs(group['massResidual'])<1e-9*(1+abs(group['mass'])) and np.max(abs(group['actionResidual']))/(1+np.max(abs(group['direct'])))<1e-10,'accepted actual FORM measure/direct-action guards')
        rows.update(zip(group['rowIndices'],group['matrices']))
    f.require(set(rows)==set(range(len(view['bound']['rows']))) and len(system['nodes'])==129,'complete baseline FORM rows and finite trial basis')
    f.atomic_pickle(base/'baseline-binding-view.pickle',view);f.atomic_pickle(base/'baseline-profile-row-matrices.pickle',rows)
    return view,case,result,system


def grade_join(case,packets,context,binding):
    r,dimensions,pencil,adapter=context;grade=case['grades']
    specs=native.grades.specs(packets['assembly']['result'],packets['factorization']['result'],binding['fieldUnits'],binding['equationUnits'],r)
    count=0
    for key,value,unit,address in specs:
        item=grade['records'][key];f.require(native.same(item['address'],address) and native.same(item['record']['ORIGINAL'],value) and tuple(item['record']['UNIT'])==tuple(unit),'actual unchanged symbolic grade source/address/unit')
        if address[0] in ('factor','source'):f.require(set(item['record']['COEFFICIENTS'])<={(0,0,0)},'actual grade-free numerical FORM operand')
        count+=1
    f.require(count==len(grade['records']),'complete unchanged independent-grade census');return count


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);args=ap.parse_args();start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False);manifest,labels=load(base)
    baseline,basecase,profile_result,system=baseline_view(base,manifest)
    signatures={};all_bindings={BASELINE:baseline}
    for row in baseline['bound']['rows']:signatures.setdefault(hash(row['original']),[]).append((BASELINE,row['index'],full_signature(row,baseline)))
    inventory={BASELINE:{'rows':len(baseline['bound']['rows']),'sources':len(baseline['jets']),'newRows':[],'reusedBaselineProfile':True,'terms':len(basecase['grades']['termJoins'])}};case_paths={};mutations=[]
    old={'domain-binding.pickle':f.unpickle(base/'accepted-domain-binding.pickle')}
    for label in labels:
        if label==BASELINE:continue
        target=base/'cases'/label;target.mkdir(parents=True);case=f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle');packets={n:f.unpickle(base/'accepted-cases'/label/(n+'.pickle')) for n in ('reduced-action','actions','assembly','factorization')}
        bound,context=native.bind_case(target,packets,old,manifest);r,dimensions,pencil,altered=context
        original=engine.NumericalReducedAction(pencil,packets['assembly']['result'],manifest['baselineInput'])
        endpoints={key:(original.input.limits[key],altered.input.limits[key]) for key in original.input.limits}
        f.atomic_pickle(target/'endpoint-pairs.pickle',endpoints);f.require(all(a==c for a,c in endpoints.values()),'actual unchanged endpoints before case end-map reuse')
        f.require(native.same(original.input.profiles,profile_result['baselineShapes']) and native.same(altered.input.profiles,profile_result['alteredShapes']),'actual profiles before original moment reuse')
        grade_count=grade_join(case,packets,context,bound);reused=[];new=[];joins=[]
        for row in bound['bound']['rows']:
            actual=full_signature(row,bound);bucket=signatures.setdefault(hash(row['original']),[]);matches=[(owner,index) for owner,index,value in bucket if native.same(actual,value)]
            if matches:reused.append({'row':row['index'],'fromCase':matches[0][0],'fromRow':matches[0][1]})
            else:new.append(row['index']);bucket.append((label,row['index'],actual))
            joins.append({'row':row['index'],'sources':tuple(v['sourceIndex'] for v in row['factors']),'reused':bool(matches)})
        actual=bound['bound']['rows'][0];changed=copy.copy(actual);limit=actual['limits'][0];changed['limits']=(sp.Tuple(limit[0],limit[1],limit[2]+1),*actual['limits'][1:])
        f.require(not native.same(full_signature(changed,bound),full_signature(actual,bound)),'actual altered physical limit rejects reuse')
        changed=copy.deepcopy(actual);fi=next(i for i,v in enumerate(changed['factors']) if v['coefficient']!=0);changed['factors'][fi]['coefficient']*=2
        f.require(not native.same(full_signature(changed,bound),full_signature(actual,bound)),'actual changed bound source coefficient rejects reuse');mutations.append(label)
        f.require(all(len(bound['bound']['rows'][index]['limits'])<=2 for index in new),'no new baseline triple FORM integration')
        f.require(len(reused)+len(new)==len(bound['bound']['rows']),'complete actual FORM row partition')
        result={'binding':bound,'grades':case['grades'],'reusedRows':reused,'newRows':new,'rowJoins':joins,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'endpoints':endpoints,'momentsReusedFrom':str(base/'accepted-profile/profile-form.pickle')}
        f.atomic_pickle(target/'case-binding.pickle',result);case_paths[label]=str(target/'case-binding.pickle');all_bindings[label]=bound
        inventory[label]={'rows':len(bound['bound']['rows']),'sources':len(bound['jets']),'terms':len(case['grades']['termJoins']),'gradeRecordsReused':grade_count,'gaussianChecks':len(bound['gaussianComparisons']),'profiles':len(bound['bound']['profileUnits']),'newRows':new,'reusedRows':reused,'unchangedEndpointCount':len(endpoints)};f.save(base/'case-inventory.json',inventory)
    f.atomic_pickle(base/'profile-case-bindings.pickle',{'cases':case_paths,'baselineView':str(base/'baseline-binding-view.pickle'),'inventory':inventory,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':manifest['scope']})
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'current/frozen profile source pre/post')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original profile input pre/post')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'copied profile input pre/post')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,'status':'COMPLETED_CASE_PROFILE_BINDINGS','cases':inventory,'casePackets':case_paths,'newUnionRows':sum(len(v['newRows']) for v in inventory.values()),'actualReuseMutations':mutations,'newQuadratureNodes':0,'newModeConstructions':0,'newScatteringSolves':0,'newGradeExtractions':0,'artifacts':artifacts,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
