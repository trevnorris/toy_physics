#!/usr/bin/env python3
"""Resume FORM bindings from copied inputs and the validated baseline view."""
import json,shutil,types
from pathlib import Path
import S11c_d_remaining_case_profile_bindings as h
f,m=h.f,h.modes
OLD=f.STORE/'s11c-remaining-case-profile-20260920/bindings/complete'
FOCUS=f.STORE/'s11c-remaining-case-profile-20260920/settings-repair'
PLAN=f.M/'S11c_d_remaining_case_profile_bindings_recovery_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_profile_bindings_settings_repair.json'


def load(base):
    previous=json.loads((OLD/'inputs.json').read_text());proof=json.loads(REPAIR.read_text())
    helper=Path(h.__file__);name=str(helper.relative_to(f.ROOT))
    f.require(previous['sourceFiles'][name]==proof['oldHelperSha256'] and f.digest(helper)==proof['newHelperSha256'],'exact recorded settings helper transition')
    f.require(proof['wholeFileReverseAstJoin'] and proof['settingsMutationRejections']==16,'actual settings repair acceptance')
    for n,v in previous['sourceFiles'].items():
        f.require(f.digest(OLD/'source'/n)==v,'original frozen profile source')
        if n!=name:f.require(f.digest(f.ROOT/n)==v,'unchanged consumed profile source')
    for n,v in previous['inputPackets'].items():f.require(f.digest(Path(n))==v,'original profile operand')
    for n,v in previous['copiedInputs'].items():f.require(f.digest(OLD/n)==v,'original copied profile operand')
    _,bc,_=h.accepted(h.BCP,'ACCEPTED_BINDINGS_AND_GRADES')
    _,pc,publication=h.accepted(h.PCP,'PUBLISHED_ANNEX_VERIFIED')
    published=f.ROOT/publication['publication']['path']
    f.require(published.is_symlink() and f.digest(published)==publication['publication']['sha256'],'original FORM publication')
    f.require(h.profile.shape_input(bc['input'])==previous['input'] and bc['input']==previous['baselineInput'],'unchanged actual original and altered inputs')
    f.require(f.digest(FOCUS/'checks.json')==proof['checksSha256'] and (FOCUS/'checks.json').read_bytes()==(FOCUS/'stdout').read_bytes() and (FOCUS/'stderr').stat().st_size==0,'clean completed baseline validation')
    for n,v in proof['artifacts'].items():f.require(f.digest(FOCUS/n)==v['sha256'],'completed settings/baseline artifact')
    f.require(m.node(helper,'baseline_view')==m.node(OLD/'source'/name,'baseline_view'),'unchanged complete baseline validator')
    pins=dict(previous['sourceFiles']);pins[name]=proof['newHelperSha256']
    for path in (Path(__file__).resolve(),PLAN,REPAIR):pins[str(path.relative_to(f.ROOT))]=f.digest(path)
    manifest=dict(previous,runDirectory=str(base),sourceFiles=pins,inputPackets=dict(previous['inputPackets']),copiedInputs={})
    for n,v in previous['copiedInputs'].items():m.retain(OLD/n,base/n,manifest,v)
    m.retain(OLD/'inputs.json',base/'original-inputs.json',manifest)
    m.retain(OLD/'source'/name,base/'original-profile-bindings-helper.py',manifest)
    m.retain(FOCUS/'checks.json',base/'baseline-validation.json',manifest,proof['checksSha256'])
    for filename,key in [('baseline-binding-view.pickle','baselineViewSha256'),('baseline-profile-row-matrices.pickle','baselineRowsSha256')]:
        m.retain(FOCUS/'baseline'/filename,base/filename,manifest,proof[key])
    m.retain(FOCUS/'actual-settings.pickle',base/'actual-settings.pickle',manifest)
    m.retain(FOCUS/'baseline-layout-evidence.pickle',base/'baseline-layout-evidence.pickle',manifest)
    manifest['settings']=h.restore_settings(previous['settings'],f.unpickle(base/'accepted-finite-system.pickle')['settings'])
    manifest['completedBaselineValidationReuse']={'directory':str(FOCUS),'checksSha256':proof['checksSha256'],'baselineValidatorAstUnchanged':True,'baselineRows':80,'fieldCoefficients':129,'originalCopies':len(previous['copiedInputs']),'newBindingsBeforeRecovery':0}
    for n,v in pins.items():
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target);f.require(f.digest(target)==v,'current recovery frozen source')
    f.save(base/'inputs.json',manifest)
    return manifest,tuple(bc['cases'])


def baseline_view(base,manifest):
    # Completed original validator and full-array checks are source/hash joined
    # in load. This only restores its two completed results and original inputs.
    view=f.unpickle(base/'baseline-binding-view.pickle')
    f.require(h.native.same(view['settings'],manifest['settings']),'restored settings in validated baseline view')
    return (view,f.unpickle(base/'accepted-bindings'/h.BASELINE/'case-binding.pickle'),
            f.unpickle(base/'accepted-profile/profile-form.pickle'),
            f.unpickle(base/'accepted-profile/ablation/finite-system.pickle'))


def coordinator():
    namespace=dict(vars(h),load=load,baseline_view=baseline_view)
    result=types.FunctionType(h.main.__code__,namespace,h.main.__name__,h.main.__defaults__,h.main.__closure__)
    f.require(result.__code__ is h.main.__code__,'unchanged complete profile binding coordinator bytecode')
    return result


if __name__=='__main__':coordinator()()
