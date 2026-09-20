#!/usr/bin/env python3
"""Finish aggregation using completed numerical validation and baseline checks."""
import argparse,json,resource,shutil,signal,time
from pathlib import Path
import S11c_d_remaining_case_response_finish as finish
f,h=finish.f,finish.h
PLAN=f.M/'S11c_d_remaining_case_response_output_plan.md'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);ap.add_argument('--resume-validation',type=Path,required=True);args=ap.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    previous=args.resume_validation.resolve();old_manifest=json.loads((previous/'inputs.json').read_text())
    validation=json.loads((previous/'operand-validation.json').read_text());origin=Path(old_manifest['originalDirectory'])
    outcome=json.loads((previous.parent/'response_finish.invocation.json').read_text())
    finish.equal(outcome,json.loads((previous.parent/'active.json').read_text()),'actual previous final outcome')
    f.require(outcome['exitCode']==1 and outcome['status']=='failed' and not (previous/'checks.json').exists(),'previous finish incomplete')
    error=(previous.parent/'response_finish.stderr').read_text()
    f.require('FileNotFoundError' in error and 'LAB_HELD__RHO4_CONSTANT/continuum/checks.json' in error,'actual missing accepted-check address')
    old,manifest=finish.load(base,origin)
    finish.equal(validation['sourceFiles'],manifest['sourceFiles'],'unchanged entire validator and consumed sources')
    finish.equal(validation['originalArtifacts'],manifest['originalArtifacts'],'same complete validated operands')
    for n,v in old_manifest['copiedInputs'].items():f.require(f.digest(previous/n)==v,('previous validated copied operand',n))
    counts=json.loads((base/'case-inventory.json').read_text());labels=(h.BASELINE,*counts)
    f.require(set(validation['cases'])==set(labels) and validation['newSolves']==validation['individualEmissionsRepeated']==0,'four completed numerical validations')
    # The missing baseline checks are the accepted original producer checks,
    # joined to all six actual copied baseline packets before copying the index.
    checkpoint=json.loads(h.RCP.read_text());baseline=Path(checkpoint['runDirectory'])
    f.require(checkpoint['status']=='PUBLISHED_ANNEX_VERIFIED' and f.digest(baseline/'checks.json')==checkpoint['checksSha256'],'actual accepted baseline checks')
    for name,value in checkpoint['artifacts'].items():
        f.require(f.digest(baseline/name)==f.digest(base/'cases'/h.BASELINE/'continuum'/name)==value['sha256'],('complete baseline artifact join',name))
    for source,target in [(baseline/'checks.json',base/'cases'/h.BASELINE/'continuum/checks.json'),
                          (previous/'operand-validation.json',base/'operand-validation.json')]+[(previous/('validation-'+label+'.json'),base/('validation-'+label+'.json')) for label in labels]:
        source_hash=f.digest(source);shutil.copyfile(source,target);f.require(f.digest(target)==source_hash,'byte-identical completed validation')
        manifest['inputPackets'][str(source)]=source_hash;manifest['copiedInputs'][str(target.relative_to(base))]=source_hash
    for label in labels:finish.equal(json.loads((base/('validation-'+label+'.json')).read_text()),validation['cases'][label],'each completed case validation')
    for name in ('inputs.json',):manifest['inputPackets'][str(previous/name)]=f.digest(previous/name)
    for name in ('active.json','response_finish.invocation.json','response_finish.stderr','response_finish.stdout'):
        p=previous.parent/name;manifest['inputPackets'][str(p)]=f.digest(p)
    for p in (Path(__file__).resolve(),PLAN):
        n=str(p.relative_to(f.ROOT));manifest['sourceFiles'][n]=f.digest(p);target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    manifest['completedValidationReuse']={'directory':str(previous),'sha256':f.digest(previous/'operand-validation.json'),'wholeValidatorUnchanged':True,'cases':list(labels),'newSolves':0,'numericalValidationRepeated':False}
    manifest['acceptedBaselineChecks']={'directory':str(baseline),'sha256':checkpoint['checksSha256'],'completeArtifactJoins':len(checkpoint['artifacts'])}
    f.save(base/'inputs.json',manifest);f.save(base/'completed-validation-reuse.json',manifest['completedValidationReuse'])
    combined=finish.aggregate(base,labels);f.save(base/'aggregation-checks.json',combined)
    f.atomic_pickle(base/'remaining-case-response.pickle',{'cases':counts,'baseline':str(base/'cases'/h.BASELINE),'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':manifest['scope']})
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,('final current/frozen source',n))
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,('final original input',n))
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,('final copied artifact',n))
    finish.equal(finish.inventory(origin),manifest['originalArtifacts'],'all original production pre/post hashes')
    for n,v in old_manifest['copiedInputs'].items():f.require(f.digest(previous/n)==v,('previous finish pre/post hash',n))
    artifacts={n:v for n,v in finish.inventory(base).items() if not n.startswith('source/') and n not in ('inputs.json','checks.json')}
    checks={**manifest,'status':'COMPLETED_FOUR_CASE_RESPONSES','cases':counts,'validation':validation['cases'],'anchoring':validation['anchoring'],
            'combinedOutput':combined,'artifacts':artifacts,'newFiniteCases':0,'newContinuumCases':0,'newQuadratureNodes':0,
            'completedNewFiniteCasesReused':3,'completedNewContinuumCasesReused':3,'individualEmissionsRepeated':0,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
