#!/usr/bin/env python3
"""Serial validated regeneration with explicit per-producer local commits.

Continue the already-running b producer through b/c1/c2/d publication. The
controller performs no physics construction: every object/check is produced by
the listed native or focused instrument. An unsuccessful command stops the
queue with its logs intact. It never pushes, edits a source, or reruns b.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
PREFIX='S11c_thickness_coordinate_'
M='_measurements/'
STATUS=ROOT/M/(PREFIX+'regeneration_status.md')
CONTRACT='## Retained user-approved solver/export contract'
CONTRACT_SHA='f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def plan(base):
    baseline=base/'baseline'
    py=sys.executable
    def native(stage):
        command=[py,M+PREFIX+'run_stage.py','--run-directory',str(base/(stage+'_full')),stage]
        if stage=='d':command+=['--','--case','ALL','--dev-symbol-cache',str(base/'d_full/symbols'),
            '--channel-input-file',M+'S11c_d_channel_preflight_input.json','--channel-input-scope','spectrum']
        return command
    def inventory(stage,publish=False):
        command=[py,M+PREFIX+'stage_inventory.py',stage,'--run-directory',str(base/(stage+'_full')),
                 '--baseline',str(baseline)]
        if stage in ('b','c2') and publish:command+=['--physical-check',str(base/(stage+'_checks/checks.json'))]
        return command+(['--publish'] if publish else [])
    def export_check(stage):
        command=[py,M+PREFIX+stage+'_export_check.py','--baseline',str(baseline),
                 '--run-directory',str(base/(stage+'_checks'))]
        if stage=='c2':command+=['--b-checkpoint',M+PREFIX+'b_export_checkpoint.json']
        return command
    def export_validate(stage):return [py,M+PREFIX+'b_export_validate.py','--stage',stage,
        '--run-directory',str(base/(stage+'_checks')),'--transcript',str(base/(stage+'_checks.out')),'--publish']
    rows=[]
    for stage,title in (('b','Regenerate the kinetic-coordinate-correct b primaries'),
                        ('c1','Regenerate c1 against the corrected b export'),
                        ('c2','Regenerate c2 and account for its imported kinetic delta')):
        steps=[] if stage=='b' else [(stage+'_producer',native(stage),None)]
        steps.append((stage+'_inventory',inventory(stage),None))
        if stage in ('b','c2'):
            steps.extend([(stage+'_export_check',export_check(stage),str(base/(stage+'_checks.out'))),
                          (stage+'_export_validate',export_validate(stage),None)])
        steps.append((stage+'_publish',inventory(stage,True),None))
        script={'b':'brane_operator','c1':'bulk_closure','c2':'selfenergy_fold'}[stage]
        ordinary=['scripts/S11c_'+stage+'_exports.py',M+PREFIX+stage+'_stage_inventory.json']
        outputs=['scripts/out/S11c_'+stage+'_'+script+'_sympy_audit.out']
        if stage in ('b','c2'):
            ordinary+=[M+PREFIX+stage+'_export_checkpoint.json']
            outputs+=['scripts/out/'+PREFIX+stage+'_export.out']
        if stage=='c2':ordinary+=[M+'S11c_c2_sympy_guard_evidence.json',M+'S11c_c2_sympy_progress.json']
        rows.append({'stage':stage,'title':title,'steps':steps,'ordinary':ordinary,'outputs':outputs,
                     'evidence':[M+PREFIX+stage+'_stage_inventory.json']+
                                ([M+PREFIX+stage+'_export_checkpoint.json'] if stage in ('b','c2') else [])})
    dplan=json.loads((ROOT/M/(PREFIX+'d_recheck_plan.json')).read_text())
    if Path(dplan['nativeProducerManifest'])!=base/'d_full/manifest.json':raise ValueError('d plan run-directory mismatch')
    rows.append({'stage':'d','title':'Regenerate the full four-case d spectrum and reduced operands after the coordinate repair',
        'steps':[
            ('d_scoped_symbols',[py,M+'S11c_c2_trace_repair_d_ends.py','--run-directory',str(base/'d_ends')],None),
            ('d_producer',native('d'),None),
            ('d_rechecks',[py,M+'S11c_c2_trace_repair_d_recheck.py','--plan',M+PREFIX+'d_recheck_plan.json',
                           '--report',M+PREFIX+'d_rechecks.json'],None),
            ('d_source_join',[py,M+'S11c_c2_trace_repair_d_source_join.py','--run-root',str(base),
                             '--report',M+PREFIX+'d_source_joins.json'],None),
            ('d_publish',[py,M+'S11c_c2_trace_repair_d_publication.py','--prefix',PREFIX+'d_','--publish'],None)],
        'ordinary':[M+PREFIX+'d_'+name+'.json' for name in ('rechecks','source_joins','full_checks','publication')]+
                    [s['result'] for s in dplan['stages'] if s['name']!='codec_expand'],
        'outputs':['scripts/out/S11c_d_mixing_scattering_sympy_audit.out'],
        'evidence':[M+PREFIX+'d_full_checks.json',M+PREFIX+'d_source_joins.json',M+PREFIX+'d_publication.json']})
    return rows


def describe(row):
    if row['stage'] in ('b','c1','c2'):
        inv=json.loads((ROOT/row['evidence'][0]).read_text())
        result={'stage':row['stage'],'exportRows':inv['exportRows'],
            'changedValueSerializations':inv['changedValueSerializations'],
            'addedKeys':inv['addedKeys'],'removedKeys':inv['removedKeys'],
            'nativeOutputBytes':inv['transcriptBytes'],'nativeOutputSha256':inv['transcriptSha256'],
            'exportSha256':inv['exportSha256'],'failures':inv['operationalOrProvenanceFailures']}
        if not inv['published'] or result['failures']:raise ValueError('upstream publication incomplete')
        if len(row['evidence'])>1:
            proof=json.loads((ROOT/row['evidence'][1]).read_text())
            result.update(residualScalars=proof['residualScalars'],nonzeroResidualScalars=proof['nonzeroResidualScalars'])
            if result['nonzeroResidualScalars']:raise ValueError('coordinate residual remains')
        return result
    inv=json.loads((ROOT/row['evidence'][0]).read_text())
    if inv['failures']:raise ValueError('d publication checks incomplete')
    return {key:inv[key] for key in ('nativeWallSeconds','nativePeakRssKiB','nativeArtifact',
                                      'spectrum','jointSheet','outstandingConstructions','failures')}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run-root',type=Path,required=True)
    parser.add_argument('--plan-only',action='store_true')
    args=parser.parse_args();base=args.run_root.resolve();rows=plan(base)
    if args.plan_only:
        for row in rows:
            for _,command,_ in row['steps']:
                if not (ROOT/command[1]).is_file():raise ValueError(('instrument missing',command[1]))
            print(json.dumps(row),flush=True)
        return
    state=base/'continuation';state.mkdir(exist_ok=False)
    plan_file=state/'plan.json';plan_file.write_text(json.dumps(rows,indent=2)+'\n')
    files={str(Path(__file__).relative_to(ROOT)),M+PREFIX+'d_recheck_plan.json'}
    files|={command[1] for row in rows for _,command,_ in row['steps']}
    pins={p:sha(ROOT/p) for p in sorted(files)}
    (state/'source_pins.json').write_text(json.dumps(pins,indent=2)+'\n')
    events=[]
    def progress(event):
        event={**event,'utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
        events.append(event)
        (state/'progress.json').write_text(json.dumps(events,indent=2)+'\n')
        print(json.dumps(event),flush=True)
    def stable():
        if pins!={p:sha(ROOT/p) for p in pins}:raise ValueError('continuation instrument changed')
        report=(ROOT/M/'S11c_d_sympy_builder_report.md').read_text()
        if sha_text(report[report.index(CONTRACT):])!=CONTRACT_SHA:raise ValueError('retained solver/export contract changed')
    progress({'stage':'waiting_for_b','manifest':str(base/'b_full/manifest.json')})
    while True:
        try:producer=json.loads((base/'b_full/manifest.json').read_text())
        except json.JSONDecodeError:
            time.sleep(1);continue
        if 'exit_code' in producer:break
        time.sleep(15)
    if producer['exit_code']!=0:raise ValueError('b producer failed; inspect preserved output')
    completed=[]
    for row in rows:
        stable();progress({'stage':row['stage'],'event':'started'})
        for name,command,transcript in row['steps']:
            stable();out=Path(transcript) if transcript else state/(name+'.stdout');err=state/(name+'.stderr')
            progress({'stage':row['stage'],'operation':name,'event':'started','command':command,
                      'stdout':str(out),'stderr':str(err)})
            with out.open('xb') as stdout,err.open('xb') as stderr:
                process=subprocess.Popen(command,cwd=ROOT,stdout=stdout,stderr=stderr)
                while process.poll() is None:
                    try:process.wait(timeout=45)
                    except subprocess.TimeoutExpired:
                        print(json.dumps({'stage':row['stage'],'operation':name,'event':'running',
                                          'stdoutBytes':stdout.tell(),'stderrBytes':stderr.tell()}),flush=True)
            progress({'stage':row['stage'],'operation':name,'event':'exited','exitCode':process.returncode,
                      'stdoutSha256':sha(out),'stderrSha256':sha(err)})
            if process.returncode or err.stat().st_size:raise ValueError(('stage failed',name,process.returncode))
        result=describe(row);stable()
        evidence={p:sha(ROOT/p) for p in row['evidence']}
        outputs={p:{'sha256':sha(ROOT/p),'bytes':(ROOT/p).stat().st_size} for p in row['outputs']}
        completed.append({'stage':row['stage'],'result':result,'evidence':evidence,'outputs':outputs})
        text='# S11c thickness-coordinate regeneration status\n\n'
        text+='Serial native producers and recorded physical/artifact checks. Each completed producer is committed locally; transcripts use DataLad/git-annex. No push.\n\n'
        for item in completed:
            text+='## Completed '+item['stage']+'\n\n'
            text+='```json\n'+json.dumps(item['result'],indent=2)+'\n```\n\n'
        pending=[v['stage'] for v in rows if v['stage'] not in {x['stage'] for x in completed}]
        text+='Remaining native producers: '+(', '.join(pending) if pending else 'none in this regeneration queue')+'.\n\n'
        text+='Fresh endpoint/reference sources, two-frequency pairing and full current/adjoint normalization remain separate next steps. Native point/path/stratum records do not establish global coverage, scattering or section 3b profile-frequency bound poles.\n'
        STATUS.write_text(text)
        message=state/(row['stage']+'_commit.txt')
        message.write_text('S11c thickness coordinate: '+row['title']+'\n\nWHAT CHANGED\n'+
            'Regenerated the native producer in dependency order against the reference-normalized thickness kinetic action. Validated source/export/output pins and the recorded checks before atomic publication.\n\n'+
            'COMPUTED CHECKPOINT\n'+json.dumps(result,indent=2)+'\n\n'+
            'Evidence digests:\n'+json.dumps(evidence,indent=2)+'\n\n'+
            'STORAGE AND BOUNDARIES\nOrdinary exports/inventories/status are in Git; every .out is saved with DataLad/git-annex. Original annex payloads remain immutable. Upstream deferred controls and supplied premises retain their prior status. This checkpoint does not establish the full S-matrix, global exceptional-domain coverage or a profile-frequency bound-pole set. No S10/Lean/authority edits, reviewer/comparator/Wolfram, downstream run or push.\n')
        ordinary=[ROOT/p for p in row['ordinary']]+[STATUS]
        annex=[ROOT/p for p in row['outputs']]
        if any(not p.is_file() for p in ordinary+annex):raise ValueError('checkpoint artifact missing')
        staged=subprocess.check_output(['git','diff','--cached','--name-only'],cwd=REPO,text=True)
        if staged.strip():raise ValueError('unrelated staged changes require inspection')
        subprocess.run(['git','add','--',*[str(p.relative_to(REPO)) for p in ordinary]],cwd=REPO,check=True)
        subprocess.run(['datalad','save','-F',str(message),'--',*[str(p.relative_to(REPO)) for p in ordinary+annex]],cwd=REPO,check=True)
        for p in annex:
            mode=subprocess.check_output(['git','ls-files','--stage','--',str(p.relative_to(REPO))],cwd=REPO,text=True).split()[0]
            if mode!='120000' or not p.is_symlink():raise ValueError('output was not stored as an annex pointer')
            if sha(p)!=outputs[str(p.relative_to(ROOT))]['sha256']:raise ValueError('annex payload changed')
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
        progress({'stage':row['stage'],'event':'committed','commit':commit,'result':result})
    progress({'stage':'native_regeneration_complete','remaining':'fresh endpoint/reference sources, two-frequency pairing and current/adjoint normalization'})


def sha_text(text):return hashlib.sha256(text.encode()).hexdigest()


if __name__=='__main__':main()
