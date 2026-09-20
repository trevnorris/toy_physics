#!/usr/bin/env python3
"""Reuse accepted derivative-control rows and integrate only the new union."""
import argparse, ast, builtins, copy, dis, gc, json, resource, shutil, signal, time, types
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_remaining_case_first_jet_bindings as binding
import S11c_d_remaining_case_profile_response as reusable

f, engine, modes = binding.f, binding.engine, binding.modes
native = reusable.matrices
BASELINE = binding.BASELINE
CP = f.M/'S11c_d_remaining_case_first_jet_bindings_checkpoint.json'
MCP = f.M/'S11c_d_remaining_case_matrices_checkpoint.json'
JCP = f.M/'S11c_d_first_jet_checkpoint.json'
FCP = f.M/'S11c_d_remaining_case_first_jet_matrices_focused.json'
PLAN = f.M/'S11c_d_remaining_case_first_jet_matrices_plan.md'
SCOPE = ('Actual one-sided first-w-derivative finite interior operators, including changed material-case cell coefficients. '
         'Original and baseline-control matrices are immutable reuse operands. No new modes, currents or response solve.')


def load(base, resume=None):
    br, bc, _ = binding.source.provenance.accepted(CP, 'ACCEPTED_CASE_FIRST_JET_BINDINGS')
    mr, mc, _ = binding.source.provenance.accepted(MCP, 'ACCEPTED_FOUR_CASE_INTERIOR_MATRICES')
    jr, jc, _ = binding.source.provenance.accepted(JCP, 'PUBLISHED_ANNEX_VERIFIED')
    pins = {}
    for checks in (bc, mc, jc):
        for n,v in checks['sourceFiles'].items():
            f.require(n not in pins or pins[n]==v, 'exact consumed source version'); pins[n]=v
    for p in (Path(__file__), PLAN, CP, MCP, JCP, Path(reusable.__file__), Path(native.__file__)):
        pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':{},'copiedInputs':{},
              'input':bc['input'],'settings':bc['settings'],'scope':SCOPE}
    for root,checks in ((br,bc),(mr,mc),(jr,jc)):
        manifest['inputPackets'][str(root/'checks.json')]=f.digest(root/'checks.json')
        for n,v in checks['inputPackets'].items():
            f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n]==v,'same actual original operand'); manifest['inputPackets'][n]=v
    labels=tuple(bc['cases'])
    if resume:
        cp=json.loads(FCP.read_text()); old=json.loads((resume/'checks.json').read_text())
        f.require(cp['status']=='ACCEPTED_FIRST_JET_MATRIX_INPUTS' and f.digest(resume/'checks.json')==cp['checksSha256'],'accepted focused matrix inputs')
        f.require(old['sourceFiles']==pins and old['mode']=='focused','unchanged focused sources')
        for n,v in old['artifacts'].items(): modes.retain(resume/n,base/n,manifest,v['sha256'])
        manifest['inputPackets'][str(resume/'checks.json')]=f.digest(resume/'checks.json')
        manifest['inputPackets'][str(FCP)]=f.digest(FCP)
        manifest['completedFocusedReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'artifacts':len(old['artifacts'])}
    else:
        for n,v in bc['artifacts'].items(): modes.retain(br/n,base/'binding-inputs'/n,manifest,v['sha256'])
        for label in labels:
            for kind in ('reduced-action','actions','assembly','factorization'):
                n='accepted-cases/'+label+'/'+kind+'.pickle'
                modes.retain(br/n,base/n,manifest,bc['artifacts'][n]['sha256'])
            n='cases/'+label+'/row-matrices.pickle'; modes.retain(mr/n,base/'original-matrices'/n,manifest,mc['artifacts'][n]['sha256'])
        for n,v in mc['artifacts'].items():
            if (n.startswith('accepted-layouts/') or '/layout-' in n) and Path(n).name.startswith('layout-'):
                modes.retain(mr/n,base/'original-matrices'/n,manifest,v['sha256'])
        modes.retain(mr/'accepted-finite-system.pickle',base/'accepted-finite-system.pickle',manifest,mc['artifacts']['accepted-finite-system.pickle']['sha256'])
        paths=[n for n in jc['inputPackets'] if n.endswith('/finite-system.pickle')]
        f.require(len(paths)==1,'actual baseline derivative-control trial basis input')
        modes.retain(Path(paths[0]),base/'first-jet-finite-basis.pickle',manifest,jc['inputPackets'][paths[0]])
    for n,v in pins.items():
        dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dest)
        f.require(f.digest(dest)==v,'frozen numerical adapter source')
    system=f.unpickle(base/'accepted-finite-system.pickle')
    manifest['settings']=binding.source.provenance.restore_settings(manifest['settings'],system['settings'])
    f.save(base/'inputs.json',manifest)
    return manifest,labels,system


def cases_and_rows(base, labels):
    cases={}; pools={}
    for label in labels:
        pools['original__'+label]=f.unpickle(base/'original-matrices/cases'/label/'row-matrices.pickle')['rows']
        if label==BASELINE:
            old=f.unpickle(base/'binding-inputs/accepted-bindings'/label/'case-binding.pickle')
            saved=f.unpickle(base/'binding-inputs/accepted-first-jet/first-jet-binding.pickle')
            grade=dict(old['grades'],records=saved['records'],termJoins=saved['termJoins'])
            cases[label]={'binding':f.unpickle(base/'binding-inputs/baseline-first-jet-binding-view.pickle'),
                          'grades':grade,'newRows':[],'reusedRows':[]}
            pools[label]=f.unpickle(base/'binding-inputs/accepted-first-jet/first-jet-rows.pickle')['rows']
        else:
            saved=f.unpickle(base/'binding-inputs/bindings'/label/'case-binding.pickle')
            routes=[dict(v,fromCase=('original__'+v['fromCase'] if v['kind']=='original' else v['fromCase'])) for v in saved['reusedRows']]
            cases[label]=dict(saved,reusedRows=routes)
    return cases,pools


def selected_context(base,label,case,specification):
    r,adapter,packets=native.context(base,label,case,specification)
    if label!=BASELINE:
        selected=f.unpickle(base/'binding-inputs/cases'/label/'term-inputs.pickle')
        f.require(binding.native.same(selected['originalCells'],packets['assembly']['result']['ROWS']) and
                  binding.native.same(selected['originalFourier'],packets['factorization']['result']), 'original independent native cell/factor source join')
        assembly=dict(packets['assembly'],result=dict(packets['assembly']['result'],ROWS=selected['selectedCells']))
        packets=dict(packets,assembly=assembly)
    return r,adapter,packets


def matrix_tail(manifest):
    complete,proof=reusable.matrix_tail()
    scope=dict(complete.__globals__,manifest=manifest,context=selected_context)
    result=types.FunctionType(complete.__code__,scope,complete.__name__,complete.__defaults__,complete.__closure__)
    missing={v.argval for v in dis.get_instructions(result) if v.opname=='LOAD_GLOBAL' and v.argval not in scope and not hasattr(builtins,v.argval)}
    f.require(not missing and result.__code__ is complete.__code__,'whole native matrix-tail bytecode and resolved namespace')
    proof.update(wholeNativeTailBytecodeUnchanged=True,actualManifestIdentity=result.__globals__['manifest'] is manifest,
                 selectedContextSha256=binding.native_body(selected_context),nativeContextSha256=binding.native_body(native.context),
                 directCellsSha256=binding.native_body(native.direct_cells),assembleCaseSha256=binding.native_body(native.assemble_case))
    return result,proof


def preflight(base,manifest,labels,system):
    cases,pools=cases_and_rows(base,labels); same=binding.native.same
    prior=f.unpickle(base/'first-jet-finite-basis.pickle')
    for name in ('nodes','derivativeMatrices','settings'):
        f.require(modes.same(system[name],prior[name]),'actual unchanged collocation and derivative basis')
    f.require(len(system['nodes'])==129,'accepted finite field coefficients')
    evidence=[]
    for directory in (base/'original-matrices',base/'binding-inputs/accepted-first-jet/changed-layouts'):
        for path in directory.rglob('layout-*.pickle'):
            group=f.unpickle(path)
            f.require(group['setting']==system['settings'],'actual reused complete quadrature settings')
            f.require(np.isfinite(group['matrices']).all() and abs(group['massResidual'])<1e-9*(1+abs(group['mass'])),'saved native layout measure evidence')
            f.require(np.max(abs(group['actionResidual']))/(1+np.max(abs(group['direct'])))<1e-10,'saved native direct-action evidence')
            evidence.append({'packet':str(path.relative_to(base)),'rowIndices':group['rowIndices'],'variables':group['variables'],
                             'mass':group['mass'],'massResidual':group['massResidual'],'actionResidual':group['actionResidual']})
    # Join each saved row array to its actual original layout, including the
    # already accepted original-case and baseline-control reuse chains.
    original={label:f.unpickle(base/'binding-inputs/accepted-bindings'/label/'case-binding.pickle') for label in labels}
    coverage={('original__'+BASELINE,i):0 for i in pools['original__'+BASELINE]}
    for path in (base/'original-matrices/accepted-layouts').glob('layout-*.pickle'):
        group=f.unpickle(path)
        for i,array in zip(group['rowIndices'],group['matrices']):
            f.require(np.array_equal(array,pools['original__'+BASELINE][i]),'original baseline actual row/layout identity');coverage['original__'+BASELINE,i]+=1
    for label in labels:
        row_packet=f.unpickle(base/'original-matrices/cases'/label/'row-matrices.pickle')
        if label==BASELINE: continue
        for route in row_packet['reusedRows']:
            f.require(np.array_equal(row_packet['rows'][route['row']],pools['original__'+route['fromCase']][route['fromRow']]),'actual original case reused array chain')
        for group in row_packet['newGroups']:
            for i,array in zip(group['rowIndices'],group['matrices']): f.require(np.array_equal(array,row_packet['rows'][i]),'actual original case new layout/row chain')
    baseline_rows=f.unpickle(base/'binding-inputs/accepted-first-jet/first-jet-rows.pickle')
    covered=set(baseline_rows['reusedRows'])
    for i in baseline_rows['reusedRows']: f.require(np.array_equal(pools[BASELINE][i],pools['original__'+BASELINE][i]),'actual unchanged baseline control row')
    for group in baseline_rows['newGroups']:
        for i,array in zip(group['rowIndices'],group['matrices']):
            f.require(np.array_equal(pools[BASELINE][i],array),'accepted derivative layout/row identity');covered.add(i)
    f.require(covered==set(range(80)) and all(v==1 for v in coverage.values()),'complete baseline original/control matrix provenance')
    f.atomic_pickle(base/'accepted-layout-evidence.pickle',evidence)
    summaries={}; controls=0; new_total=0; reuse_total=0; cell_pairs={}
    views=dict(cases)
    views.update({'original__'+label:case for label,case in original.items()})
    for label,case in cases.items():
        data=case['binding']; bound=data['bound']; grade=case['grades']
        f.require(data['settings']==system['settings'],'full actual finite settings')
        r,adapter,packets=selected_context(base,label,case,manifest['input'])
        originals=f.unpickle(base/'accepted-cases'/label/'assembly.pickle')['result']['ROWS']
        selected=packets['assembly']['result']['ROWS']; addresses={tuple(v['address']):v['record'] for v in grade['records'].values()}
        pairs=[]; changed=0; terms=0
        for old,cell in zip(originals,selected):
            f.require(same({k:v for k,v in old.items() if k!='NONLOCAL'},{k:v for k,v in cell.items() if k!='NONLOCAL'}),'all native cell fields/probes/units preserved')
            f.require(len(old['NONLOCAL'])==len(cell['NONLOCAL']),'all selected native terms present')
            for ti,((oi,ov),(ni,nv)) in enumerate(zip(old['NONLOCAL'],cell['NONLOCAL'])):
                expected=addresses['cell',cell['ROW'],cell['COLUMN'],ti]['ORIGINAL']
                f.require(same(oi,ni) and same(nv,expected),'literal changed coefficient enters independent native cell route')
                pairs.append((oi,ov,ni,nv,expected));changed+=int(not same(ov,nv));terms+=1
        f.require(terms==len(grade['termJoins']),'all native cell terms retained')
        for item in grade['records'].values():
            kind,*address=item['address']; rec=item['record']
            if kind in ('factor','source'): f.require(set(rec['COEFFICIENTS'])<={(0,0,0)},'actual grade-free quadrature support')
            if kind=='factor': f.require(same(rec['ORIGINAL'],bound['rows'][address[0]]['factors'][address[1]]['symbolicCoefficient']),'changed factor grade/binding identity')
            if kind=='source': f.require(same(rec['ORIGINAL'],bound['sources'][0,address[0]]['symbolicAmplitude']),'changed source grade/binding identity')
        for term in grade['termJoins']:
            f.require(same(term['originalIntegral'],bound['rows'][term['integralIndex']]['original']),'complete ordered source address')
            f.require(all(v['factorGrades'][1:]==((0,0,0),(0,0,0)) for factor in term['factors'] for v in factor['combinations']),'full term convolution support')
        if label!=BASELINE:
            fresh=set(case['newRows']); reused={v['row'] for v in case['reusedRows']}
            f.require(not fresh&reused and fresh|reused==set(range(len(bound['rows']))),'complete computed numerical row partition')
            for route in case['reusedRows']:
                other=views[route['fromCase']]['binding'];row=bound['rows'][route['row']];old=other['bound']['rows'][route['fromRow']]
                f.require(same(binding.row_signature(row,data),binding.row_signature(old,other)),'complete consumed source/profile/field/limit/measure row join')
                if route['fromCase'] in pools:
                    array=pools[route['fromCase']][route['fromRow']]
                    f.require(array.shape==(129,129) and np.isfinite(array).all(),'actual complete reusable row array')
                else: f.require(route['fromRow'] in cases[route['fromCase']]['newRows'],'explicit earlier newly computed row route')
            for i in fresh: f.require(len(bound['rows'][i]['limits'])<=2,'no new triple integral')
            row=bound['rows'][0];limit=row['limits'][0]
            wrong=dict(row,limits=(sp.Tuple(limit[0],limit[1],limit[2]+1),*row['limits'][1:]))
            f.require(not same(binding.row_signature(wrong,data),binding.row_signature(row,data)),'actual changed-limit control');controls+=1
            if changed:
                pair=next(v for v in pairs if not same(v[1],v[3])); f.require(not same(pair[1],pair[4]),'actual omitted cell reversal control');controls+=1
            new_total+=len(fresh);reuse_total+=len(reused)
        cell_pairs[label]=pairs
        summaries[label]={'rows':len(bound['rows']),'sources':len(data['jets']),'terms':terms,'changedCellCoefficients':changed,'newRows':case['newRows'],'reusedRows':len(case['reusedRows'])}
    f.atomic_pickle(base/'selected-cell-input-joins.pickle',cell_pairs)
    _,tail=matrix_tail(manifest)
    result={'cases':summaries,'newUnionRows':new_total,'missingCaseReuses':reuse_total,'baselineControlRows':80,
            'fieldCoefficients':129,'unknowns':645,'actualMutationControls':controls,'savedLayoutPackets':len(evidence),'nativeTail':tail,
            'newQuadratureNodes':0,'newMatrixAssemblies':0,'newModes':0,'newResponseSolves':0}
    f.require(new_total==27 and reuse_total==193,'actual accepted binding union census')
    f.save(base/'preflight.json',result);return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--focused',action='store_true');p.add_argument('--resume-focused',type=Path);args=p.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    f.require(args.focused!=(args.resume_focused is not None),'focused preparation or accepted focused reuse')
    manifest,labels,system=load(base,args.resume_focused)
    focus=json.loads((base/'preflight.json').read_text()) if args.resume_focused else preflight(base,manifest,labels,system)
    summaries={};new_nodes=0
    if not args.focused:
        cases,rows=cases_and_rows(base,labels);complete,proof=matrix_tail(manifest)
        f.require(proof==focus['nativeTail'],'exact native tail/helper joins')
        for label,case in cases.items():
            if label==BASELINE: continue
            target=base/'cases'/label;target.mkdir(parents=True,exist_ok=False)
            _,new_nodes=complete(base,label,case,system,rows,summaries,new_nodes,None,target);gc.collect()
    for n,v in manifest['sourceFiles'].items(): f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'source/frozen pre/post identity')
    for n,v in manifest['inputPackets'].items(): f.require(f.digest(Path(n))==v,'original input pre/post identity')
    for n,v in manifest['copiedInputs'].items(): f.require(f.digest(base/n)==v,'copied operand pre/post identity')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    checks={**manifest,'mode':'focused' if args.focused else 'construct','status':'PASSED_FIRST_JET_MATRIX_INPUTS' if args.focused else 'COMPLETED_CASE_FIRST_JET_INTERIORS',
            'preflight':focus,'cases':summaries,'newNodes':new_nodes,'newModeConstructions':0,'newResponseSolves':0,'artifacts':artifacts,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__': main()
