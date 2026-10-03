#!/usr/bin/env python3
"""Local part of two bounded Gaussian packet pairings. No pressure action."""
import argparse
import ast
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback

ROOT=Path('/var/projects/toy_physics')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','containment','Journal','decode','one_symbol')
G=((0,0),(1,0),(0,1),(1,1))
COMPONENTS=('NATIVE_FLAT','NATIVE_HEIGHT','NATIVE_SLOPE','NATIVE_MIXED_ITERATION','INHERITED_DIRECT_WHOLE_OFF_DIAGONAL')
ZERO={'text':'0','srepr':'Integer(0)'}


def require(v,message):
    # No bool(value): unknowns and numeric truthiness must never become proofs.
    exact_true=getattr(getattr(globals().get('sp'),'S',None),'true',None)
    if v is not True and not (exact_true is not None and v is exact_true):
        raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def save(path,value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())


def definitions(text,names):
    selected=[n for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
    require({n.name for n in selected}==set(names),'exact inert helper set')
    return ast.Module(body=selected,type_ignores=[])


def verify_gate(path,manifest_path,m):
    g=json.loads(Path(path).read_text())
    require(g['status']=='READY_FOR_ONE_PACKET_LOCAL_ACTION','gate status')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest pins')
    require(g['sourcePins']==m['sourcePins'],'source census')
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and
        g['supervisor']==str(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),'actual containment route')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source '+p)
    for key in ('sharedGuard','supervisor','launcher','buildReviewRecord','authority'):
        require(sha(g[key])==g[key+'Sha256'],'gate '+key)
    require(g['launcher']==m['launcher'] and g['buildReviewRecord']==m['reviewRecordWillBe'] and g['authority']==m['executionAuthority'],'exact document routes')
    r=json.loads(Path(g['buildReviewRecord']).read_text())
    require(r['independentBuildClearance'] is False and g['independentBuildClearance'] is False and r['allChecksPassed'] is True,'honest literal build disposition')
    require(r['reports']['claude']['literalVerdict']=='NEEDS REVISION' and
        r['reports']['grok']['literalVerdict']=='CLEAR FOR THIS PACKET-ACTION LOCAL BUILD','literal source reports')
    require(g['localToolingRepairRecord']==m['localToolingRepairRecord'] and
        sha(g['localToolingRepairRecord'])==g['localToolingRepairRecordSha256'],'exact tooling record')
    repair=json.loads(Path(g['localToolingRepairRecord']).read_text())
    require(g['localToolingExecutionAuthority'] is True and repair['localToolingExecutionAuthority'] is True and
        repair['independentBuildClearance'] is False and repair['testsPassed'] is True,'standing tested tooling authority, not author CLEAR')
    for key in ('workerSha256','manifestSha256','librarySha256','launcherSha256'):
        require(repair['reviewed'][key]==r[key] and repair['current'][key]==g[key],'review/repair/gate '+key)
    for key in ('sharedGuardSha256','supervisorSha256'):
        require(r[key]==g[key],'unchanged containment '+key)
    method=json.loads(Path(m['methodRecord']).read_text())
    require(method['jointIndependentMethodClearance'] is True and method['methodSha256']==sha(m['methodPath'])==r['methodSha256'],'exact method')
    a=json.loads(Path(g['authority']).read_text())
    require(a['scope']==m['scope']==g['scope'] and a['boundedInstrumentAuthorized'] is True and
        a['scienceExecutionsAuthorized']==g['scientificRunsAuthorized']==1 and
        a['automaticScientificRetry'] is False and a['noDeadline'] is True and g['durationLimits'] is None,'standing bounded authority')
    return g


def verify_invocation(args,g,argv):
    tail=[str(Path(__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(list(argv)==tail and g['command'][-len(tail):]==tail,'worker argv')
    require(args.out.resolve()==Path(g['outputDirectory']).resolve(),'output route')



def certificate_match_indices(coefficients,certificates):
    """Exact saved-value multiset matching; list order is not scientific content."""
    require(len(coefficients)==len(certificates),'coefficient/certificate multiplicity')
    require(all(v['finite'] is True for v in certificates),'all coefficient certificates finite')
    remaining=list(range(len(certificates)));indices=[]
    for coefficient in coefficients:
        matches=[i for i in remaining if certificates[i]['value']==coefficient['value']]
        require(bool(matches),'missing exact coefficient certificate')
        index=matches[0];remaining.remove(index);indices.append(index)
    require(not remaining,'unused certificate')
    return indices


def dimension_of(expr,dimensions):
    zero=(sp.S.Zero,)*3
    if expr.is_Number or expr==sp.I:return zero
    if expr.is_Symbol:
        require(expr.name in dimensions,'native unit symbol '+expr.name)
        return tuple(map(sp.Rational,dimensions[expr.name]))
    if expr.is_Add:
        ds=[dimension_of(v,dimensions) for v in expr.args if v!=0]
        require(bool(ds) and all(d==ds[0] for d in ds),'native additive unit mismatch');return ds[0]
    if expr.is_Mul:return tuple(map(sum,zip(*(dimension_of(v,dimensions) for v in expr.args))))
    if expr.is_Pow and expr.exp.is_Rational:return tuple(expr.exp*d for d in dimension_of(expr.base,dimensions))
    raise ValueError('unsupported dimension node '+str(expr.func))


def run_science(m,J,ns):
    D=ns['decode'];raw={};copies={}
    for alias,r in m['savedInputs'].items():
        p=Path(r['path']);q=J.out/'saved'/alias;q.parent.mkdir(parents=True,exist_ok=True)
        require(sha(p)==r['sha256'],'saved source '+alias);shutil.copyfile(p,q)
        require(sha(q)==r['sha256'],'saved copy '+alias)
        raw[alias]=json.loads(q.read_text());copies[alias]={'source':str(p),'path':str(q.relative_to(J.out)),'sha256':r['sha256'],'bytes':q.stat().st_size}
        ns['replace_json'](J.out/'saved-copy-index.json',copies)
    physical=raw['preflight/physical-plan.json'];context=raw['local/context.json']
    require(physical['frequency']==3 and D(physical['cs'])==sp.sqrt(6)/2 and D(physical['kappa'])==sp.sqrt(595)/10,'actual physical frequency/speed')
    require(D(physical['centers'])==[-sp.Rational(5,2),sp.Rational(5,2)] and physical['s']==8,'same packet')
    J.emit('new-local-tail-arguments',{'actualPacket':physical,'carrierUpper':sp.Rational(5,2),'kappaSquared':D(physical['kappa'])**2,'tailRule':'positive coefficient majorants; Gaussian moment integration by parts; assessed analytic inequality'})
    require(0<D(physical['kappa'])<sp.Rational(5,2),'actual carrier majorant domain')
    require(context['physical']==raw['physical-input.json'] and context['frequencyOverride']=={'old':'1','actual':3},'original materials and held frequency')
    require(context['physical']['parameters']['L_W']=='10' and context['effectiveSpeedOnlyInPressure'] is True,'local scale/speed scope')
    selection=raw['selected/local-cells.json'];allcells=raw['local/all-local-cells.json'];cells=selection['selected']
    require(selection['sourceSha256']==m['savedInputs']['local/all-local-cells.json']['sha256'] and selection['sourceBytes']==m['savedInputs']['local/all-local-cells.json']['bytes'],'original complete local file')
    require(cells==[v for v in allcells if v['row']=='THETA_BALANCE' and v['field']=='e_W'] and len(cells)==16,'complete selected local cells')
    require([(c['xOrder'],tuple(c['grade'])) for c in cells]==list(itertools.product(range(4),G)),'all16 local grades/orders')
    partition=raw['native/THETA-partition.json'];native=raw['native/THETA-row.json'];joins=raw['units/source-text-joins.json']
    require(partition['source']==native['source']==joins['provenance'] and hashlib.sha256(native['fullConstructor'].encode()).hexdigest()==joins['rowHashes']['THETA_BALANCE'],'native units/local same full source')
    source=Path(partition['source']['source'])
    with source.open('rb') as f:line=next(line for i,line in enumerate(f,1) if i==partition['source']['valueLine'])
    require(hashlib.sha256(line).hexdigest()==partition['source']['sourceLineSha256'],'actual immutable source line')
    units=D(raw['units/source-units.json']);rules=raw['local/native-wave-profile-contract.json'];dimensions=dict(rules['dimensions'])
    require(units['nativeRows'][3]==[-3,-1,1] and units['nativeFields'][4]==[0,0,0] and units['frame']==context['physical']['unit_frame'],'actual THETA/eW unit frame')
    for name,entries in raw['units/merged.json']['generatedEntries'].items():
        values=[D(e['unit']) for e in entries]
        require(len(values)==4 and all(v==values[0] for v in values),'inherited actual registry agreement')
        require(name not in dimensions or dimensions[name]==values[0],'native registry conflict');dimensions[name]=values[0]
    for name in ['epsilon_shape','eta_bg','sigma_W','e_W']:require(dimensions[name]==[0,0,0],'dimensionless native grade/field')
    children={}
    for alias in m['localBatches']:
        inp=raw[alias.replace('-return','-input')]
        require(inp['children']==[next(v for v in partition['children'] if v['childIndex']==e['childIndex']) for e in raw[alias]],'saved batch actual input children')
        for ev in raw[alias]:
            require(ev['completed'] is True and ev['row']=='THETA_BALANCE' and ev['childIndex'] not in children,'complete unique saved native child')
            original=next(v for v in partition['children'] if v['childIndex']==ev['childIndex'])
            require(ev['sourceConstructor']==original['constructorText'] and ev['sourceSha256']==original['sha256']==hashlib.sha256(ev['sourceConstructor'].encode()).hexdigest(),'actual original child operand')
            require(all(v['cancelled']==ZERO for v in ev['identities']),'inherited child zeros')
            children[ev['childIndex']]=ev
    require(set(children)==set(partition['localChildIndices']),'complete original local partition')
    rowunit=tuple(units['nativeRows'][3]);fieldunit=tuple(units['nativeFields'][4]);specs=[];used=[]
    for index,cell in enumerate(cells):
        n=cell['xOrder'];g=cell['grade']
        matches=[v for v in children.values() if v['fieldColumn']==4 and v['xOrder']==n]
        summands=[next(r['value'] for r in v['mappedGrades'] if r['grade']==g) for v in matches]
        J.emit('cell-'+str(index)+'-ancestry',{'cell':cell,'childIds':[v['childIndex'] for v in matches],'summands':summands})
        require(cell['sourceChildren']==[v['childIndex'] for v in matches] and cell['summands']==summands,'actual native child/cell summands')
        require(all(v['cancelled']==ZERO for v in cell['identities']),'saved whole-cell zeros')
        require(cell['identities'][0]['left']==cell['coefficient'],'saved actual coefficient identity')
        for ev in matches:
            if ev['childIndex'] in used:continue
            used.append(ev['childIndex']);jet=ev['jet'];coefficient=D(ev['nativeCoefficient']);names={s.name for s in coefficient.free_symbols}
            require(not any(n.lower().startswith(('c_s','cs_')) or n.lower()=='cs' or 'speed' in n.lower() for n in names),'actual unbound selected local speed independence')
            d0=dimension_of(coefficient,dimensions)
            waveunit=tuple(map(sp.Rational,ev['nativeJetDimension']));mapped=tuple(d0[i]+(-sum(jet['spatialOrders'][1:]) if i==0 else -jet['timeOrder'] if i==1 else 0) for i in range(3))
            expect=tuple(rowunit[i]-fieldunit[i]+(ev['xOrder'] if i==0 else 0) for i in range(3))
            J.emit('new-unit-child-'+str(ev['childIndex']),{'originalChild':ev,'rawCoefficientUnit':d0,'savedWaveUnit':waveunit,'mappedCoefficientUnit':mapped,'expected':expect,'rowUnit':rowunit,'noCoefficientRecalculation':True})
            require(tuple(d0[i]+waveunit[i] for i in range(3))==rowunit and mapped==expect,'actual native raw unit and mapped jet unit')
        poly=cell['polynomial']
        J.emit('cell-'+str(index)+'-certificate-originals',{'coefficients':poly['coefficients'],'certificates':poly['coefficientCertificates'],'matching':'exact saved values with multiplicity, never list position'})
        matching=certificate_match_indices(poly['coefficients'],poly['coefficientCertificates'])
        J.emit('cell-'+str(index)+'-certificate-match',{'certificateIndicesInCoefficientOrder':matching,'originalCertificatesUnchanged':True})
        coeff=D(poly['coefficients']);cert=D(poly['coefficientCertificates'])
        for ci,vi in enumerate(matching):
            v=cert[vi]
            J.zero('new-cell-'+str(index)+'-certificate-components-'+str(ci),v['value'],v['real']+sp.I*v['imaginary'])
            require(v['value']==coeff[ci]['value'] and v['finite'] is True,'actual decoded certificate equality')
        degree=max([c['order'] for c in coeff],default=0);vector=[sp.S.Zero]*(degree+1)
        for c in coeff:vector[c['order']]=c['value']
        pairs=[]
        for i,v in enumerate(vector):
            re,im=v.as_real_imag();require(re.is_Rational is True and im.is_Rational is True,'exact rational local transport')
            J.zero('new-cell-'+str(index)+'-transport-'+str(i),v,sp.Rational(str(re))+sp.I*sp.Rational(str(im)));pairs.append([str(re),str(im)])
        T=sp.Symbol('full_weak_tanh_variable',real=True)
        # This is a new numeric adapter join, not the old coefficient derivation.
        original=D(poly['polynomial']);symbols=original.free_symbols
        require(not symbols or len(symbols)==1 and next(iter(symbols)).name==T.name,'actual polynomial symbol')
        if symbols:T=next(iter(symbols))
        J.zero('new-cell-'+str(index)+'-vector',sum(v*T**i for i,v in enumerate(vector)),original)
        unit=tuple(rowunit[i]-fieldunit[i]+(n if i==0 else 0) for i in range(3))
        J.emit('cell-'+str(index)+'-adapter',{'cell':cell,'coefficientsAscending':pairs,'nativeSourceChildren':cell['sourceChildren'],'coefficientUnit':unit,'dxUnit':[1,0,0],'rowUnit':rowunit,'dualUnit':'U_dual_THETA (kept formal)','pairingUnit':'U_dual_THETA * M_ref * L_ref^-2 * T_ref^-1','epsilonAlreadyExtracted':True,'timeAndTangentsAlreadyAbsorbed':True})
        require(cell['epsilonPower']==(0 if poly['zero'] else 1),'epsilon exactly once or explicitzero')
        specs.append({'cellIndex':index,'xOrder':n,'grade':g,'coefficientsAscending':pairs,'zero':poly['zero']})
    J.emit('local-unit-conclusion',{'pairingUnit':'U_dual_THETA * M_ref * L_ref^-2 * T_ref^-1','nativeRowUnit':rowunit,'nativeFieldUnit':fieldunit,'checkedNativeChildren':used,'pressureSummandUnitsEstablished':False,'powerInterpretation':False})
    for name,receipt in raw['rules/extraction.json']['records'].items():
        require(receipt['sha256']==m['savedInputs']['rules/'+name+'.json']['sha256'] and receipt['record']=='rules/'+name,'original full rule bytes')
    require(raw['rules/extraction.json']['database']==raw['fourier/journal-receipt.json'],'original Fourier rule provenance')
    import importlib.util
    def module(key,name):
        p=Path(m[key]);require(sha(p)==m['sourcePins'][str(p)],'pinned runtime '+key)
        s=importlib.util.spec_from_file_location(name,p);v=importlib.util.module_from_spec(s);s.loader.exec_module(v);return v
    old=module('storageLibrary','packet_storage');inner=module('quadratureLibrary','packet_quadrature');lib=module('librarySource','packet_local')
    import mpmath
    require(str(Path(mpmath.__file__).resolve())==m['runtimeLibrary']['initPath'] and mpmath.__version__==m['runtimeLibrary']['version'],'actual runtime')
    store=old.DurableStore(J.out/'local-evidence.sqlite')
    try:
        numerical=lib.LocalEvaluator(store,old,inner,{name:raw['rules/'+name+'.json'] for name in ('A-GL24','A-GL48','B-G7-K15')})
        result=numerical.run(specs)
    finally:
        store.close();J.emit('numerical-journal-receipt',{'path':'local-evidence.sqlite','bytes':(J.out/'local-evidence.sqlite').stat().st_size,'sha256':sha(J.out/'local-evidence.sqlite')})
    return {'status':'BOUNDED_LOCAL_PACKET_ACTIONS_COMPLETE_PRESSURE_PENDING','cells':16,'packets':2,'localResult':result,'completePacketAction':False,'pressureActionEvaluated':False,'currentOrLoss':None,'scientificAcceptance':False}


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);args=p.parse_args()
    m=json.loads(args.inputs.read_text());g=verify_gate(args.gate,args.inputs,m);verify_invocation(args,g,sys.argv)
    pins={**m['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False);J=None;code=1;started=time.monotonic();result={}
    try:
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(m['helperSource']).read_text(),HELPERS),'pinned-inert-helpers','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        ns.update(sp=sp,Str=Str);J=ns['Journal'](args.out)
        result=J.stage('packet-local',{'manifestSha256':g['manifestSha256'],'buildReviewSha256':g['buildReviewRecordSha256']},lambda:run_science(m,J,ns));code=0
    except BaseException:
        result={'status':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        index=args.out/'saved-copy-index.json'
        if index.exists():
            for r in json.loads(index.read_text()).values():pins[str(args.out/r['path'])]=r['sha256']
        records={}
        for path,expected in pins.items():
            try:records[path]={'expected':expected,'actual':sha(path),'error':None}
            except OSError as e:records[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',records)
        if any(r['actual']!=r['expected'] for r in records.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False)
        save(args.out/'checks.json',result if J is None else J.encode(result));sys.stdout.write((args.out/'checks.json').read_text())
    return code


if __name__=='__main__':sys.exit(main())
