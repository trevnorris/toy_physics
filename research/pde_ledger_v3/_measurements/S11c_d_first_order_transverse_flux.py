#!/usr/bin/env python3
"""Native finite transverse end flux and linear survival diagnostic. No receiving field or leakage coefficient."""
import argparse, ast, hashlib, json, os, re, resource, shutil, sys, time, traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','containment','Journal','decode')
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
GRADES=('00','10','01')
VERDICT='CLEAR FOR THIS NATIVE TRANSVERSE END-FLUX BUILD'

def require(value,message):
    if value is not True: raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def save(path,value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())

def definitions(text):
    nodes=[n for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in HELPERS]
    require({n.name for n in nodes}==set(HELPERS),'unchanged journal/containment helper census')
    return ast.Module(body=nodes,type_ignores=[])

def verify_gate(path,manifest_path,manifest):
    g=json.loads(Path(path).read_text())
    require(g['status']=='READY_FOR_ONE_FIRST_ORDER_TRANSVERSE_FLUX_INSTRUMENT','one finite receiving prerequisite gate')
    for key,p in [('worker',__file__),('manifest',manifest_path),('launcher',manifest['launcher'])]:
        require(g[key+'Sha256']==sha(p),'actual '+key)
    require(g['sourcePins']==manifest['sourcePins'],'source census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'pin '+p)
    for key in ['buildReviewRecord','methodRecord','executionAuthority']:
        require(sha(g[key])==g[key+'Sha256'],'gate '+key)
    require(g['methodRecord']==manifest['methodRecord'],'actual method record')
    review=json.loads(Path(g['buildReviewRecord']).read_text())
    require(sha(review['report'])==review['reportSha256'],'literal build report bytes')
    require(review['reviewers']==['claude'] and review['literalVerdict']==VERDICT and review['buildAssessed'] is True,'actual scoped Claude build')
    for key in ['workerSha256','manifestSha256','launcherSha256']:
        require(review[key]==g[key],'review '+key)
    method=json.loads(Path(g['methodRecord']).read_text())
    require(sha(method['report'])==method['reportSha256'],'literal method report bytes')
    require(method['methodAssessed'] is True and method['literalVerdict']=='CLEAR FOR THIS FIRST-ORDER MATCHED RECEIVING AND FACE METHOD','actual method')
    require(sha(method['methodSource'])==method['methodSha256']==manifest['methodSha256'],'reviewed method bytes')
    auth=json.loads(Path(g['executionAuthority']).read_text())
    require(auth['scienceExecutionsAuthorized']==1 and auth['scope']==manifest['scope'] and auth['reviewers']==['claude'],'bounded authority')
    require(auth['AGENTSSha256']==sha(ROOT/'AGENTS.md'),'current authority')
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and g['supervisor']==str(M/'S11c_d_end_normalization_run.py'),'actual containment paths')
    require(sha(g['sharedGuard'])==g['guardSha256'] and sha(g['supervisor'])==g['supervisorSha256'],'actual containment bytes')
    require(g['scope']==manifest['scope'] and g['scientificRunsAuthorized']==1 and g['pooledExecution'] is True,'exact scope')
    require(manifest['resources']=={'memoryBytes':4*1024**3,'poolGiB':16,'zeroSwap':True,'cpuCount':1,'threads':1,'tasksMax':32,'hostReserveGiB':4,'durationLimits':None},'fixed resources')
    return g

def verify_invocation(args,gate,argv):
    expected=[str(Path(__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(list(argv)==expected and gate['command'][-len(expected):]==expected,'exact worker command')
    require(str(args.out.resolve())==gate['outputDirectory'],'exact output')

def scientific_work(manifest,J,decode):
    copied={};cache={};directory=J.out/'saved';directory.mkdir()
    for alias,record in manifest['savedFiles'].items():
        src=Path(record['path']);require(sha(src)==record['sha256'],'saved bytes '+alias)
        dst=directory/alias;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
        require(sha(dst)==record['sha256'],'copy bytes '+alias)
        copied[alias]={'source':str(src),'path':str(dst.relative_to(J.out)),'sha256':record['sha256']}
        tmp=J.out/'saved-copy-index.new';tmp.write_text(json.dumps(copied,indent=2)+'\n');tmp.replace(J.out/'saved-copy-index.json')
    def get(alias):
        if alias not in cache:cache[alias]=decode(json.loads((directory/alias).read_text()))
        return cache[alias]
    def clean(v):return v.applyfunc(sp.cancel) if isinstance(v,sp.MatrixBase) else sp.cancel(v)
    def zero(name,left,right):
        if isinstance(left,sp.MatrixBase):
            require(isinstance(right,sp.MatrixBase) and left.shape==right.shape,'shape '+name)
            J.emit(name+'-matrix-operands',{'left':left,'right':right})
            for i in range(left.rows):
                for j in range(left.cols):J.zero(name+'-'+str(i)+'-'+str(j),left[i,j],right[i,j])
        else:J.zero(name,left,right)
    def sym(expr,name):
        found=[s for s in expr.free_symbols if s.name==name];require(len(found)==1,'unique native symbol '+name);return found[0]
    def pairs(d):return dict(d['mappingPairs']) if set(d)=={'mappingPairs'} else dict(d)
    # Only the established saved-object decoder and structural predicate execute.
    # No end/current constructor or old validator is imported or called.
    import builtins,importlib,io,pickle
    import numpy as np
    from collections import OrderedDict
    codec={'pickle':pickle,'OrderedDict':OrderedDict,'builtins':builtins,'np':np,'sp':sp,'importlib':importlib,'io':io,'require':require}
    def inert(path,names,ns):
        nodes=[n for n in ast.parse(Path(path).read_text()).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
        require({n.name for n in nodes}==set(names),'existing definition census')
        exec(compile(ast.Module(body=nodes,type_ignores=[]),'unchanged-saved-operand-reader','exec'),ns)
    inert(manifest['codecSource'],('SavedCodec','decode'),codec)
    structural={'np':np,'sp':sp};inert(manifest['structureSource'],('exact_structure',),structural)
    opaque={n:codec['decode']((directory/n).read_bytes()) for n in manifest['opaqueAliases']}
    original=opaque['native/RIGHT-original.pickle'];right=opaque['native/RIGHT-return.pickle'];actual=opaque['native/RIGHT-input.pickle']['native']
    receipt=opaque['native/RIGHT-restore-input.pickle'];producer=get('native/RIGHT-producer-checks.json')
    J.emit('right-native-restoration-operands',{'receipts':{k:copied[k] for k in manifest['opaqueAliases']},'originalRestoreInput':receipt,'producer':producer,'actualSelector':['native'],'expectedSelector':[],'originalSelector':[],'predicate':'unchanged exact_structure on full live objects; receipts alone never decide'})
    require(receipt['pin']['sha256']==sha(directory/'native/RIGHT-original.pickle')==producer['objectsSha256'],'original RIGHT object bytes')
    joins={'actualExpected':structural['exact_structure'](actual,right),'originalExpected':structural['exact_structure'](original,right)}
    J.emit('right-native-restoration-return',joins);require(all(v is True for v in joins.values()),'full original RIGHT live-object joins')
    require(producer['provenance']['case']=='LAB_HELD_RHO4_CONSTANT' and producer['provenance']['end']=='RIGHT','original native case')
    require(producer['provenance']['engineSha256']==sha(directory/'native/frozen-engine.py') and producer['provenance']['instrumentSha256']==sha(directory/'native/frozen-producer.py'),'original source origins')
    leftRecord=get('native/LEFT-original.json');require(leftRecord['passed'] is True,'inherited LEFT full object join');left=leftRecord['actual']
    context=get('source/extended-binding-context.json');values=context['numeric'];physical=get('source/physical-input.json');eps=context['epsilon']
    require(context['physical']==physical and context['frequencyOverride']=={'old':'1','actual':3},'original physical input and held omega')
    for name,value in values.items():require(value==sp.Rational(physical['parameters'][name]) if name!='omega' else value==3,'actual physical binding '+name)
    old=get('source/incident-columns.json');U=old['columns'];S=old['chart']['matrixWeakFromUniform'];p=sp.sympify(old['selectedPoint']['mapping']['uniformNormal'])
    matched=get('matching/matched-transverse-amplitudes.json');end=get('receiving/end-matching-prerequisite.json');T1=matched['T1'];B=end['B'];l=sym(B,'receiving_block_l');delta=end['deltaP'];Bp=end['BprimeAtIncident']
    nativeInput=opaque['native/RIGHT-input.pickle'];rawSource=get('native/RIGHT-raw-source-binding.json')
    actualP=nativeInput['pairing']['result']['CLOSED_PENCIL_LEGS'][0]
    J.emit('right-current-and-end-source-family',{'currentOriginalReceipt':copied['native/RIGHT-return.pickle'],'nativeReconstructionInputReceipt':copied['native/RIGHT-input.pickle'],'actualClosedPencil':actualP,'savedOriginalEndSource':rawSource['nativeSource'],'savedBinding':rawSource,'sameRestoredNativeObject':joins,'originalFrequency':nativeInput['frequency'],'originalSpeed':nativeInput['cs'],'note':'This joins original unbound pencil/current families; the historical finite speed/root is not reused as this expansion.'})
    require(actualP==rawSource['nativeSource'],'native current and saved endpoint pencil share the actual original source family')
    require(pairs(rawSource['endpoints'])==pairs(right['profileBindings']),'actual end profiles')
    require(matched['B']==B and matched['DB']==end['DB'] and matched['deltaP']==delta and matched['Bprime']==Bp,'completed matching arguments')
    require(matched['amplitudesInOriginalIncidentCoordinates'] is True and matched['physicalFluxNormalizationApplied'] is False,'actual amplitude convention')
    require(get('receiving/incident-chart-match-matrix-operands.json')['right']==U,'original incident chart argument')
    restored=[]
    for prefix,rows,cols,folder in [('full-five-row-matched-end',5,2,'receiving'),('right-field-finite-matching',5,2,'matching')]:
        for i in range(rows):
            for j in range(cols):
                name=prefix+'-'+str(i)+'-'+str(j);op=get(folder+'/'+name+'-input.json');ret=get(folder+'/'+name+'-return.json')
                require(ret['cancelled']==0,'inherited '+name);restored.append({'name':name,'operands':op,'return':ret})
    require(get('receiving/full-five-row-matched-end-matrix-operands.json')['left']==end['fullEndFirstOrder'],'actual inherited full endpoint field')
    require(get('matching/right-field-finite-matching-matrix-operands.json')['left']==matched['rightFiniteField'],'actual inherited finite field')
    J.emit('inherited-matching-and-end',{'records':restored,'completeEnd':end,'matching':matched,'incident':old,'functionsReplayed':False})
    # Keep independent momentum legs until their physical diagonal is taken.
    km,kp,lam=sp.symbols('flux_left_momentum flux_right_momentum flux_lambda',real=True);q=sp.Symbol('flux_depth')
    eta,sigma=sp.symbols('flux_eta flux_sigma',real=True)
    h1,h2=values['s11cdTangentialMomentum1'],values['s11cdTangentialMomentum2']
    native={'LEFT':left,'RIGHT':right};bound={};mappings=[]
    def bind(expr,nat,name,normal=None,amplitudes=None):
        endpoints=pairs(nat['profileBindings']);v=expr.xreplace(endpoints);mp={}
        for a in v.free_symbols:
            n=a.name
            if n=='epsilon_shape':mp[a]=eps
            elif n=='s11cdCurrentLeftMomentum':mp[a]=km
            elif n=='s11cdCurrentRightMomentum':mp[a]=kp
            elif n=='s11cdSpectralNormalMomentum':require(normal is not None,'explicit native wave leg');mp[a]=normal
            elif n=='s11cdAcousticRightNormalMomentum':mp[a]=q
            elif n=='eta_bg':mp[a]=eta
            elif n=='sigma_W':mp[a]=sigma
            elif n.startswith('s11cdCurrentPlusAmplitude'):
                index=int(n.removeprefix('s11cdCurrentPlusAmplitude'));require(amplitudes is not None and 0<=index<5,'actual five native amplitudes');mp[a]=amplitudes[index]
            elif n=='c_s0':mp[a]=sp.sqrt(6)/2
            elif n in values:mp[a]=values[n]
            else:raise ValueError('unbound native flux input '+n)
        result=v.xreplace(mp)
        J.emit(name+'-binding',{'original':expr,'profileBindings':list(endpoints.items()),'numericMapping':list(mp.items()),'result':result,'epsilonStripped':False,'etaSigmaIndependent':True})
        require(not result.atoms(sp.Limit,sp.Integral,sp.Derivative,sp.Subs),'native current binding unresolved')
        mappings.append({'name':name,'mapping':list(mp.items())});return result
    units=opaque['native/uniform-units.pickle'];unitOriginal=opaque['native/uniform-original.pickle'];unitReceipt=opaque['native/uniform-restore-input.pickle']
    J.emit('original-current-unit-operands',{'originalReceipt':unitReceipt,'restoredUnits':units,'originalSelected':{k:unitOriginal[k] for k in ('fieldUnits','currentUnit')},'originalFile':copied['native/uniform-original.pickle']})
    require(unitReceipt['pin']['sha256']==sha(directory/'native/uniform-original.pickle'),'original uniform unit source bytes')
    require(structural['exact_structure']({k:unitOriginal[k] for k in ('fieldUnits','currentUnit')},units),'actual inherited native unit values')
    require(len(units['fieldUnits'])==5 and tuple(units['currentUnit'])==(0,-3,1),'inherited field/current unit convention')
    for side,nat in native.items():
        slab,con=nat['slab'],nat['conservative'];registry=pairs(nat['knownDimensions'])
        J.emit(side+'-native-current-source',{'slabCurrent':slab['SLAB_CURRENT_MATRIX'],'conservativeCurrent':con['SLAB_CURRENT_MATRIX'],'massCorrection':slab['MASS_RATE_CORRECTION_MATRIX'],'decompositionResidual':slab['CURRENT_DECOMPOSITION_RESIDUAL'],'harmonicNormalization':slab['HARMONIC_VARIATION_NORMALIZATION'],'knownDimensions':list(registry.items()),'units':units,'sourceBoundaryWork':con['NORMAL_BOUNDARY_WORK'],'actualTimeBoundary':slab['ACTUAL_TIME_BOUNDARY'],'chemicalTransport':slab['CHEMICAL_TRANSPORT_CURRENT']})
        require(slab['CURRENT_DECOMPOSITION_RESIDUAL']==0 and slab['POLARIZATION_RESIDUAL']==0 and slab['CORRECTION_POLARIZATION_RESIDUAL']==0,'inherited original current identities')
        require(slab['HARMONIC_VARIATION_NORMALIZATION']==sp.Rational(1,4),'native real-field harmonic factor retained')
        require(tuple(registry[sym(slab['SLAB_CURRENT_MATRIX'],'epsilon_shape')])==(0,0,0),'original epsilon units')
        amplitudeUnits=[]
        for i in range(5):
            hits=[k for k in registry if getattr(k,'name',None)=='s11cdCurrentPlusAmplitude'+str(i)];require(len(hits)==1,'native field unit identity');amplitudeUnits.append(tuple(registry[hits[0]]))
        require(amplitudeUnits==[tuple(v) for v in units['fieldUnits']],'original native field/current convention')
        J.emit(side+'-field-unit-join',{'originalAmplitudeUnits':amplitudeUnits,'uniformFieldUnits':units['fieldUnits'],'currentUnit':units['currentUnit'],'newDimensionInference':False})
        entries={}
        for tag,expr in [('slab',slab['SLAB_CURRENT_MATRIX']),('conservative',con['SLAB_CURRENT_MATRIX']),('mass',slab['MASS_RATE_CORRECTION_MATRIX'])]:
            v=bind(expr,nat,side+'-'+tag);stripped=v.applyfunc(lambda a:sp.expand(a).coeff(eps,2));zero(side+'-'+tag+'-epsilon-once',v,eps**2*stripped)
            entries[tag]=S*stripped*S.T
        zero(side+'-new-matrix-decomposition',entries['slab'],entries['conservative']+entries['mass'])
        entries['path']={k:clean(v.subs({eta:lam,sigma:lam/10})) for k,v in entries.items()}
        bound[side]=entries
        # End transverse waves have no density/thickness or face loss drivers.
        # These are new end-current argument joins, not a recomputed source.
        for sign in [1,-1]:
            mom=sign*kp;columns=S.T*B.subs(l,mom)
            mass=sp.Matrix(con['MATERIAL_CONSTRAINT_COEFFICIENTS']).T
            mc=bind(mass,nat,side+'-mass-constraint-'+str(sign),normal=mom)
            zero(side+'-mass-constraint-on-wave-'+str(sign),mc*columns,sp.zeros(1,2))
            ch=bind(slab['CHEMICAL_FIELD_ROW'],nat,side+'-chemical-row-'+str(sign),normal=mom)
            zero(side+'-chemical-on-wave-'+str(sign),ch*columns,sp.zeros(1,2))
            faces=nat['acoustic']['FACE_RECORDS'];require(len(faces)==2 and {int(f['ORIENTATION']) for f in faces}=={-1,1},'both physical faces')
            for fi,face in enumerate(faces):
                for col in range(2):
                    for key in ['OUTWARD_VELOCITY','AMPLITUDE']:
                        value=bind(face[key],nat,side+'-face-'+str(fi)+'-'+str(sign)+'-'+str(col)+'-'+key,normal=mom,amplitudes=columns[:,col])
                        zero(side+'-face-zero-'+str(fi)+'-'+str(sign)+'-'+str(col)+'-'+key,value,0)
    JL=bound['LEFT']['path']['slab'].subs(lam,0);JR=bound['RIGHT']['path']['slab']
    zero('native-baseline-current-join',JR.subs(lam,0),JL)
    Bminus=B.subs(l,-p);J0=JL.subs({km:p,kp:p});Jminus=JL.subs({km:-p,kp:-p})
    G0=clean(U.H*J0*U);Gref=clean(-Bminus.H*Jminus*Bminus)
    for name,G in [('incident',G0),('outward-reflected',Gref)]:
        zero(name+'-Hermiticity',G,G.H)
        minors=[clean(G[0,0]),clean(G.det())];J.emit(name+'-positive-weight',{'matrix':G,'leadingPrincipalMinors':minors,'signConvention':'outward reflection includes minus signed slab current' if name!='incident' else '+p incident toward defect'})
        require(all(x.is_positive is True for x in minors),'NONPOSITIVE_OR_UNKNOWN_PHYSICAL_CURRENT '+name)
    # Retain both momentum derivatives and both basis derivatives explicitly.
    point={km:p,kp:p,lam:0};explicit=JR.diff(lam).subs(point);leftLeg=delta*JR.diff(km).subs(point);rightLeg=delta*JR.diff(kp).subs(point)
    basisLeft=delta*Bp.H*J0*U;basisRight=delta*U.H*J0*Bp
    G1parts={'explicitNativeEnd':clean(U.H*explicit*U),'leftMomentum':clean(U.H*leftLeg*U),'rightMomentum':clean(U.H*rightLeg*U),'leftBasis':clean(basisLeft),'rightBasis':clean(basisRight)}
    G1=clean(sum(G1parts.values(),sp.zeros(2)));pR=p+lam*delta;BR=B.subs(l,pR);GR=BR.H*JR.subs({km:pR,kp:pR})*BR
    zero('right-end-flux-order0',GR.subs(lam,0),G0)
    zero('right-end-flux-chain-rule',GR.diff(lam).subs(lam,0),G1)
    zero('right-end-G1-Hermiticity',G1,G1.H)
    # Direct current contractions verify mass-rate correction cancellation on
    # the actual changing end basis, rather than discarding the correction.
    for side,basis,legs in [('LEFT',U,{km:p,kp:p}),('RIGHT',BR,{km:pR,kp:pR})]:
        correction=bound[side]['path']['mass'].subs(legs);zero(side+'-transverse-mass-current',clean(basis.H*correction*basis),sp.zeros(2))
    K1=clean(G1+T1.H*G0+G0*T1)
    J.emit('physical-linear-survival-operands',{'G0':G0,'Gref':Gref,'rightUnboundCurrent':JR,'leftUnboundCurrent':JL,'rightMomentum':pR,'rightBasis':BR,'G1Parts':G1parts,'G1':G1,'T1':T1,'K1':K1,'epsilonConvention':'All displayed weights multiply the same original epsilon_shape^2; harmonic 1/4 is already inside native current. No extra factor.','conditionalDeficit':'-lambda a^H K1 a/(a^H G0 a), only if the stationary expansion and other-channel asymptotics exist','reflectionSurvivalOrder2':'lambda^2 a^H R1^H Gref R1 a; unevaluated R1, not a complete order2 balance'})
    zero('K1-Hermiticity',K1,K1.H)
    direct=(sp.eye(2)+lam*T1).H*GR*(sp.eye(2)+lam*T1)
    zero('linear-survival-expansion',direct.diff(lam).subs(lam,0),K1)
    # A general infinitesimal outgoing coordinate change B_R -> B_R(I+lambda X)
    # shifts T1 by -X and G1 by X^H G0+G0 X. The physical K1 must be identical.
    xr=sp.symbols('flux_gauge_r0:4',real=True);xi=sp.symbols('flux_gauge_i0:4',real=True);X=sp.Matrix(2,2,lambda i,j:xr[2*i+j]+sp.I*xi[2*i+j])
    gaugeT=T1-X;gaugeG=G1+X.H*G0+G0*X;gaugeK=clean(gaugeG+gaugeT.H*G0+G0*gaugeT)
    gaugeBasis=BR*(sp.eye(2)+lam*X);gaugeCurrent=gaugeBasis.H*JR.subs({km:pR,kp:pR})*gaugeBasis
    zero('gauge-current-from-actual-native-map',gaugeCurrent.diff(lam).subs(lam,0),gaugeG)
    zero('gauge-field-order1',((gaugeBasis*(sp.eye(2)+lam*gaugeT)).diff(lam)).subs(lam,0),(BR*(sp.eye(2)+lam*T1)).diff(lam).subs(lam,0))
    J.emit('gauge-covariance-operands',{'outgoingCoordinateChange':X,'T1':gaugeT,'G1':gaugeG,'K1':gaugeK,'incomingCoordinatesUnchanged':True})
    zero('general-infinitesimal-outgoing-gauge',gaugeK,K1)
    shift=sp.Symbol('flux_profile_translation',real=True);phase=sp.exp(sp.I*(p-pR)*shift)
    shiftedT=T1+sp.diff(phase,lam).subs(lam,0)*sp.eye(2);shiftedK=clean(G1+shiftedT.H*G0+G0*shiftedT)
    J.emit('origin-covariance-operands',{'activeProfileTranslation':shift,'transmittedPhase':phase,'reflectedPhase':sp.exp(2*sp.I*p*shift),'translatedT1':shiftedT,'translatedK1':shiftedK,'derivation':'Translate the entire matched solution x->x-b and multiply by exp(ip b) to keep the incoming exp(ip x) coefficient fixed. This gives exp(i(p-pR)b) for transmission and exp(2ipb) for reflection. No source/profile/integral replay.'})
    zero('profile-translation-K1',shiftedK,K1)
    # Controls alter the actual complete current/survival assembly.
    J.nonzero('omit-G1-current-normalization',K1[0,0],(T1.H*G0+G0*T1)[0,0])
    J.nonzero('reverse-reflected-current-orientation',Gref[0,0],(-Gref)[0,0])
    noleft=clean(G1-G1parts['leftMomentum']+T1.H*G0+G0*T1)
    J.nonzero('omit-one-native-momentum-leg',K1[0,0],noleft[0,0])
    J.nonzero('change-T1-gauge-without-current',K1[0,0],(G1+(T1-sp.eye(2)).H*G0+G0*(T1-sp.eye(2)))[0,0])
    decision=all(v==0 for v in K1);J.emit('linear-survival-decision',{'K1':K1,'identicallyZero':decision,'ifNonzero':'STOP_ASSUMED_QUADRATIC_DEFICIT_ROUTE','ifZero':'One necessary obstruction removed, not zero leakage or complete power balance','allRealReceivingRegularity':False,'otherChannelAsymptoticsProved':False,'physicalHeldProfileWorkBalance':False})
    return {'executionStatus':'COMPLETED_NATIVE_TRANSVERSE_FLUX_PREREQUISITE','K1IdenticallyZero':decision,'quadraticDeficitRouteBlocked':not decision,'G0':G0,'Gref':Gref,'G1':G1,'K1':K1,'newControls':4,'restoredLiteralZeros':len(restored),'physicalWorkBalance':False,'allRealReceivingRegularity':False,'leakageFactor':None,'scientificAcceptance':False,'oldFunctionsCalled':False}
def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest);verify_invocation(args,gate,sys.argv)
    pins={**manifest['sourcePins'],**{r['path']:r['sha256'] for r in manifest['savedFiles'].values()},str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False);J=None;result={};code=1;start=time.monotonic()
    try:
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(manifest['helperSource']).read_text()),'unchanged-journal-and-containment','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        ns.update(sp=sp,Str=Str);J=ns['Journal'](args.out)
        result=J.stage('native-transverse-end-flux-prerequisite',{'manifestSha256':gate['manifestSha256'],'buildReviewSha256':gate['buildReviewRecordSha256']},lambda:scientific_work(manifest,J,ns['decode']));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result);sys.stderr.write(result['traceback'])
    finally:
        if (args.out/'saved-copy-index.json').exists():
            for v in json.loads((args.out/'saved-copy-index.json').read_text()).values():pins[str(args.out/v['path'])]=v['sha256']
        post={}
        for path,expected in pins.items():
            try:post[path]={'expected':expected,'actual':sha(path)}
            except OSError as e:post[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',post)
        if any(v['expected']!=v['actual'] for v in post.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-start,scientificAcceptance=False)
        save(args.out/'checks.json',J.encode(result) if J else result);sys.stdout.write((args.out/'checks.json').read_text())
    return code

if __name__=='__main__':sys.exit(main())
