#!/usr/bin/env python3
"""First-order incident-transverse source columns. No inverse or power calculation."""
import argparse, ast, hashlib, json, os, re, resource, shutil, sys, time, traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','containment','Journal','decode')
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
GRADES=('00','10','01')
VERDICT='CLEAR FOR THIS FIRST-ORDER TRANSVERSE SOURCE BUILD'

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
    require(g['status']=='READY_FOR_ONE_FIRST_ORDER_SOURCE_INSTRUMENT','one finite source gate')
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
    require(method['methodAssessed'] is True and method['literalVerdict']=='CLEAR FOR THIS FIRST-ORDER TRANSVERSE SOURCE METHOD','actual method')
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

def jet_spec(name):
    if name.startswith('grad_theta_'):name='theta_d'+name.rsplit('_',1)[1]
    match=re.fullmatch(r'(u_[123]|theta|e_W)((?:_t{1,2})?(?:_?d[123])*)',name)
    require(match is not None,'native wave jet '+name)
    channel,suffix=match.groups()
    return {'channel':channel,'timeOrder':2 if '_tt' in suffix else int('_t' in suffix),
            'spatialOrders':[len(re.findall('d'+str(i),suffix)) for i in range(1,4)]}

def scientific_work(manifest,J,decode):
    saved={};copies={};copyroot=J.out/'saved';copyroot.mkdir()
    for alias,record in manifest['savedFiles'].items():
        source=Path(record['path']);require(sha(source)==record['sha256'],'saved '+alias)
        target=copyroot/alias;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        require(sha(target)==record['sha256'],'copy '+alias)
        copies[alias]={'source':str(source),'path':str(target.relative_to(J.out)),'sha256':record['sha256']}
        save_index=J.out/'saved-copy-index.json'
        temporary=save_index.with_suffix('.new');temporary.write_text(json.dumps(copies,indent=2)+'\n');temporary.replace(save_index)
        if alias.endswith('.json'):saved[alias]=decode(json.loads(target.read_text()))
    def get(alias):return saved[alias]
    def symbol(expr,name):
        found=[s for s in expr.free_symbols if s.name==name];require(len(found)==1,'one actual symbol '+name);return found[0]
    def zero(name,left,right):J.zero(name,left,right)
    inherited_cache={}
    def inherited(prefix):
        if prefix in inherited_cache:return inherited_cache[prefix]
        operands=get(prefix+'-input.json');result=get(prefix+'-return.json')
        require(result['cancelled']==0,'inherited literal zero '+prefix)
        J.emit('inherited-'+prefix.replace('/','-'),{'operands':operands,'return':result,'functionsCalled':False})
        inherited_cache[prefix]=operands
        return operands
    physical=get('consumer/physical-input.json');context=get('local/extended-binding-context.json');binding=get('consumer/binding-context.json')
    require(context['physical']==physical==binding['physicalInput'],'same original physical input')
    require(context['frequencyOverride']=={'old':'1','actual':3},'held frequency distinct from original')
    values=context['numeric'];L=values['L_W'];omega=values['omega'];eps=context['epsilon'];eta,sigma=context['independentGrades']
    zero('held-omega',omega,sp.Integer(3));zero('native-L',L,sp.Integer(10))
    require(context['effectiveSpeedOnlyInPressure'] is True,'inherited raw local speed census')
    profiles=get('local/profile-jet-certificates.json');profile={r['name']:r for r in profiles}
    require(len(profile)==len(profiles),'unique original profile records')
    T=symbol(profile['w1_profile']['value'],'full_weak_tanh_variable')
    for rec in profiles:
        require(rec['L']==L and all(v['cancelled']==0 for v in rec['identities']),'inherited scaled profile proof')
    J.emit('restored-profile-rules',{'certificates':profiles,'rule':get('local/native-wave-profile-contract.json')['profileRule'],'newDerivativesComputed':False})
    # Existing profile jets carry L^n; only rename their original T variable.
    def profile_map(expr):
        return {s:profile[s.name]['value'] for s in expr.free_symbols if s.name in profile}
    def mapped(expr):
        result=expr.xreplace(profile_map(expr));require(result.free_symbols<={T},'complete profile map');return result
    def field_value(expr):
        replacements={}
        for t in expr.atoms(sp.tanh):
            x=symbol(t,'composition_x');zero('field-argument-'+str(len(J.artifacts)),t.args[0],x/L);replacements[t]=T
        result=expr.xreplace(replacements);require(result.free_symbols<={T},'stored field polynomial context');return result
    chart=get('incident/actual-source-chart.json');lift=get('incident/LEFT-invariant-lift.json')
    require(lift['actual']==lift['expected'],'saved exact lift arguments')
    S=chart['matrixWeakFromUniform'];N=symbol(lift['actual'],'uniformNormal')
    point=get('incident/point-01-LEFT.json');leg=next(r for r in point['momentumSigns'] if r['sign']==1)
    require(leg['mapping']['uniformFrequency']=='3' and leg['mapping']['uniformSoundSpeed']=='sqrt(6)/2' and leg['stratum']['stratum']=='EXACT_GRAZING','actual selected incident point')
    p=sp.sympify(leg['mapping']['uniformNormal']);h1=values['s11cdTangentialMomentum1'];h2=values['s11cdTangentialMomentum2'];k=(p,h1,h2)
    chartp=symbol(chart['weakCarrier'],'weak_end_p')
    for i in range(3):zero('chart-carrier-'+str(i),chart['weakCarrier'][i].subs(chartp,p),k[i])
    zero('chart-time',chart['timeRate'],-sp.I*omega)
    U=S*lift['actual'].subs(N,p);require(U.shape==(5,2),'actual two-column incident lift')
    J.emit('incident-columns',{'originalLift':lift,'chart':chart,'selectedPoint':leg,'columns':U,'fieldOrder':FIELDS,'fieldCoordinateNormalizationOnly':True})
    for c in range(2):
        zero('transverse-'+str(c),sum(k[i]*U[i,c] for i in range(3)),0)
        zero('theta-'+str(c),U[3,c],0);zero('eW-'+str(c),U[4,c],0)
    zero('incident-grazing',sp.Rational(6)-sum(v*v for v in k),0)
    zero('geometric-polarization',U[0,1],0)
    def multiplier(spec,column,derivative_sign=1):
        v=U[FIELDS.index(spec['channel']),column]*(-sp.I*omega)**spec['timeOrder']
        for ki,n in zip(k,spec['spatialOrders']):v*= (derivative_sign*sp.I*ki)**n
        return v
    def action(expr,column,derivative_sign=1):
        replacements=profile_map(expr)
        for a in expr.free_symbols:
            if a in replacements:continue
            if a.name in (eta.name,sigma.name):continue
            replacements[a]=multiplier(jet_spec(a.name),column,derivative_sign)
        result=sp.cancel(expr.xreplace(replacements))
        require(result.free_symbols<={T,eta,sigma},'complete source action');return result,replacements
    sources={};controls=[]
    cfac=(sp.Integer(30300)-101000*sp.I)/11119090000
    bracket=sp.I*(400+600*sum(v*v for v in k))*profile['m1_profile_d1']['value']+130*p*profile['m1_profile_d1d1']['value']+sp.I*(1900+2100*sum(v*v for v in k))*profile['w1_profile_d1']['value']+430*p*profile['w1_profile_d1d1']['value']
    for face in ['plus','minus']:
        for grade in GRADES:
            rec=get('sources/'+face+'-source-jets-'+grade+'.json');old=inherited('sources/'+face+'-source-jet-reconstruction-'+grade)
            require(old['left']==rec['source'] and old['right']==rec['reconstruction'],'original source/jet proof arguments')
            for col in range(2):
                name=face+'-source-'+grade+'-column-'+str(col);J.active=name
                J.emit(name+'-input',{'original':rec,'proof':old,'column':U[:,col],'profileCertificates':'restored-profile-rules','wave':k,'time':-sp.I*omega})
                result,replacements=action(rec['source'],col)
                terms=[]
                for jet in rec['jets']:
                    require(jet_spec(jet['spec']['name'])=={key:jet['spec'][key] for key in ['channel','timeOrder','spatialOrders']},'native jet selector')
                    term=field_value(jet['field'])*multiplier(jet['spec'],col)
                    terms.append({'atom':jet['atom'],'spec':jet['spec'],'coefficient':jet['coefficient'],'field':jet['field'],'multiplier':multiplier(jet['spec'],col),'term':term})
                J.emit(name+'-substitution',{'map':list(replacements.items()),'terms':terms,'fullReturn':result})
                zero(name+'-stored-jet-action',result,sum(t['term'] for t in terms))
                zero(name+'-hand-identity',result,cfac*U[0,col]*bracket if grade=='01' else 0)
                J.emit(name+'-return',{'value':result,'leftEndpoint':result.subs(T,-1),'rightEndpoint':result.subs(T,1),'newAppliedSource':True})
                sources[face,grade,col]=result
                if face=='plus' and grade=='01' and col==0:
                    # Actual term deletion and wrong incident derivative sign, through this same action.
                    live=next(t for t in terms if t['spec']['channel']=='u_1' and t['spec']['timeOrder']==0 and t['spec']['spatialOrders']==[0,0,0])
                    jet=next(j for j in rec['jets'] if j['atom']==live['atom'])
                    omitted,_=action(rec['source']-jet['coefficient']*jet['atom'],col)
                    reversed_wave,_=action(rec['source'],col,-1)
                    for label,wrong in [('omit-native-u1-source',omitted),('reverse-native-spatial-derivatives',reversed_wave)]:
                        movement=sp.cancel(wrong-result);J.emit(label+'-control',{'source':rec['source'],'column':U[:,col],'baseline':result,'mutant':wrong,'movement':movement,'physicalPoint':T,'testValue':movement.subs(T,0)})
                        require(movement.subs(T,0).is_zero is False,'actual responsive '+label);controls.append(label)
    # Apply the full already-published chemical expression; no old grade calculation.
    chemical=get('consumer/chemical-amplitude-grade-split.json')
    for col in range(2):
        value,sub=action(chemical['full'],col);baseline,sub0=action(chemical['zeroGrade'],col)
        J.emit('chemical-'+str(col)+'-input',{'original':chemical,'substitution':list(sub.items()),'column':U[:,col]})
        zero('chemical-baseline-'+str(col),baseline,0)
        zero('chemical-physical-source-'+str(col),value,-sp.I*sigma*U[0,col]*bracket/10100)
        J.emit('chemical-'+str(col)+'-return',{'fullApplied':value,'independentSigmaCoefficient':-sp.I*U[0,col]*bracket/10100,'faceOrientationNotSuppressed':True})
    for face in ['plus','minus']:
        native=get('consumer/'+face+'-source-input.json')
        v=native['nativeVelocity'];mp={a:(values['W_0'] if a.name=='W_0' else eps if a.name==eps.name else 0 if a.name=='e_W_t' else None) for a in v.free_symbols}
        require(all(b is not None for b in mp.values()),'native velocity leaves')
        J.emit(face+'-velocity-input',{'nativeSource':native,'map':list(mp.items()),'liftEW':list(U[4,:])})
        zero(face+'-held-profile-velocity',v.xreplace(mp),0)
        cc=native['chemicalCoefficient'];cm={a:(sp.S.Zero if a.name==eta.name else profile[a.name]['value'] if a.name in profile else None) for a in cc.free_symbols}
        require(all(v is not None for v in cm.values()),'actual baseline chemical normalization')
        Cchem=cc.xreplace(cm)
        for col in range(2):
            zero(face+'-chemical-to-normalized-source-'+str(col),sources[face,'01',col],Cchem*(-sp.I*U[0,col]*bracket/10100))
    # Apply saved local cells; time and edge derivatives are already in coefficients.
    cells=get('local/all-local-cells.json');require(len(cells)==400,'all original local cells')
    keys={(r['row'],r['fieldColumn'],r['xOrder'],tuple(r['grade'])) for r in cells};require(len(keys)==400,'unique original local cells')
    local={};local_terms=[]
    for row in ROWS:
        for grade in GRADES:
            selected=[v for v in cells if v['row']==row and v['grade']==[int(x) for x in grade]]
            old=inherited('local/'+row+'-new-physical-row-grade-'+grade)
            for col in range(2):
                recs=[{'cell':v,'incidentMultiplier':(sp.I*p)**v['xOrder']*U[v['fieldColumn'],col],
                       'applied':v['coefficient']*(sp.I*p)**v['xOrder']*U[v['fieldColumn'],col]} for v in selected]
                name=row+'-'+grade+'-column-'+str(col);J.emit(name+'-input',{'cells':recs,'inheritedProof':old})
                total=sp.cancel(sum(v['applied'] for v in recs));local[row,grade,col]=total
                for side in ['left','right']:
                    expr=old[side];mapping={}
                    for s in expr.free_symbols:
                        if s==T:continue
                        mt=re.fullmatch(r'local_field_(\d)_jet_(\d+)',s.name);require(mt is not None,'original local row jet')
                        mapping[s]=(sp.I*p)**int(mt[2])*U[int(mt[1]),col]
                    zero(name+'-applied-'+side,total,expr.xreplace(mapping))
                endpoints={side:sp.cancel(sum(v['cell'][side+'Endpoint']*v['incidentMultiplier'] for v in recs)) for side in ['left','right']}
                J.emit(name+'-return',{'operatorAction':total,'sourceForce':-total,'endpoints':endpoints,'nativeEpsilonPower':[v['cell']['epsilonPower'] for v in recs]})
                if grade=='00':zero(name+'-baseline',total,0)
                local_terms.append({'row':row,'grade':grade,'column':col,'action':total,'endpoints':endpoints})
    # Whole flat response and both lab-normal maps, receiving l independent of p.
    census=get('consumer/reference-response-census.json');routes=get('pressure/route-inspection.json')['routes']
    for route in routes:
        origin=manifest['routeArrays'][route['sourcePath']]
        require(manifest['savedFiles'][origin]['sha256']==route['sourceSha256'],'original route-array hash')
        require(get(origin)[int(route['pointer'].removeprefix('/'))]==route['completeAddress'],'original exact route selector')
    q=sp.Symbol('first_order_outgoing_depth');ell=sp.Symbol('first_order_receiving_momentum',real=True)
    beta=omega/(10-omega*sp.I);R=sp.Rational(3,10)/(q+beta)
    flat=census['flat'];fm={s:(omega if s.name=='reference_unrestricted_frequency' else q if s.name=='reference_qi' else None) for s in flat.free_symbols}
    require(all(v is not None for v in fm.values()),'original flat leaves');zero('flat-whole-response',flat.xreplace(fm),R)
    J.emit('receiving-sheet',{'momentum':ell,'depth':q,'depthSquared':6-ell**2-h1**2-h2**2,'outgoingConvention':'q>=0 real inside; q=+i sqrt(-q^2) outside','flatResponse':R,'beta':beta,'finiteAtGrazing':R.subs(q,0),'noForcedInverseClaim':True})
    factors={}
    for route in routes:
        a=route['completeAddress'];proof=a['fullFactorProof']['proof'];op=get('flat-proofs/'+proof+'-operands.json')
        normal=inherited('flat-proofs/'+proof+'-normal-source-join');full=inherited('flat-proofs/'+proof+'-full-mapped-residual')
        require(a['responseMap']['flatSupport']=='k=l' and a['responseInputDepth']=='q(l)' and a['responseOutputDepth']=='q(l)','actual flat support')
        require(all(op['actualResponseMap'][key]==a['responseMap'][key] for key in ['original','mapped','map','flatSupport','frequency']) and op['savedNormal']==a['normalOriginal'],'actual response/normal argument join')
        mapped_factor=a['responseCoefficient']*a['normalMultiplier']
        calls=mapped_factor.atoms(sp.Function);require(len(calls)==1,'one receiving depth call')
        call=next(iter(calls));require(call.func.__name__=='common_outgoing_q' and str(call.args[0])=='composition_l','actual output depth')
        factor=mapped_factor.xreplace({call:q});expected=R*(1 if route['slot']=='pressure' else (1 if route['face']=='plus' else -1)*sp.I*q)
        zero(route['face']+'-'+route['slot']+'-flat-factor',factor,expected)
        factors[route['face'],route['slot']]=factor
        J.emit(route['face']+'-'+route['slot']+'-flat-arguments',{'address':a,'proofOperands':op,'normalProof':normal,'fullProof':full,'factor':factor})
    # Tags are the forward transform of the entire source envelope at l-p, not values.
    pressure_rows=[]
    for row in ROWS:
        inherited('pressure/'+row+'-affine-pressure-reconstruction')
        for col in range(2):
            pieces=[]
            for face in ['plus','minus']:
                for slot,prefix in [('pressure','delta_p_'),('normal','d_w_delta_p_')]:
                    split=get('pressure/'+row+'-'+prefix+face+'-split.json')
                    C={g:mapped(sp.cancel(split['retained'][str(tuple(map(int,g)))]/eps)) for g in GRADES}
                    require(not C['00'].free_symbols,'constant baseline native consumer')
                    # Drop response10/01 only on the already verified identical zero function.
                    require(sources[face,'00',col]==0 and sources[face,'10',col]==0,'whole-function zero before response')
                    tag=sp.Symbol('source01_'+face+'_column_'+str(col)+'_hat_l_minus_p')
                    applied=0 if sources[face,'01',col]==0 else C['00']*factors[face,slot]*tag
                    pieces.append({'face':face,'slot':slot,'split':split,'epsilon':eps,'consumers':C,'responseCensus':{'height':census['heightCoefficient'],'slope':census['slopeCoefficient']},'source00':sources[face,'00',col],'source10':sources[face,'10',col],'source01Envelope':sources[face,'01',col],'transformTag':tag,'transformConvention':'integral exp(-i(l-p)x) envelope(x) dx; inverse dl/(2pi); not evaluated','flatNormalFactor':factors[face,slot],'pressure01':applied})
            total=sum(v['pressure01'] for v in pieces)
            pressure_rows.append({'row':row,'column':col,'pressure10':0,'pressure01':total,'pieces':pieces})
    J.emit('first-order-pressure-assembly',pressure_rows)
    for col in range(2):
        vector10=sp.Matrix([-local[row,'10',col] for row in ROWS]);vector01=sp.Matrix([-local[row,'01',col] for row in ROWS]);force=vector10+vector01/10
        J.emit('full-local-force-column-'+str(col),{'etaForce':vector10,'sigmaForce':vector01,'physicalPathForce':force,'path':'eta=lambda,sigma=lambda/10 AFTER independent grades','pressureAddition':'negative C00 R00 S01 transform /10; see pressure assembly, not multiplied into a local profile','receivingDotDiagnostic':sp.cancel(ell*force[0]+h1*force[1]+h2*force[2]),'receivingMomentum':ell,'incidentMomentum':p,'localLeft':force.subs(T,-1),'localRight':force.subs(T,1),'diagnosticIsPhysicalLossProjection':False})
    return {'executionStatus':'COMPLETED_FIRST_ORDER_SOURCE_SCREEN','incidentColumns':2,'nativeFaces':2,'localCellsRestored':400,'appliedGrades':GRADES,'newControls':controls,'completePhysicalFaceMaps':False,'matchedEndForcing':False,'receivingInverse':False,'grazingRegularity':False,'powerComputed':False,'leakageFactor':None,'oldFunctionsCalled':False,'scienceAcceptance':False}

def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest);verify_invocation(args,gate,sys.argv)
    pins={**manifest['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False);J=None;result={};code=1;start=time.monotonic()
    try:
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(manifest['helperSource']).read_text()),'unchanged-journal-and-containment','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        ns.update(sp=sp,Str=Str);J=ns['Journal'](args.out)
        result=J.stage('first-order-source-screen',{'manifestSha256':gate['manifestSha256'],'buildReviewSha256':gate['buildReviewRecordSha256']},lambda:scientific_work(manifest,J,ns['decode']));code=0
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
