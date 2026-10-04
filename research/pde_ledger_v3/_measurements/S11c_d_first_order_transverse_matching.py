#!/usr/bin/env python3
"""Finite matched transverse amplitudes in saved coordinates. No flux or leakage calculation."""
import argparse, ast, hashlib, json, os, re, resource, shutil, sys, time, traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','containment','Journal','decode')
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
GRADES=('00','10','01')
VERDICT='CLEAR FOR THIS FINITE MATCHED TRANSVERSE AMPLITUDE BUILD'

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
    require(g['status']=='READY_FOR_ONE_FIRST_ORDER_TRANSVERSE_MATCHING_INSTRUMENT','one finite receiving prerequisite gate')
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

def jet_spec(name):
    if name.startswith('grad_theta_'):name='theta_d'+name.rsplit('_',1)[1]
    match=re.fullmatch(r'(u_[123]|theta|e_W)((?:_t{1,2})?(?:_?d[123])*)',name)
    require(match is not None,'native wave jet '+name)
    channel,suffix=match.groups()
    return {'channel':channel,'timeOrder':2 if '_tt' in suffix else int('_t' in suffix),
            'spatialOrders':[len(re.findall('d'+str(i),suffix)) for i in range(1,4)]}


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
    def clean(expr):return expr.applyfunc(sp.cancel) if isinstance(expr,sp.MatrixBase) else sp.cancel(expr)
    def zero(name,left,right):
        if isinstance(left,sp.MatrixBase):
            require(isinstance(right,sp.MatrixBase) and left.shape==right.shape,'matrix shape '+name)
            J.emit(name+'-matrix-operands',{'left':left,'right':right})
            for i in range(left.rows):
                for j in range(left.cols):J.zero(name+'-'+str(i)+'-'+str(j),left[i,j],right[i,j])
        else:J.zero(name,left,right)
    def symbol(expr,name):
        items=[v for v in expr.free_symbols if v.name==name];require(len(items)==1,'actual symbol '+name);return items[0]
    old=get('source/incident-columns.json');U=old['columns'];p=sp.sympify(old['selectedPoint']['mapping']['uniformNormal'])
    context=get('source/extended-binding-context.json');L=context['numeric']['L_W'];h1=context['numeric']['s11cdTangentialMomentum1'];h2=context['numeric']['s11cdTangentialMomentum2'];H=sp.sqrt(h1*h1+h2*h2)
    require(context['frequencyOverride']=={'old':'1','actual':3} and context['numeric']['omega']==3 and L==10,'actual physical scale/frequency')
    require(p.is_positive is True and H.is_positive is True and L.is_positive is True,'positive outgoing and profile scales')
    end=get('receiving/end-matching-prerequisite.json');Bmat=end['B'];DB=end['DB'];l=symbol(Bmat,'receiving_block_l');delta=end['deltaP'];Bp=end['BprimeAtIncident'];DBp=end['DBprimeAtIncident']
    # The old exact checks are restored, not recalculated.
    restored=[]
    for prefix,rows,cols in [('FULL-offwave-source-correspondence',5,5),('transverse-block',2,2),('transverse-to-scalars-coupling',3,2),('scalars-to-transverse-coupling',2,3),('full-five-row-matched-end',5,2),('end-RHS-minus-operator-sign',5,2)]:
        for i in range(rows):
            for j in range(cols):
                name=prefix+'-'+str(i)+'-'+str(j);op=get('receiving/'+name+'-input.json');ret=get('receiving/'+name+'-return.json');require(ret['cancelled']==0,'inherited '+name);restored.append({'name':name,'operands':op,'return':ret})
    transformed=get('receiving/complete-transformed-operator.json');require(transformed['nonzeroCouplingPolicy']=='STOP_NO_AUTOMATIC_SCHUR_INVERSE','actual accepted block scope')
    require(get('receiving/transverse-block-matrix-operands.json')['left']==transformed['matrix'][:2,:2],'actual inherited transverse block operand')
    require(get('receiving/transverse-block-matrix-operands.json')['right']==sp.eye(2)*sp.Rational(3,2)*(l*l-p*p),'actual p and scalar propagation normalization')
    require(get('receiving/transverse-to-scalars-coupling-matrix-operands.json')['left']==transformed['matrix'][2:,:2] and get('receiving/scalars-to-transverse-coupling-matrix-operands.json')['left']==transformed['matrix'][:2,2:],'actual inherited coupling operands')
    require(get('receiving/end-RHS-minus-operator-sign-matrix-operands.json')['right']==end['rightForcing'],'actual end-sign operand')
    require(old['columns']==get('receiving/incident-chart-match-matrix-operands.json')['right'],'actual incident match operand')
    J.emit('inherited-receiving-returns',{'records':restored,'end':end,'transformed':transformed,'functionsCalled':False})
    source_assembly=get('source/first-order-pressure-assembly.json')
    for rec in source_assembly:
        if rec['row'] in ['U0','U1','U2']:
            require(all(all(v==0 for v in part['consumers'].values()) for part in rec['pieces']),'actual pressure displacement consumers vanish at every retained grade')
    J.emit('pressure-transverse-absence',{'actualCompleteAssembly':source_assembly,'claim':'All displacement-row pressure consumer grades are literal zero. No pressure term is inserted into this transverse forcing.'})
    forces=[get('source/full-local-force-column-'+str(i)+'.json') for i in range(2)]
    require(forces==end['savedForces'],'exact completed force arguments')
    T=symbol(forces[0]['etaForce'],'full_weak_tanh_variable');w=(1+T)/2;m=(1-T*T)/3;A=9*w-6*m
    Feta=sp.Matrix.hstack(*(f['etaForce'] for f in forces));Fsigma=sp.Matrix.hstack(*(f['sigmaForce'] for f in forces));F=sp.Matrix.hstack(*(f['physicalPathForce'] for f in forces))
    zero('new-height-scalar-step-identity',Feta,A*U)
    zero('new-complete-path-identity',F,Feta+Fsigma/10)
    require(all(f['path']=='eta=lambda,sigma=lambda/10 AFTER independent grades' for f in forces),'actual independent path')
    # These are NEW polynomial moment calculations on the saved forces.
    # H(x) is split at the original tanh centre x=0; no value of Ahat(0).
    fminus=F;fplus=F-9*U;intervals=[('left',fminus,-1,0),('right',fplus,0,1)];moment=sp.zeros(5,2);moment_records=[]
    for side,poly,a,b in intervals:
        factor=1+T if side=='left' else 1-T;remaining=1-T if side=='left' else 1+T
        for i in range(5):
            for j in range(2):
                # Cancel only the endpoint approached on this half-line.
                numerator=sp.Poly(poly[i,j],T);quot,rem=sp.div(numerator,sp.Poly(factor,T))
                J.emit(side+'-moment-division-'+str(i)+'-'+str(j),{'force':poly[i,j],'endpointFactor':factor,'quotient':quot.as_expr(),'remainder':rem.as_expr(),'interval':[a,b],'jacobian':L/(1-T*T)})
                zero(side+'-moment-endpoint-rem-'+str(i)+'-'+str(j),rem.as_expr(),0)
                qpoly,rpoly=sp.div(quot,sp.Poly(remaining,T));rval=rpoly.as_expr();require(not rval.has(T),'constant simple-fraction remainder')
                polynomial_primitive=sum(c*T**(power[0]+1)/sp.Integer(power[0]+1) for power,c in qpoly.terms())
                log_primitive=(-rval*sp.log(1-T) if side=='left' else rval*sp.log(1+T))
                primitive=L*(polynomial_primitive+log_primitive)
                # New elementary primitive certificate, not an old profile derivative.
                zero(side+'-moment-primitive-'+str(i)+'-'+str(j),sp.diff(primitive,T),L*poly[i,j]/(1-T*T))
                value=sp.expand(primitive.subs(T,b)-primitive.subs(T,a));require(value.is_finite is True,'finite new half-line moment')
                moment[i,j]+=value;moment_records.append({'side':side,'indices':[i,j],'force':poly[i,j],'primitive':primitive,'endpoints':[a,b],'endpointValues':[primitive.subs(T,a),primitive.subs(T,b)],'value':value})
    moment=clean(moment);J.emit('new-localized-forward-moment',{'fullForce':F,'rightStep':9*U,'halfLineRecords':moment_records,'moment':moment,'transformConvention':'Fhat(Q)=integral exp(-iQx) F(x) dx; inverse measure dl/(2pi)','momentumTransfer':0,'origin':0,'profileLength':L})
    # Elementary outgoing scalar Green jump and exact scalar-step reference.
    x=sp.Symbol('matching_x',real=True);y=sp.Symbol('matching_y',real=True);eps=sp.Symbol('matching_abel_epsilon',positive=True)
    gplus=sp.I*sp.exp(sp.I*p*x)/(3*p);gminus=sp.I*sp.exp(-sp.I*p*x)/(3*p)
    zero('outgoing-green-right-ODE',-sp.Rational(3,2)*(sp.diff(gplus,x,2)+p*p*gplus),0)
    zero('outgoing-green-left-ODE',-sp.Rational(3,2)*(sp.diff(gminus,x,2)+p*p*gminus),0)
    zero('outgoing-green-continuity',gplus.subs(x,0),gminus.subs(x,0))
    zero('outgoing-green-delta-jump',-sp.Rational(3,2)*(sp.diff(gplus,x).subs(x,0)-sp.diff(gminus,x).subs(x,0)),1)
    abelPrimitive=sp.exp((2*sp.I*p-eps)*y)/(2*sp.I*p-eps)
    zero('Abel-half-line-primitive',sp.diff(abelPrimitive,y),sp.exp((2*sp.I*p-eps)*y))
    abelValue=-1/(2*sp.I*p);scalar_step=-delta/(2*p)
    zero('scalar-step-reflection-normalization',9*sp.I*abelValue/(3*p),scalar_step)
    lam=sp.Symbol('matching_lambda',real=True);pR=p+lam*delta
    exact_r=(p-pR)/(p+pR);exact_t=2*p/(p+pR)
    zero('scalar-exact-step-R1',sp.diff(exact_r,lam).subs(lam,0),scalar_step)
    zero('scalar-exact-step-T1',sp.diff(exact_t,lam).subs(lam,0),scalar_step)
    J.emit('scalar-outgoing-matching',{'greenPositive':gplus,'greenNegative':gminus,'AbelPrimitive':abelPrimitive,'AbelUpperEndpointForEpsilonPositive':0,'AbelLowerEndpoint':abelPrimitive.subs(y,0),'AbelLimit':abelValue,'rightSecularCoefficient':sp.I*delta,'rightFinite':scalar_step,'leftReflection':scalar_step,'exactStepR':exact_r,'exactStepT':exact_t,'analyticStatement':'For eps>0 the primitive tends to0 at infinity. Its lower endpoint gives the Abel integral; eps->0 is taken only after integration. No product of delta(Q) with the outgoing pole.'})
    # The rational transverse projector introduces only exponentially decaying
    # coordinate kernels; it has no acoustic q branch or real pole of its own.
    projector=clean(Bmat*DB);projectedU=clean(projector*U);dprojected=clean(projectedU.diff(l).subs(l,p));DB0=DB.subs(l,p);DBminus=DB.subs(l,-p)
    zero('new-projector-step-polarization-split',dprojected,Bp+Bmat.subs(l,p)*DBp*U)
    helmholtz=sp.exp(-H*sp.Abs(x))/(2*H)
    # A polynomial numerator over l^2+H^2 certifies the only coordinate poles.
    pole_numerators=[]
    for v in projector:
        n=sp.cancel((l*l+H*H)*v);require(sp.denom(n).has(l) is False,'no other coordinate pole');sp.Poly(n,l);pole_numerators.append(n)
    J.emit('transverse-projection-asymptotic-argument',{'projector':projector,'coordinatePoleNumerators':pole_numerators,'coordinateDenominator':l*l+H*H,'coordinateKernel':helmholtz,'H':H,'sourceTailRate':2/L,'sourcePolynomials':{'left':fminus,'right':fplus},'outgoingScalarGreen':{'positive':gplus,'negative':gminus},'argument':'Each half-line remainder has a removable approached-end factor, hence exponential decay at rate2/L away from origin. Projector entries are polynomial numerators over l^2+H^2, whose inverse kernel is exp(-H|x|)/(2H); differentiated coordinate kernels give decaying tails plus local distributions. The incident step generates the retained secular/finite terms. Remaining coordinate tails vanish, so only the real transverse poles contribute propagating transverse end amplitudes. Exact block separation and displacement pressure zeros remove q from this projected channel. This is an assessed analytic argument, not machine measure theory or a statement about other channels.','fullCoupledFieldAsymptoticsProved':False})
    stepT=scalar_step*sp.eye(2)+delta*DBp*U
    localT=clean(sp.I*DB0*moment/(3*p));T1=clean(stepT+localT)
    rightFinite=clean(scalar_step*U+delta*dprojected+sp.I*projector.subs(l,p)*moment/(3*p))
    zero('right-field-finite-matching',rightFinite,delta*Bp+U*T1)
    # Reflection retains an exact convergent integral, not a numerical value.
    source_local_minus=DBminus*fminus;source_local_plus=DBminus*fplus
    reflection_integrals=sp.Matrix(2,2,lambda i,j:sp.Symbol('R_local_integral_'+str(i)+'_'+str(j)))
    Rstep=clean(scalar_step*DBminus*U);R1=Rstep+sp.I*reflection_integrals/(3*p)
    J.emit('matched-transverse-amplitudes',{'p':p,'deltaP':delta,'B':Bmat,'DB':DB,'Bprime':Bp,'DBprime':DBp,'moment':moment,'stepT1':stepT,'localizedT1':localT,'T1':T1,'rightFiniteField':rightFinite,'rightSecularField':sp.I*delta*U,'reflectionStep':Rstep,'reflectionLocalIntegralTags':reflection_integrals,'reflectionIntegrands':{'left':source_local_minus,'right':source_local_plus,'phase':sp.exp(2*sp.I*p*x),'profileArgument':sp.tanh(x/L),'variable':T,'split':0,'transfer':-2*p},'R1Formal':R1,'reflectionIntegralEvaluated':False,'amplitudesInOriginalIncidentCoordinates':True,'arbitraryIncidentColumnVector':True,'physicalFluxNormalizationApplied':False})
    J.nonzero('omit-dual-momentum-derivative-control',T1[1,0],(T1-delta*DBp*U)[1,0])
    J.nonzero('reverse-scalar-outgoing-control',scalar_step,-scalar_step)
    return {'executionStatus':'COMPLETED_FINITE_MATCHED_TRANSVERSE_AMPLITUDE_PREREQUISITE','T1ComputedInSavedCoordinates':True,'R1ExactIntegralRepresentation':True,'reflectionIntegralEvaluated':False,'physicalFluxNormalizationApplied':False,'G0Computed':False,'GrefComputed':False,'G1Computed':False,'K1Computed':False,'allRealCoupledRegularity':False,'physicalWorkBalance':False,'leakageFactor':None,'scienceAcceptance':False,'oldFunctionsCalled':False,'newControls':2}
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
        result=J.stage('matched-transverse-amplitude-prerequisite',{'manifestSha256':gate['manifestSha256'],'buildReviewSha256':gate['buildReviewRecordSha256']},lambda:scientific_work(manifest,J,ns['decode']));code=0
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
