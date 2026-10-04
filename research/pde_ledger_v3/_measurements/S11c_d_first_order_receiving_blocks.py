#!/usr/bin/env python3
"""Finite receiving-block, end-sign and affine-face prerequisite. No scattering/power calculation."""
import argparse, ast, hashlib, json, os, re, resource, shutil, sys, time, traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','containment','Journal','decode')
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
GRADES=('00','10','01')
VERDICT='CLEAR FOR THIS FIRST-ORDER RECEIVING-BLOCK AND FACE BUILD'

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
    require(g['status']=='READY_FOR_ONE_FIRST_ORDER_RECEIVING_BLOCKS_INSTRUMENT','one finite receiving prerequisite gate')
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
    # Complete saved files are copied; only operands used below are decoded.
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
    def sym(expr,name):
        found=[a for a in expr.free_symbols if a.name==name];require(len(found)==1,'actual symbol '+name);return found[0]
    def clean(expr):
        return expr.applyfunc(sp.cancel) if isinstance(expr,sp.MatrixBase) else sp.cancel(expr)
    def zeros(name,left,right):
        if isinstance(left,sp.MatrixBase):
            require(isinstance(right,sp.MatrixBase) and left.shape==right.shape,'matrix operands '+name)
            J.emit(name+'-matrix-operands',{'left':left,'right':right})
            for i in range(left.rows):
                for j in range(left.cols):J.zero(name+'-'+str(i)+'-'+str(j),left[i,j],right[i,j])
        else:J.zero(name,left,right)
    def substitute(name,expr,mapping,allowed):
        J.emit(name+'-substitution-input',{'expression':expr,'mapping':list(mapping.items()),'allowedSymbols':sorted(allowed,key=str)})
        result=expr.xreplace(mapping);J.emit(name+'-substitution-return',{'value':result,'remaining':sorted(result.free_symbols,key=str)})
        require(result.free_symbols<=set(allowed),'unbound '+name);return result
    def nonzero(name,baseline,mutant):
        J.nonzero(name,baseline,mutant)
    context=get('source-input/local/extended-binding-context.json');values=context['numeric'];eps=context['epsilon']
    physical=get('source-input/consumer/physical-input.json')
    require(context['physical']==physical==get('source-input/consumer/binding-context.json')['physicalInput'],'same original physical inputs')
    require(context['frequencyOverride']=={'old':'1','actual':3},'same held omega3')
    omega,L=values['omega'],values['L_W'];zeros('held-frequency',omega,sp.Integer(3));zeros('profile-length',L,sp.Integer(10))
    old=get('source-result/incident-columns.json');U=old['columns'];S=old['chart']['matrixWeakFromUniform']
    require(U.shape==(5,2) and S.shape==(5,5),'original two-column/chart dimensions')
    p=sp.sympify(old['selectedPoint']['mapping']['uniformNormal']);h1,h2=(values[n] for n in ['s11cdTangentialMomentum1','s11cdTangentialMomentum2'])
    l=sp.Symbol('receiving_block_l',real=True);q=sp.Symbol('receiving_block_q');lam=sp.Symbol('receiving_small_lambda',real=True)
    psi=sp.Matrix(sp.symbols('receiving_field_0:5'));k=(l,h1,h2)
    sheet=get('source-result/receiving-sheet.json');oldq=sheet['depth'];oldl=sheet['momentum']
    zeros('actual-sheet',sheet['depthSquared'].xreplace({oldl:l}),p**2-l**2)
    R=sheet['flatResponse'].xreplace({oldq:q});beta=sheet['beta'];rho=values['rho_m']
    zeros('flat-beta-memory',beta,omega*values['Lambda_A_0']/(rho*(1-sp.I*omega*values['tau_A'])))
    J.emit('outgoing-domain',{'qSquared':p**2-l**2,'realInterior':'q=positive sqrt(p^2-l^2)','realExterior':'q=+i sqrt(l^2-p^2)','threshold':q,'beta':beta,'ReBeta':sp.re(beta),'extension':'cancel the original rational expression for q != 0, then q=0; no change of physical frequency','fullReceivingRegularity':False})
    require(sp.re(beta).is_positive is True,'nonzero outgoing closure denominator')
    def wave(expr,label):
        mp={}
        for a in expr.free_symbols:
            spec=jet_spec(a.name);v=psi[FIELDS.index(spec['channel'])]*(-sp.I*omega)**spec['timeOrder']
            for momentum,n in zip(k,spec['spatialOrders']):v*= (sp.I*momentum)**n
            mp[a]=v
        return substitute(label,expr,mp,{l,*psi})
    def row(expr,label):
        coefficients=sp.Matrix([[sp.diff(expr,a) for a in psi]])
        zeros(label+'-linear',expr,(coefficients*psi)[0]);return coefficients
    # New receiving argument l, using the saved cell coefficients and no producer.
    cells=get('source-input/local/all-local-cells.json');require(len(cells)==400,'complete original cell inventory')
    local=sp.zeros(5);entries=[];seen=set()
    for cell in cells:
        key=(cell['row'],cell['fieldColumn'],cell['xOrder'],tuple(cell['grade']))
        require(key not in seen,'unique saved cell');seen.add(key)
        if cell['grade']!=[0,0]:continue
        require(not cell['coefficient'].free_symbols,'constant baseline coefficient')
        require(all(proof['cancelled']==0 for proof in cell['identities']),'inherited cell returns')
        require(cell['polynomial']['original']==cell['coefficient'],'actual cell polynomial operand')
        value=cell['coefficient']*(sp.I*l)**cell['xOrder'];local[ROWS.index(cell['row']),cell['fieldColumn']]+=value
        entries.append({'cell':cell,'newReceivingMultiplier':(sp.I*l)**cell['xOrder'],'term':value})
    require(len(entries)==100,'all baseline receiving cells');J.emit('local-receiving-assembly',{'entries':entries,'matrix':local,'oldFunctionsCalled':False})
    source_rows={}
    for face in ['plus','minus']:
        rec=get('source-input/sources/'+face+'-source-jets-00.json')
        prior=get('source-input/sources/'+face+'-source-jet-reconstruction-00-return.json')
        op=get('source-input/sources/'+face+'-source-jet-reconstruction-00-input.json')
        require(prior['cancelled']==0 and op['left']==rec['source'] and op['right']==rec['reconstruction'],'saved source proof arguments')
        J.emit(face+'-inherited-source00',{'record':rec,'operands':op,'return':prior})
        source_rows[face]=row(wave(rec['source'],face+'-receiving-source00'),face+'-source00')
    assembly=get('source-result/first-order-pressure-assembly.json');pressure=sp.zeros(5);pressure_entries=[]
    for name in ROWS:
        records=[v for v in assembly if v['row']==name and v['column']==0];require(len(records)==1,'one complete pressure row')
        pieces=records[0]['pieces'];require(len(pieces)==4,'all pressure slots')
        for piece in pieces:
            face,slot=piece['face'],piece['slot'];factor=piece['flatNormalFactor'].xreplace({oldq:q});consumer=piece['consumers']['00']
            require(not consumer.free_symbols,'constant baseline consumer')
            original=get('source-result/'+face+'-'+slot+'-flat-arguments.json')
            require(original['factor']==piece['flatNormalFactor'],'saved factor argument identity')
            pressure[ROWS.index(name),:]+=consumer*factor*source_rows[face]
            pressure_entries.append({'row':name,'originalPiece':piece,'factorArguments':original,'newFactor':factor,'receivingSourceRow':source_rows[face],'contribution':consumer*factor*source_rows[face]})
    require(len(pressure_entries)==20,'complete five-row/four-slot assembly')
    J.emit('pressure-receiving-assembly',{'entries':pressure_entries,'matrix':pressure,'sourceMomentum':l,'normalOutputDepth':q,'flatSupport':'input=output=l','wholeFlatResponseOnce':True})
    Pweak=clean(local+pressure);uniform=get('source-input/incident/LEFT-invariant-P.json')
    require(uniform['actual']==uniform['expected'],'original uniform operand identity')
    mp={sym(uniform['actual'],'uniformNormal'):l,sym(uniform['actual'],'uniformFrequency'):omega,sym(uniform['actual'],'uniformPhysicalDepth'):q}
    native=substitute('uniform-native-receiving',uniform['actual'],mp,{l,q});nativeWeak=S*native*S.T
    zeros('FULL-offwave-source-correspondence',Pweak,nativeWeak)
    # Independently bind the complete original native source, not just selected columns.
    rawends={}
    for side in ['LEFT','RIGHT']:
        saved=get('source-input/incident/LEFT-raw-source-binding.json' if side=='LEFT' else 'receiving/RIGHT-raw-source-binding.json')
        mp=dict(saved['map']['mappingPairs']);mp={a:b for a,b in mp.items() if a.name not in ['eta_bg','sigma_W']}
        raw=saved['nativeSource'].xreplace(mp);remaining={}
        for a in raw.free_symbols:
            if a.name=='weak_end_p':remaining[a]=l
            elif a.name=='weak_end_q':remaining[a]=q
            elif a.name=='weak_end_cs':remaining[a]=sp.sqrt(6)/2
            elif a.name=='eta_bg':remaining[a]=lam
            elif a.name=='sigma_W':remaining[a]=lam/10
            else:raise ValueError('original end mapping leaf '+str(a))
        J.emit(side+'-raw-source-arguments',{'saved':saved,'originalMapping':list(mp.items()),'newArgumentMapping':list(remaining.items()),'originalFiniteOriginNotUsed':True})
        rawends[side]=clean(S*substitute(side+'-raw-end',raw,remaining,{l,q,lam})*S.T)
    zeros('LEFT-raw-vs-receiving',rawends['LEFT'],Pweak)
    zeros('RIGHT-baseline-vs-receiving',rawends['RIGHT'].subs(lam,0),Pweak)
    H2=h1*h1+h2*h2;K2=l*l+H2;t=sp.Matrix([0,h2,-h1]);b=sp.Matrix([H2,-l*h1,-l*h2]);kv=sp.Matrix(k)
    chart=sp.diag(1,1,1,1,1);chart[:3,0]=t;chart[:3,1]=b;chart[:3,2]=kv
    dual=sp.eye(5);dual[0,:3]=t.T/H2;dual[1,:3]=b.T/(H2*K2);dual[2,:3]=kv.T/K2
    zeros('geometric-dual',dual*chart,sp.eye(5));zeros('geometric-H2',H2,sp.Rational(1,20))
    transformed=clean(dual*Pweak*chart);DT=sp.Rational(3,2)*(l*l-p*p)
    J.emit('complete-transformed-operator',{'physicalMatrix':Pweak,'chart':chart,'dual':dual,'matrix':transformed,'DTCandidate':DT,'bothCouplingDirectionsRetained':True,'nonzeroCouplingPolicy':'STOP_NO_AUTOMATIC_SCHUR_INVERSE'})
    zeros('transverse-block',transformed[:2,:2],sp.eye(2)*DT)
    zeros('transverse-to-scalars-coupling',transformed[2:,:2],sp.zeros(3,2))
    zeros('scalars-to-transverse-coupling',transformed[:2,2:],sp.zeros(2,3))
    block=transformed[2:,2:];thresholds=[]
    for sign in [-1,1]:
        name='grazing-'+('minus' if sign<0 else 'plus');point={l:sign*p,q:sp.S.Zero}
        values_at=block.subs(point);denominators=[sp.fraction(sp.cancel(v))[1] for v in block]
        evaluated=[sp.cancel(v.subs(point)) for v in denominators];det=sp.cancel(values_at.det())
        record={'point':list(point.items()),'originalBlock':block,'originalRationalDenominators':denominators,'evaluatedDenominators':evaluated,'matrix':values_at,'determinant':det,'finite':det.is_finite,'nonzero':det.is_zero is False,'scope':'local continuous outgoing threshold only; no all-real pole exclusion or field'}
        J.emit(name+'-determinant-and-domains',record);thresholds.append(record)
        require(all(v.is_finite is True and v.is_zero is False for v in evaluated),'threshold original domain '+name)
        require(det.is_finite is True and det.is_zero is False,'SINGULAR_OR_UNKNOWN_RECEIVING_THRESHOLD '+name)
    # Matched-end sign and polarization derivative; no transverse radiation integral.
    E=chart[:,:2];D=dual[:2,:];C=clean(D.subs(l,p)*U);require(C.det()!=0,'incident coordinate chart')
    B=clean(E*C);DB=clean(C.inv()*D);zeros('incident-chart-match',B.subs(l,p),U);zeros('doublet-dual',DB*B,sp.eye(2))
    P1=clean(rawends['RIGHT'].diff(lam).subs(lam,0));rightForcing=sp.zeros(5,2)
    forces=[]
    for column in range(2):
        force=get('source-result/full-local-force-column-'+str(column)+'.json');forces.append(force)
        rightForcing[:,column]=force['localRight']
        require(force['pressureAddition'].startswith('negative C00 R00 S01'),'actual force sign convention')
    zeros('end-RHS-minus-operator-sign',-P1.subs(l,p)*U,rightForcing)
    zeros('both-column-end-step',rightForcing,9*U)
    endPerturbation=clean(DB.subs(l,p)*P1.subs(l,p)*U);zeros('end-scalar-perturbation',endPerturbation,-9*sp.eye(2))
    delta=sp.cancel(-endPerturbation[0,0]/sp.diff(DT,l).subs(l,p));correction=clean(delta*B.diff(l).subs(l,p))
    zeros('delta-p-from-end-sign',delta,3/p)
    fullEndFirst=clean(P1.subs(l,p)*U+delta*(Pweak.diff(l).subs(l,p)*U+Pweak.subs(l,p)*B.diff(l).subs(l,p)))
    zeros('full-five-row-matched-end',fullEndFirst,sp.zeros(5,2))
    J.emit('end-matching-prerequisite',{'RIGHTOriginal':rawends['RIGHT'],'P1':P1,'savedForces':forces,'rightForcing':rightForcing,'C':C,'B':B,'DB':DB,'BprimeAtIncident':B.diff(l).subs(l,p),'DBprimeAtIncident':DB.diff(l).subs(l,p),'deltaP':delta,'fieldPolarizationCorrection':correction,'fullEndFirstOrder':fullEndFirst,'qRemainsIndependent':True,'T1Computed':False,'K1Computed':False})
    # Actual native face records supply homogeneous maps; direct drive is affine.
    fullnative=get('receiving/LEFT-native-1-original.json');require(fullnative['passed'] is True,'saved native join return')
    acoustic=fullnative['actual']['acoustic'];J.emit('native-face-original-operands',acoustic)
    uniformFields=S.T*psi
    def native_map(expr,label):
        mp={}
        for a in expr.free_symbols:
            if a.name in values:mp[a]=values[a.name]
            elif re.fullmatch('s11cdCurrentPlusAmplitude[0-4]',a.name):mp[a]=uniformFields[int(a.name[-1])]
            elif a.name in ['s11cdSpectralNormalMomentum','s11cdCurrentRightMomentum']:mp[a]=l
            elif a.name=='s11cdAcousticRightNormalMomentum':mp[a]=q
            else:raise ValueError('native face map leaf '+str(a))
        return clean(substitute(label,expr,mp,{l,q,*psi}))
    mu0=native_map(acoustic['CHEMICAL_AFFINITY_DRIVER'],'native-chemical-driver')
    mu_saved=wave(get('source-input/consumer/chemical-amplitude-grade-split.json')['zeroGrade'],'saved-chemical00-receiving')
    zeros('native-vs-saved-chemical00',mu0,mu_saved)
    kernels={n:native_map(v,'native-memory-'+n) for n,v in acoustic['MEMORY_KERNELS'].items()}
    zeros('native-V-kernel-zero',kernels['V'],0);zeros('native-X-kernel-zero',kernels['X'],0)
    zeros('native-R-whole-response',R,rho*omega/(q+omega*kernels['A']/rho))
    face_outputs=[];controls=[]
    for index,face in enumerate(acoustic['FACE_RECORDS']):
        sign=face['ORIENTATION'];require(sign in [-1,1],'native face orientation');name='plus' if sign==1 else 'minus'
        V0=native_map(face['OUTWARD_VELOCITY'],name+'-native-outward-velocity');lift=native_map(face['VIRTUAL_LIFT'],name+'-native-virtual-lift')
        P0=native_map(face['PRESSURE'],name+'-native-pressure');murow=row(mu0,name+'-mu-row')
        zeros(name+'-native-affine-closure-homogeneous',P0,R*(V0+kernels['A']*mu0/rho))
        zeros(name+'-native-mass-flux',native_map(face['RELATIVE_MASS_FLUX'],name+'-native-j'),q*P0/omega-rho*V0)
        zeros(name+'-native-affinity',native_map(face['AFFINITY'],name+'-native-affinity'),mu0-P0/rho)
        zeros(name+'-native-total-load',native_map(face['MECHANICAL_LOAD'],name+'-native-load'),lift*(P0+kernels['X']*(mu0-P0/rho)))
        zeros(name+'-face-mu-transverse',murow*B,sp.zeros(1,2))
        zeros(name+'-face-V-transverse',row(V0,name+'-V-row')*B,sp.zeros(1,2))
        zeros(name+'-face-P-transverse',row(P0,name+'-P-row')*B,sp.zeros(1,2))
        normal=get('source-result/'+name+'-normal-flat-arguments.json')['factor'].xreplace({oldq:q})
        zeros(name+'-original-lab-normal-join',normal,sign*sp.I*q*R)
        velocity_context=get('source-result/'+name+'-velocity-input.json')
        velocity_operands=get('source-result/'+name+'-held-profile-velocity-input.json')
        velocity_return=get('source-result/'+name+'-held-profile-velocity-return.json')
        require(velocity_context['nativeSource']==get('source-input/consumer/'+name+'-source-input.json'),'actual saved velocity source')
        require(velocity_context['liftEW']==list(U[4,:]),'actual velocity incident arguments')
        require(velocity_return['cancelled']==0 and velocity_operands['left']==velocity_operands['right']==0,'published zero direct velocity')
        V1=velocity_operands['left']
        J.emit(name+'-restored-direct-velocity',{'originalContext':velocity_context,'operands':velocity_operands,'return':velocity_return,'value':V1,'oldFunctionCalled':False})
        for column in range(2):
            chemical=get('source-result/chemical-'+str(column)+'-return.json');source=get('source-result/'+name+'-source-01-column-'+str(column)+'-return.json')
            inherited=get('source-result/'+name+'-chemical-to-normalized-source-'+str(column)+'-return.json')
            inherited_operands=get('source-result/'+name+'-chemical-to-normalized-source-'+str(column)+'-input.json')
            require(inherited['cancelled']==0 and inherited_operands['left']==source['value'],'inherited normalization actual source')
            # Join the old RIGHT operand to the newly bound native closure map;
            # do not recompute the completed old LEFT-minus-RIGHT identity.
            zeros(name+'-native-normalization-argument-'+str(column),inherited_operands['right'],kernels['A']*chemical['independentSigmaCoefficient']/rho)
            envelope=chemical['independentSigmaCoefficient']/10
            tag=sp.Symbol('receiving_mu1_column_'+str(column)+'_hat_l_minus_p') if envelope!=0 else sp.S.Zero
            mu=mu0+tag;V=V0+V1;P=R*(V+kernels['A']*mu/rho);v=q*P/(rho*omega);j=rho*(v-V);affinity=mu-P/rho;load=lift*(P+kernels['X']*affinity)
            J.emit(name+'-affine-face-column-'+str(column),{'nativeFace':face,'mu00':mu0,'V00':V0,'directChemicalEnvelope':envelope,'directChemicalTag':tag,'transformArgument':l-p,'transformEvaluated':False,'inheritedSourceNormalization':{'operands':inherited_operands,'return':inherited},'totalMu':mu,'totalV':V,'pressure':P,'bulkOutwardVelocity':v,'massFlux':j,'affinity':affinity,'mechanicalResponse':kernels['X']*affinity,'totalLoad':load,'labPressureDerivative':sign*sp.I*q*P,'outwardPressureDerivative':sp.I*q*P,'directPressureCount':1,'sourcePressureNotAddedAgain':True})
            zeros(name+'-affine-mass-equation-'+str(column),j,kernels['A']*affinity+kernels['V']*V)
            face_outputs.append((name,column))
            if name=='plus' and column==0:
                T=sym(envelope,'full_weak_tanh_variable');point={q:sp.S.One,l:sp.S.Zero,tag:envelope.subs(T,0),**{a:sp.S.Zero for a in psi}}
                nonzero('omit-direct-chemical-control',P.subs(point),P.subs(tag,0).subs(point));controls.append('omit-direct-chemical')
                point={q:sp.S.One,l:sp.S.Zero,tag:sp.S.Zero,**{a:sp.S.Zero for a in psi}};point[psi[3]]=1
                nonzero('reversed-lab-normal-control',(sign*sp.I*q*P).subs(point),(-sign*sp.I*q*P).subs(point));controls.append('reversed-lab-normal')
                actual=(murow*B)[0,0];wrong=(murow*B.subs(l,p))[0,0]
                nonzero('incoming-for-receiving-projection-control',actual.subs(l,0),wrong.subs(l,0));controls.append('incoming-for-receiving-projection')
    require(len(face_outputs)==4 and len(controls)==3,'all native face/column routes and controls')
    return {'executionStatus':'COMPLETED_RECEIVING_BLOCK_AND_FACE_PREREQUISITE','sourceCorrespondence':'full five rows at independent receiving l,q','transverseCouplingIdenticallyZero':True,'thresholdsChecked':2,'endRHSAndPolarizationChecked':True,'affineFaceColumns':4,'newControls':controls,'transverseSolvabilityOrMomentsComputed':False,'T1Computed':False,'R1Computed':False,'G1Computed':False,'K1Computed':False,'matchedAsymptoticProjectionComplete':False,'physicalReflectionWeightComputed':False,'allRealRegularity':False,'receivingInverse':False,'powerComputed':False,'leakageFactor':None,'oldFunctionsCalled':False,'scienceAcceptance':False}
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
        result=J.stage('receiving-block-and-face-prerequisite',{'manifestSha256':gate['manifestSha256'],'buildReviewSha256':gate['buildReviewRecordSha256']},lambda:scientific_work(manifest,J,ns['decode']));code=0
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
