#!/usr/bin/env python3
"""Finite source projection, cross current and real-axis C3 algebra. No field or power integral."""
import argparse, ast, hashlib, json, os, re, resource, shutil, sys, time, traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','containment','Journal','decode')
ROWS=('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')
FIELDS=('u_1','u_2','u_3','theta','e_W')
GRADES=('00','10','01')
VERDICT='CLEAR FOR THIS REAL-AXIS RECEIVING AND CROSS-FLUX BUILD'

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
    require(g['status']=='READY_FOR_ONE_FIRST_ORDER_RECEIVING_REGULAR_INSTRUMENT','one finite receiving prerequisite gate')
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
    require(method['methodAssessed'] is True and method['literalVerdict']=='CLEAR FOR THIS REAL-AXIS RECEIVING AND CROSS-FLUX METHOD','actual method')
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

def variations(signs):
    require(all(type(s) is int and s in (-1,0,1) for s in signs),'exact sign sequence')
    live=[s for s in signs if s]
    return sum(a!=b for a,b in zip(live,live[1:]))


def exact_sign(value):
    if value.is_zero is True:return 0
    if value.is_positive is True:return 1
    if value.is_negative is True:return -1
    raise ValueError('UNRESOLVED_EXACT_SIGN '+str(value))


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
    def selected(alias,keys):
        # Only a literal JSON selection before the same established decoder.
        raw=json.loads((directory/alias).read_text())
        for key in keys:raw=raw[key]
        J.emit('selected-'+alias.replace('/','-').removesuffix('.json')+'-'+'-'.join(map(str,keys)),{'file':copied[alias],'selector':keys,'selectedRawJSON':raw})
        return decode(raw)
    def clean(v):return v.applyfunc(sp.cancel) if isinstance(v,sp.MatrixBase) else sp.cancel(v)
    def zero(name,left,right):
        if isinstance(left,sp.MatrixBase):
            require(isinstance(right,sp.MatrixBase) and left.shape==right.shape,'shape '+name)
            J.emit(name+'-matrix-operands',{'left':left,'right':right})
            for i in range(left.rows):
                for j in range(left.cols):J.zero(name+'-'+str(i)+'-'+str(j),left[i,j],right[i,j])
        else:J.zero(name,left,right)
    def symbol(expr,name):
        found=[s for s in expr.free_symbols if s.name==name];require(len(found)==1,'original symbol '+name);return found[0]
    inherited=[]
    def old_scalar(folder,name,left=None,right=None):
        op=get(folder+'/'+name+'-input.json');ret=get(folder+'/'+name+'-return.json')
        J.emit('inherited-'+folder+'-'+name,{'operands':op,'return':ret,'functionsCalled':False})
        require(ret['cancelled']==0,'inherited exact zero '+name)
        if left is not None:require(op['left']==left,'inherited actual LEFT '+name)
        if right is not None:require(op['right']==right,'inherited actual RIGHT '+name)
        inherited.append(folder+'/'+name);return op
    def old_matrix(folder,name,left=None,right=None):
        op=get(folder+'/'+name+'-matrix-operands.json')
        J.emit('inherited-'+folder+'-'+name+'-full',op)
        if left is not None:require(op['left']==left,'original full LEFT '+name)
        if right is not None:require(op['right']==right,'original full RIGHT '+name)
        require(op['left'].shape==op['right'].shape,'old matrix shape')
        for i in range(op['left'].rows):
            for j in range(op['left'].cols):old_scalar(folder,name+'-'+str(i)+'-'+str(j),op['left'][i,j],op['right'][i,j])
        return op
    context=get('original/source/extended-binding-context.json');values=context['numeric'];eps=context['epsilon']
    physical=get('original/source/physical-input.json')
    require(context['physical']==physical and context['frequencyOverride']=={'old':'1','actual':3},'same physical input and held frequency')
    require(context==get('source-input/local/extended-binding-context.json'),'original local binding context')
    for name,value in values.items():require(value==sp.Rational(physical['parameters'][name]) if name!='omega' else value==3,'actual physical binding '+name)
    omega,L=values['omega'],values['L_W'];h1,h2=(values[n] for n in ['s11cdTangentialMomentum1','s11cdTangentialMomentum2'])
    incident=get('source/incident-columns.json');U=incident['columns'];S=incident['chart']['matrixWeakFromUniform']
    require(incident['fieldOrder']==list(FIELDS),'original field order')
    p=sp.sympify(incident['selectedPoint']['mapping']['uniformNormal'])
    cs=sp.sympify(incident['selectedPoint']['mapping']['uniformSoundSpeed'])
    require(cs==sp.sympify(manifest['scope']['effectiveCs']) and sp.sympify(incident['selectedPoint']['mapping']['uniformFrequency'])==omega,'actual selected speed and held frequency')
    J.emit('selected-physical-speed-context',{'selectedMapping':incident['selectedPoint']['mapping'],'effectiveCs':cs,'originalInputCs0':physical['parameters']['c_s0'],'effectiveScope':manifest['scope']['effectiveCs'],'physicalInput':physical,'notOriginalInputSpeedFallback':True})
    full=get('receiving/complete-transformed-operator.json');E5,D5,transformed=full['chart'],full['dual'],full['matrix']
    l=symbol(E5,'receiving_block_l');q=symbol(transformed,'receiving_block_q');C=transformed[2:5,2:5];R3=D5[2:5,:];E3=E5[:,2:5]
    H2=h1*h1+h2*h2;K2=l*l+H2;end=get('receiving/end-matching-prerequisite.json');B=end['B'];Bp=end['BprimeAtIncident'];delta=end['deltaP']
    matched=get('matching/matched-transverse-amplitudes.json')
    require(matched['B']==B and matched['Bprime']==Bp and matched['deltaP']==delta,'actual saved moving basis')
    J.emit('moving-basis-symbol-identities',{'B':B,'Bprime':Bp,'T1':matched['T1'],'U':U,'receivingMomentum':l,'BFreeSymbols':sorted(B.free_symbols,key=str)})
    require(B.free_symbols<={l} and not Bp.free_symbols and not U.free_symbols and not matched['T1'].free_symbols,'foreign moving-basis or incident symbol')
    dual_proof=old_matrix('receiving','geometric-dual',right=sp.eye(5))
    J.emit('live-basis-product-arguments',{'savedDualLeft':dual_proof['left'],'D5':D5,'E5':E5,'physicalMatrix':full['physicalMatrix'],'savedTransformed':transformed,'newArgumentJoinsNotOldResidualReplay':True})
    zero('new-dual-left-argument-join',dual_proof['left'],D5*E5)
    zero('new-complete-transformation-join',transformed,D5*full['physicalMatrix']*E5)
    old_matrix('receiving','transverse-to-scalars-coupling',transformed[2:5,:2],sp.zeros(3,2))
    old_matrix('receiving','scalars-to-transverse-coupling',transformed[:2,2:5],sp.zeros(2,3))
    old_matrix('receiving','transverse-block',transformed[:2,:2],sp.eye(2)*full['DTCandidate'])
    old_matrix('receiving','incident-chart-match',B.subs(l,p),U)
    correspondence=old_matrix('receiving','FULL-offwave-source-correspondence',left=full['physicalMatrix'])
    raw_proof=old_matrix('receiving','LEFT-raw-vs-receiving',right=full['physicalMatrix'])
    # Original source/row families, their actual return operands and native provenance.
    local=get('receiving/local-receiving-assembly.json');pressure=get('receiving/pressure-receiving-assembly.json')
    cells=get('source-input/local/all-local-cells.json');require(len(cells)==400 and len(local['entries'])==100,'full original local00 row census')
    J.emit('receiving-source-basis-ancestry',{'fullNativeCorrespondence':get('receiving/FULL-offwave-source-correspondence-matrix-operands.json'),'local':local,'pressure':pressure,'fullTransformed':full,'nativeSourceArguments':get('receiving/LEFT-raw-source-arguments.json'),'oldSourceFunctionsCalled':False})
    # New operand connections use saved terms/returns; no cell, source or substitution producer.
    baseline=[cell for cell in cells if cell['grade']==[0,0]]
    require([entry['cell'] for entry in local['entries']]==baseline,'complete ordered local00 ancestry')
    local_terms=sp.zeros(5)
    for index,entry in enumerate(local['entries']):
        cell=entry['cell'];require(cell['field']==FIELDS[cell['fieldColumn']],'local physical column')
        require(entry['newReceivingMultiplier']==(sp.I*l)**cell['xOrder'],'original receiving jet argument')
        # This connects the saved contribution to its declared coefficient, not its original child sum.
        zero('new-local-term-argument-'+str(index),entry['term'],cell['coefficient']*entry['newReceivingMultiplier'])
        local_terms[ROWS.index(cell['row']),cell['fieldColumn']]+=entry['term']
    zero('new-local-saved-summand-join',local['matrix'],local_terms)
    psi=sp.Matrix(sp.symbols('receiving_field_0:5'));source_rows={}
    for face in ['plus','minus']:
        inherited_source=get('receiving/'+face+'-inherited-source00.json')
        source=get('source-input/sources/'+face+'-source-jets-00.json')
        original_op=get('source-input/sources/'+face+'-source-jet-reconstruction-00-input.json')
        original_ret=get('source-input/sources/'+face+'-source-jet-reconstruction-00-return.json')
        sub_input=get('receiving/'+face+'-receiving-source00-substitution-input.json')
        sub_return=get('receiving/'+face+'-receiving-source00-substitution-return.json')
        J.emit(face+'-restored-source00-frame',{'inherited':inherited_source,'source':source,'operands':original_op,'return':original_ret,'substitutionInput':sub_input,'substitutionReturn':sub_return,'fieldOrder':list(FIELDS),'fields':psi,'momentum':l,'frequency':omega,'tangents':[h1,h2],'acceptedSubstitutionExecutionDependency':True})
        require(inherited_source=={'record':source,'operands':original_op,'return':original_ret},'original source00 receipt '+face)
        require(original_ret['cancelled']==0 and original_op['left']==source['source'] and original_op['right']==source['reconstruction'],'original source00 proof '+face)
        require(sub_input['expression']==source['source'] and set(sub_input['allowedSymbols'])=={l,*psi},'original source00 substitution input '+face)
        jet_map=dict(sub_input['mapping'])
        require(len(jet_map)==len(sub_input['mapping']) and set(jet_map)==source['source'].free_symbols,'complete source00 substitution keys '+face)
        for atom,value in jet_map.items():
            jets=[jet for jet in source['jets'] if jet['atom']==atom];require(len(jets)==1,'unique original source jet '+face)
            spec=jets[0]['spec'];powers=spec['spatialOrders']
            require(spec['name']==atom.name and len(powers)==3,'original native jet specification '+face)
            expected=psi[FIELDS.index(spec['channel'])]*(-sp.I*omega)**spec['timeOrder']
            for momentum,order in zip([l,h1,h2],powers):expected*= (sp.I*momentum)**order
            zero('new-'+face+'-jet-frame-'+atom.name,value,expected)
        require(set(sub_return['remaining'])==sub_return['value'].free_symbols and sub_return['value'].free_symbols<={l,*psi},'original source00 return symbols '+face)
        linear=old_scalar('receiving',face+'-source00-linear',left=sub_return['value'])
        rows=[entry['receivingSourceRow'] for entry in pressure['entries'] if entry['originalPiece']['face']==face]
        require(len(rows)==10 and all(row==rows[0] for row in rows),'one complete source row per face '+face)
        source_rows[face]=rows[0]
        zero('new-'+face+'-source-row-proof-argument',linear['right'],(rows[0]*psi)[0])
    require(pressure['sourceMomentum']==l and pressure['normalOutputDepth']==q and pressure['wholeFlatResponseOnce'] is True and pressure['flatSupport']=='input=output=l','actual flat receiving arguments')
    assembly=get('source/first-order-pressure-assembly.json');pressure_terms=sp.zeros(5);seen=set()
    for index,entry in enumerate(pressure['entries']):
        piece=entry['originalPiece'];face,slot=piece['face'],piece['slot'];key=(entry['row'],face,slot)
        require(key not in seen,'unique receiving pressure slot');seen.add(key)
        originals=[record for record in assembly if record['row']==entry['row'] and record['column']==0]
        require(len(originals)==1 and piece in originals[0]['pieces'],'actual full source/consumer piece')
        factor=get('source/'+face+'-'+slot+'-flat-arguments.json')
        require(entry['factorArguments']==factor and factor['factor']==piece['flatNormalFactor'],'same saved flat factor')
        olddepth=get('source/receiving-sheet.json')['depth']
        require(entry['newFactor']==factor['factor'].xreplace({olddepth:q}),'actual flat output depth mapping')
        zero('new-pressure-term-argument-'+str(index),entry['contribution'],piece['consumers']['00']*entry['newFactor']*source_rows[face])
        pressure_terms[ROWS.index(entry['row']),:]+=entry['contribution']
    require(seen=={(row,face,slot) for row in ROWS for face in ['plus','minus'] for slot in ['pressure','normal']},'complete receiving pressure census')
    zero('new-pressure-saved-summand-join',pressure['matrix'],pressure_terms)
    zero('new-physical-local-pressure-join',full['physicalMatrix'],local['matrix']+pressure['matrix'])
    # The two native LEFT operands live in the physical weak basis, before D5/E5.
    uniform=get('source-input/incident/LEFT-invariant-P.json')
    native_input=get('receiving/uniform-native-receiving-substitution-input.json');native_return=get('receiving/uniform-native-receiving-substitution-return.json')
    raw_args=get('receiving/LEFT-raw-source-arguments.json');raw_saved=get('source-input/incident/LEFT-raw-source-binding.json')
    raw_input=get('receiving/LEFT-raw-end-substitution-input.json');raw_return=get('receiving/LEFT-raw-end-substitution-return.json')
    J.emit('native-physical-basis-argument-frames',{'S':S,'uniform':uniform,'nativeInput':native_input,'nativeReturn':native_return,'rawArguments':raw_args,'rawSaved':raw_saved,'rawInput':raw_input,'rawReturn':raw_return,'fullCorrespondence':correspondence,'rawProof':raw_proof,'oldSubstitutionFunctionsCalled':False,'acceptedSubstitutionExecutionDependency':True})
    require(uniform['actual']==uniform['expected']==native_input['expression'],'original uniform source input')
    require(dict(native_input['mapping'])=={symbol(uniform['actual'],'uniformNormal'):l,symbol(uniform['actual'],'uniformFrequency'):omega,symbol(uniform['actual'],'uniformPhysicalDepth'):q},'actual native receiving substitution frame')
    require(all(raw_args['saved'][key]==raw_saved[key] for key in ['nativeSource','raw','map']) and raw_input['expression']==raw_saved['raw'],'actual raw LEFT source input')
    require(raw_args['originalFiniteOriginNotUsed'] is True,'unused finite-origin context stays separate')
    require(dict(raw_args['originalMapping'])=={a:b for a,b in raw_saved['map']['mappingPairs'] if a.name not in ['eta_bg','sigma_W']},'original raw LEFT physical bindings')
    require(dict(raw_input['mapping'])==dict(raw_args['newArgumentMapping'])=={symbol(raw_input['expression'],'weak_end_p'):l,symbol(raw_input['expression'],'weak_end_q'):q},'actual raw LEFT receiving frame')
    for label,ret in [('uniform',native_return),('raw',raw_return)]:
        require(set(ret['remaining'])==ret['value'].free_symbols and ret['value'].free_symbols<={l,q},'native return symbols '+label)
    zero('new-native-uniform-proof-right-argument',correspondence['right'],S*native_return['value']*S.T)
    zero('new-native-raw-proof-left-argument',raw_proof['left'],S*raw_return['value']*S.T)
    zero('new-chart-selection-R3',R3,sp.Matrix([[l/K2,h1/K2,h2/K2,0,0],[0,0,0,1,0],[0,0,0,0,1]]))
    zero('new-chart-selection-E3',E3,sp.Matrix([[l,0,0],[h1,0,0],[h2,0,0],[0,1,0],[0,0,1]]))
    sheet=get('source/receiving-sheet.json');oldl,oldq=sheet['momentum'],sheet['depth'];beta=sheet['beta'];rho=values['rho_m']
    def transport(name,expr,allowed):
        expr=sp.sympify(expr);mapping={oldl:l,oldq:q}
        J.emit(name+'-argument-input',{'original':expr,'mapping':list(mapping.items()),'allowed':sorted(allowed,key=str),'oldMomentumEqualsNew':oldl==l,'oldDepthEqualsNew':oldq==q})
        result=expr.xreplace(mapping)
        J.emit(name+'-argument-return',{'result':result,'freeSymbols':sorted(result.free_symbols,key=str)})
        require(result.free_symbols<=allowed,'FOREIGN_RECEIVING_ARGUMENT '+name)
        return result
    zero('new-normal-depth-law',sheet['depthSquared'].xreplace({oldl:l}),omega**2/(cs*cs)-H2-l*l)
    zero('new-normal-vs-total-momentum',p*p,omega**2/(cs*cs)-H2)
    zero('new-native-closure-beta',beta,omega*values['Lambda_A_0']/(rho*(1-sp.I*omega*values['tau_A'])))
    domain=get('receiving/outgoing-domain.json')
    require(domain['qSquared']==p*p-l*l and domain['beta']==beta,'actual outgoing domain operands')
    require(sheet['outgoingConvention']=='q>=0 real inside; q=+i sqrt(-q^2) outside','saved positive outgoing sheet')
    J.emit('full-row-selection-and-outgoing-sheet',{'incident':incident,'chart':E5,'dual':D5,'selectedRows':[2,3,4],'selectedColumns':[2,3,4],'C3':C,'R3':R3,'E3':E3,'sourceSheet':sheet,'receivingDomain':domain,'H2':H2,'K2':K2,'p':p,'L':L,'noPhysicalLossProjection':True})
    require(H2.is_positive is True and p.is_positive is True,'positive chart domain and threshold')
    # Project the COMPLETE five physical force rows before selecting C3 rows.
    forces=[get('source/full-local-force-column-'+str(col)+'.json') for col in range(2)]
    T=symbol(forces[0]['etaForce'],'full_weak_tanh_variable')
    Feta=sp.Matrix.hstack(*(f['etaForce'] for f in forces));Fsigma=sp.Matrix.hstack(*(f['sigmaForce'] for f in forces))
    A=9*(1+T)/2-2*(1-T*T)
    sourceproofs=[]
    for col in range(2):
        for row,name in enumerate(ROWS):
            for grade,target in [('10',Feta),('01',Fsigma)]:
                rec=get('source/'+name+'-'+grade+'-column-'+str(col)+'-return.json');op=get('source/'+name+'-'+grade+'-column-'+str(col)+'-input.json')
                require(rec['sourceForce']==target[row,col],'actual saved force row/grade '+name)
                for entry in op['cells']:require(entry['cell'] in cells,'full source/receiving common row cell')
                sourceproofs.append({'row':name,'column':col,'grade':grade,'operands':op,'return':rec})
        require(forces[col]['incidentMomentum']==p and forces[col]['receivingMomentum']==oldl,'original independent input/output arguments')
        require(forces[col]['pressureAddition'].startswith('negative C00 R00 S01'),'original forcing sign')
    J.emit('complete-force-row-ancestry',sourceproofs)
    projection=D5*Feta
    table=sp.Matrix([[-2*sp.I*p*A,sp.I*A],[sp.I*(sp.Rational(1,5)+4*l*p)*A/K2,0],[sp.I*(l-p)*A/(5*K2),0],[0,0],[0,0]])
    zero('new-full-five-row-eta-table',projection,table)
    right=sp.Matrix.hstack(*(f['localRight'] for f in forces))
    zero('new-right-step-projection',R3.subs(l,p)*right,sp.zeros(3,2))
    old_matrix('receiving','both-column-end-step',right,9*U)
    J.emit('new-full-projection-before-selection',{'Feta':Feta,'Fsigma':Fsigma,'D5':D5,'allEtaRows':projection,'allSigmaRows':D5*Fsigma/10,'R3':R3,'rightStep':right,'nonlocalizedTransverseRows':projection[:2,:],'physicalGradePath':'eta=lambda,sigma=lambda/10','noTransverseStepDiscarded':True})
    # NEW combined smooth derivative, not a replay of old source/profile or moment routines.
    Aprime=sp.diff(A,T)*(1-T*T)/L
    zero('new-combined-step-derivative',Aprime,(sp.Rational(9,2)+4*T)*(1-T*T)/L)
    zero('new-step-endpoint-jump',A.subs(T,1)-A.subs(T,-1),sp.Integer(9))
    profiles=get('source/restored-profile-rules.json');require(len(profiles['certificates'])==40,'original profile certificate count')
    for rec in profiles['certificates']:
        require(rec['L']==L,'same physical profile scale')
        require(all(v['cancelled']==0 for v in rec['identities']),'inherited profile derivative returns')
    J.emit('inherited-profile-arguments',profiles)
    Ah=sp.Symbol('receiving_Aprime_hat_l_minus_p');Q=l-p
    eta_projected=sp.Matrix([[Ah/(5*K2),0],[0,0],[0,0]])
    zero('new-Q-before-cusp-projection',R3*U,sp.Matrix([[sp.I*Q/(5*K2),0],[0,0],[0,0]]))
    J.emit('whole-step-distribution-transport',{'A':A,'Aprime':Aprime,'inputArgument':p,'outputArgument':l,'Q':Q,'profileRule':profiles['rule'],'L':L,'forwardConvention':'integral exp(-i Q x); inverse dl/(2pi)','stepTransform':'9*pi*delta(Q)+9*PV(1/(iQ))+hat(A-9H)','combinedIdentity':'Q*Ahat=-i*AprimeHat before applying any receiving cusp','projectedEta':eta_projected,'zeroTransferAprimeHat':A.subs(T,1)-A.subs(T,-1),'fullStepNotLocalized':True,'transformEvaluated':False})
    # Full pressure source/consumer maps, with actual both-face chemical normalization.
    assembly=get('source/first-order-pressure-assembly.json');require(len(assembly)==10 and len(pressure['entries'])==20,'full original pressure row/slot census')
    native=selected('original/native/LEFT-original.json',['actual','acoustic'])
    require(native==get('receiving/native-face-original-operands.json'),'same original native face bodies')
    pressure_force=sp.zeros(5,2);face_maps={};tag_envelopes={};face_records=[];forcing_factors=[]
    field_symbols=set(sp.symbols('receiving_field_0:5'))
    memory_input=get('receiving/native-memory-A-substitution-input.json');memory_return=get('receiving/native-memory-A-substitution-return.json')
    J.emit('inherited-native-memory-A-arguments',{'input':memory_input,'return':memory_return,'originalKernel':native['MEMORY_KERNELS']['A']})
    require(memory_input['expression']==native['MEMORY_KERNELS']['A'],'actual unbound memory A kernel')
    A_memory=memory_return['value'];require(not A_memory.free_symbols,'fully bound native memory kernel')
    zero('new-native-memory-A-physical-binding',A_memory,values['Lambda_A_0']/(1-sp.I*omega*values['tau_A']))
    coefficient=clean(A_memory/rho)
    J.emit('native-normalization-multiplier',{'savedMemory':A_memory,'rho':rho,'rawRatio':A_memory/rho,'canonicalRatio':coefficient})
    for face,sign in [('plus',1),('minus',-1)]:
        nativeface=[f for f in native['FACE_RECORDS'] if f['ORIENTATION']==sign];require(len(nativeface)==1,'original face selection')
        closure=old_scalar('receiving',face+'-native-affine-closure-homogeneous')
        normal=transport('normal-'+face,get('source/'+face+'-normal-flat-arguments.json')['factor'],{l,q})
        old_scalar('receiving',face+'-original-lab-normal-join',normal)
        for col in range(2):
            aff=get('receiving/'+face+'-affine-face-column-'+str(col)+'.json');face_maps[(face,col)]=aff
            require(aff['nativeFace']==nativeface[0],'native face identity')
            for key in ['mu00','V00','totalMu','totalV','pressure','labPressureDerivative']:
                allowed={l,q,*field_symbols}|({aff['directChemicalTag']} if aff['directChemicalTag']!=0 else set())
                aff[key]=transport(face+'-'+str(col)+'-'+key,aff[key],allowed)
            require(transport(face+'-'+str(col)+'-transfer',aff['transformArgument'],{l})==Q,'actual Fourier transfer')
            norm=aff['inheritedSourceNormalization'];source=get('source/'+face+'-source-01-column-'+str(col)+'-return.json')
            require(norm['operands']['left']==source['value'] and norm['return']['cancelled']==0,'actual source normalization inherited return')
            proof=old_scalar('receiving',face+'-native-normalization-argument-'+str(col),norm['operands']['right'])
            chemical=get('source/chemical-'+str(col)+'-return.json')
            require(aff['directChemicalEnvelope']==chemical['independentSigmaCoefficient']/10,'physical sigma to lambda chemical map')
            require(aff['totalMu']==aff['mu00']+aff['directChemicalTag'],'actual affine chemical addition')
            require(aff['totalV']==aff['V00'],'actual held-profile zero direct velocity')
            velocity=get('receiving/'+face+'-restored-direct-velocity.json');require(velocity['value']==0 and velocity['return']['cancelled']==0,'inherited held velocity')
            zero('new-'+face+'-memory-normalization-argument-'+str(col),proof['right'],coefficient*chemical['independentSigmaCoefficient'])
            # New transport from complete physical affine pressure to the source-force tag.
            require(aff['mu00'].free_symbols<={l,*field_symbols},'only explicit five field leaves and receiving momentum')
            zero_field={v:sp.S.Zero for v in field_symbols}
            direct=aff['pressure'].subs(zero_field);tag=aff['directChemicalTag'];response=transport('response-'+face+'-'+str(col),sheet['flatResponse'],{l,q})
            zero('new-'+face+'-direct-pressure-count-'+str(col),direct,response*coefficient*tag)
            face_records.append({'face':face,'column':col,'affine':aff,'source':source,'chemical':chemical,'normalizationProof':proof,'nativeClosureProof':closure,'directPressure':direct,'sameNativeUnits':get('flux/LEFT-field-unit-join.json')})
    for col in range(2):
        for ri,row in enumerate(ROWS):
            hits=[v for v in assembly if v['row']==row and v['column']==col];require(len(hits)==1,'complete saved pressure row')
            rec=hits[0];require(len(rec['pieces'])==4 and {(v['face'],v['slot']) for v in rec['pieces']}=={(f,s) for f in ['plus','minus'] for s in ['pressure','normal']},'both faces and slots')
            total=sp.S.Zero;details=[]
            for piece in rec['pieces']:
                face,slot=piece['face'],piece['slot'];aff=face_maps[(face,col)]
                baseline=[v for v in pressure['entries'] if v['row']==row and v['originalPiece']['face']==face and v['originalPiece']['slot']==slot];require(len(baseline)==1,'same baseline consumer ancestry')
                base=baseline[0]
                require(base['originalPiece']==next(v for r in assembly if r['column']==0 and r['row']==row for v in r['pieces'] if v['face']==face and v['slot']==slot),'complete original receiving piece')
                require(piece['consumers']==base['originalPiece']['consumers'] and piece['split']==base['originalPiece']['split'] and piece['epsilon']==base['originalPiece']['epsilon'],'actual shared consumer grade/epsilon basis')
                args=get('source/'+face+'-'+slot+'-flat-arguments.json');require(piece['flatNormalFactor']==args['factor']==base['factorArguments']['factor'],'whole response and output normal argument')
                source=get('source/'+face+'-source-01-column-'+str(col)+'-return.json')
                require(piece['source00']==piece['source10']==0 and piece['source01Envelope']==source['value'],'actual first-order source split')
                require(piece['epsilon']==eps,'native epsilon once')
                original_tag=piece['transformTag'];tag=original_tag if piece['source01Envelope']!=0 else sp.S.Zero
                tag_envelopes[(face,col)]=(tag,piece['source01Envelope'])
                factor=transport('factor-'+str(col)+'-'+row+'-'+face+'-'+slot,piece['flatNormalFactor'],{l,q});term=piece['consumers']['00']*factor*tag
                require(not piece['consumers']['00'].free_symbols and factor.free_symbols<={q},'constant consumer and q-only native flat factor')
                forcing_factors.append({'row':row,'rowIndex':ri,'column':col,'face':face,'slot':slot,'consumer':piece['consumers']['00'],'factor':factor,'coefficient':piece['consumers']['00']*factor/10})
                # This is a new receiving transport; old pressure source functions remain uncalled.
                require(piece['transformConvention']=='integral exp(-i(l-p)x) envelope(x) dx; inverse dl/(2pi); not evaluated','exact Fourier convention')
                zero('new-pressure-piece-transport-'+str(col)+'-'+row+'-'+face+'-'+slot,transport('piece-'+str(col)+'-'+row+'-'+face+'-'+slot,piece['pressure01'],{l,q,original_tag}),term)
                affine_direct=aff['pressure'].subs({v:sp.S.Zero for v in field_symbols})
                normal_sign=1 if face=='plus' else -1
                from_affine=piece['consumers']['00']*affine_direct*(normal_sign*sp.I*q if slot=='normal' else 1)
                normalized_tag_map={original_tag:10*coefficient*aff['directChemicalTag']}
                zero('new-affine-force-transport-'+str(col)+'-'+row+'-'+face+'-'+slot,term.xreplace(normalized_tag_map)/10,from_affine)
                total+=term;details.append({'originalPiece':piece,'receivingBaselineEntry':base,'flatArguments':args,'receivingTerm':term,'affineChemicalMap':aff['inheritedSourceNormalization']})
            zero('new-pressure-row-transport-'+str(col)+'-'+row,transport('row-'+str(col)+'-'+row,rec['pressure01'],{l,q,*[v['transformTag'] for v in rec['pieces']]}),total)
            pressure_force[ri,col]=-total/10
            J.emit('new-complete-pressure-row-'+str(col)+'-'+row,{'saved':rec,'pieces':details,'force':pressure_force[ri,col],'sign':'minus native operator action','sourceDerivativeAt':p,'normalDepthAt':q})
    J.emit('native-affine-source-transport',face_records)
    # Each sigma component has its own unevaluated Fourier tag; no hidden projection.
    sigma_tags=sp.zeros(5,2);profile_records=[]
    def localized(name,expr):
        J.emit(name+'-profile-input',{'expression':expr,'T':T,'L':L,'substitution':'T=tanh(x/L)','arguments':'Q=l-p'})
        poly=sp.Poly(expr,T,extension=[sp.I,sp.sqrt(595)])
        endpoints=[poly.eval(v) for v in [-1,1]];quotient,remainder=sp.div(poly,sp.Poly(1-T*T,T,extension=[sp.I,sp.sqrt(595)]))
        J.emit(name+'-profile-return',{'polynomial':poly.as_expr(),'coefficients':poly.all_coeffs(),'endpoints':endpoints,'factorQuotient':quotient.as_expr(),'remainder':remainder.as_expr(),'analyticInduction':'d/dx=(1-T^2)/L*d/dT preserves an endpoint-zero factor; every derivative decays exponentially','machineMeasureProof':False})
        require(all(v==0 for v in endpoints) and remainder.is_zero,'NONLOCALIZED_OR_UNKNOWN_THREE_FIELD_PROFILE '+name)
        profile_records.append(name)
    localized('combined-Aprime',Aprime)
    for col in range(2):
        for ri in range(5):
            value=Fsigma[ri,col]/10;localized('sigma-'+str(col)+'-'+str(ri),value)
            sigma_tags[ri,col]=sp.Symbol('sigma_'+str(col)+'_'+str(ri)+'_hat_l_minus_p') if value!=0 else sp.S.Zero
    for (face,col),(tag,expr) in tag_envelopes.items():localized('source-envelope-'+face+'-'+str(col),expr)
    # The exact U_B receiving cancellation acts on a common transform of the
    # saved polynomial relation; do not assign unrelated tags to proportional profiles.
    projected_sigma=R3*Fsigma/10
    zero('new-complete-UB-sigma-projection',projected_sigma[:,1],sp.zeros(3,1))
    zero('new-complete-UB-pressure-projection',(R3*pressure_force)[:,1],sp.zeros(3,1))
    zero('new-complete-UB-eta-projection',(R3*U)[:,1],sp.zeros(3,1))
    J.emit('full-three-row-source-class',{'R3':R3,'eta':eta_projected,'sigmaProfiles':Fsigma/10,'sigmaProjectionBeforeTransform':projected_sigma,'sigmaIndependentTags':sigma_tags,'pressureForce':pressure_force,'pressureProjection':R3*pressure_force,'UBThreeRowsZero':True,'tagRelations':'Apply the Fourier transform linearly to the saved profile relations; the UB independent display tags are not independent physical sources.','profileCertificates':profile_records,'transformsEvaluated':False})
    # Moving output projector on a delta-prime: f delta-prime=f(p)delta-prime-fprime(p)delta.
    Rprime=R3.diff(l).subs(l,p);moving=delta*(R3.subs(l,p)*Bp+Rprime*U)
    J.emit('moving-projector-distribution-operands',{'R3':R3,'B':B,'Bprime':Bp,'U':U,'deltaP':delta,'R3primeAtP':Rprime,'matchedT1':matched['T1'],'deltaCoefficient':moving,'deltaPrimeCoefficient':-delta*R3.subs(l,p)*U,'fourierJet':'2pi*[U delta+lambda*((deltaP Bprime+U T1)delta-deltaP U delta-prime)]','homogeneousRightEndOnly':True})
    zero('new-moving-projector-delta-prime-cancellation',moving,sp.zeros(3,2))
    zero('new-both-column-delta-prime-coefficient',-delta*R3.subs(l,p)*U,sp.zeros(3,2))
    zero('new-both-column-T1-delta-coefficient',R3.subs(l,p)*U*matched['T1'],sp.zeros(3,2))
    J.emit('declared-derivative-control-entry',{'row':0,'column':0,'projector':R3.subs(l,p),'basisDerivative':Bp,'rawEntry':(R3.subs(l,p)*Bp)[0,0],'baseline':moving[0,0],'mutant':(delta*R3.subs(l,p)*Bp)[0,0],'deltaP':delta})
    J.nonzero('drop-moving-projector-derivative-control',moving[0,0],(delta*R3.subs(l,p)*Bp)[0,0])
    # Actual opposite-leg native current; no current constructor, G0, Gref or K1 replay.
    flux=get('flux/physical-linear-survival-operands.json');JL=flux['leftUnboundCurrent']
    currentproof=old_matrix('flux','native-baseline-current-join',right=JL)
    for name,key in [('incident','G0'),('outward-reflected','Gref')]:
        weight=get('flux/'+name+'-positive-weight.json');J.emit('inherited-'+name+'-weight',weight)
        require(weight['matrix']==flux[key] and all(v.is_positive is True for v in weight['leadingPrincipalMinors']),'actual original signed current weight')
    binding=get('flux/LEFT-slab-binding.json');original_current=selected('flux/LEFT-native-current-source.json',['slabCurrent'])
    require(binding['original']==original_current==selected('original/native/LEFT-original.json',['actual','slab','SLAB_CURRENT_MATRIX']),'actual original current source body')
    epsilonproof=old_matrix('flux','LEFT-slab-epsilon-once',left=binding['result'])
    # New explicit transport join from the saved epsilon-stripped uniform-chart operand.
    zero('new-native-current-chart-join',epsilonproof['right'],eps**2*S.T*JL*S)
    units=get('flux/LEFT-field-unit-join.json');require(tuple(units['currentUnit'])==(0,-3,1),'actual inherited physical current units')
    km=symbol(JL,'flux_left_momentum');kp=symbol(JL,'flux_right_momentum')
    require(JL.free_symbols<={km,kp},'foreign native current momentum')
    Bm=B.subs(l,-p);Bp_inc=B.subs(l,p)
    require(get('receiving/incident-chart-match-matrix-operands.json')['left']==Bp_inc,'actual incoming leg restoration')
    original_flux_source=Path(manifest['nativeFluxSource']).read_text()
    original_function=next(n for n in ast.parse(original_flux_source).body if isinstance(n,ast.FunctionDef) and n.name=='scientific_work')
    requested_assignments={
        'JL':"JL=bound['LEFT']['path']['slab'].subs(lam,0)",
        'Bminus':'Bminus=B.subs(l,-p)',
        'J0':'J0=JL.subs({km:p,kp:p})',
        'Jminus':'Jminus=JL.subs({km:-p,kp:-p})',
        'G0':'G0=clean(U.H*J0*U)',
        'Gref':'Gref=clean(-Bminus.H*Jminus*Bminus)'}
    assignment_receipts=[]
    for name,expected in requested_assignments.items():
        actual=[n for n in original_function.body if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id==name]
        require(len(actual)==1 and ast.dump(actual[0],include_attributes=False)==ast.dump(ast.parse(expected).body[0],include_attributes=False),'actual completed diagonal-current source assignment '+name)
        assignment_receipts.append({'name':name,'originalSource':ast.get_source_segment(original_flux_source,actual[0]),'line':actual[0].lineno})
    # These are explicit old argument frames and saved returns, NOT another evaluation
    # of U.H*J0*U or -Bminus.H*Jminus*Bminus. Original execution provenance is inherited.
    J.emit('inherited-diagonal-current-argument-frames',{'sourcePath':manifest['nativeFluxSource'],'sourceSha256':sha(manifest['nativeFluxSource']),'assignments':assignment_receipts,'sameOriginalIncident':incident,'sameOriginalMovingBasis':end,'sameCurrentInput':currentproof,'frames':[{'name':'G0','bra':Bp_inc,'ket':Bp_inc,'current':JL,'momentumLegs':[[km,p],[kp,p]],'outwardSign':1,'savedReturn':flux['G0']},{'name':'Gref','bra':Bm,'ket':Bm,'current':JL,'momentumLegs':[[km,-p],[kp,-p]],'outwardSign':-1,'savedReturn':flux['Gref']}],'noDiagonalContractionReplayed':True,'provenance':'Exact completed source assignments plus live original operand identity and published weight returns; conditional on preserved original guarded execution, not a new numerical weight proof.'})
    minusplus=Bm.H*JL.subs({km:-p,kp:p})*Bp_inc;plusminus=Bp_inc.H*JL.subs({km:p,kp:-p})*Bm
    x=sp.Symbol('receiving_cross_flux_x',real=True)
    J.emit('opposite-leg-native-current-operands',{'JL':JL,'braMinus':Bm,'ketPlus':Bp_inc,'km':km,'kp':kp,'minusPlus':minusplus,'plusMinus':plusminus,'minusPlusPhase':sp.exp(2*sp.I*p*x),'plusMinusPhase':sp.exp(-2*sp.I*p*x),'inheritedG0':flux['G0'],'inheritedGref':flux['Gref'],'inheritedK1':flux['K1'],'harmonicAndEpsilonConvention':flux['epsilonConvention'],'units':units,'physicalPowerBalance':False})
    zero('new-cross-current-conjugate-relation',minusplus,plusminus.H)
    zero('new-cross-current-minus-plus',minusplus,sp.zeros(2))
    zero('new-cross-current-plus-minus',plusminus,sp.zeros(2))
    mutant=Bm.H*JL.subs({km:p,kp:p})*Bp_inc
    J.emit('wrong-bra-leg-control-operands',{'baseline':minusplus,'mutant':mutant,'bra':Bm,'ket':Bp_inc,'baselineLegs':[[-p,p]],'mutantLegs':[[p,p]]})
    J.nonzero('wrong-bra-momentum-sign-control',minusplus[1,1],mutant[1,1])
    # NEW global three-field determinant on the actual two outgoing rays.
    zero('new-C3-evenness',C,C.subs(l,-l))
    s=sp.Symbol('receiving_squared_l',real=True);C_s=C.xreplace({l*l:s})
    J.emit('even-block-transport',{'originalC3':C,'s':s,'Cs':C_s,'sheetSubstitution':p*p-q*q})
    require(not C_s.has(l),'UNRESOLVED_EVEN_BLOCK_REDUCTION')
    zero('new-even-block-reconstruction',C_s.subs(s,l*l),C)
    Cq=C_s.subs(s,p*p-q*q)
    def exclude_ray(name,polynomial,substitution,upper):
        r=sp.Symbol('receiving_ray_'+name,real=True);expr=sp.expand(polynomial.subs(q,substitution(r)))
        realpart,imagpart=sp.expand_complex(expr).as_real_imag();a=sp.Poly(realpart,r,domain=sp.QQ);b=sp.Poly(imagpart,r,domain=sp.QQ)
        J.emit(name+'-real-imag-input',{'originalPolynomial':polynomial,'substitution':substitution(r),'expression':expr,'real':a.as_expr(),'imaginary':b.as_expr(),'realCoefficients':a.all_coeffs(),'imaginaryCoefficients':b.all_coeffs(),'interval':[0,upper]})
        require(not(a.is_zero and b.is_zero),'IDENTICALLY_ZERO_OUTGOING_DOMAIN_POLYNOMIAL')
        # Save the exact Euclidean operands and every remainder, including zero components.
        aa,bb=a,b;euclid=[]
        while not bb.is_zero:
            quotient,remainder=aa.div(bb);euclid.append({'dividend':aa.as_expr(),'divisor':bb.as_expr(),'quotient':quotient.as_expr(),'remainder':remainder.as_expr()});aa,bb=bb,remainder
        u,v,g=sp.gcdex(a.as_expr(),b.as_expr(),r);G=sp.Poly(g,r,domain=sp.QQ)
        J.emit(name+'-gcd-return',{'euclideanSteps':euclid,'u':u,'v':v,'gcd':g,'degree':G.degree(),'nonzero':not G.is_zero})
        zero(name+'-exact-Bezout',u*a.as_expr()+v*b.as_expr(),g)
        require(not G.is_zero,'ZERO_COMMON_POLYNOMIAL')
        if G.degree()==0:
            J.emit(name+'-root-exclusion',{'gcd':g,'noCommonRealRoot':True,'reason':'nonzero constant exact Bezout identity','interval':[0,upper],'numericalRootFinder':False})
            return {'ray':name,'gcd':g,'count':0,'certificate':'constant Bezout'}
        sf=G.sqf_part().monic();sequence=[sf,sf.diff()];divisions=[]
        while not sequence[-1].is_zero:
            quotient,remainder=sequence[-2].div(sequence[-1]);divisions.append({'left':sequence[-2].as_expr(),'right':sequence[-1].as_expr(),'quotient':quotient.as_expr(),'negativeRemainder':-remainder.as_expr()})
            if remainder.is_zero:break
            sequence.append(-remainder)
        at_zero=[v.eval(0) for v in sequence]
        at_upper=[v.LC() if upper==sp.oo else v.eval(upper) for v in sequence]
        J.emit(name+'-Sturm-operands',{'gcd':g,'squareFree':sf.as_expr(),'factorization':[(f.as_expr(),m) for f,m in G.sqf_list()[1]],'sequence':[v.as_expr() for v in sequence],'divisions':divisions,'atZero':at_zero,'atUpper':at_upper,'upper':upper,'infinityRule':'positive infinity uses exact leading coefficient signs'})
        lo=[exact_sign(v) for v in at_zero];hi=[exact_sign(v) for v in at_upper]
        endpoints=[G.eval(0)]+([] if upper==sp.oo else [G.eval(upper)])
        count=variations(lo)-variations(hi)
        J.emit(name+'-root-exclusion',{'lowerSigns':lo,'upperSigns':hi,'lowerVariations':variations(lo),'upperVariations':variations(hi),'openIntervalRoots':count,'endpointGcdValues':endpoints,'upper':upper})
        require(all(v.is_zero is False for v in endpoints) and count==0,'REAL_AXIS_RECEIVING_POLE_OR_ENDPOINT_ROOT '+name)
        return {'ray':name,'gcd':g,'count':count,'certificate':'exact square-free Sturm and endpoints'}
    denominator_certificates=[]
    def domain_certificate(name,expr):
        J.emit(name+'-domain-input',{'expression':expr,'q':q,'allowedRays':['q=t,0<=t<=p','q=i*r,r>=0']})
        require(expr.free_symbols<={q},'FOREIGN_OR_UNREDUCED_DOMAIN_ARGUMENT '+name)
        joined=sp.together(expr);n,d=sp.fraction(joined)
        # No cancel: both the numerator and denominator of this original base
        # must be defined/nonzero before any later rational cancellation.
        J.emit(name+'-domain-rational-pair',{'original':expr,'joined':joined,'numerator':n,'denominator':d})
        group=[]
        for tag,poly in [('numerator',n),('denominator',d)]:
            require(sp.Poly(poly,q,extension=sp.I).as_expr()==sp.expand(poly),'exact Gaussian-rational denominator polynomial')
            for ray,sub,upper in [('propagating',lambda r:r,p),('evanescent',lambda r:sp.I*r,sp.oo)]:
                group.append(exclude_ray(name+'-'+tag+'-'+ray,poly,sub,upper))
        denominator_certificates.append({'name':name,'original':expr,'rationalPair':[n,d],'certificates':group})
    # Retain every original negative-power base even when together/cancel would
    # erase it. Domains are transported onto q with the same physical sheet.
    original_denominators=[]
    for i,entry in enumerate(C):
        bases=[]
        for position,node in enumerate(sp.preorder_traversal(entry)):
            if node.is_Pow and node.exp.is_negative is True:
                J.emit('entry-'+str(i)+'-base-'+str(len(bases))+'-input',{'entry':entry,'treePosition':position,'negativePower':node,'base':node.base,'exponent':node.exp})
                require(node.exp.is_Integer is True,'NONRATIONAL_ORIGINAL_DENOMINATOR_POWER')
                transported=node.base.xreplace({l*l:s}).subs(s,p*p-q*q)
                require(not transported.has(l),'original denominator even-sheet map')
                bases.append({'position':position,'base':node.base,'exponent':node.exp,'transported':transported})
        joined=sp.together(Cq[i]);n,d=sp.fraction(joined)
        J.emit('entry-'+str(i)+'-full-domain-operands',{'original':entry,'sheetEntry':Cq[i],'originalNegativePowers':bases,'joinedEntry':joined,'joinedNumerator':n,'joinedDenominator':d})
        for j,base in enumerate(bases):domain_certificate('entry-'+str(i)+'-base-'+str(j),base['transported'])
        domain_certificate('entry-'+str(i)+'-joined-denominator',d)
        original_denominators.append({'original':entry,'negativePowers':bases,'joinedDenominator':d})
    # The receiving chart has its own domain even if simplified C3 has lost K2.
    domain_certificate('chart-K2',p*p+H2-q*q)
    J.emit('new-determinant-input',{'C3':C,'Cq':Cq,'originalEntryDenominators':original_denominators,'K2':K2,'physicalDepthLaw':domain['qSquared'],'fullFiveFieldHasTransversePoles':True})
    rawdet=Cq.det(method='berkowitz');joined=sp.together(rawdet);rawnum,rawden=sp.fraction(joined)
    J.emit('new-raw-determinant-before-cancel',{'rawDeterminant':rawdet,'joined':joined,'rawNumerator':rawnum,'rawDenominator':rawden})
    domain_certificate('raw-determinant-denominator',rawden)
    det=sp.cancel(joined);N,D=sp.fraction(det)
    J.emit('new-determinant-return',{'rawDeterminant':rawdet,'joined':joined,'rawNumerator':rawnum,'rawDenominator':rawden,'cancelled':det,'numerator':N,'denominator':D,'numeratorCoefficients':sp.Poly(N,q,extension=sp.I).all_coeffs(),'denominatorCoefficients':sp.Poly(D,q,extension=sp.I).all_coeffs()})
    zero('new-determinant-rational-reconstruction',rawnum*D,N*rawden)
    domain_certificate('cancelled-determinant-denominator',D)
    for side,sign in [('minus',-1),('plus',1)]:
        threshold=get('receiving/grazing-'+side+'-determinant-and-domains.json')
        require(threshold['originalBlock']==C and dict(threshold['point'])=={l:sign*p,q:sp.S.Zero},'actual original threshold arguments')
        J.emit('inherited-'+side+'-threshold',threshold)
        zero('new-'+side+'-threshold-join',det.subs(q,0),threshold['determinant'])
    certificates=[exclude_ray('determinant-propagating',N,lambda r:r,p),exclude_ray('determinant-evanescent',N,lambda r:sp.I*r,sp.oo)]
    # The off-diagonal control is fixed before any outcome, never searched/fitted.
    mutant=Cq.copy();mutant[0,1]=sp.S.Zero;mutant[1,0]=sp.S.Zero
    J.emit('delete-coupling-control-input',{'baseline':Cq,'mutant':mutant,'deletedPositions':[[0,1],[1,0]],'q':0})
    mutantdet=mutant.subs(q,0).det(method='berkowitz')
    J.nonzero('delete-C3-coupling-control',det.subs(q,0),mutantdet)
    # Only after ray/domain clearance. No forcing solve, inverse field or transform evaluation.
    adj=Cq.adjugate(method='berkowitz');J.emit('adjugate-input',{'Cq':Cq,'determinant':det,'adjugate':adj,'rootCertificates':certificates})
    zero('new-right-adjugate-identity',Cq*adj,det*sp.eye(3))
    zero('new-left-adjugate-identity',adj*Cq,det*sp.eye(3))
    multiplier=clean(adj/det);r=sp.Symbol('receiving_growth_r',positive=True);growth=[]
    for i,value in enumerate(multiplier):
        ev=sp.cancel(value.subs(q,sp.I*r));n,d=sp.fraction(ev);pn,pd=sp.Poly(n,r,extension=sp.I),sp.Poly(d,r,extension=sp.I)
        entry={'index':i,'value':value,'evanescent':ev,'numerator':n,'denominator':d,'numeratorDegree':pn.degree(),'denominatorDegree':pd.degree(),'leadingNumerator':pn.LC(),'leadingDenominator':pd.LC(),'zeroNumerator':pn.is_zero}
        J.emit('large-momentum-entry-'+str(i),entry)
        require(not pd.is_zero and pd.LC().is_zero is False,'actual nonzero leading denominator')
        entry['growthPower']=0 if pn.is_zero else max(0,int(pn.degree()-pd.degree()));growth.append(entry)
    def degree_record(name,expr,var):
        require(expr.free_symbols<={var},'foreign growth argument '+name)
        value=sp.cancel(expr);n,d=sp.fraction(value);pn=sp.Poly(n,var,extension=[sp.I,sp.sqrt(595)]);pd=sp.Poly(d,var,extension=[sp.I,sp.sqrt(595)])
        record={'name':name,'actualCoefficient':expr,'variable':var,'numerator':n,'denominator':d,'numeratorDegree':pn.degree(),'denominatorDegree':pd.degree(),'leadingNumerator':pn.LC(),'leadingDenominator':pd.LC(),'zeroNumerator':pn.is_zero,'growthPower':0 if pn.is_zero else max(0,int(pn.degree()-pd.degree()))}
        J.emit('forcing-growth-'+name,record)
        require(not pd.is_zero and pd.LC().is_zero is False,'actual forcing denominator leading coefficient')
        physical_denominator=d if var==q else d.xreplace({l*l:p*p-q*q})
        domain_certificate('forcing-'+name+'-denominator',physical_denominator)
        return record
    projection_growth=[degree_record('R3-'+str(i),v,l) for i,v in enumerate(R3)]
    eta_growth=degree_record('Aprime-multiplier',1/(5*K2),l)
    pressure_growth=[]
    for i,term in enumerate(forcing_factors):
        factor_growth=degree_record('pressure-'+str(i),term['coefficient'],q)
        # Product bounds add positive exponents. No cancellation buys decay.
        row_power=max(projection_growth[j*R3.cols+term['rowIndex']]['growthPower'] for j in range(R3.rows))
        pressure_growth.append({'operand':term,'factor':factor_growth,'receivingProjectionGrowth':row_power,'productGrowth':row_power+factor_growth['growthPower']})
    chart_growth=[degree_record('E3-'+str(i),v,l) for i,v in enumerate(E3)]
    forcing_power=max([eta_growth['growthPower']]+[v['growthPower'] for v in projection_growth]+[v['productGrowth'] for v in pressure_growth])
    chart_extra=max(v['growthPower'] for v in chart_growth)
    power=max(v['growthPower'] for v in growth)+chart_extra+forcing_power
    J.emit('complete-source-and-chart-growth-ledger',{'eta':eta_growth,'sigmaProjection':projection_growth,'pressure':pressure_growth,'chart':chart_growth,'sourceGrowthPower':forcing_power,'chartPower':chart_extra,'receivingGrowthPower':max(v['growthPower'] for v in growth),'totalPower':power,'sufficientSchwartzDecayExponentForL1':power+2,'analyticRelation':'On evanescent rays r=sqrt(l^2-p^2) is asymptotic to |l|; finite compact portions are bounded by the established domains. Arbitrarily stronger Schwartz weights give all required weighted L2 classes.','domainCertificates':denominator_certificates,'noUniformNormOrRate':True})
    J.emit('conditional-three-field-response-class',{'symbolicMultiplier':multiplier,'entryGrowth':growth,'physicalChart':E3,'chartAdditionalPower':chart_extra,'forcingAdditionalPower':forcing_power,'sufficientCombinedPolynomialWeight':power,'rootCertificates':certificates,'profileCertificates':profile_records,'compactArgument':'Actual denominator and determinant nonvanishing give continuity on each compact outgoing ray including grazing; rational large-r degrees give polynomial growth. No numerical norm bound.','FourierArgument':'Aprime, every sigma and direct source envelope have endpoint-zero polynomial factors and hence Schwartz transforms. R3 and all flat normal factors have at most polynomial growth. Choose weights larger than the displayed receiving/chart growth for ordinary L1 and weighted L2.','conclusion':'Conditional outgoing weighted-L2 three-field uniqueness and Riemann-Lebesgue decay of E3 times the three-field contribution only. Analytic induction/measure arguments assessed, not machine theorem proofs.','notClaimed':['total matched field decay','decay rate','far-bulk angular flux','threshold-supported distribution uniqueness','nonuniform work balance','quadratic leakage coefficient'],'powerIdentityRequiredNext':True})
    return {'executionStatus':'COMPLETED_FINITE_SOURCE_PROJECTION_CROSS_FLUX_AND_REAL_AXIS_RECEIVING','inheritedZeros':len(inherited),'projectedUBThreeRowsZero':True,'crossCurrentZero':True,'realAxisC3PoleExclusion':certificates,'conditionalGrowthPowerIncludingForcingAndChart':power,'newControls':3,'fullFiveFieldInverse':False,'fieldComputed':False,'FourierIntegralEvaluated':False,'powerBalanceDerived':False,'leakageFactor':None,'oldFunctionsCalled':False,'scientificAcceptance':False}

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
        result=J.stage('source-projection-cross-flux-real-axis-receiving',{'manifestSha256':gate['manifestSha256'],'buildReviewSha256':gate['buildReviewRecordSha256']},lambda:scientific_work(manifest,J,ns['decode']));code=0
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
