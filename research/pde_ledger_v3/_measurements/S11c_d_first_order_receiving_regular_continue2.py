#!/usr/bin/env python3
"""Saved receiving/domain prefix; exact polynomial representation repair."""
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
    require(g['status']=='READY_FOR_ONE_FIRST_ORDER_RECEIVING_REGULAR_CONTINUATION2','one finite receiving prerequisite gate')
    for key,p in [('worker',__file__),('manifest',manifest_path),('launcher',manifest['launcher'])]:
        require(g[key+'Sha256']==sha(p),'actual '+key)
    require(g['sourcePins']==manifest['sourcePins'],'source census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'pin '+p)
    require(g['librarySourcePins']==manifest['librarySourcePins'],'exact installed polynomial library sources')
    for p,h in manifest['librarySourcePins'].items():require(sha(p)==h,'library source pin '+p)
    for key in ['buildReviewRecord','methodRecord','executionAuthority']:
        require(sha(g[key])==g[key+'Sha256'],'gate '+key)
    require(g['methodRecord']==manifest['methodRecord'],'actual method record')
    review=json.loads(Path(g['buildReviewRecord']).read_text())
    require(sha(review['report'])==review['reportSha256'],'literal build report bytes')
    require(review['reviewers']==['claude'] and review['literalVerdict']==VERDICT and review['buildAssessed'] is True,'actual scoped Claude build')
    for key in ['workerSha256','manifestSha256','launcherSha256']:
        require(review[key]==g['baseReviewPins'][key],'original review '+key)
    require(g['independentBuildClearance'] is False and g['localToolingAuthority'] is True,'honest local continuation authority')
    require(sha(manifest['baseWorker'])==g['baseReviewPins']['workerSha256'],'original reviewed source')
    require(sha(manifest['priorRecord'])==g['priorRecordSha256'],'inspected failed prefix')
    require(sha(g['toolingRepairRecord'])==g['toolingRepairRecordSha256'],'repair bytes')
    repair=json.loads(Path(g['toolingRepairRecord']).read_text())
    require(repair['workerSha256']==sha(__file__) and repair['testsPassed'] is True and repair['unchangedRemainingTailAST'] is True,'tested finite continuation')
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


def exact_gcd_triplet(a,b,variable,engine):
    """Same monic Bezout specification, with explicit one-zero inputs."""
    require(not (a.is_zero and b.is_zero),'IDENTICALLY_ZERO_OUTGOING_DOMAIN_POLYNOMIAL')
    if b.is_zero:return engine.S.One/a.LC(),engine.S.Zero,a.monic().as_expr()
    if a.is_zero:return engine.S.Zero,engine.S.One/b.LC(),b.monic().as_expr()
    return engine.gcdex(a.as_expr(),b.as_expr(),variable)

def exact_polynomial_reconstruction(name,polynomial,variable,journal,engine):
    journal.emit(name+'-polynomial-input',{'polynomial':polynomial,'variable':variable,'coefficientField':'QQ(i)'})
    native=engine.Poly(polynomial,variable,extension=engine.I)
    reconstructed=native.as_expr();expanded=engine.expand(polynomial)
    journal.emit(name+'-polynomial-reconstruction',{'polynomial':polynomial,'variable':variable,'nativePolynomial':native,'coefficients':native.all_coeffs(),'reconstructed':reconstructed,'expandedOriginal':expanded,'structurallyEqual':reconstructed==expanded})
    journal.zero(name+'-exact-polynomial-reconstruction',engine.expand(reconstructed),expanded)


def scientific_work(manifest,J,decode):
    # Entire immutable failed run, including every original input and complete return.
    priorRoot=Path(manifest['priorRoot']);directory=J.out/'prior';directory.mkdir();priorCopies={}
    require({str(p.relative_to(priorRoot)) for p in priorRoot.rglob('*') if p.is_file()}==set(manifest['priorFiles']),'complete failed tree census')
    for name,record in manifest['priorFiles'].items():
        src=priorRoot/name;require(sha(src)==record['sha256'] and src.stat().st_size==record['bytes'],'prior bytes '+name)
        dst=directory/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);require(sha(dst)==record['sha256'],'prior copy '+name)
        priorCopies[name]={'source':str(src),'path':str(dst.relative_to(J.out)),**record}
    save(J.out/'prior-copy-index.json',priorCopies)
    copied={}
    for alias,record in manifest['savedFiles'].items():
        path=directory/'complete/prior/complete/saved'/alias
        require(sha(path)==record['sha256']==sha(record['path']),'actual saved original '+alias)
        copied[alias]={'source':record['path'],'path':str(path.relative_to(J.out)),'sha256':record['sha256']}
    save(J.out/'saved-copy-index.json',copied)
    artifact=json.loads((directory/'complete/artifact-index.json').read_text());previous={};cache={}
    for name,record in artifact.items():
        pth=directory/'complete'/name;require(sha(pth)==record['sha256'] and pth.stat().st_size==record['bytes'],'complete previous artifact '+name)
        previous[name]=json.loads(pth.read_text())
    oldDirectory=directory/'complete/prior/complete';oldArtifact=json.loads((oldDirectory/'artifact-index.json').read_text())
    def prior(name):
        key='original/'+name
        if key not in cache:
            pth=oldDirectory/name;require(sha(pth)==oldArtifact[name]['sha256'],'original prefix artifact '+name)
            cache[key]=decode(json.loads(pth.read_text()))
        return cache[key]
    def completed(name):
        if name not in cache:cache[name]=decode(previous[name])
        return cache[name]
    def get(alias):
        key='saved/'+alias
        if key not in cache:cache[key]=decode(json.loads((oldDirectory/'saved'/alias).read_text()))
        return cache[key]
    def clean(v):return v.applyfunc(sp.cancel) if isinstance(v,sp.MatrixBase) else sp.cancel(v)
    def zero(name,left,right):
        if isinstance(left,sp.MatrixBase):
            require(isinstance(right,sp.MatrixBase) and left.shape==right.shape,'shape '+name)
            J.emit(name+'-matrix-operands',{'left':left,'right':right})
            for i in range(left.rows):
                for j in range(left.cols):J.zero(name+'-'+str(i)+'-'+str(j),left[i,j],right[i,j])
        else:J.zero(name,left,right)
    def symbol(expr,name):
        found=[s for s in expr.free_symbols if s.name==name];require(len(found)==1,'actual restored symbol '+name);return found[0]
    failure=json.loads((directory/'complete/failure.json').read_text())
    require(failure['traceback'].rstrip().endswith('ValueError: exact Gaussian-rational denominator polynomial'),'exact unfinished representation predicate')
    require(len(artifact)==99 and list(artifact)[-1]=='entry-3-base-0-domain-rational-pair.json','exact complete previous boundary')
    literal_zero={'text':'0','srepr':'Integer(0)'};restored_zeros=[]
    for name,value in previous.items():
        if name.endswith('-return.json') and isinstance(value,dict) and 'cancelled' in value:
            require(value['cancelled']==literal_zero,'completed polynomial proof '+name)
            stem=name.removesuffix('-return.json');require(stem+'-input.json' in artifact and stem+'-raw.json' in artifact,'full completed proof operands')
            restored_zeros.append({'input':artifact[stem+'-input.json'],'raw':artifact[stem+'-raw.json'],'return':artifact[name]})
    require(len(restored_zeros)==12,'twelve complete Bezout returns')
    inheritedPrefix=completed('restored-complete-prefix.json');context=completed('unpublished-tail-context-reconstitution.json')
    require(len(inheritedPrefix['zeroProofs'])==639 and len(inheritedPrefix['inheritedRecords'])==171 and len(inheritedPrefix['controls'])==2,'original prefix census')
    for proof in inheritedPrefix['zeroProofs']:
        for key in ['input','raw','return']:require(proof[key]==oldArtifact[proof[key]['path']],'original full proof receipt')
        require(prior(proof['return']['path'])['cancelled']==sp.S.Zero,'original literal zero return')
    for name in inheritedPrefix['inheritedRecords']:
        value=prior(name);require(value['return']['cancelled']==sp.S.Zero and value['functionsCalled'] is False,'original inherited return')
    for control in inheritedPrefix['controls']:
        for key in ['input','return']:require(control[key]==oldArtifact[control[key]['path']],'original full control receipt')
        value=prior(control['return']['path']);require(value['zero'] is False and value['finite'] is True,'original responsive control')
    state=prior('full-row-selection-and-outgoing-sheet.json');transport=prior('even-block-transport.json');classes=prior('full-three-row-source-class.json')
    C=transport['originalC3'];s=transport['s'];C_s=transport['Cs'];p=state['p'];R3=state['R3'];E3=state['E3'];H2=state['H2'];K2=state['K2'];domain=state['receivingDomain']
    l=symbol(state['chart'],'receiving_block_l');q=symbol(C,'receiving_block_q')
    require(C==state['C3'] and C==get('receiving/grazing-plus-determinant-and-domains.json')['originalBlock'],'actual original C3')
    require(inheritedPrefix['sourceClass']==classes and inheritedPrefix['crossCurrent']==prior('opposite-leg-native-current-operands.json'),'actual complete source/current prefix')
    require(context['originalCs']==C_s and context['sheetSubstitution']==transport['sheetSubstitution'],'actual published context arguments')
    Cq=context['reconstitutedCq'];forcing_factors=context['forcingFactors'];profile_records=context['profileRecordNames']
    require(len(forcing_factors)==40 and profile_records==classes['profileCertificates'],'complete published forcing context')
    # These expressions are now published returns from the prior continuation;
    # do not repeat their sheet substitution or forty coefficient products.
    for term in forcing_factors:
        record=prior('new-complete-pressure-row-'+str(term['column'])+'-'+term['row']+'.json')
        matches=[v for v in record['pieces'] if v['originalPiece']['face']==term['face'] and v['originalPiece']['slot']==term['slot']]
        require(len(matches)==1 and matches[0]['originalPiece']['consumers']['00']==term['consumer'] and matches[0]['receivingBaselineEntry']['newFactor']==term['factor'],'published forcing argument identity')
    J.emit('restored-complete-source-and-domain-prefix',{'priorFailure':failure,'originalCompletePrefix':inheritedPrefix,'publishedContext':context,'newBezoutProofsRestored':restored_zeros,'oldSourceCurrentFunctionsCalled':False,'completeContextRecomputed':False,'nonoperationAttestation':'Conditional on unchanged original control flow reaching the published failure; unpublished polynomial comparison is not accepted or restored.'})
    def exclude_ray(name,polynomial,substitution,upper):
        r=sp.Symbol('receiving_ray_'+name,real=True);expr=sp.expand(polynomial.subs(q,substitution(r)))
        realpart,imagpart=sp.expand_complex(expr).as_real_imag();a=sp.Poly(realpart,r,domain=sp.QQ);b=sp.Poly(imagpart,r,domain=sp.QQ)
        J.emit(name+'-real-imag-input',{'originalPolynomial':polynomial,'substitution':substitution(r),'expression':expr,'real':a.as_expr(),'imaginary':b.as_expr(),'realCoefficients':a.all_coeffs(),'imaginaryCoefficients':b.all_coeffs(),'interval':[0,upper]})
        require(not(a.is_zero and b.is_zero),'IDENTICALLY_ZERO_OUTGOING_DOMAIN_POLYNOMIAL')
        # Save the exact Euclidean operands and every remainder, including zero components.
        aa,bb=a,b;euclid=[]
        while not bb.is_zero:
            quotient,remainder=aa.div(bb);euclid.append({'dividend':aa.as_expr(),'divisor':bb.as_expr(),'quotient':quotient.as_expr(),'remainder':remainder.as_expr()});aa,bb=bb,remainder
        J.emit(name+'-gcd-call-arguments',{'real':a.as_expr(),'imaginary':b.as_expr(),'variable':r,'realZero':a.is_zero,'imaginaryZero':b.is_zero,'realLeadingCoefficient':a.LC(),'imaginaryLeadingCoefficient':b.LC(),'specification':'monic gcd and exact Bezout; one-zero adapter only'})
        u,v,g=exact_gcd_triplet(a,b,r,sp);G=sp.Poly(g,r,domain=sp.QQ)
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
        if name=='entry-3-base-0':
            saved=completed(name+'-domain-rational-pair.json');original=completed(name+'-domain-input.json')
            require(original['expression']==expr==saved['original'] and original['q']==q,'exact failed negative-power domain frame')
            joined,n,d=saved['joined'],saved['numerator'],saved['denominator']
            J.emit('restored-failed-domain-pair',{'input':original,'pair':saved,'receipts':[artifact[name+'-domain-input.json'],artifact[name+'-domain-rational-pair.json']],'newWorkBegins':'exact expanded polynomial reconstruction; no completed ray calculation replayed'})
        else:
            J.emit(name+'-domain-input',{'expression':expr,'q':q,'allowedRays':['q=t,0<=t<=p','q=i*r,r>=0']})
            require(expr.free_symbols<={q},'FOREIGN_OR_UNREDUCED_DOMAIN_ARGUMENT '+name)
            joined=sp.together(expr);n,d=sp.fraction(joined)
            # No cancel: both the numerator and denominator of this original base
            # must be defined/nonzero before any later rational cancellation.
            J.emit(name+'-domain-rational-pair',{'original':expr,'joined':joined,'numerator':n,'denominator':d})
        group=[]
        for tag,poly in [('numerator',n),('denominator',d)]:
            exact_polynomial_reconstruction(name+'-'+tag,poly,q,J,sp)
            for ray,sub,upper in [('propagating',lambda r:r,p),('evanescent',lambda r:sp.I*r,sp.oo)]:
                group.append(exclude_ray(name+'-'+tag+'-'+ray,poly,sub,upper))
        denominator_certificates.append({'name':name,'original':expr,'rationalPair':[n,d],'certificates':group})
    # Retain every original negative-power base even when together/cancel would
    # erase it. Domains are transported onto q with the same physical sheet.
    original_denominators=[]
    for i in range(3):
        frame=completed('restored-first-entry-domain-operands.json')['value'] if i==0 else completed('entry-'+str(i)+'-full-domain-operands.json')
        require(frame['original']==C[i] and frame['sheetEntry']==Cq[i] and frame['originalNegativePowers']==[],'actual completed entry frame')
        name='entry-'+str(i)+'-joined-denominator'
        pair=completed('restored-first-domain-pair.json')['pair'] if i==0 else completed(name+'-domain-rational-pair.json')
        require(pair['original']==frame['joinedDenominator'],'same completed denominator')
        group=[];full=[]
        for tag in ['numerator','denominator']:
            for ray in ['propagating','evanescent']:
                stem=name+'-'+tag+'-'+ray
                args=completed('restored-first-ray-polynomials.json')['saved'] if i==0 and tag=='numerator' and ray=='propagating' else completed(stem+'-real-imag-input.json')
                gcd=completed(stem+'-gcd-return.json');decision=completed(stem+'-root-exclusion.json');proof=completed(stem+'-exact-Bezout-return.json')
                require(args['originalPolynomial']==pair[tag] and args['interval']==[0,p if ray=='propagating' else sp.oo],'exact completed ray arguments')
                require(gcd['nonzero'] is True and gcd['degree']==0 and decision['gcd']==gcd['gcd'] and decision['noCommonRealRoot'] is True and proof['cancelled']==sp.S.Zero,'complete actual constant Bezout decision')
                group.append({'ray':stem,'gcd':gcd['gcd'],'count':0,'certificate':'constant Bezout'})
                full.append({'arguments':args,'gcd':gcd,'decision':decision,'proof':proof,'receipt':artifact[stem+'-root-exclusion.json']})
        denominator_certificates.append({'name':name,'original':pair['original'],'rationalPair':[pair['numerator'],pair['denominator']],'certificates':group})
        original_denominators.append({'original':frame['original'],'negativePowers':frame['originalNegativePowers'],'joinedDenominator':frame['joinedDenominator']})
        J.emit('restored-completed-entry-'+str(i),{'frame':frame,'pair':pair,'allFourRays':full,'functionsCalled':False,'summaryReconstructedFromCompletePublishedDecisions':True})
    for i in range(3,len(C)):
        entry=C[i]
        if i==3:
            frame=completed('entry-3-full-domain-operands.json');base_input=completed('entry-3-base-0-input.json')
            require(frame['original']==entry and frame['sheetEntry']==Cq[i] and len(frame['originalNegativePowers'])==1 and frame['originalNegativePowers'][0]['base']==base_input['base'],'actual partial entry3 frame')
            bases=frame['originalNegativePowers'];d=frame['joinedDenominator']
            J.emit('restored-partial-entry3',{'frame':frame,'originalBaseOperand':base_input,'nextDomain':'entry-3-base-0','completedEntryDomainsReplayed':False})
        else:
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
    return {'executionStatus':'COMPLETED_FINITE_SOURCE_PROJECTION_CROSS_FLUX_AND_REAL_AXIS_RECEIVING','inheritedZeros':171,'projectedUBThreeRowsZero':True,'crossCurrentZero':True,'realAxisC3PoleExclusion':certificates,'conditionalGrowthPowerIncludingForcingAndChart':power,'newControls':1,'restoredPrefixControls':2,'restoredNewZeros':639+len(restored_zeros),'fullFiveFieldInverse':False,'fieldComputed':False,'FourierIntegralEvaluated':False,'powerBalanceDerived':False,'leakageFactor':None,'oldFunctionsCalled':False,'scientificAcceptance':False}

def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest);verify_invocation(args,gate,sys.argv)
    pins={**manifest['sourcePins'],**manifest['librarySourcePins'],**{r['path']:r['sha256'] for r in manifest['savedFiles'].values()},str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate),**{str(Path(manifest['priorRoot'])/n):v['sha256'] for n,v in manifest['priorFiles'].items()}}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False);J=None;result={};code=1;start=time.monotonic()
    try:
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(manifest['helperSource']).read_text()),'unchanged-journal-and-containment','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        ns.update(sp=sp,Str=Str);J=ns['Journal'](args.out)
        result=J.stage('saved-prefix-real-axis-receiving-continuation',{'manifestSha256':gate['manifestSha256'],'buildReviewSha256':gate['buildReviewRecordSha256']},lambda:scientific_work(manifest,J,ns['decode']));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result);sys.stderr.write(result['traceback'])
    finally:
        if (args.out/'saved-copy-index.json').exists():
            for v in json.loads((args.out/'saved-copy-index.json').read_text()).values():pins[str(args.out/v['path'])]=v['sha256']
        if (args.out/'prior-copy-index.json').exists():
            for v in json.loads((args.out/'prior-copy-index.json').read_text()).values():pins[str(args.out/v['path'])]=v['sha256']
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
