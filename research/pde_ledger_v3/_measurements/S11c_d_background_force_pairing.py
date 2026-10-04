#!/usr/bin/env python3
"""Finite background force/area/rate comparisons; no physical work or leakage closure."""
import argparse, ast, hashlib, json, os, resource, shutil, sys, time, traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','containment','Journal','decode')
VERDICT='CLEAR FOR THIS BACKGROUND FORCE-PAIRING BUILD AND AMENDED METHOD'

def require(v,message):
    if v is not True: raise ValueError(message)

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
    require({n.name for n in nodes}==set(HELPERS),'unchanged helper census')
    return ast.Module(body=nodes,type_ignores=[])

def verify_gate(path,manifest_path,manifest):
    g=json.loads(Path(path).read_text())
    require(g['status']=='READY_FOR_ONE_BACKGROUND_FORCE_PAIRING_INSTRUMENT','one bounded gate')
    for key,p in [('worker',__file__),('manifest',manifest_path),('launcher',manifest['launcher']),('method',manifest['methodPath'])]:
        require(g[key+'Sha256']==sha(p),'actual '+key)
    require(g['sourcePins']==manifest['sourcePins'],'source census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'pin '+p)
    for key in ('buildReviewRecord','executionAuthority'):
        require(sha(g[key])==g[key+'Sha256'],'gate '+key)
    review=json.loads(Path(g['buildReviewRecord']).read_text())
    require(sha(review['report'])==review['reportSha256'],'literal report bytes')
    require(review['reviewers']==['claude'] and review['literalVerdict']==VERDICT and review['buildAssessed'] is True and review['methodAssessed'] is True,'actual combined scoped assessment')
    for key in ('workerSha256','manifestSha256','launcherSha256','methodSha256'):require(review[key]==g[key],'review '+key)
    auth=json.loads(Path(g['executionAuthority']).read_text())
    require(auth['scienceExecutionsAuthorized']==1 and auth['scope']==manifest['scope'] and auth['reviewers']==['claude'],'standing bounded authority')
    require(auth['AGENTSSha256']==sha(ROOT/'AGENTS.md'),'current authority')
    require(g['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py') and g['supervisor']==str(M/'S11c_d_end_normalization_run.py'),'actual containment paths')
    require(sha(g['sharedGuard'])==g['guardSha256'] and sha(g['supervisor'])==g['supervisorSha256'],'containment bytes')
    require(g['scope']==manifest['scope'] and g['scientificRunsAuthorized']==1 and g['pooledExecution'] is True,'exact scope')
    require(manifest['resources']=={'memoryBytes':4*1024**3,'poolGiB':16,'zeroSwap':True,'cpuCount':1,'threads':1,'tasksMax':32,'hostReserveGiB':4,'durationLimits':None},'fixed resources')
    return g

def named(value,name):
    found=[v for k,v in value if str(k)==name]
    require(len(found)==1,'unique native label '+name)
    return found[0]

def scientific_work(manifest,J,decode):
    copied={};directory=J.out/'saved';directory.mkdir()
    for alias,record in manifest['savedFiles'].items():
        src=Path(record['path']);require(sha(src)==record['sha256'],'saved '+alias)
        dst=directory/alias;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
        require(sha(dst)==record['sha256'],'copy '+alias)
        copied[alias]={'source':str(src),'path':str(dst.relative_to(J.out)),'sha256':record['sha256']}
        tmp=J.out/'saved-copy-index.next';tmp.write_text(json.dumps(copied,indent=2)+'\n');tmp.replace(J.out/'saved-copy-index.json')
    registries={};raw={}
    for alias,record in manifest['registryInputs'].items():
        view=json.loads((directory/record['savedAlias']).read_text());raw[alias]=view
        export=Path(record['originalSource']);require(sha(export)==view['sourceSha256'],'original registry '+alias)
        lines=export.read_text().splitlines()
        display=ast.literal_eval(lines[view['displaySourceLine']-1].strip().removeprefix("'display': ").removesuffix(','))
        constructor=ast.literal_eval(lines[view['constructorSourceLine']-1].strip().removeprefix("'value': _restore(").removesuffix('),'))
        require(display==view['completePublishedDisplay'] and constructor==view['completePublishedConstructorString'],'literal original export '+alias)
        J.emit('original-'+alias,{'source':str(export),'sourceSha256':view['sourceSha256'],'display':display,'constructor':constructor,'registryKey':view['registryKey'],'originalLines':[view['displaySourceLine'],view['constructorSourceLine']]})
        if alias=='legacy_residual':
            J.emit('legacy-residual-not-decoded',{'arithmeticCertificate':False,'faceLabelsNotUsed':True,'logicalConstructorsPreserved':True})
        else:
            registries[alias]=J.stage('restore-'+alias,{'constructor':constructor},lambda c=constructor,d=display:decode({'text':d,'srepr':c}))
    def case(alias,labels):
        entries=[v for k,v in registries[alias] if tuple(str(x) for x in k)==tuple(str(x) for x in labels)]
        require(len(entries)==1,'unique full case '+alias+str(labels));return entries[0]
    op=case('operator',('LAB_HELD','RHO4_CONSTANT'));support=case('support',('LAB_HELD','RHO4_CONSTANT'))
    opv=named(op,'VALUE');sv=named(support,'VALUE');body=named(opv,'BODY_FORCE');sbody=named(sv,'BODY_FORCE')
    faces=list(named(opv,'PER_FACE_TRACTION'));sfaces=list(named(sv,'PER_FACE_TRACTION'))
    J.emit('original-background-selection',{'operator':op,'support':support,'operatorFaces':faces,'supportFaces':sfaces,'legacyLabelsUsed':False})
    require([int(k) for k,v in faces]==[1,-1] and [int(k) for k,v in sfaces]==[1,-1],'original face signs and positions')
    require(len(faces)==len(sfaces)==2 and all(len(v)==4 for k,v in faces+sfaces),'complete four-vector faces')
    allsymbols=set().union(*(v.atoms(sp.Symbol) for v in registries.values()))
    def sym(name):
        selected=[v for v in allsymbols if v.name==name];require(len(selected)==1,'native symbol '+name);return selected[0]
    eps=sym('epsilon_shape');eta=sym('eta_bg');sigma=sym('sigma_W');W=sym('W_0');L=sym('L_W');kap=sym('kappa_W');w=sym('w1_profile')
    grad=sp.Matrix([sigma*sym('w1_profile_d'+str(i)) for i in range(1,4)])
    lap=sum(sym('w1_profile_d'+str(i)+'d'+str(i)) for i in range(1,4))
    ut=sp.Matrix([sym('u_'+str(i)+'_t') for i in range(1,4)]);et=sym('e_W_t');zt=sym('zeta_c_t')
    J.emit('native-symbols-domains',{'epsilon':eps,'eta':eta,'sigma':sigma,'W0':W,'L':L,'gradient':grad,'laplacianProfile':lap,'rates':[ut,et,zt],'W0Positive':W.is_positive,'LPositive':L.is_positive})
    require(W.is_positive is True and L.is_positive is True,'native positive W0/L domain')
    e=named(body,'E_W');J.zero('saved-E-explicit-source-grade',e,-W**3*kap*sigma*(1+2*eta*w)*lap/L)
    J.emit('background-force-grade-scope',{'savedE':e,'publishedHeader':named(op,'MULTIGRADE'),'literalBodyU':named(body,'U'),'unknownHigherE':'UNSUPPLIED; no physical sigma2 or higher completion','actualNewProductGradesNotHeader':True})
    for j,(sign,t) in enumerate(faces):J.zero('saved-face-E-origin-'+str(j),t[3],sign*e/W)
    # Zero dimensions do not carry evidence of homogeneity or equal physical units.
    dims=named(op,'DIMENSION_L_T_M');sdims=named(support,'DIMENSION_L_T_M')
    unit_records=[]
    def unit_slot(name,value,unit,svalue,sunit):
        rec={'slot':name,'operatorValue':value,'operatorUnit':unit,'supportValue':svalue,'supportUnit':sunit,'classification':'VACUOUS_LITERAL_ZERO_NO_UNIT_AGREEMENT' if value==0 else 'NONZERO_NATIVE_UNIT_COMPARISON'}
        J.emit('unit-'+name,rec)
        if value!=0:require(unit==sunit,'native nonzero slot unit '+name)
        unit_records.append(rec)
    for j in range(3):unit_slot('U'+str(j),named(body,'U')[j],dims[0][1][0][1][j],named(sbody,'U')[j],sdims[0][1][0][1][j])
    for j,label in [(1,'THETA'),(2,'E_W')]:unit_slot(label,named(body,label),dims[0][1][j][1],named(sbody,label),sdims[0][1][j][1])
    for j in range(2):
        for k in range(4):unit_slot('face'+str(j)+'-'+str(k),faces[j][1][k],dims[1][1][j][1][k],sfaces[j][1][k],sdims[1][1][j][1][k])
    # New comparisons on actual saved operands, never old geometry producers.
    areas=[];rate_records=[]
    def rect(name,expression):
        def work():
            coeffs=[{'grade':[i,j],'value':sp.cancel(sp.diff(expression,eta,i,sigma,j).subs({eta:0,sigma:0}))} for i in range(2) for j in range(2)]
            retained=sum(v['value']*eta**v['grade'][0]*sigma**v['grade'][1] for v in coeffs)
            return {'fullNewExpression':expression,'coefficients':coeffs,'retainedComparison':retained,'rawComplement':expression-retained,'physicalHigherOrders':'UNSUPPLIED'}
        return J.stage(name,{'expression':expression,'variables':[eta,sigma],'orders':[[0,0],[1,0],[0,1],[1,1]],'newWorkComparisonOnly':True},work)
    for j,sign in enumerate((1,-1)):
        for dof in ('DELTA_W','ZETA_C'):
            nrec=case('normal',('LAB_HELD',sign,dof));mrec=case('measure',('LAB_HELD',sign,dof));vrec=case('velocity',('LAB_HELD',sign,dof))
            nv=named(nrec,'VALUE');mv=named(mrec,'VALUE');vv=named(vrec,'VALUE')
            J.emit('geometry-'+str(j)+'-'+dof,{'normal':nrec,'measure':mrec,'velocity':vrec,'orders':{'normalEntry0':0,'normalEntry1':1,'measureEntry0':0,'measureEntry1':1,'velocity':1}})
            require(len(nv)==len(mv)==3 and nv[0].shape==(4,1),'complete saved geometry shape')
            expected=sp.Matrix([*list(-grad/2),sign])
            for k in range(4):J.zero('saved-normal-'+str(j)+'-'+dof+'-'+str(k),nv[0][k],expected[k])
            J.zero('saved-background-area-'+str(j)+'-'+dof,mv[0],1)
            J.zero('velocity-order-'+str(j)+'-'+dof,vv,eps*(W*et/2 if dof=='DELTA_W' else sign*zt))
            for k,v in enumerate(nv[1]):J.zero('normal-increment-order-'+str(j)+'-'+dof+'-'+str(k),v,eps*sp.diff(v,eps))
            J.zero('area-increment-order-'+str(j)+'-'+dof,mv[1],eps*sp.diff(mv[1],eps))
            if dof=='DELTA_W':n=nv[0]
        rawarea=sp.sqrt(sum(n[k]**2 for k in range(4)))
        ar=rect('new-c2-area-comparison-'+str(j),rawarea)
        J.zero('raw-area-square-'+str(j),rawarea**2,1+(grad.dot(grad))/4)
        J.zero('new-retained-area-'+str(j),ar['retainedComparison'],1)
        c2=J.stage('new-area-sigma2-coefficient-'+str(j),{'rawArea':rawarea,'sigma':sigma},lambda:sp.diff(rawarea,sigma,2).subs(sigma,0)/2)
        J.zero('area-debt-leading-'+str(j),c2,sum(sym('w1_profile_d'+str(k))**2 for k in range(1,4))/8)
        exactn=n/rawarea;vSaved=(named(case('velocity',('LAB_HELD',sign,'DELTA_W')),'VALUE')+named(case('velocity',('LAB_HELD',sign,'ZETA_C')),'VALUE'))/eps
        h=zt+sign*W*et/2+sign*grad.dot(ut)/2;r=sp.Matrix([*list(ut),h]);exactV=exactn.dot(r)
        hSaved=(vSaved-n[:3,0].dot(ut))/n[3];hMixed=(exactV-n[:3,0].dot(ut))/n[3]
        J.emit('rate-normalization-'+str(j),{'savedNormal':n,'sourceExactNormalCandidate':exactn,'rawNorm':rawarea,'savedV':vSaved,'exactVFromSourceRate':exactV,'rateVector':r,'savedReconstruction':hSaved,'mixedExactVelocityReconstruction':hMixed,'physicalChoiceMade':False})
        J.zero('source-unit-normal-'+str(j),exactn.dot(exactn),1)
        J.zero('saved-rate-reconstruction-'+str(j),hSaved,h)
        J.zero('exact-normal-velocity-'+str(j),exactV,(sign*zt+W*et/2)/rawarea)
        J.zero('normal-density-conversion-'+str(j),n.dot(r),rawarea*exactV)
        J.zero('mixed-rate-full-debt-'+str(j),hMixed-h,sign*(exactV-vSaved))
        rr=rect('new-mixed-rate-debt-'+str(j),hMixed-h);J.zero('mixed-rate-retained-zero-'+str(j),rr['retainedComparison'],0)
        rate_records.append({'sign':sign,'rate':r,'area':rawarea,'mixedRateDebt':hMixed-h,'savedRate':hSaved});areas.append(rawarea)
    J.zero('equal-face-raw-area',areas[0],areas[1]);A=areas[0]
    def pull(tractions,area):
        result=sp.zeros(5,1)
        for sign,t in tractions:
            for k in range(3):result[k]+=area*(t[k]+sign*t[3]*grad[k]/2)
            result[3]+=area*sign*W*t[3]/2;result[4]+=area*t[3]
        return result
    mapped=J.stage('new-formal-background-pullback',{'faces':faces,'areaCandidate':A,'gradient':grad,'W0':W,'physicalAreaConvention':'UNSUPPLIED'},lambda:pull(faces,A))
    expected=sp.Matrix([*list(A*e*grad/W),A*e,0])
    for k in range(5):J.zero('formal-pullback-identity-'+str(k),mapped[k],expected[k])
    differences=[]
    for k in range(3):
        full=mapped[k]-named(body,'U')[k];rec=rect('new-face-body-U-grade-'+str(k),full);J.zero('face-body-U-retained-'+str(k),rec['retainedComparison'],0);differences.append(rec)
    ecompare=rect('new-face-body-E-grade',mapped[3]-e);J.zero('face-body-E-retained',ecompare['retainedComparison'],0)
    hold=J.stage('new-formal-support-pullback',{'supportFaces':sfaces,'areaCandidate':A,'physicalAreaConvention':'UNSUPPLIED'},lambda:pull(sfaces,A))
    compatibility={'U':sp.Matrix(named(sbody,'U'))-hold[:3,0],'E_W':named(sbody,'E_W')-hold[3]}
    J.emit('support-and-centre-unresolved',{'operator':mapped,'support':hold,'formalFaceCentreResidual':mapped[4]-hold[4],'bodyFaceCompatibility':compatibility,'sourceAssignmentCentreZero':mapped[4],'centreZeroIndependentEvidence':False,'backgroundTractionMeasure':'TRACTION_MEASURE_UNSUPPLIED','totalR0_c':'UNSUPPLIED','incrementalSupportLaw':'UNSUPPLIED','physicalPairing':'ORIGINAL_FORCE_PAIRING_UNRESOLVED','symmetricCentreForceFixedByThickness':False,'unknownHigherE':'UNSUPPLIED sigma2 and higher; product remainders do not supply it'})
    # A literal polynomial product ledger. No second-order physical input is invented.
    aa=sp.symbols('ledger_A0:3');tt=sp.symbols('ledger_t0:3');rr=sp.symbols('ledger_r0:3');z=sp.Symbol('ledger_wave_epsilon')
    product=sum(aa[k]*z**k for k in range(3))*sum(tt[k]*z**k for k in range(3))*sum(rr[k]*z**k for k in range(3))
    power=J.stage('new-wave-product-ledger',{'areaCoefficients':aa,'tractionCoefficients':tt,'rateCoefficients':rr,'product':product,'stationaryRate0':0,'scalarComponent':True},lambda:{'order0':sp.expand(product).coeff(z,0),'order1Stationary':sp.expand(product.subs(rr[0],0)).coeff(z,1),'order2Stationary':sp.expand(product.subs(rr[0],0)).coeff(z,2),'forceOrder2':aa[0]*tt[2]+aa[1]*tt[1]+aa[2]*tt[0],'unknownPhysicalSecondOrderAssigned':False})
    J.zero('wave-order1',power['order1Stationary'],aa[0]*tt[0]*rr[1])
    J.zero('wave-order2',power['order2Stationary'],aa[0]*tt[0]*rr[2]+aa[0]*tt[1]*rr[1]+aa[1]*tt[0]*rr[1])
    # Formal controls prove coefficient sensitivity, not nonzero physical power.
    def control(name,baseline,mutant):
        def work():
            move=sp.cancel(mutant-baseline);num,den=sp.fraction(move);variables=sorted(num.free_symbols,key=str);poly=sp.Poly(num,*variables) if variables else None
            nonzero=bool(poly.is_zero is False) if poly is not None else bool(num!=0)
            return {'movement':move,'numerator':num,'denominator':den,'formalPolynomialNonzero':nonzero,'physicalNonzeroValueClaimed':False}
        v=J.stage(name,{'baseline':baseline,'mutant':mutant,'interpretation':'FORMAL_COEFFICIENT_SENSITIVITY_ONLY'},work)
        require(v['formalPolynomialNonzero'],'applicable formal control '+name)
    control('duplicate-thickness-work-control',pull(faces,sp.S.One)[3]*et,(pull(faces,sp.S.One)[3]+e)*et)
    corrupt=[(sign,sp.Tuple(*t[:3],-t[3] if sign==1 else t[3])) for sign,t in faces]
    control('reverse-upper-centre-control',pull(faces,sp.S.One)[4],pull(corrupt,sp.S.One)[4])
    corrupt_hold=[(sign,sp.Tuple(*t[:3],2*t[3] if sign==1 else t[3])) for sign,t in sfaces]
    control('support-symbol-coefficient-control',pull(sfaces,sp.S.One)[4],pull(corrupt_hold,sp.S.One)[4])
    return {'executionStatus':'COMPLETED_FORMAL_BACKGROUND_COMPARISONS_WITH_UNSUPPLIED_PHYSICAL_PAIRING','backgroundTractionMeasure':'TRACTION_MEASURE_UNSUPPLIED','totalR0_c':'UNSUPPLIED','physicalWorkClosed':False,'fullIncrementalCentreDriveComputed':False,'oldFunctionsCalled':False,'newFormalControls':3,'leakageFactor':None,'scienceAcceptance':False}

def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest)
    expected=[str(Path(__file__).resolve()),'--out',str(args.out),'--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(sys.argv==expected and gate['command'][-len(expected):]==expected,'exact invocation')
    require(str(args.out.resolve())==gate['outputDirectory'],'gate output')
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
        result=J.stage('background-force-comparisons',{'manifestSha256':gate['manifestSha256'],'buildReviewSha256':gate['buildReviewRecordSha256']},lambda:scientific_work(manifest,J,ns['decode']));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result);sys.stderr.write(result['traceback'])
    finally:
        if (args.out/'saved-copy-index.json').exists():
            for v in json.loads((args.out/'saved-copy-index.json').read_text()).values():pins[str(args.out/v['path'])]=v['sha256']
        post={}
        for path,expected_hash in pins.items():
            try:post[path]={'expected':expected_hash,'actual':sha(path)}
            except OSError as e:post[path]={'expected':expected_hash,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',post)
        if any(v['expected']!=v['actual'] for v in post.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-start,scientificAcceptance=False)
        save(args.out/'checks.json',J.encode(result) if J else result);sys.stdout.write((args.out/'checks.json').read_text())
    return code
if __name__=='__main__':sys.exit(main())
