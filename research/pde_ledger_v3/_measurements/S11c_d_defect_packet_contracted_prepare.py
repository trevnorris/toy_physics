"""Guarded new Gaussian eligibility, formula, panel and tail joins.

Old proof returns and operands are restored, not recomputed. All new identities
and actual geometry execute only from the guarded numerical worker preamble.
"""
import ast
from fractions import Fraction
import math
from pathlib import Path
import S11c_d_defect_packet_contracted_validation as V

ZERO={'text':'0','srepr':'Integer(0)'}
COMPONENTS=('NATIVE_MIXED_ITERATION','INHERITED_DIRECT_WHOLE_OFF_DIAGONAL')


def require(v,message):
    if v is not True:raise ValueError(message)


def prepare(raw,manifest,J,sp,C,G,N):
    get=raw.__getitem__;R=lambda x:C.restore_scalar(sp,x)
    def pair(z):
        re,im=sp.expand_complex(z).as_real_imag();require(re.is_Rational is True and im.is_Rational is True,'exact complex rational coefficient')
        return [str(re),str(im)]
    def exact(name,left,right,context):
        J.start(name,{'left':left,'right':right,'context':context,'newDerivation':True})
        res=left-right;J.emit(name+'-raw',{'residual':res});value=sp.cancel(sp.together(res));J.emit(name+'-decision',{'cancelled':value})
        require(value is sp.S.Zero,'exact new join '+name);J.finish({'residual':value})
    accepted=get('contraction/result-record.json')
    require(accepted['allChecksPassed'] and accepted['status']=='BOUNDED_FINITE_WINDOW_CONTRACTION_CERTIFICATE_ACCEPTED_NUMERICAL_EVALUATOR_PENDING','actual accepted contraction')
    for alias,receipt in manifest['savedInputs'].items():
        if alias.startswith('contraction/complete/'):
            old=accepted['records']['complete/'+alias.removeprefix('contraction/complete/')]
            require((old['sha256'],old['bytes'])==(receipt['sha256'],receipt['bytes']),'accepted contraction operand receipt '+alias)
    old=lambda name:get('contraction/complete/'+name+'.json')
    def inherited(name):
        args=old(name+'-input');ret=old(name+'-return')
        J.emit('restored-'+name,{'arguments':args,'return':ret,'completedFunctionCalled':False})
        require(ret.get('residual')==ZERO,'literal accepted contraction zero '+name)
        return args
    oldproof={name:inherited('new-full-factorization-'+name) for name in ('J','Dr','Dh','Dq')}
    inherited('new-wrong-root-mutant-factorization')
    contracts=get('numeric-source-contracts.json')
    for fragment in contracts['fragments'].values():C.fragment(fragment)
    J.emit('new-original-source-contracts',contracts)
    physical=get('preflight/physical-plan.json');origins=get('accepted-units/parameter-quantity-origins.json')
    J.emit('restored-physical-context',{'physical':physical,'origins':origins,'originalPhysical':get('saved/physical-input.json')})
    require(physical==get('saved/preflight/physical-plan.json')==origins['physicalPlan'],'same exact original packets')
    require(origins['context']==get('saved/native/binding.json') and origins['context']['physicalInput']==get('saved/physical-input.json'),'same original physical parameter and native bindings')
    require(physical['frequency']==3 and R(physical['cs'])==sp.sqrt(6)/2 and R(physical['kappa'])==sp.sqrt(595)/10,'same real frequency, cs and kappa')
    require(physical['s']==8 and [R(v) for v in physical['centers']]==[-sp.Rational(5,2),sp.Rational(5,2)],'same centers and width')
    require(R(origins['context']['numeric']['W_0'])==1 and R(origins['context']['numeric']['L_W'])==10,'native W/L quantities')
    class ExactParameters:
        j=sp.I
        @staticmethod
        def mpf(v):return sp.Integer(v)
    native_a=R(get('saved/inner/runtime-a-input.json')['left']);native_mu=R(get('saved/inner/runtime-mu-input.json')['left'])
    beta_value,cj_value,cd_value=N.Evaluator.params(None,ExactParameters())
    for name,new,old_value in [('beta',beta_value,native_a*native_mu),('Cj',cj_value,native_a*native_mu**2*10/4),('Cd',cd_value,(-sp.I*native_mu)*10/(4*sp.I))]:
        exact('new-implemented-parameter-'+name,new,old_value,{'actualOriginalA':get('saved/inner/runtime-a-input.json'),'actualOriginalMu':get('saved/inner/runtime-mu-input.json'),'nativeW':origins['context']['numeric']['W_0'],'nativeL':origins['context']['numeric']['L_W'],'newParametersCallableExecutedOnExactConstants':True})
    numeric_tree=ast.parse(Path(manifest['librarySource']).read_text())
    newq=[v for v in ast.walk(numeric_tree) if isinstance(v,ast.FunctionDef) and v.name=='q'];require(len(newq)==1,'one actual numerical outgoing branch')
    oldq=C.fragment(contracts['fragments']['inner-root'])
    newrad=next(v.value for v in newq[0].body if isinstance(v,ast.Assign));oldrad=next(v.value for v in oldq.body if isinstance(v,ast.Assign))
    require(ast.dump(newrad)==ast.dump(oldrad),'actual original radicand AST')
    mapped=ast.parse(ast.unparse(newq[0]).replace('rad','d'))
    require(ast.dump(next(v.value for v in mapped.body[0].body if isinstance(v,ast.Return)))==ast.dump(next(v.value for v in oldq.body if isinstance(v,ast.Return))),'same native positive outgoing branch AST')
    J.emit('new-outgoing-source-join',{'original':contracts['fragments']['inner-root'],'actualNewSource':ast.get_source_segment(Path(manifest['librarySource']).read_text(),newq[0]),'sameRadicandAndBranch':True})
    V.source_ast_joins(contracts,numeric_tree,ast.parse(Path(__file__).read_text()),C,J)
    joined=V.whole_densities(raw,manifest,J,sp,C,exact)
    fields=get('saved/pressure/fields.json');adapters=get('saved/preflight/numeric-factor-adapters.json')
    selected=get('saved/selected/pressure-addresses.json')['selected'];row=get('saved/inventory/THETA_BALANCE-ordered-addresses.json')
    require(len(fields)==34 and len(adapters['definitions'])==20 and len(row)==2652,'all inherited source and response operands')
    require(selected==[v for v in row if v['row']=='THETA_BALANCE' and v['jet']['channel']=='e_W'] and len(selected)==544,'complete original row selection')
    for fid,field in fields.items():
        prefix='saved/field/'+fid
        proof,ret=get(prefix+'-reconstruction-input.json'),get(prefix+'-reconstruction-return.json')
        J.emit('restored-field-'+fid,{'field':field,'polynomial':get(prefix+'-polynomial.json'),'arguments':proof,'return':ret})
        require(proof['left']==field['field'] and ret['cancelled']==ZERO,'old full polynomial reconstruction arguments and exact zero')
    for i in range(17):
        stem='saved/factors/address-full-factor-'+str(i);operands=get(stem+'-operands.json')
        for suffix in ('full-mapped-residual','normal-source-join'):
            args,ret=get(stem+'-'+suffix+'-input.json'),get(stem+'-'+suffix+'-return.json')
            J.emit('restored-factor-'+str(i)+'-'+suffix,{'operands':operands,'input':args,'return':ret})
            require(ret['cancelled']==ZERO,'inherited full expression exact zero; simplified sides may be identical')
    for label,adapter in adapters['definitions'].items():
        context=old('template-context-'+label);args=get('saved/preflight/numeric-factor-'+label+'-arguments.json')
        require(context['adapter']==adapter and context['actualArguments']==args and context['accepted']==get('accepted-units/complete-template-'+label+'.json'),'actual complete native template and unit join')
        J.emit('restored-template-'+label,context)
    entries=[];applicable=[];grade_locations=[];k=sp.Symbol('new_numeric_k',real=True);l=sp.Symbol('new_numeric_l',real=True)
    ku,gu,gv=sp.symbols('new_numeric_kernel new_numeric_Gu new_numeric_Gv')
    for index,a in enumerate(selected):
        ident=a['addressId'];ui=get('accepted-units/summand-'+str(ident)+'-input.json');ur=get('accepted-units/summand-'+str(ident)+'-return.json')
        ai,ar=old('address-'+str(ident)+'-input'),old('address-'+str(ident)+'-return');label='-'.join((a['face'],a['slot'],a['component']))
        operands={'address':a,'unitInput':ui,'unitReturn':ur,'contractionInput':ai,'contractionReturn':ar}
        J.emit('restored-address-'+str(ident),operands)
        require(ai['address']==ui['address']==ui['sourceTransportInput']['address']==a and ai['acceptedUnitArguments']==ui and ai['acceptedUnitReturn']==ur==ar['actualUnitsInherited'],'original source/grade/wave/unit/complete address operands')
        require(ai['adapter']==adapters['definitions'][label] and ai['originalFactorOperands']==get('saved/factors/'+a['fullFactorProof']['proof']+'-operands.json'),'native response factor and argument joins')
        require(ui['waveProof']['left']==a['waveMultiplier'] and ui['normal']==a['normalMultiplier'],'original wave and normal operands')
        require(ar['allOriginalJoins'] and ar['addressId']==ident and ar['explicitZero']==a['status'].startswith('EXACT_ZERO'),'accepted full address classification')
        require(adapters['addressJoins'][index]=={'addressId':ident,'adapter':label},'same complete adapter address map')
        for role in ('source','consumer'):require(fields[a[role+'Transform']['coefficientId']]['field']==a[role+'Field'],'actual address coefficient')
        grade_locations.append({'addressId':ident,'target':a['targetGrade'],'source':a['sourceGrade'],'response':a['responseGrade'],'consumer':a['consumerGrade'],'status':a['status']})
        if a['component'] not in COMPONENTS:continue
        applicable.append(a)
        if a['status'].startswith('EXACT_ZERO'):continue
        require(a['slot']=='pressure' and a['normalMultiplier']==a['normalOriginal']=={'text':'1','srepr':'Integer(1)'},'actual pressure-only normal selector')
        require(a['sourceGrade']==a['consumerGrade']==[0,0] and a['targetGrade']==a['responseGrade']==[1,1] and a['epsilonCount']==1,'actual independent grades/epsilon')
        require(a['deltaSupport']==['k=p','r=l'] and a['measure']=='dl dk' and a['responseMap']['flatSupport'] is None,'constant coefficient deltas leave k and l independent')
        require(a['sourceTransform']['transfer']=='k-p' and a['consumerTransform']['transfer']=='r-l' and a['responseInputDepth']=='q(k)' and a['responseOutputDepth']=='q(l)' and a['responseMap']['frequency']==3 and a['responseMap']['positiveRegulatorContinuation'] is False,'actual source/test transfer and real-frequency response arguments')
        values={};specs={}
        for role in ('source','consumer'):
            fid=a[role+'Transform']['coefficientId'];field=fields[fid];poly=get('saved/field/'+fid+'-polynomial.json');proof=get('saved/field/'+fid+'-reconstruction-input.json')
            J.emit('new-constant-'+str(ident)+'-'+role+'-input',{'actualField':field,'polynomial':poly,'oldProof':proof,'address':a})
            require(field['constant'] is True and poly['degree']==0 and len(poly['coefficients'])==1,'actual full constant polynomial eligibility')
            numerator,denominator=R(poly['polynomial']),R(poly['denominator'])
            require(denominator.is_Rational is True and denominator>0 and not numerator.free_symbols,'actual constant numerator and regular denominator')
            exact('new-constant-numerator-'+str(ident)+'-'+role,numerator,R(poly['coefficients'][0]),{'originalPolynomial':poly})
            value=numerator/denominator;values[role]=value
            exact('new-constant-quotient-'+str(ident)+'-'+role,value,R(a[role+'Field']),{'oldProof':proof,'originalAddress':a})
            exact('new-constant-value-'+str(ident)+'-'+role,value,R(a[role+'Transform']['constantValue']),{'actualField':field})
            interface=V.original_selector(ur,a,role,[pair(value)])
            J.emit('new-original-selector-'+str(ident)+'-'+role,{'actualUnitReturn':ur,'selectedInterface':interface,'address':a})
            specs[role]={'role':'X' if role=='source' else 'Y','coefficientsAscending':[pair(value)],'argumentDerivative':interface['argumentDerivative'],'originalInterface':interface,'center':'-5/2' if role=='source' else '5/2','spatialOrders':a['jet']['spatialOrders'] if role=='source' else [0,0,0],'timeOrder':a['jet']['timeOrder'] if role=='source' else 0,'argument':'k' if role=='source' else '-l','carrier':'p0' if role=='source' else '-p0','width':'8','normalization':'1/(2*pi)' if role=='source' else '1','coefficientProof':{'alias':'saved/field/'+fid+'-reconstruction-input.json','returnAlias':'saved/field/'+fid+'-reconstruction-return.json'},'nativeFieldId':fid,'unit':ur['X' if role=='source' else 'Y']}
        n,n2,n3=a['jet']['spatialOrders'];nt=a['jet']['timeOrder'];require(n in (0,2),'only actual n=0,2 subset')
        P=(-3*sp.I)**nt*(sp.I/5)**n2*(sp.I/10)**n3
        wave=R(a['waveMultiplier']);symbols=C.named_symbols(wave);require(set(symbols)<={'composition_p'},'original wave symbols')
        mapped=wave.xreplace({s:k for s in symbols.values()})
        exact('new-native-wave-on-delta-'+str(ident),mapped,P*(sp.I*k)**n,{'deltaSupport':a['deltaSupport'],'originalWave':a['waveMultiplier'],'originalUnitProof':ui['waveProof']})
        b,c=values['source'],values['consumer'];alpha=b*c*P*sp.I**n
        exact('new-consumer-epsilon-binding-'+str(ident),R(a['consumerOriginal']),R(a['epsilon'])*c,{'actualAddress':a,'unitInput':ui,'unitReturn':ur,'epsilonCount':a['epsilonCount']})
        V.response_selection(raw,J,sp,C,exact,a,ai['adapter'],ar,joined)
        exact('new-family-complete-summand-'+str(ident),b*c*mapped*gu*gv*ku/(2*sp.pi),alpha*k**n*gu*gv*ku/(2*sp.pi),{'completeTemplate':ai['adapter'],'contractionUnits':ar.get('contractionUnits'),'actualIndependentDepths':['q(k)','q(l)'],'nativeNormal':a['normalMultiplier'],'onlyOrdinaryJ':a['component']==COMPONENTS[0]})
        for tag in ('packet_J','packet_D'):
            proof=old('new-native-injection-'+label+'-'+tag+'-return');require(proof['residual']==ZERO,'actual accepted complete native injection proof')
        entry={'addressId':ident,'face':a['face'],'component':a['component'],'n':n,'b':pair(b),'c':pair(c),'Pj':pair(P),'alpha':pair(alpha),'Xmultiplier':pair(b*P*sp.I**n),'Ymultiplier':pair(c),'primitives':['J'] if a['component']==COMPONENTS[0] else ['Dr','Dh','Dq'],'grade':a['targetGrade'],'specs':specs,'unit':ur['total'],'fullOriginalOperands':{'copyAlias':'contraction/complete/address-'+str(ident)+'-input.json','returnAlias':'contraction/complete/address-'+str(ident)+'-return.json','literalSelectors':{'address':['address'],'units':['acceptedUnitReturn'],'template':['adapter'],'factor':['originalFactorOperands']}},'epsilonAlreadyExtracted':True}
        J.emit('new-eligible-'+str(ident),entry);entries.append(entry)
    require(len(applicable)==64 and len(entries)==20 and sum(v['component']==COMPONENTS[0] for v in entries)==10,'actual applicable/live census')
    require(sum(a['status']=='EXACT_ZERO_SOURCE_JET' for a in applicable)==12 and sum(a['status']=='EXACT_ZERO_CONSUMER' for a in applicable)==32,'all44 explicit zero classifications')
    ids=[e['addressId'] for e in entries];require(min(e['addressId'] for e in entries if e['component']==COMPONENTS[1])==8347 and min(e['addressId'] for e in entries if e['component']==COMPONENTS[0] and e['n']==2)==8350,'actual applicable native control selection')
    J.emit('new-complete-grade-coverage',grade_locations)
    # New Gaussian derivative/moment identities. The transform theorem and
    # integration by parts on Schwartz Gaussians are assessed analytic facts.
    w,nu,p0,s=sp.symbols('new_w new_nu new_p0 new_s',real=True);Q=sp.S.One;moments=[sp.S.One,-sp.I*s*s*nu]
    for r in range(1,2):moments.append(-sp.I*s*s*nu*moments[r]+r*s*s*moments[r-1])
    for n in range(3):
        if n:Q=sp.expand(sp.diff(Q,w)+(sp.I*p0-w/s**2)*Q)
        if n in {e['n'] for e in entries}:
            poly=sp.Poly(Q,w);moment=sum(poly.nth(i)*moments[i] for i in range(n+1))
            exact('new-centered-Gaussian-n'+str(n),moment,(sp.I*(nu+p0))**n,{'Q':Q,'momentsOverM0':moments,'w':'x-xu','nu':'k-p0','normalization':'1/(2pi)','sNonzero':True})
    exact('new-Y-phase',(p0-l)*sp.Rational(-5,2),(l-p0)*sp.Rational(5,2),{'argument':'-l','carrier':'-p0','center':'5/2','extraConjugation':False})
    # The same implemented combination is joined to each inherited full
    # contraction formula, with original assumptions and symbol objects intact.
    definitions=old('new-contraction-definitions');outer={name:R(v) for name,v in definitions['outerIntegrands'].items()};symbols=C.named_symbols(sp.Add(*outer.values()))
    z={name:symbols['contraction_'+name] for name in definitions['integrands']}
    cm,qm,beta,aa,mu,W,L=[symbols['contract_'+name] for name in ('m','qm','beta','a','mu','W','L')]
    generated=N.combine(cm,qm,beta,aa*mu**2*W*L/4,(-sp.I*mu)*W*L/(4*sp.I),z)
    for name,value in generated.items():exact('new-numeric-combine-source-'+name,value,outer[name],{'acceptedDefinitions':definitions,'oldExactProof':oldproof.get(name),'newCallable':'combine; no numerical integration yet'})
    J.emit('new-Gaussian-envelope-transport',{'oldBounds':get('preflight/Fourier-envelope-constants.json'),'actualEligibleSources':entries,'originalFourierSource':contracts['fragments']['fourier-product'],'reason':'Actual constant quotient/specs and physical centered Q identities above describe exactly the same untruncated source/test Fourier integrals as the old strip proof. The assessed Gaussian identity evaluates those integrals. Therefore their saved CX exp(-5|k-p0|), CY exp(-5|l-p0|) bounds transport by equality, not a new magnitude-only assumption. No old integral or envelope constructor is called.','carrierDomain':'|p0|<3','scope':'real omega3 and original centers/width, argumentDerivative0','sourceOfAnalyticInequalities':'accepted original strip proof, not independently rederived'})
    # Restore explicit absolute numerator envelopes so a signed-sum bound is
    # never used to bound an individual direct summand.
    pieces={name:get('absolute-bounds/'+name+'.json') for name in ('D-reflected-numerator-envelope','D-height-numerator-envelope','J-numerator-envelope','global-q-triangle-gap')}
    for name,rec in pieces.items():
        J.emit('restored-absolute-'+name,rec)
        require(all(R(v).is_nonnegative is True for v in rec['coefficients']),'inherited nonnegative numerator coefficients')
    common=C.named_symbols(R(pieces['D-reflected-numerator-envelope']['larger']))
    ak,al,at=[common['weak_abs_'+n] for n in ('k','l','t')];Pabs=1+ak+al
    expected={
      'D-reflected-numerator-envelope':(2*Pabs**2*(1+at),ak*(2*al+at)),
      'D-height-numerator-envelope':(18*Pabs**2*(1+at),ak*(2*ak+at)+16*(1+ak)**2),
      'J-numerator-envelope':(2*Pabs**2*(1+at)**2,(ak+at)*(2*ak+at)),
      'global-q-triangle-gap':((ak+4)**2,ak**2+sp.Rational(453,50))}
    for name,(larger,smaller) in expected.items():
        rec=pieces[name]
        # Only new argument/primitive joins; do not reconstruct the accepted gap.
        exact('new-absolute-larger-join-'+name,larger,R(rec['larger']),{'savedCertificate':rec,'actualOriginalDensityProofs':oldproof,'absoluteVariableContext':common})
        exact('new-absolute-smaller-join-'+name,smaller,R(rec['smaller']),{'savedCertificate':rec,'actualOriginalDensities':joined['originalDensities'],'outgoingRatioBound':'|qi/(qi+qh)|<=1; |qs+qo|>=|qs|; |qi|^2<=16(1+|k|)^2','heightAndQuadraticAddedBeforeBounding':True})
    tail_context=get('preflight/tail-bound-derivation.json');domain=get('saved/pressure/global-parameter-domain.json');profile=get('saved/pressure/profile-envelope.json')
    require(tail_context['inheritedBounds']==get('saved/pressure/whole-envelopes.json') and profile['productBound']=='121 exp(-|t|)','actual positive global envelope ancestry')
    J.emit('new-per-primitive-absolute-tail-transport',{'pieces':pieces,'domain':domain,'profile':profile,'actualDensityProofs':oldproof,'lemma':'Outgoing quadrant: |qi/(qh+qi)|<=1 and |qs+qo|>=|qs|. With |Cd|<=1 and |E|>=b*, reflected numerator <=2 P^2(1+|t|), combined height+quadratic numerator <=18 P^2(1+|t|). Each added direct term is therefore bounded by the positive SUM 18 P^2(1+|t|)*(1/|qs|+1/|qh|)*|A1 A2|/b*^2. This is a triangle majorant before integration; it is not |sum D|. The accepted shift-root integral yields the inherited D envelope. If T>=K+4 and kappa<3, both depth moduli exceed1 outside the middle window, giving 36*121/b*^2 times the degree1 exponential tail.','eachPrimitiveUsesFullPositiveDEnvelope':True,'notRecomputedGlobalProof':True})
    tails=[];base=get('preflight/tail-plan.json');byid={r['addressId']:r for r in base['allAddresses']};bounds=get('preflight/Fourier-envelope-constants.json')['bounds']
    require((base['K'],base['T'])==(27,122),'accepted base window')
    def exponential_moment(n,start,rate):
        # Exact polynomial prefactor of integral_start^infty x^n exp(-rate*x) dx.
        return sum(sp.Rational(math.factorial(n),math.factorial(j))*start**j/rate**(n-j+1) for j in range(n+1))
    def weighted_tail(degree,start,rate=5):
        # For integer rate*start, e^(-rate*start) <= 2^(-rate*start).
        require(isinstance(start,int) and start>=0 and isinstance(rate,int) and rate>0,'integer tail plan')
        rate=sp.Integer(rate)
        return 2*sp.Rational(1,2)**int(rate*start)*sum(math.comb(degree,n)*exponential_moment(n,sp.Integer(start),rate) for n in range(degree+1))
    bstar=R(tail_context['b']);F3=R(tail_context['fullExponentialMoments']['3']);ordinary=R(tail_context['ordinaryKernelUpper'])
    require(bool(bstar==sp.Rational(3000,11101) and bstar>0 and F3>0 and ordinary>0),'same positive original tail constants')
    for e in entries:
        ident=e['addressId'];cx,cy=R(bounds[str(ident)]['X']),R(bounds[str(ident)]['Y']);require(bool(cx>0 and cy>0),'positive actual source/test bound')
        for K,T in ((27,122),(29,124)):
            if K==27:
                r=byid[ident];outerbound=R(r['outer']);middlebound=R(r['middle'])
                J.emit('restored-base-tail-'+str(ident),r)
            else:
                J.start('new-enlarged-tail-'+str(ident),{'K':K,'T':T,'CX':cx,'CY':cy,'originalConstants':tail_context,'address':e,'originalSource':contracts['fragments']['tail-contributions']})
                b=bstar;E30=sp.Integer(3)**30;Cordinary=ordinary;Tlim=T
                F={2:R(tail_context['fullExponentialMoments']['2']),3:F3}
                outerbound=2*Cordinary*cx*cy*E30*F[3]*weighted_tail(3,K)
                if e['component']==COMPONENTS[0]:
                    middlebound=4*cx*cy*E30*F[3]**2*(sp.Rational(4,5)*121/b**3)*weighted_tail(2,Tlim,1)
                    middlebound+=(4/b)*cx*cy*E30*F[2]**2*sp.Rational(55,3)*sp.Rational(1,2)**Tlim
                else:
                    middlebound=4*cx*cy*E30*F[3]**2*(36*121/b**2)*weighted_tail(1,Tlim,1)
                J.emit('new-enlarged-tail-'+str(ident)+'-decision',{'outer':outerbound,'middle':middlebound,'originalDomain':T>=K+4,'includesPositiveHOvercount':e['component']==COMPONENTS[0],'sourceFunctionCalled':False})
                require(bool(T>=K+4 and outerbound>=0 and middlebound>=0),'positive new tail domain')
                J.finish({'outer':outerbound,'middle':middlebound})
            for primitive in e['primitives']:tails.append({'addressId':ident,'primitive':primitive,'K':K,'T':T,'outer':outerbound,'middle':middlebound,'total':outerbound+middlebound,'baselineOnly':True,'includesExtraHOvercount':primitive=='J'})
    for row in tails:
        epsilon=sp.Rational(1,80*10**11)
        J.emit('new-tail-address-budget-'+str(row['addressId'])+'-'+row['primitive']+'-'+str(row['K']),{'row':row,'epsilon':epsilon,'noOtherRowDonation':True})
        require(bool(row['total']>=0 and row['total']<epsilon),'actual positive per-address per-primitive tail allocation')
    for K in (27,29):
        total=sum(v['total'] for v in tails if v['K']==K);J.emit('new-positive-tail-total-'+str(K),{'rows':[v for v in tails if v['K']==K],'sum':total,'budget':sp.Rational(1,10**11),'noCancellation':True})
        require(bool(total<sp.Rational(1,10**11)),'positive all-address/primitive tail allocation')
    plans={};receipts={};kap=G.Quad(0,1,Fraction(119,20))
    for K,T in ((27,122),(29,124)):
        for carrier in (kap,kap*0):
            key=str(K)+'/'+str(carrier.b)
            label_alias='new-label-transport-'+('matching' if carrier==kap else 'zero')+'-K'+str(K)+'-square-full-arrangement'
            labels=old(label_alias)
            require(labels['originalPlan']==get('geometry/'+label_alias.removeprefix('new-label-transport-')+'.json'),'actual same window/carrier original labels')
            require({v['primitive'] for v in labels['transport']}=={'J','Dr','Dh','Dq'},'all four primitive label transports')
            for primitive in ('J','Dr','Dh','Dq'):
                J.emit('new-selected-label-transport-'+key.replace('/','-')+'-'+primitive,{'originalPlan':labels['originalPlan'],'transport':[v for v in labels['transport'] if v['primitive']==primitive],'sourceAlias':label_alias,'nativeMethod':'contracted variables; original external labels are retained as provenance, not an old quadrature request'})
            J.start('new-geometry-'+key.replace('/','-'),{'K':K,'T':T,'kappa':kap,'carrier':carrier,'width':8,'length':10,'inheritedLabels':labels})
            plan=G.mesh(kap,K,T,carrier,8,10);rec=J.emit('new-geometry-'+key.replace('/','-')+'-full',G.packed(plan));J.finish({'audit':True,'slabs':len(plan['slabs']),'fullWings':True})
            for side in (-1,1):J.emit('new-wing-control-'+key.replace('/','-')+'-'+str(side),G.packed(G.wing_control(plan,side)))
            plans[key]=plan;receipts[key]=rec
    J.emit('new-nested-allocation-lemma',{'epsilon':sp.Rational(1,80*10**11),'pointDensityTarget':'epsilon*(1+1/abs(q(m)))/(16*(3*M+4))','momentumNondimensionalizedInNativeUnit':True,'bound':'integral[-M,M]1/abs(q(m)) dm=pi+2 acosh(M/kappa)<4+M for kappa>2 and M>=149; hence integral w<3M+4','discreteGuardSeparate':'actual absolute K/G or GL weighted nested indicators <=epsilon/4','empiricalIndicatorsNotRigorous':True})
    return {'entries':entries,'geometryReceipts':receipts,'scope':manifest['scope'],'ruleReceipts':{r:manifest['savedInputs']['saved/rules/'+n+'.json'] for r,n in [('A24','A-GL24'),('A48','A-GL48'),('B50','B-G7-K15')]},'sourceManifestSha256':manifest['selfSourceIdentity'],'originalProofDependence':True,'inferredGammaUnits':True},plans,tails
