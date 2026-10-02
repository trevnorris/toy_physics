#!/usr/bin/env python3
"""Bounded real-frequency source/consumer address inventory; no integral/solve."""
import argparse
import ast
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import resource
import shutil
import sys
import time
import traceback
from types import SimpleNamespace

ROOT = Path('/var/projects/toy_physics')
THREADS = ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS = ('require','sha','save','replace_json','posthash_records','containment',
           'Journal','decode','one_symbol','named','function_source','assignment_source',
           'expanded_sinh_arguments')
TEXT_HELPERS = ('require','text_sha','literal_record','tuple_arguments','literal_key',
                'named','selected_case')
G = ((0,0),(1,0),(0,1),(1,1))
SLOTS = ('delta_p_plus','delta_p_minus','d_w_delta_p_plus','d_w_delta_p_minus')
ROWS = ('U0','U1','U2','THETA_BALANCE','E_W_BALANCE')


def require(v, message):
    if v is not True:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def save(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


def definitions(text, names):
    nodes = [n for n in ast.parse(text).body
             if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in names]
    require({n.name for n in nodes} == set(names), 'exact helper census')
    return ast.Module(body=nodes, type_ignores=[])


def verify_gate(path, manifest_path, manifest):
    gate = json.loads(Path(path).read_text())
    require(gate['status']=='READY_FOR_ONE_SOURCE_COMPOSITION_INSTRUMENT','gate status')
    require(gate['workerSha256']==sha(__file__) and
            gate['manifestSha256']==sha(manifest_path), 'worker/manifest')
    require(gate['sourcePins']==manifest['sourcePins'], 'source census')
    for p,h in gate['sourcePins'].items():
        require(sha(p)==h,'source '+p)
    for key in ('sharedGuard','supervisor','launcher','buildReviewRecord','authority'):
        require(sha(gate[key])==gate[key+'Sha256'],'gate '+key)
    require(gate['launcher']==manifest['launcher'], 'launcher route')
    review=json.loads(Path(gate['buildReviewRecord']).read_text())
    require(review['independentBuildClearance'] is True and
            review['correctedMethodAssessed'] is True and review['allChecksPassed'] is True,
            'substantive implementation and corrected-method assessment')
    for engine in ('claude','grok'):
        require(review['reports'][engine]['literalVerdict']==
                'CLEAR FOR THIS BOUNDED SOURCE-COMPOSITION BUILD','literal build verdict')
    for key in ('workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256','launcherSha256'):
        require(review[key]==gate[key], 'review/gate '+key)
    require(review['methodSha256']==sha(manifest['methodPath']), 'exact corrected method')
    authority=json.loads(Path(gate['authority']).read_text())
    require(authority['boundedInstrumentAuthorized'] is True and
            authority['automaticScientificRetry'] is False, 'standing bounded authority')
    require(gate['durationLimits'] is None and gate['scientificRunsAuthorized']==1 and
            gate['scope']==manifest['scope'], 'one bounded no-deadline job')
    return gate


def triples(target):
    return [(a,b,c) for a,b,c in itertools.product(G,repeat=3)
            if tuple(a[i]+b[i]+c[i] for i in (0,1))==target]


def jet_spec(name):
    original=name
    if name.startswith('grad_theta_'):
        name='theta_d'+name.rsplit('_',1)[1]
    match=re.fullmatch(r'(u_[123]|theta|e_W)((?:_t{1,2})?(?:_?d[123])*)', name)
    if not match:
        return None
    channel,suffix=match.groups()
    return {'name':original,'channel':channel,'timeOrder':2 if '_tt' in suffix else 1 if '_t' in suffix else 0,
            'spatialOrders':[len(re.findall('d'+str(i),suffix)) for i in (1,2,3)],
            'baseDimension':[1,0,0] if channel.startswith('u_') else [0,0,0]}


def broad_row_census(text):
    """Full source-text partition; no symbolic evaluation of the local slab part."""
    tree=ast.parse(text,mode='eval').body
    children=tree.args if isinstance(tree,ast.Call) and getattr(tree.func,'id',None)=='Add' else [tree]
    result=[]
    for index,node in enumerate(children):
        fragment=ast.get_source_segment(text,node);hits=[]
        for a in ast.walk(node):
            if (isinstance(a,ast.Call) and getattr(a.func,'id',None) in ('Symbol','Function')
                and a.args and isinstance(a.args[0],ast.Constant) and isinstance(a.args[0].value,str)
                and ('delta_p' in a.args[0].value or 'd_w_' in a.args[0].value)):
                hits.append({'constructor':a.func.id,'name':a.args[0].value})
        result.append({'childIndex':index,'constructorText':fragment,
                       'sha256':hashlib.sha256(fragment.encode()).hexdigest(),'hits':hits})
    return result


def polynomial_terms(expression, eta, sigma):
    out={}
    for term in sp.Add.make_args(sp.expand(expression)):
        pd=term.as_powers_dict();powers=[];coefficient=term
        for x in (eta,sigma):
            n=pd.get(x,sp.S.Zero)
            require(n.is_integer is True and n.is_nonnegative is True, 'polynomial grade domain')
            powers.append(int(n));coefficient=coefficient/x**n
        require(not coefficient.has(eta,sigma), 'grade-independent coefficient')
        key=tuple(powers);out[key]=out.get(key,sp.S.Zero)+coefficient
    return {k:sp.cancel(v) for k,v in out.items()}


def quotient_recurrence(numerator, denominator, grades, cancel):
    """Exact multivariate quotient recurrence; no full coefficient-domain Poly."""
    result={}
    for a,b in sorted(grades,key=lambda g:(sum(g),g)):
        remainder=numerator.get((a,b),0)
        for (i,j),value in denominator.items():
            if (i or j) and i<=a and j<=b:
                remainder-=value*result.get((a-i,b-j),0)
        result[a,b]=cancel(remainder/denominator[(0,0)])
    return result


def grade_split(J,name,full,eta,sigma,zero_saved,nonzero):
    J.emit(name+'-operands',{'full':full,'savedZero':zero_saved,'grades':[eta,sigma]})
    numerator,denominator=sp.fraction(sp.cancel(full))
    n=polynomial_terms(numerator,eta,sigma);d=polynomial_terms(denominator,eta,sigma)
    den0=d.get((0,0),sp.S.Zero)
    nonzero(J,name+'-regular-denominator',den0)
    table=quotient_recurrence(n,d,tuple(itertools.product(range(3),repeat=2)),sp.cancel)
    retained=sum(table[g]*eta**g[0]*sigma**g[1] for g in G)
    raw=numerator-denominator*retained;terms=polynomial_terms(raw,eta,sigma)
    J.emit(name+'-split',{'numerator':numerator,'denominator':denominator,'denominatorAtZero':den0,
        'retained':{str(g):table[g] for g in G},'excludedPure':{str(g):table[g] for g in ((2,0),(0,2))},
        'fullHigherRemainder':full-retained,'quotientRingNumeratorRemainder':raw,
        'zeroGradeReused':True,'nativeShapeCoefficientsCalled':False})
    J.zero(name+'-saved-zero',table[0,0],zero_saved)
    for g in G:J.zero(name+'-quotient-'+str(g[0])+str(g[1]),terms.get(g,sp.S.Zero),sp.S.Zero)
    return {g:table[g] for g in G}


def run_science(manifest,J,helpers,text_helpers,nonzero):
    cache={};copies={};used=set()
    for name,r in manifest['savedFiles'].items():
        src=Path(r['path']);require(sha(src)==r['sha256'] and src.stat().st_size==r['bytes'],'saved '+name)
        dst=J.out/'saved'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
        require(sha(dst)==r['sha256'],'copied '+name)
        copies[name]={'source':str(src),'path':str(dst.relative_to(J.out)),'sha256':r['sha256'],'bytes':r['bytes']}
    save(J.out/'saved-copy-index.json',copies)
    def load(name):
        if name not in cache:cache[name]=helpers['decode'](json.loads((J.out/'saved'/name).read_text()))
        used.add(name);return cache[name]
    one=helpers['one_symbol'];c2=Path(manifest['c2Source']).read_text()
    routes=load('native/function-routes.json')['routes']
    for r in routes:
        require(sha(ROOT/r['sourcePath'])==r['sourceSha256'],'route source')
        nodes=[n for n in ast.parse(Path(ROOT/r['sourcePath']).read_text()).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==r['name']]
        require(len(nodes)==1 and ast.get_source_segment(Path(ROOT/r['sourcePath']).read_text(),nodes[0])==r['text'],'route text '+r['name'])
    context=load('consumer/binding-context.json');eta,sigma=context['independentGrades'];eps=context['epsilon']
    require(context['frequency']==3, 'real-frequency3 binding')
    physical=json.loads(Path(manifest['physicalInput']).read_text())
    require(context['physicalInput']==physical,'actual physical input join')
    def bind(value):
        syms=value.atoms(sp.Symbol);named={s.name:s for s in syms}
        require(len(named)==len(syms),'no ambiguous raw symbol assumptions')
        dm={s:context['densityMap'][s.name] for s in syms if s.name in context['densityMap']}
        value=value.xreplace(dm)
        for _ in range(len(context['profileEqualities'])+2):
            mapping={s:context['profileEqualities'][s.name] for s in value.atoms(sp.Symbol)
                     if s.name in context['profileEqualities']}
            changed=value.xreplace(mapping)
            if changed==value:break
            value=changed
        else:raise ValueError('background substitution cycle')
        mapping={s:context['numeric'][s.name] for s in value.atoms(sp.Symbol) if s.name in context['numeric']}
        require(not any(s in mapping for s in (eta,sigma,eps)),'independent grades/epsilon')
        return value.xreplace(mapping)
    J.emit('actual-binding-context',{'restored':context,'fixedFrequency':3,'densityRestored':True,
        'regulatedCompositionClaim':False,'noSigmaEtaIdentification':True,'csOnlyInDepth':True})
    require(context['numeric']['omega']==3 and context['numeric']['W_0']==1 and
            context['numeric']['Lambda_X_0']==0,'material/frequency settings')
    native=load('native/consumer-census.json')
    for path,rec in native['sourcePins'].items():require(sha(path)==rec['sha256'],'original source pin')
    # Broad full-row metadata scan reuses text only; old symbolic row assembly is never called.
    txt,provenance=text_helpers['literal_record'](Path(manifest['bSource']),'slab_operator')
    _,selected=text_helpers['selected_case'](txt,('LAB_HELD','RHO4_CONSTANT'))
    value=text_helpers['named'](selected,'VALUE')
    require(hashlib.sha256(selected.encode()).hexdigest()==native['slabConsumers']['caseConstructorSha256'],'selected slab case source')
    row_texts={}
    for name in ('U_BODY_BALANCE','THETA_BALANCE','E_W_BALANCE'):
        row=text_helpers['named'](value,name);expanded=text_helpers['named'](row,'EXPANDED')
        if name=='U_BODY_BALANCE':
            row_texts.update({'U'+str(i):t for i,t in enumerate(text_helpers['tuple_arguments'](expanded))})
        else:row_texts[name]=expanded
    consumers={};consumer_full={};census={};unit_by_key={}
    units=load('consumer-finish/consumer-unit-joins.json')
    for u in units:unit_by_key[u['row'],str(u['slot'])]=u
    for rowname in ROWS:
        entries=broad_row_census(row_texts[rowname]);unknown=[h for e in entries for h in e['hits'] if h['constructor']!='Symbol' or h['name'] not in SLOTS]
        J.emit(rowname+'-full-native-partition',{'source':provenance,'fullConstructor':row_texts[rowname],
            'children':entries,'unknownPressureAtoms':unknown,'localChildIndices':[e['childIndex'] for e in entries if not e['hits']],
            'pressureChildIndices':[e['childIndex'] for e in entries if e['hits']],
            'localPartUnchanged':True,'localPartRestoredAsScience':False})
        require(not unknown,'broad pressure/normal derivative census')
        old=load('raw/'+rowname+'-executed-native-census.json')
        hits=[e for e in entries if e['hits']]
        require(len(entries)==old['totalChildren'],'all native children')
        require([(e['childIndex'],e['sha256']) for e in hits]==[(e['childIndex'],e['constructorSha256']) for e in old['selected']],'actual selected child identities')
        # Bind only old selected pressure children; non-pressure children stay original source data.
        row_saved=load('consumer/'+rowname+'-consumer-input.json');row_raw=row_saved['raw'];row_bound=row_saved['bound']
        actual=sum(sp.sympify(e['constructorText']) for e in hits)
        J.zero(rowname+'-pressure-source-join',actual,row_raw)
        J.zero(rowname+'-pressure-bound-join',bind(row_raw),row_bound)
        coeff=row_saved['slotCoefficients'];symbols={s.name:s for s in row_bound.atoms(sp.Symbol)}
        for label in ('plus','minus'):
            slotname='delta_p_'+label;atom=symbols.get(slotname,sp.Symbol(slotname))
            ablated=row_bound.subs(atom,0);removed=row_bound-ablated
            J.emit(rowname+'-'+label+'-native-pressure-slot-ablation',{'baseline':row_bound,'ablated':ablated,'removed':removed,'slot':atom})
            J.zero(rowname+'-'+label+'-native-pressure-ablation-join',removed,coeff[slotname]*atom)
        J.zero(rowname+'-affine-pressure-reconstruction',row_bound,sum(coeff[n]*symbols.get(n,sp.Symbol(n)) for n in SLOTS))
        consumers[rowname]={};consumer_full[rowname]=coeff;census[rowname]=entries
        for slot in SLOTS:
            expr=coeff[slot]
            require(not any(s.name in SLOTS for s in expr.atoms(sp.Symbol)), 'affine independent-slot coefficients')
            oldsplit=load('consumer/'+rowname+'-'+slot+'-grade-split.json')
            J.zero(rowname+'-'+slot+'-full-coefficient',expr,oldsplit['full'])
            table=grade_split(J,rowname+'-'+slot,expr,eta,sigma,oldsplit['zeroGrade'],nonzero)
            consumers[rowname][slot]=table
            for g,v in table.items():
                J.zero(rowname+'-'+slot+'-epsilon-'+str(g[0])+str(g[1]),v,eps*sp.cancel(v/eps))
                require(eps not in sp.cancel(v/eps).free_symbols,'epsilon exactly once in nonzero consumer')
            if rowname.startswith('U'):
                require(all(v==0 for v in table.values()),'native U pressure absence')
            else:
                u=unit_by_key[rowname,slot]
                J.zero(rowname+'-'+slot+'-native-unit-coefficient',bind(u['coefficient']),expr)
                for i in range(3):J.zero(rowname+'-'+slot+'-unit-'+str(i),u['total'][i],u['expected'][i])
    sources={};source_jets={};source_full={};x=sp.Symbol('composition_x',real=True)
    w=(1+sp.tanh(x/10))/2;m=(1-sp.tanh(x/10)**2)/3
    def profile(expr):
        mapping={}
        for s in expr.free_symbols:
            hit=re.fullmatch(r'([wm])1_profile((?:_?d[123])*)',s.name)
            if hit:
                base,suffix=hit.groups();directions=re.findall('d([123])',suffix)
                mapping[s]=sp.S.Zero if any(d!='1' for d in directions) else 10**len(directions)*sp.diff(w if base=='w' else m,x,len(directions))
        result=sp.cancel(expr.xreplace(mapping))
        J.emit('profile-map-'+str(len(J.artifacts)),{'original':expr,'map':[[a,b] for a,b in mapping.items()],
            'oneDimensional':result,'L':10,'independentSigma':True,'profileTransformEvaluated':False})
        require(result.free_symbols<={x},'coefficient field fully bound')
        return result
    for face in ('plus','minus'):
        inp=load('consumer/'+face+'-source-input.json');saved=load('consumer/'+face+'-source-grade-split.json')
        J.zero(face+'-saved-source-parts',saved['full'],inp['parts']['velocity']+inp['parts']['chemical'])
        J.zero(face+'-own-velocity-normalization',bind(inp['nativeVelocity']/eps),inp['velocityAmplitude'])
        J.zero(face+'-own-velocity-coefficient',inp['parts']['velocity'],inp['velocityCoefficient']*inp['velocityAmplitude'])
        J.zero(face+'-own-chemical-coefficient',inp['parts']['chemical'],inp['chemicalCoefficient']*inp['chemicalAmplitude'])
        J.zero(face+'-native-composed-source',inp['combined'],saved['full'])
        require(eps not in saved['full'].free_symbols,'source epsilon stripped exactly once')
        require(not any(s.name in ('omega','c_s0','W_bg','rho_br_bg_rho4_constant') for s in saved['full'].free_symbols),'numeric and live density binding complete')
        table=grade_split(J,face+'-source',saved['full'],eta,sigma,saved['zeroGrade'],nonzero)
        source_full[face]=saved['full'];sources[face]=table;source_jets[face]={}
        alljets=sorted([a for a in saved['full'].free_symbols if jet_spec(a.name)],key=str)
        for g,expr in table.items():
            jetrows=[];reconstruction=sp.S.Zero
            for atom in alljets:
                coefficient=sp.cancel(sp.diff(expr,atom));reconstruction+=coefficient*atom
                require(not any(j in coefficient.free_symbols for j in alljets),'linear wave-jet coefficient')
                field=profile(coefficient)
                spec=jet_spec(atom.name);orders=spec['spatialOrders'];dim=[spec['baseDimension'][0]-sum(orders),-spec['timeOrder'],0]
                jetrows.append({'atom':atom,'spec':spec,'coefficient':coefficient,'field':field,'jetDimension':dim,
                    'coefficientDimensionFromSource':[1-dim[0],-1-dim[1],-dim[2]],
                    'status':'EXACT_ZERO' if field==0 else 'FORMAL_COEFFICIENT_AVAILABLE'})
            J.emit(face+'-source-jets-'+str(g[0])+str(g[1]),{'grade':g,'source':expr,'jets':jetrows,'reconstruction':reconstruction,
                'dimensionScope':'Native un-restricted wave_jet dimensions and inherited normalized source dimension; material units remain in saved coefficient provenance'})
            J.zero(face+'-source-jet-reconstruction-'+str(g[0])+str(g[1]),expr,reconstruction)
            source_jets[face][g]=jetrows
    J.zero('both-face-source-equality',source_full['plus'],source_full['minus'])
    # Response pieces are copied, not reconstructed by the old response producers.
    rc=load('reference/retained-response-census.json');objects=[rc[n] for n in ('flat','heightCoefficient','slopeCoefficient','mixedIteration','taggedTotalMixed')]
    omega=one(objects,'reference_unrestricted_frequency');qi=one(objects,'reference_qi');qo=one(objects,'reference_qo')
    k=one(objects,'reference_k');l=one(objects,'reference_l');p,r=sp.symbols('composition_p composition_r',real=True)
    cs=sp.Symbol('composition_cs',positive=True);beta=omega/(10-sp.I*omega)
    wholeD=one(objects,'certified_closed_direct_whole_convolution')
    beta3=beta.subs(omega,3);R=lambda q:q/(q+beta3)
    # Distinct bare factor, closed-density denominator and whole tag. No density-minus-factor test.
    direct=load('direct/closed-density.json');O=one([direct['density']],'grazing_unrestricted_frequency')
    dq_i=one([direct['density']],'grazing_qi');dq_o=one([direct['density']],'grazing_qo')
    J.zero('beta-identity',direct['beta'].subs(O,omega),beta)
    E=1/((qi+beta3)*(qo+beta3));Rprod=R(qi)*R(qo)
    J.emit('typed-direct-objects',{'bareClosureFactor':Rprod,'factoredDensityExternalDenominator':E,
        'closedDensity':direct,'wholeTag':wholeD,'typedRelation':'Rprod = qi*qo*E; entire closed density is not Rprod',
        'barePlaceholderIsWholeTag':False,'multiplyWholeTagByResolvents':False})
    J.zero('typed-factor-relation',Rprod,qi*qo*E)
    # Inherit complete full-expression equalities; don't repeat their cancellation.
    for face,sign in (('plus',1),('minus',-1)):
        for op in ('actual-closed-reference','actual-closed-normalJet'):
            prior_in=load('direct/'+face+'-'+op+'-original-input.json')
            prior_out=load('direct/'+face+'-'+op+'-canonical-return.json')
            J.emit(face+'-'+op+'-restored-proof',{'input':prior_in,'return':prior_out,'recomputed':False,
                 'scope':'Completed full-expression identity; actual current right-operand join below'})
            require(prior_out['cancelled']==0,'literal completed density return')
            target=direct['density'].subs(O,3)*(sp.I*sign*dq_o if op.endswith('normalJet') else 1)
            J.sinh_zero(face+'-'+op+'-current-operand-join',prior_in['right'],target)
        raw=load('raw/'+face+'-closure-operands.json');before=load('raw/'+face+'-closed-before-cancel.json')
        bare=one([raw['matrix']],'{}_raw_direct'.format(face))
        mapping={s:{'q_i':qi,'q_o':qo}.get(s.name,s) for s in raw['factor'].free_symbols}
        J.zero(face+'-bare-linear-factor',sp.diff(raw['closed'][0,2],bare).xreplace(mapping),Rprod)
        J.zero(face+'-restored-closure-factor',raw['factor'].xreplace(mapping),Rprod)
        J.zero(face+'-restored-reference-factor',before['reference'].xreplace(mapping),Rprod)
        J.zero(face+'-restored-normal-factor',before['jet'].xreplace(mapping),sp.I*sign*qo*Rprod)
        # New independent-D native assignment join. The adapter is deliberately not kernel_apply.
        D=sp.Symbol(face+'_independent_bare_D');P=eta*sigma*Rprod*D
        trace=load('reference/'+face+'-new-native-trace.json');slot=load('reference/'+face+'-final-native-slot-routing.json')
        require(trace['newTrace'][0,0]==1 and trace['newTrace'][0,2]==0 and slot['valueCoefficient']==1,'own face trace diagonal/direct slot')
        T=trace['newTrace'];delta_matrix=sp.zeros(3);delta_matrix[0,2]=P
        J.zero(face+'-unchanged-native-T01-direct-subtraction',(T*delta_matrix)[0,2],P)
        H=slot['affineHeight'];solution=slot['equationSolution'];target=one([solution],face+'_physical_target');jet=one([solution],face+'_jet_slot')
        normal=sp.Symbol('composition_NORMAL',real=True);calls=[]
        def adapter(inputs,diagonal,off_diagonal,source,kout,kin,second=sp.S.Zero):
            require(diagonal==0 and second==0 and source==1,'direct-only formal adapter')
            require(kout==(l,) and kin==(k,), 'adapter ordered arguments')
            calls.append({'diagonal':diagonal,'offDiagonal':off_diagonal,'source':source,'kout':kout,'kin':kin,'second':second})
            return off_diagonal
        ns={'sp':sp,'face':sign,'qo':qo,'NORMAL':normal,'reference':sp.Rational(sign,2),
            'reference_matrix':sp.Matrix([[0,P],[0,0]]),'jet_diagonal':sp.S.Zero,'jet_second':sp.S.Zero,
            'inputs':None,'composed_source':sp.S.One,'ko':(l,),'ki':(k,),'kernel_apply':adapter,
            'trace_map':{'REFERENCE_VALUE_SOLVE':solution,'PHYSICAL_PRESSURE_TARGET':target},'pressure':P,'jet_slot':jet}
        fragments={name:helpers['assignment_source'](c2,'build_face',name) for name in ('extension','jet_transfer','normal_jet','reference_pressure')}
        J.emit(face+'-direct-adapter-input',{'D':D,'physicalInjection':P,'referenceTransferInjection':P,
            'bothDifferOnlyBeyondRetainedRectangle':True,'traceHeight':slot['savedHeight'],'traceSolution':solution,
            'nativeAssignments':fragments,'newAddress':'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL','nativeFunctionsCalled':False})
        for name in fragments:exec(fragments[name],ns)
        height=slot['savedHeight'];unprojected=ns['reference_pressure'].subs(H,height)
        pt=polynomial_terms(unprojected,eta,sigma);projected=sum(pt.get(g,sp.S.Zero)*eta**g[0]*sigma**g[1] for g in G)
        J.emit(face+'-direct-adapter-return',{'calls':calls,'jet':ns['normal_jet'],'unprojectedReference':unprojected,
            'excluded':unprojected-projected,'retainedReference':projected,'algebraicRoutingOnly':True,
            'noMiddleIntegral':True,'physicalSourceNotExcited':True})
        J.zero(face+'-direct-normal-routing',ns['normal_jet'],sp.I*sign*qo*P)
        J.zero(face+'-direct-reference-routing',projected,P)
        require(set(pt)<={(1,1),(2,1)}, 'direct trace only excluded21')
        J.zero(face+'-saved-flat-diagonal-jet',rc['jetKernels'][face]['flat'],sp.I*sign*qo*rc['flat'].subs(qi,qo))
        flat_good=rc['jetKernels'][face]['flat'].subs({omega:3,qo:sp.Rational(3,2)},simultaneous=True)
        flat_bad=(sp.I*sign*qo*rc['flat']).subs({omega:3,qo:sp.Rational(3,2),qi:2},simultaneous=True)
        J.emit(face+'-unbound-flat-pole-control',{'baseline':flat_good,'corrupt':flat_bad,'notPhysicalOnDiagonal':True,'independentDepths':[2,sp.Rational(3,2)]})
        nonzero(J,face+'-unbound-flat-pole-movement',flat_bad-flat_good)
        for key in ('height','slope'):
            J.zero(face+'-saved-'+key+'-jet',rc['jetKernels'][face][key],sp.I*sign*qo*rc[key+'Coefficient'])
        J.zero(face+'-saved-mixed-sum',rc['jetKernels'][face]['mixed'],sp.I*sign*qo*(rc['mixedIteration']+wholeD))
    J.zero('saved-total-mixed-sum',rc['taggedTotalMixed'],rc['mixedIteration']+wholeD)
    # Whole definitions and maps retain their bound variable roles; no Integral constructor is called.
    tagdefs={}
    t,td=sp.symbols('composition_middle_transfer composition_direct_transfer',real=True)
    qfun=sp.Function('common_outgoing_q')
    # The value is a declared outgoing sheet, not a mode solver or a grazing assignment.
    depth_variable=sp.Symbol('composition_depth_momentum',real=True)
    rad=9/cs**2-sp.Rational(1,20)-depth_variable**2
    saved_sheet=load('raw/physical-sheet.json')
    require(saved_sheet['frequency']==3 and saved_sheet['edge']==[sp.Rational(1,5),sp.Rational(1,10)],'actual inherited depth sheet frequency/edges')
    old_q=saved_sheet['input'];qm_names={'k':depth_variable,'increment_effective_bulk_speed':cs}
    require({a.name for a in old_q.free_symbols}==set(qm_names),'saved input sheet map census')
    qm_map={a:qm_names[a.name] for a in old_q.free_symbols}
    qdefinition=old_q.subs(qm_map,simultaneous=True)
    require(isinstance(qdefinition,sp.Piecewise) and len(qdefinition.args)==3 and qdefinition.args[-1].expr is sp.nan,'unassigned grazing branch')
    J.zero('mapped-positive-depth-branch',qdefinition.args[0].expr,sp.sqrt(rad))
    J.zero('mapped-decaying-depth-branch',qdefinition.args[1].expr,sp.I*sp.sqrt(-rad))
    J.emit('inherited-depth-map',{'saved':old_q,'map':[[a,b] for a,b in qm_map.items()],'new':qdefinition})
    new_direct_map={'grazing_unrestricted_frequency':sp.Integer(3),'grazing_output':l,'k':k,
        'grazing_transfer':td,'grazing_qi':qfun(k),'grazing_qo':qfun(l),
        'grazing_qh':qfun(k+td),'grazing_qs':qfun(l-td)}
    require({a.name for a in direct['density'].free_symbols}==set(new_direct_map),'complete direct map symbol census')
    actual_direct_map={a:new_direct_map[a.name] for a in direct['density'].free_symbols}
    mapped_direct=direct['density'].subs(actual_direct_map,simultaneous=True)
    # Join the inherited full density and argument map as complete objects.
    old_join=load('reference/restored-direct-argument-join.json')
    J.zero('saved-direct-density-provenance',old_join['savedDensity']['density'],direct['density'])
    old_mapped=direct['density'].subs(dict(old_join['map']),simultaneous=True)
    J.zero('saved-direct-argument-provenance',old_mapped,old_join['transported'])
    jd=load('reference/right-height-PV-operands.json')['density']
    new_iter_map={'reference_unrestricted_frequency':sp.Integer(3),'reference_k':k,'reference_l':l,
        'reference_t':t,'reference_qi':qfun(k),'reference_qo':qfun(l),'reference_qm':qfun(k+t)}
    require({a.name for a in jd.free_symbols}==set(new_iter_map),'complete iteration density symbol census')
    actual_iter_map={a:new_iter_map[a.name] for a in jd.free_symbols}
    J.emit('actual-whole-density-argument-maps',{'direct':{'old':direct['density'],'map':[[a,b] for a,b in actual_direct_map.items()],
        'mapped':mapped_direct,'boundVariable':td,'assumptions':[[a,a.assumptions0] for a in actual_direct_map]},
        'iteration':{'old':jd,'map':[[a,b] for a,b in actual_iter_map.items()],'mapped':jd.subs(actual_iter_map,simultaneous=True),'boundVariable':t},
        'depthFunction':qfun(depth_variable),'depthDefinition':qdefinition,'cs':cs,
        'noValueAtGrazing':True,'noGlobalComposedDomainClaim':True})
    for name,path in [('H','reference/left-height-subtracted-PV.json'),('Jwhole','reference/right-height-PV-operands.json'),('Dwhole','direct/closed-density.json')]:
        tagdefs[name]={'savedDefinition':load(path),'sha256':copies[path]['sha256'],
            'freeMomenta':['l','k'],'boundVariable':'td' if name=='Dwhole' else 't',
            'middleMomentum':'k+td' if name=='Dwhole' else 'k+t','reflectedMomentum':'l-td' if name=='Dwhole' else None,
            'valueEvaluated':False,'globalDomain':'UNRESOLVED outside certified response k/l[-3,3]'}
    tagdefs['H']['freeMomenta']=['l-k'];tagdefs['H']['middleMomentum']=None
    J.emit('whole-tag-definitions',tagdefs)
    direct_map=load('reference/restored-direct-argument-join.json')
    J.emit('inherited-direct-map',{'saved':direct_map,'realFrequency':3,'argumentMapRequired':True,'newConvolution':False})
    contract=load('raw/native-fourier-contract.json');edge_contract=load('raw/edge-delta-reduction.json')
    require(contract['profileForwardPower']==-3 and contract['sourceInversePower']==-3 and contract['sourceForwardPower']==0,'native Fourier normalization')
    require(edge_contract['plainReducedMiddleMeasure'] is True,'saved plain edge-reduced measure')
    J.emit('fourier-and-unit-provenance',{'contract':contract,'edgeReduction':edge_contract,
        'sourceUnit':load('consumer-finish/units-and-scope.json')['normalizedSourceDimension'],
        'consumerUnitJoins':units,'nativeWaveJetSource':next(a['text'] for a in routes if a['name']=='wave_jet'),
        'nativeInputsSource':next(a['text'] for a in routes if a['name']=='Inputs'),
        'oneDimensionalForwardNormalization':'1/(2*pi)','measures':['dl','dl dk'],'evaluatedTransform':False})
    # Fourier routing keeps whole coefficient functions as tags, not products of transforms.
    ft={}
    def transform(field,transfer):
        key=hashlib.sha256(sp.srepr(field).encode()).hexdigest()
        if key not in ft:ft[key]={'field':field,'constant':x not in field.free_symbols,
            'definition':'hat b(s)=(1/(2*pi))*integral exp(-i*s*x)*b(x) dx; distribution not evaluated',
            'constantRule':'b*delta(s)','productTransformedAsWhole':True}
        return {'coefficientId':key,'transfer':transfer,'constantValue':field if x not in field.free_symbols else None}
    def wave(spec,momentum):
        out=(-sp.I*3)**spec['timeOrder']
        for n,mom in zip(spec['spatialOrders'],(momentum,sp.Rational(1,5),sp.Rational(1,10))):out*= (sp.I*mom)**n
        return out
    wholeH=one(objects,'whole_height_slope_convolution');wholeJ=one(objects,'whole_iterated_density_integral')
    Htag=sp.Function('Hwhole')(l-k,sp.Integer(1),sp.Integer(10))
    base_components={ (0,0):[('NATIVE_FLAT',rc['flat'].subs({omega:3,qi:qo},simultaneous=True))],
        (1,0):[('NATIVE_HEIGHT',rc['heightCoefficient'].subs(omega,3))],
        (0,1):[('NATIVE_SLOPE',rc['slopeCoefficient'].subs(omega,3))] }
    def components_for(face,grade):
        if grade!=(1,1):return base_components[grade]
        signature=(l,k,sp.Integer(3),cs,sp.Rational(1,5),sp.Rational(1,10),sp.Integer(1),sp.Integer(10))
        Jtag=sp.Function('Jwhole_'+face)(*signature);Dtag=sp.Function('Dwhole_'+face)(*signature)
        iteration=rc['mixedIteration'].subs({omega:3,wholeH:Htag,wholeJ:Jtag},simultaneous=True)
        return [('NATIVE_MIXED_ITERATION',iteration),('INHERITED_DIRECT_WHOLE_OFF_DIAGONAL',Dtag)]
    coverage=[];addresses=[];consumer_fields={}
    for rowname in ROWS:
        for face,sign in (('plus',1),('minus',-1)):
            for slotkind,prefix in (('pressure','delta_p_'),('normal','d_w_delta_p_')):
                slot=prefix+face;consumer_fields[rowname,slot]={g:profile(sp.cancel(c/eps)) for g,c in consumers[rowname][slot].items()}
                for target_grade in G:
                    for a,b,c in triples(target_grade):
                        cv=consumers[rowname][slot][a];cf=consumer_fields[rowname,slot][a];source_list=source_jets[face][c]
                        core={'row':rowname,'face':face,'slot':slotkind,'targetGrade':target_grade,'consumerGrade':a,'responseGrade':b,'sourceGrade':c}
                        zero='EXACT_ZERO_CONSUMER' if cf==0 else 'EXACT_ZERO_SOURCE' if all(v['field']==0 for v in source_list) else None
                        coverage.append({**core,'status':zero or 'ADDRESSED','sourceJets':len(source_list)})
                        for item in source_list:
                            for component,rv in components_for(face,b):
                                reason=zero or ('EXACT_ZERO_SOURCE_JET' if item['field']==0 else None)
                                reduced=b==(0,0)
                                addr={**core,'component':component,'jet':item['spec'],'sourceField':item['field'],
                                    'consumerField':cf,'epsilonCount':1,'epsilon':eps,'consumerOriginal':cv,
                                    'sourceTransform':transform(item['field'],'l-p' if reduced else 'k-p'),
                                    'consumerTransform':transform(cf,'r-l'),'waveMultiplier':wave(item['spec'],p),
                                    'responseCoefficient':rv,'normalMultiplier':sp.I*sign*qo if slotkind=='normal' else sp.S.One,
                                    'responseInputDepth':'q(l)' if reduced else 'q(k)','responseOutputDepth':'q(l)',
                                    'freeMomenta':['p','l','r'] if reduced else ['p','k','l','r'],
                                    'measure':'dl' if reduced else 'dl dk','deltaSupport':['k=l'] if reduced else [],
                                    'status':reason or 'FORMAL_ADDRESS_AVAILABLE','globalComposition':'UNRESOLVED',
                                    'nativeRowChildHashes':[e['sha256'] for e in census[rowname] if any(h['name']==slot for h in e['hits'])],
                                    'sourceOperandSha256':copies['consumer/'+face+'-source-grade-split.json']['sha256'],
                                    'responseCensusSha256':copies['reference/retained-response-census.json']['sha256'],
                                    'wholeValueEvaluated':False}
                                if x not in item['field'].free_symbols:addr['deltaSupport'].append('l=p' if reduced else 'k=p')
                                if x not in cf.free_symbols:addr['deltaSupport'].append('r=l')
                                addresses.append(addr)
        J.emit(rowname+'-ordered-addresses',[a for a in addresses if a['row']==rowname])
    require(len(coverage)==len(ROWS)*2*2*16,'complete 16-triple coverage per face/slot/row')
    J.emit('grade-coverage',coverage);J.emit('whole-coefficient-transforms',ft)
    # Independent formal expansion of ordered multiplication, with coefficient fields noncommuting.
    C={g:sp.Symbol('Mconsumer'+str(g),commutative=False) for g in G}
    F={g:sp.Symbol('Fresponse'+str(g),commutative=False) for g in G}
    S={g:sp.Symbol('MsourceWave'+str(g),commutative=False) for g in G}
    formal=sp.expand(sp.prod(sum(v*eta**g[0]*sigma**g[1] for g,v in tab.items()) for tab in (C,F,S)))
    for g in G:
        expected=sum(C[a]*F[b]*S[c] for a,b,c in triples(g))
        actual=formal.coeff(eta,g[0]).coeff(sigma,g[1]);residual=sp.expand(actual-expected)
        J.emit('ordered-grade-algebra-'+str(g[0])+str(g[1]),{'actual':actual,'expected':expected,'expandedResidual':residual})
        require(residual==0,'noncommuting ordered grade reconstruction')
    J.emit('ordered-grade-meaning',{'newAlgebraNotNativeExecution':True,'multiplicationOrder':['consumer','response','source then wave derivative'],
        'actualAddressesFileSuffix':'-ordered-addresses.json','nativeAffineRowJoinedBeforePlaceholders':True})
    controls=[]
    qpoint={sp.Rational(3,2):sp.Integer(2),sp.Integer(2):sp.Rational(3,2),sp.Rational(30,13):sp.Rational(25,26)}
    for momentum,depth in qpoint.items():J.zero('physical-control-depth-'+str(len(J.artifacts)),depth**2,9/sp.Rational(10,7)-sp.Rational(1,20)-momentum**2)
    def sensitivity(name,baseline,corrupt,context_record):
        movement=sp.cancel(corrupt-baseline)
        J.emit(name+'-control-operands',{'baselineCoefficient':baseline,'corruptCoefficient':corrupt,'movement':movement,
             'context':context_record,'formalTagCoefficientOnly':True,'physicalConvolutionNotEvaluated':True})
        cert=nonzero(J,name+'-movement',movement);controls.append({'name':name,'movement':movement,'certificate':cert,'formalTagCoefficientOnly':True})
    for rowname in ('THETA_BALANCE','E_W_BALANCE'):
        for face,sign in (('plus',1),('minus',-1)):
            cp=consumer_fields[rowname,'delta_p_'+face][0,0];cn=consumer_fields[rowname,'d_w_delta_p_'+face][1,0]
            require(x not in cp.free_symbols,'constant pressure consumer')
            nonzero(J,rowname+'-'+face+'-pressure-applicability',cp)
            nonzero(J,rowname+'-'+face+'-normal-applicability',cn.subs(x,0))
            src10=next((a for a in source_jets[face][1,0] if a['field']!=0 and a['spec']['spatialOrders'][0]>0),None)
            require(src10 is not None,'actual nonzero source10 spatial jet')
            nonzero(J,rowname+'-'+face+'-source10-profile-applicability',src10['field'].subs(x,0))
            # Entire source transform is retained as a tag. Baseline/corruption are its coefficients.
            flat3=sp.Rational(3,10)/(sp.Rational(3,2)+beta3)
            wp=wave(src10['spec'],sp.Rational(3,2));wk=wave(src10['spec'],sp.Integer(2))
            context10={'sourceJet':src10,'sourceTransform':transform(src10['field'],'2-3/2'),'support':'k=l=r=2,p=3/2','csSquared':sp.Rational(10,7)}
            sensitivity(rowname+'-'+face+'-source-p-to-k',cp*flat3*wp,cp*flat3*wk,context10)
            sensitivity(rowname+'-'+face+'-omit-source10',cp*flat3*wp,sp.S.Zero,context10)
            src00=next((a for a in source_jets[face][0,0] if a['field']!=0 and x not in a['field'].free_symbols and wave(a['spec'],sp.Rational(3,2))!=0),None)
            require(src00 is not None,'nonzero constant source00')
            sourceval=src00['field']*wave(src00['spec'],sp.Rational(3,2));nonzero(J,rowname+'-'+face+'-source00-applicability',sourceval)
            slope=(sp.Rational(3,10)*sp.Rational(3,2))/((sp.Rational(3,2)+beta3)*(2+beta3))
            tags={'consumerTransform':transform(cn,'30/13-2'),'slopeProfile':'j(2-3/2)','support':'p=k=3/2,l=2,r=30/13','sourceJet':src00,'csSquared':sp.Rational(10,7)}
            baseline=sp.I*sign*sp.Rational(3,2)*slope*sourceval
            sensitivity(rowname+'-'+face+'-normal-q-l-to-r',baseline,sp.I*sign*sp.Rational(25,26)*slope*sourceval,tags)
            sensitivity(rowname+'-'+face+'-omit-consumer10',baseline,sp.S.Zero,tags)
            directcoef=cp*sourceval
            for label,factor in (('omit-direct',0),('double-direct',2)):
                sensitivity(rowname+'-'+face+'-'+label,directcoef,factor*directcoef,
                    {'tag':'Dwhole_'+face+'(2,3/2)','sourceJet':src00,'grades':[[0,0],[1,1],[0,0]],'externalResolventMultiplier':1})
            if face=='minus':
                normflat=sp.I*sign*2*sp.Rational(3,10)/(2+beta3)*sourceval
                sensitivity(rowname+'-lower-jet-sign',normflat,-normflat,
                    {'consumerTransform':transform(cn,'30/13-3/2'),'support':'p=k=l=3/2,r=30/13','grade':[1,0]})
    J.emit('responsive-formal-controls',controls)
    J.emit('applicability-obligations',{'globalComposition':'UNRESOLVED','globalTestSpace':'UNRESOLVED','composedGrazingLimit':'UNRESOLVED',
        'sourceRegulatorContinuation':'NOT_CLAIMED','internalCertifiedDomain':{'cs':[1,2],'k/l':[-3,3]},
        'hiddenCutoff':False,'oldCutoffs4and6Reused':False,'currentOrPowerComputed':False,'savedWorkReplayed':False})
    J.emit('consumed-source-index',{'files':sorted(used),'copies':copies})
    return {'executionStatus':'COMPLETED_BOUNDED_SOURCE_COMPOSITION_INVENTORY','bothFaces':True,'gradeTriplesPerRowFaceSlot':16,
        'coverageEntries':len(coverage),'addressEntries':len(addresses),'formalResponsiveControls':len(controls),
        'realFrequency':3,'scientificAcceptance':False,'integralsEvaluated':False,'finiteSolves':0,'productionChanges':False,
        'globalComposition':'UNRESOLVED','lossClaim':False,'independentPhysicsValidationByControls':False}


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest)
    pins={**manifest['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False)
    J=None;result={};code=1;started=time.monotonic()
    try:
        ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(definitions(Path(manifest['helperSource']).read_text(),HELPERS),'unchanged-inert-helpers','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        ns.update(sp=sp,Str=Str);J=ns['Journal'](args.out)
        tex={'ast':ast,'hashlib':hashlib,'Path':Path}
        exec(compile(definitions(Path(manifest['textHelperSource']).read_text(),TEXT_HELPERS),'unchanged-text-helpers','exec'),tex)
        exact={'sp':sp,'require':require}
        exec(compile(definitions(Path(manifest['exactHelperSource']).read_text(),('exact_nonzero_number',)),'unchanged-exact-nonzero','exec'),exact)
        result=J.stage('source-consumer-composition',{'manifestSha256':gate['manifestSha256'],'reviewSha256':gate['buildReviewRecordSha256']},
            lambda:run_science(manifest,J,ns,tex,exact['exact_nonzero_number']));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False}
        save(args.out/'failure.json',result)
    finally:
        if (args.out/'saved-copy-index.json').exists():
            copies=json.loads((args.out/'saved-copy-index.json').read_text());pins.update({str(args.out/v['path']):v['sha256'] for v in copies.values()})
        records={}
        for path,expected in pins.items():
            try:records[path]={'expected':expected,'actual':sha(path),'error':None}
            except OSError as e:records[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',records)
        if any(v['expected']!=v['actual'] for v in records.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code


if __name__=='__main__':sys.exit(main())
