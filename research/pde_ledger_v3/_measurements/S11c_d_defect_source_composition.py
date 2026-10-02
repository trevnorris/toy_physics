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


def verify_helper_paths(gate):
    require(gate['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py'), 'actual guard route')
    require(gate['supervisor']==str(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),
            'actual supervisor route')


def verify_invocation(args, gate, argv):
    require(args.out.resolve()==Path(gate['outputDirectory']).resolve(),'gate output route')
    expected=[str(Path(__file__).resolve()),'--out',str(args.out),
              '--inputs',str(args.inputs),'--gate',str(args.gate)]
    require(list(argv)==expected,'actual worker argv')
    require(gate['command'][-len(expected):]==expected,'gate worker command tail')


def verify_gate(path, manifest_path, manifest):
    gate = json.loads(Path(path).read_text())
    verify_helper_paths(gate)
    require(gate['status']=='READY_FOR_ONE_SOURCE_COMPOSITION_INSTRUMENT','gate status')
    require(gate['workerSha256']==sha(__file__) and
            gate['manifestSha256']==sha(manifest_path), 'worker/manifest')
    require(gate['sourcePins']==manifest['sourcePins'], 'source census')
    for p,h in gate['sourcePins'].items():
        require(sha(p)==h,'source '+p)
    for key in ('sharedGuard','supervisor','launcher','buildReviewRecord','authority'):
        require(sha(gate[key])==gate[key+'Sha256'],'gate '+key)
    require(gate['launcher']==manifest['launcher'], 'launcher route')
    require(gate['buildReviewRecord']==manifest['reviewRecordWillBe'],'manifest review route')
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


def join_source_input(J, face, inp, chemical, normalization, bind, one, epsilon):
    """New argument joins on restored operands; no prior source function call."""
    velocity_atom=one([inp['raw']], 's11cc1_V_lab_held_'+face)
    chemical_atom=one([inp['raw']], 's11cc1_mu_theta_lab_held_'+face)
    density_atom=one([inp['raw']], 'rho_br_bg_rho4_constant')
    stage2={velocity_atom:inp['velocityAmplitude'], chemical_atom:chemical['amplitude']}
    J.emit(face+'-native-source-join-operands', {'saved':inp,'nativeChemical':chemical,
        'inheritedNormalization':normalization,'epsilon':epsilon,
        'stage2Map':[[a,b] for a,b in stage2.items()],'liveDensityAtom':density_atom,
        'previousSourceFunctionCalled':False})
    identified=inp['raw'].subs(stage2, simultaneous=True)
    bound=bind(identified)
    J.emit(face+'-native-source-join-bound', {'identified':identified,'bound':bound,
        'target':inp['combined'],'densityBindingContext':'actual-binding-context'})
    J.zero(face+'-native-chemical-amplitude-join',inp['chemicalAmplitude'],chemical['amplitude'])
    J.zero(face+'-native-chemical-epsilon-join',bind(chemical['raw'][1]/epsilon),chemical['amplitude'])
    J.zero(face+'-inherited-velocity-normalization-join',inp['velocityCoefficient'],normalization['velocityCoefficient'])
    J.zero(face+'-raw-stage2-live-density-join',bound,inp['combined'])


def speed_symbol_inventory(operands):
    records=[]
    for address,expression in sorted(operands.items()):
        names=sorted({s.name for s in expression.free_symbols})
        hits=[n for n in names if n.lower().startswith('c_s') or 'speed' in n.lower()
              or n.lower()=='cs' or n.lower().startswith('cs_')]
        records.append({'address':address,'freeSymbolNames':names,'speedSymbols':hits})
    return records


def native_profile_scale_rule(text):
    classes=[n for n in ast.parse(text).body if isinstance(n,ast.ClassDef) and n.name=='Inputs']
    require(len(classes)==1,'native Inputs source identity')
    methods=[n for n in classes[0].body if isinstance(n,ast.FunctionDef) and n.name=='at_source']
    require(len(methods)==1,'native at_source source identity')
    expected=ast.parse("if base in ('w1_profile','m1_profile'):\n    value*=self.values['L_W']**len(indices)").body[0]
    matches=[n for n in ast.walk(methods[0]) if isinstance(n,ast.If)
             and ast.dump(n.test)==ast.dump(expected.test)]
    require(len(matches)==1 and ast.dump(matches[0])==ast.dump(expected),'native profile jet scale rule')
    return {'source':ast.get_source_segment(text,matches[0]),
            'rule':'L_W**number_of_native_spatial_indices','functionExecuted':False}


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
            'spatialOrders':[len(re.findall('d'+str(i),suffix)) for i in (1,2,3)]}


def address_sum_for(addresses,row,face,grade):
    selected=[a for a in addresses if a['row']==row and a['face']==face and a['targetGrade']==grade]
    total=sum(a['consumerOriginal']*a['responsePlaceholder']*a['sourceOriginal']*a['sourceAtom'] for a in selected)
    return selected,total


def complete_factor_map(response_symbols,normal_symbols,roles,recorded_pairs):
    expected={a:b for a,b in roles.items() if a in response_symbols}
    require(len(recorded_pairs)==len(expected) and dict(recorded_pairs)==expected,'complete source-bound response map roles')
    return {a:b for a,b in roles.items() if a in response_symbols | normal_symbols}


def select_certified_candidate(candidates):
    return next((v for v in candidates if v['eligible'] is True and
                 v['zero'] is False and v['finite'] is True),None)


def epsilon_power(exact_zero):
    return 0 if exact_zero else 1


def native_jet_dimensions(text):
    """Read the unrestricted field dimension rule, without executing wave_jet."""
    calls=[n for n in ast.walk(ast.parse(text)) if isinstance(n,ast.Call)
           and isinstance(n.func,ast.Name) and n.func.id=='field' and len(n.args)==3
           and isinstance(n.args[0],ast.BinOp)]
    require(len(calls)==1, 'unique unrestricted field rule')
    rule=calls[0].args[2]
    expected=ast.parse("(1,0,0) if base.startswith('u_') else (0,0,0)",mode='eval').body
    require(ast.dump(rule)==ast.dump(expected),'native field dimension decision')
    return {'velocity':list(ast.literal_eval(rule.body)),
            'scalar':list(ast.literal_eval(rule.orelse)),
            'source':ast.get_source_segment(text,calls[0]),
            'sourceSha256':hashlib.sha256(text.encode()).hexdigest()}


def broad_row_census(text):
    """Partition literal constructor text and refuse hidden/dynamic name routes."""
    tree=ast.parse(text,mode='eval').body
    children=tree.args if isinstance(tree,ast.Call) and getattr(tree.func,'id',None)=='Add' else [tree]
    result=[]
    for index,node in enumerate(children):
        fragment=ast.get_source_segment(text,node);hits=[];unsupported=[]
        for call in ast.walk(node):
            if not isinstance(call,ast.Call):continue
            name=getattr(call.func,'id',getattr(call.func,'attr',None))
            if name in ('Symbol','Function','symbols'):
                literal=bool(call.args and isinstance(call.args[0],ast.Constant) and isinstance(call.args[0].value,str))
                if not isinstance(call.func,ast.Name) or name=='symbols' or not literal:
                    unsupported.append({'text':ast.get_source_segment(text,call),'reason':'unsupported constructor spelling or nonliteral name'})
                elif 'delta_p' in call.args[0].value or 'd_w_' in call.args[0].value:
                    hits.append({'constructor':name,'name':call.args[0].value})
            if name=='Add' and not isinstance(call.func,ast.Name):
                unsupported.append({'text':ast.get_source_segment(text,call),'reason':'nonliteral Add partition'})
        text_count=sum(fragment.count(v) for v in ('delta_p','d_w_'))
        hit_count=sum(h['name'].count(v) for h in hits for v in ('delta_p','d_w_'))
        result.append({'childIndex':index,'constructorText':fragment,
            'sha256':hashlib.sha256(fragment.encode()).hexdigest(),'hits':hits,
            'rawPressureSubstringCount':text_count,'accountedSubstringCount':hit_count,
            'unsupportedConstructors':unsupported,'completeNameCoverage':text_count==hit_count and not unsupported})
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
    # Inventory native operands before any density/profile/numeric binding.
    speed_operands={'chemical/raw':load('consumer/native-chemical-amplitude.json')['raw'][1]}
    for face in ('plus','minus'):
        source_input=load('consumer/'+face+'-source-input.json')
        for key in ('raw','nativeVelocity'):
            speed_operands[face+'/source/'+key]=source_input[key]
    for rowname in ROWS:
        speed_operands[rowname+'/consumer/raw']=load('consumer/'+rowname+'-consumer-input.json')['raw']
    speed_inventory=speed_symbol_inventory(speed_operands)
    speed_absent=all(not item['speedSymbols'] for item in speed_inventory)
    J.emit('prebinding-native-speed-inventory',{'operands':speed_operands,'inventory':speed_inventory,
        'noNativeSpeedSymbols':speed_absent,'nameRule':'case-insensitive c_s*, cs, cs_*, *speed*',
        'scope':'Raw source, chemical, velocity and pressure-consumer operands; before binding',
        'legacyPhysicalCs':physical['parameters']['c_s0'],'legacyCsRetuned':False})
    require(speed_absent,'unaccounted native source/consumer speed dependency')
    profile_length=context['numeric']['L_W']
    scale_rule=native_profile_scale_rule(c2)
    J.emit('native-profile-scale-join',{'savedLength':profile_length,
        'physicalLength':physical['parameters']['L_W'],'nativeRule':scale_rule,'declaredLength':10})
    require(profile_length==10 and profile_length==sp.Rational(physical['parameters']['L_W']),
            'native saved physical profile length')
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
    J.emit('actual-binding-context',{'restored':context,'fixedFrequency':3,'densityStatus':'RESTORED_PENDING_SOURCE_JOINS',
        'regulatedCompositionClaim':False,'noSigmaEtaIdentification':True,'csOnlyInDepth':speed_absent,
        'speedScopeEvidence':'prebinding-native-speed-inventory'})
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
    consumers={};consumer_full={};native_rows={};census={};unit_by_key={}
    units=load('consumer-finish/consumer-unit-joins.json')
    for u in units:unit_by_key[u['row'],str(u['slot'])]=u
    for rowname in ROWS:
        entries=broad_row_census(row_texts[rowname]);unknown=[h for e in entries for h in e['hits'] if h['constructor']!='Symbol' or h['name'] not in SLOTS]
        J.emit(rowname+'-full-native-partition',{'source':provenance,'fullConstructor':row_texts[rowname],
            'children':entries,'unknownPressureAtoms':unknown,'localChildIndices':[e['childIndex'] for e in entries if not e['hits']],
            'pressureChildIndices':[e['childIndex'] for e in entries if e['hits']],
            'localPartUnchanged':True,'localPartRestoredAsScience':False})
        require(not unknown and all(e['completeNameCoverage'] for e in entries),'broad pressure/normal derivative census')
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
        for slotname in SLOTS:
            atom=symbols.get(slotname,sp.Symbol(slotname))
            ablated=row_bound.subs(atom,0);removed=row_bound-ablated
            J.emit(rowname+'-'+slotname+'-native-slot-ablation',{'baseline':row_bound,'ablated':ablated,'removed':removed,'slot':atom})
            J.zero(rowname+'-'+slotname+'-native-ablation-join',removed,coeff[slotname]*atom)
        J.zero(rowname+'-affine-pressure-reconstruction',row_bound,sum(coeff[n]*symbols.get(n,sp.Symbol(n)) for n in SLOTS))
        consumers[rowname]={};consumer_full[rowname]=coeff;native_rows[rowname]=row_bound;census[rowname]=entries
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
    unit_scope=load('consumer-finish/units-and-scope.json')
    source_dim=unit_scope['normalizedSourceDimension']
    require(list(source_dim)==[1,-1,0], 'inherited normalized source dimension')
    native_dim=native_jet_dimensions(next(a['text'] for a in routes if a['name']=='wave_jet'))
    J.emit('native-unrestricted-jet-dimension-rule',native_dim)
    historical=load('consumer-finish/native-pressure-consumer-omission.json')
    J.emit('historical-pressure-ablation-scope',{'records':historical,
        'lowerRecordInPriorEvidence':False,'newLowerAblations':'own native rows above',
        'excludedGrade':[2,1],'includedAsRetainedMixed':False})
    require(len(historical)==2 and {h['row'] for h in historical}=={'THETA_BALANCE','E_W_BALANCE'}, 'historical upper records')
    for h in historical:
        require(str(h['removedSlot'])=='delta_p_plus' and h['mixedPerSource']==0 and h['ablatedRowPerD']==0,'old mixed ablation zero')
        hg=polynomial_terms(h['selectedIncrement'],eta,sigma)
        require(set(hg)=={(2,1)},'old unprojected21 must be present')
        coefficient=hg[2,1]
        formal_factors=[one([coefficient],n) for n in ('epsilon_shape','w1_profile','inherited_whole_bare_mixed_kernel')]
        scalar=sp.cancel(coefficient/sp.prod(formal_factors))
        J.emit('historical-'+h['row']+'-excluded21-formal-nonzero',{'coefficient':coefficient,
            'formalFactors':formal_factors,'scalar':scalar,'nonzeroPhysicalFieldAsserted':False})
        require(not scalar.free_symbols,'old21 coefficient after independent formal factors is constant')
        J.zero('historical-'+h['row']+'-excluded21-reconstruction',coefficient,scalar*sp.prod(formal_factors))
        nonzero(J,'historical-'+h['row']+'-excluded21-scalar',scalar)
        J.zero('historical-'+h['row']+'-native-join',bind(h['native']),native_rows[h['row']])
    sources={};source_jets={};source_full={};x=sp.Symbol('composition_x',real=True)
    w=(1+sp.tanh(x/profile_length))/2;m=(1-sp.tanh(x/profile_length)**2)/3
    def profile(expr):
        mapping={}
        for s in expr.free_symbols:
            hit=re.fullmatch(r'([wm])1_profile((?:_?d[123])*)',s.name)
            if hit:
                base,suffix=hit.groups();directions=re.findall('d([123])',suffix)
                mapping[s]=sp.S.Zero if any(d!='1' for d in directions) else profile_length**len(directions)*sp.diff(w if base=='w' else m,x,len(directions))
        result=sp.cancel(expr.xreplace(mapping))
        J.emit('profile-map-'+str(len(J.artifacts)),{'original':expr,'map':[[a,b] for a,b in mapping.items()],
            'oneDimensional':result,'L':profile_length,'independentSigma':True,'profileTransformEvaluated':False})
        require(result.free_symbols<={x},'coefficient field fully bound')
        return result
    for face in ('plus','minus'):
        inp=load('consumer/'+face+'-source-input.json');saved=load('consumer/'+face+'-source-grade-split.json')
        chemical=load('consumer/native-chemical-amplitude.json')
        chemical_domain=load('consumer/chemical-amplitude-domain.json')
        normalization=load('consumer/inherited-source-normalization.json')
        # The complete native records are already saved; restore scalar leaves only.
        mu_text=text_helpers['tuple_arguments'](native['chemicalSource']['valueConstructorText'])[1]
        density_text=text_helpers['tuple_arguments'](native['geometry']['background_density_map']['cases'][0]['valueConstructorText'])[1]
        native_mu=sp.sympify(mu_text);native_density=sp.sympify(density_text)
        J.emit(face+'-chemical-density-source-operands',{'nativeChemicalConstructor':mu_text,
            'nativeChemical':native_mu,'savedChemical':chemical,'chemicalDomain':chemical_domain,
            'nativeDensityConstructor':density_text,'nativeDensity':native_density,
            'savedDensity':context['density'],'densityMap':context['densityMap']})
        J.zero(face+'-native-chemical-raw-join',native_mu,chemical['raw'][1])
        J.zero(face+'-native-density-raw-join',native_density,context['density'][1])
        J.zero(face+'-native-density-map-join',context['densityMap']['rho_br_bg_rho4_constant'],native_density)
        J.zero(face+'-chemical-domain-original-join',chemical_domain['original'],chemical['amplitude'])
        J.zero(face+'-chemical-domain-reduced-join',chemical_domain['reduced'],chemical['amplitude'])
        J.zero(face+'-chemical-domain-fraction-join',chemical_domain['numerator'],chemical_domain['denominator']*chemical_domain['reduced'])
        J.zero(face+'-chemical-domain-zero-grade-join',chemical_domain['denominator'].subs({eta:0,sigma:0},simultaneous=True),chemical_domain['denominatorAtZero'])
        nonzero(J,face+'-chemical-domain-denominator',chemical_domain['denominatorAtZero'])
        join_source_input(J,face,inp,chemical,normalization,bind,one,eps)
        J.emit(face+'-native-source-join-status',{'densityRestoredAndJoined':True,
            'chemicalAmplitudeJoined':True,'inheritedNormalizationJoined':True,
            'scope':'Exact restored-operand joins, not independent revalidation of the source producer'})
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
                spec=jet_spec(atom.name);orders=spec['spatialOrders'];base_dim=native_dim['velocity' if spec['channel'].startswith('u_') else 'scalar'];dim=[base_dim[0]-sum(orders),base_dim[1]-spec['timeOrder'],base_dim[2]]
                jetrows.append({'atom':atom,'spec':spec,'coefficient':coefficient,'field':field,'jetDimension':dim,
                    'requiredCoefficientDimension':[source_dim[i]-dim[i] for i in range(3)],
                    'coefficientDimensionIndependentlyVerified':False,'nativeDimensionRuleSha256':native_dim['sourceSha256'],
                    'status':'EXACT_ZERO' if field==0 else 'FORMAL_COEFFICIENT_AVAILABLE'})
            J.emit(face+'-source-jets-'+str(g[0])+str(g[1]),{'grade':g,'source':expr,'jets':jetrows,'reconstruction':reconstruction,
                'dimensionScope':'Native unrestricted jet rule checked; normalized source unit inherited. Required coefficient units are inferred expectations, not a new dimensional verification after numeric binding'})
            J.zero(face+'-source-jet-reconstruction-'+str(g[0])+str(g[1]),expr,reconstruction)
            source_jets[face][g]=jetrows
    J.zero('both-face-source-equality',source_full['plus'],source_full['minus'])
    selected_fourier=load('consumer/selected-fourier-contraction.json')
    J.emit('inherited-selected-fourier-context',{'saved':selected_fourier,
        'sourceZeroGrades':{f:sources[f][0,0] for f in ('plus','minus')},
        'scope':'Historical upper off-shell source context only; plane-wave contraction is not replayed or promoted to a general Fourier map'})
    J.zero('selected-fourier-frequency-context',selected_fourier['omega'],context['frequency'])
    for axis,key in ((1,'s11cdTangentialMomentum1'),(2,'s11cdTangentialMomentum2')):
        for label in ('kin','kout'):
            J.zero('selected-fourier-'+label+'-edge-'+str(axis),selected_fourier[label][axis],sp.Rational(physical['parameters'][key]))
    require(selected_fourier['incidentTransverseMode'] is False and
            selected_fourier['lowerFaceCorrectionConstructed'] is False,'historical selected Fourier scope')
    # Response pieces are copied, not reconstructed by the old response producers.
    rc=load('reference/retained-response-census.json');objects=[rc[n] for n in ('flat','heightCoefficient','slopeCoefficient','mixedIteration','taggedTotalMixed')]
    omega=one(objects,'reference_unrestricted_frequency');qi=one(objects,'reference_qi');qo=one(objects,'reference_qo')
    k=one(objects,'reference_k');l=one(objects,'reference_l');p,r=sp.symbols('composition_p composition_r',real=True)
    K,L=sp.symbols('composition_k composition_l',real=True)
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
        fi=load('direct/'+face+'-reference-factor-input.json');fo=load('direct/'+face+'-reference-factor-return.json')
        fmap={s:{'grazing_qi':qi,'grazing_qo':qo}[s.name] for s in fi['right'].free_symbols}
        J.emit(face+'-inherited-isolated-factor',{'input':fi,'return':fo,'map':[[a,b] for a,b in fmap.items()],
            'completeDensityIdentity':False,'recomputedPriorProof':False})
        require(fo['cancelled']==0,'saved isolated factor literal zero')
        J.zero(face+'-isolated-factor-current-right',fi['right'].xreplace(fmap),Rprod)
        oldfactor_map={a:{'q_i':dq_i,'q_o':dq_o}[a.name] for a in before['reference'].free_symbols}
        J.zero(face+'-isolated-factor-current-left',before['reference'].xreplace(oldfactor_map),fi['left'])
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
        Tdelta=T*delta_matrix
        J.emit(face+'-whole-native-trace-direct-injection',{'T':T,'delta':delta_matrix,'product':Tdelta,'T01':T[0,1],
            'newReference12':delta_matrix[1,2],'T01TimesNew12':T[0,1]*delta_matrix[1,2]})
        for i in range(3):
            J.zero(face+'-native-trace-unit-diagonal-'+str(i),T[i,i],sp.S.One)
            for j in range(i):J.zero(face+'-native-trace-lower-'+str(i)+str(j),T[i,j],sp.S.Zero)
            for j in range(3):J.zero(face+'-unchanged-native-Tdelta-'+str(i)+str(j),Tdelta[i,j],delta_matrix[i,j])
        H=slot['affineHeight'];solution=slot['equationSolution'];target=one([solution],face+'_physical_target');jet=one([solution],face+'_jet_slot')
        normal=sp.Symbol('composition_NORMAL',real=True);calls=[]
        def adapter(inputs,diagonal,off_diagonal,source,kout,kin,second=sp.S.Zero):
            require(diagonal==0 and second==0 and source==1,'direct-only formal adapter')
            require(kout==(l,) and kin==(k,), 'adapter ordered arguments')
            calls.append({'diagonal':diagonal,'offDiagonal':off_diagonal,'source':source,'kout':kout,'kin':kin,'second':second})
            return off_diagonal
        ns={'sp':sp,'face':sign,'qo':qo,'NORMAL':normal,'reference':sign*context['numeric']['W_0']/2,
            'reference_matrix':sp.Matrix([[0,P],[0,0]]),'jet_diagonal':sp.S.Zero,'jet_second':sp.S.Zero,
            'inputs':SimpleNamespace(values=context['numeric']),'composed_source':sp.S.One,'ko':(l,),'ki':(k,),'kernel_apply':adapter,
            'trace_map':{'REFERENCE_VALUE_SOLVE':solution,'PHYSICAL_PRESSURE_TARGET':target},'pressure':P,'jet_slot':jet}
        fragments={name:helpers['assignment_source'](c2,'build_face',name) for name in ('reference','extension','jet_transfer','normal_jet','reference_pressure')}
        J.emit(face+'-direct-adapter-input',{'D':D,'physicalInjection':P,'referenceTransferInjection':P,
            'bothDifferOnlyBeyondRetainedRectangle':True,'traceHeight':slot['savedHeight'],'traceSolution':solution,
            'nativeAssignments':fragments,'newAddress':'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL','nativeFunctionsCalled':False})
        for name in fragments:exec(fragments[name],ns)
        J.zero(face+'-native-reference-location-join',ns['reference'],sign*context['numeric']['W_0']/2)
        height=slot['savedHeight'];unprojected=ns['reference_pressure'].subs(H,height)
        pt=polynomial_terms(unprojected,eta,sigma);projected=sum(pt.get(g,sp.S.Zero)*eta**g[0]*sigma**g[1] for g in G)
        J.emit(face+'-direct-adapter-return',{'calls':calls,'jet':ns['normal_jet'],'unprojectedReference':unprojected,
            'excluded':unprojected-projected,'retainedReference':projected,'algebraicRoutingOnly':True,
            'noMiddleIntegral':True,'physicalSourceNotExcited':True})
        J.zero(face+'-direct-normal-routing',ns['normal_jet'],sp.I*sign*qo*P)
        J.zero(face+'-direct-reference-routing',projected,P)
        require(set(pt)<={(1,1),(2,1)}, 'direct trace only excluded21')
        J.zero(face+'-excluded21-trace-correction',unprojected-projected,-height*sp.I*sign*qo*P)
        J.emit(face+'-native-reference-location',{'W0':context['numeric']['W_0'],'reference':ns['reference'],
            'face':sign,'nativeReferenceAssignment':helpers['assignment_source'](c2,'build_face','reference')})
        J.zero(face+'-saved-flat-diagonal-jet',rc['jetKernels'][face]['flat'],sp.I*sign*qo*rc['flat'].subs(qi,qo))
        flat_good=rc['jetKernels'][face]['flat'].subs({omega:3,qo:sp.Rational(3,2)},simultaneous=True)
        flat_bad=(sp.I*sign*qo*rc['flat']).subs({omega:3,qo:sp.Rational(3,2),qi:2},simultaneous=True)
        J.emit(face+'-unbound-flat-pole-control',{'baseline':flat_good,'corrupt':flat_bad,'notPhysicalOnDiagonal':True,'independentDepths':[2,sp.Rational(3,2)]})
        nonzero(J,face+'-unbound-flat-pole-movement',flat_bad-flat_good)
        for key in ('height','slope'):
            J.zero(face+'-saved-'+key+'-jet',rc['jetKernels'][face][key],sp.I*sign*qo*rc[key+'Coefficient'])
        J.zero(face+'-saved-mixed-sum',rc['jetKernels'][face]['mixed'],sp.I*sign*qo*(rc['mixedIteration']+wholeD))
    J.zero('saved-total-mixed-sum',rc['taggedTotalMixed'],rc['mixedIteration']+wholeD)
    # Join inherited raw kernels and completed row-density returns, without redoing their proofs.
    density_symbols={a.name:a for a in direct['density'].free_symbols}
    raw_density_targets={'q_i':dq_i,'q_o':dq_o,'q_h':density_symbols['grazing_qh'],
        'q_s':density_symbols['grazing_qs'],'increment_transfer':density_symbols['grazing_transfer'],
        'increment_difference':density_symbols['grazing_output']-density_symbols['k'],
        'k':density_symbols['k']}
    for rowname in ROWS:
        rawrow=load('raw/'+rowname+'-retained-increment.json')
        for face in ('plus','minus'):
            kernel=rawrow['rawKernelPlus' if face=='plus' else 'rawKernelMinus']
            before=load('raw/'+face+'-closed-before-cancel.json')
            J.zero(rowname+'-'+face+'-raw-kernel-operand-join',kernel,before['raw'])
            J.zero(rowname+'-'+face+'-raw-source-zero-join',rawrow['sourcePlus' if face=='plus' else 'sourceMinus'],sources[face][0,0])
        mixed=rawrow['mixedCoefficient']
        mixmap={a:{'q_i':qi,'q_o':qo}.get(a.name,a) for a in mixed.free_symbols}
        bare_terms={face:sp.Symbol('composition_bare_'+face) for face in ('plus','minus')}
        for a in mixed.free_symbols:
            if a.name in ('increment_raw_plus','increment_raw_minus'):mixmap[a]=bare_terms[a.name.rsplit('_',1)[1]]
        joined=sum(consumers[rowname]['delta_p_'+f][0,0]*sources[f][0,0]*Rprod*bare_terms[f] for f in ('plus','minus'))
        J.emit(rowname+'-raw-direct-row-map',{'original':mixed,'map':[[a,b] for a,b in mixmap.items()],
            'currentBareDirectRow':joined,'closedWholeUsedHere':False})
        J.zero(rowname+'-raw-direct-row-join',mixed.xreplace(mixmap),joined)
        if rowname.startswith('U'):continue
        old_in=load('direct/'+rowname+'-actual-closed-row-density-original-input.json')
        old_out=load('direct/'+rowname+'-actual-closed-row-density-canonical-return.json')
        rawmap={a:raw_density_targets[a.name] for a in rawrow['rawRowDensity'].free_symbols if a.name in raw_density_targets}
        current=direct['density'].subs(O,3)*sum(consumers[rowname]['delta_p_'+f][0,0]*sources[f][0,0] for f in ('plus','minus'))
        J.emit(rowname+'-inherited-full-row-density-proof',{'originalInput':old_in,'return':old_out,
            'rawMap':[[a,b] for a,b in rawmap.items()],'currentDirect11Row':current,
            'recomputedPriorCancellation':False,'perFaceDirectWholeMultiplicity':1})
        require(old_out['cancelled']==0,'inherited full row density literal zero')
        J.sinh_zero(rowname+'-full-row-density-left-join',old_in['left'],rawrow['rawRowDensity'].subs(rawmap,simultaneous=True))
        J.sinh_zero(rowname+'-full-row-density-right-join',old_in['right'],current)
    left_join=load('reference/height-left-plus-trace-input.json')
    left_return=load('reference/height-left-plus-trace-return.json')
    right_join=load('reference/right-height-PV-operands.json')
    require(left_return['cancelled']==0,'inherited left-assignment plus trace literal zero')
    native_iteration=[]
    for face in ('plus','minus'):
        before=load('reference/'+face+'-reference-before-guards.json')
        parts=before['uncombined'];require(len(parts)==3 and before['noDirectInNativeSecond'] is True,'two native assignments plus trace only')
        middle=one([parts[1]],'reference_m');transfer=one([right_join['originalCoefficient']],'reference_t')
        J.emit(face+'-mixed-native-assignment-coverage',{'savedUncombined':parts,
            'roles':['left-height/right-slope','left-slope/right-height','native trace subtraction'],
            'leftInput':left_join,'leftReturn':left_return,'rightInput':right_join,
            'rightArgumentMap':[[middle,k+transfer]],'oldFunctionsCalled':False,'priorCancellationReplayed':False})
        J.zero(face+'-left-trace-input-join',parts[0]+parts[2],left_join['left'])
        J.zero(face+'-right-assignment-input-join',parts[1].subs(middle,k+transfer),right_join['originalCoefficient'])
        native_iteration.append({'face':face,'leftCoefficient':left_join['right'],'rightDensity':right_join['density']})
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
    expected_conditions=(sp.Lt(depth_variable**2-9/cs**2,-sp.Rational(1,20)),
                         sp.Gt(depth_variable**2-9/cs**2,-sp.Rational(1,20)),sp.true)
    J.emit('physical-depth-branch-conditions',{'actual':[a.cond for a in qdefinition.args],'expected':expected_conditions})
    require(tuple(a.cond for a in qdefinition.args)==expected_conditions,'physical branch conditions')
    J.zero('mapped-positive-depth-branch',qdefinition.args[0].expr,sp.sqrt(rad))
    J.zero('mapped-decaying-depth-branch',qdefinition.args[1].expr,sp.I*sp.sqrt(-rad))
    J.emit('inherited-depth-map',{'saved':old_q,'map':[[a,b] for a,b in qm_map.items()],'new':qdefinition})
    new_direct_map={'grazing_unrestricted_frequency':sp.Integer(3),'grazing_output':L,'k':K,
        'grazing_transfer':td,'grazing_qi':qfun(K),'grazing_qo':qfun(L),
        'grazing_qh':qfun(K+td),'grazing_qs':qfun(L-td)}
    require({a.name for a in direct['density'].free_symbols}==set(new_direct_map),'complete direct map symbol census')
    actual_direct_map={a:new_direct_map[a.name] for a in direct['density'].free_symbols}
    mapped_direct=direct['density'].subs(actual_direct_map,simultaneous=True)
    # Join the inherited full density and argument map as complete objects.
    old_join=load('reference/restored-direct-argument-join.json')
    J.zero('saved-direct-density-provenance',old_join['savedDensity']['density'],direct['density'])
    old_mapped=direct['density'].subs(dict(old_join['map']),simultaneous=True)
    J.zero('saved-direct-argument-provenance',old_mapped,old_join['transported'])
    jd=load('reference/right-height-PV-operands.json')['density']
    new_iter_map={'reference_unrestricted_frequency':sp.Integer(3),'reference_k':K,'reference_l':L,
        'reference_t':t,'reference_qi':qfun(K),'reference_qo':qfun(L),'reference_qm':qfun(K+t)}
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
    hdef=tagdefs['H']['savedDefinition'];hs=one([hdef['ordinarySubtractedIntegrand']],'reference_left_height_transfer')
    J.zero('actual-H-bound-variable-map',hdef['variableChange'][1],l-hs)
    tagdefs['H'].update(freeMomenta=['l-k'],boundVariable=hs,middleMomentum=None,
        originalChangeOfVariable=hdef['variableChange'],boundVariableAssumptions=hs.assumptions0,
        integrationNotSameAsJwhole=True)
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
    for joined in native_iteration:
        J.zero(joined['face']+'-mixed-census-from-inherited-assignments',rc['mixedIteration'],joined['leftCoefficient']*wholeH+wholeJ)
    Htag=sp.Function('Hwhole')(L-K,sp.Integer(1),sp.Integer(10))
    base_components={ (0,0):[('NATIVE_FLAT',rc['flat'].subs({omega:3,qi:qo},simultaneous=True))],
        (1,0):[('NATIVE_HEIGHT',rc['heightCoefficient'].subs(omega,3))],
        (0,1):[('NATIVE_SLOPE',rc['slopeCoefficient'].subs(omega,3))] }
    def components_for(face,grade):
        if grade!=(1,1):return base_components[grade]
        signature=(L,K,sp.Integer(3),cs,sp.Rational(1,5),sp.Rational(1,10),sp.Integer(1),sp.Integer(10))
        Jtag=sp.Function('Jwhole_'+face)(*signature);Dtag=sp.Function('Dwhole_'+face)(*signature)
        iteration=rc['mixedIteration'].subs({omega:3,wholeH:Htag,wholeJ:Jtag},simultaneous=True)
        return [('NATIVE_MIXED_ITERATION',iteration),('INHERITED_DIRECT_WHOLE_OFF_DIAGONAL',Dtag)]
    component_maps={};component_cache={}
    for face in ('plus','minus'):
        for grade in G:
            for component,rv in components_for(face,grade):
                mapping={omega:sp.Integer(3),qi:qfun(L if grade==(0,0) else K),qo:qfun(L),k:K,l:L}
                actual_map={a:b for a,b in mapping.items() if a in rv.free_symbols}
                mapped=rv.subs(actual_map,simultaneous=True)
                require(not (rv.free_symbols & {omega,qi,qo,k,l}) - set(actual_map) and
                        not any(a.name.startswith('reference_') and a not in actual_map for a in rv.free_symbols),'complete response argument map')
                key=face+'-'+str(grade)+'-'+component
                digest=hashlib.sha256(sp.srepr(sp.Tuple(*[sp.Tuple(a,b) for a,b in actual_map.items()])).encode()).hexdigest()
                component_maps[face,grade,component]={'id':key,'sha256':digest,'original':rv,'mapped':mapped,
                    'map':[[a,b] for a,b in actual_map.items()], 'symbolAssumptions':[[a,a.assumptions0] for a in actual_map],
                    'flatSupport':'k=l' if grade==(0,0) else None,'frequency':3,'positiveRegulatorContinuation':False}
                component_cache[face,grade,component]=mapped
                J.emit('component-map-'+str(len(component_maps)),component_maps[face,grade,component])
        assembled=sum(v for _,v in components_for(face,(1,1)))
        signature=(L,K,sp.Integer(3),cs,sp.Rational(1,5),sp.Rational(1,10),sp.Integer(1),sp.Integer(10))
        back={Htag:wholeH,sp.Function('Jwhole_'+face)(*signature):wholeJ,sp.Function('Dwhole_'+face)(*signature):wholeD}
        transported=assembled.xreplace(back)
        J.emit(face+'-assembled-mixed-components',{'components':components_for(face,(1,1)),'sum':assembled,
            'backMap':[[a,b] for a,b in back.items()],'transported':transported})
        J.zero(face+'-assembled-mixed-reference',transported,rc['taggedTotalMixed'].subs(omega,3))
        for grade in G:
            for component,rv in components_for(face,grade):
                tags=rv.atoms(sp.Function)
                require(all(v.func.__name__ in ('reference_height_hat','reference_slope_hat','Hwhole','Jwhole_'+face,'Dwhole_'+face) for v in tags),'zero-profile tag census')
                zero_map={v:sp.S.Zero for v in tags}
                reduced=rv.xreplace(zero_map)
                J.emit(face+'-'+component+'-zero-profile-response',{'original':rv,'formalProfileMap':[[a,b] for a,b in zero_map.items()],
                    'reduced':reduced,'distributionValueEvaluated':False})
                J.zero(face+'-'+component+'-zero-profile-response-identity',reduced,rv if grade==(0,0) else sp.S.Zero)
        J.zero(face+'-assembled-mixed-normal',sp.I*(1 if face=='plus' else -1)*qo*transported,
               rc['jetKernels'][face]['mixed'].subs(omega,3))
    coverage=[];addresses=[];consumer_fields={};response_tokens={}
    def response_token(face,slot,b,component,c,atom):
        key=(face,slot,b,component,c,str(atom))
        if key not in response_tokens:response_tokens[key]=sp.Symbol('ordered_response_'+str(len(response_tokens)))
        return response_tokens[key]

    factor_proofs={}
    def join_address_factor(addr):
        face=addr['face'];flat=addr['responseGrade']==(0,0)
        normal=rc['normalJet'][face] if addr['slot']=='normal' else sp.S.One
        original=addr['responseOriginal']*normal
        roles={omega:sp.Integer(3),qi:qfun(L if flat else K),qo:qfun(L),k:K,l:L}
        mapping={a:b for a,b in roles.items() if a in original.free_symbols}
        response_map={a:b for a,b in roles.items() if a in addr['responseOriginal'].free_symbols}
        actual=addr['responseCoefficient']*addr['normalMultiplier']
        operands=(addr['responseOriginal'],addr['normalOriginal'],normal,tuple(mapping.items()),
                  tuple(addr['responseMap']['map'][i] for i in range(len(addr['responseMap']['map']))),actual)
        digest=hashlib.sha256(repr(tuple(sp.srepr(v) if isinstance(v,sp.Basic) else repr(v) for v in operands)).encode()).hexdigest()
        if digest not in factor_proofs:
            name='address-full-factor-'+str(len(factor_proofs))
            J.emit(name+'-operands',{'addressId':addr['addressId'],'original':original,'savedNormal':normal,
                'addressNormalOriginal':addr['normalOriginal'],'requiredMap':[[a,b] for a,b in mapping.items()],
                'actualResponseMap':addr['responseMap'],'mappedAddressFactor':actual,'flatSupport':flat})
            verified_map=complete_factor_map(addr['responseOriginal'].free_symbols,normal.free_symbols,roles,addr['responseMap']['map'])
            require(verified_map==mapping,'complete normal and response role map')
            J.zero(name+'-normal-source-join',addr['normalOriginal'],normal)
            J.zero(name+'-full-mapped-residual',original.subs(mapping,simultaneous=True),actual)
            factor_proofs[digest]=(operands,name)
        else:require(factor_proofs[digest][0]==operands,'identical factor operands before local proof reuse')
        return {'proof':factor_proofs[digest][1],'operandSha256':digest,'completeNormalMap':[[a,b] for a,b in mapping.items()]}

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
                                reason=zero or ('EXACT_ZERO_SOURCE_JET' if item['field']==0 else 'EXACT_ZERO_RESPONSE' if rv==0 else None)
                                reduced=b==(0,0)
                                addr={**core,'addressId':len(addresses),'component':component,'jet':item['spec'],'sourceField':item['field'],
                                    'responsePlaceholder':response_token(face,slot,b,component,c,item['atom']),
                                    'consumerField':cf,'epsilonCount':epsilon_power(bool(reason)),'epsilon':eps,'consumerOriginal':cv,
                                    'epsilonMeaning':'exact zero has no power; otherwise explicit homogeneous factor, not proof of nonzero field',
                                    'sourceOriginal':item['coefficient'],'sourceAtom':item['atom'],
                                    'sourceTransform':transform(item['field'],'l-p' if reduced else 'k-p'),
                                    'consumerTransform':transform(cf,'r-l'),'waveMultiplier':wave(item['spec'],p),
                                    'responseCoefficient':component_cache[face,b,component], 'responseOriginal':rv,
                                    'responseMap':component_maps[face,b,component],
                                    'normalMultiplier':sp.I*sign*qfun(L) if slotkind=='normal' else sp.S.One,
                                    'normalOriginal':sp.I*sign*qo if slotkind=='normal' else sp.S.One,
                                    'responseInputDepth':'q(l)' if reduced else 'q(k)','responseOutputDepth':'q(l)',
                                    'freeMomenta':['p','l','r'] if reduced else ['p','k','l','r'],
                                    'measure':'dl' if reduced else 'dl dk','deltaSupport':['k=l'] if reduced else [],
                                    'status':reason or 'FORMAL_ADDRESS_AVAILABLE_NONZERO_NOT_ASSERTED','globalComposition':'UNRESOLVED',
                                    'nativeRowChildHashes':[e['sha256'] for e in census[rowname] if any(h['name']==slot for h in e['hits'])],
                                    'sourceOperandSha256':copies['consumer/'+face+'-source-grade-split.json']['sha256'],
                                    'responseCensusSha256':copies['reference/retained-response-census.json']['sha256'],
                                    'wholeValueEvaluated':False}
                                if x not in item['field'].free_symbols:addr['deltaSupport'].append('l=p' if reduced else 'k=p')
                                if x not in cf.free_symbols:addr['deltaSupport'].append('r=l')
                                addr['fullFactorProof']=join_address_factor(addr)
                                addresses.append(addr)
        J.emit(rowname+'-ordered-addresses',[a for a in addresses if a['row']==rowname])
    require(len(coverage)==len(ROWS)*2*2*16,'complete 16-triple coverage per face/slot/row')
    J.emit('grade-coverage',coverage);J.emit('whole-coefficient-transforms',ft)
    # Actual native affine rows, with independent source/response addresses substituted.
    # Tokens retain the response/source-grade/jet identity; they have no evaluated field value.
    for rowname in ROWS:
        for face in ('plus','minus'):
            substitutions={a:sp.S.Zero for a in native_rows[rowname].free_symbols if a.name in SLOTS}
            slot_inputs={}
            for prefix in ('delta_p_','d_w_delta_p_'):
                slot=prefix+face;signal=sp.S.Zero
                for b,c in itertools.product(G,repeat=2):
                    if any(b[i]+c[i]>1 for i in (0,1)):continue
                    for component,_ in components_for(face,b):
                        for item in source_jets[face][c]:
                            token=response_token(face,slot,b,component,c,item['atom'])
                            signal+=eta**(b[0]+c[0])*sigma**(b[1]+c[1])*token*item['coefficient']*item['atom']
                slot_inputs[slot]=signal
                for atom in substitutions:
                    if atom.name==slot:substitutions[atom]=signal
            substituted=native_rows[rowname].subs(substitutions,simultaneous=True)
            numerator,denominator=sp.fraction(sp.cancel(substituted))
            n=polynomial_terms(numerator,eta,sigma);d=polynomial_terms(denominator,eta,sigma)
            nonzero(J,rowname+'-'+face+'-formal-row-denominator',d.get((0,0),sp.S.Zero))
            actual=quotient_recurrence(n,d,G,sp.cancel)
            J.emit(rowname+'-'+face+'-actual-row-placeholder-input',{'nativeRow':native_rows[rowname],
                'substitution':[[a,b] for a,b in substitutions.items()],'slotSignals':slot_inputs,
                'substituted':substituted,'tokenMeaning':'Response component acting on the specified source grade and jet; no kernel evaluation',
                'operatorOrder':['consumer multiplication','response','source multiplication then wave derivative']})
            for g in G:
                selected,address_sum=address_sum_for(addresses,rowname,face,g)
                J.emit(rowname+'-'+face+'-actual-address-sum-'+str(g),{'addressIds':[a['addressId'] for a in selected],
                    'sum':address_sum,'directRowCoefficient':actual[g]})
                J.zero(rowname+'-'+face+'-actual-address-reconstruction-'+str(g),actual[g],address_sum)
    J.emit('response-placeholder-index',[{'key':key,'token':v,
        'actualAddresses':[{'id':a['addressId'],'mappedResponse':a['responseCoefficient']*a['normalMultiplier'],
                           'operand':a['responseOriginal']*a['normalOriginal'],'map':a['responseMap']['sha256']}
                          for a in addresses if a['responsePlaceholder']==v]}
        for key,v in response_tokens.items()])
    # Zero-profile and constant-coefficient reductions do not evaluate response convolutions.
    def profile_limit(expr,constant):
        mapping={}
        for a in expr.free_symbols:
            match=re.fullmatch(r'([wm])1_profile((?:_?d[123])*)',a.name)
            if match:mapping[a]=constant[match.group(1)] if not match.group(2) else sp.S.Zero
        return expr.xreplace(mapping),mapping
    const_profile={'w':sp.Symbol('constant_height',real=True),'m':sp.Symbol('constant_modulation',real=True)}
    for face in ('plus','minus'):
        for g in G:
            reduced,mapping=profile_limit(sources[face][g],{'w':sp.S.Zero,'m':sp.S.Zero})
            constant,cmap=profile_limit(sources[face][g],const_profile)
            J.emit(face+'-source-profile-limits-'+str(g),{'source':sources[face][g],
                'zeroMap':[[a,b] for a,b in mapping.items()],'zeroProfile':reduced,
                'constantMap':[[a,b] for a,b in cmap.items()],'constantCoefficientSource':constant})
            J.zero(face+'-zero-profile-source-'+str(g),reduced,sources[face][0,0] if g==(0,0) else sp.S.Zero)
    for rowname in ROWS:
        for slot in SLOTS:
            for g in G:
                reduced,mapping=profile_limit(consumers[rowname][slot][g],{'w':sp.S.Zero,'m':sp.S.Zero})
                J.emit(rowname+'-'+slot+'-zero-profile-'+str(g),{'map':[[a,b] for a,b in mapping.items()],'value':reduced})
                J.zero(rowname+'-'+slot+'-zero-profile-consumer-'+str(g),reduced,consumers[rowname][slot][0,0] if g==(0,0) else sp.S.Zero)
    J.emit('constant-end-and-zero-profile-scope',{'sourceConsumerConstantMultipliers':'constant*delta of own transfer',
        'flatResponseSupport':'k=l; complete input pole and normal prefactor use q(l)',
        'zeroProfileNonflatTags':'set height/slope/H/Jwhole/Dwhole to zero as formal profile-degree reduction',
        'nonzeroConstantHeightResponseReduction':'NOT_COMPUTED; no new distribution product or constant-end response proof',
        'fullOperatorConstantEndClaim':False})
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
    for momentum,depth in qpoint.items():
        radpoint=rad.subs({depth_variable:momentum,cs:sp.sqrt(sp.Rational(10,7))},simultaneous=True)
        selected_depth=qdefinition.subs({depth_variable:momentum,cs:sp.sqrt(sp.Rational(10,7))},simultaneous=True)
        J.emit('physical-control-depth-input-'+str(len(J.artifacts)),{'momentum':momentum,'depth':depth,'radicand':radpoint,
            'selectedSheetDepth':selected_depth,'positiveDepth':depth.is_positive,'positiveRadicand':radpoint.is_positive})
        require(depth.is_positive is True and radpoint.is_positive is True,'control on positive outgoing branch')
        J.zero('physical-control-depth-'+str(len(J.artifacts)),depth,selected_depth)
    def sensitivity(name,baseline,corrupt,context_record):
        movement=sp.cancel(corrupt-baseline)
        J.emit(name+'-control-operands',{'baselineCoefficient':baseline,'corruptCoefficient':corrupt,'movement':movement,
             'context':context_record,'formalTagCoefficientOnly':True,'physicalConvolutionNotEvaluated':True,
             'actualTransformAtTransfer':'NOT_EVALUATED_OR_CERTIFIED_NONZERO'})
        cert=nonzero(J,name+'-movement',movement);controls.append({'name':name,'movement':movement,'certificate':cert,'formalTagCoefficientOnly':True})
    def nonconstant_profile(name,field):
        # Exact polynomial in the already-bound tanh profile; no Fourier value is inferred.
        tau=sp.Symbol('composition_tanh_variable',real=True)
        rational=sp.cancel(field.xreplace({sp.tanh(x/profile_length):tau}))
        numerator,denominator=sp.fraction(rational)
        initial={'field':field,'profileVariable':tau,'map':[[sp.tanh(x/profile_length),tau]],
            'rational':rational,'numerator':numerator,'denominator':denominator,
            'transformValue':'NOT_EVALUATED_OR_CERTIFIED_NONZERO','pointwiseProfileProbeUsed':False}
        J.emit(name+'-input',initial)
        if x not in field.free_symbols or not numerator.free_symbols<={tau} or denominator.free_symbols:
            result={'status':'NOT_CERTIFIED_NONCONSTANT','reason':'constant or outside polynomial tanh coefficient class'}
            J.emit(name+'-decision',result);return result
        terms=polynomial_terms(numerator,tau,sp.Symbol('unused_profile_grade'))
        coefficients={g[0]:v for g,v in terms.items() if g[1]==0}
        require(len(coefficients)==len(terms),'one independent profile polynomial variable')
        eligible=[(n,v) for n,v in sorted(coefficients.items()) if n>0 and v.is_zero is False and v.is_finite is True]
        J.emit(name+'-polynomial',{'coefficients':{str(n):v for n,v in coefficients.items()},
            'eligibleNonconstantCoefficients':eligible,'constantDenominator':denominator})
        if not eligible or not all(v.is_finite is True and not v.free_symbols for v in coefficients.values()):
            result={'status':'NOT_CERTIFIED_NONCONSTANT','reason':'no exact nonzero positive-degree coefficient or unresolved finiteness'}
            J.emit(name+'-decision',result);return result
        nonzero(J,name+'-denominator',denominator)
        power,coefficient=eligible[-1];nonzero(J,name+'-positive-degree-coefficient',coefficient)
        J.zero(name+'-field-reconstruction',field,(sum(v*tau**n for n,v in coefficients.items())/denominator).subs(tau,sp.tanh(x/profile_length)))
        result={'status':'NONCONSTANT_POLYNOMIAL_IN_TANH_CERTIFIED','positivePower':power,
            'nonzeroCoefficient':coefficient,'constantDenominator':denominator,
            'argument':'tanh(x/10) ranges over an interval; a polynomial with a nonzero positive-degree coefficient is nonconstant',
            'transformValue':'NOT_EVALUATED_OR_CERTIFIED_NONZERO','deltaOnlyConstantField':False}
        J.emit(name+'-decision',result);return result
    def choose_source(face,grade,spatial):
        candidates=[]
        for index,item in enumerate(source_jets[face][grade]):
            if spatial:
                applicable=item['spec']['spatialOrders'][0]>0 and x in item['field'].free_symbols
                certificate=nonconstant_profile(face+'-source10-candidate-'+str(index),item['field']) if applicable else None
                certified=certificate is not None and certificate['status']=='NONCONSTANT_POLYNOMIAL_IN_TANH_CERTIFIED'
                val=wave(item['spec'],sp.Rational(3,2))
                candidates.append({'jet':item,'formalWaveCoefficient':val,'eligible':applicable and certified,
                    'zero':val.is_zero,'finite':val.is_finite,'nonconstantCertificate':certificate,
                    'transformValue':'NOT_EVALUATED_OR_CERTIFIED_NONZERO'})
            else:
                applicable=x not in item['field'].free_symbols
                val=sp.cancel(item['field']*wave(item['spec'],sp.Rational(3,2))) if applicable else None
                candidates.append({'jet':item,'constantFieldTimesWave':val,'eligible':applicable,
                    'zero':None if val is None else val.is_zero,'finite':None if val is None else val.is_finite})
        selected=select_certified_candidate(candidates)
        chosen=None if selected is None else selected['jet']
        J.emit(face+'-source-control-candidates-'+str(grade),{'candidates':candidates,
            'selected':selected,'noneApplicable':chosen is None,'unknownTreatedAsApplicable':False,
            'pointwiseProfileNonzeroUsed':False})
        require(chosen is not None,'no certified applicable source control; evidence preserved')
        return {**chosen,'controlSelectionEvidence':selected}
    selected_sources={(f,g):choose_source(f,g,g==(1,0)) for f in ('plus','minus') for g in ((1,0),(0,0))}
    def control_address(row,face,slot,a,b,c,jet,component):
        matches=[v for v in addresses if v['row']==row and v['face']==face and v['slot']==slot
                 and v['consumerGrade']==a and v['responseGrade']==b and v['sourceGrade']==c
                 and v['sourceAtom']==jet['atom'] and v['component']==component]
        J.emit('control-address-selection-'+str(len(J.artifacts)),{'requested':[row,face,slot,a,b,c,jet['atom'],component],
            'matches':matches,'count':len(matches)})
        require(len(matches)==1 and matches[0]['epsilonCount']==1,'one actual nonzero-candidate address')
        return matches[0]
    def response_scalar(address,kin,kout,strip_slope=False):
        original=address['responseOriginal']*address['normalOriginal']
        substitutions={omega:sp.Integer(3),k:kin,l:kout,qi:qpoint[kin],qo:qpoint[kout]}
        if strip_slope:
            slopes=[v for v in original.atoms(sp.Function) if v.func.__name__=='reference_slope_hat']
            require(len(slopes)==1,'one actual saved slope profile factor')
            substitutions[slopes[0]]=sp.S.One
        value=original.subs(substitutions,simultaneous=True)
        J.emit('control-response-map-'+str(len(J.artifacts)),{'addressId':address['addressId'],'original':original,
            'map':[[a,b] for a,b in substitutions.items()],'value':value,
            'slopeTagCoefficientOnly':strip_slope,'responseMapSha256':address['responseMap']['sha256']})
        require(not value.free_symbols and not value.atoms(sp.Function),'fully bound scalar control coefficient')
        return value
    for rowname in ('THETA_BALANCE','E_W_BALANCE'):
        for face,sign in (('plus',1),('minus',-1)):
            cp=consumer_fields[rowname,'delta_p_'+face][0,0];cn=consumer_fields[rowname,'d_w_delta_p_'+face][1,0]
            cncert=nonconstant_profile(rowname+'-'+face+'-consumer10-nonconstant',cn)
            J.emit(rowname+'-'+face+'-consumer-control-applicability',{'pressure':cp,'normal':cn,
                'nonconstantCertificate':cncert,'offDeltaTransfers':['30/13-2','30/13-3/2'],
                'transformValue':'NOT_EVALUATED_OR_CERTIFIED_NONZERO'})
            require(x not in cp.free_symbols,'constant pressure consumer')
            nonzero(J,rowname+'-'+face+'-pressure-applicability',cp)
            require(cncert['status']=='NONCONSTANT_POLYNOMIAL_IN_TANH_CERTIFIED','off-delta consumer needs nonconstant coefficient field')
            src10=selected_sources[face,(1,0)];src00=selected_sources[face,(0,0)]
            aflat=control_address(rowname,face,'pressure',(0,0),(0,0),(1,0),src10,'NATIVE_FLAT')
            flat3=response_scalar(aflat,sp.Integer(2),sp.Integer(2))
            J.zero(rowname+'-'+face+'-flat-control-saved-join',aflat['responseOriginal'],rc['flat'].subs({omega:3,qi:qo},simultaneous=True))
            J.zero(rowname+'-'+face+'-flat-control-scalar',flat3,sp.Rational(3,10)/(sp.Rational(3,2)+beta3))
            wp=wave(src10['spec'],sp.Rational(3,2));wk=wave(src10['spec'],sp.Integer(2))
            context10={'addressId':aflat['addressId'],'sourceJet':src10,'sourceTransform':transform(src10['field'],'2-3/2'),
                'support':'k=l=r=2,p=3/2','csSquared':sp.Rational(10,7),'actualResponseScalar':flat3}
            sensitivity(rowname+'-'+face+'-source-p-to-k',cp*flat3*wp,cp*flat3*wk,context10)
            sensitivity(rowname+'-'+face+'-omit-source10',cp*flat3*wp,sp.S.Zero,context10)
            sourceval=src00['field']*wave(src00['spec'],sp.Rational(3,2));nonzero(J,rowname+'-'+face+'-source00-applicability',sourceval)
            aslope=control_address(rowname,face,'normal',(1,0),(0,1),(0,0),src00,'NATIVE_SLOPE')
            slope_normal=response_scalar(aslope,sp.Rational(3,2),sp.Integer(2),True)
            J.zero(rowname+'-'+face+'-slope-control-saved-join',aslope['responseOriginal']*aslope['normalOriginal'],rc['jetKernels'][face]['slope'].subs(omega,3))
            tags={'addressId':aslope['addressId'],'consumerTransform':transform(cn,'30/13-2'),
                'slopeProfile':'j(2-3/2)','support':'p=k=3/2,l=2,r=30/13','sourceJet':src00,'csSquared':sp.Rational(10,7)}
            baseline=slope_normal*sourceval
            sensitivity(rowname+'-'+face+'-normal-q-l-to-r',baseline,baseline*qpoint[sp.Rational(30,13)]/qpoint[sp.Integer(2)],tags)
            sensitivity(rowname+'-'+face+'-omit-consumer10',baseline,sp.S.Zero,tags)
            adirect=control_address(rowname,face,'pressure',(0,0),(1,1),(0,0),src00,'INHERITED_DIRECT_WHOLE_OFF_DIAGONAL')
            dtag=adirect['responseCoefficient'];require(dtag.func.__name__=='Dwhole_'+face,'actual direct tag')
            pointtag=dtag.subs({K:sp.Rational(3,2),cs:sp.sqrt(sp.Rational(10,7))},simultaneous=True)
            J.zero(rowname+'-'+face+'-direct-whole-unit-multiplier',adirect['normalOriginal'],sp.S.One)
            directcoef=cp*sourceval
            for label,factor in (('omit-direct',0),('double-direct',2)):
                sensitivity(rowname+'-'+face+'-'+label,directcoef,factor*directcoef,
                    {'addressId':adirect['addressId'],'tag':pointtag,'fullSignature':list(pointtag.args),
                     'support':['k=p=3/2','r=l'],'freeExternalMomentum':L,'sourceJet':src00,
                     'grades':[[0,0],[1,1],[0,0]],'externalResolventMultiplier':1})
            if face=='minus':
                anormal=control_address(rowname,face,'normal',(1,0),(0,0),(0,0),src00,'NATIVE_FLAT')
                normflat=response_scalar(anormal,sp.Rational(3,2),sp.Rational(3,2))*sourceval
                J.zero(rowname+'-lower-flat-control-saved-join',anormal['responseOriginal']*anormal['normalOriginal'],rc['jetKernels'][face]['flat'].subs(omega,3))
                sensitivity(rowname+'-lower-jet-sign',normflat,-normflat,
                    {'addressId':anormal['addressId'],'consumerTransform':transform(cn,'30/13-3/2'),
                     'support':'p=k=l=3/2,r=30/13','grade':[1,0]})
    J.emit('responsive-formal-controls',controls)
    J.emit('applicability-obligations',{'globalComposition':'UNRESOLVED','globalTestSpace':'UNRESOLVED','composedGrazingLimit':'UNRESOLVED',
        'sourceRegulatorContinuation':'NOT_CLAIMED','sourceCoefficientDimensionAudit':'NOT_RECONSTRUCTED_AFTER_NUMERICAL_BINDING; required units recorded, source unit inherited','internalCertifiedDomain':{'cs':[1,2],'k/l':[-3,3]},
        'hiddenCutoff':False,'oldCutoffs4and6Reused':False,'currentOrPowerComputed':False,'savedWorkReplayed':False})
    J.emit('consumed-source-index',{'files':sorted(used),'copies':copies})
    return {'executionStatus':'COMPLETED_BOUNDED_SOURCE_COMPOSITION_INVENTORY','bothFaces':True,'gradeTriplesPerRowFaceSlot':16,
        'coverageEntries':len(coverage),'addressEntries':len(addresses),'formalResponsiveControls':len(controls),
        'realFrequency':3,'scientificAcceptance':False,'integralsEvaluated':False,'finiteSolves':0,'productionChanges':False,
        'globalComposition':'UNRESOLVED','lossClaim':False,'independentPhysicsValidationByControls':False}


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest)
    verify_invocation(args,gate,sys.argv)
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
