"""New separate-tail checks and restoration of the completed failed prefix.

Import is inert. All source restoration and new arithmetic run under containment.
No old preamble, exact identity, tail-expression or numerical-bank function runs.
"""
import ast,copy,hashlib,json
from pathlib import Path
from S11c_d_defect_packet_contracted_continue_resume import (
    require,ZERO,UNITS,COMPONENTS,canonical_slots,check_census,group_indices,
    exact_constant,guard_inventory,prefix_return,join_base_census,bound_source,independent_slots,tail_formula_join)
from S11c_d_defect_packet_contracted_continue_tail import continue_geometry


def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canonical(v):return json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False)


def restore_prefix(raw,manifest,J):
    prior=J.out/'prior';root=prior/'complete';get=lambda n:read(root/(n+'.json'))
    failed=get('checks');require(failed['status']=='FAILED_PRESERVED' and failed['activeOperation'] is None,'actual original failed stage')
    require(failed['failure']==(prior/'defect_packet_contracted_numeric.stderr').read_text(),'strict original failure identity')
    require('actual positive per-address per-primitive tail allocation' in failed['failure'],'declared original refusal')
    require((prior/'defect_packet_contracted_numeric.stdout').read_bytes()==(root/'checks.json').read_bytes(),'original stdout checks bytes')
    require(not (root/'numerical.sqlite').exists(),'no original numerical action')
    records=[json.loads(v) for v in (root/'evidence-chain.jsonl').read_text().splitlines()]
    names={};previous='0'*64
    for i,record in enumerate(records):
        body={k:v for k,v in record.items() if k!='chainSha256'};p=root/(body['name']+'.json')
        require(body['sequence']==i and body['previous']==previous and body['name'] not in names,'original append-only sequence')
        require(hashlib.sha256(canonical(body).encode()).hexdigest()==record['chainSha256'] and sha(p)==body['sha256'] and p.stat().st_size==body['bytes'],'original chain operands')
        names[body['name']]=i;previous=record['chainSha256']
    require(names['new-tail-address-budget-8360-J-27']+1==names['failure'],'exact failed checkpoint')
    budget_names=[n for n in names if n.startswith('new-tail-address-budget-')]
    budget_records=[{'name':n,'chainSequence':names[n],'record':get(n),'originalControlFlowDisposition':'refused' if i==32 else 'passed before the next checkpoint; inherited control-flow attestation only'} for i,n in enumerate(budget_names)]
    J.emit('restored-all-original-budget-checkpoints',{'records':budget_records,'independentInequalityReevaluation':False,'sourceGuard':{'source':manifest['attestedGuardSources']['prepare'],'line':216}})
    require(len(budget_names)==33 and budget_names[-1]=='new-tail-address-budget-8360-J-27','exact32 prior checkpoints and33rd refusal')
    for rec in budget_records:
        row=rec['record']['row']
        require(rec['name']=='new-tail-address-budget-'+str(row['addressId'])+'-'+row['primitive']+'-'+str(row['K']) and rec['record']['epsilon']=={'text':'1/8000000000000','srepr':'Rational(1, 8000000000000)'} and rec['record']['noOtherRowDonation'] is True,'original full budget address/epsilon arguments')
    require(all(names[a]<names[b] for a,b in zip(budget_names,budget_names[1:])),'original checkpoint order; no claim of new arithmetic pass')
    require(len(records)==5007 and get('journal-result')==failed,'complete original chained result')
    require(not any(n.startswith(('new-positive-tail-total','new-geometry-','numeric-')) for n in names),'only unfinished science ahead')
    old_manifest=read(prior/'source/research/pde_ledger_v3/_measurements/S11c_d_defect_packet_contracted_numeric_inputs.json')
    require(old_manifest['savedInputs']==manifest['savedInputs'],'same3142 original source operands')
    require(get('actual-containment')['durationLimits'] is None,'original contained prefix')
    for alias,rec in manifest['savedInputs'].items():
        require(sha(root/'saved'/alias)==rec['sha256'],'prior argument copy '+alias)
    for source,digest in old_manifest['sourcePins'].items():
        require(sha(source)==digest,'unchanged prior source '+source)
    units=get('posthashes');require(units['copiesIntact'] and failed['copiesIntact'] and failed['sourcePosthashesIntact'] and failed['identityCopiesIntact'],'original integrity flags')
    for source,rec in units['sources'].items():require(rec['intact'] and sha(source)==rec['expected']==rec['actual'],'actual prior source posthash')
    operations=failed['completedOperations'];require(len(operations)==len(set(operations))==263,'all263 complete operations')
    counts={'inherited-new-zero':0,'enlarged-tail-observation':0}
    for name in operations:
        kinds=['input','decision','return']+([] if name.startswith('new-enlarged-tail-') else ['raw'])
        value={k:get(name+'-'+k) for k in kinds}
        J.emit('restored-prefix-'+name,{'operation':name,'published':value,'sourceReceipt':manifest['priorFiles']['complete/'+name+'-return.json'],'originalFunctionCalled':False})
        kind=prefix_return(name,value);counts[kind]+=1
        require(names[name+'-input']<names[name+'-decision']<names[name+'-return'],'original operation order')
    require(counts=={'inherited-new-zero':243,'enlarged-tail-observation':20},'exact prefix census')
    # Every old guard is inventoried, not all branch-dependent calls are claimed
    # executed. Reaching the fixed failed checkpoint attests passed executed
    # instances. Current source/JSON/hash checks join the actual evidence anew.
    inventory={}
    for label,path in manifest['attestedGuardSources'].items():
        text=Path(path).read_text();rows=guard_inventory(text,216 if label=='prepare' else None)
        for row in rows:
            if row['function'] not in manifest['attestedGuardFunctions'][label]:row['disposition']='not on the attested executed prefix call path'
            if label=='restore_library' and row['function'] in ('restore_scalar','walk'):
                row['disposition']='structural parser guards reused only to materialize saved exact rational tail operands; no old scientific calculation'
                row['originalPredicateFunctionReexecuted']=True
        inventory[label]={'sha256':sha(path),'executedCallPathFunctions':manifest['attestedGuardFunctions'][label],'guards':rows}
    require(any(v['disposition']=='failed' and 'per-address per-primitive' in v['source'] for v in inventory['prepare']['guards']),'original guard location')
    J.emit('restored-nonoperation-guard-attestation',{'sources':inventory,'originalFailure':failed['failure'],'chainHead':previous,
        'originalScientificPredicateFunctionsReexecuted':False,'structuralRationalParserGuardsReused':True,'metadataJoinsRerun':['source and copy hashes','chain sequence/full input-return records','actual native address/field/unit/selector maps','source AST contracts and untouched numerical tail'],
        'scientificPredicates':'Successful executed instances depend on original pinned source control flow and published prefix; no new independent arithmetic verification. Branch-dependent unexecuted predicates are not claimed passed.'})
    for name in ('new-original-source-contracts','new-source-ast-join-input','new-profile-tail-source-joins','new-outgoing-source-join','restored-physical-context','restored-original-density-operands','new-Gaussian-envelope-transport','new-per-primitive-absolute-tail-transport'):
        J.emit('restored-attested-'+name,get(name))
    contracts=raw['numeric-source-contracts.json']
    require(get('new-original-source-contracts')==contracts,'actual original source contracts')
    for rec in contracts['fragments'].values():
        text=Path(rec['source']).read_text();nodes=[n for n in ast.walk(ast.parse(text)) if getattr(n,'lineno',None)==rec['line'] and getattr(n,'end_lineno',None)==rec['endLine'] and ast.get_source_segment(text,n)==rec['text']]
        require(len(nodes)==1 and sha(rec['source'])==rec['sourceSha256'],'actual original fragment and hash')
    return get,operations,counts


def prepare(raw,manifest,J,sp,C,G,N):
    get,operations,counts=restore_prefix(raw,manifest,J)
    original=raw.__getitem__;entries=[];all_addresses=original('saved/selected/pressure-addresses.json')['selected']
    require(all_addresses==[a for a in original('saved/inventory/THETA_BALANCE-ordered-addresses.json') if a['row']=='THETA_BALANCE' and a['jet']['channel']=='e_W'],'complete original native row selection')
    fields=original('saved/pressure/fields.json');adapters=original('saved/preflight/numeric-factor-adapters.json')
    for a in all_addresses:
        ident=str(a['addressId']);v=get('restored-address-'+ident)
        ui=original('accepted-units/summand-'+ident+'-input.json');ur=original('accepted-units/summand-'+ident+'-return.json')
        ai=original('contraction/complete/address-'+ident+'-input.json');ar=original('contraction/complete/address-'+ident+'-return.json')
        require(v=={'address':a,'unitInput':ui,'unitReturn':ur,'contractionInput':ai,'contractionReturn':ar},'actual old native address operands')
        require(ai['address']==ui['address']==ui['sourceTransportInput']['address']==a and ai['acceptedUnitReturn']==ur==ar['actualUnitsInherited'],'same source/grade/unit ancestry')
        label='-'.join((a['face'],a['slot'],a['component']))
        require(ai['adapter']==adapters['definitions'][label] and ui['waveProof']['left']==a['waveMultiplier'] and ui['normal']==a['normalMultiplier'],'actual response/wave/normal')
        J.emit('restored-address-context-'+ident,v)
        if a['component'] not in COMPONENTS or a['status'].startswith('EXACT_ZERO'):continue
        e=get('new-eligible-'+ident)
        require(e['addressId']==a['addressId'] and e['face']==a['face'] and e['component']==a['component'] and e['grade']==a['targetGrade'],'eligible original identity')
        require(e['unit']==ur['total']==UNITS and all(v==UNITS for v in ar['totals'].values()),'actual common native summand units before any sum')
        for role in ('source','consumer'):
            fid=a[role+'Transform']['coefficientId'];p=get('new-constant-'+ident+'-'+role+'-input');s=get('new-original-selector-'+ident+'-'+role)
            require(p=={'actualField':fields[fid],'polynomial':original('saved/field/'+fid+'-polynomial.json'),'oldProof':original('saved/field/'+fid+'-reconstruction-input.json'),'address':a},'complete constant field proof arguments')
            require(p['actualField']['constant'] is True and p['polynomial']['degree']==0 and s['actualUnitReturn']==ur and s['address']==a,'inherited constant eligibility')
            require(e['specs'][role]['originalInterface']==s['selectedInterface'] and e['specs'][role]['nativeFieldId']==fid,'actual resumed numerical selector')
        J.emit('restored-live-numerical-entry-'+ident,{'entry':e,'originalAddress':a,'unitArguments':ui,'unitReturn':ur,'noNewCoefficientArithmetic':True});entries.append(e)
    require(len(all_addresses)==544 and len(entries)==20,'full selected/live census')
    J.emit('new-tail-unit-scope',{'entries':entries,'commonUnit':UNITS,'coordinates':'Original native reference-unit magnitudes; analytic1e-11, numerical1e-11/80 and absolutecomparison1e-9 carry this same summand unit.',
        'analyticCeilingIsJDirectSubsetOnly':True,'pendingComponentsRequireSeparateFutureBudgetAndFullLedger':True,'noGlobalBudgetAlreadyEstablished':True})
    baseplan=original('preflight/tail-plan.json');physical=original('preflight/physical-plan.json')
    J.emit('new-independent-tail-domain-input',{'basePlan':baseplan,'physicalPlan':physical,'plannedCarriers':manifest['scope']['plannedCarriers'],'originalContext':get('restored-physical-context'),'nativeAddresses':all_addresses,'entries':entries})
    require((baseplan['K'],baseplan['T'])==(27,122) and baseplan['coversBothCarriers'] is True,'actual original base window and shared carrier envelope')
    require(physical==get('restored-physical-context')['physical']==original('saved/preflight/physical-plan.json')==original('accepted-units/parameter-quantity-origins.json')['physicalPlan'],'same saved physical carrier plan')
    require(manifest['scope']['plannedCarriers']==[v['text'] for v in physical['carriers']]==['sqrt(595)/10','0'] and physical['carriers']==[physical['kappa'],ZERO],'actual two carrier bindings')
    kappa=C.restore_scalar(sp,physical['kappa']);carriers=[C.restore_scalar(sp,v) for v in physical['carriers']]
    domain_flags={'positiveKappa':bool(kappa>0),'kappaBelow3':bool(kappa<3),**{'carrierBelow3/'+str(i):bool(v.is_real is True and abs(v)<3) for i,v in enumerate(carriers)}}
    J.emit('new-native-carrier-domain-decision',{'kappa':kappa,'carriers':carriers,'flags':domain_flags,'newMissingDomainCheckOnly':True})
    require(all(domain_flags.values()),'actual kappa and carrier envelope domain')
    selected=join_base_census(baseplan,all_addresses,entries)
    J.emit('new-independent-native-base-census',{'selectedBaseRecords':selected,'actualEntries':entries,'faceCounts':{'plus':10,'minus':10},'completeNativeRowsJoined':True})
    majorant=get('new-per-primitive-absolute-tail-transport');density={name:original('contraction/complete/new-full-factorization-'+name+'-input.json') for name in ('J','Dr','Dh','Dq')}
    pieces={name:original('absolute-bounds/'+name+'.json') for name in ('D-reflected-numerator-envelope','D-height-numerator-envelope','J-numerator-envelope','global-q-triangle-gap')}
    J.emit('new-inherited-primitive-majorant-input',{'record':majorant,'originalPieces':pieces,'actualDensityInputs':density,'source':manifest['attestedGuardSources']['prepare'],'priorReceipt':manifest['priorFiles']['complete/new-per-primitive-absolute-tail-transport.json']})
    require(majorant['eachPrimitiveUsesFullPositiveDEnvelope'] is True and majorant['notRecomputedGlobalProof'] is True,'actual inherited per-primitive positive D lemma')
    require(majorant['pieces']==pieces and majorant['actualDensityProofs']==density and majorant['domain']==original('saved/pressure/global-parameter-domain.json') and majorant['profile']==original('saved/pressure/profile-envelope.json'),'complete primitive-majorant source operands')
    for name in density:
        require(get('restored-new-full-factorization-'+name)['arguments']==density[name] and get('restored-new-full-factorization-'+name)['return']=={'residual':ZERO},'same inherited factorization input and return')
    J.emit('new-inherited-primitive-majorant-dependency',{'accepted':True,'perPrimitiveFullDOuterAndMiddleDominance':'Inherited triangle/shift-root analytic lemma, not a new magnitude calculation or independent proof.','combinedDEnvelopeUsedThreeTimes':True})
    formula_input={'executedSource':Path(manifest['attestedGuardSources']['prepare']).read_text(),'originalBase':original('numeric-source-contracts.json')['fragments']['tail-contributions']}
    J.emit('new-tail-formula-provenance-input',formula_input)
    formula_join=tail_formula_join(formula_input['executedSource'],formula_input['originalBase']['text'])
    J.emit('new-tail-formula-provenance-decision',formula_join)
    require(all(v['same'] for v in formula_join.values()),'actual original base/enlarged outer J H D AST identity')
    const=lambda record:exact_constant(sp,C,record)
    baseplan=original('preflight/tail-plan.json');base={v['addressId']:v for v in baseplan['allAddresses']}
    context=original('preflight/tail-bound-derivation.json');env=original('preflight/Fourier-envelope-constants.json')['bounds']
    J.emit('new-original-tail-domains-input',{'context':context,'physical':original('preflight/physical-plan.json'),'globalDomain':original('saved/pressure/global-parameter-domain.json'),'profile':original('saved/pressure/profile-envelope.json')})
    require(bool(const(context['b'])==sp.Rational(3000,11101) and const(context['b'])>0),'same positive beta domain')
    bounds={};tail_arguments={};enlarged={};source_records={}
    for e in entries:
        ident=e['addressId'];saved=get('new-enlarged-tail-'+str(ident)+'-input');ret=get('new-enlarged-tail-'+str(ident)+'-return');dec=get('new-enlarged-tail-'+str(ident)+'-decision')
        require(get('restored-base-tail-'+str(ident))==base[ident] and base[ident]['heightQ']==ZERO,'exact original selected base tail')
        require(saved['address']==e and saved['originalConstants']==context and saved['CX']==env[str(ident)]['X'] and saved['CY']==env[str(ident)]['Y'],'actual enlarged arguments to same entry/constants')
        require(saved['originalSource']==original('numeric-source-contracts.json')['fragments']['tail-contributions'],'original positive expression provenance')
        require(dec['includesPositiveHOvercount']==(e['component']==COMPONENTS[0]),'saved H overcount label')
        require(bool(const(saved['CX'])>0 and const(saved['CY'])>0),'positive native envelope arguments')
        J.emit('new-restored-tail-arguments-'+str(ident),{'base':base[ident],'enlargedInput':saved,'enlargedReturn':ret,'baselineOnlyMeaning':'unmutated kernel; both windows identified explicitly','reEvaluatedClosedForm':False})
        enlarged[ident]={'input':saved,'return':ret};source_records[ident]={27:base[ident],29:ret}
        pair={27:{k:const(base[ident][k]) for k in ('outer','middle')},29:{k:const(ret[k]) for k in ('outer','middle')}}
        flags={str(K)+'/'+k:bool(v>=0) for K,x in pair.items() for k,v in x.items()}
        mono={k:bool(pair[29][k]<=pair[27][k]) for k in ('outer','middle')}
        J.emit('new-tail-component-regression-'+str(ident),{'actualPairs':{str(k):v for k,v in pair.items()},'nonnegative':flags,'componentwiseMonotonicity':mono,'classification':'argument/provenance regression only; not independent arithmetic or physics certificate'})
        require(all(flags.values()) and all(mono.values()),'saved nonnegative component monotonicity')
        require(baseplan['T']>=baseplan['K']+4 and saved['T']>=saved['K']+4 and (saved['K'],saved['T'])==tuple(manifest['scope']['originalWindowPairs'][1])==(29,124) and tuple(manifest['scope']['originalWindowPairs'][0])==(baseplan['K'],baseplan['T']),'both actual saved window domains')
        bounds[ident]=pair;tail_arguments[str(ident)]=saved
    # The separate H addend was not persisted. Derive only its SOURCE-BOUND
    # scaling relation from the original exponential factor, not its magnitude
    # or either already completed combined-middle value.
    fragment=original('numeric-source-contracts.json')['fragments']['tail-contributions']
    tree=ast.parse(fragment['text']);terms=[n.value for n in ast.walk(tree) if isinstance(n,ast.AugAssign) and isinstance(n.target,ast.Name) and n.target.id=='middle']
    require(len(terms)==1,'one original H overcount source term')
    reference=ast.parse('(4/b)*cx*cy*E30*F[2]**2*sp.Rational(55,3)*sp.Rational(1,2)**Tlim',mode='eval').body
    require(ast.dump(terms[0])==ast.dump(reference),'exact separate H exponential source form')
    J.start('new-H-overcount-scaling',{'originalSource':fragment,'HExpression':ast.unparse(terms[0]),'actualTailArguments':tail_arguments,'baseT':122,'enlargedT':124,'missingMagnitudeNotRecomputed':True})
    h=sp.Symbol('unchanged_positive_H_prefactor',positive=True);left=h*sp.Rational(1,2)**124;right=h*sp.Rational(1,2)**122/4;residual=sp.cancel(left-right)
    J.emit('new-H-overcount-scaling-decision',{'left':left,'right':right,'residual':residual,'ratio':sp.Rational(1,4),'positiveCoefficientPremise':'same4/b,CX,CY,3^30,F2^2,55/3; exact saved positive constants checked','notIndependentTailReevaluation':True})
    require(bool(residual is sp.S.Zero and const(context['fullExponentialMoments']['2'])>0),'separate H-source scaling and positive common factor');J.finish({'residual':residual,'classification':'new source-bound regression relation only'})
    ceiling=sp.Rational(1,10**11);tails=[];all_groups=[]
    encoded=lambda v:{'text':str(v),'srepr':sp.srepr(v)}
    for K,T in ((27,122),(29,124)):
        for carrier in manifest['scope']['plannedCarriers']:
            expected=independent_slots(baseplan,enlarged,all_addresses,K,T,carrier);rows=[]
            for spec in canonical_slots(entries,K,T,carrier):
                ident=spec['addressId'];pair=bounds[ident][K]
                rows.append({**spec,**pair,'boundSource':bound_source(ident,K,T,source_records[ident][K]),'outerOperand':encoded(pair['outer']),'middleOperand':encoded(pair['middle']),'total':pair['outer']+pair['middle'],'baselineOnly':True,'includesExtraHOvercount':spec['primitive']=='J'})
            label=str(K)+'-'+('zero' if carrier=='0' else 'matching')
            J.emit('new-analytic-tail-census-'+label,{'expected':expected,'actual':rows,'numericalEpsilonUnchanged':sp.Rational(1,80*10**11)})
            check_census(rows,expected)
            # Same census path; mutations have no second baseline arithmetic.
            mutations={'drop':copy.deepcopy(rows[:-1]),'duplicate':copy.deepcopy(rows+[rows[0]])}
            for tag,key,value in [('wrong-T','T',T+1),('wrong-K','K',29 if K==27 else 27),('wrong-carrier','carrier',manifest['scope']['plannedCarriers'][1] if carrier==manifest['scope']['plannedCarriers'][0] else manifest['scope']['plannedCarriers'][0]),('wrong-face','face','minus' if rows[0]['face']=='plus' else 'plus')]:
                v=copy.deepcopy(rows);v[0][key]=value;mutations[tag]=v
            value_mutant=copy.deepcopy(rows);other=bounds[rows[0]['addressId']][29 if K==27 else 27]
            value_mutant[0].update(other,outerOperand=encoded(other['outer']),middleOperand=encoded(other['middle']),total=other['outer']+other['middle'])
            mutations['swapped-window-values']=value_mutant
            address_mutant=copy.deepcopy(rows);other_row=next(v for v in rows if v['addressId']!=rows[0]['addressId'] and (v['outerOperand'],v['middleOperand'])!=(rows[0]['outerOperand'],rows[0]['middleOperand']))
            for key in ('outer','middle','total','outerOperand','middleOperand'):address_mutant[0][key]=copy.deepcopy(other_row[key])
            mutations['swapped-address-values']=address_mutant
            for tag,mutant in mutations.items():
                J.emit('new-census-control-'+label+'-'+tag+'-input',{'actualExpected':expected,'mutant':mutant})
                error=None
                try:check_census(mutant,expected)
                except ValueError as ex:error=str(ex)
                J.emit('new-census-control-'+label+'-'+tag+'-return',{'refused':error is not None,'reason':error,'noNumericResponseClaim':True})
                require(error is not None,'responsive same-path coverage control')
            groups=group_indices(rows);sums={name:sum((rows[i]['total'] for i in indices),sp.S.Zero) for name,indices in groups.items()}
            flags={name:bool(value>=0 and value<ceiling) for name,value in sums.items()}
            J.emit('new-separate-analytic-tail-sums-'+label,{'rows':rows,'groups':groups,'sums':sums,'ceiling':ceiling,'unit':UNITS,'decisions':flags,'carriersNotAdded':True,'windowsNotAdded':True,'eachDCountedThreeTimes':True,'noDonationOrCancellation':True,'subsetOnly':True})
            all_groups.append({'label':label,'flags':flags,'sumRecord':'new-separate-analytic-tail-sums-'+label})
            tails.extend(rows)
    J.emit('new-all-analytic-tail-group-decisions',{'groups':all_groups,'carrierGroupsUseSameUniformEnvelopeNotIndependentEvidence':True,'allFourRecordsPersistedBeforeBudgetGuard':True})
    require(len(all_groups)==4 and all(all(v['flags'].values()) for v in all_groups),'actual separate analytic subset aggregate budget')
    J.emit('new-analytic-tail-acceptance',{'scope':'Only fixed J/Dr/Dh/Dq subset','fourWindowCarrierGroupsPassed':True,'tailRows':len(tails),'wholePacketBudget':None,'pendingComponents':'H/flat/heightPV/slope require future independent positive allocations and complete error ledger before full action acceptance','empiricalErrorNotRigorous':True})
    return continue_geometry(raw,manifest,J,sp,C,G,N,entries,tails)
