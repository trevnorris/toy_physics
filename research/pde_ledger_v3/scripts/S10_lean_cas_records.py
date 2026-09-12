"""Bind full coincidence records and spectral status metadata to certified objects.

Association keys are unique and exhaustive. Every repeated field is translated
from its own payload. Unresolved CAS signs remain unresolved observations even
when an independent Lean theorem determines the mathematical sign.
"""
import re
import S10_lean_cas_coincidence as co
import S10_lean_cas_roots as roots
import S10_lean_cas_reruns as reruns

CASES = roots.CASES
FIELDS = ('Pair','Equation','Locus','Allowed','IntersectionOutcome','Intersects')


def aggregate_suffixes(case):
    if case=='generic':
        return ['Q3_ROOT_COINCIDENCE_LOCI',*[f'ROOT{i}_Q8_ROOT_COINCIDENCE_LOCI' for i in range(1,4)]]
    return [reruns.prefix('WL',case)+'Q3_STRATUM_ROOT_COINCIDENCE_RECORDS']


def associations(bridge,tree,case):
    fields=FIELDS if case=='generic' else ('Pair','Equation','ParameterLocus','AllowedParameterRegion','IntersectionOutcome','Intersects')
    result=[]
    for obj in co.seq(bridge,tree,1 if case=='parallel' else 3):
        bridge.require(obj[0]=='assoc','expected a coincidence association')
        values={}
        for rule in obj[1:]:
            bridge.require(rule[0]=='rule' and rule[1][0]=='key','expected a named association field')
            key=rule[1][1]
            bridge.require(key in fields and key not in values,'unknown or duplicate association field')
            values[key]=rule[2]
        bridge.require(set(values)==set(fields),'missing coincidence association field')
        pair=co.seq(bridge,values['Pair'],2)
        bridge.require(all(x[0]=='num' and 1<=x[1]<=3 for x in pair),'pair indices must be integers in 1..3')
        eq=values['Equation']
        bridge.require(eq[0]=='eq' and eq[2]==['num',0],'aggregate equation must compare to literal zero')
        co.node(bridge,eq[1],'WL')
        co.prop(bridge,co.locus(bridge,values[fields[2]]))
        co.prop(bridge,values[fields[3]])
        bridge.require(values['IntersectionOutcome'] in (['symbol','decidedEmpty'],['symbol','decidedNonempty']),
                       'unsupported aggregate outcome')
        bridge.require(values['Intersects'] in (['symbol','True'],['symbol','False']),'expected aggregate Boolean decision')
        result.append(values)
    return result


def build(bridge,engine,path,case):
    generic=case=='generic';point='k' if generic else engine+'Loci.'+case+'Point'
    dims=bridge.DIMS if generic else {**bridge.DIMS,'z':(2,-2,0),**{f'k{i}':(0,0,0) for i in range(3)}}
    b=bridge.Builder(engine,path,namespace=engine+'Records'+case.title(),support='RecordSupport',
                     tactic='cas_eval' if generic else 'rerun_eval',units='units' if generic else 'coordinateUnits',dims=dims)
    if engine=='WL':
        pairs=co.PAIRS[:1] if case=='parallel' else co.PAIRS
        cells=[]
        for i,j in pairs:
            cells += [(str(i+1),(0,0,0),[]),(str(j+1),(0,0,0),[]),
                      (f'referenceRoot rho mu sigma {point} {i} - referenceRoot rho mu sigma {point} {j}',
                       (0,-2,0) if generic else (2,-2,0),['rho','sigma'])]
        for suffix in aggregate_suffixes(case):
            def extract(raw,e):
                values=associations(bridge,co.parse(bridge,raw,e),case)
                return [co.node(bridge,t,e) for v in values for t in [*co.seq(bridge,v['Pair'],2),v['Equation'][1]]]
            b.record(suffix,lambda x:x,cells,parse_payload=extract)
            b.records[-1]['record_semantics']={'kind':'aggregate','case':case,'field_order':list(FIELDS if generic else ('Pair','Equation','ParameterLocus','AllowedParameterRegion','IntersectionOutcome','Intersects')),
                                             'cell_order':'two pair indices followed by root difference, per pair'}
    else:
        def operand(raw,e):
            parsed=bridge.Parser(raw,e).parse()
            bridge.require(isinstance(parsed,list) and len(parsed)==2 and parsed[1]==bridge.Node('atom',('z',)),
                           'spectrum operands require one expression and the exact omegaSquared variable')
            return parsed
        suffix=('' if generic else reruns.prefix(engine,case))+'Q3_SPECTRUM_SOLVE_CONDITION_OPERANDS'
        b.record(suffix,lambda x:x,
                 [(f'Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma z {point})',(-9,-6,3) if generic else (-3,-6,3),[]),
                  ('z',(0,-2,0) if generic else (2,-2,0),[])],parse_payload=operand)
        b.records[-1]['record_semantics']={'kind':'spectrum_operands','case':case,'variable':'omegaSquared'}
    return b


def field_proofs(cell):
    domain=cell['raw_denominator_factors']+[['atom',x] for x in cell['reference_domain_symbols']]
    return ' '.join('(by positivity)' for _ in dict.fromkeys(map(str,domain)))


def reported_sign(bridge,raw):
    table={'0':'zero','1':'positive','-1':'negative','undecided_under_joint_assumptions':'undecided'}
    bridge.require(raw in table,'unsupported reported sign')
    return table[raw]


def generate(bridge):
    outputs,engines,records,audits={},{},[],[]
    builders={(e,case):build(bridge,e,p,case) for e,p in bridge.INPUTS.items() for case in CASES}
    for (e,case),b in builders.items():
        outputs[bridge.GENERATED/(b.namespace+'.lean')]=b.output()
        engines[b.namespace]={'path':str(b.path.relative_to(bridge.BASE)),'sha256':bridge.sha(b.path.read_bytes()),
                              'records':b.records,'unique_trees':len(b.nodes)}
        audits += [b.namespace+'.'+a for a in b.audits]
    header=['import S10Audit.CAS.RecordSupport',*['import S10Audit.CAS.'+b.namespace for b in builders.values()],
            '', 'set_option backward.isDefEq.respectTransparency false','','namespace S10Audit.CAS.Records',
            'open S10Pilot S10Anisotropic Polynomial','noncomputable section','']
    lines=header.copy();status_lines=header.copy()
    params='(rho mu sigma lam : ℝ) (a k : Vec 3)';args='rho mu sigma lam a k'
    hd='(hd : coincidenceDomain rho mu sigma lam)'
    def entry(engine,suffix,kind):
        tag=bridge.PREFIX+suffix;rows=builders[engine,'generic'].rows
        bridge.require(tag in rows,'missing metadata record '+tag)
        line,raw=rows[tag]
        rec={'tag':engine+'_'+tag,'line':line,'payload':raw,'payload_sha256':bridge.sha(raw.encode()),
             'kind':kind,'components':[],'lean_semantic_claims':[],
             'domain':'D=3; rho,mu,sigma,lambdaScale>0; sigma != 1; real k,a; explicit nonzero k when needed'}
        records.append(rec)
        return raw,rec
    def claim(rec,name):
        full='Records.'+name;rec['lean_semantic_claims'].append(full);audits.append(full)
    def translated(rec,name,tree,primary,expected):
        lines.extend([f'def {name} {params} : Prop := {co.prop(bridge,tree)}','',
                      f'theorem {name}_primary {params} : {name} {args} ↔ Coincidence.{primary} {args} := by',
                      '  rfl','',f'theorem {name}_correct {params} {hd} : {name} {args} ↔ {expected} := by',
                      f'  exact ({name}_primary rho mu sigma lam a k).trans',
                      f'    (Coincidence.{primary}_correct rho mu sigma lam a k hd)',''])
        rec['components'].append({'predicate_ast':tree,'lean_predicate':'Records.'+name,'primary':primary,'expected':expected})
        for suffix in ('primary','correct'): claim(rec,name+'_'+suffix)
    for root in range(1,4):
        for suffix,kind in [('LOCI','locus'),('ALLOWED_OPERANDS','allowed'),('ALLOWED_TESTS','status')]:
            raw,rec=entry('PY',f'ROOT{root}_Q8_ROOT_COINCIDENCE_'+suffix,'q8_duplicate')
            tree=co.parse(bridge,raw,'PY');rec['mixed_ast']=tree
            for index,item in enumerate(co.seq(bridge,tree,3)):
                if kind=='locus': pt=co.locus(bridge,item)
                elif kind=='status': pt=item
                else:
                    loc,positive,assumptions=co.seq(bridge,item,3)
                    pt=co.junction('and',[co.locus(bridge,loc),positive,assumptions])
                expected=f'coincidenceGeometry {index} k' if kind=='locus' else f'allowedCoincidence {index} k' if kind=='allowed' else f'(∃ w : Vec 3, allowedCoincidence {index} w)'
                translated(rec,f'PY_q8root{root}_{kind}{index}',pt,f'PY_generic_{kind}{index}',expected)
    for case in CASES:
        generic=case=='generic';b=builders['WL',case];point='k' if generic else 'WLLoci.'+case+'Point'
        for number,suffix in enumerate(aggregate_suffixes(case)):
            raw,rec=entry('WL',suffix,'coincidence_aggregate');tree=co.parse(bridge,raw,'WL');rec['mixed_ast']=tree
            fields=associations(bridge,tree,case);ar=b.records[number]
            base=f'WL_{case}_aggregate{number}'
            for index,v in enumerate(fields):
                for key,kind in [('Locus' if generic else 'ParameterLocus','locus'),
                                 ('Allowed' if generic else 'AllowedParameterRegion','allowed'),
                                 ('IntersectionOutcome','outcome'),('Intersects','status')]:
                    pt=co.locus(bridge,v[key]) if kind=='locus' else v[key]
                    if kind=='outcome': pt=['symbol','True' if v[key][1]=='decidedNonempty' else 'False']
                    expected=(f'coincidenceGeometry {index} k' if kind=='locus' else f'allowedCoincidence {index} k' if kind=='allowed' else f'(∃ w : Vec 3, allowedCoincidence {index} w)') if generic else 'False'
                    translated(rec,base+f'_{kind}{index}',pt,f'WL_{case}_{kind}{index}',expected)
                c=ar['cells'][index*3+2];name=base+f'_equation{index}'
                lines += [f'def {name} (rho mu sigma z : ℝ) (k : Vec 3) : Prop :=',
                          f'  Expr.eval (values rho mu sigma z k) {b.namespace}.{c["lean_tree"]} = 0','',
                          f'theorem {name}_locus (rho mu sigma lam z : ℝ) (a k : Vec 3) {hd} :',
                          f'    {name} rho mu sigma z k ↔ {base}_locus{index} {args} := by',
                          '  have hr := hd.1', '  have hs := hd.2.2.1',
                          f'  rw [{name}, {b.namespace}.{c["lean_claim"]} rho mu sigma z k {field_proofs(c)},',
                          f'    sub_eq_zero, Coincidence.WL_{case}_pair{index}_roots rho mu sigma lam k hd,',
                          f'    {base}_locus{index}_correct rho mu sigma lam a k hd]','']
                claim(rec,name+'_locus');ar['record_semantics'].setdefault('lean_semantic_claims',[]).append('Records.'+name+'_locus')
                rec['components'].append({'field':'Equation','arithmetic_cell':index*3+2,'lean_predicate':'Records.'+name})
                for kind in ('status','outcome'):
                    name=base+f'_{kind}{index}_decides'
                    lines += [f'theorem {name} {params} {hd} :',
                              f'    {base}_{kind}{index} {args} ↔ ∃ w : Vec 3, {base}_allowed{index} rho mu sigma lam a w := by',
                              f'  simp_rw [{base}_{kind}{index}_correct rho mu sigma lam a k hd,',
                              f'    {base}_allowed{index}_correct rho mu sigma lam a _ hd]']
                    if not generic: lines += ['  simp']
                    lines += [''];claim(rec,name)
            pairs=[co.seq(bridge,v['Pair'],2) for v in fields]
            literal='['+', '.join(f'({p[0][1]}, {p[1][1]})' for p in pairs)+']'
            reference='[(1, 2)]' if case=='parallel' else '[(1, 2), (1, 3), (2, 3)]'
            lines += [f'def {base}_pairs : List (ℕ × ℕ) := {literal}',
                      f'theorem {base}_pairs_complete : {base}_pairs = {reference} := by decide','']
            claim(rec,base+'_pairs_complete')
            rec['components'].append({'field':'Pair','lean_pairs':'Records.'+base+'_pairs','expected_pairs':[[i+1,j+1] for i,j in (co.PAIRS[:1] if case=='parallel' else co.PAIRS)]})
    # Sign observations and empty solver-condition lists are bound to each actual root expression.
    for engine,path in bridge.INPUTS.items():
        for case in CASES:
            generic=case=='generic';pre='' if generic else reruns.prefix(engine,case)
            point='k' if generic else engine+'Loci.'+case+'Point'
            rb=roots.build(bridge,engine,path,case);rn=rb.namespace
            root_cells=rb.records[1]['cells'];nonzero='hk' if generic else point+'_nonzero'
            for index,c in enumerate(root_cells):
                name=f'{engine}_{case}_root{index}_sign'
                raw,rec=entry(engine,pre+f'ROOT{index+1}_Q3_SIGN','root_sign')
                sign=reported_sign(bridge,raw);rec['reported_sign']=sign;rec['computed_sign']='zero' if index==0 else 'positive'
                root=f'Expr.eval (values rho mu sigma z k) {rn}.{c["lean_tree"]}'
                status_lines += [f'def {name} : ReportedSign := .{sign}',
                                 f'theorem {name}_observed : {name} = .{sign} := rfl','',
                                 f'theorem {name}_computed (rho mu sigma lam z : ℝ) (k : Vec 3) {hd} '+('(hk : k ≠ 0)' if generic else '')+' :',
                                 f'    (classifiedSign {index}).Holds ({root}) := by',
                                 '  have hr := hd.1', '  have hm := hd.2.1', '  have hs := hd.2.2.1',
                                 f'  rw [{rn}.{c["lean_claim"]} rho mu sigma z k {field_proofs(c)}]',
                                 f'  exact reference_root_sign rho mu sigma {point} {index} hr hm hs {nonzero}','',
                                 f'theorem {name}_sound (rho mu sigma lam z : ℝ) (k : Vec 3) {hd} '+('(hk : k ≠ 0)' if generic else '')+' :',
                                 f'    {name}.Holds ({root}) := by']
                if sign=='undecided': status_lines += ['  trivial','']
                else: status_lines += [f'  exact {name}_computed rho mu sigma lam z k hd'+(' hk' if generic else ''),'']
                for kind in ('observed','computed','sound'): claim(rec,name+'_'+kind)
                if engine=='WL':
                    raw,rec=entry(engine,pre+f'ROOT{index+1}_Q3_SOLVER_CONDITIONS','solver_conditions')
                    t=co.parse(bridge,raw,engine);co.seq(bridge,t,0);rec['mixed_ast']=t
                    sn=f'{engine}_{case}_root{index}_conditions'
                    status_lines += [f'def {sn} : List Prop := []',
                                     f'theorem {sn}_empty : {sn} = [] := rfl','',
                                     f'theorem {sn}_root (rho mu sigma lam z : ℝ) (k : Vec 3) {hd} :',
                                     f'    Matrix.det ((1/2 : ℝ) • referenceMatrix rho mu sigma ({root}) {point}) = 0 := by',
                                     '  have hr := hd.1','  have hs := hd.2.2.1',
                                     f'  rw [{rn}.{c["lean_claim"]} rho mu sigma z k {field_proofs(c)}]',
                                     f'  exact reference_root_isRoot rho mu sigma {point} {index} (ne_of_gt hr) (ne_of_gt hs)','']
                    claim(rec,sn+'_empty');claim(rec,sn+'_root')
            if engine=='PY':
                raw,rec=entry(engine,pre+'Q3_SPECTRUM_SOLVE_CONDITION','spectrum_status')
                bridge.require(raw=='roots_returned','unsupported spectrum solve disposition')
                b=builders[engine,case];ar=b.records[0];c=ar['cells'][0]
                name=f'PY_{case}_spectrum_returned'
                status_lines += [f'def {name} : Bool := true',
                                 f'theorem {name}_nonempty (rho mu sigma lam z : ℝ) (k : Vec 3) {hd} :',
                                 f'    {name} = true ↔ {rn}.candidates rho mu sigma z k ≠ [] := by',
                                 f'  rw [{rn}.candidates_reference rho mu sigma z k (ne_of_gt hd.1) (ne_of_gt hd.2.2.1)]',
                                 f'  simp [{name}]','',
                                 f'theorem {name}_complete (rho mu sigma lam z : ℝ) (k : Vec 3) {hd} :',
                                 f'    Expr.eval (values rho mu sigma z k) {b.namespace}.{c["lean_tree"]} = 0 ↔',
                                 f'      z ∈ {rn}.candidates rho mu sigma z k := by',
                                 f'  rw [{b.namespace}.{c["lean_claim"]} rho mu sigma z k]']
                if generic:
                    status_lines += ['  rw [rootPolynomial_complete rho mu sigma z k (ne_of_gt hd.1) (ne_of_gt hd.2.2.1),',
                                     '    rootPolynomial_roots rho mu sigma k (ne_of_gt hd.1) (ne_of_gt hd.2.2.1),',
                                     f'    {rn}.candidates_reference rho mu sigma z k (ne_of_gt hd.1) (ne_of_gt hd.2.2.1)]',
                                     '  simp [referenceRoot]','']
                else:
                    status_lines += [f'  exact {rn}.candidates_complete rho mu sigma z k (ne_of_gt hd.1) (ne_of_gt hd.2.1)',
                                     '    hd.2.2.1 hd.2.2.2.1 z','']
                for kind in ('nonempty','complete'): claim(rec,name+'_'+kind)
                ar['record_semantics']['lean_semantic_claims']=['Records.'+name+'_complete']
    for case,index in [('parallel',1),('perpendicular',2)]:
        raw,rec=entry('PY',f'Q8_STRATUM{index}_SKIP_STATUS','stratum_status')
        bridge.require(raw=='not_skipped_allowed_branch','unsupported stratum skip disposition')
        name='PY_'+case+'_not_skipped';point='PYLoci.'+case+'Point'
        status_lines += [f'def {name} : Bool := true',f'theorem {name}_allowed :',
                         f'    {name} = true ∧ {point} ≠ 0 ∧ PYLoci.extraTransverse {point} :=',
                         f'  ⟨rfl, {point}_nonzero, {point}_target⟩','']
        claim(rec,name+'_allowed')
    # Record and justify the origin stratum that the SymPy audit deliberately skips.
    raw,rec=entry('PY','Q8_SKIPPED_STRATUM1_BRANCH','skipped_branch')
    tree=co.parse(bridge,raw,'PY');rec['mixed_ast']=tree
    pt=co.junction('and',[co.assignment(bridge,x) for x in co.seq(bridge,tree)])
    status_lines += [f'def PY_skipped_branch (k : Vec 3) : Prop := {co.prop(bridge,pt)}',
                     'theorem PY_skipped_branch_geometry (k : Vec 3) : PY_skipped_branch k ↔ k = 0 := by',
                     '  simp only [PY_skipped_branch, wavevector_zero_iff]',
                     '  tauto','']
    claim(rec,'PY_skipped_branch_geometry')
    rec['components'].append({'predicate_ast':pt,'lean_predicate':'Records.PY_skipped_branch'})
    raw,rec=entry('PY','Q8_SKIPPED_STRATUM1_ALLOWED_TEST','skipped_decision')
    bridge.require(raw in ('True','False'),'expected skipped-stratum Boolean')
    status_lines += [f'def PY_skipped_test : Bool := {raw.lower()}',
                     'theorem PY_skipped_test_decides : PY_skipped_test = true ↔',
                     '    ∃ k : Vec 3, PY_skipped_branch k ∧ 0 < normSq k := by',
                     '  simp [PY_skipped_test, PY_skipped_branch_geometry, normSq_positive_iff]','']
    claim(rec,'PY_skipped_test_decides')
    raw,rec=entry('PY','Q8_SKIPPED_STRATUM1_REASON','skipped_reason')
    reason='locus_conflicts_with_positive_wavevector_norm'
    bridge.require(raw==reason,'unsupported skipped-stratum reason')
    status_lines += ['theorem PY_skipped_reason_sound (k : Vec 3) (h : PY_skipped_branch k) : normSq k = 0 := by',
                     '  rw [PY_skipped_branch_geometry] at h', '  subst k','  simp [normSq, dot]','']
    claim(rec,'PY_skipped_reason_sound')
    raw,rec=entry('PY','Q8_SKIPPED_STRATA','skipped_aggregate')
    tree=co.Parser(bridge,raw,'PY',extra_symbols=(reason,)).parse();rec['mixed_ast']=tree
    branch,test,why=co.seq(bridge,co.seq(bridge,tree,1)[0],3)
    bridge.require(test in (['symbol','True'],['symbol','False']) and why==['symbol',reason],'unsupported skipped aggregate disposition')
    pt=co.junction('and',[co.assignment(bridge,x) for x in co.seq(bridge,branch)])
    status_lines += [f'def PY_skipped_aggregate_branch (k : Vec 3) : Prop := {co.prop(bridge,pt)}',
                     f'def PY_skipped_aggregate_test : Bool := {test[1].lower()}',
                     'theorem PY_skipped_aggregate_branch_primary (k : Vec 3) :',
                     '    PY_skipped_aggregate_branch k ↔ PY_skipped_branch k := Iff.rfl',
                     'theorem PY_skipped_aggregate_test_decides : PY_skipped_aggregate_test = true ↔',
                     '    ∃ k : Vec 3, PY_skipped_aggregate_branch k ∧ 0 < normSq k := by',
                     '  simp [PY_skipped_aggregate_test, PY_skipped_aggregate_branch_primary,',
                     '    PY_skipped_branch_geometry, normSq_positive_iff]',
                     'theorem PY_skipped_aggregate_reason_sound (k : Vec 3) (h : PY_skipped_aggregate_branch k) : normSq k = 0 :=',
                     '  PY_skipped_reason_sound k ((PY_skipped_aggregate_branch_primary k).mp h)','']
    for kind in ('branch_primary','test_decides','reason_sound'): claim(rec,'PY_skipped_aggregate_'+kind)
    for module,body in [('RecordBindings',lines),('StatusBindings',status_lines)]:
        body += ['end','end S10Audit.CAS.Records','']
        outputs[bridge.GENERATED/(module+'.lean')]=re.sub(r'\b(rho|mu|sigma|lam|a|k|hd|hk)\b',r'_\1','\n'.join(body))
    audits += ['reference_root_sign','reference_root_isRoot']
    return outputs,engines,records,audits
