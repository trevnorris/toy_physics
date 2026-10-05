"""Preserve coincidence predicates, guards, pair indices and allowed witnesses.

The small mixed grammar is non-evaluating. It specializes the transcript's
dimension variable to D=3 and represents real membership by the Lean types.
All arithmetic, logical branches and conditional-rule guards remain explicit.
"""
import re
import S10_lean_cas_reruns as reruns

PAIRS = [(0,1),(0,2),(1,2)]


class Parser:
    token = re.compile(r'\s+|<\||\|>|->|\*\*|&&|\|\||==|!=|>=|<=|"[A-Za-z]+"|[A-Za-z_][A-Za-z_0-9]*(?:\.[a-z]+)?|[0-9]+|[+*/^(),{}\[\]<>&|=-]')
    precedence = {'->':0,'||':1,'|':1,'&&':2,'&':2,'==':3,'!=':3,'>':3,'<':3,'>=':3,'<=':3,
                  '+':4,'-':4,'*':5,'/':5,'^':7,'**':7}
    operations = {'->':'rule','||':'or','|':'or','&&':'and','&':'and','==':'eq','!=':'ne',
                  '>':'gt','<':'lt','>=':'ge','<=':'le','+':'add','-':'sub','*':'mul','/':'div','^':'pow','**':'pow'}

    def __init__(self, bridge, raw, engine, *, extra_symbols=()):
        self.b,self.engine,self.tokens,self.pos = bridge,engine,[],0
        self.extra_symbols = frozenset(extra_symbols)
        bridge.require(engine in ('PY','WL') and len(raw)<100000,'coincidence parser domain')
        end=0
        for m in self.token.finditer(raw):
            bridge.require(m.start()==end,'unsupported coincidence syntax')
            end=m.end()
            if not m.group().isspace(): self.tokens.append(m.group())
        bridge.require(end==len(raw),'unsupported coincidence trailing syntax')

    def peek(self): return self.tokens[self.pos] if self.pos<len(self.tokens) else None

    def take(self, expected=None):
        t=self.peek()
        self.b.require(t is not None and (expected is None or t==expected),f'expected {expected}, got {t}')
        self.pos+=1
        return t

    def sequence(self, close):
        items=[]
        while self.peek()!=close:
            if items:
                self.take(',')
                if self.peek()==close: break
            items.append(self.expression())
        self.take(close)
        return items

    def expression(self, minimum=0):
        t=self.take()
        if t in ('+','-'):
            a=self.expression(6)
            left=a if t=='+' else ['neg',a]
        elif t in ('(', '[','{','<|'):
            self.b.require(t!='[' or self.engine=='PY','unexpected Wolfram square bracket')
            self.b.require(t not in ('{','<|') or self.engine=='WL','unexpected Python container')
            close={'(':')','[':']','{':'}','<|':'|>'}[t]
            if t=='(':
                if self.peek()==')': self.take(')'); left=['seq']
                else:
                    first=self.expression()
                    if self.peek()==')': self.take(')'); left=first
                    else:
                        self.take(','); left=['seq',first,*self.sequence(')')]
            else: left=['assoc' if t=='<|' else 'seq',*self.sequence(close)]
        elif t.isdigit():
            self.b.require(len(t)<=30,'oversized coincidence integer')
            left=['num',int(t)]
        elif t.startswith('"'):
            self.b.require(t[1:-1] in ('Pair','Equation','Locus','Allowed','IntersectionOutcome','Intersects',
                                      'ParameterLocus','AllowedParameterRegion'),'unsupported association key')
            left=['key',t[1:-1]]
        elif self.peek()==('(' if self.engine=='PY' else '['):
            allowed=('Eq','Ne','Q.integer','Q.real','Q.positive') if self.engine=='PY' else ('Element','Inequality','Sqrt','ConditionalExpression')
            self.b.require(t in allowed,'unsupported coincidence function '+t)
            self.take()
            left=['call',t,*self.sequence(')' if self.engine=='PY' else ']')]
        else:
            names={'rho_br','mu_R','s_rho','rhoBr','muR','sRho','lambdaScale','D','braneDimension',
                   'k1','k2','k3','a1','a2','a3','True','False','EmptySet','Reals','Integers',
                   'Less','LessEqual','Greater','GreaterEqual','decidedEmpty','decidedNonempty'}
            self.b.require(t in names | self.extra_symbols,'unsupported coincidence symbol '+t)
            left=['symbol',t]
        while self.peek() in self.precedence and self.precedence[self.peek()]>=minimum:
            op=self.take(); prec=self.precedence[op]
            if op in ('**','^'):
                self.b.require(op==('**' if self.engine=='PY' else '^'),'wrong power syntax')
            if op=='->': self.b.require(self.engine=='WL','unexpected Python rule')
            if op in ('&&','||','&','|'):
                self.b.require((op in ('&&','||'))==(self.engine=='WL'),'wrong logical operator syntax')
            right=self.expression(prec if op in ('**','^') else prec+1)
            left=[self.operations[op],left,right]
        return left

    def parse(self):
        result=self.expression()
        self.b.require(self.peek() is None,'unconsumed coincidence tokens')
        return result


def seq(bridge, tree, length=None):
    bridge.require(tree[0]=='seq' and (length is None or len(tree)==length+1),'wrong coincidence sequence shape')
    return tree[1:]


def term(bridge, t):
    op,*a=t
    if op=='num': return str(a[0])
    if op=='symbol':
        names={'rho_br':'rho','rhoBr':'rho','mu_R':'mu','muR':'mu','s_rho':'sigma','sRho':'sigma',
               'lambdaScale':'lam','D':'(3 : ℝ)','braneDimension':'(3 : ℝ)',
               **{f'k{i+1}':f'k {i}' for i in range(3)},**{f'a{i+1}':f'a {i}' for i in range(3)}}
        bridge.require(a[0] in names,'non-arithmetic coincidence symbol')
        return names[a[0]]
    if op=='neg': return '(-'+term(bridge,a[0])+')'
    if op=='call':
        bridge.require(a[0]=='Sqrt' and len(a)==2,'non-arithmetic coincidence function')
        return '(Real.sqrt '+term(bridge,a[1])+')'
    bridge.require(op in ('add','sub','mul','div','pow'),'non-arithmetic coincidence tree')
    if op=='pow': bridge.require(a[1][0]=='num' and 0<=a[1][1]<=32,'unbounded coincidence power')
    if op=='div':
        denominator=node(bridge,a[1],None)
        bridge.numeric(bridge.Node('div',(bridge.Node('num',(1,)),denominator)))
        for factor in bridge.denominators(bridge.Node('div',(bridge.Node('num',(1,)),denominator))):
            bridge.require(factor.op=='atom' and factor.args[0] in ('rho','mu','sigma'),'denominator outside declared coefficient domain')
    return '('+term(bridge,a[0])+{'add':' + ','sub':' - ','mul':' * ','div':' / ','pow':' ^ '}[op]+term(bridge,a[1])+')'


def prop(bridge,t):
    op,*a=t
    if op=='symbol':
        bridge.require(a[0] in ('True','False','EmptySet'),'non-logical coincidence symbol')
        return 'True' if a[0]=='True' else 'False'
    if op in ('and','or'): return '('+prop(bridge,a[0])+(' ∧ ' if op=='and' else ' ∨ ')+prop(bridge,a[1])+')'
    if op in ('eq','ne','lt','gt','le','ge'):
        return '('+term(bridge,a[0])+{'eq':' = ','ne':' ≠ ','lt':' < ','gt':' > ','le':' ≤ ','ge':' ≥ '}[op]+term(bridge,a[1])+')'
    if op=='call':
        f,*xs=a
        if f in ('Eq','Ne'):
            bridge.require(len(xs)==2,'wrong equality arity')
            return prop(bridge,['eq' if f=='Eq' else 'ne',*xs])
        if f in ('Q.real','Q.integer','Q.positive','Element'):
            bridge.require(len(xs)==(2 if f=='Element' else 1),'wrong membership arity')
            term(bridge,xs[0])
            sort = xs[1] if f=='Element' else ['symbol','Integers' if f=='Q.integer' else 'Reals']
            if f=='Q.positive': return '(0 < '+term(bridge,xs[0])+')'
            bridge.require(sort in (['symbol','Reals'],['symbol','Integers']),'unsupported membership sort')
            if sort==['symbol','Integers']:
                bridge.require(xs[0] in (['symbol','D'],['symbol','braneDimension']),'integer membership only for fixed D=3')
            return 'True'
        if f=='Inequality':
            bridge.require(len(xs)==5 and xs[1]==xs[3]==['symbol','Less'],'unsupported chained inequality')
            return prop(bridge,['and',['lt',xs[0],xs[2]],['lt',xs[2],xs[4]]])
    raise bridge.BridgeError('unsupported logical tree '+repr(t))


def junction(kind, items):
    if not items: return ['symbol','True' if kind=='and' else 'False']
    result=items[0]
    for item in items[1:]: result=[kind,result,item]
    return result


def assignment(bridge,t):
    if t[0]=='call' and t[1]=='Eq':
        bridge.require(len(t)==4,'wrong coordinate equality')
        bridge.require(t[2] in (['symbol','k1'],['symbol','k2'],['symbol','k3']),'unsupported Python locus assignment target')
        term(bridge,t[2]);term(bridge,t[3])
        return ['eq',t[2],t[3]]
    bridge.require(t[0]=='rule' and t[1][0]=='symbol','expected guarded assignment')
    target,value=t[1:]
    bridge.require(target[1] in ('k1','k2','k3','muR','sRho'),'unsupported locus assignment target')
    if value[0]=='call' and value[1]=='ConditionalExpression':
        bridge.require(len(value)==4,'wrong conditional assignment arity')
        term(bridge,value[2]);prop(bridge,value[3])
        return ['and',['eq',target,value[2]],value[3]]
    term(bridge,value)
    return ['eq',target,value]


def locus(bridge,t):
    if t==['symbol','EmptySet']: return ['symbol','False']
    return junction('or',[junction('and',[assignment(bridge,a) for a in seq(bridge,b)]) for b in seq(bridge,t)])


def node(bridge,t,engine):
    op,*a=t
    if op=='num': return bridge.Node('num',tuple(a))
    if op=='symbol':
        names=bridge.NAMES[engine] if engine else {**bridge.NAMES['PY'],**bridge.NAMES['WL']}
        bridge.require(a[0] in names,'unexpected arithmetic operand symbol')
        return bridge.Node('atom',(names[a[0]],))
    if op=='neg': return bridge.Node('mul',(bridge.Node('num',(-1,)),node(bridge,a[0],engine)))
    bridge.require(op in ('add','sub','mul','div','pow'),'unsupported arithmetic operand')
    if op=='pow':
        bridge.require(a[1][0]=='num' and 0<=a[1][1]<=32,'unbounded arithmetic operand power')
        return bridge.Node('pow',(node(bridge,a[0],engine),a[1][1]))
    return bridge.Node(op,tuple(node(bridge,x,engine) for x in a))


def parse(bridge,raw,engine): return Parser(bridge,raw,engine).parse()

CASES = ('generic','parallel','perpendicular')


def build(bridge,engine,path,case):
    generic=case=='generic'
    pre='' if generic else reruns.prefix(engine,case)
    ns=engine+'CoincidenceArithmetic'+case.title()
    dims=bridge.DIMS if generic else {**bridge.DIMS,'z':(2,-2,0),**{f'k{i}':(0,0,0) for i in range(3)}}
    b=bridge.Builder(engine,path,namespace=ns,support='CoincidenceSupport',
                     tactic='cas_eval' if generic else 'rerun_eval',
                     units='units' if generic else 'coordinateUnits',dims=dims)
    point='k' if generic else engine+'Loci.'+case+'Point'
    pairs=PAIRS[:1] if case=='parallel' else PAIRS
    cells=[(f'referenceRoot rho mu sigma {point} {i} - referenceRoot rho mu sigma {point} {j}',
            (0,-2,0) if generic else (2,-2,0),['rho','sigma']) for i,j in pairs]
    if engine=='PY':
        b.record(pre+'Q3_ROOT_COINCIDENCE_EQUATIONS',lambda x:x,cells,
                 parse_payload=lambda raw,e:[node(bridge,t,e) for t in seq(bridge,parse(bridge,raw,e),len(pairs))])
        b.record(pre+'Q3_ROOT_COINCIDENCE_PAIR_INDICES',lambda x:x,
                 [(str(x+1),(0,0,0),[]) for pair in pairs for x in pair],
                 parse_payload=lambda raw,e:[node(bridge,t,e) for p in seq(bridge,parse(bridge,raw,e),len(pairs)) for t in seq(bridge,p,2)])
    else:
        for pair,cell in zip(pairs,cells):
            i,j=pair
            def operand(raw,e):
                equation,premise=seq(bridge,parse(bridge,raw,e),2)
                bridge.require(equation[0]=='eq' and equation[2]==['num',0],'expected equation to literal zero')
                prop(bridge,premise)
                return [node(bridge,equation[1],e)]
            b.record(pre+f'ROOT{i+1}_ROOT{j+1}_Q3_COINCIDENCE_OPERANDS',lambda x:x,[cell],parse_payload=operand)
    for r in b.records: r['coincidence_arithmetic']={'case':case,'kind':'pair_indices' if 'PAIR_INDICES' in r['tag'] else 'equation'}
    return b


def generate(bridge):
    outputs,engines,audits,records={},{},[],[]
    builders={}
    for engine,path in bridge.INPUTS.items():
        for case in CASES:
            b=build(bridge,engine,path,case);builders[engine,case]=b
            outputs[bridge.GENERATED/(b.namespace+'.lean')]=b.output()
            engines[b.namespace]={'path':str(path.relative_to(bridge.BASE)),'sha256':bridge.sha(path.read_bytes()),
                                  'records':b.records,'unique_trees':len(b.nodes)}
            audits += [b.namespace+'.'+a for a in b.audits]
    lines=['import S10Audit.CAS.CoincidenceSupport',
           *['import S10Audit.CAS.'+b.namespace for b in builders.values()], '',
           'set_option backward.isDefEq.respectTransparency false', '',
           'namespace S10Audit.CAS.Coincidence', 'open S10Pilot S10Anisotropic', 'noncomputable section','']
    params='(rho mu sigma lam : ℝ) (a k : Vec 3)'
    args='rho mu sigma lam a k'
    hypothesis='(hd : coincidenceDomain rho mu sigma lam)'
    def entry(engine,suffix):
        rows=builders[engine,'generic'].rows;tag=bridge.PREFIX+suffix
        bridge.require(tag in rows,'missing coincidence record '+tag)
        line,raw=rows[tag];tree=parse(bridge,raw,engine)
        rec={'tag':engine+'_'+tag,'line':line,'payload':raw,'payload_sha256':bridge.sha(raw.encode()),
             'mixed_ast':tree,'components':[],'lean_semantic_claims':[],
             'coefficient_domain':['0 < rho','0 < mu','0 < sigma','sigma != 1','0 < lambdaScale'],
             'dimension_specialization':3,'coordinates':'arbitrary real k and a; no generic-chart exclusion'}
        records.append(rec)
        return tree,rec
    def emit(rec,label,t,expected,kind='predicate',proof=None):
        name=rec['tag'][0:2]+'_'+label
        lines.extend([f'def {name} {params} : Prop := {prop(bridge,t)}','',
                      f'theorem {name}_correct {params} {hypothesis} :',
                      f'    {name} {args} ↔ ({expected}) := by',
                      proof or ('  rcases hd with ⟨hr, hm, hs, hs1, hl⟩\n'
                      '  have hrn : ¬ rho < 0 := not_lt_of_gt hr\n'
                      '  have hsn : ¬ sigma < 0 := not_lt_of_gt hs\n'
                      '  have hmn : mu ≠ 0 := ne_of_gt hm\n'
                      '  have hsplit := positive_sigma_split sigma hs hs1\n'
                      f'  coincidence_eval {name}'),''])
        claim='Coincidence.'+name+'_correct'
        audits.append(claim);rec['lean_semantic_claims'].append(claim)
        rec['components'].append({'kind':kind,'predicate_ast':t,'lean_predicate':'Coincidence.'+name,
                                  'expected':expected})
        return name
    def allowed(t):
        loc,positive,assumptions=seq(bridge,t,3)
        return junction('and',[locus(bridge,loc),positive,assumptions])
    for engine,path in bridge.INPUTS.items():
        for case in CASES:
            generic=case=='generic';pre='' if generic else reruns.prefix(engine,case)
            pairs=PAIRS[:1] if case=='parallel' else PAIRS
            point='k' if generic else engine+'Loci.'+case+'Point'
            for index,(i,j) in enumerate(pairs):
                # A proof at each actual point also rules out coincidence of distinct stratum roots.
                name=f'{engine}_{case}_pair{index}_roots'
                lines += [f'theorem {name} (rho mu sigma lam : ℝ) (k : Vec 3) {hypothesis} :',
                          f'    referenceRoot rho mu sigma {point} {i} = referenceRoot rho mu sigma {point} {j} ↔',
                          (f'      coincidenceGeometry {index} k' if generic else '      False')+' := by',
                          '  rcases hd with ⟨hr, hm, hs, hs1, _hl⟩',
                          f'  have h := coincidence_pair_geometry rho mu sigma {point} {index} (ne_of_gt hr) (ne_of_gt hm) hs hs1']
                if generic:
                    lines += ['  exact h','']
                else:
                    lines += [f'  norm_num [pairLeft, pairRight, coincidenceGeometry, {point},',
                              '    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, Fin.ext_iff, Fin.coe_ofNat_eq_mod] at h',
                              f'  norm_num [{point}]', '  exact h','']
                audits.append('Coincidence.'+name)
            if engine=='PY':
                for suffix,kind in [('LOCI','locus'),('ALLOWED_OPERANDS','allowed'),('ALLOWED_TESTS','status')]:
                    t,rec=entry(engine,pre+'Q3_ROOT_COINCIDENCE_'+suffix)
                    for index,item in enumerate(seq(bridge,t,len(pairs))):
                        expected=(f'coincidenceGeometry {index} k' if kind=='locus' else
                                  f'allowedCoincidence {index} k' if kind=='allowed' else
                                  f'(∃ w : Vec 3, allowedCoincidence {index} w)') if generic else 'False'
                        pt=locus(bridge,item) if kind=='locus' else allowed(item) if kind=='allowed' else item
                        proof=None
                        if kind=='status' and generic:
                            proof=f'  rw [allowedCoincidence_exists]\n  norm_num [{engine}_{case}_{kind}{index}, Fin.ext_iff, Fin.coe_ofNat_eq_mod]'
                        emit(rec,f'{case}_{kind}{index}',pt,expected,kind,proof)
                # A witness list certifies soundness and existence; it does not cover the entire axis.
                t,rec=entry(engine,pre+'Q3_ROOT_COINCIDENCE_ALLOWED_WITNESSES')
                for index,item in enumerate(seq(bridge,t,len(pairs))):
                    pt=locus(bridge,item)
                    target='(k 0 = 1 ∧ k 1 = 0 ∧ k 2 = 0)' if generic and index==2 else 'False'
                    name=emit(rec,f'{case}_witness{index}',pt,target,'witness')
                    expected=f'allowedCoincidence {index} k' if generic else 'False'
                    lines += [f'theorem {name}_sound {params} {hypothesis}',
                              f'    (hw : {name} {args}) : {expected} := by',
                              f'  rw [{name}_correct rho mu sigma lam a k hd] at hw']
                    if generic and index==2:
                        lines += ['  rw [allowedCoincidence_geometry]', '  rcases hw with ⟨h0, h1, h2⟩', '  simp [h0, h1, h2]','']
                    else: lines += ['  exact hw.elim','']
                    claim='Coincidence.'+name+'_sound';audits.append(claim);rec['lean_semantic_claims'].append(claim)
                    # For the nonempty case, the printed point itself is admissible.
                    if generic and index==2:
                        lines += [f'theorem {name}_exists (rho mu sigma lam : ℝ) (a : Vec 3) {hypothesis} :',
                                  f'    ∃ k : Vec 3, {name} {args} := by',
                                  '  refine ⟨![1, 0, 0], ?_⟩',f'  rw [{name}_correct rho mu sigma lam a _ hd]',
                                  '  norm_num [Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]','']
                        claim='Coincidence.'+name+'_exists';audits.append(claim);rec['lean_semantic_claims'].append(claim)
                for suffix,kind in [('BRANCH_ALLOWED_OPERANDS','branch_allowed'),('BRANCH_ALLOWED_TESTS','branch_status'),
                                    ('BRANCH_ALLOWED_WITNESSES','branch_witness')]:
                    t,rec=entry(engine,pre+'Q3_ROOT_COINCIDENCE_'+suffix)
                    for index,group in enumerate(seq(bridge,t,len(pairs))):
                        branches=seq(bridge,group,1 if generic else 0)
                        if not branches:
                            emit(rec,f'{case}_{kind}{index}',['symbol','False'],'False',kind)
                        for branch in branches:
                            if kind=='branch_allowed':
                                loc,positive,assumptions=seq(bridge,branch,3)
                                pt=junction('and',[junction('and',[assignment(bridge,x) for x in seq(bridge,loc)]),positive,assumptions])
                                expected=f'allowedCoincidence {index} k';proof=None
                            elif kind=='branch_status':
                                pt=branch;expected=f'(∃ w : Vec 3, allowedCoincidence {index} w)'
                                proof=f'  rw [allowedCoincidence_exists]\n  norm_num [{engine}_{case}_{kind}{index}, Fin.ext_iff, Fin.coe_ofNat_eq_mod]'
                            else:
                                pt=junction('and',[assignment(bridge,x) for x in seq(bridge,branch)]) if seq(bridge,branch) else ['symbol','False']
                                expected='k 0 = 1 ∧ k 1 = 0 ∧ k 2 = 0' if index==2 else 'False';proof=None
                            emit(rec,f'{case}_{kind}{index}',pt,expected,kind,proof)
            else:
                for index,(i,j) in enumerate(pairs):
                    pp=pre+f'ROOT{i+1}_ROOT{j+1}_Q3_COINCIDENCE_'
                    t,rec=entry(engine,pp+'OPERANDS')
                    _,premise=seq(bridge,t,2)
                    emit(rec,f'{case}_premise{index}',premise,'k ≠ 0' if generic else 'True','premise')
                    for suffix,kind in [(('LOCUS' if generic else 'PARAMETER_LOCUS'),'locus'),
                                        (('ALLOWED_INTERSECTION' if generic else 'ALLOWED_PARAMETER_REGION'),'allowed'),
                                        (('INTERSECTION_TEST' if generic else 'PARAMETER_INTERSECTION_TEST'),'status'),
                                        (('INTERSECTION_OUTCOME' if generic else 'PARAMETER_INTERSECTION_OUTCOME'),'outcome')]:
                        t,rec=entry(engine,pp+suffix)
                        if kind=='outcome':
                            bridge.require(t in (['symbol','decidedEmpty'],['symbol','decidedNonempty']),'unsupported coincidence outcome')
                            pt=['symbol','True' if t[1]=='decidedNonempty' else 'False']
                        else: pt=locus(bridge,t) if kind=='locus' else t
                        expected=(f'coincidenceGeometry {index} k' if kind=='locus' else
                                  f'allowedCoincidence {index} k' if kind=='allowed' else
                                  f'(∃ w : Vec 3, allowedCoincidence {index} w)') if generic else 'False'
                        proof=None
                        if kind in ('status','outcome') and generic:
                            proof=f'  rw [allowedCoincidence_exists]\n  norm_num [{engine}_{case}_{kind}{index}, Fin.ext_iff, Fin.coe_ofNat_eq_mod]'
                        emit(rec,f'{case}_{kind}{index}',pt,expected,kind,proof)
    # Bind the actual scalar differences and decision records to their loci.
    for (engine,case),b in builders.items():
        generic=case=='generic';pairs=PAIRS[:1] if case=='parallel' else PAIRS
        arithmetic=b.records[:1] if engine=='PY' else b.records
        for index,_ in enumerate(pairs):
            r=arithmetic[0] if engine=='PY' else arithmetic[index]
            c=r['cells'][index if engine=='PY' else 0]
            # Each emitted denominator obligation is discharged from the fixed coefficient domain.
            domain=c['raw_denominator_factors']+[ ['atom',n] for n in c['reference_domain_symbols'] ]
            domain=list(dict.fromkeys(str(d) for d in domain))
            # Use inferred proof arguments, in the same order as Builder.record.
            hs=' '.join('(by positivity)' if d not in ("[\'atom\', \'rho\']","[\'atom\', \'sigma\']") else
                        '(ne_of_gt hr)' if d=="['atom', 'rho']" else '(ne_of_gt hs)' for d in domain)
            name=f'{engine}_{case}_equation{index}_locus'
            lines += [f'theorem {name} (rho mu sigma lam z : ℝ) (a k : Vec 3) {hypothesis} :',
                      f'    Expr.eval (values rho mu sigma z k) {b.namespace}.{c["lean_tree"]} = 0 ↔',
                      f'      {engine}_{case}_locus{index} {args} := by',
                      '  have hr := hd.1', '  have hs := hd.2.2.1',
                      f'  rw [{b.namespace}.{c["lean_claim"]} rho mu sigma z k {hs}, sub_eq_zero,',
                      f'    {engine}_{case}_pair{index}_roots rho mu sigma lam k hd,',
                      f'    {engine}_{case}_locus{index}_correct rho mu sigma lam a k hd]','']
            claim='Coincidence.'+name;audits.append(claim)
            r['coincidence_arithmetic'].setdefault('lean_semantic_claims',[]).append(claim)
            for kind in ('status','outcome') if engine=='WL' else ('status',):
                name=f'{engine}_{case}_{kind}{index}_decides'
                lines += [f'theorem {name} {params} {hypothesis} :',
                          f'    {engine}_{case}_{kind}{index} {args} ↔',
                          f'      ∃ w : Vec 3, {engine}_{case}_allowed{index} rho mu sigma lam a w := by',
                          f'  simp_rw [{engine}_{case}_{kind}{index}_correct rho mu sigma lam a k hd,',
                          f'    {engine}_{case}_allowed{index}_correct rho mu sigma lam a _ hd]']
                if not generic: lines += ['  simp']
                lines += [''];audits.append('Coincidence.'+name)
                rec=next(r for r in records if any(c['lean_predicate']==f'Coincidence.{engine}_{case}_{kind}{index}' for c in r['components']))
                rec['lean_semantic_claims'].append('Coincidence.'+name)
    for index in range(3):
        for kind in ('locus','allowed','status'):
            name=f'generic_{kind}{index}_cross_engine'
            lines += [f'theorem {name} {params} {hypothesis} :',
                      f'    PY_generic_{kind}{index} {args} ↔ WL_generic_{kind}{index} {args} := by',
                      f'  rw [PY_generic_{kind}{index}_correct rho mu sigma lam a k hd,',
                      f'    WL_generic_{kind}{index}_correct rho mu sigma lam a k hd]','']
            audits.append('Coincidence.'+name)
    lines += ['end','end S10Audit.CAS.Coincidence','']
    # Uniform argument names preserve correspondence even when a record ignores a coordinate.
    text=re.sub(r'\b(rho|mu|sigma|lam|a|k|hd)\b',r'_\1','\n'.join(lines))
    outputs[bridge.GENERATED/'CoincidenceBindings.lean']=text
    audits += ['normSq_positive_iff','sumSquares_positive_iff','coincidence_pair_geometry',
               'allowedCoincidence_geometry','allowedCoincidence_exists','positive_sigma_split']
    return outputs,engines,records,audits
