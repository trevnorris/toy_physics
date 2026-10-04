"""Inert metadata helpers for an exact saved-prefix continuation.

No scientific imports. Scientific numeric objects are supplied only by the
contained worker. These functions do not call any completed producer or proof.
"""
import ast,json
from collections import Counter
from pathlib import Path

ZERO={'text':'0','srepr':'Integer(0)'}
UNITS=['-2','-1','1']
COMPONENTS=('NATIVE_MIXED_ITERATION','INHERITED_DIRECT_WHOLE_OFF_DIAGONAL')


def require(v,message):
    if v is not True:raise ValueError(message)


def canonical_slots(entries,K,T,carrier):
    rows=[]
    for e in entries:
        expected=['J'] if e['component']==COMPONENTS[0] else ['Dr','Dh','Dq']
        require(e['component'] in COMPONENTS and e['primitives']==expected,'native primitive membership')
        require(e['face'] in ('plus','minus') and e['grade']==[1,1] and e['unit']==UNITS,'native face grade unit')
        for p in expected:
            rows.append({'addressId':e['addressId'],'primitive':p,'face':e['face'],'grade':e['grade'],
                         'K':K,'T':T,'carrier':carrier,'unit':e['unit'],'purpose':'baseline-unmutated'})
    require(len({e['addressId'] for e in entries})==len(entries)==20,'twenty unique live native addresses')
    require(Counter(e['face'] for e in entries)=={'plus':10,'minus':10},'both native faces in actual rows')
    require(len(rows)==40 and Counter(v['primitive'] for v in rows)=={'J':10,'Dr':10,'Dh':10,'Dq':10},'forty exact primitive occurrences')
    return rows


def check_census(rows,expected):
    fields=('addressId','primitive','face','grade','K','T','carrier','unit','purpose','boundSource','outerOperand','middleOperand')
    def key(v):
        require(set(v)>=set(fields),'full census operands')
        return tuple(json.dumps(v[f],sort_keys=True,separators=(',',':')) for f in fields)
    require(Counter(map(key,rows))==Counter(map(key,expected)),'exact native primitive/window/face census')
    require(len(rows)==len(expected)==40,'no missing or duplicate primitive')
    return True


def group_indices(rows):
    groups={'full':list(range(len(rows)))}
    for i,r in enumerate(rows):
        for key in ['face/'+r['face'],'grade/'+','.join(map(str,r['grade'])),
                    'component/'+('J' if r['primitive']=='J' else 'D'),
                    'primitive/'+r['primitive'],'address/'+str(r['addressId'])]:
            groups.setdefault(key,[]).append(i)
    return groups


def exact_constant(sp,C,record):
    # Materialization only inside containment; no parser extension or evaluation
    # of a source string. Tail scalars must be actual exact Rational instances.
    value=C.restore_scalar(sp,record)
    require(value.is_Rational is True and value.is_finite is True,'finite exact rational tail operand')
    return value


def guard_inventory(text,stop_line=None):
    tree=ast.parse(text);out=[]
    for node in ast.walk(tree):
        if not (isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='require'):continue
        owners=[f for f in ast.walk(tree) if isinstance(f,ast.FunctionDef) and f.lineno<=node.lineno<=f.end_lineno]
        owner=min(owners,key=lambda n:n.end_lineno-n.lineno).name
        position='conditional executed instances attested from original control flow'
        if stop_line is not None and owner=='prepare':
            position='failed' if node.lineno==stop_line else ('unreached' if node.lineno>stop_line else position)
        out.append({'function':owner,'line':node.lineno,'endLine':node.end_lineno,'source':ast.get_source_segment(text,node),
                    'ast':ast.dump(node,include_attributes=False),'disposition':position,
                    'originalPredicateFunctionReexecuted':False})
    return sorted(out,key=lambda x:x['line'])


def prefix_return(name,record):
    require(set(record)>= {'input','return'},'full original input and return')
    if name.startswith('new-enlarged-tail-'):
        decision=record['decision'];ret=record['return']
        require(ret=={k:decision[k] for k in ('outer','middle')} and decision['originalDomain'] is True and decision['sourceFunctionCalled'] is False,'saved enlarged tail return')
        require(record['input']['K']==29 and record['input']['T']==124,'original enlarged arguments')
        return 'enlarged-tail-observation'
    require(record['return']=={'residual':ZERO} and record['decision']=={'cancelled':ZERO},'literal exact saved zero return')
    require(record['input']['newDerivation'] is True and set(record['input'])>={'left','right','context'},'original complete identity arguments')
    require(set(record['raw'])=={'residual'},'original raw residual retained')
    return 'inherited-new-zero'


def tree_inventory(root):
    root=Path(root)
    return {str(p.relative_to(root)):{'bytes':p.stat().st_size,'symlink':p.is_symlink()}
            for p in sorted(root.rglob('*')) if p.is_file() or p.is_symlink()}


def check_tree(observed,expected,count=8271,total=158628313):
    require(set(observed)==set(expected),'exact prior file tree')
    require(len(observed)==len(expected)==count,'exact prior file count')
    require(sum(v['bytes'] for v in observed.values())==sum(v['bytes'] for v in expected.values())==total,'exact prior byte total')
    require(all(not observed[n]['symlink'] and observed[n]['bytes']==v['bytes'] for n,v in expected.items()),'regular complete prior files')
    return True


def join_base_census(baseplan,addresses,entries):
    native={a['addressId']:a for a in addresses}
    base=baseplan['allAddresses'];ids=[v['addressId'] for v in base]
    require(len(native)==len(addresses)==544 and len(ids)==len(set(ids)),'unique complete native/base ids')
    require(set(ids)=={i for i,a in native.items() if not a['status'].startswith('EXACT_ZERO')},'base list is complete native live row')
    require(all(v['component']==native[v['addressId']]['component'] and v['grade']==native[v['addressId']]['targetGrade'] for v in base),'base/native component and grade')
    selected=[v for v in base if v['component'] in COMPONENTS]
    require({v['addressId'] for v in selected}=={e['addressId'] for e in entries} and len(selected)==len(entries)==20,'independent saved preflight/live entry census')
    require(Counter(native[v['addressId']]['face'] for v in selected)=={'plus':10,'minus':10},'both native faces ten addresses each')
    require(Counter((native[v['addressId']]['face'],v['component']) for v in selected)=={(f,c):5 for f in ('plus','minus') for c in COMPONENTS},'five J and five direct per face')
    return selected


def bound_source(ident,K,T,value):
    return {'addressId':ident,'K':K,'T':T,'alias':'preflight/tail-plan.json/allAddresses/'+str(ident) if K==27 else 'prior/complete/new-enlarged-tail-'+str(ident)+'-return.json',
            'outer':value['outer'],'middle':value['middle']}


def independent_slots(baseplan,enlarged,addresses,K,T,carrier):
    # Expected identities and operands come from the original full preflight
    # list/native rows and full published returns, not the continuation entries.
    require((baseplan['K'],baseplan['T'])==(27,122) and (K,T) in ((27,122),(29,124)),'actual saved/selected windows')
    native={a['addressId']:a for a in addresses};result=[]
    for v in baseplan['allAddresses']:
        if v['component'] not in COMPONENTS:continue
        ident=v['addressId'];a=native[ident]
        require(not a['status'].startswith('EXACT_ZERO'),'base member is actually live')
        if K==27:value=v
        else:
            pair=enlarged[ident];require(pair['input']['K']==K and pair['input']['T']==T and pair['input']['address']['addressId']==ident,'actual enlarged source window/address');value=pair['return']
        for p in (['J'] if v['component']==COMPONENTS[0] else ['Dr','Dh','Dq']):
            result.append({'addressId':ident,'primitive':p,'face':a['face'],'grade':v['grade'],'K':K,'T':T,'carrier':carrier,'unit':UNITS,'purpose':'baseline-unmutated',
                           'boundSource':bound_source(ident,K,T,value),'outerOperand':value['outer'],'middleOperand':value['middle']})
    require(len(result)==40,'independent forty primitive values')
    return result


def tail_formula_join(executed_text,base_text):
    # Constructor/source AST only. No formula or completed callable evaluation.
    executed=ast.parse(executed_text);base=ast.parse(base_text)
    funcs=[v for v in ast.walk(base) if isinstance(v,ast.FunctionDef) and v.name=='contributions'];require(len(funcs)==1,'one original contributions source')
    fn=funcs[0]
    def branch_eq(node,left,right):return isinstance(node,ast.If) and ast.dump(node.test)==ast.dump(ast.parse(left+' == '+right,mode='eval').body)
    window=[n for n in ast.walk(executed) if branch_eq(n,'K','27')];require(len(window)==1,'one executed enlarged branch')
    enlarged=window[0].orelse
    selection=[n for n in enlarged if branch_eq(n,"e['component']",'COMPONENTS[0]')];require(len(selection)==1,'actual enlarged J/direct selection')
    baseJ=[n for n in ast.walk(fn) if branch_eq(n,'comp',repr(COMPONENTS[0]))];require(len(baseJ)==1,'actual base J selection')
    require(len(baseJ[0].orelse)==1 and branch_eq(baseJ[0].orelse[0],'comp',repr(COMPONENTS[1])),'actual base direct selection')
    baseOuter=[n for n in ast.walk(fn) if isinstance(n,ast.If) and branch_eq(n,'comp',"'NATIVE_FLAT'")];require(len(baseOuter)==1 and len(baseOuter[0].orelse)==1,'actual base ordinary outer branch')
    outer=[n for n in enlarged if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id=='outerbound'];require(len(outer)==1,'actual enlarged outer assignment')
    class Names(ast.NodeTransformer):
        def visit_Name(self,node):return ast.copy_location(ast.Name(id={'outerbound':'outer','middlebound':'middle'}.get(node.id,node.id),ctx=node.ctx),node)
    import copy
    pairs={'outer':(baseOuter[0].orelse[0],outer[0]),'J':(baseJ[0].body[0],selection[0].body[0]),'H-overcount':(baseJ[0].body[1],selection[0].body[1]),'D':(baseJ[0].orelse[0].body[0],selection[0].orelse[0])}
    require(len(baseJ[0].body)==len(selection[0].body)==2 and len(selection[0].orelse)==1,'complete J/H/D assignment census')
    return {k:{'baseSource':ast.get_source_segment(base_text,b),'executedSource':ast.get_source_segment(executed_text,e),
               'baseAST':ast.dump(b),'executedMappedAST':ast.dump(Names().visit(copy.deepcopy(e))),
               'same':ast.dump(b)==ast.dump(Names().visit(copy.deepcopy(e))),
               'renamesOnly':{'outerbound':'outer','middlebound':'middle'},'KandTlim':'same symbolic arguments; actual saved windows are joined separately'} for k,(b,e) in pairs.items()}


def checked_totals(rows,expected,materialize,encode):
    check_census(rows,expected)
    totals=[]
    for row in rows:
        # The encoded operands have already joined the independent saved census.
        outer=materialize(row['outerOperand']);middle=materialize(row['middleOperand'])
        require(encode(row['outer'])==row['outerOperand'] and encode(row['middle'])==row['middleOperand'],'live outer/middle objects equal validated operands')
        total=outer+middle
        require(encode(row['total'])==encode(total),'live total equals validated operand sum')
        totals.append(total)
    return totals


def continuation_source_joins(worker,original_worker,tail,original_prepare):
    def fn(text,name):
        nodes=[v for v in ast.parse(text).body if isinstance(v,ast.FunctionDef) and v.name==name]
        require(len(nodes)==1,'one source function '+name);return nodes[0]
    newrun=fn(worker,'run');oldrun=fn(original_worker,'run');old=fn(original_prepare,'prepare');new=fn(tail,'continue_geometry')
    starts=[i for i,n in enumerate(old.body) if isinstance(n,ast.Assign) and any(isinstance(v,ast.Name) and v.id=='plans' for v in n.targets)]
    require(len(starts)==1,'one original unfinished geometry start')
    oldtail=ast.Module(body=old.body[starts[0]:],type_ignores=[]);newtail=ast.Module(body=new.body[2:],type_ignores=[])
    return {'run':{'original':ast.unparse(oldrun),'new':ast.unparse(newrun),'same':ast.dump(newrun)==ast.dump(oldrun)},
            'geometryTail':{'original':ast.unparse(oldtail),'new':ast.unparse(newtail),'same':ast.dump(oldtail)==ast.dump(newtail)}}


def tail_binding_joins(base_text,executed_text):
    base=ast.parse(base_text);exe=ast.parse(executed_text);result={}
    def bind(tree,target):
        nodes=[n for n in ast.walk(tree) if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id==target]
        require(len(nodes)==1,'one actual binding '+target);return nodes[0].value
    expected_base={'b':'sp.Rational(3000,11101)','E15':'sp.Integer(3)**15','E30':'E15**2','F':'{n:weighted_tail(n,0) for n in range(4)}',
                   'KJ':'sp.Rational(4,5)*121*Cq/b**3','KD':'18*121*2*Cq/b**2','Cq':'sp.Integer(36)','Cordinary':'1+4*KJ+4*KD+400/b**2'}
    expected_new={'tail_context':"get('preflight/tail-bound-derivation.json')",'bstar':"R(tail_context['b'])",'ordinary':"R(tail_context['ordinaryKernelUpper'])",'F3':"R(tail_context['fullExponentialMoments']['3'])",'b':'bstar','Cordinary':'ordinary','E30':'sp.Integer(3)**30','F':"{2:R(tail_context['fullExponentialMoments']['2']),3:F3}"}
    for label,tree,specs in [('base',base,expected_base),('enlarged',exe,expected_new)]:
        for name,spec in specs.items():
            actual=bind(tree,name);ref=ast.parse(spec,mode='eval').body
            result[label+'/'+name]={'actualSource':ast.unparse(actual),'requiredSource':spec,'same':ast.dump(actual)==ast.dump(ref)}
    for name in ('exponential_moment','weighted_tail'):
        old=[n for n in ast.walk(base) if isinstance(n,ast.FunctionDef) and n.name==name];new=[n for n in ast.walk(exe) if isinstance(n,ast.FunctionDef) and n.name==name]
        require(len(old)==len(new)==1,'one actual tail helper '+name)
        result['helper/'+name]={'baseSource':ast.unparse(old[0]),'enlargedSource':ast.unparse(new[0]),'same':ast.dump(old[0])==ast.dump(new[0])}
    emits=[n for n in ast.walk(base) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='emit' and n.args and isinstance(n.args[0],ast.Constant) and n.args[0].value=='tail-bound-derivation']
    require(len(emits)==1 and isinstance(emits[0].args[1],ast.Dict),'actual base saved binding output')
    payload={k.value:v for k,v in zip(emits[0].args[1].keys,emits[0].args[1].values)}
    for key,expression in [('b','b'),('ordinaryKernelUpper','Cordinary'),('fullExponentialMoments','{str(k):v for k,v in F.items()}')]:
        result['base-output/'+key]={'actualSource':ast.unparse(payload[key]),'requiredSource':expression,'same':ast.dump(payload[key])==ast.dump(ast.parse(expression,mode='eval').body)}
    # Exact constructor pattern plus integer-power law; no saved bound or
    # large constant is reevaluated. The individual source patterns above gate it.
    result['E30-identity']={'base':'(3^15)^2','enlarged':'3^30','law':'Integer power composition on positive base3; 15*2=30. Formal source-binding identity only.',
                           'same':all(result[k]['same'] for k in ('base/E15','base/E30','enlarged/E30'))}
    return result
