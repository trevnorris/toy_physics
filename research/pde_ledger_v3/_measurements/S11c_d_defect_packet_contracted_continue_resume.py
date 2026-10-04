"""Inert metadata helpers for an exact saved-prefix continuation.

No scientific imports. Scientific numeric objects are supplied only by the
contained worker. These functions do not call any completed producer or proof.
"""
import ast
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
    require(len(rows)==40 and Counter(v['primitive'] for v in rows)=={'J':10,'Dr':10,'Dh':10,'Dq':10},'forty exact primitive occurrences')
    return rows


def check_census(rows,expected):
    fields=('addressId','primitive','face','grade','K','T','carrier','unit','purpose')
    def key(v):
        require(set(v)>=set(fields),'full census operands')
        return tuple(tuple(v[f]) if isinstance(v[f],list) else v[f] for f in fields)
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
