import json, sympy as sp
from sympy.core.symbol import Str
IDX='/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_c2_sympy_stdout_index.json'
OUT='/tmp/S11c_c2_selfenergy_fold_sympy_audit.out'
_idx=json.load(open(IDX)); _tags=_idx['tags']

def _slice(tag):
    e=_tags[tag]
    with open(OUT,'rb') as f:
        f.seek(e['byte_offset']); return f.read(e['byte_count']).decode('utf-8','replace')

def _balanced(s,start):
    depth=0
    for i in range(start,len(s)):
        if s[i]=='(':depth+=1
        elif s[i]==')':
            depth-=1
            if depth==0:return s[start:i+1]
    return s[start:]

def _field_src(payload,name):
    key="Str('%s')"%name
    j=payload.find(key)
    if j<0:return None
    k=payload.find(',',j+len(key)); m=k+1
    while payload[m]==' ':m+=1
    if payload[m]=='(':return _balanced(payload,m)
    depth=0
    for i in range(m,len(payload)):
        if payload[i]=='(':depth+=1
        elif payload[i]==')':
            if depth==0:return payload[m:i]
            depth-=1

_ns={'Str':Str, **vars(sp)}
from sympy.functions.elementary.piecewise import ExprCondPair
_ns['ExprCondPair']=ExprCondPair
def restore(src):
    return eval(src,{'__builtins__':{}}, _ns)

def value_of(tag):
    p=_slice(tag)
    # strip leading "TAG: " prefix
    j=p.find(': '); p=p[j+2:]
    src=_field_src(p,'VALUE')
    return restore(src)

def to_dict(t):
    # nested Tuple of (Str(name), value) -> dict, recursively
    if isinstance(t,sp.Tuple):
        items=list(t)
        if items and all(isinstance(x,sp.Tuple) and len(x)==2 and isinstance(x[0],Str) for x in items):
            return {str(k):to_dict(v) for k,v in items}
        return [to_dict(x) for x in items]
    return t

if __name__=='__main__':
    import sys
    v=value_of(sys.argv[1])
    d=to_dict(v)
    print(type(d), list(d.keys()) if isinstance(d,dict) else 'list/expr')
