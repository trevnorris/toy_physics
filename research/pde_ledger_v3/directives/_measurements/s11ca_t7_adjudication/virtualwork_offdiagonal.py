#!/usr/bin/env python3
"""S11c-a T7 step-1 VIRTUAL-WORK off-diagonal test: are WL's physical!=virtual DOF pairings nonzero?
Prints the computed SHAPE_DERIVATIVE.EXPRESSION per (physicalDOF,virtualDOF); asserts nothing (rule 2)."""
import re
WL="/var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S11c_a_interface_geometry_mathematica_audit.out"
line=[l for l in open(WL) if l.startswith('WL_S11CA_VIRTUAL_WORK_SHAPE_DERIV:')][0].split(': ',1)[1]
def split_top(s):
    a=[];d=0;i=0;st=0;ins=None
    while i<len(s):
        ch=s[i]
        if ins:
            if ch==ins:ins=None
            i+=1;continue
        if ch in "'\"":ins=ch;i+=1;continue
        two=s[i:i+2]
        if two=='<|':d+=1;i+=2;continue
        if two=='|>':d-=1;i+=2;continue
        if ch in '[{(':d+=1
        elif ch in ']})':d-=1
        elif ch==',' and d==0:a.append(s[st:i]);st=i+1
        i+=1
    a.append(s[st:]);return[x.strip() for x in a if x.strip()]
def entries(v):
    v=v.strip();inner=v[2:-2] if v.startswith('<|') and v.endswith('|>') else v
    return [(m.group(1),m.group(2)) for e in split_top(inner) if (m:=re.match(r'^"([^"]*)"\s*->\s*(.*)$',e.strip(),re.S))]
def field(v,f):
    for k,val in entries(v):
        if k==f:return val
    return None
for k,v in entries(line.strip()):
    pd=re.search(r'\|DOF_(DELTA_W|ZETA_C)\|',k+'|'); vd=re.search(r'VIRTUAL_DOF_(DELTA_W|ZETA_C)',k)
    pd=pd.group(1) if pd else '?';vd=vd.group(1) if vd else '?'
    sd=field(v,'SHAPE_DERIVATIVE'); ex=field(sd,'EXPRESSION') if sd else None
    z = (ex is not None and ex.strip() in ('0','{0}'))
    print(f"{'OFF' if pd!=vd else 'DIAG':4s} phys={pd:8s} virt={vd:8s} zero={z}")
