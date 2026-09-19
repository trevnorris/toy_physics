#!/usr/bin/env python3
"""Resume saved end pairs; retain unsimplified current derivative operands."""
import ast
import copy
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import sympy as sp
import S11c_d_continuum_boundary as c

f=c.f
ORIGINAL=f.STORE/'s11c-continuum-boundary-20260919/production/complete'
PLAN=f.M/'S11c_d_continuum_boundary_recovery_plan.md'
RESTORED={}
CONSUMED={}


def equal(a,b):
    if isinstance(a,np.ndarray):return isinstance(b,np.ndarray) and a.dtype==b.dtype and np.array_equal(a,b)
    if isinstance(a,dict):return isinstance(b,dict) and list(a)==list(b) and all(equal(v,b[k]) for k,v in a.items())
    if isinstance(a,(tuple,list)):return type(a) is type(b) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return bool(a==b)


def check_sources():
    inputs=json.loads((ORIGINAL/'inputs.json').read_text());name=str(Path(c.__file__).resolve().relative_to(f.ROOT))
    for source,sha in inputs['sourceFiles'].items():
        f.require(f.digest(ORIGINAL/'source'/source)==sha,('original frozen source',source))
        if source!=name:f.require(f.digest(f.ROOT/source)==sha,('unchanged original source',source))
    for source,sha in inputs['inputPackets'].items():f.require(f.digest(Path(source))==sha,('original operand',source))
    before=ast.parse((ORIGINAL/'source'/name).read_text());after=ast.parse(Path(c.__file__).read_text())
    old_table=next(n for n in before.body if getattr(n,'name',None)=='current_tables')
    replacements=[]
    class WithoutCancel(ast.NodeTransformer):
        def visit_Call(self,n):
            self.generic_visit(n)
            if isinstance(n.func,ast.Attribute) and n.func.attr=='applyfunc' and len(n.args)==1 and ast.unparse(n.args[0])=='sp.cancel':
                replacements.append(ast.dump(n));return n.func.value
            return n
    WithoutCancel().visit(old_table)
    f.require(len(replacements)==1 and ast.dump(before)==ast.dump(after),'whole-constructor join: remove only final current-table cancel')
    artifacts={p.name:{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in ORIGINAL.glob('*.pickle')}
    f.require(set(artifacts)=={'left-pencil.pickle',*(f'left-mode-{i}.pickle' for i in (1,5,7,16,17))},'exact saved left-end construction inventory')
    return {'originalDirectory':str(ORIGINAL),'originalInputsSha256':f.digest(ORIGINAL/'inputs.json'),
            'removedCalls':replacements,'artifacts':artifacts,'wholeCheckerJoin':True}


def restore_jet(end,*args):
    path=ORIGINAL/(end.lower()+'-pencil.pickle')
    if not path.exists():return c.J(*args)
    saved=f.unpickle(path);pencil,wave,k,q,eta,sigma,origin=args
    f.require(pencil==saved['pencil'] and wave==saved['waves'][0],'saved pencil and wave source identity')
    obj=c.J.__new__(c.J);obj.pencil_coefficients=saved['pencilCoefficients'];obj.radical_coefficients=saved['radicalCoefficients']
    obj.evaluate=sp.lambdify((k,q),tuple(obj.pencil_coefficients.values()),'numpy',cse=True)
    RESTORED[path.name]={'sourceSha256':f.digest(path),'pencilJoin':True,'waveJoin':True};return obj


def restore_pair(end,coefficients,basis):
    matches=[]
    for path in ORIGINAL.glob(end.lower()+'-mode-*.pickle'):
        packet=f.unpickle(path)
        if equal(packet['coefficients'],coefficients) and equal(packet['R'][(0,0)],basis):matches.append((path,packet))
    if not matches:
        f.require(end!='LEFT','every saved left mode must match');return c.J.pair(coefficients,basis)
    f.require(len(matches)==1,'unique complete saved subspace')
    path,packet=matches[0];modes=packet['R'];shift={g:v for g,v in packet['K'].items() if g!=(0,0)}
    computed=c.J.equation(coefficients,modes,shift)
    f.require(equal(computed[(0,0)],packet['diagnostics']['BASE_EQUATION']),'saved full base equation')
    for g in c.J.grades:
        f.require(equal(computed[g],packet['diagnostics']['COEFFICIENT_RESIDUALS'][g]['EQUATION']) and
                  equal(basis.conj().T@modes[g],packet['diagnostics']['COEFFICIENT_RESIDUALS'][g]['GAUGE']),
                  'saved full coefficient/gauge residual replay')
    RESTORED[path.name]={'sourceSha256':f.digest(path),'coefficientsJoin':True,'basisJoin':True,'fullResidualReplay':True}
    return modes,shift,packet['diagnostics'],packet['jacobian']


def write_artifact(path,value):
    old=ORIGINAL/path.name
    if old.exists():
        f.require(equal(value,f.unpickle(old)),'original/copy complete construction packet identity')
        f.require(not path.exists(),'new recovery artifact');shutil.copyfile(old,path)
        f.require(f.digest(path)==f.digest(old),'byte-identical recovered packet')
    else:f.atomic_pickle(path,value)


def save_table_entry(base,end,name,index,value):
    directory=base/'current-tables';directory.mkdir(exist_ok=True)
    path=directory/(end.lower()+'-'+name+'-'+''.join(map(str,index))+'.pickle')
    f.atomic_pickle(path,{'index':index,'value':value})
    CONSUMED[str(path.relative_to(base))]={'bytes':path.stat().st_size,'sha256':f.digest(path)}
    f.save(base/'current-table-inventory.json',CONSUMED)


def derived_constructors():
    source=ast.parse(Path(c.__file__).read_text())
    construct=copy.deepcopy(next(n for n in source.body if getattr(n,'name',None)=='construct_end'))
    table=copy.deepcopy(next(n for n in source.body if getattr(n,'name',None)=='current_tables'))
    original_construct=ast.dump(construct);original_table=ast.dump(table)
    calls=[]
    class Hooks(ast.NodeTransformer):
        def visit_Call(self,n):
            self.generic_visit(n);key=ast.unparse(n.func)
            mapping={'J':'restore_jet','J.pair':'restore_pair','f.atomic_pickle':'write_artifact','current_tables':'tracked_tables'}
            if key in mapping:
                calls.append(key);n.func=ast.Name(mapping[key],ast.Load())
                if key in ('J','J.pair'):n.args.insert(0,ast.Name('end',ast.Load()))
                if key=='current_tables':n.args=[ast.Name(v,ast.Load()) for v in ('base','end','name')]+n.args
            return n
    Hooks().visit(construct)
    f.require(sorted(calls)==sorted(['J','J.pair','f.atomic_pickle','f.atomic_pickle','f.atomic_pickle','current_tables']),
              'only constructor cache/artifact/current checkpoint hooks')
    class Undo(ast.NodeTransformer):
        def visit_Call(self,n):
            self.generic_visit(n);key=ast.unparse(n.func)
            mapping={'restore_jet':'J','restore_pair':'J.pair','write_artifact':'f.atomic_pickle','tracked_tables':'current_tables'}
            if key in mapping:
                n.func=ast.parse(mapping[key],mode='eval').body
                if key in ('restore_jet','restore_pair'):n.args=n.args[1:]
                if key=='tracked_tables':n.args=n.args[3:]
            return n
    f.require(ast.dump(Undo().visit(copy.deepcopy(construct)))==original_construct,'whole numerical constructor reverse AST join')
    table.name='tracked_tables';table.args.args=[ast.arg(v) for v in ('base','end','name')]+table.args.args
    additions=[]
    class TableCheckpoint(ast.NodeTransformer):
        def visit_Assign(self,n):
            if any(isinstance(t,ast.Subscript) and isinstance(t.value,ast.Name) and t.value.id=='table' for t in n.targets):
                additions.append(ast.dump(n));index=copy.deepcopy(n.targets[0].slice)
                call=ast.Expr(ast.Call(ast.Name('save_table_entry',ast.Load()),[ast.Name(v,ast.Load()) for v in ('base','end','name')]+[index,ast.Subscript(ast.Name('table',ast.Load()),copy.deepcopy(index),ast.Load())],[]))
                return [n,call]
            return n
    TableCheckpoint().visit(table);f.require(len(additions)==1,'one per-entry durable table checkpoint')
    stripped=copy.deepcopy(table);stripped.name='current_tables';stripped.args.args=stripped.args.args[3:]
    class RemoveHook(ast.NodeTransformer):
        def visit_Expr(self,n):
            if isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='save_table_entry':return None
            return n
    RemoveHook().visit(stripped);f.require(ast.dump(stripped)==original_table,'whole derivative constructor reverse AST join')
    namespace=dict(vars(c));namespace.update(restore_jet=restore_jet,restore_pair=restore_pair,
        write_artifact=write_artifact,save_table_entry=save_table_entry)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[table,construct],type_ignores=[])),str(Path(__file__)), 'exec'),namespace)
    return namespace['construct_end'],{'constructEndReverseAstJoin':True,'currentTablesReverseAstJoin':True,'hooks':calls}


def main():
    original_check=check_sources();construct,joins=derived_constructors();old_load=c.load
    def load(base):
        values=old_load(base);pins,operands=values[-2:]
        for p in (Path(__file__),PLAN):
            name=str(p.resolve().relative_to(f.ROOT));pins[name]=f.digest(p)
            dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
        for name,item in original_check['artifacts'].items():operands[str(ORIGINAL/name)]=item['sha256']
        operands[str(ORIGINAL/'inputs.json')]=original_check['originalInputsSha256']
        inputs=json.loads((base/'inputs.json').read_text());inputs.update(sourceFiles=pins,inputPackets=operands,
            savedOperandRecovery=original_check,instrumentJoins=joins);f.save(base/'inputs.json',inputs)
        return values
    c.load=load;c.construct_end=construct
    c.main()
    # The original main's complete validation has finished; retain the explicit
    # cache/derivative inventories without changing its checks/stdout contract.
    import sys
    base=Path(sys.argv[sys.argv.index('--run-directory')+1])
    f.save(base/'recovery.json',{'original':original_check,'joins':joins,'restored':RESTORED,'currentTables':CONSUMED})
    f.require(len(RESTORED)==6,'all saved left pencil and complete modes resumed')


if __name__=='__main__':main()
