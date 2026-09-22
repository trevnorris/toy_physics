#!/usr/bin/env python3
"""Continue unstarted material source calculations across a JSON settings boundary."""
import ast
import copy
import json
from pathlib import Path
import S11c_d_remaining_case_frequency_material_sources as original

NAME='S11c_d_remaining_case_frequency_material_sources_recover'
OLD_SHA='7a9b4804868c38992d8b5895a242a38f440b4dc1fe49483346d3c5036ddc4e2b'
ROOT=original.F/'material-sources-recovery-01'
INVENTORY_SHA='4f89011c80748828e73f9c5301f2174a44d8a47b9cc79d3267c75b65590e0d0c'


def reuse_metadata(reader,journal,name):
    inventory=reader.json(ROOT/'failed-file-inventory.json',INVENTORY_SHA)
    old=original.F/'material-sources/complete'/name;record=inventory[str(old)]
    r=reader.retain(old,record['sha256']);target=journal.base/name
    original.require(not target.exists() and not target.is_symlink(),'fresh completed metadata reference')
    target.symlink_to(r['canonical'])
    journal.artifacts[name]={'path':str(target),'sha256':r['sha256'],'bytes':r['bytes']}
    reader.retain(target,r['sha256'])


def settings_boundary(journal,folder,common,context,raw,view):
    # Native settings keep their types. JSON is a separately labelled metadata view.
    flags={'contextPair':original.same(*common['contextPair']),
           'basisPair':original.same(*common['basisPair']),
           'nativeSettings':original.same(context['settings'],raw['settings']),
           'metadataView':original.same(json.loads(json.dumps(raw['settings'])),view['settings'])}
    journal.json(folder+'/native-settings-metadata-boundary.json',{
        'nativeSource':view['nativeSource'],'physicalContext':view['physicalRoutes']['context'],
        'nativeSettingsTypes':{k:type(v).__name__ for k,v in raw['settings'].items()},
        'metadataSettingsTypes':{k:type(v).__name__ for k,v in view['settings'].items()},
        'checks':flags,'physicalSettingsCoerced':False})
    return all(flags.values())


def provenance(reader,journal):
    inventory=reader.json(ROOT/'failed-file-inventory.json',INVENTORY_SHA)
    for path,r in inventory.items():reader.retain(path,r['sha256'])
    journal.json('material-source-settings-recovery.json',{
        'originalHelper':reader.retain(original.__file__,OLD_SHA),'recoveryHelper':reader.retain(__file__),
        'failedInventory':reader.retain(ROOT/'failed-file-inventory.json',INVENTORY_SHA),
        'failureInspection':reader.retain(ROOT/'failure-inspection.json'),
        'wholeMainReverseAST':True,'completedMetadataReferences':2,'originalNumericalCalls':0,
        'change':'Strict raw context/basis/settings equality plus separate native JSON metadata view; all computational settings retain native types. Numerical method and controls unchanged.'})


def adapted(text):
    fn=next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name=='main');before=copy.deepcopy(fn)
    replacements=[]
    for parent in ast.walk(fn):
        for field,value in ast.iter_fields(parent):
            children=list(enumerate(value)) if isinstance(value,list) else [(None,value)]
            for i,node in children:
                replacement=None
                if isinstance(node,ast.Expr) and isinstance(node.value,ast.Call):
                    call=node.value
                    if ast.unparse(call.func)=='journal.json' and isinstance(call.args[0],ast.Constant) and call.args[0].value in ('source-preparation-scope.json','native-and-numerical-method-joins.json'):
                        replacement=ast.parse('reuse_metadata(reader,journal,'+repr(call.args[0].value)+')').body[0]
                    elif ast.unparse(call.func)=='require' and len(call.args)>1 and isinstance(call.args[1],ast.Constant) and call.args[1].value=='full own context/unit basis/settings':
                        replacement=ast.parse("require(settings_boundary(journal,folder,common,context,raw,view),'full own context/unit basis/settings')").body[0]
                elif isinstance(node,ast.Subscript) and ast.unparse(node)=="view['settings']":
                    replacement=ast.parse("raw['settings']",mode='eval').body
                if replacement is not None:
                    replacements.append((parent,field,i,copy.deepcopy(node)))
                    if i is None:setattr(parent,field,replacement)
                    else:value[i]=replacement
    # The whole require is replaced; mutations of its detached children are irrelevant.
    reverse=copy.deepcopy(fn)
    for parent,field,i,old in reversed(replacements):
        if i is None:setattr(parent,field,old)
        else:getattr(parent,field)[i]=old
    original.require(ast.dump(fn)==ast.dump(before),'entire original material source main reverses')
    changed=reverse
    original.require(sum(isinstance(n,ast.Call) and ast.unparse(n.func)=='reuse_metadata' for n in ast.walk(changed))==2,'two completed metadata references')
    original.require(sum(isinstance(n,ast.Call) and ast.unparse(n.func)=='settings_boundary' for n in ast.walk(changed))==1,'one failed metadata boundary')
    original.require(not any(isinstance(n,ast.Subscript) and ast.unparse(n)=="view['settings']" for n in ast.walk(changed)),'computational settings stay native')
    at=next(i for i,n in enumerate(changed.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='reader.postcheck')
    changed.body.insert(at,ast.parse('provenance(reader,journal)').body[0])
    return ast.fix_missing_locations(ast.Module(body=[changed],type_ignores=[]))


def main():
    original.require(original.saved.digest(original.__file__)==OLD_SHA,'immutable failed material source helper')
    env=dict(vars(original),NAME=NAME,__file__=str(Path(__file__).resolve()),reuse_metadata=reuse_metadata,settings_boundary=settings_boundary,provenance=provenance)
    exec(compile(adapted(Path(original.__file__).read_text()),'<native settings material source continuation>','exec'),env);env['main']()


if __name__=='__main__':main()
