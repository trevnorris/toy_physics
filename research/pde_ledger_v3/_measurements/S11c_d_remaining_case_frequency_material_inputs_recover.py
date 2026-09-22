#!/usr/bin/env python3
"""Continue material metadata with exact canonical source-value receipt identity."""
import ast
import copy
from pathlib import Path
import S11c_d_remaining_case_frequency_material_inputs as original
NAME='S11c_d_remaining_case_frequency_material_inputs_recover'
OLD_SHA='b6a3e4411142822b037b436f4a8c6e381862593c5225062a1c8e2aab36843e20'
ROOT=original.F/'material-inputs-recovery-01'
INVENTORY_SHA='2235a35a6c54f365120cb58445ac39544c9e363208d6943bb5e49ba669641f19'


def reuse_native_caller(reader,journal):
    inventory=reader.json(ROOT/'failed-file-inventory.json',INVENTORY_SHA)
    path=original.F/'material-inputs/complete/native-callers-and-prior-preparation.json';r=inventory[str(path)]
    rec=reader.retain(path,r['sha256']);target=journal.base/path.name
    original.require(not target.exists() and not target.is_symlink(),'fresh completed metadata reference')
    target.symlink_to(rec['canonical']);journal.artifacts[path.name]={'path':str(target),'sha256':rec['sha256'],'bytes':rec['bytes']};reader.retain(target,rec['sha256'])


def canonical_value_equal(reader,journal,a,b):
    left=reader.retain(a['path'],a['sha256']);right=reader.retain(b['path'],b['sha256'])
    matched=all(left[k]==right[k] for k in ('canonical','sha256','bytes')) and a['bytes']==b['bytes']==left['bytes']
    journal.json('source22-original-and-accepted-value-routes.json',{'originalReceipt':a,'acceptedManifest':b,'originalRoute':left,'acceptedRoute':right,'sameCanonicalBytes':matched,
        'physicalValuesChanged':False,'readOrReconstructSourceValue':False})
    return matched


def provenance(reader,journal):
    inventory=reader.json(ROOT/'failed-file-inventory.json',INVENTORY_SHA)
    for path,rec in inventory.items():reader.retain(path,rec['sha256'])
    journal.json('material-reader-recovery-join.json',{'originalHelper':reader.retain(original.__file__,OLD_SHA),'newHelper':reader.retain(__file__),
        'failedInventory':reader.retain(ROOT/'failed-file-inventory.json',INVENTORY_SHA),'failureInspection':reader.retain(ROOT/'failure-inspection.json'),
        'wholeMainReverseAST':True,'completedMetadataReused':1,'physicalWriterOrScientificCalls':0,'repair':'Only source22 original/linked accepted value addresses compare canonical file/hash/size; both raw/logical/link routes retained and postchecked.'})


def adapted(text):
    fn=next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name=='main');before=copy.deepcopy(fn)
    changes=[]
    for i,node in enumerate(fn.body):
        if not isinstance(node,ast.Expr) or not isinstance(node.value,ast.Call):continue
        call=node.value
        if ast.unparse(call.func)=='journal.json' and isinstance(call.args[0],ast.Constant) and call.args[0].value=='native-callers-and-prior-preparation.json':
            changes.append((i,copy.deepcopy(node)));fn.body[i]=ast.parse('reuse_native_caller(reader,journal)').body[0]
        elif ast.unparse(call.func)=='require' and len(call.args)>1 and isinstance(call.args[1],ast.Constant) and call.args[1].value=='actual completed source22 recurrence':
            changes.append((i,copy.deepcopy(node)));call.args[0]=ast.parse("canonical_value_equal(reader,journal,fr['value'],fa['source-action/value.pickle'])",mode='eval').body
    original.require(len(changes)==2,'exact one completed metadata route and failed source22 address guard')
    at=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='reader.postcheck')
    fn.body.insert(at,ast.parse('provenance(reader,journal)').body[0])
    reverse=copy.deepcopy(fn);reverse.body.pop(at)
    for i,node in changes:reverse.body[i]=node
    original.require(ast.dump(reverse)==ast.dump(before),'entire original reader reverse AST')
    return ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))


def main():
    original.require(original.saved.digest(original.__file__)==OLD_SHA,'immutable original material reader')
    env=dict(vars(original),NAME=NAME,__file__=str(Path(__file__).resolve()),reuse_native_caller=reuse_native_caller,canonical_value_equal=canonical_value_equal,provenance=provenance)
    exec(compile(adapted(Path(original.__file__).read_text()),'<canonical saved-source metadata recovery>','exec'),env);env['main']()


if __name__=='__main__':main()
