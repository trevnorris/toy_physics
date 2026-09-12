#!/usr/bin/env python3
"""Classify decoded baseline differences without weakening algebraic equality.

Only a globally injective renaming of printed Dummy identities, documented
unordered mappings/closure edges, and explicitly identified run-provenance
records receive special treatment. Existing native values use their separate
field-preservation inventory; every other difference remains unclassified.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from ledger_fold import _restore
from S11c_d_output_codec import decoded_lines
from S11c_d_joint_sheet_inventory import association,binding_order_comparison,sp

DUMMY=re.compile(r"Dummy\(([^()]*)dummy_index=(-?[0-9]+)\)")


def inspect(raw,native):
    before,after=Path(raw['baseline']),Path(raw['transcript'])
    if native['sha256']!=raw['sha256'] or native['nativePreservation']['baselineSha256']!=raw['baselineSha256']:
        raise ValueError('native inventory transcript mismatch')
    changed=set(raw['changedPayloads']);unchanged_dummy_ids=set()
    def collect(path):
        result={}
        for line in decoded_lines(path):
            tag,_,body=line.rstrip('\n').partition(': ')
            if tag in changed:result[tag]=body
            else:unchanged_dummy_ids.update(match[2] for match in DUMMY.finditer(body))
        return result
    a,b=collect(before),collect(after);classes=defaultdict(list)
    preserved=set(native['nativePreservation']['associationOrderOrAddedFieldOnlyTags'])
    for name in sorted(changed & preserved):classes['nativeAddedFieldsOrAssociationOrder'].append(name)
    operations={'PY_S11CD_'+v for v in ('BUILD_INPUT_DIGESTS','IMPLEMENTATION_CHECKPOINT','RESOURCE_MEASUREMENTS','EMISSION_LINES','PROCESS_COMPLETION')}
    for name in sorted(changed & operations):classes['producerProvenanceOrCheckpoint'].append(name)
    accounted=preserved|operations
    names=[k.removeprefix('PY_S11CD_')for k in sorted(changed-accounted)if k.endswith('_BOUND_CARRIERS')]
    for k in binding_order_comparison(before,after,names):
        name='PY_S11CD_'+k;classes['boundCarrierAssociationOrder'].append(name);accounted.add(name)
    for name in sorted(changed-accounted):
        if name=='PY_S11CD_IMPORT_CLOSURE':
            left,right=association(_restore(a[name])),association(_restore(b[name]))
            wanted={'import_keys','closure','symbol_edges','dimension_edges'}
            if set(left)==set(right)==wanted and all(Counter(map(sp.srepr,left[k]))==Counter(map(sp.srepr,right[k]))for k in wanted):
                classes['closureSetAndEdgeOrder'].append(name);accounted.add(name)
        elif name=='PY_S11CD_INFERRED_INPUT_DIMENSIONS':
            left,right=association(_restore(a[name])),association(_restore(b[name]))
            if left==right:
                classes['inferredDimensionAssociationOrder'].append(name);accounted.add(name)
        elif name.startswith('PY_S11CD_CONSERVATIVE_SOURCE_PARAMETER_ALIGNMENT_'):
            left,right=_restore(a[name]),_restore(b[name]);meta=name.replace('PY_S11CD_','PY_S11CD_METADATA_',1)
            if Counter(map(sp.srepr,left))!=Counter(map(sp.srepr,right)):continue
            def aligned(body,parent):
                entries=_restore(body);values={}
                for path,descriptor in entries:
                    key=(sp.srepr(parent[int(path[0])]),tuple(path[1:]))
                    if key in values:raise ValueError('noninjective alignment metadata key')
                    values[key]=sp.srepr(descriptor)
                return values
            if meta in a and meta in b and aligned(a[meta],left)==aligned(b[meta],right):
                classes['sourceAlignmentOrderWithReindexedMetadata'].extend((name,meta));accounted.update((name,meta))
    forward={};reverse={};dummy_conflicts=[];candidates=[]
    for name in sorted(changed-accounted):
        x,y=list(DUMMY.finditer(a[name])),list(DUMMY.finditer(b[name]))
        if not x or len(x)!=len(y):continue
        if DUMMY.sub(lambda m:'Dummy('+m[1]+'dummy_index=IDENTITY)',a[name])!=DUMMY.sub(lambda m:'Dummy('+m[1]+'dummy_index=IDENTITY)',b[name]):continue
        for p,q in zip(x,y):
            if forward.get(p[2],q[2])!=q[2] or reverse.get(q[2],p[2])!=p[2]:dummy_conflicts.append(name)
            if p[2]!=q[2] and (p[2] in unchanged_dummy_ids or q[2] in unchanged_dummy_ids):dummy_conflicts.append(name)
            forward[p[2]]=q[2];reverse[q[2]]=p[2]
        candidates.append(name)
    if not dummy_conflicts:
        for name in candidates:
            renamed=DUMMY.sub(lambda m:'Dummy('+m[1]+'dummy_index='+forward[m[2]]+')',a[name])
            if renamed==b[name]:classes['injectiveDummyIdentityRenaming'].append(name);accounted.add(name)
    unclassified=sorted(changed-accounted)
    return {'baseline':str(before),'baselineSha256':raw['baselineSha256'],'transcript':str(after),'sha256':raw['sha256'],
        'baselineTagCount':raw['baselineTagCount'],'newTagCount':raw['newTagCount'],
        'unchangedPayloadCount':raw['unchangedPayloadCount'],'addedTagCount':len(raw['addedTags']),
        'rawChangedPayloadCount':len(changed),'missingTags':raw['missingTags'],
        'classificationCounts':{k:len(v) for k,v in classes.items()},'classifications':dict(classes),
        'dummyIdentityMap':forward,'dummyConflicts':dummy_conflicts,'unclassifiedChanges':unclassified,
        'nativeExistingValueDifferences':native['nativePreservation']['differentExistingValues'],
        'nativeMissingTags':native['nativePreservation']['missingTags']}


def run():
    parser=argparse.ArgumentParser();parser.add_argument('raw_inventory',type=Path);parser.add_argument('native_inventory',type=Path);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    raw=json.loads(args.raw_inventory.read_text());native=json.loads(args.native_inventory.read_text())['results'][0]
    result=inspect(raw,native)
    result.update(instrumentSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),wallSeconds=time.monotonic()-started,
        inputInventorySha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.raw_inventory,args.native_inventory)})
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='classifications'},indent=2))
    if any(result[k] for k in ('missingTags','dummyConflicts','unclassifiedChanges','nativeExistingValueDifferences','nativeMissingTags')):raise SystemExit(1)


if __name__=='__main__':run()
