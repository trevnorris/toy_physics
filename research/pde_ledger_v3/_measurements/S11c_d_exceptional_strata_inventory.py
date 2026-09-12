#!/usr/bin/env python3
"""Inventory computed exceptional loci, exact spaces and regularity criteria."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import resource
import time

from S11c_d_joint_sheet_inventory import association,portable,sp
from S11c_d_joint_sheet_check import engine,_restore
from S11c_d_output_codec import decoded_lines,restore_emission_index


def native_preservation(before,after):
    def collect(path):
        result={}
        for line in decoded_lines(path):
            tag,_,body=line.rstrip('\n').partition(': ')
            if tag.startswith(('PY_S11CD_END_SPECTRUM_','PY_S11CD_METADATA_END_SPECTRUM_')):
                result[tag]=body
        return result
    a,b=collect(before),collect(after);changed=[];additive=[];missing=[];different=[]
    def metadata(body):
        result={}
        for item in _restore(body):
            entry=association(item)
            descriptor=tuple(sp.srepr(entry[k])for k in ('DIMENSION_L_T_M','MULTIGRADE','EPSILON_LAMBDA_SUPPORT'))
            for path in entry['PATHS']:result[sp.srepr(path)]=descriptor
        return result
    for tag,body in a.items():
        if tag not in b:missing.append(tag);continue
        if body==b[tag]:continue
        changed.append(tag)
        if tag.startswith('PY_S11CD_METADATA_'):
            old,new=metadata(body),metadata(b[tag])
        else:
            old,new=association(_restore(body)),association(_restore(b[tag]))
            if not old or not new:different.append(tag);continue
            old={k:sp.srepr(v)for k,v in old.items()};new={k:sp.srepr(v)for k,v in new.items()}
        if all(new.get(k)==v for k,v in old.items()):additive.append(tag)
        else:different.append(tag)
    return {'baseline':str(before.resolve()),'baselineSha256':hashlib.sha256(before.read_bytes()).hexdigest(),
        'baselineNativeTags':len(a),'newNativeTags':len(b),'missingTags':missing,'rawChangedTags':changed,
        'associationOrderOrAddedFieldOnlyTags':additive,'differentExistingValues':different}


def inspect(path):
    tags=[];tag_set=set();duplicates=[];pending={};gaps=[];nonfinite=[];packets={};native={};constraints={};index_checks=[];residuals={}
    for line in decoded_lines(path):
        full_tag,_,payload=line.rstrip('\n').partition(': ')
        tag=full_tag.removeprefix('PY_S11CD_')
        if full_tag in tag_set:duplicates.append(full_tag)
        if tag=='EMISSION_LINES':
            value=association(_restore(payload))
            if 's11cdIndexedTagCount' in value:
                restored=restore_emission_index(value,tags)
                index_checks.append({'indexedTags':len(restored),'missingAssignments':len(tags)-len(restored)})
        tags.append(full_tag);tag_set.add(full_tag)
        if tag.endswith('DIMENSION_CONSTRAINTS')or tag in ('REDUCED_DIMENSION_CONSTRAINT_RESIDUALS',
            'REDUCED_DIMENSION_UNRESOLVED','PENCIL_DIMENSION_CONSTRAINT_RESIDUALS'):
            constraints[tag]=portable(_restore(payload))
        if tag.startswith('END_SPECTRUM_') and tag.endswith(('_REGULARITY_CRITERIA','_SUMMARY')):
            native[tag]=portable(_restore(payload))
        if tag.startswith(('METADATA_END_EXCEPTIONAL_SLICE_','METADATA_BULK_EXCEPTIONAL_SLICE_')):
            name=tag.removeprefix('METADATA_');groups=_restore(payload);units={}
            for group in groups:
                record=association(group)
                if not {'PATHS','DIMENSION_L_T_M','MULTIGRADE','EPSILON_LAMBDA_SUPPORT'}<=set(record):
                    gaps.append([name,'metadata fields']);continue
                for p in record['PATHS']:
                    key=tuple(str(v) if isinstance(v,sp.core.symbol.Str)else int(v) for v in p)
                    if key in units:gaps.append([name,'duplicate path',str(key)])
                    units[key]=tuple(record['DIMENSION_L_T_M'])
            value=pending.pop(name,None)
            if value is None:gaps.append([name,'missing object'])
            elif {'OBJECT_SHA256','OBJECT_SHA_AND_NUMERIC_PIT'}&set(association(value)):
                for p,v in engine.leaves(value):
                    if not isinstance(v,sp.core.symbol.Str) and v.has(sp.nan,sp.zoo,sp.oo,-sp.oo):
                        nonfinite.append([name,'fingerprint',str(p),str(v)])
            elif not {'OBJECT_SHA256','OBJECT_SHA_AND_NUMERIC_PIT'}&set(association(value)):
                for p,v in engine.leaves(value):
                    if isinstance(v,sp.core.symbol.Str):continue
                    if p not in units:gaps.append([name,'missing dimension',str(p)])
                    if v.has(sp.nan,sp.zoo,sp.oo,-sp.oo):nonfinite.append([name,str(p),str(v)])
                    for category in ('RIGHT_INVERSE_RESIDUAL','LEFT_INVERSE_RESIDUAL','MATRIX_PRECISION_REFINEMENT',
                                     'INVERSE_PRECISION_REFINEMENT','RADICAL_RESIDUAL','RIGHT_RESIDUAL','LEFT_RESIDUAL'):
                        if name.endswith('_'+category) and getattr(v,'is_number',False):
                            dimension=tuple(map(int,units.get(p,())))
                            key=category+'|'+str(dimension);magnitude=abs(complex(v))
                            if key not in residuals or magnitude>residuals[key]['maximum']:
                                residuals[key]={'maximum':magnitude,'tag':name,'dimension':dimension}
            continue
        if not tag.startswith(('END_EXCEPTIONAL_SLICE_','BULK_EXCEPTIONAL_SLICE_')):continue
        value=_restore(payload);pending[tag]=value
        prefix,_,label=tag.partition('_CONSTANT_')
        packet=packets.setdefault(prefix,{'records':{},'statuses':Counter()})
        data=portable(value)
        if label.endswith('STATUS'):packet['statuses'][str(value)]+=1
        if 'STATUS' in association(value):packet['statuses'][str(association(value)['STATUS'])]+=1
        if not ('OBJECT_SHA256' in association(value) or 'OBJECT_SHA_AND_NUMERIC_PIT' in association(value)):
            packet['records'][label]=data
    for p in packets.values():p['statuses']=dict(p['statuses'])
    return {'transcript':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'bytes':path.stat().st_size,'uniqueTags':len(tag_set),'duplicates':duplicates,'metadataGaps':gaps,
        'nonfiniteObjects':nonfinite,'unmatchedObjects':sorted(pending),'dimensionConstraints':constraints,
        'packetCount':len(packets),'packets':packets,'nativeRegularity':native,'sourceIndexChecks':index_checks,
        'residualMaximaByDimension':residuals,
        'completionMarkers':sum(t.endswith(('_EXCEPTIONAL_PREFLIGHT_PROCESS_COMPLETION','_PROCESS_COMPLETION'))
            and not t.startswith('PY_S11CD_METADATA_') for t in tag_set)}


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('transcripts',type=Path,nargs='+')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--baseline',type=Path)
    args=parser.parse_args();start=time.monotonic()
    records=[inspect(p)for p in args.transcripts]
    if args.baseline:
        for path,record in zip(args.transcripts,records):record['nativePreservation']=native_preservation(args.baseline,path)
    result={'instrumentSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'results':records,'wallSeconds':time.monotonic()-start,
            'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':run()
