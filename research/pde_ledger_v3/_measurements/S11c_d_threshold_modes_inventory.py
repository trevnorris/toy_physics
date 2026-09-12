#!/usr/bin/env python3
"""Inventory full threshold chains, local connections, units and residuals."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

from S11c_d_joint_sheet_inventory import association,portable,sp
from S11c_d_joint_sheet_check import engine,_restore
from S11c_d_output_codec import decoded_lines,restore_emission_index


def inspect(path):
    tags=[];seen=set();pending={};duplicates=[];gaps=[];nonfinite=[];maxima={};records={};counts=Counter();index=[]
    secants={}
    for line in decoded_lines(path):
        name,_,payload=line.rstrip('\n').partition(': ')
        tag=name.removeprefix('PY_S11CD_')
        if name in seen:duplicates.append(name)
        if tag=='EMISSION_LINES':
            value=association(_restore(payload))
            if 's11cdIndexedTagCount' in value:
                restored=restore_emission_index(value,tags)
                index.append({'indexedTags':len(restored),'missingAssignments':len(tags)-len(restored)})
        tags.append(name);seen.add(name)
        if tag in ('PROCESS_COMPLETION','EXCEPTIONAL_PREFLIGHT_PROCESS_COMPLETION'):counts['completion']+=1
        if tag.startswith('THRESHOLD_MODE_'):
            value=_restore(payload);pending[tag]=value;counts['objects']+=1
            heavy=bool({'OBJECT_SHA256','OBJECT_SHA_AND_NUMERIC_PIT'}&set(association(value)))
            if heavy:counts['fingerprints']+=1
            if not heavy and tag.endswith(('_CENSUS','_SUMMARY','_DOMAIN','_DOMAINS','_REAL_AXIS_SHEET_JOIN','_LOCAL_VALUATIONS')):
                records[tag]=portable(value)
            if tag.endswith('_RIGHT_CENSUS'):counts['rightChainSpaces']+=1
            if tag.endswith('_LEFT_CENSUS'):counts['leftChainSpaces']+=1
            if '_CHAIN_' in tag and tag.endswith('_LENGTH'):counts['individualChains']+=1
            if '_NODE_' in tag and tag.endswith('_RANK_DOMAIN'):counts['pathModeNodes']+=1
            if '_APPROACH_' in tag and tag.endswith('_RANK_DOMAIN'):counts['approachModeNodes']+=1
            if tag.endswith('_REAL_AXIS_SHEET_JOIN'):
                counts['sourceSheetMatches' if association(value)['MATCHES_REDUCED_SEED']==sp.true else 'oppositeSheetPoints']+=1
            if tag.endswith('_NORMAL_PATH'):
                r=association(value);counts['normalPaths']+=1;counts['definedNormalPaths']+=int(r['PATH_DEFINED']==sp.true)
                if '_LOOP_' in tag:records[tag]=portable(value)
            if tag.endswith('_BULK_PATH'):
                r=association(value);counts['bulkPaths']+=1;counts['definedBulkPaths']+=int(r['PATH_DEFINED']==sp.true)
            continue
        if not tag.startswith('METADATA_THRESHOLD_MODE_'):continue
        object_tag=tag.removeprefix('METADATA_')
        groups=_restore(payload);units={}
        for group in groups:
            r=association(group)
            if not {'PATHS','DIMENSION_L_T_M','MULTIGRADE','EPSILON_LAMBDA_SUPPORT'}<=set(r):
                gaps.append([object_tag,'metadata fields']);continue
            for path_entry in r['PATHS']:
                key=tuple(str(v) if isinstance(v,sp.core.symbol.Str) else int(v) for v in path_entry)
                if key in units:gaps.append([object_tag,'duplicate path',str(key)])
                units[key]=tuple(map(int,r['DIMENSION_L_T_M']))
        value=pending.pop(object_tag,None)
        if value is None:gaps.append([object_tag,'missing object']);continue
        heavy=bool({'OBJECT_SHA256','OBJECT_SHA_AND_NUMERIC_PIT'}&set(association(value)))
        for p,v in engine.leaves(value):
            if isinstance(v,sp.core.symbol.Str):continue
            if v.has(sp.nan,sp.zoo,sp.oo,-sp.oo):nonfinite.append([object_tag,str(p),str(v)])
            if heavy:continue
            if p not in units:gaps.append([object_tag,'missing dimension',str(p)])
            if not getattr(v,'is_number',False):continue
            category=None
            if 'RESIDUAL' in object_tag:
                for candidate in ('PLANE_EQUATION_RESIDUAL','EQUATION_RESIDUAL','RADICAL_TAYLOR_RESIDUAL',
                    'DETERMINANT_CHAIN_MULTIPLICITY_RESIDUAL','NORMAL_ENDPOINT_RESIDUAL','NODE_BULK_ENDPOINT_RESIDUAL',
                    'BULK_PATH_ENDPOINT_RESIDUAL','BULK_SEED_RESIDUAL','RADICAL_RESIDUAL','NORMAL_SQUARE_RESIDUAL',
                    'GAUGE_RESIDUAL','SOURCE_ROOT_RESIDUAL','THRESHOLD_RESIDUAL'):
                    if candidate in object_tag:category=candidate;break
                if category is None:category='OTHER_RESIDUAL'
            elif 'MATRIX_REFINEMENT' in object_tag:category='MATRIX_REFINEMENT'
            if category:
                dimension=units.get(p,());key=category+'|'+str(dimension);magnitude=abs(complex(v))
                if key not in maxima or magnitude>maxima[key]['maximum']:
                    maxima[key]={'maximum':magnitude,'dimension':dimension,'tag':object_tag}
                if '_APPROACH_' in object_tag and object_tag.endswith(('_RIGHT_RESIDUAL','_LEFT_RESIDUAL')):
                    record=secants.setdefault(object_tag,{'maximum':0.,'dimensions':set()})
                    record['maximum']=max(record['maximum'],magnitude);record['dimensions'].add(dimension)
    for record in secants.values():record['dimensions']=sorted(record['dimensions'])
    return {'transcript':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'bytes':path.stat().st_size,'uniqueTags':len(seen),'duplicates':duplicates,'metadataGaps':gaps,
        'nonfiniteObjects':nonfinite,'unmatchedObjects':sorted(pending),'counts':dict(counts),
        'records':records,'residualMaximaByDimension':maxima,'secantResiduals':secants,'sourceIndexChecks':index}


def run():
    parser=argparse.ArgumentParser();parser.add_argument('transcripts',type=Path,nargs='+')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();start=time.monotonic()
    data={'instrumentSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'results':[inspect(path) for path in args.transcripts],'wallSeconds':time.monotonic()-start}
    args.output.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':run()
