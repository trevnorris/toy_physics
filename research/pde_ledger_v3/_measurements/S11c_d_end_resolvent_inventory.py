#!/usr/bin/env python3
"""Inventory end-resolvent objects, dimensions and contour/bank residuals."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import resource
import time

from S11c_d_joint_sheet_inventory import association,portable,binding_order_comparison
from S11c_d_joint_sheet_check import engine,_restore
from S11c_d_output_codec import decoded_lines
import sympy as sp


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect(path):
    tags=set();duplicates=[];pending={};gaps=[];nonfinite=[];packets={};constraints={}
    residuals={};old={};metadata_counts={}
    for line in decoded_lines(path):
        tag,sep,payload=line.rstrip('\n').partition(': ')
        if not sep:continue
        tag=tag.removeprefix('PY_S11CD_')
        if tag in tags:duplicates.append(tag)
        tags.add(tag)
        if tag.startswith(('JOINT_SHEET_','METADATA_JOINT_SHEET_')):
            old[tag]=hashlib.sha256(payload.encode()).hexdigest()
        if tag in ('REDUCED_DIMENSION_CONSTRAINT_RESIDUALS','REDUCED_DIMENSION_UNRESOLVED',
                   'PENCIL_DIMENSION_CONSTRAINT_RESIDUALS','END_RESOLVENT_PREFLIGHT_DIMENSION_CONSTRAINTS'):
            constraints[tag]=portable(_restore(payload))
        if tag.startswith('METADATA_END_RESOLVENT_'):
            name=tag.removeprefix('METADATA_');groups=_restore(payload);units={}
            axis=association(groups)
            if str(axis.get('REPRESENTATION',''))=='MATRIX_AXIS_DIMENSIONS':
                rows,columns=map(int,axis['SHAPE'])
                r,c=axis['ROW_DIMENSIONS_L_T_M'],axis['COLUMN_DIMENSION_OFFSETS_L_T_M']
                if len(r)!=rows or len(c)!=columns:gaps.append([name,'axis shape mismatch'])
                if any(len(v)!=3 for v in (*r,*c)):gaps.append([name,'dimension vector shape'])
                if any(any(v)for v in axis['AXIS_ENCODING_RESIDUAL']):gaps.append([name,'axis encoding residual'])
                if not {'MULTIGRADE','EPSILON_LAMBDA_SUPPORT'}<=set(axis):gaps.append([name,'missing grade'])
                for i in range(rows):
                    for j in range(columns):units[(i*columns+j,)]=tuple(int(a+b)for a,b in zip(r[i],c[j]))
                value=pending.pop(name,None)
                if value is None:gaps.append([name,'missing object'])
                elif not {'OBJECT_SHA256','OBJECT_SHA_AND_NUMERIC_PIT'} & set(association(value)):
                    gaps.append([name,'axis metadata without tensor fingerprint'])
                metadata_counts[name]=len(units)
                continue
            for group in groups:
                g=association(group)
                if not {'PATHS','DIMENSION_L_T_M','MULTIGRADE','EPSILON_LAMBDA_SUPPORT'}<=set(g):
                    gaps.append([name,'missing metadata field']);continue
                for p in g['PATHS']:
                    key=tuple(str(v)if isinstance(v,sp.core.symbol.Str)else int(v) for v in p)
                    if key in units:gaps.append([name,str(key),'duplicate path'])
                    units[key]=tuple(int(v)for v in g['DIMENSION_L_T_M'])
            value=pending.pop(name,None)
            if value is None:gaps.append([name,'missing object'])
            elif not {'OBJECT_SHA256','OBJECT_SHA_AND_NUMERIC_PIT'} & set(association(value)):
                for p,v in engine.leaves(value):
                    if isinstance(v,sp.core.symbol.Str):continue
                    if p not in units:gaps.append([name,str(p),'missing dimension'])
                    if v.has(sp.nan,sp.zoo,sp.oo,-sp.oo):nonfinite.append([name,str(p),str(v)])
            metadata_counts[name]=len(units)
            continue
        if not tag.startswith('END_RESOLVENT_'):continue
        value=_restore(payload);pending[tag]=value
        match=re.match(r'(END_RESOLVENT_(?:INPUT|PIT)_.+?_0)_(.*)',tag)
        if not match:continue
        prefix,label=match.groups()
        packet=packets.setdefault(prefix,{'poles':{},'contours':{},'banks':{},'scalarResiduals':{}})
        data=association(value)
        pole=re.fullmatch(r'POLE_(\d+)_RECORD',label)
        contour=re.fullmatch(r'POLE_(\d+)_CONTOUR_(\d+)_(\d+)_RECORD',label)
        if pole:packet['poles'][pole[1]]=portable(value)
        elif contour:packet['contours'][label]=portable(value)
        elif label.endswith('_BANK_RECORDS'):packet['banks'][label]=portable(value)
        elif label in ('SUMMARY','BANK_SUMMARY','SCOPE','UNIT_FRAME','GRADE_ORIGIN','BOUND_CARRIERS','FREQUENCY',
                       'BRANCH_POINTS','DENOMINATOR_POINTS','INPUT_BINDING','COEFFICIENT_COORDINATE_UNITS'):
            packet[label]=portable(value)
        elif label.endswith('_MAXIMA'):
            category=re.sub(r'^POLE_\d+_(?:CONTOUR_\d+_\d+_)?','',label.removesuffix('_MAXIMA'))
            category=re.sub(r'^(FREQUENCY|MOMENTUM)_CUT_BANK_\d+_\d+_(?:M1_|1_)?',r'\1_BANK_',category)
            for group in value:
                g=association(group);dim=tuple(int(v)for v in g['DIMENSION_L_T_M']);v=float(g['MAXIMUM_ABSOLUTE'])
                key=category+'|'+str(dim)
                if key not in residuals or v>residuals[key]['maximum']:
                    residuals[key]={'maximum':v,'tag':tag,'dimension':dim}
        elif ('RESIDUAL'in label) and not data:
            numbers=[abs(complex(v))for _,v in engine.leaves(value)if getattr(v,'is_number',False)]
            if numbers:packet['scalarResiduals'][label]=max(numbers)
    statuses=Counter(p.get('STATUS','UNRESOLVED')for packet in packets.values()for p in packet['poles'].values())
    summary_keys=('NATIVE_CANDIDATE_COUNT','RESIDUE_COUNT','CONTOUR_CANDIDATE_COUNT','NULLITY_DIFFERENCE_COUNT',
                  'BANK_PAIR_COUNT','BANK_INVERSE_COUNT','SUBTRACTION_UNRESOLVED_COUNT','ORIGINAL_UNRESOLVED_CANDIDATE_COUNT')
    totals={key:sum(p.get('SUMMARY',{}).get(key,0)for p in packets.values())for key in summary_keys}
    contours=[c for p in packets.values()for c in p['contours'].values()]
    contour_diagnostics={key:max((abs(complex(c[key]['real'],c[key]['imag']))if isinstance(c[key],dict)else abs(c[key])
        for c in contours if key in c),default=0)for key in ('MAXIMUM_RADICAL_RESIDUAL','MAXIMUM_COEFFICIENT_CONDITION',
        'MAXIMUM_COEFFICIENT_INVERSE_RESIDUAL','ROOT_CLOSURE_RESIDUAL','SEED_REFINEMENT_DIFFERENCE','SEED_ODE_DIFFERENCE')}
    return {'transcript':str(path.resolve()),'bytes':path.stat().st_size,'sha256':digest(path),'uniqueTags':len(tags),
        'duplicates':duplicates,'metadataGaps':gaps,'unmatchedObjects':sorted(pending),'nonfiniteObjects':nonfinite,
        'dimensionalConstraints':constraints,'packetCount':len(packets),'packets':packets,'totals':totals,
        'poleStatuses':dict(statuses),'contourCount':len(contours),'contourDiagnosticMaxima':contour_diagnostics,
        'residualMaximaByDimension':residuals,'oldContinuationPayloadDigests':old,
        'completionMarkers':int('PROCESS_COMPLETION'in tags),'metadataObjectCount':len(metadata_counts)}


def run():
    started=time.monotonic()
    parser=argparse.ArgumentParser()
    parser.add_argument('transcripts',nargs='+',type=Path)
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    results=[inspect(p)for p in args.transcripts]
    if args.baseline:
        baseline=inspect(args.baseline)
        for result in results:
            differences=[k for k,v in baseline['oldContinuationPayloadDigests'].items()
                         if result['oldContinuationPayloadDigests'].get(k)!=v]
            ordering=binding_order_comparison(args.baseline,Path(result['transcript']),differences)
            result['baselineComparison']={'path':str(args.baseline.resolve()),'sha256':baseline['sha256'],
                'rawDifferences':differences,'associationOrderOnlyDifferences':ordering,
                'semanticDifferences':sorted(set(differences)-set(ordering))}
    record={'instrumentSha256':digest(Path(__file__)),'results':results,'wallSeconds':time.monotonic()-started,
            'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'transcripts':len(results),'packets':sum(v['packetCount']for v in results),
        'metadataGaps':sum(len(v['metadataGaps'])for v in results),'nonfiniteObjects':sum(len(v['nonfiniteObjects'])for v in results),
        'wallSeconds':record['wallSeconds']}))


if __name__=='__main__':run()
