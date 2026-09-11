#!/usr/bin/env python3
"""Inventory emitted joint-path operands, dimensions, residuals and banks."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import resource
import sys
import time

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from ledger_fold import _restore
from S11c_d_output_codec import decoded_lines
from S11c_d_mixing_scattering_sympy_audit import leaves
from S11c_d_end_spectrum_inventory import input_json


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def binding_order_comparison(left,right,names):
    """Compare association keys and metadata path sets; retain raw differences."""
    def collect(path):
        values={}
        for line in decoded_lines(path):
            tag,separator,payload=line.rstrip('\n').partition(': ')
            tag=tag.removeprefix('PY_S11CD_')
            if separator and tag in names:values[tag]=payload
        return values
    a,b=collect(left),collect(right)
    def canonical(tag,payload):
        value=_restore(payload)
        if tag.startswith('METADATA_'):
            groups=[]
            for group in value:
                entry=association(group)
                if set(entry)!={'PATHS','DIMENSION_L_T_M','MULTIGRADE','EPSILON_LAMBDA_SUPPORT'}:
                    return sp.srepr(value)
                groups.append((tuple(sorted(sp.srepr(p) for p in entry['PATHS'])),
                    sp.srepr(entry['DIMENSION_L_T_M']),
                    tuple(sorted(sp.srepr(g) for g in entry['MULTIGRADE'])),
                    tuple(sorted(sp.srepr(g) for g in entry['EPSILON_LAMBDA_SUPPORT']))))
            return tuple(sorted(groups))
        entry=association(value)
        if len(entry)!=len(value):return sp.srepr(value)
        return tuple(sorted((key,sp.srepr(v)) for key,v in entry.items()))
    return [name for name in names if name.endswith('_BOUND_CARRIERS') and name in a and name in b
            and canonical(name,a[name])==canonical(name,b[name])]


def association(value):
    if not isinstance(value,sp.Tuple):return {}
    if not all(isinstance(v,sp.Tuple) and len(v)==2 and isinstance(v[0],sp.core.symbol.Str) for v in value):return {}
    return {str(k):v for k,v in value}


def portable(value):
    if isinstance(value,dict):return {str(k):portable(v) for k,v in value.items()}
    assoc=association(value)
    if assoc:return portable(assoc)
    if isinstance(value,(tuple,list,sp.Tuple)):return [portable(v) for v in value]
    if isinstance(value,sp.logic.boolalg.BooleanAtom):return bool(value)
    if isinstance(value,sp.Integer):return int(value)
    if isinstance(value,sp.Float):return float(value)
    if isinstance(value,sp.core.symbol.Str):return str(value)
    if isinstance(value,sp.Basic) and value.is_number and value.is_finite:
        z=complex(value)
        return float(z.real) if z.imag==0 else {'real':z.real,'imag':z.imag}
    return str(value)


def inspect(path):
    tags=set();duplicates=[];pending={};metadata=set();gaps=[];packets={};native={};legacy={};pins=None
    residuals={};constraints={};nonfinite=[];specification=None
    for line in decoded_lines(path):
        tag,separator,payload=line.rstrip('\n').partition(': ')
        if not separator:continue
        tag=tag.removeprefix('PY_S11CD_')
        if tag in tags:duplicates.append(tag)
        tags.add(tag)
        if tag.startswith(('END_SPECTRUM_INPUT_','END_SPECTRUM_PIT_','METADATA_END_SPECTRUM_INPUT_','METADATA_END_SPECTRUM_PIT_')):
            native[tag]=hashlib.sha256(payload.encode()).hexdigest()
        if tag.startswith(('FULL_SECTOR_PENCIL_SYMBOL_','FULL_PENCIL_MODE_PIT_')):
            legacy[tag]=hashlib.sha256(payload.encode()).hexdigest()
        if tag=='BUILD_INPUT_DIGESTS':pins=portable(_restore(payload))
        if tag=='CHANNEL_INPUT_SPECIFICATION':specification=input_json(_restore(payload))
        if tag in ('REDUCED_DIMENSION_CONSTRAINT_RESIDUALS','REDUCED_DIMENSION_UNRESOLVED','PENCIL_DIMENSION_CONSTRAINT_RESIDUALS'):
            constraints[tag]=portable(_restore(payload))
        if tag.startswith('METADATA_JOINT_SHEET_'):
            base=tag.removeprefix('METADATA_');metadata.add(base)
            groups=_restore(payload);units={}
            for group in groups:
                data=association(group)
                for p in data['PATHS']:
                    key=tuple(str(v) if isinstance(v,sp.core.symbol.Str) else int(v) for v in p)
                    if key in units:gaps.append((base,str(key),'duplicate metadata path'))
                    units[key]=tuple(int(v) for v in data['DIMENSION_L_T_M'])
                if 'MULTIGRADE' not in data or 'EPSILON_LAMBDA_SUPPORT' not in data:gaps.append((base,'missing grade'))
            value=pending.pop(base,None)
            if value is None:gaps.append((base,'missing object'))
            elif 'OBJECT_SHA256' not in association(value):
                for p,x in leaves(value):
                    if isinstance(x,sp.core.symbol.Str):continue
                    if p not in units:gaps.append((base,str(p),'missing dimension'))
                    if x.has(sp.nan,sp.zoo):nonfinite.append((base,str(p),str(x)))
                    if ('RESIDUAL' in base or 'REFINEMENT' in base) and x.is_number and x.is_finite:
                        key=base.split('_0_',1)[-1]+'|'+str(units.get(p))
                        residuals[key]=max(residuals.get(key,0.),abs(complex(x)))
            continue
        if not tag.startswith('JOINT_SHEET_'):continue
        value=_restore(payload);pending[tag]=value;data=association(value)
        if tag.endswith('_DIMENSION_CONSTRAINTS'):constraints[tag]=portable(value)
        match=re.match(r'(JOINT_SHEET_(?:INPUT|PIT)_.+?_0)_(.*)',tag)
        if not match:continue
        packet_name,label=match.groups()
        packet=packets.setdefault(packet_name,{'paths':{},'banks':{},'literalMatrixJoins':[],'matrixJumpFingerprints':{}})
        if 'VERTICES' in data:
            record={key:portable(data[key]) for key in ('STATUS','PATH_DEFINED','SEED_Q','END_Q','ODE_END_Q',
                'REFINEMENT_DIFFERENCE','ODE_DIFFERENCE','ODE_MAXIMUM_RADICAL_RESIDUAL','CLOSED_COORDINATE_PATH',
                'END_TO_SEED_RATIO','RADICAND_ARGUMENT_TURNS') if key in data}
            record['vertices']=portable(data['VERTICES'])
            record['cutEncounters']=sum(len(association(s)['DOWNWARD_FREQUENCY_CUT_ENCOUNTERS']) for s in data.get('SEGMENTS',()))
            record['branchIntersections']=sum(bool(association(s)['BRANCH_INTERSECTION']) for s in data.get('SEGMENTS',()))
            packet['paths'][label]=record
        elif '_BANK_OPERANDS_' in label:packet['banks'][label]=portable(value)
        elif label.startswith(('REAL_AXIS_MATRIX_RESIDUAL_','REAL_AXIS_REFINED_MATRIX_RESIDUAL_')):
            packet['literalMatrixJoins'].append({'label':label,'entries':len(value),
                'maximumAbsoluteResidual':max((abs(complex(v)) for v in value),default=0.),
                'nonzero':[(i,str(v)) for i,v in enumerate(value) if v!=0]})
        elif label.startswith('REAL_AXIS_JOIN_PRECISION_'):packet[label]=portable(value)
        elif '_MATRIX_JUMP_' in label:packet['matrixJumpFingerprints'][label]=portable(value)
        elif label in ('SUMMARY','UNIT_FRAME','BOUND_CARRIERS','GRADE_ORIGIN','SOURCE_GRADE_SUPPORT',
                       'INPUT_BINDING','MOMENTUM_BANK_TARGET_SOURCE','MOMENTUM_BANK_TARGETS','ORIGINAL_UNRESOLVED_CANDIDATES'):
            packet[label]=portable(value)
    input_digest=hashlib.sha256(json.dumps(specification,sort_keys=True,separators=(',',':')).encode()).hexdigest() if specification else None
    statuses=Counter(p['STATUS'] for v in packets.values() for p in v['paths'].values())
    return {'transcript':str(path.resolve()),'bytes':path.stat().st_size,'sha256':digest(path),
        'uniqueTags':len(tags),'duplicates':duplicates,'metadataGaps':gaps,'unmatchedObjects':sorted(pending),
        'nonfiniteObjects':nonfinite,'packets':packets,'pathStatuses':dict(statuses),
        'residualMaximumByDimension':residuals,'dimensionalConstraints':constraints,'sourcePins':pins,
        'canonicalInputSha256':input_digest,'completionMarkers':int('PROCESS_COMPLETION' in tags),
        'nativeSpectrumPayloadDigests':native,'legacyPayloadDigests':legacy}


def run():
    started=time.monotonic()
    parser=argparse.ArgumentParser()
    parser.add_argument('transcripts',type=Path,nargs='+')
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    results=[inspect(p) for p in args.transcripts]
    if args.baseline:
        baseline=inspect(args.baseline)
        for result in results:
            differences=[k for k,v in baseline['nativeSpectrumPayloadDigests'].items()
                         if result['nativeSpectrumPayloadDigests'].get(k)!=v]
            ordering=binding_order_comparison(args.baseline,Path(result['transcript']),differences)
            result['baselineComparison']={'path':str(args.baseline.resolve()),'sha256':baseline['sha256'],
                'nativeDifferences':differences,'nativeAssociationOrderOnlyDifferences':ordering,
                'nativeSemanticDifferences':sorted(set(differences)-set(ordering)),
                'legacyDifferences':[k for k,v in baseline['legacyPayloadDigests'].items()
                                     if result['legacyPayloadDigests'].get(k)!=v]}
    record={'instrumentSha256':digest(Path(__file__)),'results':results,'wallSeconds':time.monotonic()-started,
            'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    args.output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'transcripts':len(results),'packets':sum(len(r['packets']) for r in results),
        'metadataGaps':sum(len(r['metadataGaps']) for r in results),'unmatchedObjects':sum(len(r['unmatchedObjects']) for r in results),
        'nonfiniteObjects':sum(len(r['nonfiniteObjects']) for r in results),'wallSeconds':record['wallSeconds']}))


if __name__=='__main__':run()
