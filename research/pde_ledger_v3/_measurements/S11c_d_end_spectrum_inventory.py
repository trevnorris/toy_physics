#!/usr/bin/env python3
"""Read emitted end-spectrum evidence without supplying a physical disposition."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys
import time
import resource

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from ledger_fold import _restore
from S11c_d_mixing_scattering_sympy_audit import leaves


def association(value):
    return {str(k):v for k,v in value}


def portable(value):
    if isinstance(value,dict):return {str(k):portable(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,sp.Tuple)):return [portable(v) for v in value]
    if isinstance(value,sp.logic.boolalg.BooleanAtom):return bool(value)
    if isinstance(value,sp.Integer):return int(value)
    if isinstance(value,sp.Basic):return str(value)
    return value


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def input_json(value):
    if isinstance(value,sp.Tuple):
        if value and all(isinstance(v,sp.Tuple) and len(v)==2 and isinstance(v[0],sp.core.symbol.Str) for v in value):
            return {str(k):input_json(v) for k,v in value}
        return [input_json(v) for v in value]
    return portable(value)


def read(path):
    tags=set();duplicates=[];packets={};metadata=set();new_tags=set();pending={};residuals={}
    legacy_symbols={};legacy_packets={};coverage={};source_pins=None;constraints={}
    metadata_gaps=[]
    specification=None;emitted_input_digest=None
    for line in path.open():
        tag,separator,payload=line.rstrip('\n').partition(': ')
        if not separator:continue
        tag=tag.removeprefix('PY_S11CD_')
        if tag in tags:duplicates.append(tag)
        tags.add(tag)
        if tag in ('CHANNEL_INPUT_SPECIFICATION','END_SPECTRUM_PREFLIGHT_INPUT'):
            specification=input_json(_restore(payload))
        if tag=='CHANNEL_INPUT_SHA256':emitted_input_digest=str(_restore(payload))
        if tag.startswith('FULL_SECTOR_PENCIL_SYMBOL_'):
            legacy_symbols[tag]=hashlib.sha256(payload.encode()).hexdigest()
        if tag.startswith('FULL_PENCIL_MODE_PIT_'):
            legacy_packets[tag]={'digest':hashlib.sha256(payload.encode()).hexdigest(),
                'roots':[(complex(association(v)['K']),complex(association(v)['Q']),
                          bool(association(v).get('HELMHOLTZ_CHART_DOMAIN_LIMITATION',False))) for v in _restore(payload)]}
        if tag=='BUILD_INPUT_DIGESTS':source_pins=portable(association(_restore(payload)))
        if tag in ('REDUCED_DIMENSION_CONSTRAINT_RESIDUALS','REDUCED_DIMENSION_UNRESOLVED','PENCIL_DIMENSION_CONSTRAINT_RESIDUALS'):
            constraints[tag]=portable(_restore(payload))
        if tag.startswith('METADATA_END_SPECTRUM_'):
            base=tag.removeprefix('METADATA_');metadata.add(base)
            groups=_restore(payload)
            covered=[];units={}
            for group in groups:
                entry=association(group)
                paths=[tuple(str(x) if isinstance(x,sp.core.symbol.Str) else int(x) for x in p) for p in entry['PATHS']]
                dimension=tuple(int(x) for x in entry['DIMENSION_L_T_M'])
                covered.extend(paths);units.update({p:dimension for p in paths})
                if 'MULTIGRADE' not in entry or 'EPSILON_LAMBDA_SUPPORT' not in entry:
                    metadata_gaps.append((base,'missing grades'))
            if base in pending:
                value=pending.pop(base)
                kind=re.sub(r'^.*_MODE_\d+_','',base)
                if kind not in ('RIGHT_RESIDUAL','LEFT_RESIDUAL','PROJECTOR_RESIDUAL'):
                    kind=base.rsplit('_',2)[-2]+'_RESIDUAL'
                for p,x in leaves(value):
                    if isinstance(x,sp.core.symbol.Str):continue
                    if p not in units:metadata_gaps.append((base,str(p)))
                    dimension=units.get(p,('MISSING',))
                    key=kind+'|'+str(dimension)
                    residuals[key]=max(residuals.get(key,0.),abs(complex(x)))
            continue
        if not tag.startswith('END_SPECTRUM_'):continue
        new_tags.add(tag)
        value=_restore(payload)
        if tag.endswith('_RESIDUAL') and isinstance(value,(sp.MatrixBase,sp.Number)):
            pending[tag]=value
        match=re.match(r'(END_SPECTRUM_(?:INPUT|PIT)_.+_\d+)_MODE_(\d+)_RECORD$',tag)
        if match:
            key,index=match.groups();record=association(value)
            packet=packets.setdefault(key,{'modes':[]})
            packet['modes'].append({'index':int(index),'k':[complex(record['K']).real,complex(record['K']).imag],
                'q':[complex(record['Q']).real,complex(record['Q']).imag], 'nullity':int(record.get('NULLITY',0)),
                'multiplicity':int(record['MULTIPLICITY']),
                'multiplicityDifference':int(record.get('ALGEBRAIC_GEOMETRIC_MULTIPLICITY_DIFFERENCE',-1)),
                'sheet':str(record['FIXED_FREQUENCY_SHEET_MEMBERSHIP']),
                'classifier':str(record.get('CLASSIFIER_STATUS','UNRESOLVED')),
                'radicalResidual':abs(complex(record['RADICAL_RESIDUAL']))})
        elif tag.endswith('_ROOT_COVERAGE'):
            key=tag.removesuffix('_ROOT_COVERAGE');certificate=association(value)
            disks=[association(v) for v in certificate['ROOT_DISKS']]
            coverage[key]={'polynomialDegree':int(certificate['DEGREE']),
                'distinctRadicalRoots':int(certificate['DISTINCT_ROOT_COUNT']),
                'countWithMultiplicity':int(certificate['COUNT_WITH_MULTIPLICITY']),
                'degreeCountResidual':int(certificate['DEGREE_COUNT_RESIDUAL']),
                'allDisksDisjoint':bool(certificate['ALL_DISKS_DISJOINT']),
                'allDisksIsolated':all(bool(v['ONE_ROOT_DISK']) for v in disks),
                'boundSigns':dict(Counter(str(v['EXACT_BOUND_DIFFERENCE_SIGN']) for v in disks)),
                'precisionRefinementMaximum':max((abs(complex(v['PRECISION_REFINEMENT_DIFFERENCE'])) for v in disks),default=0.),
                'finitePolynomialRootCoverage':bool(certificate['FINITE_POLYNOMIAL_ROOT_COVERAGE'])}
        elif tag.endswith('_SUMMARY'):
            packets.setdefault(tag.removesuffix('_SUMMARY'),{'modes':[]})['summary']=portable(association(value))
        elif tag.endswith('_EXCEPTION_DEGREES'):
            packets.setdefault(tag.removesuffix('_EXCEPTION_DEGREES'),{'modes':[]})['exceptionDegrees']=portable(association(value))
        elif tag.endswith('_POLYNOMIAL_DEGREES'):
            packets.setdefault(tag.removesuffix('_POLYNOMIAL_DEGREES'),{'modes':[]})['polynomialDegrees']=portable(association(value))
        elif tag.startswith(('END_SPECTRUM_INPUT_','END_SPECTRUM_PIT_')):
            for ending,label in (('_GRADE_ORIGIN','gradeOrigin'),('_BOUND_CARRIERS','boundCarriers'),('_UNIT_FRAME','unitFrame')):
                if tag.endswith(ending):
                    packet=packets.setdefault(tag.removesuffix(ending),{'modes':[]})
                    packet[label]=([str(v) for v in value] if label=='unitFrame' else
                                   {str(k):str(v) for k,v in value})
                    break
    comparisons=[]
    input_digest=(hashlib.sha256(json.dumps(specification,sort_keys=True,separators=(',',':')).encode()).hexdigest()
                  if specification is not None else None)
    for key,packet in packets.items():
        if key.startswith('END_SPECTRUM_INPUT_') and specification is not None:
            packet['inputBinding']={'profiles':specification['profiles'],'unitFrame':specification['unit_frame'],
                'parameters':specification['parameters'],'canonicalInputSha256':input_digest}
        if key.startswith('END_SPECTRUM_PIT_') and '_REFERENCE_' not in key:
            legacy='FULL_PENCIL_MODE_PIT_'+key.removeprefix('END_SPECTRUM_PIT_')
            if legacy in legacy_packets:
                old=[(k,q) for k,q,chart in legacy_packets[legacy]['roots'] if not chart]
                new=[(complex(*v['k']),complex(*v['q'])) for v in packet['modes']]
                gaps=[min((max(abs(a-c),abs(b-d)) for c,d in old),default=float('inf')) for a,b in new]
                comparisons.append({'packet':key,'legacyRegularCount':len(old),'nativeCount':len(new),
                                    'maximumNearestRootDifference':max(gaps,default=0.)})
    return {'path':str(path.resolve()),'bytes':path.stat().st_size,'sha256':digest(path),
        'uniqueTags':len(tags),'duplicates':duplicates,'newSpectrumTags':len(new_tags),
        'missingMetadata':sorted(new_tags-metadata),'metadataGaps':metadata_gaps,
        'unmatchedResidualMetadata':sorted(pending),'dimensionConstraints':constraints,
        'sourcePins':source_pins,'packets':packets,'coverage':coverage,
        'emittedInputSha256':emitted_input_digest,'recomputedInputSha256':input_digest,
        'inputDigestMismatch':emitted_input_digest is not None and emitted_input_digest!=input_digest,
        'maximumResidualByDimension':residuals,'legacyRegularRootComparisons':comparisons,
        'legacySymbolDigests':legacy_symbols,'legacyPacketDigests':{k:v['digest'] for k,v in legacy_packets.items()},
        'completionMarkers':int('PROCESS_COMPLETION' in tags)}


def run():
    started=time.monotonic()
    p=argparse.ArgumentParser()
    p.add_argument('transcripts',type=Path,nargs='+');p.add_argument('--baseline',type=Path)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    results=[read(path) for path in args.transcripts]
    if args.baseline:
        baseline=read(args.baseline)
        for result in results:
            result['baselineComparison']={'path':str(args.baseline.resolve()),'sha256':baseline['sha256'],
                'symbolDigestDifferences':[k for k,v in baseline['legacySymbolDigests'].items()
                                           if result['legacySymbolDigests'].get(k)!=v],
                'packetDigestDifferences':[k for k,v in baseline['legacyPacketDigests'].items()
                                           if result['legacyPacketDigests'].get(k)!=v]}
    report={'instrumentSha256':digest(Path(__file__)),'results':results,
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'transcripts':len(results),'packets':sum(len(v['packets']) for v in results),
                     'missingMetadata':sum(len(v['missingMetadata']) for v in results),
                     'wallSeconds':report['wallSeconds'],'peakRssKiB':report['peakRssKiB']}))

if __name__=='__main__':run()
