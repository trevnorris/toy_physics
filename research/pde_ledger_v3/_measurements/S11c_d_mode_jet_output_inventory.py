#!/usr/bin/env python3
"""Inventory the full transcript and its stored rectangular-jet residuals."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys
from sympy import Integer
from sympy.core.symbol import Str

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from ledger_fold import _restore


def association(value):
    return {str(k):v for k,v in value}


def path_component(value):
    return str(value) if isinstance(value,Str) else int(value) if isinstance(value,Integer) else value


def spectrum(record,packet,index):
    return {'packet':packet,'candidate':index,'k':str(record['K']),
            'q':str(record['Q']),'nullity':int(record.get('NULLITY',0)),
            'sheet':str(record['PHYSICAL_BULK_SHEET'])}


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('transcript',type=Path)
    parser.add_argument('--baseline',type=Path)
    args=parser.parse_args()
    baseline=subprocess.run([sys.executable,str(ROOT/'_measurements/S11c_d_sheet_output_inventory.py'),
                             str(args.transcript)],capture_output=True,text=True)
    result=json.loads(baseline.stdout)
    pending={}
    states=Counter()
    clusters=Counter()
    residuals={}
    missing=[]
    spectra=[]
    symbol_payloads={}
    fingerprint_residuals={}
    tags=set()
    fingerprints=[]
    for line in args.transcript.open():
        tag,separator,value=line.partition(': ')
        if not separator:
            continue
        tags.add(tag)
        if tag.startswith('PY_S11CD_FULL_SECTOR_PENCIL_SYMBOL_'):
            symbol_payloads[tag]=value.strip()
        prefixes=('PY_S11CD_FIRST_JET_ROUTE_RIGHT_RESIDUAL_',
                  'PY_S11CD_FIRST_JET_ROUTE_LEFT_RESIDUAL_',
                  'PY_S11CD_FIRST_JET_ROUTE_ROOT_RESIDUAL_',
                  'PY_S11CD_CLASSIFIER_PROJECTOR_RECTANGLE_RESIDUAL_')
        for prefix in prefixes:
            if tag.startswith(prefix):
                fingerprints.append(tag)
                body=association(_restore(value))
                samples=body['NUMERIC_UNIT_FRAME_TENSOR_PIT']
                entry=fingerprint_residuals.setdefault(prefix,{'objects':0,'sample_values':0,
                    'maximum_absolute_numeric_frame_projection':0.})
                entry['objects']+=1
                entry['sample_values']+=len(samples)
                entry['maximum_absolute_numeric_frame_projection']=max(
                    entry['maximum_absolute_numeric_frame_projection'],*(abs(complex(v)) for v in samples))
        if tag.startswith('PY_S11CD_FULL_PENCIL_MODE_'):
            pending[tag]=[association(v) for v in _restore(value)]
        if tag.startswith('PY_S11CD_METADATA_FULL_PENCIL_MODE_'):
            mode_tag=tag.replace('PY_S11CD_METADATA_','PY_S11CD_',1)
            packet=pending.pop(mode_tag)
            metadata=[association(v) for v in _restore(value)]
            units={tuple(map(path_component,p)):tuple(v['DIMENSION_L_T_M'])
                   for v in metadata for p in v['PATHS']}
            for i,record in enumerate(packet):
                spectra.append(spectrum(record,mode_tag,i))
                if 'NULLITY' not in record or int(record['NULLITY'])==0:
                    states['NO_NUMERICAL_NULLSPACE']+=1
                    continue
                status=str(record.get('RECTANGULAR_MODE_JET_DEFINED','MISSING'))
                states[status]+=1
                clusters[str(record['NULLITY'])]+=1
                if status=='MISSING':
                    missing.append((mode_tag,i,'definition_status'))
                if status!='True':
                    continue
                data=association(record['RECTANGULAR_EQUATION_COEFFICIENT_RESIDUALS'])
                for side,grades in data.items():
                    for grade,matrix in association(grades).items():
                        for j,v in enumerate(matrix):
                            path=(i,'RECTANGULAR_EQUATION_COEFFICIENT_RESIDUALS',side,grade,j)
                            if path not in units:
                                missing.append((mode_tag,str(path),'residual_dimension'))
                                continue
                            unit=units[path]
                            key=side+'_'+grade+'_'+','.join(map(str,unit))
                            measured=abs(complex(v))
                            entry=residuals.setdefault(key,{'side':side,'grade':grade,
                                'dimension_l_t_m':list(map(str,unit)),
                                'maximum_absolute_coefficient':0.,'entries':0})
                            entry['entries']+=1
                            entry['maximum_absolute_coefficient']=max(entry['maximum_absolute_coefficient'],measured)
    if args.baseline:
        old=[]
        old_symbols={}
        for line in args.baseline.open():
            tag,separator,value=line.partition(': ')
            if separator and tag.startswith('PY_S11CD_FULL_SECTOR_PENCIL_SYMBOL_'):
                old_symbols[tag]=value.strip()
            if separator and tag.startswith('PY_S11CD_FULL_PENCIL_MODE_'):
                old.extend(spectrum(association(record),tag,i) for i,record in enumerate(_restore(value)))
        result['baseline_spectrum_inventory']={
            'transcript':str(args.baseline.resolve()),'baseline_records':len(old),
            'current_records':len(spectra),'identical_literal_records':old==spectra,
            'different_record_indices':[i for i,(a,b) in enumerate(zip(old,spectra)) if a!=b]}
        result['baseline_symbol_inventory']={
            'baseline_symbols':len(old_symbols),'current_symbols':len(symbol_payloads),
            'identical_literal_payloads':old_symbols==symbol_payloads,
            'different_symbols':[tag for tag in sorted(old_symbols.keys()|symbol_payloads.keys())
                                 if old_symbols.get(tag)!=symbol_payloads.get(tag)]}
    missing.extend((tag,'fingerprint_metadata') for tag in fingerprints
                   if tag.replace('PY_S11CD_','PY_S11CD_METADATA_',1) not in tags)
    result.update(rectangular_jet_states=dict(states),nullspace_dimensions=dict(clusters),
                  rectangular_residuals_by_dimension=list(residuals.values()),
                  rectangular_coverage_gaps=missing,baseline_inventory_exit_code=baseline.returncode,
                  baseline_inventory_stderr=baseline.stderr,
                  fingerprint_residual_projections=fingerprint_residuals,
                  spectral_records=spectra)
    print(json.dumps(result,indent=2))
    raise SystemExit(bool(baseline.returncode or missing or pending))


if __name__=='__main__':
    run()
