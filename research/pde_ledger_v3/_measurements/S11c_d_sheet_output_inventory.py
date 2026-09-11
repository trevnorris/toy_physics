#!/usr/bin/env python3
"""Read tag/payload shapes, stored labels, source hashes, and file metadata."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from ledger_fold import _restore


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('transcript',type=Path)
    args=parser.parse_args()
    tags=[]
    modes=[]
    labels=Counter()
    dimensions={}
    completion=[]
    source_pins={}
    outstanding=[]
    symbols=[]
    execution={}
    missing_paths=[]
    for line in args.transcript.open():
        tag,separator,value=line.partition(': ')
        if not separator or not tag.startswith('PY_S11CD_'):
            continue
        tags.append(tag)
        if tag=='PY_S11CD_BUILD_INPUT_DIGESTS':
            source_pins={str(k):str(v) for k,v in _restore(value)}
        if tag.startswith('PY_S11CD_FULL_SECTOR_PENCIL_SYMBOL_'):
            symbols.append(tag)
        if tag.startswith('PY_S11CD_FULL_PENCIL_MODE_'):
            records=_restore(value)
            counts=Counter()
            path_count=0
            for index,record in enumerate(records):
                slots={str(k):v for k,v in record}
                label=str(slots.get('PHYSICAL_BULK_SHEET','UNCOMPUTED'))
                counts[label]+=1
                labels[label]+=1
                if 'BULK_SHEET_PATH' in slots:
                    path_count+=1
                elif str(slots.get('FINITE_PENCIL'))!='False':
                    missing_paths.append({'tag':tag,'candidate_index':index})
            modes.append({'tag':tag,'candidates':len(records),'labels':dict(counts),
                          'branch_path_records':path_count})
        if tag in ('PY_S11CD_REDUCED_DIMENSION_CONSTRAINT_RESIDUALS',
                   'PY_S11CD_REDUCED_DIMENSION_UNRESOLVED',
                   'PY_S11CD_PENCIL_DIMENSION_CONSTRAINT_RESIDUALS'):
            dimensions[tag]=value.strip()
        if tag=='PY_S11CD_OUTSTANDING_CONSTRUCTIONS':
            outstanding=list(map(str,_restore(value)))
        if tag=='PY_S11CD_PROCESS_COMPLETION':
            completion.append(value.strip())
        if tag in ('PY_S11CD_IMPLEMENTATION_CHECKPOINT','PY_S11CD_CHANNEL_INPUT_EXECUTION',
                   'PY_S11CD_RESOURCE_MEASUREMENTS'):
            execution[tag]=value.strip()
    counts=Counter(tags)
    cases=[a+'_'+r for a in ('LAB_HELD','MATERIAL_ADVECTED')
           for r in ('RHO4_CONSTANT','RHOBR_CONSTANT')]
    expected_symbols={'PY_S11CD_FULL_SECTOR_PENCIL_SYMBOL_'+end+'_'+case
                      for case in cases for end in ('REFERENCE','LEFT','RIGHT')}
    expected_modes={'PY_S11CD_FULL_PENCIL_MODE_PIT_'+end+'_'+case+'_'+str(sample)
                    for case in cases for end in ('LEFT','RIGHT') for sample in range(3)}
    actual_modes={p['tag'] for p in modes}
    coverage={'missing_symbols':sorted(expected_symbols-set(symbols)),
              'extra_symbols':sorted(set(symbols)-expected_symbols),
              'missing_mode_packets':sorted(expected_modes-actual_modes),
              'extra_mode_packets':sorted(actual_modes-expected_modes),
              'missing_mode_metadata':sorted('PY_S11CD_METADATA_'+tag[len('PY_S11CD_'):]
                  for tag in actual_modes if 'PY_S11CD_METADATA_'+tag[len('PY_S11CD_'):] not in counts),
              'missing_branch_paths':missing_paths}
    changed=[path for path,h in source_pins.items() if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=h]
    result={'transcript':str(args.transcript.resolve()),'bytes':args.transcript.stat().st_size,
            'sha256':hashlib.sha256(args.transcript.read_bytes()).hexdigest(),
            'tag_count':len(tags),'duplicate_tags':{k:v for k,v in counts.items() if v>1},
            'source_pins':source_pins,'changed_source_pins':changed,
            'full_symbol_tags':symbols,'mode_packets':modes,'sheet_label_counts':dict(labels),
            'dimension_records':dimensions,'outstanding_constructions':outstanding,
            'completion_records':completion,'execution_records':execution,'coverage_gaps':coverage}
    print(json.dumps(result,indent=2))
    raise SystemExit(bool(changed or result['duplicate_tags'] or len(completion)!=1 or any(coverage.values())))


if __name__=='__main__':
    run()
