#!/usr/bin/env python3
"""Source-pinned constant-end resolvent and Laurent-contour construction."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import time

from S11c_d_joint_sheet_check import load,engine,_restore
from S11c_d_output_codec import decoded_lines


def run():
    started=time.monotonic()
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--end',choices=('REFERENCE','LEFT','RIGHT'),default='REFERENCE')
    parser.add_argument('--pit',action='store_true')
    args=parser.parse_args()
    modes,strong,units,bindings,inputs,provenance=load(args)
    provenance['resolventInstrumentSha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dims=engine.PHYSICAL_METADATA.dimensions
    engine.emit('END_RESOLVENT_PREFLIGHT_PROVENANCE',provenance)
    engine.emit('METADATA_END_RESOLVENT_PREFLIGHT_PROVENANCE',modes.numeric_metadata(engine.cas(provenance),lambda p:dims.zero))
    manifest=json.loads(args.manifest.read_text())
    prefix='PY_S11CD_END_SPECTRUM_'+('PIT_' if args.pit else 'INPUT_')+args.end+'_LAB_HELD_RHO4_CONSTANT_0_MODE_'
    records=[]
    for line in decoded_lines(Path(manifest['run_directory'])/'full.out'):
        tag,_,payload=line.partition(': ')
        if tag.startswith(prefix) and tag.endswith('_RECORD'):
            records.append({str(k):v for k,v in _restore(payload)})
    if not records:raise ValueError('missing native spectrum in pinned producer')
    audit=engine.EndResolventAudit(modes,strong,units,bindings)
    audit.construct(args.end+'_LAB_HELD_RHO4_CONSTANT',channel_input=None if args.pit else inputs,
                    reference=args.end=='REFERENCE',spectrum={'RECORDS':records})
    constraints=tuple(dims.constraints)
    engine.emit('END_RESOLVENT_PREFLIGHT_DIMENSION_CONSTRAINTS',constraints)
    engine.emit('METADATA_END_RESOLVENT_PREFLIGHT_DIMENSION_CONSTRAINTS',modes.numeric_metadata(engine.cas(constraints),lambda p:dims.zero))
    resources={'WALL_SECONDS':time.monotonic()-started,'PEAK_RSS_KIB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    engine.emit('END_RESOLVENT_PREFLIGHT_RESOURCES',resources)
    engine.emit('METADATA_END_RESOLVENT_PREFLIGHT_RESOURCES',modes.numeric_metadata(engine.cas(resources),lambda p:(0,1,0)if p[-1]=='WALL_SECONDS'else dims.zero))


if __name__=='__main__':run()
