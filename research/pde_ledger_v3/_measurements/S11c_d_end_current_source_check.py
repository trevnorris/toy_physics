#!/usr/bin/env python3
"""Pinned end-energy and closed-face source reconstruction before mode use."""
import argparse
import ast
from copy import copy
import json
from pathlib import Path
import pickle
import resource
import time

import sympy as sp

from S11c_d_joint_sheet_check import load, engine
from S11c_d_modal_current_check import digest, source_node


def main():
    parser=argparse.ArgumentParser()
    for name in ('manifest','input','run-directory'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--end',choices=('LEFT','RIGHT'),required=True)
    args=parser.parse_args()
    destination=args.run_directory
    destination.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    def progress(stage,**values):
        with (destination/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps({'stage':stage,'elapsedSeconds':time.monotonic()-started,**values})+'\n')
    modes,source_strong,units,branches,inputs,provenance=load(args)
    manifest=json.loads(args.manifest.read_text());base=Path(manifest['run_directory'])
    frozen=base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    if source_node(frozen.read_text(),'UniformSlabCurrent')!=source_node(Path(engine.__file__).read_text(),'UniformSlabCurrent'):
        raise ValueError('end energy constructor changed from native producer')
    energy=pickle.loads((base/'symbols'/f'{args.end}_LAB_HELD_RHO4_CONSTANT.pickle').read_bytes())[5]
    r=copy(modes.r)
    r.x=tuple(r.symbols['s11cc2X'+str(i)] for i in (1,2,3));r.t=r.symbols['s11cc2Time']
    r.z=r.symbols['s11cdNormalPosition']
    r.end_values={key:inputs.limits[value] for key,value in r.end_values.items()}
    endpoint=-sp.oo if args.end=='LEFT' else sp.oo
    strong=source_strong.xreplace(inputs.limits)
    ends=engine.ConstantEndPencil.__new__(engine.ConstantEndPencil);ends.r,ends.kn=r,modes.k
    modes=engine.FullPencilModes(ends,modes.curl,modes.units)
    current=engine.UniformSlabCurrent(r,{'value':energy},ends,strong[3,:])
    balance=engine.SlabEnergyBalance(current)
    acoustic=engine.ClosedAcousticEnergy(balance,modes,strong)
    provenance.update({'instrumentSha256':digest(Path(__file__)),
        'profileLimits':{str(key):str(value) for key,value in inputs.limits.items()},
        'endpoint':str(endpoint),'sourceEnergySha256':__import__('hashlib').sha256(sp.srepr(energy).encode()).hexdigest()})
    root=Path(engine.__file__).resolve().parents[1]
    for name in ('scripts/S11c_d_mixing_scattering_sympy_audit.py','_measurements/S11c_d_end_current_source_check.py'):
        target=destination/'source'/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((root/name).read_bytes())
    progress('loaded',provenance=provenance)
    source=current.construct('LAB_HELD',endpoint);progress('conservative_constructed')
    slab=balance.construct('LAB_HELD',endpoint);progress('balance_constructed')
    bulk=acoustic.construct('LAB_HELD',endpoint);progress('acoustic_constructed')
    payload={'conservative':source,'slab':slab,'acoustic':bulk,'strong':strong,
        'profileBindings':inputs.limits,'knownDimensions':engine.PHYSICAL_METADATA.dimensions.known}
    objects=destination/'objects.pickle';objects.write_bytes(pickle.dumps(payload,protocol=5))
    progress('objects_saved',sha256=digest(objects))
    prefix='END_CURRENT_SOURCE_'+args.end+'_LAB_HELD_RHO4_CONSTANT'
    def output(name,value,unit=lambda path:(0,0,0),heavy=False):
        body=engine.cas(value)
        engine.emit(prefix+'_'+name,engine.carrier_fingerprint(body) if heavy else body)
        engine.emit('METADATA_'+prefix+'_'+name,modes.numeric_metadata(body,unit))
    output('PROVENANCE',provenance)
    output('SOURCE_BOUND_PENCIL',strong,lambda path:units[path],True)
    current.emit('LAB_HELD',endpoint,prefix)
    balance.emit('LAB_HELD',endpoint,prefix)
    acoustic.emit('LAB_HELD',endpoint,prefix)
    progress('source_emitted')
    # Grade projection is performed on actual residuals with eta/sigma live.
    # The raw source residuals above remain in the transcript.
    checks={}
    for group,values in (('conservative',source),('slab',slab),('acoustic',bulk)):
        for key,value in values.items():
            if not key.endswith('_RESIDUAL'):continue
            entries=[v for _,v in engine.leaves(engine.cas(value)) if not isinstance(v,engine.Str)]
            projected=[current.retained(v) for v in entries]
            checks[group+'_'+key]={'scalars':len(entries),
                'rawNonzeroScalars':sum(v!=0 for v in entries),
                'retainedNonzeroScalars':sum(v!=0 for v in projected),
                'retainedGrades':sorted({grade for v in projected for grade in engine.PHYSICAL_METADATA.coefficients(v)})}
    output('RESIDUAL_GRADE_CENSUS',checks)
    output('DIMENSION_CONSTRAINTS',tuple(engine.PHYSICAL_METADATA.dimensions.constraints))
    index=engine.emission_index(engine.EMISSION_LINES)
    engine.emit(prefix+'_EMISSION_LINES',index)
    engine.emit('METADATA_'+prefix+'_EMISSION_LINES',modes.numeric_metadata(engine.cas(index),lambda path:(0,0,0)))
    summary={'provenance':provenance,'residualInventory':checks,'objectsSha256':digest(objects),
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'retainedNonzeroScalars':sum(v['retainedNonzeroScalars'] for v in checks.values())}
    (destination/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    progress('completed',summary=summary)
    if summary['retainedNonzeroScalars']:
        raise ValueError('end source retained-grade residual; inspect emitted operands')


if __name__=='__main__':main()
