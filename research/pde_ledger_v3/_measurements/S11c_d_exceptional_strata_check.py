#!/usr/bin/env python3
"""Pinned one-case exceptional-slice and native-regularity construction."""
import argparse
import hashlib
import inspect
from pathlib import Path
import pickle
import resource
import time
from types import SimpleNamespace

from S11c_d_joint_sheet_check import load,engine


def run():
    started=time.monotonic()
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--end',choices=('REFERENCE','LEFT','RIGHT'),default='REFERENCE')
    parser.add_argument('--native',action='store_true')
    parser.add_argument('--end-only',action='store_true')
    parser.add_argument('--threshold-only',action='store_true')
    parser.add_argument('--analysis-cache',type=Path)
    args=parser.parse_args()
    modes,strong,units,bindings,inputs,provenance=load(args)
    sp=engine.sp;dims=engine.PHYSICAL_METADATA.dimensions;r=modes.r
    base=Path(engine.json.loads(args.manifest.read_text())['run_directory'])
    full=pickle.loads((base/'symbols'/f'{args.end}_LAB_HELD_RHO4_CONSTANT.pickle').read_bytes())[0]
    # Re-enter the existing field ansatz to construct the physical lift.
    r.z=sp.Symbol('s11cdExceptionalCheckNormalCoordinate',real=True)
    pencil=engine.ReducedPencil.__new__(engine.ReducedPencil);pencil.r=r
    pencil.fields=tuple(sp.Function('s11cdReducedField'+name) for name in ('u1','u2','u3','theta','eW'))
    coordinates=sp.symbols('s11cdExceptionalCheckTangent0:2',real=True)
    plane=sp.exp(sp.I*sum(k*x for k,x in zip(r.tangents,coordinates)))
    pencil.tangent_derivatives=tuple(sp.diff(plane,x)/plane for x in coordinates)
    ends=engine.ConstantEndPencil.__new__(engine.ConstantEndPencil);ends.r=r;ends.kn=modes.k
    lift=ends.field_lift(pencil)
    spectrum=engine.EndSpectrumCoverage(modes,strong,full,lift,units)
    if args.analysis_cache:
        args.analysis_cache.mkdir(parents=True,exist_ok=True)
        def cached_method(method,dependencies=()):
            source=inspect.getsource(method)+''.join(inspect.getsource(f) for f in dependencies)
            def cached(*operands):
                key=hashlib.sha256((source+sp.srepr(engine.cas(operands))).encode()).hexdigest()
                path=args.analysis_cache/(key+'.pickle')
                if path.exists():return pickle.loads(path.read_bytes())
                value=method(*operands)
                path.write_bytes(pickle.dumps(value))
                return value
            return cached
        for name in ('analyze','threshold_data'):
            setattr(engine.EndExceptionalSlice,name,staticmethod(cached_method(getattr(engine.EndExceptionalSlice,name))))
        engine.BulkExceptionalSlice.analyze=staticmethod(cached_method(engine.BulkExceptionalSlice.analyze,
            (engine.BulkExceptionalSlice.real_normal_projection,)))
    provenance['exceptionalInstrumentSha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    provenance['codecSha256']=hashlib.sha256((engine.ROOT/'scripts/S11c_d_output_codec.py').read_bytes()).hexdigest()
    def output(tag,value,unit=None):
        spectrum.emit('EXCEPTIONAL_PREFLIGHT_'+tag,value,unit)
    output('PROVENANCE',provenance)
    if args.native:
        spectrum.construct(args.end+'_LAB_HELD_RHO4_CONSTANT',-1 if args.end=='LEFT' else 1,
            channel_input=inputs,reference=args.end=='REFERENCE')
    else:
        data=engine.EndExceptionalSlice(spectrum).construct(args.end+'_LAB_HELD_RHO4_CONSTANT',inputs,reference=args.end=='REFERENCE')
        if args.threshold_only:
            engine.ThresholdModeAudit(spectrum,bindings).construct(args.end+'_LAB_HELD_RHO4_CONSTANT',inputs,
                reference=args.end=='REFERENCE',end_data=data)
        elif not args.end_only:
            engine.BulkExceptionalSlice(spectrum,bindings).construct(args.end+'_LAB_HELD_RHO4_CONSTANT',inputs,
                reference=args.end=='REFERENCE',end_data=data)
    output('DIMENSION_CONSTRAINTS',tuple(dims.constraints))
    output('RESOURCES',{'WALL_SECONDS':time.monotonic()-started,'PEAK_RSS_KIB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
           lambda p:(0,1,0) if p[-1]=='WALL_SECONDS' else dims.zero)
    engine.emit('EMISSION_LINES',engine.emission_index(engine.EMISSION_LINES))
    output('PROCESS_COMPLETION',True)


if __name__=='__main__':run()
