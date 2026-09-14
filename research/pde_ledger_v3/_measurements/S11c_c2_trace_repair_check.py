#!/usr/bin/env python3
"""Actual inherited trace, retained kernel solve, and source-row factorization."""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import sys
import time
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from ledger_fold import load_model,check_consumer,assert_lookups_equal_manifest
import S11c_c2_selfenergy_fold_sympy_audit as c
import S11c_d_mixing_scattering_sympy_audit as d


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class TraceMetadata(d.PhysicalMetadata):
    @lru_cache(maxsize=None)
    def coefficients(self,value):
        # Exact Laurent monomials account for the inherited wave-amplitude
        # normalization. No Taylor replacement of a rational sum is allowed.
        if value.is_Pow and value.base in self.generators and value.exp.is_Integer:
            grade=[0]*len(self.generators);grade[self.generators.index(value.base)]=int(value.exp)
            return {tuple(grade):sp.S.One}
        return super().coefficients(value)


def canonical_power(expression,progress=lambda *a,**k:None):
    carriers=sorted(expression.atoms(sp.Integral,sp.Derivative),key=sp.default_sort_key)
    mapping={v:sp.Dummy('traceRepairPowerCarrier') for v in carriers}
    shielded=expression.xreplace(mapping)
    active=tuple(v for v in mapping.values() if v in shielded.free_symbols)
    if not active:return sp.cancel(shielded)
    terms=sp.Poly(shielded,*active).terms();result=sp.S.Zero
    for index,(powers,coefficient) in enumerate(terms):
        progress('power_coefficient',index=index+1,count=len(terms))
        result+=sp.cancel(coefficient)*sp.Mul(*(v**p for v,p in zip(active,powers)))
    return result.xreplace({v:k for k,v in mapping.items()})


def run():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--inputs-cache',type=Path);parser.add_argument('--all-cases',action='store_true')
    parser.add_argument('--build-case',action='store_true')
    parser.add_argument('--model-cache',type=Path);parser.add_argument('--power-cache',type=Path)
    args=parser.parse_args()
    if (args.model_cache or args.power_cache) and (not args.build_case or args.all_cases):
        raise ValueError('one-case build required for model/power replay')
    if args.power_cache and not args.model_cache:raise ValueError('power replay requires pinned model')
    base=args.run_directory;base.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    paths=[Path(__file__),ROOT/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py',ROOT/'scripts/S11c_d_mixing_scattering_sympy_audit.py',
        ROOT/'scripts/S11c_b_exports.py',ROOT/'scripts/S11c_c1_exports.py',ROOT/'scripts/ledger_fold.py',
        ROOT/'directives/S11c_c2_SHARED_PHYSICS.md',ROOT/'directives/S11c_c2_sympy_build_directive.md']
    pins={str(p.relative_to(ROOT)):digest(p) for p in paths}
    for p in paths:
        target=base/'source'/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    def progress(stage,**extra):
        with (base/'progress.jsonl').open('a') as f:f.write(json.dumps({'stage':stage,'seconds':time.monotonic()-started,**extra})+'\n')
    progress('load')
    sourcepins={name:pin for name,pin in pins.items() if name in ('scripts/S11c_b_exports.py','scripts/S11c_c1_exports.py','scripts/ledger_fold.py')}
    if args.inputs_cache:
        meta=json.loads(args.inputs_cache.with_suffix('.json').read_text())
        if meta['sourcePins']!=sourcepins or digest(args.inputs_cache)!=meta['sha256']:raise ValueError('input cache pins')
        values,audit=pickle.loads(args.inputs_cache.read_bytes());inputs=c.Inputs(values)
    else:
        fold,audit=load_model(ROOT/'scripts/S11c_b_exports.py',ROOT/'scripts/S11c_c1_exports.py')
        closure=check_consumer(fold,c.IMPORT_KEYS);lookup=assert_lookups_equal_manifest(c.bind_inputs,fold,c.IMPORT_KEYS)
        inputs=lookup['result'];audit|={'lookups':sorted(lookup['lookups']),'closure':sorted(closure['closure'])}
        cache=base/'inputs.pickle';cache.write_bytes(pickle.dumps((inputs.values,audit),protocol=5))
        cache.with_suffix('.json').write_text(json.dumps({'sourcePins':sourcepins,'sha256':digest(cache)},indent=2)+'\n')
    progress('inputs')
    records=[];writekeys=set();total=nonzero=0;results={}
    metadata=TraceMetadata.__new__(TraceMetadata);metadata.generators=(inputs.eps,inputs.eta,inputs.sigma)
    adopted={}
    if args.model_cache:
        old=args.model_cache.parent/'source'
        for name in ('scripts/S11c_c2_selfenergy_fold_sympy_audit.py','scripts/S11c_b_exports.py',
                     'scripts/S11c_c1_exports.py','scripts/ledger_fold.py'):
            if digest(old/name)!=pins[name]:raise ValueError(('model cache source',name))
        adopted['model']={'path':str(args.model_cache),'sha256':digest(args.model_cache)}
    def output(case,face,name,value,units=None,residual=False):
        nonlocal total,nonzero
        value=c.cas(value);units=c.cas(c.dimensions(value) if units is None else units)
        if units.has(sp.nan,sp.zoo):raise ValueError(('unresolved units',case,face,name))
        flatunits=d.payload_units(value,units)
        for path,leaf in d.leaves(value):
            if isinstance(leaf,c.Str):continue
            if leaf!=0 and sp.ImmutableMatrix(c.dimension(leaf))!=sp.ImmutableMatrix(flatunits[path]):
                raise ValueError(('dimension',case,face,name,path))
        key='s11cc2TraceRepair'+''.join(v.title().replace('_','') for v in (*case,'Plus' if face==1 else 'Minus',name))
        if key in writekeys or key in inputs.values:raise ValueError('trace write key')
        writekeys.add(key)
        grade_support=set();lambda_support=set()
        for _,leaf in d.leaves(value):
            if isinstance(leaf,c.Str):continue
            coefficients=metadata.coefficients(leaf.xreplace(inputs.profiles));grade_support.update(coefficients)
            homotopy_coefficients={}
            for (e,a,b),coefficient in coefficients.items():
                grade=(e,a+b)
                homotopy_coefficients[grade]=homotopy_coefficients.get(grade,sp.S.Zero)+coefficient*(inputs.values['W_0']/inputs.values['L_W'])**b
            lambda_support.update(g for g,v in homotopy_coefficients.items() if v!=0)
        grades=sp.Tuple(*(sp.Tuple(*g) for g in sorted(grade_support)))
        homotopy=sp.Tuple(*(sp.Tuple(*g) for g in sorted(lambda_support)))
        body={'writeKey':key,'value':d.carrier_fingerprint(value) if d.dag_size(value)>1200 else value,
              'multigrade':grades,'epsilonLambdaSupport':homotopy,'dimensionLTM':units}
        print(sp.srepr(c.cas(body)),flush=True)
        records.append({'case':case,'face':face,'name':name,'value':value,'units':units,'payload':body})
        if residual:
            scalars=[v for _,v in d.leaves(value) if not isinstance(v,c.Str)]
            total+=len(scalars);nonzero+=sum(v!=0 for v in scalars)
    cases=[(a,rho) for a in c.ANCHORINGS for rho in c.DENSITIES] if args.all_cases else [('LAB_HELD','RHO4_CONSTANT')]
    for case in cases:
        anchoring,density=case;results[case]={}
        for face in c.FACES:
            progress('face_started',case=case,face=face)
            response,z,resolvent,zm,inv_operand,inv,pmat,ko,ki,qo,second=c.kernel_bridge(inputs,*case[:1],face,density,{})
            reference,reference_second,trace=c.reference_pressure_kernels(inputs,anchoring,face,density,pmat,second,ko,ki,qo,{})
            traceunits=c.reference_trace_dimensions(inputs,trace)
            for key,value in trace.items():output(case,face,key,value,traceunits[key],key.endswith('_RESIDUAL'))
            label='plus' if face==1 else 'minus';p=inputs.a('delta_p_'+label);j=inputs.a('d_w_delta_p_'+label)
            row=c.expanded_rows(inputs.slab[case])
            expressions=[*row['U'],row['THETA'],row['E_W'],
                *inputs.geometry['traction'][(*case[:1],face,c.REPRESENTATION,density)],
                inputs.geometry['closure_shape_deriv'][(*case[:1],face,c.REPRESENTATION,density)]]
            joins=[]
            for index,source in enumerate(expressions):
                value_coef=sp.diff(source,p);jet_coef=sp.diff(source,j)
                factorized=trace['NORMAL_JET_COEFFICIENT']*value_coef/trace['VALUE_COEFFICIENT']
                difference=sp.cancel(c.retained_shape(jet_coef-factorized,inputs))
                unit=c.dimension(source);jetunit=c.dimension(j)
                unit=sp.ImmutableMatrix(tuple(a-b for a,b in zip(unit,jetunit)))
                output(case,face,'SourceJetCoefficient'+str(index),jet_coef,unit)
                output(case,face,'TraceFactorizedCoefficient'+str(index),factorized,unit)
                output(case,face,'SourceTraceResidual'+str(index),difference,unit,True)
                joins.append((jet_coef,factorized,difference))
            mixed=sp.diff(reference_second,inputs.eta,inputs.sigma).subs({inputs.eta:0,inputs.sigma:0})
            output(case,face,'ReferenceMixedKernel',mixed)
            results[case][face]={'trace':trace,'sourceJoins':joins,'reference':reference,'referenceSecond':reference_second}
            progress('face_completed',case=case,face=face)
        if args.build_case:
            progress('case_build_started',case=case)
            model=pickle.loads(args.model_cache.read_bytes()) if args.model_cache else c.build_case(inputs,*case)
            target=base/('model_'+'_'.join(case)+'.pickle');target.write_bytes(pickle.dumps(model,protocol=5))
            progress('case_model_saved',case=case,modelSha256=digest(target))
            if args.power_cache:
                if args.power_cache.parent!=args.model_cache.parent:raise ValueError('power/model cache origin')
                cache_record=json.loads(args.power_cache.with_suffix('.json').read_text())
                if (cache_record['modelSha256']!=digest(args.model_cache) or
                    cache_record['sha256']!=digest(args.power_cache) or
                    cache_record['nativeSourceSha256']!=pins['scripts/S11c_c2_selfenergy_fold_sympy_audit.py']):
                    raise ValueError('power cache provenance')
                saved=pickle.loads(args.power_cache.read_bytes())
                c.NEW_DIMENSIONS.update(saved['dimensions']);c.DIMENSION_SCHEMA.update(saved['schema']);c.dimension.cache_clear()
                covectors,pairing,residual=saved['pairing']
                adopted['power']={'path':str(args.power_cache),'sha256':digest(args.power_cache)}
            else:
                covectors,pairing,residual=c.traction_pairing(inputs,case,model)
            if args.model_cache:
                # Re-enter the native declaration sites skipped by loading a
                # model pickle. Units come from the unchanged source bindings.
                c.outgoing_spectral(inputs,tuple(ko),tuple(ki))
                for fmap in model['maps'].values():
                    for _,source in fmap['IDENTIFICATIONS']:
                        inputs.at_source(source.xreplace(inputs.profiles))
                replay_open=c.tree(c.extract(model['open'],inputs),lambda e:c.retained_shape(e,inputs))
                # Compare only after printing below; this also restores the
                # native trial/test declarations used by the closed kernel.
                open_replay_residual=c.difference(replay_open,model['open_kernel'])
                c.dimension.cache_clear()
            power_target=base/'power.pickle'
            power_target.write_bytes(pickle.dumps({'pairing':(covectors,pairing,residual),
                'dimensions':c.NEW_DIMENSIONS,'schema':c.DIMENSION_SCHEMA},protocol=5))
            power_target.with_suffix('.json').write_text(json.dumps({'modelSha256':digest(target),
                'sha256':digest(power_target),'nativeSourceSha256':pins['scripts/S11c_c2_selfenergy_fold_sympy_audit.py']},indent=2)+'\n')
            progress('power_saved',sha256=digest(power_target))
            output(case,1,'ClosedSlabOperator',model['closed'])
            output(case,1,'ClosedCouplingKernel',model['closed_kernel'])
            if args.model_cache:
                output(case,1,'OpenKernelReplayOperand',replay_open)
                output(case,1,'OpenKernelCacheOperand',model['open_kernel'])
                output(case,1,'OpenKernelReplayCanonicalResidual',open_replay_residual,c.dimensions(replay_open),True)
            powers=sorted({v for v in c.cas((covectors,pairing,residual)).atoms(sp.Pow)
                if v.exp.is_negative and v.has(*metadata.generators)},key=sp.default_sort_key)
            output(case,1,'AmplitudeNormalizationReciprocals',sp.Tuple(*powers))
            output(case,1,'AmplitudeNormalizationExcludedZeroBases',sp.Tuple(*(v.base for v in powers)))
            output(case,1,'TractionPowerOperands',pairing)
            output(case,1,'TractionPowerResidual',residual,c.dimensions(pairing['SLAB_POWER']))
            normalized=canonical_power(residual,progress)
            output(case,1,'TractionPowerCanonicalResidual',normalized,c.dimensions(pairing['SLAB_POWER']),True)
            results[case]['power']={'raw':residual,'canonical':normalized}
            progress('case_build_completed',case=case,modelSha256=digest(target))
    payload={'records':records,'results':results,'metadataGenerators':metadata.generators,'profiles':inputs.profiles,'homotopyRatio':inputs.values['W_0']/inputs.values['L_W'],'dimensions':c.NEW_DIMENSIONS,'schema':c.DIMENSION_SCHEMA}
    (base/'objects.pickle').write_bytes(pickle.dumps(payload,protocol=5))
    summary={'sourcePins':pins,'sourcePinsAfter':{str(p.relative_to(ROOT)):digest(p) for p in paths},'inputAudit':audit,'adoptedCaches':adopted,
        'objectsSha256':digest(base/'objects.pickle'),'objects':len(records),'residualScalars':total,'nonzeroResidualScalars':nonzero,
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (base/'checks.json').write_text(json.dumps(summary,indent=2)+'\n');progress('completed')
    if nonzero or summary['sourcePins']!=summary['sourcePinsAfter']:raise ValueError('emitted trace or source reconstruction discrepancy')


if __name__=='__main__':run()
