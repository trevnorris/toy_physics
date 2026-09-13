#!/usr/bin/env python3
"""Source-pinned chemical/density and inherited pressure-jet provenance operands."""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import re
import resource
import shutil
import tempfile
import time
from types import SimpleNamespace
import sympy as sp
from S11c_d_joint_sheet_check import engine, digest, decoded_lines, _restore
from S11c_inertia_artifact_audit import export_data
from S11c_d_output_codec import restore_emission_index

ROOT=Path(__file__).resolve().parents[1]


def named(body,key):
    return next(v for k,v in body if str(k)==key)


def cases(body):
    return {tuple(map(str,k)):named(v,'VALUE') for k,v in body}


def metadata_context(data):
    known=dict(data['knownDimensions'])
    symbols={s.name:s for s in known if isinstance(s,sp.Symbol)}
    symbols.update({s.name:s for s in data['strong'].free_symbols})
    dimensions=engine.DimensionAnalysis.__new__(engine.DimensionAnalysis)
    dimensions.known,dimensions.unknown,dimensions.constraints,dimensions.solution=known,{},set(),{}
    dimensions.zero=(sp.S.Zero,)*3
    r=SimpleNamespace(symbols=symbols,ell=symbols['L_W'])
    engine.PHYSICAL_METADATA=engine.PhysicalMetadata(dimensions,r)
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r
    return symbols,dimensions,modes


def compute(args):
    started=time.monotonic();base=args.run_directory
    base.mkdir(parents=True,exist_ok=False)
    paths=[Path(__file__),ROOT/'scripts/S11c_d_mixing_scattering_sympy_audit.py',
        ROOT/'scripts/S11c_b_exports.py',ROOT/'scripts/S11c_c1_exports.py',ROOT/'scripts/S11c_c2_exports.py',
        ROOT/'scripts/S11c_a_interface_geometry_sympy_audit.py',ROOT/'scripts/S11c_b_brane_operator_sympy_audit.py',
        ROOT/'scripts/S11c_c1_bulk_closure_sympy_audit.py',ROOT/'scripts/S11c_c2_selfenergy_fold_sympy_audit.py',
        ROOT/'_measurements/S11c_d_joint_sheet_check.py',ROOT/'_measurements/S11c_inertia_artifact_audit.py',
        ROOT/'scripts/ledger_fold.py',ROOT/'scripts/S11c_d_output_codec.py',
        ROOT/'directives/S11c_d_SHARED_PHYSICS.md',ROOT/'directives/S11c_c1_SHARED_PHYSICS.md']
    pins={str(p):digest(p) for p in paths}
    for p in paths:
        target=base/'source'/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    values,_,_=export_data(ROOT/'scripts/S11c_b_exports.py')
    inputs={k:_restore(values[k]) for k in ('mu_theta_operator','slab_operator','background_density_map')}
    case=('LAB_HELD','RHO4_CONSTANT')
    chemical=cases(inputs['mu_theta_operator'])[case][1]
    operator=cases(inputs['slab_operator'])[case]
    density=cases(inputs['background_density_map'])[('RHO4_CONSTANT',)][1]
    rows={key:named(named(operator,slot),'EXPANDED') for key,slot in
          [('MASS','THETA_BALANCE'),('MECHANICAL','E_W_BALANCE')]}
    packets=[];results={};run_pins={};proof_scalars=0;nonzero_proof=0
    for end,run in [('LEFT',args.left_run),('RIGHT',args.right_run)]:
        checks=json.loads((run/'checks.json').read_text())
        if digest(run/'objects.pickle')!=checks['objectsSha256']:raise ValueError('saved end object pin')
        if digest(ROOT/'scripts/S11c_d_mixing_scattering_sympy_audit.py')!=checks['provenance']['engineSha256']:
            raise ValueError('native source changed from completed end packet')
        producer=checks['provenance']['producerSources']
        for name in ('scripts/S11c_b_exports.py','scripts/S11c_c1_exports.py','scripts/S11c_c2_exports.py'):
            if digest(ROOT/name)!=producer[name]:raise ValueError(('producer input pin',name))
        run_pins[end]={str(run/name):digest(run/name) for name in ('objects.pickle','checks.json','full.out')}
        data=pickle.loads((run/'objects.pickle').read_bytes())
        symbols,dimensions,modes=metadata_context(data)
        eps,eta,sigma=(symbols[n] for n in ('epsilon_shape','eta_bg','sigma_W'))
        amps=[symbols['s11cdCurrentPlusAmplitude'+str(i)] for i in range(5)]
        field_units=[dimensions.measure(a) for a in amps]
        x=tuple(symbols['s11cc2X'+str(i)] for i in (1,2))+(symbols['s11cdNormalPosition'],)
        k=tuple(symbols['s11cdTangentialMomentum'+str(i)] for i in (1,2))+(symbols['s11cdSpectralNormalMomentum'],)
        t=symbols['s11cc2Time'];phase=sp.exp(sp.I*(sum(a*b for a,b in zip(k,x))-symbols['omega']*t))
        fields=[a*phase for a in amps]
        def end_map(value,waves=False):
            mapping={}
            for atom in value.atoms(sp.Symbol):
                name=re.sub(r'^grad_theta_([123])$',r'theta_d\1',atom.name)
                match=re.fullmatch(r'(u_[123]|theta|e_W)((?:_t{1,2})?(?:_?d[123])*)',name)
                if match and waves:
                    field,suffix=match.groups();v=fields[{'u_1':0,'u_2':1,'u_3':2,'theta':3,'e_W':4}[field]]
                    for i in re.findall(r'd([123])',suffix):v=sp.diff(v,x[int(i)-1])
                    if '_tt' in suffix:v=sp.diff(v,t,2)
                    elif '_t' in suffix:v=sp.diff(v,t)
                    mapping[atom]=v
                elif re.fullmatch(r'([wm])1_profile((?:_d[123](?:d[123])*)?)',atom.name):
                    profile='w' if atom.name.startswith('w') else 'm'
                    # Apply the already supplied, pinned profile limit; derivative jets vanish on this constant end.
                    limit=next(v for key,v in data['profileBindings'].items() if
                        isinstance(key,sp.Limit) and str(key.args[0].func)=='s11cd'+profile.upper()+'Profile'
                        and key.args[2]==(-sp.oo if end=='LEFT' else sp.oo))
                    mapping[atom]=sp.diff(limit,symbols['s11cdProfileCoordinate'],len(re.findall(r'd[123]',atom.name)))
                elif atom.name in symbols:mapping[atom]=symbols[atom.name]
            return value.xreplace(mapping),mapping
        mapped,mapping=end_map(chemical,True);mapped=sp.expand(mapped/eps/phase)
        row=sp.ImmutableMatrix(1,5,[sp.diff(mapped,a) for a in amps])
        energy_row=data['slab']['CHEMICAL_FIELD_ROW']
        chemical_difference=(row-energy_row).applyfunc(sp.cancel)
        density_end,density_map=end_map(density)
        density_source=data['acoustic']['SLAB_SURFACE_DENSITY']
        objects={'IMPORTED_CHEMICAL_ROW':row,'ENERGY_CHEMICAL_ROW':energy_row,
                 'CHEMICAL_JOIN_RESIDUAL':chemical_difference,'IMPORTED_SURFACE_DENSITY':density_end,
                 'MASS_RATE_SURFACE_DENSITY':density_source,'DENSITY_JOIN_RESIDUAL':sp.cancel(density_end-density_source)}
        chemical_units={(i,):dimensions.measure(v) for i,v in enumerate(energy_row)}
        unit_maps={name:chemical_units for name in objects if 'CHEMICAL' in name}
        unit_maps.update({name:{():dimensions.measure(density_source)} for name in objects if 'DENSITY' in name})
        proof_names=['CHEMICAL_JOIN_RESIDUAL','DENSITY_JOIN_RESIDUAL'];jet_operands={}
        for name,source_row in rows.items():
            terms=[];jet_operands[name]=[]
            for face in data['acoustic']['FACE_RECORDS']:
                sign=int(face['ORIENTATION']);label='plus' if sign==1 else 'minus'
                slot=next(s for s in source_row.free_symbols if s.name=='d_w_delta_p_'+label)
                source_eps=next(s for s in source_row.free_symbols if s.name=='epsilon_shape')
                coefficient,coefficient_map=end_map(sp.diff(source_row,slot)/source_eps)
                coefficient=sp.cancel(coefficient)
                pressure=face['PRESSURE'].subs({eta:0,sigma:0})
                # Outgoing bulk continuation ansatz, with the lab-normal derivative computed before taking its trace.
                normal=sp.Symbol('s11cdTraceDiagnosticBulkPosition',real=True)
                dimensions.known[normal]=dimensions.measure(x[2])
                continuation=pressure*sp.exp(sign*sp.I*symbols['s11cdAcousticRightNormalMomentum']*normal)
                jet=sp.diff(continuation,normal).subs(normal,0)
                terms.append(coefficient*jet)
                jet_operands[name].append({'sourceCoefficient':sp.diff(source_row,slot)/source_eps,
                    'alignment':coefficient_map,'coefficient':coefficient,'pressure':pressure,'continuation':continuation,'jet':jet})
                for suffix,value in [('ROW_JET_COEFFICIENT',coefficient),('REFERENCE_PRESSURE',pressure),('REFERENCE_NORMAL_JET',jet)]:
                    key=name+'_'+label.upper()+'_'+suffix;objects[key]=value
                    if value!=0:unit=dimensions.measure(value)
                    else:
                        row_index=3 if name=='MASS' else 4
                        j=next(j for j in range(5) if data['strong'][row_index,j]!=0)
                        row_unit=tuple(a+b for a,b in zip(dimensions.measure(data['strong'][row_index,j]),field_units[j]))
                        unit=tuple(a-b for a,b in zip(row_unit,dimensions.measure(jet)))
                    unit_maps[key]={():unit}
            radical_map={symbols['s11cdAcousticRightNormalMomentum']:
                         data['acoustic']['ACOUSTIC_RADICAL_SCALE']*symbols['s11cdBulkRadical']}
            contribution=sp.Add(*terms).subs(radical_map)
            jet_row=sp.ImmutableMatrix(1,5,[sp.cancel(sp.diff(contribution,a)) for a in amps])
            discrepancy=data['acoustic']['REDUCED_'+name+'_FACE_JOIN_RESIDUAL']
            difference=(discrepancy-jet_row).applyfunc(lambda v:sp.cancel(sp.together(v)))
            index=3 if name=='MASS' else 4
            j=next(j for j in range(5) if data['strong'][index,j]!=0)
            row_unit=tuple(a+b for a,b in zip(dimensions.measure(data['strong'][index,j]),field_units[j]))
            units={(i,):tuple(a-b for a,b in zip(row_unit,field_units[i])) for i in range(5)}
            for suffix,value in [('SAVED_SOURCE_DISCREPANCY',discrepancy),('INHERITED_JET_CONTRIBUTION',jet_row),('JET_ACCOUNTING_RESIDUAL',difference)]:
                key=name+'_'+suffix;objects[key]=value;unit_maps[key]=units
            proof_names.append(name+'_JET_ACCOUNTING_RESIDUAL')
        denominators=tuple(sorted({p.base for value in objects.values() for p in value.atoms(sp.Pow)
                                   if p.exp.is_negative},key=sp.default_sort_key))
        objects['RECIPROCAL_DOMAIN_OPERANDS']=denominators
        unit_maps['RECIPROCAL_DOMAIN_OPERANDS']={(i,):dimensions.measure(v) for i,v in enumerate(denominators)}
        for name,value in objects.items():
            body=engine.cas(value);units=unit_maps[name]
            for path,leaf in engine.leaves(body):
                if leaf!=0 and dimensions.measure(leaf)!=units[path]:raise ValueError(('source dimension',end,name,path))
            metadata=modes.numeric_metadata(body,lambda path:units[path])
            tag='PY_S11CD_RIGHT_SOURCE_TRACE_'+end+'_'+name
            heavy=engine.dag_size(body)>1200
            write_key='s11cdRightSourceTrace'+''.join(word.title() for word in (end+'_'+name).split('_'))
            if write_key in values or write_key in engine.IMPORT_KEYS:raise ValueError('imported write key')
            packets.append({'tag':tag,'value':body,'metadata':metadata,'heavy':heavy,'units':units,'writeKey':write_key})
        proof=[v for name in proof_names for _,v in engine.leaves(engine.cas(objects[name]))]
        proof_scalars+=len(proof);nonzero_proof+=sum(v!=0 for v in proof)
        results[end]={'objects':objects,'chemicalAlignment':mapping,'densityAlignment':density_map,'jetOperands':jet_operands,'knownDimensions':dict(dimensions.known),'strong':data['strong'],
                      'proofScalars':len(proof),'nonzeroProofScalars':sum(v!=0 for v in proof),
                      'sourceDiscrepancyScalars':sum(v!=0 for name in ('MASS','MECHANICAL') for v in objects[name+'_SAVED_SOURCE_DISCREPANCY'])}
        if dimensions.constraints:raise ValueError('source trace dimension constraints')
    payload={'inputs':inputs,'packets':packets,'results':results}
    (base/'objects.pickle').write_bytes(pickle.dumps(payload,protocol=5))
    for packet in packets:
        engine.emit(packet['tag'].removeprefix('PY_S11CD_'),engine.carrier_fingerprint(packet['value']) if packet['heavy'] else packet['value'])
        engine.emit('METADATA_'+packet['tag'].removeprefix('PY_S11CD_'),packet['metadata'])
    index=engine.emission_index(engine.EMISSION_LINES)
    engine.emit('RIGHT_SOURCE_TRACE_EMISSION_LINES',index)
    engine.emit('METADATA_RIGHT_SOURCE_TRACE_EMISSION_LINES',modes.numeric_metadata(engine.cas(index),lambda path:dimensions.zero))
    summary={'sourcePins':pins,'sourcePinsAfter':{str(p):digest(p) for p in paths},'inputPackets':run_pins,
             'objectsSha256':digest(base/'objects.pickle'),'proofScalars':proof_scalars,'nonzeroProofScalars':nonzero_proof,
             'cases':{k:{n:v[n] for n in ('proofScalars','nonzeroProofScalars','sourceDiscrepancyScalars')} for k,v in results.items()},
             'objects':len(packets),'writeKeys':{v['tag']:v['writeKey'] for v in packets},'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (base/'checks.json').write_text(json.dumps(summary,indent=2)+'\n')
    if nonzero_proof or summary['sourcePins']!=summary['sourcePinsAfter']:raise ValueError('emitted source-trace reconstruction or pin discrepancy')


def validate(args):
    base=args.run_directory;summary=json.loads((base/'checks.json').read_text())
    if digest(base/'objects.pickle')!=summary['objectsSha256']:raise ValueError('trace packet pin')
    for group in [summary['sourcePins'],*summary['inputPackets'].values()]:
        for path,pin in group.items():
            if digest(Path(path))!=pin:raise ValueError(('trace source pin',path))
    payload=pickle.loads((base/'objects.pickle').read_bytes());records={};keys=set()
    for line in decoded_lines(base/'full.out'):
        tag,sep,value=line.rstrip('\n').partition(': ')
        if not sep or tag in records:raise ValueError('trace tag framing')
        records[tag]=_restore(value)
    tags=list(records)
    index_tag='PY_S11CD_RIGHT_SOURCE_TRACE_EMISSION_LINES'
    index={str(k):v for k,v in records.pop(index_tag)}
    restore_emission_index(index,tags[:tags.index(index_tag)])
    symbols,dimensions,modes=metadata_context(payload['results']['RIGHT'])
    index_metadata=modes.numeric_metadata(engine.cas(index),lambda path:dimensions.zero)
    if records.pop('PY_S11CD_METADATA_RIGHT_SOURCE_TRACE_EMISSION_LINES')!=index_metadata:raise ValueError('trace index metadata')
    for packet in payload['packets']:
        metadata=modes.numeric_metadata(packet['value'],lambda path:packet['units'][path])
        if metadata!=packet['metadata']:raise ValueError('trace metadata reconstruction')
        expected=engine.carrier_fingerprint(packet['value']) if packet['heavy'] else packet['value']
        if records.pop(packet['tag'])!=expected:raise ValueError('trace payload')
        if records.pop('PY_S11CD_METADATA_'+packet['tag'].removeprefix('PY_S11CD_'))!=packet['metadata']:raise ValueError('trace metadata')
        key=packet['writeKey']
        if not re.fullmatch(r'[a-z][A-Za-z0-9]*',key) or key in engine.IMPORT_KEYS:raise ValueError('trace write key')
        if key in keys:raise ValueError('trace write-key collision')
        keys.add(key)
    if records:raise ValueError('extra trace tags')
    summary['transcript']={'bytes':(base/'full.out').stat().st_size,'sha256':digest(base/'full.out')}
    summary['runDirectory']=str(base);summary['validatedTags']=len(tags)
    if args.publish:
        target=ROOT/'scripts/out/S11c_d_right_source_trace.out'
        if target.exists() or target.is_symlink():raise FileExistsError(target)
        with tempfile.NamedTemporaryFile(dir=target.parent,delete=False) as stream:
            staging=Path(stream.name);stream.write((base/'full.out').read_bytes())
        staging.replace(target);summary['publication']=str(target)
        (ROOT/'_measurements/S11c_d_right_source_trace_checkpoint.json').write_text(json.dumps(summary,indent=2)+'\n')
    (base/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ('proofScalars','nonzeroProofScalars','cases','objects','wallSeconds','peakRssKiB','validatedTags','transcript')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--left-run',type=Path);parser.add_argument('--right-run',type=Path)
    parser.add_argument('--validate',action='store_true');parser.add_argument('--publish',action='store_true')
    args=parser.parse_args()
    validate(args) if args.validate else compute(args)
