#!/usr/bin/env python3
"""Recompute three constant-background pencils from accepted reduced rows."""
import argparse,ast,contextlib,json,resource,signal,time
from pathlib import Path
import sympy as sp
from sympy.core.function import AppliedUndef
import S11c_d_continuum_grades as grades
import S11c_d_continuum_boundary as boundary
from S11c_d_reduced_action_source_check import restore_context

f=grades.f;engine=f.engine
PLAN=f.M/'S11c_d_uniform_source_plan.md'
PREFIX='UNIFORM_SOURCE_LAB_HELD_RHO4_CONSTANT'
PRODUCER=f.STORE/'s11c-thickness-coordinate-20260914/d_full'
CASE=('LAB_HELD','RHO4_CONSTANT')


def class_ast(source,name):
    return ast.dump(next(n for n in ast.parse(source).body if getattr(n,'name',None)==name))


def load(base):
    source,cp,path,checkpoint=grades.packet('S11c_d_reduced_action_source_checkpoint.json','reduced-action.pickle')
    r,dimensions=restore_context(source)
    f.require(source['schema']==1 and {(k,tuple(map(str,c))) for k,c in source['payloads']}=={(k,CASE) for k in engine.CLOSED_KEYS},'complete native reduced row pair')
    engine_name=str(engine.HERE.relative_to(f.ROOT))
    source_engine=Path(cp['runDirectory'])/'source'/engine_name
    f.require(f.digest(source_engine)==source['sourceDigests'][engine_name],'actual frozen reduction engine')
    for name,h in source['sourceDigests'].items():
        if name!=engine_name:f.require(f.digest(f.ROOT/name)==h,('unchanged reduction physics',name))
    manifest=json.loads((PRODUCER/'manifest.json').read_text())
    f.require(manifest['exit_code']==0 and manifest['source_hashes_before']==manifest['source_hashes_after'],'accepted stable end producer')
    frozen=PRODUCER/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    f.require(f.digest(frozen)==manifest['source_hashes_after']['scripts/S11c_d_mixing_scattering_sympy_audit.py'],'original constant-end constructor snapshot')
    current=ast.parse(engine.HERE.read_text());definitions={n.name:n for n in current.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
    helpers={'ReducedPencil','ConstantEndPencil','EdgeReduction','DimensionAnalysis','PhysicalMetadata','ChannelInput'}
    while True:
        expanded=helpers|{n.id for key in helpers for n in ast.walk(definitions[key]) if isinstance(n,ast.Name) and n.id in definitions}
        if expanded==helpers:break
        helpers=expanded
    joins={label:{name:class_ast(old.read_text(),name)==ast.dump(definitions[name]) for name in sorted(helpers)}
           for label,old in (('reduction',source_engine),('endSymbols',frozen))}
    f.require(all(all(v.values()) for v in joins.values()),('unchanged complete consumed native definition closure',joins))
    values={key:engine.named(payload,'VALUE') for (key,_),payload in source['payloads'].items()}
    pencil=engine.ReducedPencil(*(values[k] for k in engine.CLOSED_KEYS),r)
    inputs=engine.ChannelInput(r,json.loads((f.M/'S11c_d_variable_profile_development_input.json').read_text()))
    cached={};operands={str(path):f.digest(path),str(PRODUCER/'manifest.json'):f.digest(PRODUCER/'manifest.json'),str(frozen):f.digest(frozen),str(source_engine):f.digest(source_engine)}
    for end in ('REFERENCE','LEFT','RIGHT'):
        relative='symbols/'+end+'_LAB_HELD_RHO4_CONSTANT.pickle';p=PRODUCER/relative
        f.require(f.digest(p)==manifest['artifacts'][relative]['sha256'],('accepted end symbol source',end))
        cached[end]=f.unpickle(p);operands[str(p)]=f.digest(p)
    pins=dict(source['sourceDigests'])
    for p in (Path(__file__),PLAN,engine.HERE,Path(grades.__file__),Path(boundary.__file__),checkpoint,
              f.M/'S11c_d_variable_profile_development_input.json',f.ACCEPTANCE,
              f.ROOT/'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md'):
        pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for name,h in pins.items():
        f.require(f.digest(f.ROOT/name)==h,'current source');dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes((f.ROOT/name).read_bytes())
    f.save(base/'inputs.json',{'sourceFiles':pins,'inputPackets':operands,'nativeHelperJoins':joins,
        'input':inputs.specification,'case':CASE,'profileEndpointValues':{str(k):str(v) for k,v in inputs.limits.items()},
        'constantBackgrounds':{'REFERENCE':'eta=sigma=0','LEFT':'computed left values everywhere; live zero-jet contrast','RIGHT':'computed right values everywhere; live zero-jet contrast'},
        'scope':'Fresh reduced constant-symbol and coupling sources under native distributional prescription; no new mode/current or scattering result.'})
    return r,pencil,inputs,cached,pins,operands


def residual(a,b):
    if isinstance(a,sp.MatrixBase):return (a-b).applyfunc(sp.expand)
    return sp.expand(a-b)


def construct(base,r,pencil,inputs,cached,progress):
    ends=engine.ConstantEndPencil(r)
    blocks=pencil.sector_blocks()
    units={'strong':ends.strong_matrix_units(pencil),'weak':ends.weak_matrix_units(blocks),
           'action':{p:engine.PHYSICAL_METADATA.dimensions.measure(v) for p,v in engine.leaves(pencil.strong)},
           'kernel':{p:engine.PHYSICAL_METADATA.dimensions.measure(v) for p,v in engine.leaves(pencil.kernel)}}
    curl,gauge,gauge_residual=ends.curl_gauge_operands(pencil)
    f.atomic_pickle(base/'common.pickle',{'strong':pencil.strong,'kernel':pencil.kernel,'blocks':blocks,'units':units,
        'curl':curl,'gauge':gauge,'gaugeResidual':gauge_residual,'constantFourierMass':ends.constant_mass})
    progress('common_saved')
    records={}
    for label,end in (('REFERENCE',None),('LEFT',-sp.oo),('RIGHT',sp.oo)):
        progress(label+'_background_started')
        strong=ends.background(pencil.strong,end);kernel=ends.background(pencil.kernel,end);weak=ends.background(blocks,end)
        background={'strong':strong,'kernel':kernel,'blocks':weak,'endpointValues':inputs.limits,
          'profileLimitOperands':() if end is None else ends.profile_limit_operands(end)}
        f.atomic_pickle(base/(label.lower()+'-background.pickle'),background);progress(label+'_background_saved')
        strong_symbol=ends.strong_matrix(pencil,strong)
        f.atomic_pickle(base/(label.lower()+'-strong.pickle'),{'symbol':strong_symbol});progress(label+'_strong_saved')
        weak_symbol=ends.weak_matrix(weak)
        f.atomic_pickle(base/(label.lower()+'-weak.pickle'),{'symbol':weak_symbol});progress(label+'_weak_saved')
        previous=cached[label]
        pairs={'strong':(strong_symbol,previous[4]),'weak':(weak_symbol,previous[0]),'curl':(curl,previous[1])}
        f.atomic_pickle(base/(label.lower()+'-source-pairs.pickle'),pairs)
        replay={k:residual(a,b) for k,(a,b) in pairs.items()}
        f.atomic_pickle(base/(label.lower()+'-source-residuals.pickle'),replay)
        f.require(all(all(v==0 for v in a) for a in replay.values()),('fresh/accepted complete constant symbols',label))
        f.require(units['weak']==previous[2],'same inherited weak block units')
        # End values are computed from the actual approved profiles. Keep the
        # material symbols and independent zero-jet contrast live.
        actual_strong=strong_symbol.xreplace(inputs.limits);actual_weak=weak_symbol.xreplace(inputs.limits)
        coupling={name:actual_weak.extract(list(rows),list(cols)) for name,rows,cols in
                  (('TH',range(3),range(3,6)),('HT',range(3,6),range(3)))}
        profile_occurrences=set().union(*(v.atoms(AppliedUndef) for v in actual_strong))
        profile_occurrences={v for v in profile_occurrences if v.func in r.profiles.values()}
        unresolved={'integrals':tuple(sorted(actual_strong.atoms(sp.Integral)|actual_weak.atoms(sp.Integral),key=sp.default_sort_key)),
          'normalCoordinates':tuple(sorted((actual_strong.free_symbols|actual_weak.free_symbols)&{r.z,r.zp,ends.shift},key=sp.default_sort_key)),
          'regulator':actual_strong.has(r.regulator) or actual_weak.has(r.regulator),'profiles':tuple(profile_occurrences)}
        record={'background':background,'strong':actual_strong,'weak':actual_weak,'coupling':coupling,'sourceResiduals':replay,
                'unresolved':unresolved,'unitJoins':units['weak']==previous[2]}
        f.atomic_pickle(base/(label.lower()+'.pickle'),record);progress(label+'_saved')
        f.require(not any(unresolved.values()),('resolved constant symbol',label,unresolved))
        records[label]=record
    differences={end:{name:records[end]['coupling'][name]-records['REFERENCE']['coupling'][name] for name in ('TH','HT')}
                 for end in ('LEFT','RIGHT')}
    return {'records':records,'differences':differences,'units':units,'curl':curl,'gaugeResidual':gauge_residual,
        'constantFourierMass':ends.constant_mass,'profileEndpoints':inputs.limits,
        'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))}


def emit_result(result):
    units=result['units']
    for label,record in result['records'].items():
        for name,key in (('STRONG','strong'),('WEAK','weak')):
            engine.fingerprinted(PREFIX+'_'+label+'_'+name,record[key],units[key])
        for direction,rows,cols in (('TH',range(3),range(3,6)),('HT',range(3,6),range(3))):
            u={(3*i+j,):units['weak'][(6*ri+cj,)] for i,ri in enumerate(rows) for j,cj in enumerate(cols)}
            if label=='REFERENCE':
                engine.fingerprinted(PREFIX+'_REFERENCE_K_'+direction,record['coupling'][direction],u)
            else:
                for name,value in (('BASE',result['records']['REFERENCE']['coupling'][direction]),('OPERAND',record['coupling'][direction]),('RESIDUAL',result['differences'][label][direction])):
                    engine.fingerprinted(PREFIX+'_'+label+'_K_'+direction+'_'+name,value,u)
        for name in ('strong','weak'):
            engine.fingerprinted(PREFIX+'_'+label+'_SOURCE_'+name+'_RESIDUAL',record['sourceResiduals'][name],units[name])
    # Structural metadata contains only dimensionless manifest values.
    boundary.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],
        'profileEndpoints':{str(k):str(v) for k,v in result['profileEndpoints'].items()},'scope':result['scope']})


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);args=p.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));started=time.monotonic()
    def timeout(*_):raise TimeoutError('uniform source budget; preserve all completed operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    def progress(stage):
        with (base/'progress.jsonl').open('a') as stream:stream.write(json.dumps({'stage':stage,'wallSeconds':time.monotonic()-started})+'\n')
    r,pencil,inputs,cached,pins,operands=load(base);progress('sources_joined')
    result=construct(base,r,pencil,inputs,cached,progress)
    result.update(sourceFiles=pins,inputPackets=operands,scope='Three separate constant reduced operator/coupling sources. Native distributional constant-background prescription; no new modal/current/scattering construction or global decoupling claim.')
    f.atomic_pickle(base/'uniform-source.pickle',result);before=f.digest(base/'uniform-source.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result);keys={tag:'s11cdUniformSource'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')}
        boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique uniform tag');entries[tag]=grades._restore(body)
    old_emit=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('uniform full payload replay',tag));seen.add(tag)
    engine.emit=replay
    try:
        emit_result(result);boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=old_emit
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'full uniform payload/key census')
    paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        for path,fields in body:
            info={str(k):v for k,v in fields}
            dimension=info['DIMENSION_L_T_M']
            f.require(isinstance(dimension,sp.Tuple) and len(dimension)==3 and all(not v.free_symbols for v in dimension),'restored uniform dimensions')
            f.require('MULTIGRADE' in info and 'EPSILON_LAMBDA_SUPPORT' in info,'uniform independent/homotopy grades')
            paths+=1
    last='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[last]},list(entries)[:list(entries).index(last)])
    f.require(before==f.digest(base/'uniform-source.pickle') and not engine.PHYSICAL_METADATA.dimensions.constraints,'uniform source packet and dimension closure')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()) and all(f.digest(Path(n))==h for n,h in operands.items()),'unchanged sources and operands')
    checks={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':operands,'sourceResidualScalars':sum(len(a) for rec in result['records'].values() for a in rec['sourceResiduals'].values()),
       'unresolved':{k:{n:bool(v) for n,v in rec['unresolved'].items()} for k,rec in result['records'].items()},
       'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':paths,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'uniform-source.pickle'),
       'artifacts':{str(p.relative_to(base)):{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.rglob('*') if p.suffix in ('.pickle','.out') and 'source' not in p.relative_to(base).parts},
       'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',checks);progress('complete');signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
