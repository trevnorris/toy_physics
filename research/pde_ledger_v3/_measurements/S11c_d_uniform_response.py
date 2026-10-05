#!/usr/bin/env python3
"""Full-subspace current and constant-background matching controls."""
import argparse,contextlib,json,resource,signal,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import sympy as sp
import S11c_d_continuum_boundary as b
f=b.f;engine=f.engine;grades=b.grades
PLAN=f.M/'S11c_d_uniform_response_plan.md';SOURCE=f.M/'S11c_d_uniform_source_checkpoint.json'
PREFIX='UNIFORM_RESPONSE_LAB_HELD_RHO4_CONSTANT'


def load(base):
    r,inp,modal,pairing,_,fields,current,pins,operands=b.load(base)
    source,checkpoint,path=f.accepted_packet(SOURCE,'uniform-source.pickle');operands[str(path)]=f.digest(path)
    for label in modal:
        cp=json.loads((f.M/('S11c_d_end_normalization_'+label.lower()+'_thickness_repair_checkpoint.json')).read_text())
        original=next(Path(n) for n in checkpoint['checks']['inputPackets'] if Path(n).name==label+'_LAB_HELD_RHO4_CONSTANT.pickle')
        f.require(f.digest(original)==cp['provenance']['cacheSha256'],'fresh uniform and isolated-root native source identity')
        frozen=(Path(cp['runDirectory'])/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py').read_text()
        f.require(grades.method(frozen,'FullPencilModes','analytic')==grades.method(engine.HERE.read_text(),'FullPencilModes','analytic'),'unchanged actual branch/radical continuation helper')
        original_physics=cp['provenance']['producerSources']
        for n in ('scripts/S11c_b_exports.py','scripts/S11c_c1_exports.py','scripts/S11c_c2_exports.py','directives/S11c_d_SHARED_PHYSICS.md'):
            f.require(f.digest(f.ROOT/n)==original_physics[n],'same isolated-root physical source')
    for p in (Path(__file__),PLAN,SOURCE,Path(f.__file__),Path(b.__file__)):
        pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for n,h in pins.items():
        f.require(f.digest(f.ROOT/n)==h,'current source');target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((f.ROOT/n).read_bytes())
    inputs=json.loads((base/'inputs.json').read_text());inputs.update(sourceFiles=pins,inputPackets=operands,
        scope='Three matched constant backgrounds; reused isolated roots after exact fresh-symbol/source joins; newly evaluated full subspaces, currents and response.')
    f.save(base/'inputs.json',inputs)
    return r,inp,modal,pairing,source,fields,current,pins,operands


def construct_modes(base,label,r,inp,original,pairing,source):
    ends=engine.ConstantEndPencil.__new__(engine.ConstantEndPencil);ends.r=r;ends.kn=sp.Symbol('s11cdSpectralNormalMomentum',real=True)
    modes=engine.FullPencilModes(ends,source['curl'],source['units']['weak'])
    algebraic,relation,branches=modes.analytic(source['records'][label]['strong'])
    origin={s:0 if label=='REFERENCE' else v for s,v in inp.origin.items()}
    mapping=inp.mapping(algebraic,relation,(modes.k,modes.q,*origin));mapping.update(origin)
    physical=algebraic.xreplace(mapping);curve=relation.xreplace(mapping)
    f.require(not(physical.free_symbols|curve.free_symbols)-{modes.k,modes.q},'complete fresh uniform material binding')
    f.require(all(v==0 for v in branches),'native real-axis branch join')
    functions={'fresh':sp.lambdify((modes.k,modes.q),physical,'numpy',cse=True),
               'curve':sp.lambdify((modes.k,modes.q),curve,'numpy',cse=True)}
    variables=tuple(pairing['FREQUENCY_LEGS'])+tuple(pairing['NORMAL_LEGS'])+tuple(pairing['BULK_LEGS'])
    symbolic={name:original['SYMBOLIC_OPERANDS'][name] for name in ('PENCIL_PLUS','NORMAL_PENCIL_PLUS','FREQUENCY_PENCIL_PLUS','CURRENT_SLAB','CURRENT_BULK')}
    symbolic['INFINITE_DEPTH_INTEGRAL']=original['SCALAR_OPERANDS']['INFINITE_DEPTH_INTEGRAL']
    bound={name:b.bind(value,inp,(*variables,*origin)).subs(origin) for name,value in symbolic.items()}
    evaluate={name:sp.lambdify(variables,value,'numpy',cse=True) for name,value in bound.items()}
    f.atomic_pickle(base/(label.lower()+'-symbolic.pickle'),{'freshPencil':physical,'curve':curve,'branchResiduals':branches,'variables':variables,'bound':bound,'origin':origin})
    records=[]
    for old in original['RECORDS']:
        i=old['INDEX'];k=complex(old['K']);q=complex(old['Q']);pq=complex(old['PHYSICAL_Q']);w=float(inp.parameters['omega'])
        point=(w,w,k.conjugate(),k,pq.conjugate(),pq)
        values={n:np.asarray(fn(*point),complex) for n,fn in evaluate.items() if n!='INFINITE_DEPTH_INTEGRAL'}
        matrix=np.asarray(functions['fresh'](k,q),complex);matrix_join=matrix-values['PENCIL_PLUS']
        accepted_join=matrix-old['OPERANDS']['PENCIL_PLUS'];scale=1+b.norm(matrix)
        u,s,vh=np.linalg.svd(matrix);n=int(np.sum(s<1e-8*max(1.,s[0])))
        f.require(n==old['NULLITY'] and n>0,'complete native root nullity')
        raw_right=vh.conj().T[:,-n:];left=u[:,-n:]
        # Procrustes is only a complete-subspace coordinate alignment.
        uu,_,vv=np.linalg.svd(raw_right.conj().T@old['FORMS']['RIGHT']);right=raw_right@(uu@vv)
        normal=left.conj().T@values['NORMAL_PENCIL_PLUS']@right
        frequency=left.conj().T@values['FREQUENCY_PENCIL_PLUS']@right
        normal_rank=np.linalg.matrix_rank(normal,tol=1e-9*max(1.,np.linalg.norm(normal,2)))
        frequency_rank=np.linalg.matrix_rank(frequency,tol=1e-9*max(1.,np.linalg.norm(frequency,2)))
        residuals={'freshPairingPencil':matrix_join/scale,'acceptedPencil':accepted_join/scale,
          'rightKernel':matrix@right/scale,'leftKernel':matrix.conj().T@left/scale,
          'rightFrameGram':right.conj().T@right-np.eye(n),'leftFrameGram':left.conj().T@left-np.eye(n),
          'rightSubspace':right@right.conj().T-old['FORMS']['RIGHT_COORDINATE_PROJECTOR'],
          'leftSubspace':left@left.conj().T-old['FORMS']['LEFT_COORDINATE_PROJECTOR']}
        for key in ('CURRENT_SLAB','CURRENT_BULK'):
            residuals[key+'SourceReplay']=(values[key]-old['OPERANDS'][key])/(1+b.norm(old['OPERANDS'][key]))
        info={key:old[key] for key in ('INDEX','ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','K','Q','PHYSICAL_Q','OMEGA','NULLITY','SHEET_MEMBERSHIP','EXACT_REAL_NORMAL','BULK_DECAY_DISK_CERTIFIED')}
        native=original['NATIVE_RECORDS'][i];info['CLASSIFIER_STATUS']=native['CLASSIFIER_STATUS']
        record={'info':info,'right':right,'rawRight':raw_right,'left':left,'singularValues':s,'pencil':matrix,
          'normalPairing':normal,'frequencyPairing':frequency,'normalRank':int(normal_rank),'frequencyRank':int(frequency_rank),
          'currentOperands':{n:values[n] for n in ('CURRENT_SLAB','CURRENT_BULK')},'residuals':residuals,
          'waveResidual':complex(functions['curve'](k,q)),'currentDefined':False}
        if old['BULK_DECAY_DISK_CERTIFIED']:
            depth=complex(evaluate['INFINITE_DEPTH_INTEGRAL'](*point));current=values['CURRENT_SLAB']+depth*values['CURRENT_BULK'];gram=right.conj().T@current@right
            vals,rotation=np.linalg.eigh((gram+gram.conj().T)/2)
            record.update(depthIntegral=depth,currentMatrix=current,currentGram=gram,currentEigenvalues=vals)
            residuals['currentHermitian']=(gram-gram.conj().T)/(1+b.norm(gram))
            eligible=bool(old['SHEET_MEMBERSHIP'] and old['EXACT_REAL_NORMAL'] and normal_rank==n and frequency_rank==n and np.all(abs(vals)>1e-9*max(1.,np.linalg.norm(gram,2))))
            record['currentDefined']=eligible
            if eligible:
                transform=rotation@np.diag(1/np.sqrt(abs(vals)));flux=right@transform;record.update(fluxRight=flux,fieldToFlux=transform,signedCurrent=np.diag(np.sign(vals)))
                residuals['fluxNormalization']=flux.conj().T@current@flux-record['signedCurrent']
        info['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED']=record['currentDefined']
        f.atomic_pickle(base/(label.lower()+f'-mode-{i}.pickle'),record)
        f.require(b.norm(residuals)<1e-8 and abs(record['waveResidual'])<1e-7*(1+abs(q)**2+abs(k)**2),'full fresh candidate and current replay')
        f.require(normal_rank==old['NORMAL_PAIRING_RANK'] and frequency_rank==old['FREQUENCY_PAIRING_RANK'],'full derivative pairing ranks')
        f.require(record['currentDefined']==bool(old['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED']),'current-domain classification')
        records.append(record)
    f.require(len(records)==18 and {(r['info']['ROOT_DISK_INDEX'],r['info']['NORMAL_LIFT_SIGN']) for r in records}=={(i,s) for i in range(9) for s in (-1,1)},'complete root/lift census')
    return records,evaluate


def match(base,label,records,evaluate):
    by_index={r['info']['INDEX']:r for r in records};pair_cache={}
    def pair(i,j):
        if (i,j) not in pair_cache:
            a,c=by_index[i],by_index[j];ai,ci=a['info'],c['info'];qa,qc=complex(ai['PHYSICAL_Q']),complex(ci['PHYSICAL_Q'])
            f.require(qa.imag>0 and qc.imag>0,'actual current depth-decay domain')
            point=(complex(ai['OMEGA']),complex(ci['OMEGA']),complex(ai['K']).conjugate(),complex(ci['K']),qa.conjugate(),qc)
            slab=np.asarray(evaluate['CURRENT_SLAB'](*point),complex);bulk=np.asarray(evaluate['CURRENT_BULK'](*point),complex);depth=complex(evaluate['INFINITE_DEPTH_INTEGRAL'](*point))
            item={'slab':slab,'bulk':bulk,'depthIntegral':depth,'total':slab+depth*bulk,'records':(i,j)}
            f.atomic_pickle(base/(label.lower()+f'-current-{i}-{j}.pickle'),item);pair_cache[i,j]=item
        return pair_cache[i,j]['total']
    eligible=[r for r in records if r['currentDefined']];channels={};candidates=[]
    for r in records:
        bases={'RIGHT':r['right']}
        if r['currentDefined']:bases['FLUX_RIGHT']=r['fluxRight']
        candidates.append({'INFO':r['info'],'BASES':bases})
    open_entries=[(r,j) for r in eligible for j in range(r['info']['NULLITY'])]
    flux=np.asarray([[a['fluxRight'][:,i].conj()@pair(a['info']['INDEX'],c['info']['INDEX'])@c['fluxRight'][:,j] for c,j in open_entries] for a,i in open_entries])
    for end,orientation in (('LEFT',-1),('RIGHT',1)):
        entries=[]
        for column,(r,j) in enumerate(open_entries):
            sign=orientation*r['signedCurrent'][j,j].real
            f.require(sign!=0,'nonzero computed open direction')
            entries.append({'END':end,'RECORD_INDEX':r['info']['INDEX'],'BASIS_COLUMN':j,'MATRIX_COLUMN':column,'DIRECTION':'OUTGOING' if sign>0 else 'INCOMING'})
        channels[end]=f.boundary_map({'CANDIDATES':candidates,'CHANNELS':entries,'OUTWARD_CURRENT':orientation*flux},orientation)
    global_modes=channels['LEFT']['outgoing']+channels['RIGHT']['outgoing'];f.require(len(global_modes)==10,'complete physical homogeneous solution space')
    anchors=np.asarray([64. if v['k'].imag<0 else -64. if v['k'].imag>0 else 0. for v in global_modes]);ks=np.asarray([v['k'] for v in global_modes]);vectors=np.column_stack([v['vector'] for v in global_modes])
    def at(x):
        phase=np.exp(1j*ks*(x-anchors));return vectors*phase[None,:],vectors*(1j*ks*phase)[None,:]
    a=[];rhs=np.zeros((10,4),complex);values={}
    for block,(end,x) in enumerate((('LEFT',-64.),('RIGHT',64.))):
        v,d=at(x);a.append(d-channels[end]['traceMap']@v);rhs[5*block:5*block+5,2*block:2*block+2]=channels[end]['incomingBoundaryData'];values[end]=(v,d)
    matrix=np.vstack(a);rs=np.max(abs(matrix),axis=1);cs=np.linalg.norm(matrix/rs[:,None],axis=0);balanced=matrix/rs[:,None]/cs[None,:]
    coefficient=np.linalg.solve(balanced,rhs/rs[:,None])/cs[:,None]
    independent,_,rank,singular=np.linalg.lstsq(balanced,rhs/rs[:,None],rcond=None);independent/=cs[:,None]
    f.atomic_pickle(base/(label.lower()+'-matching-system.pickle'),{'matrix':matrix,'rhs':rhs,'anchors':anchors,'momenta':ks,'vectors':vectors,'channels':channels,'values':values,'rowScale':rs,'columnScale':cs})
    residual=(matrix@coefficient-rhs)/rs[:,None];mutation=rhs.copy();mutation[:5,:2]*=-1;mutation_residual=(matrix@coefficient-mutation)/rs[:,None]
    modal={};out=[];labels=[];full_currents={};boundary_checks={}
    for block,end in enumerate(('LEFT','RIGHT')):
        v,d=values[end];field=v@coefficient;trace=field.copy();trace[:,2*block:2*block+2]-=channels[end]['incomingValues']
        amp=np.linalg.solve(channels[end]['right'],trace);modal[end]=amp
        boundary_checks[end]=(d@coefficient-channels[end]['traceMap']@field-rhs[5*block:5*block+5])/rs[5*block:5*block+5,None]
        selected=channels[end]['outgoing']+channels[end]['incoming'];amplitudes=np.vstack((amp,np.eye(4)[2*block:2*block+2]))
        current=np.asarray([[left['vector'].conj()@pair(left['RECORD_INDEX'],right['RECORD_INDEX'])@right['vector'] for right in selected] for left in selected])
        full_currents[end]={'matrix':current,'amplitudes':amplitudes,'quadratic':amplitudes.conj().T@current@amplitudes,'outward':channels[end]['orientation']*(amplitudes.conj().T@current@amplitudes),'labels':selected}
        for j,item in enumerate(channels[end]['outgoing']):
            if item['kind']=='open':out.append(amp[j]);labels.append((end,item))
    scattering=np.vstack(out);din=np.asarray([np.exp(1j*v['k']*x) for end,x in (('LEFT',-64.),('RIGHT',64.)) for v in channels[end]['incoming']]);dout=np.asarray([np.exp(-1j*v['k']*(-64. if end=='LEFT' else 64.)) for end,v in labels]);origin=dout[:,None]*scattering*din[None,:]
    result={'coefficients':coefficient,'independentDifference':coefficient-independent,'rank':int(rank),'condition':float(singular[0]/singular[-1]),'residual':residual,'boundaryResiduals':boundary_checks,'mutationResidual':mutation_residual,
      'channels':channels,'modal':modal,'scattering':scattering,'originScattering':origin,'incomingPhases':din,'outgoingPhases':dout,'fullEndCurrents':full_currents,'openCurrent':flux,'globalModes':global_modes,'currentPairs':pair_cache}
    f.atomic_pickle(base/(label.lower()+'-response.pickle'),result)
    f.require(rank==10 and b.norm(residual)<1e-9 and b.norm(coefficient-independent)/(1+b.norm(coefficient))<1e-9 and b.norm(boundary_checks)<1e-9,'uniform matching rank and independent/direct equations')
    f.require(b.norm(mutation_residual)>0,'actual one-sided incident sign response')
    return result


def emit_result(result,r):
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r;modes.eta=r.symbols['eta_bg'];modes.sigma=r.symbols['sigma_W'];eps=r.symbols['epsilon_shape']
    def tensor(name,value,unit=lambda p:(0,0,0),power=0):
        a=np.asarray(value,complex)
        if a.ndim==1:a=a.reshape(-1,1)
        body=eps**power*sp.ImmutableMatrix(*a.shape,[modes.number(v) for v in a.ravel()]);engine.emit(PREFIX+'_'+name,modes.compact_fingerprint(body));engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(body,unit))
    for label,data in result['backgrounds'].items():
        for rec in data['modes']:
            i=rec['info']['INDEX'];n=rec['info']['NULLITY'];tag=label+'_MODE_'+str(i)
            tensor(tag+'_RIGHT_FIELD_BASIS',rec['right'],lambda p:result['fieldUnits'][p[0]//n],1)
            tensor(tag+'_SINGULAR_VALUES',rec['singularValues'])
            for name,value in rec['residuals'].items():tensor(tag+'_FRAME_RESIDUAL_'+name,value)
            if rec['currentDefined']:tensor(tag+'_SIGNED_CURRENT',rec['signedCurrent'],power=2)
        response=data['response'];tensor(label+'_BOUNDARY_SCATTERING',response['scattering']);tensor(label+'_ORIGIN_SCATTERING',response['originScattering'])
        tensor(label+'_MATCHING_SCALED_RESIDUAL',response['residual']);tensor(label+'_INCIDENT_SIGN_MUTATION',response['mutationResidual'])
        for end,current in response['fullEndCurrents'].items():tensor(label+'_'+end+'_FULL_NORMAL_CURRENT_PER_INCIDENT_UNIT',current['quadratic'],power=2)
        boundary_info=[{'index':v['info']['INDEX'],'rootDisk':v['info']['ROOT_DISK_INDEX'],'lift':v['info']['NORMAL_LIFT_SIGN'],'nullity':v['info']['NULLITY'],'sheet':v['info']['SHEET_MEMBERSHIP'],'bulkDecay':v['info']['BULK_DECAY_DISK_CERTIFIED'],'exactRealNormal':v['info']['EXACT_REAL_NORMAL'],'openCurrentDefined':v['currentDefined']} for v in data['modes']]
        b.structural_flags(PREFIX+'_'+label+'_CENSUS',boundary_info)
    b.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],'scope':result['scope'],
      'numericalGradeConvention':'Finite evaluations at each declared grade origin; frame residuals are dimensionless coefficients. Mode/current amplitude powers remain explicit.'})


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);args=p.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));start=time.monotonic()
    def timeout(*_):raise TimeoutError('uniform response budget; preserve completed modes/current/matching operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    r,inp,modal,pairing,source,fields,current,pins,operands=load(base);backgrounds={}
    for label in ('REFERENCE','LEFT','RIGHT'):
        modes,evaluate=construct_modes(base,label,r,inp,modal[label],pairing[label],source);response=match(base,label,modes,evaluate);backgrounds[label]={'modes':modes,'response':response}
    result={'backgrounds':backgrounds,'fieldUnits':fields,'currentUnit':current,'sourceFiles':pins,'inputPackets':operands,
      'scope':'Three constant-background finite-contrast physical selected-mode matching controls. All root candidates and complete subspaces retained. No new variable-profile transparency, continuum expansion, global exceptional coverage or frequency-pole claim.'}
    f.atomic_pickle(base/'uniform-response.pickle',result);before=f.digest(base/'uniform-response.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result,r);keys={tag:'s11cdUniformResponse'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')};b.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);b.structural_flags(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique uniform response tags');entries[tag]=grades._restore(body)
    original=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('full uniform response emission',tag));seen.add(tag)
    engine.emit=replay
    try:emit_result(result,r);b.structural_flags(PREFIX+'_WRITE_KEYS',keys);b.structural_flags(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=original
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'complete uniform response key/metadata census')
    paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        structural=tag.endswith(('_CENSUS','_MANIFEST','_WRITE_KEYS','_EMISSION_LINES'))
        for item in body:
            v={str(k):x for k,x in (item[1] if structural else item)};d=v['DIMENSION_L_T_M'];f.require(len(d)==3 and all(not x.free_symbols for x in d) and 'MULTIGRADE' in v and 'EPSILON_LAMBDA_SUPPORT' in v,'uniform response resolved units/grades');paths+=1 if structural else len(v['PATHS'])
    last='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[last]},list(entries)[:list(entries).index(last)])
    f.require(before==f.digest(base/'uniform-response.pickle') and not engine.PHYSICAL_METADATA.dimensions.constraints,'uniform response packet and dimension closure')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()) and all(f.digest(Path(n))==h for n,h in operands.items()),'uniform unchanged sources/operands')
    summary={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':operands,'backgrounds':{name:{'candidates':len(v['modes']),'basisDirections':sum(a['info']['NULLITY'] for a in v['modes']),'modeResidual':max(b.norm(a['residuals']) for a in v['modes']),'currentPairs':len(v['response']['currentPairs']),'matchingRank':v['response']['rank'],'matchingCondition':v['response']['condition'],'matchingResidual':b.norm(v['response']['residual']),'independentDifference':b.norm(v['response']['independentDifference'])} for name,v in backgrounds.items()},
      'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':paths,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'uniform-response.pickle'),
      'artifacts':{str(p.relative_to(base)):{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.rglob('*') if p.suffix in ('.pickle','.out') and 'source' not in p.relative_to(base).parts},'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
