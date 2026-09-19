#!/usr/bin/env python3
"""Construct continuum end subspace, boundary and physical-current coefficients."""
import argparse
import ast
import contextlib
import itertools
import json
import math
from pathlib import Path
import resource
import signal
import time
from types import SimpleNamespace

import numpy as np
import sympy as sp
from scipy.linalg import block_diag
import S11c_d_continuum_matrices as matrices

grades=matrices.grades; f=matrices.f; engine=f.engine
J=engine.RectangularModeJets
G=((0,0),*J.grades)
PLAN=f.M/'S11c_d_continuum_boundary_plan.md'
PREFIX='CONTINUUM_BOUNDARY_LAB_HELD_RHO4_CONSTANT'


def inverse(series):
    result={(0,0):np.linalg.inv(series[(0,0)])}
    for g in J.grades:
        known=J.multiply(series,result).get(g,np.zeros_like(result[(0,0)]))
        result[g]=-result[(0,0)]@known
    residual=J.multiply(series,result);residual[(0,0)]-=np.eye(len(result[(0,0)]))
    return result,residual


def subtract(a,b):
    return {g:a.get(g,np.zeros_like(next(iter(b.values()))))-b.get(g,np.zeros_like(next(iter(a.values())))) for g in set(a)|set(b)}


def norm(value):
    if isinstance(value,dict):return max([0.]+[norm(v) for v in value.values()])
    if isinstance(value,(list,tuple)):return max([0.]+[norm(v) for v in value])
    return float(np.max(np.abs(value))) if np.size(value) else 0.


def evaluate(series,eta,sigma):
    return sum((eta**a*sigma**b*v for (a,b),v in series.items()),np.zeros_like(series[(0,0)]))


def node(source,name):
    return ast.dump(next(n for n in ast.parse(source).body if getattr(n,'name',None)==name))


def load(base):
    _,matrix_checkpoint,mp=f.accepted_packet(f.M/'S11c_d_continuum_matrix_checkpoint.json','bound-coefficients.pickle')
    grade_packet,gcp,gpath=f.accepted_packet(matrices.GRADE_CHECKPOINT,'continuum-grades.pickle')
    rp=next(Path(p) for p in grade_packet['inputPackets'] if Path(p).name=='reduced-action.pickle')
    f.require(f.digest(rp)==grade_packet['inputPackets'][str(rp)],'reduction context identity')
    r,dimensions=f.prior.domain.momentum.source.native.source.restore_context(f.unpickle(rp))
    dimensions.__dict__.update(grade_packet['dimensionState'])
    specification=json.loads((f.M/'S11c_d_variable_profile_development_input.json').read_text())
    input_=engine.ChannelInput(r,specification)
    modal={};pairing={};pins={};operands={str(mp):f.digest(mp),str(gpath):f.digest(gpath),str(rp):f.digest(rp)};joins={}
    for end in ('REFERENCE','LEFT','RIGHT'):
        (value,known),cp,path,checkpoint=grades.packet('S11c_d_end_normalization_'+end.lower()+'_thickness_repair_checkpoint.json','modal.pickle')
        f.require(not cp['unaccountedResidualNormsAboveDiagnosticThreshold'] and cp['recordCount']==18,'accepted complete end normalization')
        dimensions.known.update(known);modal[end]=value
        frozen=(Path(cp['runDirectory'])/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py').read_text()
        names=('ModalCurrentSubspaces','RectangularModeJets','polynomial_terms')
        joins[end]={n:node(frozen,n)==node(engine.HERE.read_text(),n) for n in names}
        f.require(all(joins[end].values()),'unchanged consumed modal and invariant-pair constructors')
        args=json.loads((Path(cp['runDirectory'])/'arguments.json').read_text());pc=Path(args['pairing_checkpoint'])
        pcp=json.loads(pc.read_text());pp=(Path(pcp['runDirectory'])/'complete.pickle').resolve()
        f.require(f.digest(pp)==pcp['artifacts']['complete.pickle']['sha256'],'accepted original pairing packet')
        pairing[end]=f.unpickle(pp)[0]['result']
        f.require(cp['provenance']['pairingCacheSha256']==f.digest(pp),'modal/original pairing identity')
        native=value['NATIVE_RECORDS'];pairs={(int(v['ROOT_DISK_INDEX']),int(v['NORMAL_LIFT_SIGN'])) for v in native}
        f.require(len(native)==len(value['RECORDS'])==18 and pairs=={(i,s) for i in range(9) for s in (-1,1)},'complete isolated-root and lift census')
        for rec,old in zip(value['RECORDS'],native):
            f.require(all(rec[k]==int(old[k]) for k in ('ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','NULLITY')),'full-subspace native record join')
        f.require(value['SYMBOLIC_OPERANDS']['PENCIL_PLUS']==pairing[end]['CLOSED_PENCIL_LEGS'][0],'actual closed pencil source')
        for name,key in [('CURRENT_SLAB','SLAB_CURRENT_MATRIX'),('CURRENT_BULK','BULK_NORMAL_CURRENT_DENSITY_MATRIX')]:
            extracted=pairing[end][key].applyfunc(lambda v:dict(engine.polynomial_terms(v,(r.symbols['epsilon_shape'],))).get((2,),sp.S.Zero))
            f.require(extracted==value['SYMBOLIC_OPERANDS'][name],'actual polarized current source')
        for p in (path,pp,Path(cp['runDirectory'])/'arguments.json'):operands[str(p)]=f.digest(p)
        for p in (checkpoint,pc):pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for p in (Path(__file__),PLAN,engine.HERE,Path(f.__file__),Path(grades.__file__),Path(matrices.__file__),f.ACCEPTANCE,
              matrices.GRADE_CHECKPOINT,f.M/'S11c_d_continuum_matrix_checkpoint.json',f.M/'S11c_d_variable_profile_development_input.json',
              f.ROOT/'directives/S11c_d_SHARED_PHYSICS.md',f.ROOT/'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md'):
        pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for name,h in gcp['sourceFiles'].items():
        if name.startswith('scripts/') or name.startswith('directives/'):
            f.require(f.digest(f.ROOT/name)==h,('accepted physical source',name));pins[name]=h
    channels={}
    for end in ('LEFT','RIGHT'):
        channels[end],_,p=f.accepted_packet(f.CHANNELS,end.lower()+'.pickle');operands[str(p)]=f.digest(p)
    pins[str(f.CHANNELS.relative_to(f.ROOT))]=f.digest(f.CHANNELS)
    field_units=grade_packet['fieldUnits']
    # Derive the common physical-current unit from each nonzero matrix entry.
    units={tuple(dimensions.measure(v)[d]+field_units[i][d]+field_units[j][d] for d in range(3))
           for i in range(5) for j in range(5) if (v:=modal['REFERENCE']['SYMBOLIC_OPERANDS']['CURRENT_SLAB'][i,j])!=0}
    f.require(len(units)==1,'homogeneous source current unit');current_unit=next(iter(units))
    row_units=[]
    for i in range(5):
        row={tuple(dimensions.measure(v)[d]+field_units[j][d] for d in range(3)) for j in range(5)
             if (v:=modal['REFERENCE']['SYMBOLIC_OPERANDS']['PENCIL_PLUS'][i,j])!=0}
        f.require(len(row)==1,'homogeneous closed pencil row');row_units.append(next(iter(row)))
    for name in pins:
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((f.ROOT/name).read_bytes())
    f.save(base/'inputs.json',{'sourceFiles':pins,'inputPackets':operands,'nativeHelperJoins':joins,'input':specification,
        'fieldUnits':[list(map(str,u)) for u in field_units],'rowUnits':[list(map(str,u)) for u in row_units],
        'currentUnit':list(map(str,current_unit)),
        'scope':'Reference-anchored analytic channel gauge; current coefficients retain normalization dependence.'})
    return r,input_,modal,pairing,channels,field_units,current_unit,pins,operands


def bind(value,input_,live):
    frequencies={'s11cdPairingLeftFrequency','s11cdPairingRightFrequency'}
    fixed={v:input_.parameters['omega'] if v.name in frequencies else input_.parameters[v.name]
           for v in value.free_symbols-set(live)}
    result=value.xreplace(fixed)
    f.require(not (result.free_symbols-set(live)),'complete material binding');return result


def current_tables(expression,waves,variables,eta,sigma):
    kl,kr,ql,qr=variables
    transports={kl:sp.cancel(-sp.diff(waves[1],kl)/sp.diff(waves[1],ql)),
                kr:sp.cancel(-sp.diff(waves[0],kr)/sp.diff(waves[0],qr))}
    parameter_proofs={ (i,v):sp.diff(wave,v) for i,wave in enumerate(waves) for v in (eta,sigma)}
    f.require(all(v==0 for v in parameter_proofs.values()),'computed parameter-independent end acoustic relation')
    def derivative(value,k,q):return value.diff(k)+value.diff(q)*transports[k]
    table={}
    for a,b in G:
        for jl in range(3-a-b):
            for jr in range(3-a-b-jl):
                value=expression.diff(eta,a).diff(sigma,b)
                for _ in range(jl):value=derivative(value,kl,ql)
                for _ in range(jr):value=derivative(value,kr,qr)
                table[(a,b,jl,jr)]=(value.subs({eta:0,sigma:0})/(math.factorial(jl)*math.factorial(jr))).applyfunc(sp.cancel)
    return table,parameter_proofs


def current_pair(table,functions,left,right):
    lr={g:v.conj().T for g,v in left['R'].items()};rr=right['R']
    ls={g:v.conj().T for g,v in left['K'].items() if g!=(0,0)};rs={g:v for g,v in right['K'].items() if g!=(0,0)}
    lp=[{(0,0):np.eye(lr[(0,0)].shape[0])},ls,J.multiply(ls,ls)]
    rp=[{(0,0):np.eye(rr[(0,0)].shape[1])},rs,J.multiply(rs,rs)]
    point=(left['k'].conjugate(),right['k'],left['q'].conjugate(),right['q'])
    result={g:np.zeros((lr[(0,0)].shape[0],rr[(0,0)].shape[1]),complex) for g in G}
    for (a,b,jl,jr) in table:
        coefficient=np.asarray(functions[(a,b,jl,jr)](*point),complex)
        l=J.multiply(lp[jl],lr);r=J.multiply(rr,rp[jr])
        for (i,j),lv in l.items():
            for (m,n),rv in r.items():
                g=(a+i+m,b+j+n)
                if max(g)<=1:result[g]+=lv@coefficient@rv
    return result


def concatenate(records,key):return {g:np.column_stack([v[key][g] for v in records]) for g in G}


def phase(cluster,position):
    base=cluster['K'][(0,0)];n=len(base)
    f.require(np.array_equal(base,cluster['k']*np.eye(n)),'scalar reference momentum per isolated cluster')
    shift={g:v for g,v in cluster['K'].items() if g!=(0,0)}
    powers=[{(0,0):np.eye(n)},shift,J.multiply(shift,shift)]
    result={g:np.zeros_like(base) for g in G}
    for j,power in enumerate(powers):
        for g,v in power.items():result[g]+=np.exp(1j*cluster['k']*position)*(1j*position)**j/math.factorial(j)*v
    return result


def construct_end(end,base,r,input_,modal,pairing,accepted,current_unit):
    eta,sigma=(r.symbols[n] for n in ('eta_bg','sigma_W'));origin={eta:0,sigma:0}
    symbolic=modal[end]['SYMBOLIC_OPERANDS'];scalars=modal[end]['SCALAR_OPERANDS'];original=pairing[end]
    waves=original['ACOUSTIC_WAVE_ROWS']
    kr=original['NORMAL_LEGS'][1];kl=original['NORMAL_LEGS'][0]
    ql,qr=original['BULK_LEGS'];variables=(kl,kr,ql,qr);live=(*variables,eta,sigma)
    pencil=bind(symbolic['PENCIL_PLUS'],input_,live);waves=tuple(bind(v,input_,live) for v in waves)
    reference=bind(modal['REFERENCE']['SYMBOLIC_OPERANDS']['PENCIL_PLUS'],input_,live)
    reference_join=(pencil.subs(origin)-reference).applyfunc(sp.cancel)
    f.require(reference_join==sp.zeros(5),'zero-grade end/reference pencil join')
    jet=J(pencil,waves[0],kr,qr,eta,sigma,origin)
    f.atomic_pickle(base/(end.lower()+'-pencil.pickle'),{'pencil':pencil,'waves':waves,'referenceJoin':reference_join,
        'pencilCoefficients':jet.pencil_coefficients,'radicalCoefficients':jet.radical_coefficients})
    evaluator=sp.lambdify((kr,qr,eta,sigma),pencil,'numpy',cse=True)
    finite_residuals=[]
    for record in modal[end]['RECORDS']:
        point=(complex(record['K']),complex(record['PHYSICAL_Q']),float(input_.origin[eta]),float(input_.origin[sigma]))
        finite_residuals.append(np.asarray(evaluator(*point),complex)-record['OPERANDS']['PENCIL_PLUS'])
    f.require(norm(finite_residuals)<1e-8,'accepted finite end pencil evaluation replay')
    orientation={'LEFT':-1,'RIGHT':1}[end];clusters=[];census=[];mode_residuals={};mutation={}
    for record in modal['REFERENCE']['RECORDS']:
        i=record['INDEX'];k=complex(record['K']);q=complex(record['PHYSICAL_Q']);n=record['NULLITY']
        info={v:record[v] for v in ('INDEX','ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','NULLITY','SHEET_MEMBERSHIP','BULK_DECAY_DISK_CERTIFIED','EXACT_REAL_NORMAL')}
        direction='excluded';reason='sheet_or_bulk_decay';kind='excluded'
        if record['SHEET_MEMBERSHIP'] and record['BULK_DECAY_DISK_CERTIFIED']:
            if record['EXACT_REAL_NORMAL']:
                signed=np.real(np.diag(record['FORMS']['SIGNED_CURRENT']))*orientation
                f.require(np.all(signed>0) or np.all(signed<0),'uniform current orientation in complete open cluster')
                direction='outgoing' if np.all(signed>0) else 'incoming';kind='open';reason='physical_open'
            elif orientation*k.imag>0:direction='outgoing';kind='evanescent';reason='outward_decay'
            else:reason='outward_growth'
        info.update(direction=direction,kind=kind,reason=reason);census.append(info)
        if direction=='excluded':continue
        basis=record['FORMS']['FLUX_RIGHT' if kind=='open' else 'RIGHT']
        f.require(basis.shape==(5,n) and np.linalg.matrix_rank(basis)==n,'full reference mode subspace')
        coefficients={index:np.asarray(v,complex) for index,v in zip(jet.pencil_coefficients,jet.evaluate(k,q))}
        modes,shift,diagnostics,jacobian=J.pair(coefficients,basis)
        packet={'info':info,'k':k,'q':q,'R':modes,'K':None if shift is None else {(0,0):k*np.eye(n),**shift},
                'diagnostics':diagnostics,'jacobian':jacobian,'coefficients':coefficients,
                'amplitudeUnit':tuple(v/2 for v in current_unit) if kind=='open' else (0,0,0)}
        f.atomic_pickle(base/(end.lower()+f'-mode-{i}.pickle'),packet)
        f.require(modes is not None,'complete nonsingular invariant pair')
        mode_residuals[i]={'base':diagnostics['BASE_EQUATION'],**diagnostics['COEFFICIENT_RESIDUALS']}
        f.require(norm(mode_residuals[i])<1e-8,'full invariant-pair and gauge residuals')
        # Remove the actual eta forcing coefficient on one side only.
        changed={**coefficients,(1,0,0):np.zeros_like(coefficients[(1,0,0)])}
        mutation[i]=J.equation(changed,modes,shift)[(1,0)]-J.equation(coefficients,modes,shift)[(1,0)]
        packet['D']={g:1j*v for g,v in J.multiply(modes,packet['K']).items()}
        clusters.append(packet)
    outgoing=[v for v in clusters if v['info']['direction']=='outgoing'];incoming=[v for v in clusters if v['info']['direction']=='incoming']
    R=concatenate(outgoing,'R');D=concatenate(outgoing,'D');I=concatenate(incoming,'R');DI=concatenate(incoming,'D')
    f.require(R[(0,0)].shape==(5,5) and I[(0,0)].shape==(5,2),'complete outgoing and incident traces')
    inv,inv_residual=inverse(R);trace=J.multiply(D,inv);insertion=subtract(DI,J.multiply(trace,I))
    boundary_residual=subtract(J.multiply(trace,R),D)
    f.require(norm(inv_residual)<1e-8 and norm(boundary_residual)<1e-8,'coefficient inverse and boundary residuals')
    position=orientation*64.
    incoming_phase={g:block_diag(*(phase(c,position)[g] for c in incoming)) for g in G}
    outgoing_phase={g:block_diag(*(phase(c,-position)[g] for c in outgoing if c['info']['kind']=='open')) for g in G}
    # Every current pair is evaluated, including evanescent matching amplitudes.
    ordered=outgoing+incoming;offsets=np.cumsum([0]+[v['R'][(0,0)].shape[1] for v in ordered]);currents={};tables={};proofs={}
    expressions={'slab':symbolic['CURRENT_SLAB'],'bulk':symbolic['CURRENT_BULK']*scalars['INFINITE_DEPTH_INTEGRAL']}
    for name,expression in expressions.items():
        bound=bind(expression,input_,live);tables[name],proofs[name]=current_tables(bound,waves,variables,eta,sigma)
        functions={g:sp.lambdify(variables,v,'numpy',cse=True) for g,v in tables[name].items()}
        result={g:np.zeros((offsets[-1],offsets[-1]),complex) for g in G}
        for i,left in enumerate(ordered):
            for j,right in enumerate(ordered):
                f.require(left['q'].imag>0 and right['q'].imag>0,'actual convergent bulk-current pair')
                values=current_pair(tables[name],functions,left,right)
                for g,v in values.items():result[g][offsets[i]:offsets[i+1],offsets[j]:offsets[j+1]]=v
        currents[name]=result
    currents['total']={g:currents['slab'][g]+currents['bulk'][g] for g in G}
    hermitian={name:{g:v-v.conj().T for g,v in series.items()} for name,series in currents.items()}
    f.require(all(norm(hermitian[name])/(1+norm(series))<1e-9 for name,series in currents.items()),'full current Hermitian coefficient residuals')
    source_current_residuals={}
    for i,c in enumerate(ordered):
        if c['info']['kind']=='open':
            actual=currents['total'][(0,0)][offsets[i]:offsets[i+1],offsets[i]:offsets[i+1]]
            source_current_residuals[c['info']['INDEX']]=actual-modal['REFERENCE']['RECORDS'][c['info']['INDEX']]['FORMS']['SIGNED_CURRENT']
    f.require(norm(source_current_residuals)<1e-8,'reference physical signed-current replay')
    finite=f.boundary_map(accepted,orientation)
    finite_comparison=evaluate(trace,float(input_.origin[eta]),float(input_.origin[sigma]))-finite['traceMap']
    f.require(end!='RIGHT' or norm(mutation)>1e-10,'actual eta forcing mutation response')
    result={'end':end,'orientation':orientation,'census':census,'clusters':ordered,'trace':trace,'outgoingInverse':inv,
        'outgoing':R,'outgoingDerivative':D,'incoming':I,'incomingDerivative':DI,'insertion':insertion,
        'incomingOriginPhase':incoming_phase,'outgoingOriginPhase':outgoing_phase,
        'currents':currents,'currentTables':tables,'acousticParameterProofs':proofs,'offsets':offsets,
        'residuals':{'invariantPair':mode_residuals,'inverse':inv_residual,'boundary':boundary_residual,
            'currentHermitian':hermitian,'referenceSignedCurrent':source_current_residuals,'finitePencilReplay':finite_residuals},
        'etaForcingMutation':mutation,'finiteTraceTruncationDifference':finite_comparison,
        'traceCondition':float(np.linalg.cond(R[(0,0)])),'coordinateGauge':'reference anchored, flux-normalized at zero only'}
    f.atomic_pickle(base/(end.lower()+'-boundary.pickle'),result);return result


def emit_result(result,r):
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r
    eta,sigma=(r.symbols[n] for n in ('eta_bg','sigma_W'));eps=r.symbols['epsilon_shape'];fields=result['fieldUnits'];cu=result['currentUnit']
    zero=(0,0,0);length=(1,0,0)
    def tensor(name,array,units,grade=(0,0),quadratic=False,literal=False):
        # Scalar entries are small here; fingerprints retain the complete tensor
        # while entry metadata is replayed on the actual coefficient components.
        weight=eta**grade[0]*sigma**grade[1]*(eps**2 if quadratic else 1)
        array=np.asarray(array,complex);body=sp.ImmutableMatrix(*array.shape,[engine.FullPencilModes.number(v)*weight for v in array.ravel()])
        engine.emit(PREFIX+'_'+name,body if literal else modes.compact_fingerprint(body))
        engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(body,lambda p:units[p[0]//array.shape[1]][p[0]%array.shape[1]]))
    def unit_matrix(rows,cols,extra=zero):return [[tuple(a-b+c for a,b,c in zip(u,v,extra)) for v in cols] for u in rows]
    for end,data in result['ends'].items():
        amplitude=[c['amplitudeUnit'] for c in data['clusters'] for _ in range(c['R'][(0,0)].shape[1])]
        outgoing_units=amplitude[:5];incoming_units=amplitude[5:]
        for key,rows,columns,extra in [('trace',fields,fields,(-1,0,0)),('outgoingInverse',outgoing_units,fields,zero),
            ('outgoing',fields,outgoing_units,zero),('outgoingDerivative',fields,outgoing_units,(-1,0,0)),
            ('incoming',fields,incoming_units,zero),('incomingDerivative',fields,incoming_units,(-1,0,0)),
            ('insertion',fields,incoming_units,(-1,0,0)),('incomingOriginPhase',incoming_units,incoming_units,zero),
            ('outgoingOriginPhase',[u for u in outgoing_units if u==tuple(v/2 for v in cu)],[u for u in outgoing_units if u==tuple(v/2 for v in cu)],zero)]:
            units=unit_matrix(rows,columns,extra)
            for g,value in data[key].items():tensor(end+'_'+key+'_'+''.join(map(str,g)),value,units,g)
        current_units=[[tuple(c-a-b for c,a,b in zip(cu,u,v)) for v in amplitude] for u in amplitude]
        for name,series in data['currents'].items():
            for g,value in series.items():tensor(end+'_CURRENT_'+name+'_'+''.join(map(str,g)),value,current_units,g,True)
        for name,series in data['residuals']['currentHermitian'].items():
            for g,value in series.items():tensor(end+'_RESIDUAL_CURRENT_'+name+'_'+''.join(map(str,g)),value,current_units,g,True,True)
        for cluster in data['clusters']:
            i=cluster['info']['INDEX'];n=cluster['R'][(0,0)].shape[1];amps=[cluster['amplitudeUnit']]*n
            for key,rows,cols,extra in [('R',fields,amps,zero),('K',amps,amps,(-1,0,0))]:
                for g,array in cluster[key].items():tensor(end+f'_MODE_{i}_{key}_'+''.join(map(str,g)),array,unit_matrix(rows,cols,extra),g)
            equation_units=unit_matrix(result['rowUnits'],amps)
            diagnostics=cluster['diagnostics']
            tensor(end+f'_MODE_{i}_BASE_RESIDUAL',diagnostics['BASE_EQUATION'],equation_units,literal=True)
            for g,residuals in diagnostics['COEFFICIENT_RESIDUALS'].items():
                tensor(end+f'_MODE_{i}_EQUATION_RESIDUAL_'+''.join(map(str,g)),residuals['EQUATION'],equation_units,g,literal=True)
                # The Euclidean gauge and mixed-unit solve are numerical unit-frame
                # diagnostics, separate from the dimensioned equation residual.
                for label in ('GAUGE','LINEAR_SYSTEM'):
                    a=np.asarray(residuals[label]).reshape(-1,1)
                    tensor(end+f'_MODE_{i}_{label}_FRAME_RESIDUAL_'+''.join(map(str,g)),a,[[zero] for _ in a],g,literal=True)
            tensor(end+f'_MODE_{i}_ETA_FORCING_MUTATION',data['etaForcingMutation'][i],equation_units,(1,0),literal=True)
            if i in data['residuals']['referenceSignedCurrent']:
                tensor(end+f'_MODE_{i}_SIGNED_CURRENT_RESIDUAL',data['residuals']['referenceSignedCurrent'][i],
                       [[zero]*n for _ in range(n)],quadratic=True,literal=True)
        tensor(end+'_FINITE_TRACE_TRUNCATION_DIFFERENCE',data['finiteTraceTruncationDifference'],unit_matrix(fields,fields,(-1,0,0)),literal=True)
        # Literal small residuals, separated by grade and source equation unit.
        for name,series,units in [('inverse',data['residuals']['inverse'],unit_matrix(fields,fields)),
            ('boundary',data['residuals']['boundary'],unit_matrix(fields,outgoing_units,(-1,0,0)))]:
            for g,array in series.items():
                weight=eta**g[0]*sigma**g[1];body=sp.ImmutableMatrix(*array.shape,[engine.FullPencilModes.number(v)*weight for v in array.ravel()])
                tag=PREFIX+'_'+end+'_RESIDUAL_'+name+'_'+''.join(map(str,g));engine.emit(tag,body)
                engine.emit('METADATA_'+tag,modes.numeric_metadata(body,lambda p,u=units,n=array.shape[1]:u[p[0]//n][p[0]%n]))
        grades.structural(PREFIX+'_'+end+'_CENSUS',{'candidates':data['census'],'traceCondition':data['traceCondition'],
            'gauge':data['coordinateGauge'],'currentPairCount':len(data['clusters'])**2})
    grades.structural(PREFIX+'_SOURCE_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],
        'scope':result['scope'],'grades':G})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('continuum boundary budget; preserve completed packets')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);started=time.monotonic()
    r,input_,modal,pairing,channels,fields,current,pins,operands=load(base)
    ends={}
    for end in ('LEFT','RIGHT'):
        ends[end]=construct_end(end,base,r,input_,modal,pairing,channels[end],current)
        f.save(base/'record-inventory.json',{p.name:{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.glob('*.pickle')})
    inputs=json.loads((base/'inputs.json').read_text())
    result={'ends':ends,'fieldUnits':fields,'rowUnits':[tuple(map(sp.Rational,u)) for u in inputs['rowUnits']],
            'currentUnit':current,'sourceFiles':pins,'inputPackets':operands,
            'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions)),
            'scope':'Complete finite modal-boundary and current coefficients in reference-anchored channel coordinates; no continuum solve yet.'}
    f.atomic_pickle(base/'continuum-boundary.pickle',result);before=f.digest(base/'continuum-boundary.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result,r)
        keys={tag:'s11cdContinuumBoundary'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')}
        grades.structural(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);grades.structural(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique boundary tag');entries[tag]=grades._restore(body)
    original=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('boundary emission replay',tag));seen.add(tag)
    engine.emit=replay
    try:
        emit_result(result,r);grades.structural(PREFIX+'_WRITE_KEYS',keys);grades.structural(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=original
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'complete boundary emission/key replay')
    metadata_paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        structural=tag.endswith(('_CENSUS','_SOURCE_MANIFEST','_WRITE_KEYS','_EMISSION_LINES'))
        for item in body:
            fields_={str(k):v for k,v in (item[1] if structural else item)}
            f.require(len(fields_['DIMENSION_L_T_M'])==3 and all(not v.free_symbols for v in fields_['DIMENSION_L_T_M']),'resolved boundary unit')
            f.require('MULTIGRADE' in fields_ and 'EPSILON_LAMBDA_SUPPORT' in fields_,'boundary grades')
            metadata_paths+=1 if structural else len(fields_['PATHS'])
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    f.require(before==f.digest(base/'continuum-boundary.pickle'),'boundary packet pre/post hash')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()) and all(f.digest(Path(n))==h for n,h in operands.items()),'source/input post hashes')
    summary={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':operands,'ends':{end:{'candidates':len(v['census']),
        'clusters':len(v['clusters']),'directions':sum(c['R'][(0,0)].shape[1] for c in v['clusters']),
        'residualMaxima':{k:norm(a) for k,a in v['residuals'].items()},'etaMutationMaximum':norm(v['etaForcingMutation']),
        'finiteTraceTruncationMaximum':norm(v['finiteTraceTruncationDifference']),'traceCondition':v['traceCondition']} for end,v in ends.items()},
        'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':metadata_paths,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'continuum-boundary.pickle'),
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
