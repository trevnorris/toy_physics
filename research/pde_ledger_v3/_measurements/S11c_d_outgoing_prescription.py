#!/usr/bin/env python3
"""Bounded spatial outgoing-prescription candidate from accepted saved inputs.

No producer import, root search, mode/current constructor or profile solve.
Independent method review and whole-job containment are required before use.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import time
import traceback

ROOT = Path(__file__).resolve().parents[3]
STORE = ROOT/'_scratch/s11c'
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


def require(test, message):
    if not test:
        raise ValueError(message)


def route(path):
    p = Path(path)
    return {'path': str(p), 'canonicalPath': str(p.resolve(strict=True)),
            'bytes': p.stat().st_size, 'sha256': digest(p)}


def containment():
    group = next(v[3:] for v in Path('/proc/self/cgroup').read_text().splitlines() if v.startswith('0::'))
    root = Path('/sys/fs/cgroup')/group.lstrip('/')
    r = {n:(root/n).read_text().strip() for n in ('memory.max','memory.swap.max','pids.max')}
    r.update(nice=os.getpriority(os.PRIO_PROCESS,0), affinity=sorted(os.sched_getaffinity(0)),
             threads={n:os.environ.get(n) for n in THREADS})
    require(r['memory.max']==str(2*1024**3) and r['memory.swap.max']=='0' and r['pids.max']=='32'
        and r['nice']>=15 and len(r['affinity'])==1 and all(v=='1' for v in r['threads'].values()),
        'required containment missing; no scientific library imported')
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    def timeout(*_):
        raise TimeoutError('900-second native limit; preserve completed work')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    return r


class SavedCodec(pickle.Unpickler):
    def find_class(self,module,name):
        if not (module.startswith(('sympy.','numpy.')) or module in ('sympy','numpy','builtins','collections')):
            raise pickle.UnpicklingError((module,name))
        return super().find_class(module,name)


class Journal:
    def __init__(self,out):
        self.out,self.records=out,[]

    def value(self,name,value):
        p=self.out/name;p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('xb') as f:
            pickle.dump(value,f,protocol=4);f.flush();os.fsync(f.fileno())
        return route(p)

    def op(self,name,function,*args):
        folder='operations/%03d-%s'%(len(self.records),name)
        r={'name':name,'input':self.value(folder+'/input.pickle',args),
           'startedUtc':datetime.now(timezone.utc).isoformat()}
        save(self.out/folder/'started.json',r)
        result=function(*args)
        r.update(value=self.value(folder+'/value.pickle',result),finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json',r);self.records.append(r)
        return result


def construct(spec,out):
    import sympy as sp
    import numpy as np
    j=Journal(out)
    actual={key:route(r['path']) for key,r in spec['inputs'].items()}
    save(out/'consumed.json',actual)
    require(actual==spec['inputs'],'frozen consumed routes/bytes changed')
    def read_json(key):
        return json.loads(Path(actual[key]['path']).read_text())
    def read_pickle(key):
        with Path(actual[key]['path']).open('rb') as f:
            return SavedCodec(f).load()
    kernel_cp=read_json('kernelCheckpoint');modal_cp=read_json('modalCheckpoint');pair_cp=read_json('pairingCheckpoint')
    require(kernel_cp['status']=='ACCEPTED_REFERENCE_INVERSE_AND_REGULAR_SPECTRAL_DENSITY_ONLY','accepted kernel ingredient')
    require(modal_cp['end']==pair_cp['end']=='REFERENCE' and modal_cp['case']=='LAB_HELD_RHO4_CONSTANT','own reference case')
    require(not modal_cp['unaccountedResidualNormsAboveDiagnosticThreshold'] and pair_cp['retainedNonzeroScalars']==0,'accepted modal/pairing scope')
    require(modal_cp['artifacts']['modal.pickle']['sha256']==actual['modalPacket']['sha256'],'modal checkpoint join')
    require(modal_cp['provenance']['inputSha256']==actual['originalPhysicalInput']['sha256'],'modal original physical input join')
    require(pair_cp['artifacts']['complete.pickle']['sha256']==actual['pairingPacket']['sha256']
        ==modal_cp['provenance']['pairingCacheSha256'],'preserved relocated pairing identity')
    for key in spec['kernelArtifactKeys']:
        rel=spec['kernelArtifactKeys'][key]
        require(kernel_cp['artifactSha256'][rel]==actual[key]['sha256'],'kernel checkpoint artifact join: '+key)
    selected=j.op('restore-selected',read_pickle,'selectedInputs')
    context=j.op('restore-kernel-context',read_pickle,'kernelContext')
    units=j.op('restore-kernel-units',read_pickle,'kernelUnits')
    modal,known=j.op('restore-modal',read_pickle,'modalPacket')
    pairing,unused_known=j.op('restore-pairing',read_pickle,'pairingPacket')
    pairing=pairing['result']
    source=selected['reference']['strong'];kn=context['normalMomentum']
    physical=read_json('physicalInput');original=read_json('originalPhysicalInput')
    require(physical==selected['physicalInput'],'same physical binding as accepted kernel')
    require(original['unit_frame']==physical['unit_frame'] and all(physical['parameters'].get(k)==v for k,v in original['parameters'].items()),'original/development parameter join')
    require(modal['SYMBOLIC_OPERANDS']['PENCIL_PLUS']==pairing['CLOSED_PENCIL_LEGS'][0],'modal/pairing exact pencil identity')
    kl,kr=pairing['NORMAL_LEGS'];ql,qr=pairing['BULK_LEGS'];wl,wr=pairing['FREQUENCY_LEGS']
    parameters={k:sp.Rational(v) for k,v in physical['parameters'].items()}
    def binding(expression,live=()):
        mapping={}
        for s in expression.free_symbols-set(live):
            if s.name in ('eta_bg','sigma_W'):
                mapping[s]=sp.S.Zero
            else:
                require(s.name in parameters,'missing material coefficient '+s.name)
                mapping[s]=parameters[s.name]
        return mapping
    source_bind=binding(source,(kn,))
    j.value('reference-physical-binding.pickle',{
        'suppliedPhysicalInput':physical,'effectiveReferenceGrades':{'eta_bg':sp.S.Zero,'sigma_W':sp.S.Zero},
        'actualSourceBinding':source_bind,'unchangedOtherPhysicalParameters':parameters})
    fixed_source=j.op('bind-accepted-symbol',lambda a,b:a.xreplace(b),source,source_bind)
    wave=pairing['ACOUSTIC_WAVE_ROWS'][0]
    wave_map={**binding(wave,(kr,qr,wr)),kr:kn,wr:parameters['omega']}
    fixed_wave=j.op('bind-physical-wave',lambda a,b:sp.cancel(a.xreplace(b)),wave,wave_map)
    polynomial=sp.Poly(fixed_wave,qr)
    require(polynomial.degree()==2 and polynomial.nth(1)==0,'quadratic physical acoustic branch')
    square=j.op('physical-q-square',lambda a,b:sp.cancel(-a/b),polynomial.nth(0),polynomial.nth(2))
    positive=sp.Poly(-square,kn)
    branch_evidence={'minusSquare':positive.as_expr(),'degree':positive.degree(),
        'linear':positive.nth(1),'quadratic':positive.nth(2),'constant':positive.nth(0),
        'normalUnits':tuple(known[kr]),'physicalBulkUnits':tuple(known[qr]),
        'kernelMomentumUnits':units['spectralMeasureUnit']}
    j.value('branch-and-coordinate-evidence.pickle',branch_evidence)
    require(positive.degree()==2 and positive.nth(1)==0 and positive.nth(2).is_positive is True
        and positive.nth(0).is_positive is True,'approved slice must be strictly evanescent on the entire real normal axis')
    require(tuple(known[kr])==tuple(known[qr])==tuple(units['spectralMeasureUnit']),'physical normal/bulk coordinate units')
    q_out=sp.I*sp.sqrt(-square)
    pencil=modal['SYMBOLIC_OPERANDS']['PENCIL_PLUS']
    modal_map={**binding(pencil,(kr,qr,wr)),kr:kn,qr:q_out,wr:parameters['omega']}
    fixed_modal=j.op('bind-saved-modal-symbol',lambda a,b:a.xreplace(b),pencil,modal_map)
    raw=fixed_modal-fixed_source
    symbol_join=j.op('actual-symbol-join',lambda a:a.applyfunc(sp.simplify),raw)
    require(symbol_join==sp.zeros(5),'actual full modal/kernel symbol join')
    field_names=('U1','U2','U3','Theta','E')
    field_units=[];row_units=[]
    for i,name in enumerate(field_names):
        field=[v for s,v in known.items() if getattr(s,'__name__',None)=='s11cdCurrentPlus'+name]
        row=[v for s,v in known.items() if getattr(s,'name',None)=='s11cdPairingPlusRowResidual'+str(i)]
        require(len(field)==len(row)==1,'unique physical field and row units')
        field_units.append(tuple(field[0]));row_units.append(tuple(row[0]))
    j.value('row-field-unit-join.pickle',(field_units,row_units,units))
    require(field_units==list(map(tuple,units['fieldUnits'])) and row_units==list(map(tuple,units['rowUnits'])),'field and row orientation')
    inverse=j.op('restore-inverse',read_pickle,'inverse')
    determinant=j.op('restore-determinant',read_pickle,'determinant')
    cofactors=[[j.op('restore-cofactor-%d-%d'%(i,k),read_pickle,'cofactor_%d_%d'%(i,k)) for k in range(5)] for i in range(5)]
    fixed_inverse=j.op('bind-saved-inverse',lambda a,b:a.xreplace(b),inverse,source_bind)
    fixed_det=j.op('bind-saved-determinant',lambda a,b:a.xreplace(b),determinant,source_bind)
    fixed_adj=j.op('bind-saved-cofactors',lambda a,b:sp.ImmutableMatrix(a).xreplace(b),cofactors,source_bind)
    inverse_join=j.op('inverse-determinant-adjugate-join',
        lambda a,d,c:(a*d-c).applyfunc(sp.simplify),fixed_inverse,fixed_det,fixed_adj)
    require(inverse_join==sp.zeros(5),'actual integrated inverse and residue operands joined')

    # Use a positive algebraic coordinate for the already selected real-axis
    # radical. This rewrites accepted expressions, without a new determinant.
    radical=sp.Symbol('s11cdOutgoingPositiveRadical',positive=True)
    radicand=-square
    def algebraic_chart(value,variable,radicand,radical):
        prepared=sp.simplify(value);replacements={}
        for atom in prepared.atoms(sp.Pow):
            if atom.exp.is_Rational and atom.exp.q==2 and atom.base.has(variable):
                ratio=sp.cancel(atom.base/radicand)
                if not ratio.has(variable) and ratio.is_positive is True:
                    replacements[atom]=ratio**atom.exp*radical**(2*atom.exp)
        rational=sp.cancel(sp.together(prepared.xreplace(replacements)))
        return {'prepared':prepared,'radicalReplacements':replacements,'rational':rational,
            'rationalInCoordinates':rational.is_rational_function(variable,radical) is True and
                rational.free_symbols<=set((variable,radical)),
            'backJoin':sp.simplify(rational.xreplace({radical:sp.sqrt(radicand)})-value)}
    def denominator_certificate(rational,variable,radical,constant,quadratic):
        denominator=sp.fraction(sp.together(rational))[1]
        # Along r^2=a+b*k^2, an allowed denominator reduces to c*r^n.
        relation=sp.Poly(variable**2-(radical**2-constant)/quadratic,variable,domain='EX')
        remainder=sp.rem(sp.Poly(denominator,variable,domain='EX'),relation).as_expr()
        terms=sp.Poly(remainder,radical,domain='EX').terms()
        return {'denominator':denominator,'reducedOnBranch':remainder,'terms':terms,
            'allowed':not remainder.has(variable) and len(terms)==1 and terms[0][1].is_zero is False}
    denominator_summaries=[]
    for name,value in [('determinant',fixed_det)]+[
            ('adjugate-%d-%d'%(i,k),fixed_adj[i,k]) for i in range(5) for k in range(5)]:
        chart=j.op(name+'-algebraic-chart',algebraic_chart,value,kn,radicand,radical)
        require(chart['backJoin']==0 and chart['rationalInCoordinates'],'actual rational denominator chart identity: '+name)
        certificate=j.op(name+'-denominator-certificate',denominator_certificate,
            chart['rational'],kn,radical,positive.nth(0),positive.nth(2))
        summary={'name':name,'allowed':certificate['allowed'],
                 'reducedOnBranch':str(certificate['reducedOnBranch'])}
        save(out/(name+'-denominator.json'),summary);denominator_summaries.append(summary)
        require(certificate['allowed'],'only nonzero radical powers allowed in source denominators: '+name)

    # Both tails are Laurent expansions in positive t=1/abs(k). The radical
    # sqrt(a*t^2+b) is analytic near t=0 since b>0. Positive leading powers
    # imply decaying entries and integrable tail derivatives. Nondecay stops
    # this bounded separated-point construction before any kernel is emitted.
    reciprocal=sp.Symbol('s11cdOutgoingReciprocalMomentum',positive=True)
    def tail_data(rational,variable,radical,t,side,constant,quadratic):
        transformed=sp.cancel(rational.subs(
            {variable:side/t,radical:sp.sqrt(constant*t**2+quadratic)/t},simultaneous=True))
        leading=transformed.as_leading_term(t,cdir=1)
        coefficient,power=leading.as_coeff_exponent(t)
        series=sp.series(transformed,t,0,max(1,int(power)+1)) if leading!=0 and power.is_Integer else transformed
        return {'transformed':transformed,'leading':leading,'coefficient':coefficient,
            'power':power,'series':series,'identicallyZero':transformed==0,
            'decays':transformed==0 or (not coefficient.has(t) and power.is_Integer and power>0)}
    tails=[]
    for i in range(5):
        for k in range(5):
            chart=j.op('inverse-%d-%d-algebraic-chart'%(i,k),algebraic_chart,
                fixed_inverse[i,k],kn,radicand,radical)
            require(chart['backJoin']==0 and chart['rationalInCoordinates'],'actual rational inverse tail chart identity')
            for side in (-1,1):
                label='inverse-%d-%d-tail-%s'%(i,k,'minus' if side<0 else 'plus')
                tail=j.op(label,tail_data,chart['rational'],kn,radical,reciprocal,side,
                    positive.nth(0),positive.nth(2))
                summary={'row':i,'column':k,'side':side,'identicallyZero':bool(tail['identicallyZero']),
                    'powerOfReciprocalMomentum':str(tail['power']),'leadingCoefficient':str(tail['coefficient']),
                    'decays':bool(tail['decays'])}
                save(out/(label+'.json'),summary);tails.append(summary)
                require(tail['decays'],'nondecaying inverse tail needs explicit contact/distribution treatment; stop here')
    save(out/'large-momentum-summary.json',tails)
    records=modal['RECORDS'];native=modal['NATIVE_RECORDS']
    coverage=modal['NATIVE_COVERAGE'];reality=modal['NORMAL_REALITY_COVERAGE']
    j.value('inherited-real-axis-coverage.pickle',{
        'nativeCoverage':coverage,'normalReality':reality,
        'exactSpectrumResiduals':modal['EXACT_SPECTRUM_RECONSTRUCTION_RESIDUALS']})
    require(coverage['FINITE_POLYNOMIAL_ROOT_COVERAGE'] in (True,sp.true),'inherited finite polynomial root coverage')
    require(all(v==0 for v in reality['CHECKS'].values()) and
        all(v==0 for v in modal['EXACT_SPECTRUM_RECONSTRUCTION_RESIDUALS'].values()),'saved literal coverage residuals')
    disks=reality['DISKS']
    require(len(disks)==len(coverage['ROOT_DISKS'])==modal_cp['provenance']['rootDiskCount'],'complete inherited disk inventory')
    require(all(d['CERTIFIED_DISK_INPUT'] and d['DENOMINATOR_EXCLUDED'] and
        d['NORMAL_REALITY_STATUS'] in ('PROVED_REAL','PROVED_NONREAL') for d in disks),'no unresolved inherited normal-reality classification')
    require(len(records)==len(native)==modal_cp['recordCount'],'saved full candidate list')
    lifted={(int(r['ROOT_DISK_INDEX']),int(r['NORMAL_LIFT_SIGN'])) for r in records}
    require(len(lifted)==len(records) and lifted=={(i,s) for i in range(len(disks)) for s in (-1,1)},'complete unique inherited normal lifts')
    census=[];chosen=[]
    for record,old in zip(records,native):
        index=int(record['INDEX'])
        j.value('candidate-%d-inherited-operands.pickle'%index,{'record':record,'native':old})
        entry={'index':index,'onSheet':str(record['SHEET_MEMBERSHIP']),
            'realNormal':bool(record['EXACT_REAL_NORMAL']),'nullity':int(record['NULLITY'])}
        save(out/('candidate-%d-classification.json'%index),entry)
        require(all(int(record[k])==int(old[k]) for k in ('ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','NULLITY')),'saved candidate correspondence')
        certificate=disks[int(record['ROOT_DISK_INDEX'])]
        require(record['NORMAL_REALITY_CERTIFICATE']==certificate and
            bool(record['EXACT_REAL_NORMAL'])==(certificate['NORMAL_REALITY_STATUS']=='PROVED_REAL'),'candidate reality certificate join')
        require(record['SHEET_MEMBERSHIP'] in (True,False,sp.true,sp.false),'resolved inherited sheet classification')
        on_sheet=record['SHEET_MEMBERSHIP'] in (True,sp.true)
        real=bool(record['EXACT_REAL_NORMAL'])
        selected_mode=on_sheet and real
        census.append({'index':int(record['INDEX']),'onSheet':on_sheet,'realNormal':real,
                       'selected':selected_mode,'nullity':int(record['NULLITY'])})
        if real and not on_sheet:
            exact=record['EXACT_LOW_DEGREE_LIFT']
            require('K' in exact and 'Q' in exact and exact['REAL_K'] in (True,sp.true),'exact excluded real lift')
            k0=exact['K'];q0=j.op('excluded-%d-physical-q'%index,
                lambda a,b:sp.simplify(a.subs(kn,b)),q_out,k0)
            scale=modal['SCALAR_OPERANDS']['RADICAL_SCALE']
            scale_map={**binding(scale,(kr,qr,wr,kl,ql,wl)),kr:k0,kl:k0,qr:q0,
                ql:sp.conjugate(q0),wr:parameters['omega'],wl:parameters['omega']}
            scale0=j.op('excluded-%d-radical-scale'%index,
                lambda a,b:sp.simplify(a.xreplace(b)),scale,scale_map)
            opposite=j.op('excluded-%d-opposite-sheet-join'%index,
                lambda a,b,c:sp.simplify(a+b*c),q0,scale0,exact['Q'])
            require(opposite==0,'excluded real lift is exactly the opposite acoustic branch')
        if selected_mode:
            require(record['BULK_DECAY_DISK_CERTIFIED'] and record['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED'],'supported outgoing real block')
            exact=record['EXACT_LOW_DEGREE_LIFT']
            require('K' in exact and exact['REAL_K'] in (True,sp.true),'saved exact real pole coordinate')
            require(int(exact['MULTIPLICITY'])==int(record['NULLITY']) in (1,2),'bounded semisimple block multiplicity')
            chosen.append(record)
    save(out/'candidate-census.json',census)
    require(chosen,'nonempty supported real-axis block set')
    chosen.sort(key=lambda r:float(sp.re(r['EXACT_LOW_DEGREE_LIFT']['K'])))
    points=[r['EXACT_LOW_DEGREE_LIFT']['K'] for r in chosen]
    require(len(set(points))==len(points),'distinct exact block coordinates')
    order=max(int(r['NULLITY']) for r in chosen)
    det_derivatives=[fixed_det]+[j.op('determinant-derivative-%d'%n,lambda a,n:sp.diff(a,kn,n),fixed_det,n) for n in range(1,order+1)]
    adj_derivatives=[fixed_adj]+[j.op('cofactor-derivative-%d'%n,lambda a,n:a.diff(kn,n),fixed_adj,n) for n in range(1,order)]
    def numeric(a):
        return np.asarray(sp.ImmutableMatrix(a).evalf(50).tolist(),dtype=complex)
    def norm(a):
        return float(np.linalg.norm(a))
    summaries=[];blocks=[]
    for record,k0 in zip(chosen,points):
        index=int(record['INDEX']);n=int(record['NULLITY']);label='block-%d-'%index
        exact=record['EXACT_LOW_DEGREE_LIFT']
        j.value(label+'accepted-operands.pickle',record)
        dvalues=[j.op(label+'det-value-%d'%r,lambda a,b:sp.simplify(a.subs(kn,b)),a,k0) for r,a in enumerate(det_derivatives[:n+1])]
        avalues=[j.op(label+'adj-value-%d'%r,lambda a,b:a.subs(kn,b).applyfunc(sp.simplify),a,k0) for r,a in enumerate(adj_derivatives[:n])]
        require(all(v==0 for v in dvalues[:n]) and dvalues[n]!=0,'exact determinant zero order')
        require(all(v==sp.zeros(5) for v in avalues[:n-1]),'no higher-order inverse pole')
        residue=j.op(label+'exact-residue',lambda a,b:a.applyfunc(lambda v:sp.cancel(n*v/b)),avalues[n-1],dvalues[n])
        p0=j.op(label+'source-at-pole',lambda a,b:a.subs(kn,b).applyfunc(sp.simplify),fixed_source,k0)
        q0=j.op(label+'physical-q-at-pole',lambda a,b:sp.simplify(a.subs(kn,b)),q_out,k0)
        scale=modal['SCALAR_OPERANDS']['RADICAL_SCALE']
        scale_map={**binding(scale,(kr,qr,wr,kl,ql,wl)),kr:k0,kl:k0,qr:q0,ql:sp.conjugate(q0),wr:parameters['omega'],wl:parameters['omega']}
        scale0=j.op(label+'saved-radical-scale',lambda a,b:sp.simplify(a.xreplace(b)),scale,scale_map)
        native_join=j.op(label+'native-physical-q-join',sp.simplify,q0-scale0*exact['Q'])
        require(native_join==0,'exact native/physical acoustic coordinate join')
        forms=record['FORMS'];operands=record['OPERANDS']
        right,left=np.asarray(forms['RIGHT'],complex),np.asarray(forms['LEFT'],complex)
        require(right.shape==left.shape==(5,n),'complete saved left/right mode blocks')
        dk=np.asarray(operands['NORMAL_PENCIL_PLUS'],complex)
        dw=np.asarray(operands['FREQUENCY_PENCIL_PLUS'],complex)
        c=np.asarray(forms['N_NORMAL'],complex)
        sv=np.linalg.svd(c,compute_uv=False)
        projected=left.conj().T@dk@right
        j.value(label+'projected-operands.pickle',(right,left,dk,dw,c,projected,sv))
        require(c.shape==(n,n) and sv[-1]>1e-10*max(1.,sv[0]),'regular full projected derivative block')
        # A new small block contraction, not a repeated full-symbol LU solve.
        cinv=np.array([[1/c[0,0]]]) if n==1 else np.array([[c[1,1],-c[0,1]],[-c[1,0],c[0,0]]])/np.linalg.det(c)
        independent=right@cinv@left.conj().T
        a=numeric(residue);p=numeric(p0)
        flux_right=np.asarray(forms['FLUX_RIGHT'],complex);flux_left=np.asarray(forms['FLUX_LEFT'],complex)
        current=np.asarray(forms['SIGNED_CURRENT'],complex)
        diagonal=np.real(np.diag(current))
        require(np.all(diagonal>0) or np.all(diagonal<0),'uniform physical direction in whole block')
        direction=1 if np.all(diagonal>0) else -1
        projector=a@dk
        residuals={'normalCoordinateJoin':complex(sp.N(k0-record['K'],50)),
            'physicalBulkCoordinateJoin':complex(sp.N(q0-record['PHYSICAL_Q'],50)),
            'frequencyCoordinateJoin':complex(sp.N(parameters['omega']-record['OMEGA'],50)),
            'actualPencilJoin':p-np.asarray(operands['PENCIL_PLUS'],complex),
            'projectedDerivativeJoin':projected-c,'cofactorVsBlockResidue':a-independent,
            'leftLaurent':p@a,'rightLaurent':a@p,'derivativeProjector':projector@projector-projector,
            'fluxFrequencyNormalization':flux_left.conj().T@dw@flux_right-np.eye(n),
            'outgoingCurrentSign':direction*current-np.eye(n)}
        norms={key:norm(v) for key,v in residuals.items()}
        norms['cofactorVsBlockResidue']/=max(1.,norm(independent))
        mutation=(-direction)*current-np.eye(n)
        sign_delta=sp.I*sp.pi*direction*residue
        changed_delta=-sign_delta
        sign_movement=norm(numeric(changed_delta-sign_delta))
        j.value(label+'checks-and-mutation.pickle',{'residuals':residuals,'norms':norms,
            'exactResidue':residue,'independentBlockResidue':independent,'signedCurrent':current,
            'direction':direction,'mutatedDirection':-direction,'mutatedCurrentResidual':mutation,
            'deltaCoefficient':sign_delta,'mutatedDeltaCoefficient':changed_delta,'deltaMovementNorm':sign_movement})
        save(out/(label+'checks.json'),{'index':index,'nullity':n,'direction':direction,'norms':norms,
            'mutationNorm':norm(mutation),'deltaMovementNorm':sign_movement})
        require(all(v<1e-8 for v in norms.values()),'actual full-block residue/current joins')
        require(norm(mutation)>1 and sign_movement>1e-8,'one-sided direction mutation response')
        blocks.append({'index':index,'k':k0,'residue':residue,'direction':direction,'deltaCoefficient':sign_delta})
        summaries.append({'index':index,'nullity':n,'direction':direction,'norms':norms,
                          'deltaMovementNorm':sign_movement})
    # Explicit symmetric exclusions define the candidate principal-value part.
    # The branch is the established real-axis germ; no complex-frequency path
    # or global half-plane branch choice is manufactured here.
    state=selected['reductionState'];distance=state['z']-state['zp'];fourier=context['fourierMass']
    phase=sp.exp(sp.I*kn*distance)
    saved_density=j.op('restore-saved-spectral-density',read_pickle,'spectralDensity')
    density=j.op('bind-saved-spectral-density',lambda a,b:a.xreplace(b),saved_density,source_bind)
    density_join=j.op('spectral-density-phase-join',
        lambda a,b,c,d:(a-b.applyfunc(lambda v:c*v/d)).applyfunc(sp.simplify),density,fixed_inverse,phase,fourier)
    require(density_join==sp.zeros(5),'saved density with actual inverse, phase and Fourier measure')
    pointwise_domain=sp.Ne(state['z'],state['zp'],evaluate=False)
    domain={'pointwiseDomain':pointwise_domain,'normalMomentumDomain':sp.S.Reals,
        'radicalUnit':units['spectralMeasureUnit'],
        'reciprocalMomentumUnit':tuple(-v for v in units['spectralMeasureUnit']),
        'largeMomentumTails':tails,'denominatorCertificates':denominator_summaries,
        'polynomialContactPart':sp.zeros(5),'contactPartReason':'All 25 entries decay on both real tails; no polynomial large-momentum part.',
        'diagonalExcluded':True,'diagonalDistributionExtensionClaimed':False,
        'tailProof':'Rational functions of k and sqrt(a+b*k^2), a,b>0, have the saved analytic Laurent tails in t=1/abs(k); positive powers give decay and integrable derivatives.'}
    j.value('prescription-domain.pickle',domain)
    save(out/'prescription-domain.json',{'pointwiseDomain':str(pointwise_domain),'diagonalExcluded':True,
        'diagonalDistributionExtensionClaimed':False,'tailChecks':len(tails),'allTailsDecay':True,
        'denominatorChecks':len(denominator_summaries),'polynomialContactPart':'zero, from the checked decaying tails'})
    exclusion=sp.Symbol('s11cdOutgoingMomentumExclusion',positive=True)
    lower=[-sp.oo]+[k+exclusion for k in points]
    upper=[k-exclusion for k in points]+[sp.oo]
    intervals=tuple(zip(lower,upper))
    pv=sp.ImmutableMatrix(5,5,lambda i,k:sp.Limit(sp.Add(*(
        sp.Integral(density[i,k],(kn,a,b)) for a,b in intervals)),exclusion,0,dir='+'))
    correction=sp.zeros(5)
    distribution_correction=sp.zeros(5)
    for block in blocks:
        correction+=sp.exp(sp.I*block['k']*distance)*block['deltaCoefficient']/fourier
        distribution_correction+=block['deltaCoefficient']*sp.DiracDelta(kn-block['k'])
    # Residue-theorem orientation for exp(+ik x): upper closure at x>0,
    # lower closure with reversed orientation at x<0. Exclude x=0 contact terms.
    support=[]
    for block in blocks:
        for side in (-1,1):
            from_pv_and_delta=sp.I*sp.pi*(side+block['direction'])/fourier
            from_outgoing_closure=2*sp.I*sp.pi*side/fourier if side==block['direction'] else sp.S.Zero
            support.append({'block':block['index'],'spatialSide':side,
                'pvPlusDelta':from_pv_and_delta,'outgoingClosure':from_outgoing_closure,
                'residual':sp.simplify(from_pv_and_delta-from_outgoing_closure)})
    j.value('spatial-sign-identity.pickle',support)
    require(all(r['residual']==0 for r in support),'Fourier contour orientation identity')
    j.value('outgoing-prescription-candidate.pickle',{
        'fixedPhysicalInput':physical,'fixedSymbol':fixed_source,'fixedInverse':fixed_inverse,
        'regularDensity':density,'physicalBulkBranch':q_out,'bulkWaveRelation':fixed_wave,
        'blocks':blocks,'momentumDeltaCorrection':distribution_correction,
        'principalValueExclusionIntervals':intervals,'exclusionMomentumUnit':units['spectralMeasureUnit'],
        'principalValueKernel':pv,'singularKernelCorrection':correction,
        'candidateKernel':pv+correction,'kernelEntryUnits':units['kernelEntryUnits'],
        'pointwiseDomain':pointwise_domain,'domainEvidence':domain,
        'referenceBinding':source_bind,'effectiveReferenceGrades':{'eta_bg':sp.S.Zero,'sigma_W':sp.S.Zero},
        'residueEntryUnits':units['kernelEntryUnits'],'profileAbelRegulator':context['profileRegulator'],
        'candidateScope':'Fixed positive-real-frequency input, selected physical sheet, real normal line; symmetric PV plus current-oriented local pole contributions.',
        'domainReviewOutstanding':['Acceptance of actual saved residuals and tail/domain evidence',
            'Coincident-position distributional extension remains outside this separated-point result'],
        'complexFrequencyRetardedEquivalenceClaimed':False,'fullFORM':False})
    status={'status':'FIXED_INPUT_OUTGOING_PRESCRIPTION_CANDIDATE_BUILT',
        'case':'LAB_HELD__RHO4_CONSTANT','blocks':summaries,'completedOperations':len(j.records),
        'pointwiseDomain':str(pointwise_domain),'tailChecks':len(tails),'allTailsDecay':True,
        'denominatorChecks':len(denominator_summaries),'inverseAdjugateJoinZero':True,'savedSpectralDensityJoinZero':True,
        'newModeSolves':0,'newCurrentConstructions':0,'newRootSearches':0,
        'realAxisBranch':'Source-joined evanescent branch at the unchanged physical input',
        'outgoingGreenOperatorAccepted':False,'complexFrequencyRetardedEquivalenceClaimed':False,
        'formCompleted':False,'radiatingCoverage':False,'resultMethodAcceptancePending':True}
    post={key:route(r['path']) for key,r in actual.items()}
    save(out/'posthashes.json',post);require(post==actual,'consumed inputs unchanged')
    save(out/'operation-index.json',j.records)
    return status


def main():
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--input-manifest',type=Path,required=True)
    p.add_argument('--gate-receipt',type=Path,required=True)
    p.add_argument('--run-directory',type=Path,required=True)
    args=p.parse_args();observed=containment()
    spec=json.loads(args.input_manifest.read_text());gate=json.loads(args.gate_receipt.read_text())
    require(gate['status']=='READY_FOR_ONE_GUARDED_OUTGOING_PRESCRIPTION_STAGE','method gate incomplete')
    require(gate['workerSha256']==digest(Path(__file__)) and gate['inputManifestSha256']==digest(args.input_manifest),'method gate source/input pins')
    require(gate['scopeAuthorization']=='Continue','bounded continuation authorization')
    out=args.run_directory.resolve();out.relative_to(STORE);out.mkdir(parents=True,exist_ok=False)
    save(out/'native-containment.json',observed);save(out/'input-manifest.json',spec);save(out/'gate-receipt.json',gate)
    started=time.monotonic()
    try:
        result=construct(spec,out)
    except BaseException:
        save(out/'failure.json',{'traceback':traceback.format_exc(),'automaticRetry':False,'wallSeconds':time.monotonic()-started})
        raise
    result['wallSeconds']=time.monotonic()-started
    save(out/'checks.json',result);print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
