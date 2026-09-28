#!/usr/bin/env python3
"""Guarded saved-source threshold/channel/face pilot. No Born or producer import.

This is an unexecuted build-review candidate. Exact loci are not automatically
physical cutoffs; grid classifications remain local. Every scientific operation
saves its input and complete return before acceptance checks.
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

ROOT = Path('/var/projects/toy_physics')
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def require(test, reason):
    if not test:
        raise ValueError(reason)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def route(path):
    p = Path(path)
    return dict(path=str(p), canonicalPath=str(p.resolve(strict=True)),
                bytes=p.stat().st_size, sha256=sha(p))


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


class OperationBudget(Exception):
    pass


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.')) or
                module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError((module, name))
        return super().find_class(module, name)


class Journal:
    def __init__(self, out, deadline):
        self.out, self.deadline, self.records, self.active = out, deadline, [], None

    def op(self, name, function, *args, seconds=None):
        p = self.out/'operations'/('%04d-%s' % (len(self.records), name))
        p.mkdir(parents=True, exist_ok=False)
        with (p/'input.pickle').open('xb') as f:
            pickle.dump(args, f, protocol=4); f.flush(); os.fsync(f.fileno())
        self.active = dict(name=name, input=route(p/'input.pickle'),
                           startedUtc=datetime.now(timezone.utc).isoformat())
        save(p/'started.json', self.active)
        remaining = self.deadline-time.monotonic()
        require(remaining > 0, 'native whole-job deadline reached')
        if seconds is not None:
            signal.setitimer(signal.ITIMER_REAL, min(seconds, remaining))
        try:
            value = function(*args)
            with (p/'value.pickle').open('xb') as f:
                pickle.dump(value, f, protocol=4); f.flush(); os.fsync(f.fileno())
            self.active.update(status='COMPLETE', value=route(p/'value.pickle'))
            save(p/'completed.json', self.active)
            self.records.append(self.active); self.active = None
            return value
        except OperationBudget:
            if seconds is None or time.monotonic() >= self.deadline:
                raise
            self.active.update(status='OPTIONAL_LOCUS_BUDGET_EXHAUSTED_NO_RETRY')
            save(p/'budget-stop.json', self.active)
            self.records.append(self.active); self.active = None
            return {'status': 'LOCUS_UNRESOLVED_WITHIN_BUDGET'}
        finally:
            signal.setitimer(signal.ITIMER_REAL, max(.001, self.deadline-time.monotonic()))


def containment():
    group = next(x[3:] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::'))
    c = Path('/sys/fs/cgroup')/group.lstrip('/')
    limits = {n: (c/n).read_text().strip() for n in ('memory.max', 'memory.swap.max', 'pids.max')}
    limits.update(nice=os.getpriority(os.PRIO_PROCESS, 0), affinity=sorted(os.sched_getaffinity(0)),
                  threads={n: os.environ.get(n) for n in THREADS})
    require(limits['memory.max']=='2147483648' and limits['memory.swap.max']=='0'
            and limits['pids.max']=='32' and limits['nice']>=15 and len(limits['affinity'])==1
            and all(v=='1' for v in limits['threads'].values()), 'missing containment; no science imported')
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def timeout(*_):
        raise OperationBudget('bounded operation/native deadline; preserve work')
    signal.signal(signal.SIGALRM, timeout); signal.setitimer(signal.ITIMER_REAL, 840)
    return limits


def construct(spec, out, journal):
    import sympy as sp
    import numpy as np

    w = sp.Symbol('mapFrequency', positive=True)
    v = sp.Symbol('mapTangentialRayCoordinate', nonnegative=True)
    k = sp.Symbol('mapNormalMomentum', real=True)
    q, qb = sp.symbols('mapPhysicalDepthMomentum mapLeftDepthMomentum', complex=True)
    params = {n: sp.Rational(a) for n, a in spec['physicalInput']['parameters'].items()}

    def show(x):
        if isinstance(x, dict): return {str(a): show(b) for a,b in x.items()}
        if isinstance(x, (list, tuple)): return [show(a) for a in x]
        if isinstance(x, sp.MatrixBase): return [[str(a) for a in row] for row in x.tolist()]
        if isinstance(x, np.ndarray): return show(x.tolist())
        if isinstance(x, complex): return dict(real=x.real, imag=x.imag)
        if isinstance(x, (np.integer, np.floating, np.bool_)): return x.item()
        if isinstance(x, sp.Basic): return str(x)
        return x

    def op(name, fn, *args, seconds=None):
        result = journal.op(name, fn, *args, seconds=seconds)
        save(out/(name+'.json'), show(result))
        return result

    def restore(name):
        record = spec['inputs'][name]
        require(route(record['path'])==record, 'input changed: '+name)
        with Path(record['path']).open('rb') as f: return SavedCodec(f).load()

    def flat(x):
        if isinstance(x, sp.MatrixBase): return list(x)
        if isinstance(x, (tuple, list)): return sum((flat(a) for a in x), [])
        return [x]

    def zeros(x): return all(sp.cancel(a)==0 for a in flat(x))

    def bind(x, label, extra=None):
        mapping = dict(extra or {})
        for s in x.free_symbols-set(mapping):
            if s.name in ('omega', 's11cdFrequency'): mapping[s]=w
            elif s.name=='s11cdTangentialMomentum1': mapping[s]=2*v
            elif s.name=='s11cdTangentialMomentum2': mapping[s]=v
            elif s.name=='s11cdSpectralNormalMomentum': mapping[s]=k
            elif s.name=='eta_bg': mapping[s]=sp.S.Zero if label=='REFERENCE' else params['eta_bg']
            elif s.name=='sigma_W': mapping[s]=sp.S.Zero if label=='REFERENCE' else params['eta_bg']*params['W_0']/params['L_W']
            elif s.name in params: mapping[s]=params[s.name]
            elif s in (w,v,k,q,qb): pass
            else: raise ValueError('unbound source symbol '+s.name)
        return x.xreplace(mapping)

    def coefficient(x):
        eps = [s for s in x.free_symbols if s.name=='epsilon_shape']
        require(len(eps)==1, 'saved quadratic shape-amplitude carrier')
        e=eps[0]; y=x.applyfunc(lambda a: sp.Poly(a,e).nth(2))
        return y, x-e**2*y

    uniform=op('restore-uniform', restore, 'uniform')
    reduction=op('restore-reduced-source-context',restore,'reduction')
    branch_state=reduction['reductionState']
    branch_inventory={'sourceEquations':branch_state['branch_equations'],
                      'sourceMap':branch_state['branch_map'],
                      'momentumGroups':branch_state['momentum_groups'],
                      'tangents':branch_state['tangents'],
                      'normalMap':branch_state['normal_map']}
    op('source-branch-binding-inventory',lambda x:x,branch_inventory)
    face_index=json.loads(Path(spec['inputs']['faceIndex']['path']).read_text())
    units=op('restore-unit-context', restore, 'unitContext')
    unit_record={'unitFrame':spec['physicalInput']['unit_frame'],
                 'fieldReferenceUnits':units['fieldUnits'], 'currentUnit':units['currentUnit'],
                 'sourceStrongEntryUnits':uniform['units']['strong'],
                 'convention':'Numeric matrices are coefficients in these saved field/row reference units; current is not a field norm.'}
    save(out/'unit-context.json', show(unit_record))

    def bound_end(label, packet, pair_tuple):
        pair, known=pair_tuple; result=pair['result']
        wl,wr=result['FREQUENCY_LEGS']; kl,kr=result['NORMAL_LEGS']; ql,qr=result['BULK_LEGS']
        source_map={wl:w,wr:w,kl:k,kr:k,ql:qb,qr:q}
        M=bind(result['CLOSED_PENCIL_LEGS'][0],label,source_map).applyfunc(sp.cancel)
        wave=sp.expand(bind(result['ACOUSTIC_WAVE_ROWS'][0],label,source_map))
        wavepoly=sp.Poly(wave,q)
        supported=wavepoly.degree()==2 and wavepoly.nth(1)==0 and wavepoly.nth(2)!=0
        if not supported: return {'supported':False,'matrix':M,'wave':wave}
        q2=sp.cancel(-wavepoly.nth(0)/wavepoly.nth(2))
        oldq=next(s for s in packet['originalRelation'].free_symbols if s.name=='s11cdBulkRadical')
        O=bind(packet['originalAlgebraic'],label,{oldq:oldq})
        W=bind(packet['originalRelation'],label,{oldq:oldq})
        scale=sp.sqrt(sp.cancel(sp.Poly(W,oldq).nth(2)/wavepoly.nth(2)))
        join=(M.subs(q,scale*oldq)-O).applyfunc(sp.cancel)
        wavejoin=sp.cancel(wave.subs(q,scale*oldq)-W)
        C=bind(uniform['curl'],label)
        # The source trial ansatz is curl(A), grad(Phi), Theta and E.
        T=sp.ImmutableMatrix.vstack(C[:,1:3],sp.zeros(2,2))
        H=sp.ImmutableMatrix([[0,0,sp.I*2*v],[0,0,sp.I*v],[0,0,sp.I*k],[1,0,0],[0,1,0]])
        B=sp.ImmutableMatrix.hstack(T,H)
        # Kinematic coordinate inverse only, never the end-pencil inverse.
        dual=B.inv().applyfunc(sp.cancel)
        separated=(dual*M*B).applyfunc(sp.cancel)
        slab,slab_grade=coefficient(result['SLAB_CURRENT_MATRIX'])
        bulk,bulk_grade=coefficient(result['BULK_NORMAL_CURRENT_DENSITY_MATRIX'])
        slab=bind(slab,label,source_map); bulk=bind(bulk,label,source_map)
        amps=sorted({s for s in known if getattr(s,'name','').startswith('s11cdCurrentPlusAmplitude')},key=lambda s:s.name)
        require(len(amps)==5, 'five actual face amplitude carriers')
        e=next(s for s in result['OPEN_BULK_CURRENT_COEFFICIENTS'].free_symbols if s.name=='epsilon_shape')
        free_bulk_amplitudes={s:sp.S.One for s in result['OPEN_BULK_CURRENT_COEFFICIENTS'].free_symbols if s.name in ('s11cdAcousticLeftAmplitude','s11cdAcousticRightAmplitude')}
        open_depth=sp.Poly(result['OPEN_BULK_CURRENT_COEFFICIENTS'][3],e).nth(2)
        open_depth=bind(open_depth,label,{**source_map,**free_bulk_amplitudes})
        face_rows=[]
        for face in result['FACE_LEG_OBJECTS']:
            rows={}
            for key in ('AMPLITUDE','PRESSURE','OUTWARD_VELOCITY','RELATIVE_MASS_FLUX','AFFINITY','BULK_VELOCITY'):
                expr=face[0][key]
                row=sp.ImmutableMatrix(1,5,lambda i,j:sp.diff(expr,amps[j]))
                rows[key]=bind(row,label,source_map)
                rows[key+'_linearResidual']=bind(expr-(row*sp.ImmutableMatrix(amps))[0],label,{**source_map,**{a:a for a in amps}})
            face_rows.append(rows)
        coupling=[uniform['records'][label]['coupling'][name] for name in ('TH','HT')]
        q2poly=sp.Poly(q2,k)
        bulkmax=q2poly.nth(0) if q2poly.degree()==2 and q2poly.nth(1)==0 and q2poly.nth(2).is_negative is True else None
        return {'supported':True,'onlyDeclaredLiveSymbols':not (M.free_symbols|wave.free_symbols)-{w,v,k,q},
                'liveFrequencyAndTangential':{w,v}.issubset(M.free_symbols|wave.free_symbols),'end':label,'matrix':M,'wave':wave,'depthSquared':q2,'bulkMaximum':bulkmax,
                'radicalScale':scale,'frequencySourceJoin':join,'physicalWaveJoin':wavejoin,
                'sourceBranchJoins':packet['branchResiduals'],'originalSourceCoupling':coupling,
                'curl':C,'lift':B,'dual':dual,'chartDeterminant':sp.factor(B.det()),
                'chartJoin':dual*B-sp.eye(5),'separated':separated,
                'TH':separated[:2,2:],'HT':separated[2:,:2],
                'slabCurrent':slab,'bulkCurrentDensity':bulk,'shapeGradeResiduals':(slab_grade,bulk_grade),
                'faceRows':face_rows,'faceAmplitudes':amps,'openDepthCurrentCoefficient':open_depth,'sourceKnownDimensions':known,
                'sourceStrongUnits':uniform['units']['strong']}

    def face_test(b):
        T=b['lift'][:,:2]; rows=[]
        literal=[]
        for item in face_index['records']:
            if item['tag']=='PY_S11CC2_FOLD_SYMBOL_MAP_LAB_HELD_RHO4_CONSTANT':
                for ident in item['velocityIdentifications']:
                    expr=sp.sympify(ident['savedValueSrepr'])
                    field=next(s for s in expr.free_symbols if s.name=='e_W_t')
                    velocity=bind(expr,b['end'],{field:-sp.I*w})
                    literal.append({'face':item['face'],'source':expr,'coefficient':velocity,'sourceLiteralSha256':ident['literalSha256']})
        for i,face in enumerate(b['faceRows']):
            contractions={name:(face[name]*T).applyfunc(sp.cancel) for name in ('AMPLITUDE','PRESSURE','OUTWARD_VELOCITY','RELATIVE_MASS_FLUX','AFFINITY','BULK_VELOCITY')}
            original=face['OUTWARD_VELOCITY']; changed=original.copy(); changed[0,4]=0
            probe=sp.ImmutableMatrix([0,0,0,0,1])
            before=original*probe; after=changed*probe; difference=(after-before).applyfunc(sp.cancel)
            rows.append({'faceIndex':i,'nativeRows':face,'transverseLift':T,'contractions':contractions,
                         'referenceDriveZero':zeros(tuple(contractions.values())),
                         'savedC2VelocityIdentifications':literal,
                         'c2ReferenceVelocityJoins':[sp.cancel(original[0,4]-x['coefficient']) for x in literal] if b['end']=='REFERENCE' else [],
                         'control':{'omittedNativeSlot':'OUTWARD_VELOCITY/eW','probe':probe,'originalRow':original,'changedRow':changed,
                                    'baseline':before,'mutated':after,'difference':difference,
                                    'nonzeroFormalResponse':not zeros(difference)}})
        return {'end':b['end'],'faces':rows,'scope':'Transverse-sector face drive of this named end; REFERENCE supplies the requested reference test. No first-order map/tilt/induced-field or bulk-power conclusion.'}

    def determinant(matrix):
        ds=[sp.lcm([sp.denom(sp.cancel(a)) for a in matrix.row(i)]) for i in range(matrix.rows)]
        cleared=sp.ImmutableMatrix(matrix.rows,matrix.cols,lambda i,j:sp.cancel(ds[i]*matrix[i,j]))
        numerator,denominator=sp.fraction(sp.cancel(cleared.det(method='domain-ge')/sp.prod(ds)))
        return {'matrix':matrix,'cleared':cleared,'rowDenominators':ds,'numerator':numerator,'denominator':denominator}

    def eliminate(det,b):
        result=sp.expand(sp.resultant(det['numerator'],b['wave'],q))
        re,im=map(sp.expand,result.as_real_imag())
        field=sp.QQ.frac_field(w,v)
        R=sp.Poly(re,k,domain=field); I=sp.Poly(im,k,domain=field)
        if R.is_zero and I.is_zero:
            return {'supported':False,'resultant':result,'real':re,'imag':im}
        if I.is_zero: G=R.monic(); A=sp.Poly(1/R.LC(),k,domain=field); D=sp.Poly(0,k,domain=field)
        elif R.is_zero: G=I.monic(); A=sp.Poly(0,k,domain=field); D=sp.Poly(1/I.LC(),k,domain=field)
        else: A,D,G=sp.gcdex(R,I)
        qr,rr=sp.div(R,G);qi,ri=sp.div(I,G)
        witness=sp.cancel((A*R+D*I-G).as_expr())
        denoms=[sp.denom(x) for P in (R,I,A,D,G,qr,qi) for x in P.all_coeffs()]
        return {'supported':True,'resultant':result,'real':re,'imag':im,'gcd':G.as_expr(),
                'bezoutA':A.as_expr(),'bezoutB':D.as_expr(),'bezoutResidual':witness,
                'divisionReal':rr.as_expr(),'divisionImag':ri.as_expr(),
                'genericSpecializationExclusions':denoms,
                'scope':'Generic gcd is not used to discard exceptional parameter curves; per-point gcd is recomputed from the original real/imag polynomials.'}

    def loci(elim,b,det):
        G=sp.Poly(elim['gcd'],k).sqf_part();g=G.as_expr()
        named={'zeroNormalCandidate':sp.together(g.subs(k,0)),
               'genericLeadingCoefficient':G.LC(),
               'branchCollision':sp.resultant(g,sp.fraction(b['depthSquared'])[0],k) if G.degree()>0 else sp.S.One,
               'rootCollision':sp.discriminant(g,k) if G.degree()>1 else sp.S.One}
        R=sp.Poly(elim['real'],k);I=sp.Poly(elim['imag'],k)
        r=sp.cancel(R.as_expr()/elim['gcd']);i=sp.cancel(I.as_expr()/elim['gcd'])
        if r!=0 and i!=0: named['exceptionalRealImagCommonRoot']=sp.resultant(sp.fraction(r)[0],sp.fraction(i)[0],k)
        source_den=sp.prod([det['denominator'],*det['rowDenominators']])
        den_locus=sp.resultant(sp.fraction(sp.cancel(source_den))[0],b['wave'],q)
        if G.degree()>0:
            named['sourcePoleCollision']=sp.resultant(g,den_locus,k)
            named['chartCollision']=sp.resultant(g,b['chartDeterminant'],k)
        named.update({'gcdWitnessDenominator'+str(i):a for i,a in enumerate(elim['genericSpecializationExclusions']) if a!=1})
        facts={n:sp.factor(sp.fraction(sp.cancel(x))[0]) for n,x in named.items()}
        return {'status':'EXACT_CANDIDATE_LOCUS_EQUATIONS_NOT_PHYSICAL_CUTOFF_CLEARANCE',
                'equations':facts,'raw':named,'coordinateRelation':'kappa^2=5*v^2; v>=0',
                'interpretation':'Physical-sheet/current tests and coverage determine which candidates bound actual open cells. No blanket cell-absence claim.'}

    def root_census(elim,b,w0,v0):
        sub={w:w0,v:v0};R=sp.Poly(elim['real'].subs(sub),k,domain=sp.QQ);I=sp.Poly(elim['imag'].subs(sub),k,domain=sp.QQ)
        if R.is_zero and I.is_zero:return {'status':'DEGENERATE_ELIMINATION_UNRESOLVED','roots':[]}
        G=sp.gcd(R,I).sqf_part()
        roots=G.intervals(eps=sp.Rational(1,10**35)) if G.degree()>0 else []
        return {'status':'COMPLETE_NECESSARY_REAL_ROOT_CENSUS','real':R.as_expr(),'imag':I.as_expr(),'gcd':G.as_expr(),
                'realRemainder':sp.rem(R,G).as_expr(),'imagRemainder':sp.rem(I,G).as_expr(),
                'roots':roots,'count':G.count_roots(-sp.oo,sp.oo),'frequency':w0,'rayCoordinate':v0,
                'depthSquared':b['depthSquared'].subs(sub)}

    def numeric(x,sub):
        return np.array(x.subs(sub).evalf(40).tolist(),dtype=complex)

    def inspect_root(b,det,census,interval,sector):
        (lo,hi),multiplicity=interval;km=(lo+hi)/2
        row={'interval':interval,'multiplicity':multiplicity,'sector':sector,'status':'UNRESOLVED','currentAccepted':False}
        q2p=sp.Poly(sp.fraction(census['depthSquared'])[0],k,domain=sp.QQ)
        if lo<=0<=hi or q2p.count_roots(lo,hi):return {**row,'status':'THRESHOLD_OR_BRANCH_ENDPOINT'}
        q2=sp.cancel(census['depthSquared'].subs(k,km))
        if q2.is_real is not True or q2==0:return {**row,'status':'DEPTH_ROOT_UNRESOLVED'}
        q0=sp.sqrt(q2); sub={w:census['frequency'],v:census['rayCoordinate'],k:km,q:q0,qb:sp.conjugate(q0)}
        denominators=[x.subs(sub) for x in [det['denominator'],*det['rowDenominators'],b['chartDeterminant']]]
        if any(abs(complex(a.evalf(40)))<1e-12 for a in denominators):return {**row,'status':'POLE_OR_CHART_UNRESOLVED','denominators':denominators}
        M=numeric(b['matrix'],sub);scale=np.maximum(np.linalg.norm(M,axis=1),1.0);A=M/scale[:,None]
        U,S,Vh=np.linalg.svd(A);null=S<1e-10
        depth_flux=b['openDepthCurrentCoefficient'].subs(sub)
        row.update(outgoingDepthCurrentCoefficient=depth_flux,physicalDepthMomentum=q0,depthSquared=q2,originalMatrix=M,rowScale=scale,singularValues=S,denominators=denominators)
        if q2.is_positive is True and not (abs(complex(depth_flux.evalf(40)).imag)<1e-10 and complex(depth_flux.evalf(40)).real>0):
            return {**row,'status':'OUTGOING_DEPTH_CURRENT_GERM_UNRESOLVED'}
        if not null.any():return {**row,'status':'REJECTED_SELECTED_SHEET'}
        if np.any((S>=1e-10)&(S<1e-7)):return {**row,'status':'AMBIGUOUS_RANK'}
        V=Vh.conj().T[:,null];left=U[:,null]/scale[:,None]
        # Source-derived lift/dual distinguishes the full T/H subspaces.
        dual=numeric(b['dual'],sub);coords=dual@V
        wanted=coords[:2] if sector=='T' else coords[2:]
        unwanted=coords[2:] if sector=='T' else coords[:2]
        residual=A@V; adjoint=M.conj().T@left
        row.update(rightBasis=V,leftBasis=left,fieldCoordinates=coords,rightResidual=residual,leftResidual=adjoint)
        if np.linalg.norm(unwanted)>1e-7*max(1.,np.linalg.norm(coords)):
            return {**row,'status':'MIXED_OR_COINCIDENT_SECTORS_UNRESOLVED'}
        if np.linalg.norm(residual)>1e-8 or np.linalg.norm(adjoint)>1e-8*max(1.,np.linalg.norm(M)):
            return {**row,'status':'FAILED_ORIGINAL_MATRIX_RESIDUAL'}
        dqdw=sp.cancel(sp.diff(b['depthSquared'],w)/(2*q))
        dqdk=sp.cancel(sp.diff(b['depthSquared'],k)/(2*q))
        N=left.conj().T@numeric(b['matrix'].diff(w)+b['matrix'].diff(q)*dqdw,sub)@V
        Nk=left.conj().T@numeric(b['matrix'].diff(k)+b['matrix'].diff(q)*dqdk,sub)@V
        row.update(frequencyPairing=N,normalPairing=Nk)
        if np.linalg.matrix_rank(N,tol=1e-9)<V.shape[1] or np.linalg.matrix_rank(Nk,tol=1e-9)<V.shape[1]:
            return {**row,'status':'SINGULAR_MODAL_PAIRING_UNRESOLVED'}
        face=[{n:numeric(f[n],sub)@V for n in ('AMPLITUDE','PRESSURE','OUTWARD_VELOCITY')} for f in b['faceRows']]
        Js=numeric(b['slabCurrent'],sub);Jb=numeric(b['bulkCurrentDensity'],sub)
        Jbulk=V.conj().T@Jb@V
        if q2.is_negative is True:
            depth=sp.I/(q0-sp.conjugate(q0));J=V.conj().T@(Js+complex(depth.evalf(40))*Jb)@V
            domain='DECAYING_EXTERIOR_WITH_INFINITE_DEPTH_CURRENT'
        elif sector=='T' and b['transverseFaceZero']:
            depth=None;J=V.conj().T@Js@V;domain='EXACTLY_UNDRIVEN_REFERENCE_SECTOR_BULK_TERM_ZERO'
            if np.linalg.norm(Jbulk)>1e-8:return {**row,'status':'ZERO_DRIVE_CURRENT_JOIN_FAILED','face':face,'bulkCurrent':Jbulk}
        else:
            return {**row,'status':'RADIATING_UNIFORM_MODE_CURRENT_UNSUPPORTED','face':face,'bulkCurrent':Jbulk}
        hermitian=J-J.conj().T;evals,rotation=np.linalg.eigh((J+J.conj().T)/2)
        frozen=b['matrix'].copy();frozen[0,0]=frozen[0,0].subs(w,0)
        Mbad=numeric(frozen,sub);bad=(Mbad/scale[:,None])@V
        shift=np.linalg.norm(bad-residual)
        tangent_frozen=b['matrix'].copy();tangent_frozen[0,0]=tangent_frozen[0,0].subs(v,0)
        tangent_bad=numeric(tangent_frozen,sub);tangent_res=(tangent_bad/scale[:,None])@V
        tangent_shift=np.linalg.norm(tangent_res-residual)
        row.update(face=face,slabCurrentMatrix=Js,bulkCurrentMatrix=Jb,depthIntegral=depth,currentForm=J,
                   currentHermiticityResidual=hermitian,currentEigenvalues=evals,currentRotation=rotation,
                   currentDomain=domain,mutation={'type':'OMIT_ACTUAL_FREQUENCY_DEPENDENCE_IN_SOURCE_ENTRY_0_0',
                   'tangentialFrozenMatrix':tangent_bad,'tangentialFrozenResidual':tangent_res,'tangentialMovement':float(tangent_shift),
                   'baselineMatrix':M,'mutatedMatrix':Mbad,'baselineResidual':residual,'mutatedResidual':bad,'movement':float(shift)})
        if np.linalg.norm(hermitian)>1e-8*max(1.,np.linalg.norm(J)) or np.min(np.abs(evals))<1e-9*max(1.,np.linalg.norm(J)):
            return {**row,'status':'CURRENT_RANK_OR_REALITY_UNRESOLVED'}
        # The selected physical source mutation must respond, not a +1 wiring change.
        if sector=='T' and min(shift,tangent_shift)<=1e-8:return {**row,'status':'SOURCE_FORM_CONTROL_UNRESPONSIVE'}
        return {**row,'status':'OPEN_'+sector+'_CURRENT_CHANNEL','currentAccepted':True,
                'positiveCurrentRank':int(sum(evals>0)),'negativeCurrentRank':int(sum(evals<0))}

    ends={}; spectral={}; curves={}; faces={}
    for label in ('REFERENCE','LEFT','RIGHT'):
        packet=op(label+'-restore-frequency-source',restore,label+'Frequency')
        pair=op(label+'-restore-pairing-source',restore,label+'Pairing')
        b=op(label+'-bind-and-source-join',bound_end,label,packet,pair)
        require(b['supported'] and b['onlyDeclaredLiveSymbols'] and b['liveFrequencyAndTangential'] and zeros(b['frequencySourceJoin']) and zeros(b['physicalWaveJoin'])
                and zeros(b['sourceBranchJoins']) and zeros(b['originalSourceCoupling'])
                and zeros(b['chartJoin']) and zeros(b['TH']) and zeros(b['HT'])
                and zeros(b['shapeGradeResiduals'])
                and all(zeros(f[n+'_linearResidual']) for f in b['faceRows'] for n in ('AMPLITUDE','PRESSURE','OUTWARD_VELOCITY','RELATIVE_MASS_FLUX','AFFINITY','BULK_VELOCITY')), 'actual source, branch, sector and current-grade joins')
        ft=op(label+'-face-test',face_test,b)
        require(all(x['control']['nonzeroFormalResponse'] and zeros(x['c2ReferenceVelocityJoins']) for x in ft['faces']), 'actual native face-term control')
        b['transverseFaceZero']=all(x['referenceDriveZero'] for x in ft['faces']);faces[label]=ft
        bulk={'sourceWave':b['wave'],'depthSquared':b['depthSquared'],'maximumOverRealNormal':b['bulkMaximum'],
              'thresholdEquation':sp.factor(b['bulkMaximum']) if b['bulkMaximum'] is not None else None,
              'normalMomentumDomain':sp.Gt(b['depthSquared'],0),'coordinateRelation':'kappa^2=5*v^2',
              'status':'EXACT_BULK_KINEMATIC_BOUNDARY' if b['bulkMaximum'] is not None else 'BULK_BOUNDARY_UNRESOLVED',
              'powerComputed':False}
        op(label+'-bulk-boundary',lambda x:x,bulk)
        ends[label]=b;spectral[label]={};curves[label]={}
        for sector,block in [('T',b['separated'][:2,:2]),('H',b['separated'][2:,2:])]:
            det=op(label+'-'+sector+'-determinant',determinant,block)
            elim=op(label+'-'+sector+'-elimination',eliminate,det,b)
            require(elim['supported'] and zeros([elim['bezoutResidual'],elim['divisionReal'],elim['divisionImag']]),'exact elimination witnesses')
            lc=op(label+'-'+sector+'-candidate-loci',loci,elim,b,det,seconds=45)
            spectral[label][sector]=(det,elim);curves[label][sector]=lc
    # Main result is exact loci plus qualifications; grid is a bounded cross-check.
    save(out/'threshold-loci.json',show(curves));save(out/'reference-face-drive.json',show(faces['REFERENCE']))
    frequencies=[sp.Rational(x) for x in spec['grid']['frequencies']]
    kappas=[sp.sympify(x) for x in spec['grid']['magnitudes']]
    fixed=sp.sqrt(sp.Rational(1,20))
    points=[(x,fixed) for x in frequencies]+[(x,y) for x in frequencies for y in kappas if y!=fixed]
    summaries=[];timings=[];grid_start=time.monotonic();stopped=None
    for point_index,(w0,kappa) in enumerate(points):
        if time.monotonic()>journal.deadline-75:
            stopped='WHOLE_JOB_BUDGET_SAVED_PREFIX';break
        if point_index>=12 and timings and (sum(timings)/len(timings))*(len(points)-point_index)>journal.deadline-time.monotonic()-75:
            if kappa!=fixed:
                stopped='FIXED_MOMENTUM_FALLBACK_ONLY_WITHIN_BUDGET';break
        begin=time.monotonic();v0=sp.cancel(kappa/sp.sqrt(5));entry={'omega':str(w0),'kappa':str(kappa),'ends':{}}
        if kappa==0:
            entry['status']='ZERO_TANGENTIAL_CHART_UNRESOLVED';summaries.append(entry);continue
        if w0==params['omega'] and kappa==fixed:
            entry['status']='REUSED_SAVED_SEED_METADATA_NO_MODE_REPLAY';entry['seed']=spec['savedSeedSummary'];summaries.append(entry);continue
        for label,b in ends.items():
            per={}
            for sector,(det,elim) in spectral[label].items():
                name='point-%03d-%s-%s'%(point_index,label,sector)
                census=op(name+'-roots',root_census,elim,b,w0,v0)
                results=[]
                if census['status']=='COMPLETE_NECESSARY_REAL_ROOT_CENSUS':
                    require(census['count']==sum(n for _,n in census['roots']) and zeros([census['realRemainder'],census['imagRemainder']]),'real root coverage')
                    for i,interval in enumerate(census['roots']):
                        result=op(name+'-root-%02d'%i,inspect_root,b,det,census,interval,sector);results.append(result)
                        require(result['status'] not in ('FAILED_ORIGINAL_MATRIX_RESIDUAL','ZERO_DRIVE_CURRENT_JOIN_FAILED','SOURCE_FORM_CONTROL_UNRESPONSIVE'), 'scientific control failed; no automatic retry')
                per[sector]={'census':census['status'],'candidateCount':len(census['roots']),
                             'statuses':[r['status'] for r in results],
                             'positiveCurrentRank':sum(r.get('positiveCurrentRank',0) for r in results),
                             'negativeCurrentRank':sum(r.get('negativeCurrentRank',0) for r in results)}
            per['bulkMaximum']=str(b['bulkMaximum'].subs({w:w0,v:v0})) if b['bulkMaximum'] is not None else None
            entry['ends'][label]=per
        entry['incidentTransverseAvailable']=(entry['ends']['LEFT']['T']['positiveCurrentRank']>0 or entry['ends']['RIGHT']['T']['negativeCurrentRank']>0)
        entry['status']='SAMPLED_CLASSIFICATIONS_NOT_CELL_COVERAGE_PROOF';summaries.append(entry);timings.append(time.monotonic()-begin)
        save(out/('point-summary-%03d.json'%point_index),entry)
    save(out/'grid-summary.json',summaries)
    return {'status':'MAP_FACE_PILOT_SAVED_FOR_INSPECTION','case':spec['case'],
            'physicalLightRegime':'NOT_ESTABLISHED_IN_REFERENCE_UNITS; light names the model transverse sector only',
            'azimuth':'saved 2:1 ray only','thresholdStatus':'Exact bulk kinematic boundary where computed; end loci are candidate equations pending physical/coverage interpretation',
            'gridPointsRecorded':len(summaries),'gridSeconds':time.monotonic()-grid_start,'gridStop':stopped,
            'completePlanePartition':False,'referenceFaceDrive':show(faces['REFERENCE']),
            'currentUnits':show(unit_record),'noBorn':True,'noBulkPower':True,'noNewAnchorSolve':True,
            'next':'Inspect saved actual results and, if needed, one separately guarded saved-output validator; return to user before Born.'}


def main():
    p=argparse.ArgumentParser(__doc__)
    for n in ('input-manifest','gate-receipt','run-directory'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();spec=json.loads(a.input_manifest.read_text());gate=json.loads(a.gate_receipt.read_text())
    require(gate['status']=='READY_FOR_ONE_GUARDED_CHANNEL_MAP_FACE_JOB','independent build gate required')
    require(gate['workerSha256']==sha(__file__) and gate['inputManifestSha256']==sha(a.input_manifest),'worker/input pins')
    require(gate['independentBuildClearance'] is True and gate['scienceJobOrdinal']==1,'review and one-constructor scope')
    require(gate['seconds']==900 and gate['nativeSeconds']==840 and gate['automaticRetry'] is False,'ordinary bounded run')
    require(gate['guardSha256']==sha(ROOT/'scripts/s11c_guarded_run.py') and gate['supervisorSha256']==sha(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'),'guard/supervisor pins')
    for rec in spec['inputs'].values():require(route(rec['path'])==rec,'source hash before import')
    require(json.loads(Path(spec['inputs']['physicalInputFile']['path']).read_text())==spec['physicalInput'],'unchanged material/profile specification')
    for join in spec['checkpointJoins']:
        checkpoint=json.loads(Path(join['checkpoint']).read_text())
        require(checkpoint['artifacts'][join['artifact']]['sha256']==spec['inputs'][join['inputKey']]['sha256'],'accepted artifact checkpoint join')
    limits=containment();start=time.monotonic();out=a.run_directory.resolve();out.relative_to(ROOT/'_scratch/s11c');out.mkdir(parents=True,exist_ok=False)
    j=Journal(out,start+840);save(out/'native-limits.json',limits);save(out/'manifest.json',spec);save(out/'gate.json',gate)
    try:
        result=construct(spec,out,j)
        post={n:route(r['path']) for n,r in spec['inputs'].items()};save(out/'posthashes.json',post);require(post==spec['inputs'],'source posthashes')
        result.update(wallSeconds=time.monotonic()-start,operations=len(j.records),allSourcesUnchanged=True)
        save(out/'operation-index.json',j.records);save(out/'checks.json',result);print(json.dumps(result,indent=2,allow_nan=False))
    except BaseException:
        signal.setitimer(signal.ITIMER_REAL,0)
        save(out/'failure.json',dict(traceback=traceback.format_exc(),incompleteOperation=j.active,wallSeconds=time.monotonic()-start,noRetry=True))
        save(out/'operation-index.json',j.records);save(out/'failure-posthashes.json',{n:route(r['path']) for n,r in spec['inputs'].items()});raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        save(out/'artifact-index.json',{str(p.relative_to(out)):route(p) for p in sorted(out.rglob('*')) if p.is_file()})


if __name__=='__main__':main()
