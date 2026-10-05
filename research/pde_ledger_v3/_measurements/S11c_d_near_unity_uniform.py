#!/usr/bin/env python3
"""Selected constant-end near-unity check from pinned original operands.

Source-only preparation is inert. Scientific imports/restoration occur only
inside the fresh method gate and pooled memory containment. No deadlines,
producer replay, numerical-mode reuse, full inverse or profile integration.
"""
import argparse
import ast
import builtins
from collections import OrderedDict
import hashlib
import importlib
import io
import json
import math
import os
from pathlib import Path
import pickle
import resource
import sys
import time
import traceback
from types import SimpleNamespace
from S11c_d_numerical_radiating_blob_store import BlobStore

ROOT = Path('/var/projects/toy_physics')
PLAN_SHA = '0d1eb443689f07558526bdf06f64d02072778376b0754096e60e6337ec2336a8'
INPUT_SHA = '660b2b0b54494dce964c7f5fbd722a76a7e894d8ed92b18b5eb6290c3c10f1d5'
PACKETS = ('uniformSource','uniformResponse','uniformCommon','LEFTFrequency','RIGHTFrequency',
           'LEFTPairing','RIGHTPairing','LEFTNativeSource','RIGHTNativeSource')
THREADS = ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
PORTS = ('AMPLITUDE','PRESSURE','OUTWARD_VELOCITY','RELATIVE_MASS_FLUX','AFFINITY',
         'CHEMICAL_POTENTIAL','BULK_VELOCITY','MECHANICAL_RESPONSE')
FORMS = ('SLAB_CURRENT_MATRIX','BULK_NORMAL_CURRENT_DENSITY_MATRIX','BULK_DEPTH_CURRENT_MATRIX',
         'INTERFACE_POWER_MATRIX')


class IntegrityError(RuntimeError):
    pass


class ScientificIssue(Exception):
    def __init__(self, reason, evidence=None, status='UNRESOLVED'):
        super().__init__(reason)
        self.evidence, self.status = evidence, status


def require(ok, reason):
    if not ok:
        raise IntegrityError(reason)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''): h.update(block)
    return h.hexdigest()


def verify(pin):
    path=Path(pin['path'])
    require(path.is_file() and digest(path)==pin['sha256'],'changed pinned file: '+str(path))
    if 'bytes' in pin: require(path.stat().st_size==pin['bytes'],'changed byte count: '+str(path))


def save(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    raw=(json.dumps(value,indent=2,allow_nan=False)+'\n').encode()
    with path.open('xb') as stream:
        stream.write(raw);stream.flush();os.fsync(stream.fileno())


def readable(value):
    if isinstance(value,dict): return {str(key):readable(item) for key,item in value.items()}
    if isinstance(value,(tuple,list,set,frozenset)): return [readable(item) for item in value]
    if isinstance(value,np.ndarray): return readable(value.tolist())
    if isinstance(value,(complex,np.complexfloating)):
        return dict(real=readable(float(value.real)),imag=readable(float(value.imag)))
    if isinstance(value,(float,np.floating)):
        item=float(value);return item if math.isfinite(item) else str(value)
    if isinstance(value,(np.integer,np.bool_)): return value.item()
    if isinstance(value,sp.MatrixBase): return readable(value.tolist())
    if isinstance(value,sp.Basic): return str(value)
    if value is None or isinstance(value,(str,int,bool)): return value
    return dict(type=type(value).__name__,representation=str(value))


class SavedCodec(pickle.Unpickler):
    def find_class(self,module,name):
        safe={('collections','OrderedDict'):OrderedDict,('builtins','complex'):builtins.complex,
              ('builtins','set'):builtins.set,('builtins','frozenset'):builtins.frozenset,
              ('builtins','slice'):builtins.slice,('numpy','ndarray'):np.ndarray,('numpy','dtype'):np.dtype,
              ('numpy.core.multiarray','_reconstruct'):np.core.multiarray._reconstruct,
              ('numpy.core.multiarray','scalar'):np.core.multiarray.scalar,
              ('numpy.core.numeric','_frombuffer'):np.core.numeric._frombuffer}
        if (module,name) in safe: return safe[module,name]
        if module.startswith('sympy.') and not name.startswith('__'):
            item=getattr(importlib.import_module(module),name)
            if isinstance(item,type) and issubclass(item,(sp.Basic,sp.MatrixBase)): return item
            if (module,name)==('sympy.core.function','_rebuild_undef'): return item
        raise pickle.UnpicklingError('unsupported saved class '+module+'.'+name)

    def persistent_load(self,identifier):
        raise pickle.UnpicklingError('persistent ID not supported')


def decode(raw):
    stream=io.BytesIO(raw);value=SavedCodec(stream).load()
    require(not stream.read(1),'trailing bytes in pinned source object')
    return value


class Journal:
    """Immutable full arguments/returns plus event log; no computation timer."""
    def __init__(self,out):
        self.out=out;self.store=BlobStore(out/'operations.sqlite',create=True)
        self.count=0;self.active=[];self.refs={};self.failures=[]

    def event(self,record):
        with (self.out/'operation-index.jsonl').open('a') as stream:
            stream.write(json.dumps(record,allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno())

    def emit(self,name,value):
        ref=self.store.put(name,pickle.dumps(value,protocol=5));self.refs[name]=ref
        self.event(dict(name=name,status='COMPLETE_EVIDENCE',value=ref))
        return ref

    def check(self,name,ok,evidence,status='FAILED_CHECK'):
        ref=self.emit(name,evidence)
        if not bool(ok): raise ScientificIssue(name,dict(artifact=ref),status)

    def call(self,name,args,fn,soft=False):
        inp=self.emit(name+'/input.pickle',args);self.active.append(dict(name=name,input=inp))
        self.event(dict(name=name,status='STARTED',input=inp));started=time.monotonic()
        try:
            value=fn()
        except ScientificIssue as error:
            value=dict(status=error.status,reason=str(error),evidence=error.evidence,
                       traceback=traceback.format_exc(),operation=name)
            failed=self.emit(name+'/failure.pickle',value);self.failures.append(value)
            self.event(dict(name=name,status=error.status,input=inp,failure=failed,seconds=time.monotonic()-started))
            self.active.pop()
            if not soft: raise
            return value
        except BaseException as error:
            self.emit(name+'/exception.pickle',dict(type=type(error).__name__,message=str(error),
                 traceback=traceback.format_exc(),active=list(self.active)))
            raise
        result=self.emit(name+'/return.pickle',value);self.count+=1
        self.event(dict(name=name,status='COMPLETE',input=inp,result=result,seconds=time.monotonic()-started))
        self.active.pop();return value

    def report(self,name,value):
        save(self.out/(name+'.json'),readable(value))


def leaves(value):
    if isinstance(value,sp.MatrixBase): return list(value)
    if isinstance(value,dict): return [item for part in value.values() for item in leaves(part)]
    if isinstance(value,(tuple,list)): return [item for part in value for item in leaves(part)]
    return [value]


def clean(value):
    if isinstance(value,sp.MatrixBase): return sp.ImmutableMatrix(value).applyfunc(sp.cancel)
    return sp.cancel(value)


def zeros(value):
    return all(sp.cancel(item)==0 for item in leaves(value))


def norm(value):
    return float(np.max(np.abs(value),initial=0.))


def numeric(value):
    a=np.asarray(sp.Matrix(value).evalf(40).tolist(),dtype=complex)
    if not np.isfinite(a).all(): raise ScientificIssue('nonfinite numerical coefficient',dict(source=value,numerical=a))
    return a


def factors(value):
    result=[]
    for item in leaves(value):
        if not isinstance(item,sp.Basic): continue
        result.extend(power.base for power in item.atoms(sp.Pow) if power.exp.is_negative is True)
        result.append(sp.fraction(sp.together(item))[1])
    return tuple(dict.fromkeys(result))


def carrier_degrees(value,epsilon,cache):
    if value in cache:return cache[value]
    if value==0: result=set()
    elif not value.has(epsilon): result={0}
    elif value==epsilon: result={1}
    elif value.is_Add: result=set().union(*(carrier_degrees(child,epsilon,cache) for child in value.args))
    elif value.is_Mul:
        result={0}
        for child in value.args: result={a+b for a in result for b in carrier_degrees(child,epsilon,cache)}
    elif value.is_Pow and value.exp.is_Integer and value.exp>=0:
        result={0}
        for _ in range(int(value.exp)): result={a+b for a in result for b in carrier_degrees(value.base,epsilon,cache)}
    else: raise ScientificIssue('unsupported epsilon dependence',value)
    cache[value]=result;return result


def source_method(tree,class_name,name):
    cls=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name==class_name)
    return next(node for node in cls.body if isinstance(node,ast.FunctionDef) and node.name==name)


def selected_lift(uniform,common,params,k,engine_tree,J):
    """Evaluate only the pinned source trial ansatz, never a producer."""
    method=source_method(engine_tree,'ReducedPencil','trial_ansatz')
    require(not method.decorator_list,'unadorned source ansatz method')
    namespace={'sp':sp};exec(compile(ast.Module(body=[method],type_ignores=[]),'<pinned-trial-ansatz>','exec'),namespace)
    z,x1,x2=sp.symbols('uniformTrialNormal uniformTrialTangent1 uniformTrialTangent2',real=True)
    t1,t2=(params['s11cdTangentialMomentum'+str(i)] for i in (1,2))
    phase=sp.exp(sp.I*(t1*x1+t2*x2+k*z))
    field_names=('u1','u2','u3','theta','eW')
    # The source ansatz uses these exact physical field labels; saved unit
    # registries below supply the actual dimension correspondence.
    fields=tuple(sp.Function('s11cdReducedField'+name) for name in field_names)
    context=SimpleNamespace(fields=fields,tangent_derivatives=tuple(sp.cancel(sp.diff(phase,x)/phase) for x in (x1,x2)))
    labels=('A0','A1','A2','Theta','E','Phi');amplitudes=sp.symbols('uniformLiftAmplitude0:6')
    trials={sp.Function('s11cdReducedTrial'+label)(z):amplitudes[i]*sp.exp(sp.I*k*z) for i,label in enumerate(labels)}
    ansatz=[namespace['trial_ansatz'](context,sector) for sector in ('TRANSVERSE','THETA','E_W','LONGITUDINAL')]
    field_ansatz=sp.ImmutableMatrix([sum(item[field](z) for item in ansatz) for field in fields])
    evaluated=field_ansatz.subs(trials,simultaneous=True).doit()
    full=sp.ImmutableMatrix(5,6,lambda i,j:sp.cancel(sp.diff(evaluated[i],amplitudes[j])/sp.exp(sp.I*k*z)))
    def curl_bind(value):
        mapping={symbol:(k if symbol.name=='s11cdSpectralNormalMomentum' else params[symbol.name])
                 for symbol in value.free_symbols}
        return value.xreplace(mapping)
    curl=curl_bind(uniform['curl']);other=curl_bind(common['curl'])
    lift=sp.ImmutableMatrix.vstack(curl[:,1:3],sp.zeros(2,2))
    gauge=tuple(curl_bind(item) for item in common['gauge'])
    residual=dict(fieldLift=full[:,1:3]-lift,curl=other-curl,gauge=tuple(curl*item for item in gauge))
    gram=clean(lift.H*lift)
    evidence=dict(sourceMethod=ast.unparse(method),sourceFieldLift=ast.unparse(source_method(engine_tree,'ConstantEndPencil','field_lift')),
                  trialAnsatz=field_ansatz,evaluatedAnsatz=evaluated,fieldLift=full,selectedLift=lift,
                  curl=curl,gauge=gauge,gram=gram,residuals=residual,fieldOrder=field_names,
                  interpretation='FRESH_SOURCE_ANSATZ_COEFFICIENTS_NOT_RESTORED_FINGERPRINT')
    J.check('lift/source-joins',zeros(residual) and gram.det()!=0,evidence)
    return dict(lift=lift,full=full,curl=curl,gram=gram,gauge=gauge,fields=fields,source=evidence)


class EndBinding:
    def __init__(self,label,restored,params,w,cs,k,km,q,qm,Q,J):
        self.label,self.J=label,J;self.params=params;self.fixed_frequency=None;self.w,self.cs,self.k,self.km,self.q,self.qm,self.Q=w,cs,k,km,q,qm,Q
        self.uniform=restored['uniformSource'];self.frequency=restored[label+'Frequency']
        self.packet,self.dimensions=restored[label+'Pairing'];self.r=self.packet['result'];self.checks=self.packet['checks']
        self.native=restored[label+'NativeSource'];self.a=self.native['acoustic'];self.slab=self.native['slab']
        self.endpoints=self.uniform['profileEndpoints'];self.origin=self.frequency['origin'];self.binding_records=[]
        require(set(symbol.name for symbol in self.origin)=={'eta_bg','sigma_W'},'two original physical grades')
        self.wl,self.wr=self.r['FREQUENCY_LEGS'];self.kl,self.kr=self.r['NORMAL_LEGS'];self.ql,self.qr=self.r['BULK_LEGS']
        self.base={self.wl:w,self.wr:w,self.kl:km,self.kr:k,self.ql:qm,self.qr:q,
                   self.frequency['momentum']:k,self.frequency['radical']:Q,self.frequency['frequency']:w}
        self.amplitudes=tuple(tuple(sorted((symbol for symbol in self.dimensions
             if getattr(symbol,'name','').startswith('s11cdCurrent'+side+'Amplitude')),key=str)) for side in ('Plus','Minus'))
        require(all(len(group)==5 for group in self.amplitudes),'five actual field amplitudes per harmonic leg')
        self.epsilon=next(symbol for symbol in self.r['SLAB_CURRENT_MATRIX'].free_symbols if symbol.name=='epsilon_shape')

    def bind(self,value,carry=(),extra=None,grades=True,label=None):
        if isinstance(value,dict):return {key:self.bind(item,carry,extra,grades,label) for key,item in value.items()}
        if isinstance(value,(list,tuple)):return tuple(self.bind(item,carry,extra,grades,label) for item in value)
        if not isinstance(value,(sp.Basic,sp.MatrixBase)):return value
        source=value;value=value.xreplace(self.endpoints)
        mapping=dict(self.base);mapping.update(extra or {})
        for symbol in value.free_symbols-set(mapping):
            name=symbol.name
            if symbol in carry or (not grades and symbol in self.origin):continue
            if symbol in self.origin:mapping[symbol]=self.origin[symbol]
            elif name in ('omega','s11cdFrequency'):mapping[symbol]=self.w
            elif name=='c_s0':mapping[symbol]=self.cs
            elif name in self.params:mapping[symbol]=self.params[name]
            elif symbol in (self.w,self.cs,self.k,self.km,self.q,self.qm,self.Q):continue
            else:raise ScientificIssue('unexpected unbound source input',dict(end=self.label,symbol=symbol,source=source))
        result=value.xreplace(mapping)
        if self.fixed_frequency is not None:result=result.xreplace({self.w:self.fixed_frequency})
        if result.atoms(sp.Integral,sp.Limit,sp.Derivative,sp.Subs):
            raise ScientificIssue('unresolved source operator after fixed binding',dict(source=source,mapping=mapping,result=result))
        self.binding_records.append(dict(label=label,sourceSymbols=tuple(sorted(source.free_symbols,key=str)),mapping=mapping,
                                         remainingSymbols=tuple(sorted(result.free_symbols,key=str)),fixedFrequency=self.fixed_frequency))
        return result

    def native_leg(self,value,leg,carry=(),grades=True):
        sign=1 if leg==0 else -1
        normal=self.k if leg==0 else self.km;depth=self.q if leg==0 else self.qm
        actual={}
        for item in leaves(value):
            if not isinstance(item,(sp.Basic,sp.MatrixBase)):continue
            for symbol in item.free_symbols:
                if symbol.name in ('omega','s11cdFrequency'):actual[symbol]=sign*self.w
                elif symbol==self.frequency['momentum']:actual[symbol]=sign*normal
                elif symbol==self.qr:actual[symbol]=sign*depth
                elif symbol.name in ('s11cdTangentialMomentum1','s11cdTangentialMomentum2'):
                    actual[symbol]=sign*self.params[symbol.name]
        actual.update(dict(zip(self.amplitudes[0],self.amplitudes[leg])))
        return self.bind(value,carry=carry,extra=actual,grades=grades)

    def retained(self,value):
        eta,sigma=(next(symbol for symbol in self.origin if symbol.name==name) for name in ('eta_bg','sigma_W'))
        def scalar(item):
            return sp.expand(sum(sp.diff(item,eta,a,sigma,b).subs({eta:0,sigma:0})*eta**a*sigma**b
                                 for a in range(2) for b in range(2)))
        return value.applyfunc(scalar) if isinstance(value,sp.MatrixBase) else scalar(value)


def on_wave(value,wave,q):
    def scalar(item):
        item=sp.cancel(item)
        if item==0:return item
        numerator,denominator=sp.fraction(item)
        remainder=sp.rem(sp.Poly(numerator,q,domain='EX'),sp.Poly(wave,q,domain='EX')).as_expr()
        return sp.cancel(remainder/denominator)
    return value.applyfunc(scalar) if isinstance(value,sp.MatrixBase) else scalar(value)


def end_restriction(B,lift,J):
    label=B.label;prefix=label+'/restriction';k,w,q,Q,cs=B.k,B.w,B.q,B.Q,B.cs
    P=B.bind(B.r['CLOSED_PENCIL_LEGS'][0],label='physical right five-field P')
    Pminus=B.bind(B.r['CLOSED_PENCIL_LEGS'][1],label='physical left five-field P')
    wave=B.bind(B.r['ACOUSTIC_WAVE_ROWS'][0]);wave_minus=B.bind(B.r['ACOUSTIC_WAVE_ROWS'][1])
    algebraic=B.bind(B.frequency['originalAlgebraic']);relation=B.bind(B.frequency['originalRelation'])
    original_omega=tuple(s for s in B.frequency['originalAlgebraic'].free_symbols if s.name=='omega')
    require(len(original_omega)==1,'actual original frequency symbol correspondence for signed source leg')
    original_minus=B.bind(B.frequency['originalAlgebraic'],extra={B.frequency['momentum']:-B.km,
         B.frequency['radical']:-Q,**{s:-w for s in B.frequency['originalAlgebraic'].free_symbols if s.name=='omega'},
         **{s:-B.params[s.name] for s in B.frequency['originalAlgebraic'].free_symbols if s.name in ('s11cdTangentialMomentum1','s11cdTangentialMomentum2')}})
    wp,rp=sp.Poly(wave,q),sp.Poly(relation,Q)
    J.check(prefix+'/wave-schema',wp.degree()==rp.degree()==2 and wp.nth(1)==rp.nth(1)==0,
            dict(wave=wave,relation=relation,physicalPolynomial=wp.as_expr(),algebraicPolynomial=rp.as_expr()))
    conversion=sp.sqrt(sp.cancel(wp.nth(2)/rp.nth(2)))
    # Q = conversion * q; the source's ACOUSTIC_RADICAL_SCALE is the inverse.
    native_inverse=B.bind(B.a['ACOUSTIC_RADICAL_SCALE'])
    wave_join=clean(relation.subs(Q,conversion*q)-wave)
    plus_join=clean(algebraic.subs(Q,conversion*q)-P)
    minus_join=clean(original_minus.subs(Q,conversion*B.qm)-Pminus)
    wave_minus_join=clean(wave_minus-wave.xreplace({k:B.km,q:B.qm}))
    units=B.dimensions|B.native['knownDimensions']
    dimensional=dict(Q=units.get(B.frequency['radical']),q=units.get(B.qr),qMinus=units.get(B.ql),
                     frequency=units.get(B.wr),normal=units.get(B.kr),fieldUnits=lift['fieldUnits'],
                     rowEntryUnits=B.uniform['units']['strong'])
    J.check(prefix+'/physical-algebraic-joins',zeros((wave_join,plus_join,minus_join,wave_minus_join,conversion*native_inverse-1)),
        dict(P=P,Pminus=Pminus,algebraic=algebraic,algebraicMinus=original_minus,relation=relation,wave=wave,
             waveMinus=wave_minus,conversion=conversion,nativeInverseScale=native_inverse,units=dimensional,
             residuals=(wave_join,plus_join,minus_join,wave_minus_join,conversion*native_inverse-1)))
    J.check(prefix+'/units',tuple(dimensional['Q'] or ())==(0,-1,0) and tuple(dimensional['q'] or ())==(-1,0,0),dimensional,'UNRESOLVED')
    L=lift['lift'];gram=lift['gram'];D=clean(gram.inv()*L.H*P*L);invariant=clean(P*L-L*D)
    reduced=on_wave(invariant,wave,q)
    J.check(prefix+'/full-five-row-invariance',zeros(reduced),dict(P=P,lift=L,gram=gram,D=D,
        rawResidual=invariant,onWaveResidual=reduced,denominatorDomain=factors((P,L,D,invariant)),TH=B.bind(B.uniform['records'][label]['coupling']['TH']),
        HT=B.bind(B.uniform['records'][label]['coupling']['HT'])))
    if any(item.has(q,B.qm,Q) for item in D):
        raise ScientificIssue('selected restriction retains acoustic radical',dict(D=D,wave=wave))
    scalar=D[0,0];degeneracy=clean(D-sp.eye(2)*scalar)
    J.check(prefix+'/doublet-dispersion',zeros(degeneracy),dict(D=D,scalar=scalar,residual=degeneracy),'UNRESOLVED')
    polynomial=sp.Poly(sp.fraction(sp.cancel(scalar))[0],k)
    J.check(prefix+'/normal-dispersion-degree',polynomial.degree()==2 and polynomial.nth(1)==0,
            dict(dispersion=scalar,polynomial=polynomial.as_expr()),'UNRESOLVED')
    k2=sp.cancel(-polynomial.nth(0)/polynomial.nth(2))
    kt2=sum(B.params['s11cdTangentialMomentum'+str(i)]**2 for i in (1,2))
    cT=sp.simplify(w/sp.sqrt(sp.cancel(k2+kt2)))
    if cT.has(cs):
        fallback=sp.solve((sp.fraction(scalar.subs({q:0,w:3}))[0],wave.subs({q:0,w:3})),(k,cs),dict=True)
        J.emit(prefix+'/small-match-fallback',dict(smallSystem=(scalar,wave),solutions=fallback))
        raise ScientificIssue('sound-speed-independent selected branch unavailable',dict(cT=cT,fallback=fallback))
    cT3=sp.simplify(cT.subs(w,3));k23=sp.simplify(k2.subs(w,3))
    direction=sp.cancel(-sp.diff(scalar,k)/sp.diff(scalar,w))
    onroot=sp.simplify(scalar.subs(k,sp.sqrt(k2)))
    J.check(prefix+'/selected-branch',cT3.is_positive is True and k23.is_positive is True and onroot==0,
        dict(dispersion=scalar,normalSquared=k2,tangentSquared=kt2,cTDefinition=w/sp.sqrt(k2+kt2),
             cT=cT,cTAt3=cT3,normalSquaredAt3=k23,direction=direction,dispersionRootResidual=onroot),'UNRESOLVED')
    return dict(P=P,Pminus=Pminus,wave=wave,waveMinus=wave_minus,algebraic=algebraic,algebraicMinus=original_minus,
                relation=relation,conversion=conversion,nativeInverse=native_inverse,D=D,dispersion=scalar,
                k2=k2,k23=k23,kt2=kt2,cT=cT3,direction=direction,lift=L,gram=gram,
                invariant=invariant,onWaveInvariant=reduced,reductionDomains=factors((P,L,D,invariant)),units=dimensional)


def native_expression(node,env):
    """Scalar closure AST only; no producer call or replacement constitutive law."""
    if isinstance(node,ast.Constant):return node.value
    if isinstance(node,ast.Name):return env[node.id]
    if isinstance(node,ast.Subscript):return native_expression(node.value,env)[native_expression(node.slice,env)]
    if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -native_expression(node.operand,env)
    if isinstance(node,ast.BinOp):
        a,b=native_expression(node.left,env),native_expression(node.right,env)
        if isinstance(node.op,ast.Add):return a+b
        if isinstance(node.op,ast.Sub):return a-b
        if isinstance(node.op,ast.Mult):return a*b
        if isinstance(node.op,ast.Div):return a/b
    raise IntegrityError('unsupported pinned scalar closure AST: '+ast.dump(node))


def native_reconstruction(B,restriction,engine_tree,J):
    prefix=B.label+'/native';r,c,a=B.r,B.checks,B.a;e=B.epsilon;amps=B.amplitudes
    allamps=(*amps[0],*amps[1]);rho=B.params['rho_m'];forms={};homogeneity={}
    for name in FORMS:
        original=B.bind(r[name],carry=(e,),label=name)
        degrees=[sorted(carrier_degrees(item,e,{})) for item in original]
        coefficient=original.subs(e,1)
        scaling={str(scale):clean(original.subs(e,scale)-scale**2*coefficient)
                 for scale in (sp.Rational(1,2),sp.Integer(2))}
        J.check(prefix+'/'+name+'/carrier',all(degree in ([],[2]) for degree in degrees) and zeros(scaling),
            dict(original=original,coefficient=coefficient,degrees=degrees,scalingResiduals=scaling))
        forms[name]=coefficient;homogeneity[name]=dict(degrees=degrees,rescaling=scaling)
    native_method=source_method(engine_tree,'ClosedAcousticEnergy','construct')
    face_loop=next(node for node in native_method.body if isinstance(node,ast.For) and isinstance(node.target,ast.Name) and node.target.id=='sign')
    closure_nodes={name:next(node.value for node in face_loop.body if isinstance(node,ast.Assign)
                    and any(isinstance(target,ast.Name) and target.id==name for target in node.targets))
                   for name in ('relative','affinity','closure')}
    faces=[];closures=[];chemical=[];coefficient_residuals=[]
    require(len(r['FACE_LEG_OBJECTS'])==len(c['FACE_AMPLITUDE_ROWS'])==len(c['FACE_PORT_ROWS'])==len(a['FACE_RECORDS'])==2,
            'complete two-face native census')
    for leg in range(2):
        raw_row=B.native_leg(B.slab['CHEMICAL_FIELD_ROW'],leg,grades=False)
        density=B.native_leg(a['SLAB_SURFACE_DENSITY'],leg,grades=False)
        raw_driver=(raw_row*sp.ImmutableMatrix(amps[leg]))[0]/density
        retained=B.retained(raw_driver);driver=retained.xreplace(B.origin)
        saved=B.native_leg(a['CHEMICAL_AFFINITY_DRIVER'],leg,carry=allamps)
        residual=clean(driver-saved)
        J.check(prefix+'/chemical-leg-'+str(leg),residual==0,
            dict(rawChemicalRow=raw_row,surfaceDensity=density,rawQuotient=raw_driver,
                 retainedDriver=retained,physicalDriver=driver,savedDriver=saved,residual=residual,
                 projection='NATIVE_RECTANGULAR_ETA_0_1_SIGMA_0_1_BEFORE_PHYSICAL_ORIGIN'))
        chemical.append(driver)
    pressure_scales=B.bind(r['OPEN_BULK_PRESSURE_SCALES'])
    velocity_scales=B.bind(r['OPEN_BULK_VELOCITY_SCALES'])
    memory=[B.native_leg(a['MEMORY_KERNELS'],leg) for leg in range(2)]
    for face in range(2):
        require(len(r['FACE_LEG_OBJECTS'][face])==len(c['FACE_AMPLITUDE_ROWS'][face])==len(c['FACE_PORT_ROWS'][face])==2,
                'both native harmonic legs')
        legrows=[];face_closures=[]
        for leg in range(2):
            rawleg=r['FACE_LEG_OBJECTS'][face][leg]
            expressions=B.bind(rawleg,carry=allamps)
            rows=B.bind({'AMPLITUDE':c['FACE_AMPLITUDE_ROWS'][face][leg],**c['FACE_PORT_ROWS'][face][leg]})
            require(set(rows)==set(PORTS) and all(row.shape==(1,5) for row in rows.values()),'complete native five-column port rows')
            reconstruction={name:clean(expressions[name]-(rows[name]*sp.ImmutableMatrix(amps[leg]))[0]) for name in PORTS}
            identities=dict(pressure=expressions['PRESSURE']-pressure_scales[leg]*expressions['AMPLITUDE'],
                bulkVelocity=expressions['BULK_VELOCITY']-velocity_scales[leg]*expressions['AMPLITUDE'],
                affinity=expressions['AFFINITY']-expressions['CHEMICAL_POTENTIAL']+expressions['PRESSURE']/rho,
                relativeMass=expressions['RELATIVE_MASS_FLUX']-rho*(expressions['BULK_VELOCITY']-expressions['OUTWARD_VELOCITY']),
                memory=expressions['RELATIVE_MASS_FLUX']-expressions['MASS_RESPONSE_A']*expressions['AFFINITY']-
                       expressions['MASS_RESPONSE_V']*expressions['OUTWARD_VELOCITY'],
                chemical=expressions['CHEMICAL_POTENTIAL']-chemical[leg],
                responseA=expressions['MASS_RESPONSE_A']-memory[leg]['A'],
                responseV=expressions['MASS_RESPONSE_V']-memory[leg]['V'],
                mechanical=expressions['MECHANICAL_RESPONSE']-memory[leg]['X']*expressions['AFFINITY'])
            identities={name:clean(value) for name,value in identities.items()}
            A,p,v=sp.symbols('closureAmplitude closurePressure closureBulkVelocity')
            env=dict(rho=rho,face_v=velocity_scales[leg]*A,face_p=pressure_scales[leg]*A,
                     outward_velocity=expressions['OUTWARD_VELOCITY'],mus=chemical[leg],kernels=memory[leg])
            for name in ('relative','affinity','closure'):env[name]=native_expression(closure_nodes[name],env)
            coefficient=sp.diff(env['closure'],A);constant=env['closure'].subs(A,0)
            native_coefficient=B.native_leg(a['FACE_RECORDS'][face]['AMPLITUDE_EQUATION_COEFFICIENT'],leg)
            closure_reconstruction=clean(env['closure']-coefficient*A-constant)
            native_coefficient_join=clean(coefficient-native_coefficient)
            closure_residual=clean(env['closure'].subs(A,expressions['AMPLITUDE']))
            env2=dict(env,face_v=v,face_p=p)
            for name in ('relative','affinity','closure'):env2[name]=native_expression(closure_nodes[name],env2)
            acoustic_eliminant=sp.resultant(pressure_scales[leg]*A-p,velocity_scales[leg]*A-v,A)
            closure_matrix=sp.ImmutableMatrix([acoustic_eliminant,env2['closure']]).jacobian((p,v))
            closure_det=clean(closure_matrix.det())
            raw=dict(face=face,orientation=a['FACE_RECORDS'][face]['ORIENTATION'],leg=leg,rows=rows,
                expressions=expressions,amplitudes=amps[leg],rowReconstruction=reconstruction,identities=identities,
                coefficient=coefficient,nativeCoefficient=native_coefficient,coefficientResidual=native_coefficient_join,
                constant=constant,closure=env['closure'],closureReconstruction=closure_reconstruction,
                closureResidual=closure_residual,pressureVelocityMatrix=closure_matrix,
                pressureVelocityDeterminant=closure_det,pressureScale=pressure_scales[leg],velocityScale=velocity_scales[leg],
                impedance=sp.cancel(pressure_scales[leg]/velocity_scales[leg]),
                closureSource={name:ast.unparse(node) for name,node in closure_nodes.items()})
            J.check(prefix+'/face-%d-leg-%d/reconstruction'%(face,leg),
                zeros((reconstruction,identities,native_coefficient_join,closure_reconstruction,closure_residual)),raw)
            legrows.append(rows);face_closures.append(raw);coefficient_residuals.append(native_coefficient_join)
        faces.append(legrows);closures.append(face_closures)
    # Recompute the polarization coefficient from the original row-power
    # derivatives before the equal-frequency binding removes their variables.
    pairweight_raw=sp.cancel(sp.diff(r['PLUS_ROW_POWER_MAP'][0,0],B.wl)/sp.diff(r['BEAT_RATES'][0],B.wl))
    pairweight=B.bind(pairweight_raw,carry=(e,));savedweight=B.bind(c['PORT_BILINEAR_WEIGHT'],carry=(e,))
    J.check(prefix+'/polarization-weight',clean(pairweight-savedweight)==0,
            dict(rawComputed=pairweight_raw,computed=pairweight,saved=savedweight,residual=clean(pairweight-savedweight)))
    # Small harmonic source-work ansatz with actual chemical driver; source
    # assignment bodies are pinned and evaluated, not a new constitutive rule.
    current_method=source_method(engine_tree,'ClosedCurrentPairing','construct')
    source_assignments=[node for node in current_method.body if
        (isinstance(node,ast.Assign) and any(isinstance(target,ast.Name) and target.id=='source_ansatz' for target in node.targets)) or
        (isinstance(node,ast.AugAssign) and isinstance(node.target,ast.Name) and node.target.id=='source_ansatz')]
    require(len(source_assignments)==2,'native source-work ansatz statement schema')
    t=sp.Symbol('uniformPowerTime',real=True);h=sp.Symbol('uniformPowerHarmonic',nonzero=True)
    phase=h*sp.exp(-sp.I*B.fixed_frequency*t);residual_amps=(sp.symbols('uniformPlusRowResidual0:5'),sp.symbols('uniformMinusRowResidual0:5'))
    real_field=lambda plus,minus:e*(plus*phase+minus/phase)/2
    env={'sp':sp,'r':SimpleNamespace(t=t),'field_ansatz':tuple(real_field(x,y) for x,y in zip(*amps)),
         'residual_fields':tuple(real_field(x,y) for x,y in zip(*residual_amps)),
         'chemical_field':real_field(*chemical)}
    exec(compile(ast.Module(body=source_assignments,type_ignores=[]),'<pinned-source-work-ansatz>','exec'),env)
    average=sp.expand(env['source_ansatz']).coeff(h,0)
    plus=sp.ImmutableMatrix(5,5,lambda i,j:sp.diff(average,amps[1][i],residual_amps[0][j]))
    minus=sp.ImmutableMatrix(5,5,lambda i,j:sp.diff(average,residual_amps[1][i],amps[0][j]))
    savedplus=B.bind(r['PLUS_ROW_POWER_MAP'],carry=(e,));savedminus=B.bind(r['MINUS_ROW_POWER_MAP'],carry=(e,))
    power=B.bind(r['SOURCE_POWER_MATRIX'],carry=(e,))
    rowpower_residuals=dict(plus=clean(plus-savedplus),minus=clean(minus-savedminus),
        sourcePower=clean(savedplus*restriction['P'].subs(B.w,B.fixed_frequency)+restriction['Pminus'].subs(B.w,B.fixed_frequency).T*savedminus-power))
    J.check(prefix+'/row-power',zeros(rowpower_residuals),dict(sourceAssignments=[ast.unparse(node) for node in source_assignments],
        sourceAnsatz=env['source_ansatz'],harmonicAverage=average,plus=plus,minus=minus,savedPlus=savedplus,
        savedMinus=savedminus,savedPower=power,residuals=rowpower_residuals))
    open_names=('OPEN_BULK_ENERGY_COEFFICIENT','OPEN_BULK_CURRENT_COEFFICIENTS')
    open_values=(r[open_names[0]],r[open_names[1]][2],r[open_names[1]][3]);open_coefficients=[];open_residuals=[]
    for ordinal,value in enumerate(open_values):
        acoustic_amps=tuple(sorted((symbol for symbol in value.free_symbols if symbol.name in
            ('s11cdAcousticLeftAmplitude','s11cdAcousticRightAmplitude')),key=str))
        require(len(acoustic_amps)==2,'original acoustic amplitude pair')
        bound=B.bind(value,carry=(*acoustic_amps,e))
        coefficient=sp.diff(bound,*acoustic_amps)
        residual=clean(bound-coefficient*sp.prod(acoustic_amps))
        saved=B.bind(c['OPEN_BULK_PAIR_COEFFICIENTS'][ordinal],carry=(e,))
        open_coefficients.append(coefficient);open_residuals.append((residual,clean(coefficient-saved)))
    recon=[]
    for face,rows in enumerate(faces):
        column=lambda name:tuple(row[name] for row in rows)
        def bilinear(left,right):return pairweight*(left[1].T*right[0]+right[1].T*left[0])
        p,v,u,m,aff,mu,response=map(column,('PRESSURE','BULK_VELOCITY','OUTWARD_VELOCITY','RELATIVE_MASS_FLUX',
                                         'AFFINITY','CHEMICAL_POTENTIAL','MECHANICAL_RESPONSE'))
        constructed={'PORT_POWER_MATRIX':bilinear(p,u)+bilinear(response,u)+bilinear(mu,m),
                     'SPLIT_PORT_POWER_MATRIX':bilinear(p,v)+bilinear(aff,m)+bilinear(response,u),
                     'INTERFACE_POWER_MATRIX':bilinear(aff,m)+bilinear(response,u)}
        outer=rows[1]['AMPLITUDE'].T*rows[0]['AMPLITUDE']
        constructed.update({name:coefficient*outer for name,coefficient in zip(
             ('ENERGY_DENSITY_MATRIX','NORMAL_CURRENT_DENSITY_MATRIX','DEPTH_CURRENT_MATRIX'),open_coefficients)})
        actual={name:B.bind(r['FACE_RECORDS'][face][name],carry=(e,)) for name in constructed}
        residual={name:clean(actual[name]-value) for name,value in constructed.items()}
        J.check(prefix+'/face-%d/bilinears'%face,zeros(residual),dict(actual=actual,reconstructed=constructed,residuals=residual))
        recon.append(dict(actual=actual,reconstructed=constructed,residuals=residual))
    J.check(prefix+'/open-pair-coefficients',zeros(open_residuals),dict(coefficients=open_coefficients,residuals=open_residuals))
    raw_objects=dict(P=restriction['P'],Pminus=restriction['Pminus'],lift=restriction['lift'],faces=faces,forms=forms,
                     closureCoefficients=[[item['coefficient'] for item in face] for face in closures],
                     closureDeterminants=[[item['pressureVelocityDeterminant'] for item in face] for face in closures])
    denominator_inventory={name:factors(value) for name,value in raw_objects.items()}
    denominator_inventory['selectedRestriction']=restriction['reductionDomains']
    equal_waves=tuple(value.subs({B.w:B.fixed_frequency,B.qm:B.q}) for value in (restriction['wave'],restriction['waveMinus']))
    equal_eliminant=sp.rem(equal_waves[1],equal_waves[0],B.q)
    saved_equal=B.bind(c['EQUAL_DEPTH_WAVE_ROWS']);saved_eliminant=B.bind(c['EQUAL_DEPTH_WAVE_ELIMINANT'])
    diagonal=r['EQUAL_DEPTH_INTEGRAL'];height_symbols=tuple(diagonal.free_symbols)
    require(len(height_symbols)==1,'native diagonal depth has one physical height variable')
    height=height_symbols[0];depth_symbols=r['DEPTH_PHASE_FACTOR'].free_symbols-{B.ql,B.qr}
    require(len(depth_symbols)==1,'native depth-phase coordinate correspondence')
    depth=next(iter(depth_symbols))
    endpoint_integrand=r['DEPTH_PHASE_FACTOR'].subs({B.ql:B.qr,depth:height})
    depth_residuals=(sp.diff(diagonal,height)-endpoint_integrand,diagonal.subs(height,0))
    equal_depth=dict(waves=equal_waves,eliminant=equal_eliminant,savedWaves=saved_equal,savedEliminant=saved_eliminant,
        diagonalIntegral=diagonal,diagonalDomain=r['EQUAL_DEPTH_DOMAIN'],genericDomain=r['GENERIC_DEPTH_DOMAIN'],
        genericIntegral=r['GENERIC_DEPTH_INTEGRAL'],height=height,heightUnit=B.dimensions.get(height),
        endpointIntegrand=endpoint_integrand,derivativeAndEndpointResiduals=depth_residuals,
        usage='DIAGONAL_SOURCE_OPERAND_ONLY_NO_DEPTH_INTEGRATION_OR_TOTAL_CURRENT_CLAIM')
    J.check(prefix+'/equal-depth-admitted-source',zeros((tuple(x-y for x,y in zip(equal_waves,saved_equal)),
        equal_eliminant-saved_eliminant,depth_residuals)),equal_depth)
    J.emit(prefix+'/unprojected-domain-inventory',dict(objects=raw_objects,factors=denominator_inventory,
        dependence={name:tuple(sorted(set().union(*(item.free_symbols for item in leaves(value))),key=str))
                    for name,value in raw_objects.items()}))
    J.emit(prefix+'/actual-bindings',B.binding_records)
    return dict(faces=faces,forms=forms,closures=closures,chemical=chemical,memory=memory,
                pressureScales=pressure_scales,velocityScales=velocity_scales,pairWeight=pairweight,
                sourceRowPower=(savedplus,savedminus),homogeneity=homogeneity,
                openCoefficients=open_coefficients,bilinearReconstruction=recon,domains=denominator_inventory,
                rawObjects=raw_objects,coefficientResiduals=coefficient_residuals,equalDepth=equal_depth)


def transform(value,fn):
    if isinstance(value,dict):return {key:transform(item,fn) for key,item in value.items()}
    if isinstance(value,(tuple,list)):return tuple(transform(item,fn) for item in value)
    if isinstance(value,sp.MatrixBase):return value.applyfunc(fn)
    if isinstance(value,sp.Basic):return fn(value)
    return value


def finite(value):
    return not any(isinstance(item,sp.Basic) and item.has(sp.oo,-sp.oo,sp.zoo,sp.nan,sp.Limit)
                   for item in leaves(value))


def projected_objects(R,N):
    L=R['lift'];left=L.xreplace({next(s for s in L.free_symbols if s.name=='uniformNormal'): 
                              next(s for s in R['Pminus'].free_symbols if s.name=='uniformMinusNormal')})
    return dict(pencil=R['P']*L,minusPencil=R['Pminus']*sp.conjugate(left),
        forms={name:left.H*value*L for name,value in N['forms'].items()},
        faces=[[{name:row*(L if leg==0 else sp.conjugate(left)) for name,row in rows.items()}
                 for leg,rows in enumerate(face)] for face in N['faces']])


def physical_substitution(B,R,speed,sign):
    k=sign*sp.sqrt(R['k23']);basic={B.w:sp.Integer(3),B.cs:speed,B.k:k,B.km:k}
    polynomial=sp.Poly(R['wave'].subs(basic),B.q)
    q2=sp.simplify(-polynomial.nth(0)/polynomial.nth(2))
    if q2==0:depth=sp.Integer(0);stratum='EXACT_GRAZING'
    elif q2.is_positive is True:depth=sp.sqrt(q2);stratum='RADIATING'
    elif q2.is_negative is True:depth=sp.I*sp.sqrt(-q2);stratum='EVANESCENT'
    else:raise ScientificIssue('unresolved exact depth stratum',dict(speed=speed,k=k,depthSquared=q2))
    mapping={**basic,B.q:depth,B.qm:sp.conjugate(depth)}
    mapping[B.Q]=sp.simplify(R['conversion'].subs(basic)*depth)
    return mapping,dict(normal=k,depth=depth,depthSquared=q2,stratum=stratum)


def source_controls(B,N,J):
    """Remove the actual native eW input, then repeat the saved row route."""
    omission=[]
    for face in range(2):
        pair=[]
        for leg in range(2):
            original=B.r['FACE_LEG_OBJECTS'][face][leg]['OUTWARD_VELOCITY']
            omitted=original.xreplace({B.amplitudes[leg][4]:sp.Integer(0)})
            source_row=sp.ImmutableMatrix(1,5,lambda i,j:sp.diff(omitted,B.amplitudes[leg][j]))
            changed=B.bind(source_row);baseline=N['faces'][face][leg]['OUTWARD_VELOCITY']
            item=dict(source=original,changedSource=omitted,removedInput=B.amplitudes[leg][4],
                      rederivedRow=source_row,changed=changed,baseline=baseline)
            J.emit(B.label+'/controls/face-%d-leg-%d/ew-source-omission'%(face,leg),item)
            pair.append(item)
        omission.append(pair)
    return omission


def uniform_applicability(B,manifest,J):
    record=B.uniform['records'][B.label];background=record['background']
    operands=background['profileLimitOperands'];endpoints=B.endpoints
    source=B.bind(record['strong']);native=B.native['strong']
    if isinstance(native,sp.MatrixBase) and native.shape==source.shape:native=B.bind(native)
    # The original native strong object is an operator mapping in some source
    # versions. Its exact operand is retained; no unsupported mapping is guessed.
    source_join=clean(native-source) if isinstance(native,sp.MatrixBase) and native.shape==source.shape else None
    unresolved=record['unresolved'];actual_atoms=source.atoms(sp.Integral,sp.Derivative,sp.Limit)
    profile_evidence=[]
    for item in operands:
        profile_evidence.append(transform(item,lambda expr:expr.xreplace(endpoints)))
    actual=dict(profileLimitOperands=operands,specializedProfileOperands=profile_evidence,
        profileEndpoints=endpoints,nativeProfileBindings=B.native['profileBindings'],
        backgroundStrong=background['strong'],physicalConstantStrong=source,nativeStrong=native,
        strongCorrespondenceResidual=source_join,sourceUnresolved=unresolved,remainingOperators=actual_atoms,
        dependencySnapshot=manifest['dependencySnapshot'],sourcePins=manifest['sourcePins'],
        applicability='CONDITIONAL_ON_UNSUPPLIED_NONUNIFORM_AND_DIRECT_MIXED_GRADE_COMPOSITION',
        absentCandidateIsNotVanishingProof=True,restBulkVelocity=sp.Integer(0))
    profile_join=B.native['profileBindings']==endpoints
    actual['nativeProfileJoin']=profile_join
    J.check(B.label+'/applicability/source-correspondence',source_join is not None and zeros(source_join) and profile_join,actual)
    J.check(B.label+'/applicability/actual-constant-source',not actual_atoms and not any(unresolved.values()),actual,'UNRESOLVED')
    return actual


def exact_limits(B,R,N,sign,J):
    """Source-local wave-constrained paths, with scalar numerator/order evidence."""
    prefix=B.label+'/limits/sign-'+str(sign);t=sp.Symbol('uniformDepthApproach',positive=True)
    k=sign*sp.sqrt(R['k23']);K2=sp.simplify(k*k+R['kt2'])
    selected=projected_objects(R,N);results={};all_orders={}
    raw=dict(selected=selected,**N['rawObjects'],nativeDomainFactors=N['domains'],restrictionDomains=R['reductionDomains'],
        impedance=[[item['impedance'] for item in face] for face in N['closures']])
    for name in ('RADIATING','EVANESCENT'):
        depth=t if name=='RADIATING' else sp.I*t
        speed=3/sp.sqrt(K2+(t*t if name=='RADIATING' else -t*t))
        mapping={B.w:sp.Integer(3),B.cs:speed,B.k:k,B.km:k,B.q:depth,B.qm:sp.conjugate(depth)}
        mapping[B.Q]=R['conversion'].subs(mapping)*depth
        wave_residuals=transform((R['wave'],R['waveMinus']),lambda expr:sp.simplify(expr.subs(mapping)))
        J.check(prefix+'/'+name+'/wave-path',zeros(wave_residuals),dict(parameter=t,cs=speed,mapping=mapping,
            domain=(sp.Gt(t,0),sp.Lt(t,sp.sqrt(K2))) if name=='EVANESCENT' else (sp.Gt(t,0),),
            actualQ=mapping[B.Q],waveResiduals=wave_residuals,cT=R['cT']))
        cache={};orders=[]
        def limit_scalar(expr):
            raw_path=expr.subs(mapping)
            if raw_path in cache:return cache[raw_path]
            raw_numerator,raw_denominator=sp.fraction(sp.together(raw_path))
            J.emit(prefix+'/'+name+'/scalar-%04d/raw'%len(orders),dict(source=expr,rawPath=raw_path,
                rawNumerator=raw_numerator,rawDenominator=raw_denominator,originalFactors=factors(expr)))
            raw_leading_n=raw_numerator.as_leading_term(t);raw_leading_d=raw_denominator.as_leading_term(t)
            raw_norder=sp.oo if raw_numerator==0 else raw_leading_n.as_powers_dict().get(t,sp.Integer(0))
            raw_dorder=raw_leading_d.as_powers_dict().get(t,sp.Integer(0))
            value=sp.cancel(raw_path)
            numerator,denominator=sp.fraction(sp.together(value))
            if value==0:
                leading_n=sp.Integer(0);leading_d=denominator.as_leading_term(t)
                norder=sp.oo;dorder=leading_d.as_powers_dict().get(t,sp.Integer(0));result=sp.Integer(0)
            else:
                leading_n=numerator.as_leading_term(t);leading_d=denominator.as_leading_term(t)
                norder=leading_n.as_powers_dict().get(t,sp.Integer(0));dorder=leading_d.as_powers_dict().get(t,sp.Integer(0))
                result=sp.limit(value,t,0,dir='+')
            evidence=dict(source=expr,rawPath=raw_path,rawNumerator=raw_numerator,rawDenominator=raw_denominator,
                rawNumeratorOrder=raw_norder,rawDenominatorOrder=raw_dorder,pathValue=value,numerator=numerator,denominator=denominator,
                numeratorLeading=leading_n,denominatorLeading=leading_d,numeratorOrder=norder,
                denominatorOrder=dorder,limit=result)
            J.emit(prefix+'/'+name+'/scalar-%04d'%len(orders),evidence)
            orders.append(evidence);cache[raw_path]=result;return result
        result=transform(raw,limit_scalar);results[name]=result;all_orders[name]=orders
        # Full-P finiteness is diagnostic. Only the source-selected restriction,
        # face/current projections and unique native closure are support gates.
        selected_finite=finite(result['selected'])
        closure_values=(result['closureCoefficients'],result['closureDeterminants'])
        J.check(prefix+'/'+name+'/selected-finite-unique',selected_finite and finite(closure_values)
                and all(sp.simplify(item).is_zero is False for item in leaves(closure_values)),
            dict(selected=result['selected'],closureValues=closure_values,
                 fullPFinite=finite(result['P']),fullFormFinite=finite(result['forms']),
                 impedance=result['impedance'],orders=orders),'UNRESOLVED')
    difference=transform(results['RADIATING']['selected'],lambda expr:expr)
    def subtract(a,b):
        if isinstance(a,dict):return {key:subtract(value,b[key]) for key,value in a.items()}
        if isinstance(a,(tuple,list)):return tuple(subtract(x,y) for x,y in zip(a,b))
        return clean(a-b)
    difference=subtract(results['RADIATING']['selected'],results['EVANESCENT']['selected'])
    J.check(prefix+'/two-route-selected-join',zeros(difference),dict(radiating=results['RADIATING'],
        evanescent=results['EVANESCENT'],selectedResiduals=difference,
        evidenceGrade='SOURCE_INTERNAL_LIMIT_CONSISTENCY'),'UNRESOLVED')
    result=dict(status='COMPUTED',normal=k,paths=results,orders=all_orders,selectedJoin=difference)
    J.report(prefix.replace('/','-'),dict(normal=k,selectedJoin=difference,
        fullPFinite={name:finite(value['P']) for name,value in results.items()},
        scalarLimits={name:len(value) for name,value in all_orders.items()}))
    return result


def domain_at_point(N,mapping,J,prefix,grazing=False):
    evidence={};bad=[];magnitudes=[]
    for name,values in N['domains'].items():
        rows=[]
        for value in values:
            bound=sp.simplify(value.subs(mapping));row=dict(factor=value,value=bound)
            if not finite(bound) or bound.is_zero is not False:bad.append((name,row))
            elif not bound.free_symbols:
                magnitude=abs(complex(bound.evalf(40)));magnitudes.append(magnitude);row['magnitude']=magnitude
            rows.append(row)
        evidence[name]=rows
    J.emit(prefix+'/raw-source-domain',dict(factors=evidence,excludedAtDirectSubstitution=bad,
        sourceLocalLimitUsed=grazing,minimumNonzeroMagnitude=min(magnitudes,default=None)))
    if bad and not grazing:raise ScientificIssue('nongrazing source denominator outside admitted domain',bad)
    return dict(minimumNonzeroMagnitude=min(magnitudes,default=None),directSubstitutionExclusions=bad)


def evaluate_sign(B,R,N,speed,sign,limits,omissions,J,prefix,probe_sheet):
    mapping,stratum=physical_substitution(B,R,speed,sign);grazing=stratum['stratum']=='EXACT_GRAZING'
    domain=domain_at_point(N,mapping,J,prefix,grazing)
    actual=limits['paths']['RADIATING'] if grazing else transform(N['rawObjects'],lambda expr:sp.cancel(expr.subs(mapping)))
    selected=(actual['selected'] if grazing else transform(projected_objects(R,N),lambda expr:sp.cancel(expr.subs(mapping))))
    L=sp.ImmutableMatrix(R['lift'].subs(mapping));Ln=numeric(L);rank=int(np.linalg.matrix_rank(Ln,tol=1e-9*max(1,norm(Ln))))
    wave=transform((R['wave'],R['waveMinus']),lambda expr:sp.simplify(expr.subs(mapping)))
    J.check(prefix+'/selected-rank',rank==2,dict(lift=L,rank=rank),'UNRESOLVED')
    J.check(prefix+'/selected-source-residuals',zeros((selected['pencil'],selected['minusPencil'],wave)),
        dict(mapping=mapping,stratum=stratum,lift=L,rank=rank,fullResiduals=(selected['pencil'],selected['minusPencil']),
             waveResiduals=wave,units=R['units'],equalDepthSource=N['equalDepth'] if stratum['stratum']!='EVANESCENT' else None))
    G=numeric(selected['forms']['SLAB_CURRENT_MATRIX']);hermiticity=norm(G-G.conj().T)/max(1,norm(G))
    eigenvalues,U=np.linalg.eigh((G+G.conj().T)/2);floor=1e-9*max(1,norm(G))
    J.check(prefix+'/physical-current-hermiticity',hermiticity<1e-9,dict(Gram=G,hermiticity=hermiticity))
    J.check(prefix+'/physical-current-spectrum',all(abs(value)>floor for value in eigenvalues),
        dict(Gram=G,eigenvalues=eigenvalues,eigenvectors=U,hermiticity=hermiticity,eigenvalueFloor=floor),'UNRESOLVED')
    T=U@np.diag(1/np.sqrt(np.abs(eigenvalues)));normalized=Ln@T
    normalizedGram=T.conj().T@G@T;normalizationResidual=normalizedGram-np.diag(np.sign(eigenvalues))
    direction=sp.simplify(R['direction'].subs(mapping));direction_numeric=complex(direction.evalf(40))
    orientation=-1 if B.label=='LEFT' else 1
    J.check(prefix+'/normalization-and-direction',norm(normalizationResidual)<1e-9 and
        abs(direction_numeric.imag)<1e-12 and direction_numeric.real!=0 and
        all(np.sign(eigenvalues)==np.sign(direction_numeric.real)),
        dict(Gram=G,basis=normalized,fieldToCurrent=T,normalizedGram=normalizedGram,
            normalizationResidual=normalizationResidual,direction=direction,endOrientation=orientation,
            outwardCurrentSigns=orientation*np.sign(eigenvalues),outwardGroupSign=orientation*np.sign(direction_numeric.real)))
    face_records=[];face_ok=True;conjugacy_ok=True
    for face,rows in enumerate(actual['faces']):
        pair=[]
        for leg,rowmap in enumerate(rows):
            values={}
            for name,row in rowmap.items():
                projected=numeric(selected['faces'][face][leg][name])@(T if leg==0 else T.conj())
                row_finite=finite(row)
                scale=max(1,(norm(numeric(row))*norm(normalized) if row_finite else norm(projected)))
                relative=norm(projected)/scale
                values[name]=dict(exact=selected['faces'][face][leg][name],exactZero=zeros(selected['faces'][face][leg][name]),
                    normalizedProjection=projected,scale=scale,relative=relative,fullRowFinite=row_finite,
                    scaleSource='FULL_NATIVE_ROW' if row_finite else 'FINITE_SELECTED_LIMIT')
                face_ok=face_ok and relative<1e-8
            pair.append(values)
        conjugacy={name:pair[1][name]['normalizedProjection']-pair[0][name]['normalizedProjection'].conj() for name in PORTS}
        conjugacy_scales={name:max(1,norm(pair[0][name]['normalizedProjection']),norm(pair[1][name]['normalizedProjection'])) for name in PORTS}
        full_conjugacy={name:numeric(rows[1][name])-numeric(rows[0][name]).conj() if finite((rows[0][name],rows[1][name])) else None for name in PORTS}
        conjugacy_ok=conjugacy_ok and all(norm(value)/conjugacy_scales[name]<1e-9 for name,value in conjugacy.items())
        conjugacy_ok=conjugacy_ok and all(value is None or norm(value)/max(1,norm(numeric(rows[0][name])),norm(numeric(rows[1][name])))<1e-9
            for name,value in full_conjugacy.items())
        face_records.append(dict(face=face,legs=pair,selectedConjugacy=conjugacy,fullRowConjugacy=full_conjugacy,conjugacyScales=conjugacy_scales))
    projections={name:T.conj().T@numeric(value)@T for name,value in selected['forms'].items() if name!='SLAB_CURRENT_MATRIX'}
    projection_scales={name:max(1,norm(numeric(actual['forms'][name]))*norm(normalized)**2) if finite(actual['forms'][name])
        else max(1,norm(numeric(selected['forms'][name]))*norm(T)**2) for name in projections}
    J.check(prefix+'/faces-and-separated-forms',face_ok and conjugacy_ok and
        all(norm(value)/projection_scales[name]<1e-8 for name,value in projections.items()),
        dict(faces=face_records,forms=projections,scales=projection_scales,
             roles={'SLAB_CURRENT_MATRIX':'normal slab current','BULK_NORMAL_CURRENT_DENSITY_MATRIX':'bulk normal density',
                    'BULK_DEPTH_CURRENT_MATRIX':'bulk depth current','INTERFACE_POWER_MATRIX':'interface power'}))
    fullP=numeric(actual['P']) if finite(actual['P']) else None
    finite_coefficients=[complex(value.evalf(40)) for value in actual['P'] if finite(value)]
    fullResidual=numeric(selected['pencil'])@T
    residualScale=max(1,max((abs(value) for value in finite_coefficients),default=0)*norm(normalized))
    J.check(prefix+'/normalized-original-pencil',norm(fullResidual)/residualScale<1e-9,
            dict(P=fullP,fields=normalized,residual=fullResidual,scale=residualScale))
    # Address a native coefficient, omit it, and repeat the original matrix route.
    pcontrols=[]
    for i in range(5):
        for j in range(5):
            coefficient=actual['P'][i,j]
            if coefficient==0 or not finite(coefficient):continue
            delta=np.zeros_like(fullResidual);delta[i,:]=-complex(coefficient.evalf(40))*normalized[j,:]
            changed=None if fullP is None else fullP.copy()
            if changed is not None:changed[i,j]=0
            movement=norm(delta)/residualScale
            pcontrols.append(dict(index=(i,j),source=R['P'][i,j],original=coefficient,changed=changed,
                route='OMIT_FINITE_NATIVE_COEFFICIENT_THEN_REPROJECT_FULL_FIVE_ROW_SOURCE',
                residual=fullResidual+delta,movement=movement))
    best_p=max(pcontrols,key=lambda item:item['movement'],default=dict(movement=0))
    M=numeric(actual['forms']['SLAB_CURRENT_MATRIX']) if finite(actual['forms']['SLAB_CURRENT_MATRIX']) else None
    current_controls=[]
    for i in range(5):
        for j in range(i+1,5):
            entries=(actual['forms']['SLAB_CURRENT_MATRIX'][i,j],actual['forms']['SLAB_CURRENT_MATRIX'][j,i])
            if not finite(entries) or all(value==0 for value in entries):continue
            contribution=np.zeros((5,5),dtype=complex);contribution[i,j],contribution[j,i]=(complex(value.evalf(40)) for value in entries)
            changed=None if M is None else M.copy()
            if changed is not None:changed[i,j]=changed[j,i]=0
            changedGram=normalizedGram-normalized.conj().T@contribution@normalized
            movement=norm(changedGram-normalizedGram)/max(1,norm(normalizedGram))
            current_controls.append(dict(indices=((i,j),(j,i)),source=(N['forms']['SLAB_CURRENT_MATRIX'][i,j],
                N['forms']['SLAB_CURRENT_MATRIX'][j,i]),changed=changed,Gram=changedGram,movement=movement))
    best_current=max(current_controls,key=lambda item:item['movement'],default=dict(movement=0))
    reversed_signs=-orientation*np.sign(eigenvalues);independent_expected=orientation*np.sign(direction_numeric.real)
    orientation_control=dict(assigned=orientation,reversed=-orientation,baseline=orientation*np.sign(eigenvalues),
        reversedCurrentSigns=reversed_signs,independentGroupSign=independent_expected,
        disagreement=bool(np.all(reversed_signs!=independent_expected)))
    ew=[]
    for face in range(2):
        for leg in range(2):
            item=omissions[face][leg];before=numeric(item['baseline'].subs(mapping));after=numeric(item['changed'].subs(mapping))
            movement=norm((after-before)[:,4:5])/max(1,norm(before))
            ew.append(dict(face=face,leg=leg,baseline=before,omitted=after,unitProbeColumn=4,movement=movement))
    controls=dict(pencil=best_p,current=best_current,orientation=orientation_control,outwardVelocity=ew)
    J.check(prefix+'/responsive-controls',best_p['movement']>1e-10 and best_current['movement']>1e-10 and
        orientation_control['disagreement'] and all(item['movement']>1e-10 for item in ew),controls,'UNRESOLVED')
    if probe_sheet:
        require(not grazing,'sheet probe requires a declared nongrazing point')
        reversed_mapping={**mapping,B.q:-mapping[B.q],B.qm:-mapping[B.qm],B.Q:-mapping[B.Q]}
        sheet=[]
        for face in range(2):
            for leg in range(2):
                rows=N['faces'][face][leg]
                for column in (3,4):
                    probes={name:dict(baseline=numeric(row.subs(mapping))[:,column:column+1],
                        reversed=numeric(row.subs(reversed_mapping))[:,column:column+1]) for name,row in rows.items()}
                    movement=max(norm(item['baseline']-item['reversed'])/max(1,norm(item['baseline'])) for item in probes.values())
                    original_closure=N['closures'][face][leg]['closure']
                    changed_amplitude=(rows['AMPLITUDE'].subs(reversed_mapping)*sp.ImmutableMatrix(B.amplitudes[leg]))[0]
                    raw_closure=original_closure.subs(reversed_mapping).subs(sp.Symbol('closureAmplitude'),changed_amplitude)
                    closure=clean(raw_closure)
                    sheet.append(dict(face=face,leg=leg,column=column,rows=probes,movement=movement,
                        sourceClosure=original_closure,changedAmplitude=changed_amplitude,rawClosure=raw_closure,closureResidual=closure))
        J.check(prefix+'/nongrazing-native-sheet-probes',max(item['movement'] for item in sheet)>1e-10 and
            zeros([item['closureResidual'] for item in sheet]),dict(mapping=mapping,reversedMapping=reversed_mapping,probes=sheet),'UNRESOLVED')
        controls['sheet']=sheet
    singular=None if fullP is None else np.linalg.svd(fullP,compute_uv=False)
    rank_tolerance=None if fullP is None else 1e-9*max(1,norm(fullP))
    fullrank=None if fullP is None else int(sum(singular>rank_tolerance))
    limit_selected=limits['paths']['RADIATING']['selected'] if limits is not None else None
    approach=({name:clean(selected['forms'][name]-limit_selected['forms'][name]) for name in N['forms']}
              if limit_selected is not None else dict(status='UNRESOLVED',reason='exact limit unavailable'))
    return dict(status='SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE',end=B.label,speed=speed,sign=sign,
        mapping=mapping,stratum=stratum,rank=rank,fullRankEstimate=fullrank,fullNullityEstimate=None if fullrank is None else 5-fullrank,
        fullRankTolerance=rank_tolerance,fullSingularValues=singular,isolatedDoubletDiagnostic=fullrank==3,
        completeOutgoingBasisClaim=False,currentGram=G,currentEigenvalues=eigenvalues,normalizedBasis=normalized,
        normalizedGram=normalizedGram,direction=direction,faces=face_records,separateForms=projections,controls=controls,
        sourceDomain=domain,liftCondition=float(np.linalg.cond(Ln)),currentCondition=float(np.linalg.cond(G)),
        approachToExactLimit=approach,conditionalApplicability=True)


def containment(manifest):
    policy=manifest['resources'];memory=policy['memoryBytes']
    require(memory==8*1024**3 and policy['cpuCount']==1 and policy['nativeThreads']==1 and
            policy['tasksMax']==32 and policy['durationLimits'] is None,'approved native resource schema')
    group=next(line[3:] for line in Path('/proc/self/cgroup').read_text().splitlines() if line.startswith('0::'))
    cgroup=Path('/sys/fs/cgroup')/group.lstrip('/')
    values={name:(cgroup/name).read_text().strip() for name in ('memory.max','memory.swap.max','pids.max')}
    values.update(affinity=sorted(os.sched_getaffinity(0)),threads={name:os.environ.get(name) for name in THREADS},
                  cgroup=str(cgroup),nativeAddressSpace=memory,durationLimits=None)
    require(values['memory.max']==str(memory) and values['memory.swap.max']=='0' and values['pids.max']=='32'
        and len(values['affinity'])==1 and all(value=='1' for value in values['threads'].values()),
        'pooled whole-job memory/cpu/task/thread containment missing; scientific imports refused')
    resource.setrlimit(resource.RLIMIT_AS,(memory,memory));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    require(resource.getrlimit(resource.RLIMIT_CPU)==(resource.RLIM_INFINITY,resource.RLIM_INFINITY),
            'unexpected inherited CPU deadline')
    return values


def gate_check(args):
    manifest=json.loads(args.inputs.read_text());gate=json.loads(args.gate.read_text())
    require(gate['status']=='READY_FOR_SELECTED_NEAR_UNITY_UNIFORM' and gate['methodClearance'] is True and
            gate['pooledExecution'] is True,'fresh selected-method launch gate required')
    require(gate['workerSha256']==digest(__file__) and gate['manifestSha256']==digest(args.inputs)
            and gate['planSha256']==PLAN_SHA,'exact worker/manifest/plan gate')
    require(gate['outputDirectory']==str(args.out.resolve()) and gate['sourcePins']==manifest['sourcePins'],
            'gate output and full source correspondence')
    require(gate['methodReviewRecord']==manifest['reviewRecord'],'method clearance record gate')
    require(set(manifest['packets'])==set(PACKETS) and manifest['maximumSpeeds']==12 and
            manifest['maximumEndEvaluations']==24,'original nine input routes and bounded schedule')
    require(manifest['plan']['sha256']==PLAN_SHA and manifest['physicalInput']['sha256']==INPUT_SHA,'approved plan and physical input')
    require(manifest['sourcePins'].get(str(Path(__file__).resolve()))==digest(__file__),'worker included in source pins')
    for pin in (*manifest['packets'].values(),manifest['physicalInput'],manifest['plan'],manifest['reviewRecord']):verify(pin)
    for path,sha in manifest['sourcePins'].items():require(digest(path)==sha,'changed source '+path)
    review=json.loads(Path(manifest['reviewRecord']['path']).read_text())
    require(review['status']=='INDEPENDENT_METHOD_CLEARANCE_BOTH_REVIEWERS_NO_RUNTIME_RESULT' and
            review['allChecksPassed'] is True,'completed exact method record')
    args.out.resolve().relative_to(ROOT/'_scratch/s11c')
    require(not args.out.exists(),'fresh one-run output required; no automatic retry')
    return manifest,gate


def science(manifest,J):
    restored={}
    for name in PACKETS:
        pin=manifest['packets'][name];raw=Path(pin['path']).read_bytes()
        rawref=J.store.put('original-inputs/'+name+'.pickle',raw)
        J.event(dict(name=name,status='ORIGINAL_BYTES_SAVED',value=rawref,provenance=pin))
        def restore():
            value=decode(raw)
            if name=='uniformResponse':value={key:value[key] for key in ('fieldUnits','currentUnit')}
            return value
        restored[name]=J.call('restore/'+name,dict(pin=pin,originalBytes=rawref),restore)
    specification=json.loads(Path(manifest['physicalInput']['path']).read_text())
    params={key:sp.Rational(value) for key,value in specification['parameters'].items()}
    require(params['omega']==1 and params['c_s0']==10 and params['s11cdTangentialMomentum1']==sp.Rational(1,5)
        and params['s11cdTangentialMomentum2']==sp.Rational(1,10),'original undeformed material/input correspondence')
    params['omega']=sp.Integer(3)
    w,cs=sp.symbols('uniformFrequency uniformSoundSpeed',positive=True)
    k,km=sp.symbols('uniformNormal uniformMinusNormal',real=True)
    q,qm,Q=sp.symbols('uniformPhysicalDepth uniformMinusPhysicalDepth uniformAlgebraicRadical')
    engine_path=next(Path(path) for path in manifest['sourcePins'] if path.endswith('/S11c_d_mixing_scattering_sympy_audit.py'))
    tree=ast.parse(engine_path.read_text())
    L=J.call('selected-lift',dict(uniform=restored['uniformSource']['curl'],common=restored['uniformCommon'],
            parameters=params,source=manifest['sourcePins'][str(engine_path)]),
            lambda:selected_lift(restored['uniformSource'],restored['uniformCommon'],params,k,tree,J))
    L['fieldUnits']=restored['uniformResponse']['fieldUnits'];L['currentUnit']=restored['uniformResponse']['currentUnit']
    J.emit('physical-units',dict(field=L['fieldUnits'],current=L['currentUnit'],frame=specification['unit_frame']))
    # Recompute approved profile endpoints and jets directly from the physical
    # input, then join each actual native end-limit placeholder.
    xi=sp.Symbol('uniformProfileCoordinate',real=True)
    profiles={name:sp.sympify(expression,locals={'xi':xi,'tanh':sp.tanh}) for name,expression in specification['profiles'].items()}
    profile_joins=[]
    for label,end in (('LEFT',-sp.oo),('RIGHT',sp.oo)):
        operands=restored['uniformSource']['records'][label]['background']['profileLimitOperands']
        for original,limit,endpoint in operands:
            name=next(name for name in profiles if original.func.__name__=='s11cd'+name.upper()+'Profile')
            evaluated=sp.limit(profiles[name],xi,end)
            jets=tuple(sp.limit(sp.diff(profiles[name],xi,order),xi,end) for order in (1,2))
            supplied=restored['uniformSource']['profileEndpoints'][endpoint]
            profile_joins.append(dict(end=label,source=(original,limit,endpoint),definition=profiles[name],
                actualEnd=evaluated,supplied=supplied,residual=evaluated-supplied,jets=jets))
    J.check('profiles/constant-end-and-zero-jets',all(item['residual']==0 and zeros(item['jets']) for item in profile_joins),
        dict(profileJoins=profile_joins,translationDerivativeAtConstantEnds=tuple(item['jets'][0] for item in profile_joins),
             missingDirectMixedGradeCandidate='NO_VANISHING_INFERENCE_FROM_ABSENCE'))
    ends={}
    for label in ('LEFT','RIGHT'):
        B=EndBinding(label,restored,params,w,cs,k,km,q,qm,Q,J)
        R=J.call(label+'/restriction',dict(actualSource=B.r['CLOSED_PENCIL_LEGS'],frequency=B.frequency,lift=L,parameters=params),
                 lambda:end_restriction(B,L,J))
        B.fixed_frequency=sp.Integer(3)
        applicability=J.call(label+'/applicability',dict(source=B.uniform['records'][label],native=B.native['profileBindings']),
                             lambda:uniform_applicability(B,manifest,J))
        N=J.call(label+'/native-reconstruction',dict(pairing=B.packet,native=B.native,parameters=params,frequency=3,cs=cs),
                 lambda:native_reconstruction(B,R,tree,J))
        omissions=J.call(label+'/native-source-controls',dict(nativeFaces=B.r['FACE_LEG_OBJECTS']),lambda:source_controls(B,N,J))
        ends[label]=dict(B=B,R=R,N=N,omissions=omissions,applicability=applicability,limits={})
    speeds=[]
    for label in ('LEFT','RIGHT'):
        ratios=tuple(map(sp.Rational,('199/200','1','201/200','1999/2000','2001/2000')))
        if label=='LEFT':ratios+=tuple(map(sp.Rational,('49/50','51/50')))
        for ratio in ratios:
            value=sp.simplify(ends[label]['R']['cT']/ratio)
            if not any(sp.simplify(value-item)==0 for item in speeds):speeds.append(value)
    require(len(speeds)<=manifest['maximumSpeeds'] and 2*len(speeds)<=manifest['maximumEndEvaluations'],'bounded exact speed union')
    schedule=[dict(index=i,speed=value,modalRatios={label:sp.simplify(item['R']['cT']/value) for label,item in ends.items()},
                   bareRatio=sp.simplify(1/value)) for i,value in enumerate(speeds)]
    declared_sheet_probes={label:next(entry['index'] for entry in schedule if sp.simplify(entry['speed']-item['R']['cT'])!=0)
                           for label,item in ends.items()}
    J.emit('exact-speed-schedule',dict(schedule=schedule,declaredSheetProbes=dict(indices=declared_sheet_probes,normalSign=1,columns=(3,4)),
        dependencySnapshot=manifest['dependencySnapshot'],
        physicalInput=manifest['physicalInput'],endEvaluations=2*len(speeds),csBinding='DIRECT_EXACT_POSITIVE_ALGEBRAIC'))
    for label,item in ends.items():
        for sign in (-1,1):
            item['limits'][sign]=J.call(label+'/exact-limits/'+str(sign),dict(restriction=item['R'],native=item['N'],sign=sign),
                lambda item=item,sign=sign:exact_limits(item['B'],item['R'],item['N'],sign,J),soft=True)
    rows=[];probed=set()
    for entry in schedule:
        for label,item in ends.items():
            results=[];speed=entry['speed'];own_match=sp.simplify(speed-item['R']['cT'])==0
            for sign in (-1,1):
                prefix='points/%02d/%s/%s'%(entry['index'],label,sign)
                limit=item['limits'][sign]
                if own_match and limit.get('status')!='COMPUTED':
                    result=dict(status='UNRESOLVED',reason='exact selected limit unavailable',limitStatus=limit)
                    J.emit(prefix+'/unresolved-dependency',result)
                else:
                    probe=entry['index']==declared_sheet_probes[label] and sign==1
                    result=J.call(prefix,dict(schedule=entry,sign=sign,limit=limit,actualSource=item['N'],probeSheet=probe),
                        lambda:evaluate_sign(item['B'],item['R'],item['N'],speed,sign,
                            limit if limit.get('status')=='COMPUTED' else None,item['omissions'],J,prefix,probe),soft=True)
                    if probe and result.get('status')=='SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE':probed.add(label)
                results.append(result)
            statuses=[result['status'] for result in results]
            status=('FAILED_CHECK' if 'FAILED_CHECK' in statuses else 'SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE'
                    if statuses==['SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE']*2 else 'UNRESOLVED')
            row=dict(schedule=entry,end=label,status=status,momentumSigns=results,
                applicability=item['applicability']['applicability'],dependencySnapshot=manifest['dependencySnapshot'])
            rows.append(row);J.report('point-%02d-%s'%(entry['index'],label),row)
    statuses=[row['status'] for row in rows]
    overall=('FAILED_CHECK' if 'FAILED_CHECK' in statuses else 'SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE'
             if all(value=='SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE' for value in statuses) and probed=={'LEFT','RIGHT'} else 'UNRESOLVED')
    result=dict(status=overall,exactSpeeds=schedule,endEvaluations=len(rows),rows=rows,
        currentUnit=L['currentUnit'],fieldUnits=L['fieldUnits'],nativeSheetProbeEnds=sorted(probed),
        dependencySnapshot=manifest['dependencySnapshot'],scope='SELECTED_CONSTANT_END_TRANSVERSE_DOUBLETS_ONLY',
        applicability='CONDITIONAL_ON_UNSUPPLIED_NONUNIFORM_AND_DIRECT_MIXED_GRADE_COMPOSITION',
        evidenceGrade='ORIGINAL_SOURCE_RECONSTRUCTIONS_AND_SELECTED_NUMERICAL_CHECKS',
        restBulkVelocity=0,completeOutgoingBasisClaim=False,leakageOrLossClaim=False,physicalCalibrationClaim=False)
    J.emit('selected-uniform-result',result);return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--inputs',type=Path,required=True);parser.add_argument('--gate',type=Path,required=True)
    args=parser.parse_args();manifest,gate=gate_check(args);initial_gate_hash=digest(args.gate);limits=containment(manifest)
    global sp,np
    import sympy as sp
    import numpy as np
    args.out.mkdir(parents=True,exist_ok=False);J=Journal(args.out)
    save(args.out/'invocation.json',dict(arguments=vars(args)|{'out':str(args.out),'inputs':str(args.inputs),'gate':str(args.gate)},
        workerSha256=digest(__file__),manifestSha256=digest(args.inputs),gateSha256=digest(args.gate),containment=limits))
    result=None;exit_code=0
    try:
        result=J.call('uniform-science',dict(manifest=manifest,gate=gate,containment=limits),lambda:science(manifest,J),soft=True)
        if result.get('status')!='SUPPORTED_ON_SELECTED_UNIFORM_SUBSPACE':exit_code=2
    except BaseException as error:
        result=dict(status='FAILED_CHECK',reason=str(error),exceptionType=type(error).__name__,traceback=traceback.format_exc(),
                    activeOperations=list(J.active));J.emit('fatal-failure',result);exit_code=2
    finally:
        post={path:digest(path)==sha for path,sha in manifest['sourcePins'].items()}
        post.update({pin['path']:digest(pin['path'])==pin['sha256'] for pin in manifest['packets'].values()})
        post[str(args.inputs.resolve())]=digest(args.inputs)==gate['manifestSha256']
        post[str(args.gate.resolve())]=digest(args.gate)==initial_gate_hash
        J.emit('posthashes',post);J.store.integrity_check();count=J.store.count();J.store.close()
        if not all(post.values()):result=dict(status='FAILED_CHECK',reason='source/input posthash changed',computed=result);exit_code=2
        checks=dict(status=result['status'],result=result,posthashes=post,blobCount=count,completedOperations=J.count,
                    scientificFailures=J.failures,containment=limits,manifestSha256=digest(args.inputs),workerSha256=digest(__file__))
        rendered=(json.dumps(readable(checks),indent=2,allow_nan=False)+'\n')
        with (args.out/'checks.json').open('x') as stream:stream.write(rendered);stream.flush();os.fsync(stream.fileno())
        sys.stdout.write(rendered);sys.stdout.flush()
    return exit_code


if __name__=='__main__':
    raise SystemExit(main())
