#!/usr/bin/env python3
"""Numerical pilot's first required end-continuation gate; no producer replay.

Reads four pinned accepted packets inside the shared resource guard. Implements
only contrast/frequency continuation and its numerical checks. Field assembly,
face/current checks and finite deficits are later pilot work, never inferred
from this gate. Full inputs/returns are saved in an append-only zip journal.
"""
import argparse
import ast
import hashlib
import io
import json
import math
import os
from pathlib import Path
import pickle
import resource
import shutil
import signal
import sys
import time
import traceback
import warnings
import zipfile


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def save_json(p, value):
    tmp=p.with_suffix(p.suffix+'.new')
    with tmp.open('w') as out:
        json.dump(value,out,indent=2,allow_nan=False);out.write('\n')
        out.flush();os.fsync(out.fileno())
    tmp.replace(p)


def numeric_json(value):
    if isinstance(value,dict):return {str(k):numeric_json(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [numeric_json(v) for v in value]
    if isinstance(value,np.ndarray):return numeric_json(value.tolist())
    if isinstance(value,np.generic):return numeric_json(value.item())
    if isinstance(value,complex):return {'real':numeric_json(value.real),'imag':numeric_json(value.imag)}
    if isinstance(value,float) and not math.isfinite(value):return {'nonfinite':repr(value)}
    return value


def norm(value):return float(np.max(np.abs(value),initial=0.))


class Journal:
    def __init__(self,base):
        self.base=base;self.count=0;self.incomplete=None
        self.path=base/'operations.zip';self.names=set()
    def blob(self,name,value):
        require(name not in self.names,'unique saved operand: '+name)
        raw=pickle.dumps(value,protocol=5)
        with zipfile.ZipFile(self.path,'a',compression=zipfile.ZIP_STORED) as z:
            z.writestr(name,raw)
        self.names.add(name)
        return {'member':name,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
    def call(self,name,args,fn):
        self.incomplete=name;inp=self.blob(name+'/input.pickle',args)
        save_json(self.base/'active-operation.json',{'name':name,'input':inp,'status':'STARTED'})
        start=time.monotonic();value=fn();out=self.blob(name+'/return.pickle',value)
        receipt={'name':name,'input':inp,'return':out,'status':'COMPLETE','seconds':time.monotonic()-start}
        with (self.base/'operation-index.jsonl').open('a') as stream:
            stream.write(json.dumps(receipt)+'\n');stream.flush();os.fsync(stream.fileno())
        self.count+=1;self.incomplete=None
        save_json(self.base/'active-operation.json',receipt)
        return value
    def report(self,name,data):save_json(self.base/(name+'.json'),numeric_json(data))


def polynomial(terms,K,Q,dK=None,dQ=None):
    # Same ordered K^p Q^t rational-matrix action as the saved Pair equation.
    I=np.eye(len(K),dtype=complex);Z=np.zeros_like(I)
    value=Z.copy();delta=Z.copy()
    for (p,t),coefficient in terms:
        kp=np.linalg.matrix_power(K,p);qt=np.linalg.matrix_power(Q,t)
        dk=Z if dK is None or p==0 else sum((np.linalg.matrix_power(K,j)@dK@np.linalg.matrix_power(K,p-1-j) for j in range(p)),start=Z.copy())
        dq=Z if dQ is None or t==0 else sum((np.linalg.matrix_power(Q,j)@dQ@np.linalg.matrix_power(Q,t-1-j) for j in range(t)),start=Z.copy())
        value+=coefficient*(kp@qt);delta+=coefficient*(dk@qt+kp@dq)
    return value,delta


class Table:
    def __init__(self,pencil,wave,parameter,k,q):
        self.pencil=pencil;self.wave=wave;self.parameter=parameter;self.k=k;self.q=q
        require(pencil.shape==(5,5),'full 5x5 source pencil')
        require(not(pencil.free_symbols|wave.free_symbols)-{parameter,k,q},'only one declared live path parameter and two algebraic momenta')
        self.rows=[];expressions=[]
        for i in range(5):
            for j in range(5):
                numerator,denominator=sp.fraction(sp.cancel(pencil[i,j]))
                nt=sp.Poly(numerator,k,q).terms();dt=sp.Poly(denominator,k,q).terms()
                self.rows.append((i,j,[p for p,c in nt],[p for p,c in dt],len(expressions)))
                expressions.extend(c for p,c in nt+dt)
        wt=sp.Poly(wave,k,q).terms();self.wave_start=len(expressions);self.wave_powers=[p for p,c in wt]
        expressions.extend(c for p,c in wt)
        require(all(not c.free_symbols-{parameter} for c in expressions),'numeric material binding before coefficient compilation')
        self.fn=sp.lambdify(parameter,expressions,'numpy',cse=True)
        self.native=sp.lambdify((parameter,k,q),pencil,'numpy',cse=True)
        self.relation=sp.lambdify((parameter,k,q),wave,'numpy',cse=True)
        self.square=sp.solve(wave,q**2)
        require(len(self.square)==1,'one quadratic radical square')
        self.square=self.square[0]
        require(sp.Poly(self.square,k).degree()==2,'quadratic bound wave relation')
        self.evidence={'pencil':pencil,'wave':wave,'parameter':parameter,'k':k,'q':q,'entries':self.rows,'coefficientExpressions':expressions,'wavePowers':self.wave_powers,'radicalSquare':self.square}
        self.last=None
    def coefficients(self,z):
        if self.last is None or self.last[0]!=z:
            a=np.asarray(self.fn(z),complex).ravel();require(np.isfinite(a).all(),'finite bound coefficients')
            rows=[]
            for i,j,n,d,start in self.rows:
                rows.append((i,j,list(zip(n,a[start:start+len(n)])),list(zip(d,a[start+len(n):start+len(n)+len(d)]))))
            self.last=(z,rows,list(zip(self.wave_powers,a[self.wave_start:])))
        return self.last[1:]


class Pair:
    """Saved invariant-pair equations and exact K/Q Frechet derivatives.

    No frequency derivatives are constructed; continuation uses a secant
    predictor after its first saved-state bootstrap. Whole nullity-two clusters keep one fixed gauge.
    """
    def __init__(self,table,seed,z):
        self.table=table;self.n=seed['R'].shape[1];self.I=np.eye(self.n,dtype=complex)
        self.gauge=np.linalg.solve(seed['R'].conj().T@seed['R'],seed['R'].conj().T)
        k=complex(np.trace(seed['K'])/self.n);q=complex(np.trace(seed['Q'])/self.n)
        matrix=np.asarray(table.native(z,k,q),complex)
        self.row_scale=1+np.max(abs(matrix),axis=1)
        self.wave_scale=1+abs(q)**2+100*abs(k)**2
        self.scale=np.concatenate((np.full(5*self.n,max(1.,norm(seed['R']))),np.full(self.n**2,max(1.,abs(k))),np.full(self.n**2,max(1.,abs(q)))))
    def pack(self,R,K,Q):return np.concatenate((R.ravel(),K.ravel(),Q.ravel()))
    def unpack(self,x):
        n=self.n;return x[:5*n].reshape(5,n),x[5*n:5*n+n*n].reshape(n,n),x[5*n+n*n:].reshape(n,n)
    def equation(self,z,x,dx=None):
        R,K,Q=self.unpack(x);zero=np.zeros_like(K)
        dR,dK,dQ=(np.zeros_like(R),zero,zero) if dx is None else self.unpack(dx)
        value=np.zeros_like(R);delta=np.zeros_like(R);minimum=float('inf')
        rows,wave=self.table.coefficients(z)
        for i,j,nt,dt in rows:
            N,dN=polynomial(nt,K,Q,dK,dQ);D,dD=polynomial(dt,K,Q,dK,dQ)
            minimum=min(minimum,float(np.linalg.svd(D,compute_uv=False)[-1]))
            inv=np.linalg.solve(D,self.I);H=N@inv;dH=dN@inv-H@dD@inv
            value[i]+=R[j]@H;delta[i]+=dR[j]@H+R[j]@dH
        W,dW=polynomial(wave,K,Q,dK,dQ)
        return self.pack(value/self.row_scale[:,None],self.gauge@R-self.I,W/self.wave_scale),self.pack(delta/self.row_scale[:,None],self.gauge@dR,dW/self.wave_scale),minimum
    def jacobian(self,z,x):
        return np.column_stack([self.equation(z,x,e)[1] for e in np.eye(len(x),dtype=complex)])
    def correct(self,z,previous):
        x=self.pack(previous['R'],previous['K'],previous['Q']);history=[]
        for iteration in range(15):
            residual,_,margin=self.equation(z,x);J=self.jacobian(z,x)
            condition=float(np.linalg.cond(J*self.scale[None,:]))
            history.append({'iteration':iteration,'residual':norm(residual),'condition':condition,'denominatorMargin':margin})
            require(np.isfinite(x).all() and np.isfinite(condition),'finite invariant-pair state and Jacobian')
            if norm(residual)<2e-12:break
            step=np.linalg.solve(J*self.scale[None,:],-residual)*self.scale
            for damping in (1.,.5,.25,.125,.0625):
                trial=x+damping*step
                if norm(self.equation(z,trial)[0])<norm(residual):x=trial;break
            else:raise ValueError('BOUNDARY_UNRESOLVED: invariant-pair Newton decrease')
        residual,_,margin=self.equation(z,x);J=self.jacobian(z,x)
        R,K,Q=self.unpack(x);comm=K@Q-Q@K
        k=complex(np.trace(K)/self.n);q=complex(np.trace(Q)/self.n)
        native=np.asarray(self.table.native(z,k,q),complex);s=np.linalg.svd(native,compute_uv=False)
        nullity=int(np.sum(s<1e-8*max(1.,s[0])))
        return {'parameter':z,'R':R,'K':K,'Q':Q,'x':x,'residual':residual,'jacobian':J,'history':history,'rowScale':self.row_scale,'waveScale':self.wave_scale,'gauge':self.gauge,'denominatorMargin':margin,'commutator':comm,'scalarKResidual':K-k*self.I,'scalarQResidual':Q-q*self.I,'nativePencil':native,'nativeKernelResidual':native@R/(1+norm(native)*norm(R)),'nativeWaveResidual':complex(self.table.relation(z,k,q)),'singularValues':s,'nullity':nullity,'k':k,'q':q}
    def validate(self,state):
        require(norm(state['residual'])<2e-12 and state['denominatorMargin']>1e-12,'BOUNDARY_UNRESOLVED: corrected full pair or denominator')
        require(norm(state['commutator'])<1e-9*(1+norm(state['K'])*norm(state['Q'])),'BOUNDARY_UNRESOLVED: noncommuting pair')
        require(norm(state['scalarKResidual'])<1e-8*(1+abs(state['k'])) and norm(state['scalarQResidual'])<1e-8*(1+abs(state['q'])),'BOUNDARY_UNRESOLVED: ambiguous split doublet')
        require(state['nullity']==self.n,'BOUNDARY_UNRESOLVED: changed full-pencil nullity')
        require(norm(state['nativeKernelResidual'])<1e-9 and abs(state['nativeWaveResidual'])<1e-9*self.wave_scale,'BOUNDARY_UNRESOLVED: original-pencil/wave residual')


def projector(R):return R@np.linalg.solve(R.conj().T@R,R.conj().T)


def compare(a,b):
    return {'k':abs(a['k']-b['k'])/(1+abs(a['k'])),'q':abs(a['q']-b['q'])/(1+abs(a['q'])),'projector':norm(projector(a['R'])-projector(b['R']))}


def bind_source(source):
    # Keep the actual saved material map; alter only the authorized contrast
    # multiplier and frequency, never the physical permeability or its memory.
    k,q,w=source['momentum'],source['radical'],source['frequency']
    P=source['originalAlgebraic'].xreplace(source['mapping'])
    wave=source['originalRelation'].xreplace(source['mapping'])
    remaining=(P.free_symbols|wave.free_symbols)-{k,q}-set(source['origin'])
    require(len(remaining)==1,'one actual unbound original frequency')
    oldw=next(iter(remaining));a=sp.Symbol('radiatingContrastMultiplier',real=True)
    mapping={s:a*v for s,v in source['origin'].items()}
    P=P.subs(mapping,simultaneous=True).xreplace({oldw:w})
    wave=wave.subs(mapping,simultaneous=True).xreplace({oldw:w})
    return P,wave,a,w,k,q


def trace_path(J,prefix,pair,initial,start,end,step,joint,frequency_of):
    count=max(1,int(math.ceil(abs(end-start)/step)));state=initial;zold=start;previous=None
    for i in range(1,count+1):
        z=start+(end-start)*i/count;name=prefix+f'/point-{i:04d}'
        predictor=state
        if previous is not None:
            previous_z,previous_state=previous
            factor=(z-zold)/(zold-previous_z)
            predictor={key:state[key]+factor*(state[key]-previous_state[key]) for key in ('R','K','Q')}
        inputs={'startParameter':zold,'targetParameter':z,'previousCompleteState':state,'predictor':predictor,'predictorKind':'secant' if previous is not None else 'saved-state bootstrap','maximumStep':step}
        fresh=J.call(name,inputs,lambda:pair.correct(z,predictor))
        J.report('last-pair-check',{'operation':name,'parameter':z,'k':fresh['k'],'q':fresh['q'],'R':fresh['R'],'K':fresh['K'],'Q':fresh['Q'],'nativePencil':fresh['nativePencil'],'nativeKernelResidual':fresh['nativeKernelResidual'],'nativeWaveResidual':fresh['nativeWaveResidual'],'invariantPairResidual':fresh['residual'],'nullity':fresh['nullity'],'expectedNullity':pair.n,'singularValues':fresh['singularValues'],'denominatorMargin':fresh['denominatorMargin'],'commutator':fresh['commutator'],'scalarKResidual':fresh['scalarKResidual'],'scalarQResidual':fresh['scalarQResidual'],'history':fresh['history']})
        pair.validate(fresh)
        # Radical transport is independent of Newton's Q, with an explicit
        # saved seed for complex normal momenta. It follows the joint path.
        if joint is not None:
            vertices=[(frequency_of(zold),state['k']),(frequency_of(z),fresh['k'])]
            seed_q=state.get('transportedQ',state['q'])
            lift=J.call(name+'-joint-sheet',{'vertices':vertices,'seedQ':seed_q},lambda:joint.trace(vertices,seed=seed_q))
            gap=abs(lift.get('END_Q',complex('nan'))-fresh['q'])/(1+abs(fresh['q']))
            J.report('last-sheet-check',{'name':name,'status':lift['STATUS'],'relativeQDifference':gap})
            require(lift['PATH_DEFINED'] and gap<1e-8 and abs(lift.get('ODE_DIFFERENCE',complex('nan')))<1e-8*(1+abs(fresh['q'])),'BOUNDARY_UNRESOLVED: independent joint sheet transport')
            fresh=dict(fresh,transportedQ=lift['END_Q'])
        previous=(zold,state);state=fresh;zold=z
        J.report('last-completed-point',{'operation':name,'parameter':z,'k':fresh['k'],'q':fresh['q'],'nullity':fresh['nullity'],'residual':norm(fresh['residual']),'denominatorMargin':fresh['denominatorMargin']})
    return state


def compile_table(J,name,pencil,wave,parameter,k,q):
    holder=[]
    def construct():
        holder.append(Table(pencil,wave,parameter,k,q));return holder[0].evidence
    J.call(name,{'pencil':pencil,'wave':wave,'parameter':parameter,'k':k,'q':q},construct)
    return holder[0]


def original_seed(mode):
    info=mode['info'];R=np.asarray(mode['right'],complex);n=info['NULLITY']
    require(R.shape==(5,n),'complete saved right basis')
    k=complex(info['K']);q=complex(info['Q'])
    return {'R':R,'K':k*np.eye(n,dtype=complex),'Q':q*np.eye(n,dtype=complex),'k':k,'q':q,'parameter':1.,'sourceInfo':info}


def synthetic_check(J):
    k,q,z=sp.symbols('instrument_k instrument_q instrument_parameter')
    P=sp.diag((k-z)/(q+2),(k-z)/(q+2),1,2,3);wave=q*q-k*k-1
    table=Table(P,wave,z,k,q);seed={'R':np.eye(5,dtype=complex)[:,:2],'K':np.eye(2,dtype=complex),'Q':np.sqrt(2)*np.eye(2,dtype=complex)}
    pair=Pair(table,seed,1.);x=pair.pack(seed['R'],seed['K'],seed['Q'])
    direction=np.arange(1,len(x)+1,dtype=float)/len(x)+.3j
    exact=pair.equation(1.,x,direction)[1];h=1e-6
    finite=(pair.equation(1.,x+h*direction)[0]-pair.equation(1.,x-h*direction)[0])/(2*h)
    result=pair.correct(1.01,seed);pair.validate(result)
    data={'derivativeResidual':exact-finite,'state':result,'expectedK':1.01,'expectedQ':np.sqrt(1+1.01**2)}
    require(norm(data['derivativeResidual'])<1e-7 and abs(result['k']-1.01)<1e-10 and abs(result['q']-data['expectedQ'])<1e-10,'synthetic equation/derivative/whole-doublet instrument')
    return data


def run(base,manifest,J):
    J.call('instrument/synthetic',{'role':'no physical input, exact diagonal pencil'},lambda:synthetic_check(J))
    inputs={}
    for key,item in manifest['packets'].items():
        target=base/'inputs'/Path(item['path']).name;target.parent.mkdir(exist_ok=True)
        require(digest(item['path'])==item['sha256'],'pinned source packet '+key)
        shutil.copyfile(item['path'],target);require(digest(target)==item['sha256'],'byte-identical packet copy')
        inputs[key]=J.call('restore/'+key,{'path':str(target),'sha256':item['sha256']},lambda p=target:pickle.loads(p.read_bytes()))
    # Reuse precisely the saved quadratic-path helper definition, not a
    # producer module or its top-level constructors.
    source=Path(manifest['jointHelper']['path']).read_text();tree=ast.parse(source)
    node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='JointBulkSheetPath')
    namespace={'np':np,'sp':sp};exec(compile(ast.Module(body=[node],type_ignores=[]),manifest['jointHelper']['path'],'exec'),namespace)
    Joint=namespace['JointBulkSheetPath']
    uniform=inputs['uniform'];results={};physical_scales={}
    ref=inputs['reference'];ref_native=sp.lambdify((ref['frequency'],ref['momentum'],ref['radical']),ref['livePencil'],'numpy',cse=True)
    for label in ('LEFT','RIGHT'):
        modes=uniform['backgrounds'][label]['modes'];require(len(modes)==18,'all eighteen saved candidates')
        end=inputs[label.lower()];P,wave,a,w,k,q=J.call(label+'/source-binding',end,lambda:bind_source(end))
        J.blob(label+'/symbolic-unit-context.pickle',{'source':end,'fieldUnits':uniform['fieldUnits'],'currentUnit':uniform['currentUnit']})
        contrast_table=compile_table(J,label+'/contrast-table',P.subs(w,1),wave.subs(w,1),a,k,q)
        full_wave=wave.subs(a,1)
        require(not full_wave.has(a) and sp.expand(wave-full_wave)==0,'source wave independent of contrast')
        square=sp.solve(full_wave,q**2)[0]
        contrast_joint=Joint(full_wave,w,k,q,sp.sqrt(square))
        seeds={m['info']['INDEX']:original_seed(m) for m in modes}
        channels=uniform['backgrounds'][label]['response']['channels'][label]
        selected={v['RECORD_INDEX'] for v in channels['incoming']+channels['outgoing']}
        expected_out={v['RECORD_INDEX'] for v in channels['outgoing']}
        require(len(selected)==5 and len(channels['incoming'])==2 and len(channels['outgoing'])==5,'saved selected direction census')
        scale=[];joins=[]
        for mode in modes:
            index=mode['info']['INDEX'];s=seeds[index]
            raw=np.asarray(contrast_table.native(1.,s['k'],s['q']),complex)
            gap=norm(raw-mode['pencil'])/(1+norm(raw))
            kernel=norm(raw@s['R'])/(1+norm(raw)*norm(s['R']))
            ratio=s['q']/complex(mode['info']['PHYSICAL_Q']);scale.append(ratio)
            joins.append({'index':index,'sourceMatrixResidual':gap,'sourceKernelResidual':kernel,'algebraicToPhysicalScale':ratio})
        J.report(label+'-seed-joins',joins)
        require(all(v['sourceMatrixResidual']<1e-10 and v['sourceKernelResidual']<1e-9 for v in joins),'actual accepted pencil/basis seed joins')
        require(max(abs(v-scale[0]) for v in scale)<1e-10 and abs(scale[0].imag)<1e-12 and scale[0].real>0,'source physical radical scale')
        physical_scales[label]=scale[0].real
        contrast_states={1.:seeds};contrast_correspondence={}
        # Contrast paths preserve all candidates, including those not selected.
        # An ambiguous or failed candidate is a gate failure, never deletion.
        for contrast in (.5,.25,0.):
            at={};controls=[]
            for index,seed in seeds.items():
                pair=Pair(contrast_table,seed,1.)
                prefix=label+f'/a-{contrast:g}/candidate-{index}'
                if not contrast_table.pencil.has(a) and not contrast_table.wave.has(a):
                    # The LEFT pencil can be literally contrast-independent.
                    # Preserve the accepted return, not a redundant Newton solve.
                    fine=J.call(prefix+'/unchanged-source-reuse',{'sourcePencil':contrast_table.pencil,'wave':contrast_table.wave,'seed':seed,'contrast':contrast},lambda:dict(seed,nativePencil=np.asarray(contrast_table.native(contrast,seed['k'],seed['q']),complex),parameter=contrast,reuse='ACCEPTED_MODE_UNCHANGED_SOURCE'))
                    coarse=fine
                else:
                    coarse=trace_path(J,prefix+'/coarse',pair,seed,1.,contrast,.05,contrast_joint,lambda _:1.)
                    fine=trace_path(J,prefix+'/fine',pair,seed,1.,contrast,.025,contrast_joint,lambda _:1.)
                c=compare(coarse,fine);controls.append({'index':index,**c})
                require(max(c.values())<1e-8,'BOUNDARY_UNRESOLVED: contrast path refinement')
                # At fixed omega the wave relation is contrast-independent;
                # actual Q continuity and refinement retain its seed sheet.
                if contrast==0:
                    reference=np.asarray(ref_native(1.,fine['k'],fine['q']),complex)
                    join=norm(reference-fine['nativePencil'])/(1+norm(reference))
                    candidates=[m for m in uniform['backgrounds']['REFERENCE']['modes'] if m['info']['NULLITY']==pair.n]
                    scores=[abs(complex(m['info']['K'])-fine['k'])/(1+abs(fine['k']))+abs(complex(m['info']['Q'])-fine['q'])/(1+abs(fine['q'])) for m in candidates]
                    order=np.argsort(scores);match=candidates[order[0]]
                    sub=norm(projector(fine['R'])-projector(np.asarray(match['right'],complex)))
                    control={'index':index,'referenceIndex':match['info']['INDEX'],'pencilJoin':join,'rootDistance':scores[order[0]],'projectorJoin':sub,'secondDistance':scores[order[1]] if len(order)>1 else None}
                    J.report(label+f'-reference-join-{index}',control)
                    require(join<1e-10 and scores[order[0]]<1e-8 and sub<1e-8 and (len(order)==1 or scores[order[1]]>1e-6),'BOUNDARY_UNRESOLVED: independent reference correspondence')
                    contrast_correspondence[index]=match['info']['INDEX']
                at[index]=fine
            J.report(label+f'-contrast-{contrast:g}',controls);contrast_states[contrast]=at
        require(len(set(contrast_correspondence.values()))==18,'BOUNDARY_UNRESOLVED: one-to-one reference candidate census')
        target=3.;target_results={}
        for contrast,starting in contrast_states.items():
            if contrast!=1. and not P.has(a) and not wave.has(a):
                target_results[contrast]=J.call(label+f'/omega3-a-{contrast:g}-unchanged-source-reuse',{'pencil':P,'wave':wave,'contrast':contrast,'priorContrast':1.,'priorCompleteReturn':target_results[1.]},lambda:target_results[1.])
                continue
            table=compile_table(J,label+f'/frequency-table-a-{contrast:g}',P.subs(a,sp.Rational(str(contrast))),wave.subs(a,sp.Rational(str(contrast))),w,k,q)
            joint=Joint(table.wave,w,k,q,sp.sqrt(table.square))
            continued={};rows=[]
            for index,seed in starting.items():
                pair=Pair(table,seed,1.);ends=[]
                for height,step,down in ((.05,.05,.01),(.025,.025,.005)):
                    prefix=label+f'/omega-3/a-{contrast:g}/candidate-{index}/height-{height:g}'
                    state=seed;start=1.+0j
                    for segment,end,maximum in [('up',1.+1j*height,down),('across',target+1j*height,step),('down',complex(target),down)]:
                        state=trace_path(J,prefix+'/'+segment,pair,state,start,end,maximum,joint,lambda v:v);start=end
                    ends.append(state)
                check=compare(*ends);require(max(check.values())<1e-8,'BOUNDARY_UNRESOLVED: upper-half path dependence')
                state=ends[-1];pq=state['q']/physical_scales[label];kk=state['k'];tol=1e-8*(1+abs(pq))
                depth=bool(pq.imag>tol or (abs(pq.imag)<=tol and pq.real>tol))
                spatial=bool((-1 if label=='LEFT' else 1)*kk.imag>1e-8*(1+abs(kk)))
                real=bool(abs(kk.imag)<=1e-8*(1+abs(kk)))
                row={'index':index,'nullity':pair.n,'selected':index in selected,'selectedOutgoing':index in expected_out,'k':kk,'physicalQ':pq,'properDepthSheet':depth,'spatialOutwardDecay':spatial,'realNormalMomentum':real,'pathComparison':check,'residual':norm(state['residual']),'currentSign':'NOT_YET_EVALUATED'}
                rows.append(row);continued[index]=state;J.report(label+f'-omega3-a-{contrast:g}-census',rows)
                if index in selected:require(depth,'BOUNDARY_UNRESOLVED: selected state lacks physical-depth sheet')
                if index in expected_out and not real:require(spatial,'BOUNDARY_UNRESOLVED: selected state grows towards end')
                if index not in selected and depth and (spatial or real):
                    raise ValueError('BOUNDARY_UNRESOLVED: additional eligible candidate changes saved span; no current classification assumed')
            target_results[contrast]={'states':continued,'census':rows}
        results[label]=target_results
        J.blob(label+'/completed-candidate-continuations.pickle',target_results)
    return {'status':'CENTRAL_CANDIDATES_CONTINUED_CURRENT_FACE_AND_ASSEMBLY_PENDING','ends':list(results),'contrasts':[1,.5,.25,0],'frequency':3,'finiteSolves':0,'physicalRadicalScales':physical_scales,'scope':'Candidate continuation only. No current/face premise, complete end map, finite field or loss result follows from these data.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--gate',type=Path,required=True);parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    manifest=json.loads(args.manifest.read_text());gate=json.loads(args.gate.read_text())
    require(gate['status']=='READY_FOR_USER_AUTHORIZED_NUMERICAL_BOUNDARY_GATE','actual readiness gate')
    require(gate['workerSha256']==digest(__file__) and gate['manifestSha256']==digest(args.manifest),'worker/manifest gate pins')
    for path,h in gate['sourcePins'].items():require(digest(path)==h,'source/helper/authority pin '+path)
    require(gate['independentMethodClearance'] is False and gate['proceedAuthority']=='EXPLICIT_USER_NO_FURTHER_GROK_MINOR_CORRECTIONS_AND_RUN','honest review status and user exception')
    base=args.run_directory.resolve();require(base==Path(manifest['resultDirectory']),'fixed result directory');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def alarm(*_):raise TimeoutError('840-second native pilot stage limit; preserve completed returns, no retry')
    signal.signal(signal.SIGALRM,alarm);signal.alarm(840)
    started=time.monotonic();J=Journal(base);code=0;result={}
    try:
        global np,sp
        import numpy as np
        import sympy as sp
        warnings.simplefilter('error',RuntimeWarning)
        result=run(base,manifest,J)
    except BaseException as error:
        code=1;result={'status':'STOPPED_UNRESOLVED','exceptionType':type(error).__name__,'message':str(error),'traceback':traceback.format_exc(),'incompleteOperation':J.incomplete,'lastOperationReceipt':json.loads((base/'active-operation.json').read_text()) if (base/'active-operation.json').exists() else None,'finiteSolves':0}
        save_json(base/'failure.json',result)
    # Keep native alarm armed while writing results and posthashes.
    try:
        post={path:{'expected':h,'actual':digest(path)} for path,h in gate['sourcePins'].items()}
        post.update({v['path']:{'expected':v['sha256'],'actual':digest(v['path'])} for v in manifest['packets'].values()})
        intact=all(v['expected']==v['actual'] for v in post.values());save_json(base/'posthashes.json',post)
        if not intact:code=1;result['status']='INTEGRITY_FAILURE'
        result.update(completeOperations=J.count,wallSeconds=time.monotonic()-started,peakRssKiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,posthashesIntact=intact,independentMethodClearance=False,physicalLossInterpretation='RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED',automaticRetry=False)
        if J.path.exists():
            result['journal']={'path':str(J.path),'sha256':digest(J.path),'bytes':J.path.stat().st_size,'members':len(J.names)}
        text=json.dumps(numeric_json(result),indent=2,allow_nan=False)+'\n'
        (base/'checks.json').write_text(text);sys.stdout.write(text);sys.stdout.flush();signal.alarm(0)
    except BaseException:
        save_json(base/'persistence-failure.json',{'traceback':traceback.format_exc(),'resultStatus':result.get('status'),'incompleteOperation':J.incomplete});return 1
    return code


if __name__=='__main__':sys.exit(main())
