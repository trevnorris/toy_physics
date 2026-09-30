#!/usr/bin/env python3
"""Numerically bind saved physical currents/faces to continued end states.

No producer, roots, frequency paths, epsilon Poly extraction or depth integral.
No computation timer. Runs only with pinned inputs under the shared guard.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import pickle
import resource
import shutil
import sqlite3
import sys
import time
import traceback
import warnings
from S11c_d_numerical_radiating_blob_store import BlobStore


def require(ok,why,evidence=None):
    if not ok:
        error=ValueError(why);error.evidence=evidence;raise error

def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def save(p,value):
    tmp=p.with_suffix(p.suffix+'.new')
    with tmp.open('w') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(p)

def readable(v):
    if isinstance(v,dict):return {str(k):readable(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [readable(x) for x in v]
    if isinstance(v,np.ndarray):return readable(v.tolist())
    if isinstance(v,np.generic):return readable(v.item())
    if isinstance(v,complex):return {'real':v.real,'imag':v.imag}
    if isinstance(v,sp.Basic):return str(v)
    return v

def norm(v):return float(np.max(np.abs(v),initial=0.))

class SourceReader(pickle.Unpickler):
    def find_class(self,module,name):
        require(module in ('builtins','collections') or module.startswith(('sympy.','numpy.')),'unsupported saved storage class '+module+'.'+name)
        return super().find_class(module,name)

def decode(raw):return SourceReader(io.BytesIO(raw)).load()

class Journal:
    def __init__(self,base):
        self.base=base;self.store=BlobStore(base/'operations.sqlite',create=True);self.count=0;self.active=None;self.refs={}
    def blob(self,name,value):return self.store.put(name,pickle.dumps(value,protocol=5))
    def call(self,name,args,fn):
        self.active=name;inp=self.blob(name+'/input.pickle',args)
        save(self.base/'active-operation.json',{'name':name,'input':inp,'status':'STARTED'})
        start=time.monotonic()
        try:result=fn()
        except BaseException as error:
            if getattr(error,'evidence',None) is not None:
                ref=self.blob(name+'/failed-check-evidence.pickle',error.evidence)
                save(self.base/'failed-check-evidence.json',{'operation':name,'message':str(error),'evidence':ref})
            raise
        out=self.blob(name+'/return.pickle',result);self.refs[name]=out
        record={'name':name,'input':inp,'return':out,'seconds':time.monotonic()-start,'status':'COMPLETE'}
        with (self.base/'operation-index.jsonl').open('a') as f:
            f.write(json.dumps(record)+'\n');f.flush();os.fsync(f.fileno())
        self.count+=1;self.active=None;save(self.base/'active-operation.json',record)
        return result
    def report(self,name,value):save(self.base/(name+'.json'),readable(value))

def leaves(v):
    if isinstance(v,sp.MatrixBase):return list(v)
    if isinstance(v,dict):return [s for x in v.values() for s in leaves(x)]
    if isinstance(v,(tuple,list)):return [s for x in v for s in leaves(x)]
    return [v]

def carrier_degrees(value,epsilon,cache):
    """Homogeneity only; independent coefficient expressions stay opaque."""
    if value in cache:return cache[value]
    if value==0:result=set()
    elif not value.has(epsilon):result={0}
    elif value==epsilon:result={1}
    elif value.is_Add:
        result=set().union(*(carrier_degrees(x,epsilon,cache) for x in value.args))
    elif value.is_Mul:
        result={0}
        for x in value.args:result={a+b for a in result for b in carrier_degrees(x,epsilon,cache)}
    elif value.is_Pow and value.exp.is_Integer and value.exp>=0:
        result={0}
        for _ in range(int(value.exp)):result={a+b for a in result for b in carrier_degrees(value.base,epsilon,cache)}
    else:raise ValueError('unsupported epsilon dependence in supplied current')
    cache[value]=result;return result


def bind_end(pair_tuple,uniform,frequency,contrast,parameters):
    packet,known=pair_tuple;r=packet['result'];c=packet['checks']
    wl,wr=r['FREQUENCY_LEGS'];kl,kr=r['NORMAL_LEGS'];ql,qr=r['BULK_LEGS']
    variables=(kl,kr,ql,qr)
    params={k:sp.Rational(v) for k,v in parameters.items()};params['omega']=sp.Integer(3)
    grades={s:sp.Rational(str(contrast))*v for s,v in frequency['origin'].items()}
    require({s.name for s in grades}=={'eta_bg','sigma_W'},'two actual saved contrast origins')
    endpoints=uniform['profileEndpoints'];bindings={wl:sp.Integer(3),wr:sp.Integer(3),**grades}
    evidence=[]
    def bind(value,epsilon=None):
        value=value.xreplace(endpoints);mapping=dict(bindings)
        for symbol in value.free_symbols-set(mapping)-set(variables):
            if symbol==epsilon:continue
            require(symbol.name in params,'unsupported source parameter '+symbol.name)
            mapping[symbol]=params[symbol.name]
        result=value.xreplace(mapping)
        require(not result.atoms(sp.Integral,sp.Limit,sp.Derivative,sp.Subs),'unresolved bound source operator')
        require(not result.free_symbols-set(variables)-({epsilon} if epsilon is not None else set()),'all material inputs bound')
        evidence.append(mapping);return result
    proof_names=('FACE_LINEAR_LIFT_RESIDUAL','FACE_PORT_LINEAR_LIFT_RESIDUAL','FACE_LEG_JOIN_RESIDUAL','FACE_BULK_RECONSTRUCTION_RESIDUAL')
    require(all(all(v==0 for v in leaves(c[name])) for name in proof_names),'saved unprojected face linearity/closure/reconstruction residuals')
    faces=[]
    require(len(c['FACE_AMPLITUDE_ROWS'])==len(c['FACE_PORT_ROWS'])==len(r['FACE_LEG_OBJECTS'])==2,'two actual faces')
    for face in range(2):
        require(len(c['FACE_AMPLITUDE_ROWS'][face])==len(c['FACE_PORT_ROWS'][face])==2,'two actual harmonic legs')
        legs=[]
        for leg in range(2):
            rows={'AMPLITUDE':c['FACE_AMPLITUDE_ROWS'][face][leg],**c['FACE_PORT_ROWS'][face][leg]}
            require(set(rows)=={'AMPLITUDE','PRESSURE','OUTWARD_VELOCITY','RELATIVE_MASS_FLUX','AFFINITY','CHEMICAL_POTENTIAL','BULK_VELOCITY','MECHANICAL_RESPONSE'},'complete supplied face row inventory')
            require(all(row.shape==(1,5) for row in rows.values()),'five physical field columns')
            legs.append({name:bind(row) for name,row in rows.items()})
        faces.append(legs)
    currents={};proofs={};epsilon=None
    for name in ('SLAB_CURRENT_MATRIX','BULK_NORMAL_CURRENT_DENSITY_MATRIX','INTERFACE_POWER_MATRIX'):
        original=r[name];eps=[s for s in original.free_symbols if s.name=='epsilon_shape']
        require(len(eps)==1 or all(x==0 for x in original),'one supplied current carrier')
        e=eps[0] if eps else epsilon
        require(e is not None,'current carrier provenance')
        if epsilon is None:epsilon=e
        require(e==epsilon,'same current carrier')
        degree=[sorted(carrier_degrees(x,e,{})) for x in original]
        require(all(v in ([],[2]) for v in degree),'source epsilon-squared homogeneity')
        currents[name]=bind(original,e);proofs[name]=degree
    pencil=bind(r['CLOSED_PENCIL_LEGS'][0])
    require(not pencil.free_symbols-{kr,qr},'right single-leg source pencil')
    # Source-specific units are retained separately, never add port power to current.
    return {'variables':variables,'epsilon':epsilon,'currents':currents,'faces':faces,'pencil':pencil,
            'carrierDegreeProofs':proofs,'savedFaceResiduals':{n:c[n] for n in proof_names},
            'sourceBranchJoins':r['SOURCE_BRANCH_JOINS'],'physicalGradeBindings':grades,
            'actualBindings':evidence,'sourceRows':{'amplitude':c['FACE_AMPLITUDE_ROWS'],'port':c['FACE_PORT_ROWS']},
            'originalFaceExpressions':r['FACE_LEG_OBJECTS'],'sourceDimensionRegistry':known}


def compile_bound(bound):
    variables=bound['variables'];eps=bound['epsilon']
    def fn(expr,args):return sp.lambdify(args,expr,'numpy',cse=True,docstring_limit=0)
    return {'current':{n:fn(v,(*variables,eps)) for n,v in bound['currents'].items()},
            'faces':[[{n:fn(v,variables) for n,v in leg.items()} for leg in face] for face in bound['faces']],
            'pencil':fn(bound['pencil'],variables[1::2])}


def evaluated_state(state,functions,physical_scale):
    k=complex(state['k']);q=complex(state['q'])/physical_scale;R=state['R']
    point=(k.conjugate(),k,q.conjugate(),q)
    current={};homogeneity={}
    for n,fn in functions['current'].items():
        current[n]=np.asarray(fn(*point,1.),complex)
        require(np.isfinite(current[n]).all(),'finite physical current/port matrix')
        homogeneity[n]={str(e):np.asarray(fn(*point,e),complex)-e**2*current[n] for e in (.5,2.)}
        require(max(norm(x) for x in homogeneity[n].values())<1e-10*(1+norm(current[n])),'current carrier rescaling')
    pencil=np.asarray(functions['pencil'](k,q),complex)
    pencil_join=(pencil-state['nativePencil'])/(1+norm(pencil))
    require(norm(pencil_join)<1e-9,'pairing source/continued original pencil join',{'boundPairingPencil':pencil,'continuedPencil':state['nativePencil'],'residual':pencil_join,'point':point})
    face_values=[[{n:np.asarray(fn(*point),complex) for n,fn in leg.items()} for leg in face] for face in functions['faces']]
    require(all(np.isfinite(v).all() for face in face_values for leg in face for v in leg.values()),'finite physical face rows')
    return {'point':point,'R':R,'k':k,'physicalQ':q,'currents':current,'homogeneityResiduals':homogeneity,
            'pencilJoin':pencil_join,'faces':face_values}


def selected_control(data,orientation,transverse):
    R=data['R'];M=data['currents']['SLAB_CURRENT_MATRIX'];G=R.conj().T@M@R
    hermitian=(G-G.conj().T)/(1+norm(G));require(norm(hermitian)<1e-9,'physical slab current Hermitian on cluster',{'gram':G,'residual':hermitian})
    result={'slabGram':G,'hermitianResidual':hermitian,'physicalQ':data['physicalQ']}
    if not transverse:return dict(result,basis=R,transform=np.eye(R.shape[1]),direction='outgoing',kind='evanescent')
    require(norm(R[3:,:])<1e-8*(1+norm(R)),'selected transverse theta/thickness components')
    projections=[];conjugacy=[]
    for face in data['faces']:
        for n in face[0]:
            projections.extend([face[0][n]@R,face[1][n]@R.conjugate()])
            conjugacy.append(face[1][n]-face[0][n].conjugate())
    face_scale=1+max(norm(v) for face in data['faces'] for leg in face for v in leg.values())*norm(R)
    require(max(norm(v) for v in projections)<1e-8*face_scale,'transverse loss-side face drives vanish',{'projections':projections,'scale':face_scale})
    require(max(norm(v) for v in conjugacy)<1e-8*face_scale,'actual harmonic face-leg conjugacy',{'residuals':conjugacy,'scale':face_scale})
    bulk=R.conj().T@data['currents']['BULK_NORMAL_CURRENT_DENSITY_MATRIX']@R
    power=R.conj().T@data['currents']['INTERFACE_POWER_MATRIX']@R
    require(max(norm(bulk),norm(power))<1e-8*(1+norm(R)**2*max(norm(v) for v in data['currents'].values())),'transverse bulk-density and interface-power projections vanish')
    eigenvalues,rotation=np.linalg.eigh((G+G.conj().T)/2)
    require(np.all(abs(eigenvalues)>1e-9*max(1.,norm(G))),'resolved transverse current')
    signed=orientation*eigenvalues;require(np.all(signed>0) or np.all(signed<0),'uniform direction in complete transverse doublet')
    transform=rotation@np.diag(1/np.sqrt(abs(eigenvalues)));basis=R@transform
    flux_residual=basis.conj().T@M@basis-np.diag(np.sign(eigenvalues))
    require(norm(flux_residual)<1e-9,'actual slab-current normalization')
    # Addressed native row/column control: omit the physical eW velocity term.
    omissions=[]
    for i,face in enumerate(data['faces']):
        row=face[0]['OUTWARD_VELOCITY'];mut=row.copy();mut[0,4]=0
        delta=(row-mut)@np.eye(5,dtype=complex)[:,4]
        omissions.append({'face':i,'column':4,'original':row,'mutated':mut,'eWProbeMovement':delta})
    require(max(norm(x['eWProbeMovement']) for x in omissions)>1e-10,'native eW velocity omission responds')
    controls=[]
    for i in range(5):
        for j in range(i+1,5):
            delta=np.zeros((5,5),complex);delta[i,j]=M[i,j];delta[j,i]=M[j,i]
            movement=basis.conj().T@delta@basis
            controls.append({'address':(i,j),'removedEntries':delta,'movement':movement})
    offdiag=max(controls,key=lambda x:norm(x['movement']))
    require(norm(offdiag['movement'])>1e-10,'source off-diagonal current omission responds',{'controls':controls,'basis':basis,'matrix':M})
    return dict(result,basis=basis,transform=transform,direction='outgoing' if np.all(signed>0) else 'incoming',kind='transverse',eigenvalues=eigenvalues,fluxNormalizationResidual=flux_residual,
                faceProjections=projections,faceConjugacyResiduals=conjugacy,bulkProjection=bulk,interfacePowerProjection=power,
                nativeVelocityOmissions=omissions,currentOffDiagonalOmission=offdiag,orientationFlipMovement=-2*orientation*(basis.conj().T@M@basis))


def assemble_end(states,census,functions,physical_scale,orientation,J,prefix):
    selected=[r for r in census if r['selected']];checked={};raw={}
    for row in selected:
        i=row['index'];state=states[i]
        data=J.call(prefix+f'/candidate-{i}/current-face-values',{'state':state,'physicalRadicalScale':physical_scale,'boundSource':J.refs[prefix+'/bound-source']},lambda s=state:evaluated_state(s,functions,physical_scale))
        real=row['realNormalMomentum'];require(not real or row['nullity']==2,'real selected branch is saved transverse doublet')
        result=J.call(prefix+f'/candidate-{i}/checks',{'data':data,'orientation':orientation,'transverse':real},lambda d=data,r=real:selected_control(d,orientation,r))
        require((result['direction']=='outgoing')==row['selectedOutgoing'],'actual current direction agrees with continued selection')
        raw[i]=data;checked[i]=result
        J.report(prefix.replace('/','-')+f'-candidate-{i}',{'index':i,'kind':result['kind'],'direction':result['direction'],'k':data['k'],'physicalQ':data['physicalQ'],'currentEigenvalues':result.get('eigenvalues'),'slabGram':result['slabGram'],'maxFaceProjection':max([0.]+[norm(x) for x in result.get('faceProjections',[])])})
    records=[]
    for row in selected:
        i=row['index'];d=checked[i]
        for j in range(d['basis'].shape[1]):records.append({'recordIndex':i,'column':j,'kind':d['kind'],'direction':d['direction'],'k':states[i]['k'],'q':states[i]['q']/physical_scale,'vector':d['basis'][:,j]})
    def cross():
        forms={name:np.zeros((len(records),len(records)),complex) for name in functions['current']}
        for a,l in enumerate(records):
            for b,r in enumerate(records):
                point=(l['k'].conjugate(),r['k'],l['q'].conjugate(),r['q'])
                for name,fn in functions['current'].items():
                    matrix=np.asarray(fn(*point,1.),complex)
                    require(np.isfinite(matrix).all(),'finite cross-current coefficients')
                    forms[name][a,b]=l['vector'].conj()@matrix@r['vector']
        transverse=[i for i,r in enumerate(records) if r['kind']=='transverse']
        residual={n:(v-v.conj().T)/(1+norm(v)) for n,v in forms.items()}
        require(max(norm(x) for x in residual.values())<1e-9,'complete selected current/power congruence Hermitian',{'forms':forms,'residuals':residual})
        mixed={n:{'rows':v[transverse,:],'columns':v[:,transverse]} for n,v in forms.items() if n!='SLAB_CURRENT_MATRIX'}
        require(max(norm(x) for v in mixed.values() for x in v.values())<1e-8*(1+max(norm(x) for x in forms.values())),'zero transverse loss-side mixed contractions',{'forms':forms,'transverseIndices':transverse,'mixed':mixed})
        return {'forms':forms,'hermitianResiduals':residual,'transverseIndices':transverse,'lossSideMixed':mixed,'anchoring':'Local end amplitudes at the finite matching position; phases belong to the later solved amplitudes. No infinite-depth integral.'}
    currents=J.call(prefix+'/full-cross-current',{'records':records,'boundSource':J.refs[prefix+'/bound-source']},cross)
    outgoing=[i for i,r in enumerate(records) if r['direction']=='outgoing'];incoming=[i for i,r in enumerate(records) if r['direction']=='incoming']
    require(len(outgoing)==5 and len(incoming)==2,'five outgoing/two incoming trace directions')
    def trace_map():
        right=np.column_stack([records[i]['vector'] for i in outgoing]);derivative=right@np.diag([1j*records[i]['k'] for i in outgoing])
        inc=np.column_stack([records[i]['vector'] for i in incoming]);inc_derivative=inc@np.diag([1j*records[i]['k'] for i in incoming])
        require(np.linalg.matrix_rank(right)==5,'full outgoing trace rank')
        matrix=np.linalg.solve(right.T,derivative.T).T;residual=(matrix@right-derivative)/(1+norm(derivative))
        condition=float(np.linalg.cond(right));require(np.isfinite(condition) and norm(residual)<1e-9,'finite trace solve and scaled residual')
        return {'outgoing':[records[i] for i in outgoing],'incoming':[records[i] for i in incoming],'right':right,'derivative':derivative,'traceMap':matrix,'residual':residual,'condition':condition,'incomingValues':inc,'incomingDerivative':inc_derivative,'incomingBoundaryData':inc_derivative-matrix@inc,'current':orientation*currents['forms']['SLAB_CURRENT_MATRIX'],'currentRecordOrder':records,'outgoingIndices':outgoing,'incomingIndices':incoming,'orientation':orientation}
    trace=J.call(prefix+'/trace-map',{'records':records,'outgoing':outgoing,'incoming':incoming,'currents':currents},trace_map)
    return {'states':states,'records':records,'checks':checked,'currents':currents,'traceMap':trace}


def run(manifest,J):
    restored={}
    for name,item in manifest['packets'].items():
        require(digest(item['path'])==item['sha256'],'source packet '+name)
        target=J.base/'inputs'/name;target.parent.mkdir(exist_ok=True);shutil.copyfile(item['path'],target)
        restored[name]=J.call('restore/'+name,item,lambda p=target:decode(p.read_bytes()))
    cp=manifest['boundary'];require(digest(cp['database'])==cp['sha256'],'completed boundary database pin')
    continued={}
    with sqlite3.connect(Path(cp['database']).as_uri()+'?mode=ro',uri=True) as db:
        for label,ref in cp['results'].items():
            def get(ref=ref):
                raw,h,size=db.execute('SELECT payload,sha256,bytes FROM blobs WHERE name=?',(ref['member'],)).fetchone()
                require(h==ref['sha256'] and size==ref['bytes'] and hashlib.sha256(raw).hexdigest()==h,'saved completed candidate return')
                return decode(raw)
            continued[label]=J.call('restore/completed-boundary-'+label,ref,get)
    physical=json.loads(Path(manifest['physicalInput']['path']).read_text());units=restored['uniformResponse']
    J.report('source-units',{'fieldUnits':units['fieldUnits'],'currentUnit':units['currentUnit'],'interpretation':'Slab transverse current only; port power never added. Physical frequency calibration remains open.'})
    results={};summary=[]
    for label in ('LEFT','RIGHT'):
        results[label]={}
        previous=None
        for contrast in (1.,.5,.25,0.):
            prefix=label+f'/a-{contrast:g}';source=restored[label+'Frequency'];pair=restored[label+'Pairing']
            bound=J.call(prefix+'/bound-source',{'pairing':pair,'profileEndpoints':restored['uniformSource']['profileEndpoints'],'frequencySource':source,'contrast':contrast,'parameters':physical['parameters']},lambda:bind_end(pair,restored['uniformSource'],source,contrast,physical['parameters']))
            target=continued[label][contrast]
            # Source-independent contrasts can share finished numerical work;
            # joins use actual expressions and serialized saved state, never text.
            signature={n:bound[n] for n in ('variables','epsilon','currents','faces','pencil','carrierDegreeProofs')}
            identical=previous is not None and signature==previous['signature'] and pickle.dumps(target,protocol=5)==previous['targetBytes']
            if identical:
                result=J.call(prefix+'/unchanged-end-map-reuse',{'boundSource':J.refs[prefix+'/bound-source'],'previousPrefix':previous['prefix'],'actualBoundExpressionsEqual':True,'actualSavedStateBytesEqual':True},lambda:previous['result'])
            else:
                functions=compile_bound(bound)
                result=assemble_end(target['states'],target['census'],functions,cp['physicalRadicalScales'][label],-1 if label=='LEFT' else 1,J,prefix)
                previous={'signature':signature,'targetBytes':pickle.dumps(target,protocol=5),'result':result,'prefix':prefix}
            results[label][contrast]=result
            summary.append({'end':label,'contrast':contrast,'status':'NUMERICAL_TRANSVERSE_CURRENT_FACE_TRACE_SUPPORTED','traceCondition':result['traceMap']['condition'],'traceResidual':norm(result['traceMap']['residual']),'incomingDirections':2,'outgoingDirections':5,'transverseDirections':4,'reusedUnchangedCompletedMap':identical})
            J.report('end-summary',summary)
    artifact=J.blob('completed-end-maps.pickle',{'ends':results,'summary':summary,'fieldUnits':units['fieldUnits'],'currentUnit':units['currentUnit']});J.report('result-artifacts',{'completedEndMaps':artifact})
    return {'status':'CENTRAL_CURRENT_FACE_TRACE_MAPS_SUPPORTED_INTEGRATION_ASSEMBLY_PENDING','frequency':3,'case':'LAB_HELD/RHO4_CONSTANT','summary':summary,'finiteSolves':0,'scope':'Finite numerical end conditions only. No field, deficit or physical leakage claim.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--run-directory',type=Path,required=True);a=p.parse_args()
    gate=json.loads(a.gate.read_text());manifest=json.loads(a.manifest.read_text())
    require(gate['status']=='READY_FOR_AUTHORIZED_NUMERICAL_END_MAPS','actual readiness')
    require(gate['workerSha256']==digest(__file__) and gate['manifestSha256']==digest(a.manifest),'worker and manifest pins')
    for path,h in gate['sourcePins'].items():require(digest(path)==h,'helper/source pin '+path)
    require(gate['wallDeadlineSeconds'] is None and gate['nativeDeadlineSeconds'] is None,'standing unlimited runtime')
    require(gate['independentMethodClearance'] is False and gate['proceedAuthority']=='EXPLICIT_USER_NO_FURTHER_GROK_MINOR_CORRECTIONS_AND_RUN','honest method status and execution authority')
    require(digest(manifest['physicalInput']['path'])==manifest['physicalInput']['sha256'],'fixed physical input pin')
    base=a.run_directory.resolve();require(str(base)==manifest['resultDirectory'],'declared output');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));start=time.monotonic();code=0;J=None
    try:
        global np,sp
        import numpy as np
        import sympy as sp
        warnings.simplefilter('error',RuntimeWarning)
        J=Journal(base);result=run(manifest,J)
    except BaseException as e:
        code=1;result={'status':'STOPPED_UNRESOLVED','exceptionType':type(e).__name__,'message':str(e),'traceback':traceback.format_exc(),'incompleteOperation':J.active if J else 'initialization','finiteSolves':0};save(base/'failure.json',result)
    post={p:{'expected':h,'actual':digest(p)} for p,h in gate['sourcePins'].items()}
    post.update({v['path']:{'expected':v['sha256'],'actual':digest(v['path'])} for v in manifest['packets'].values()})
    post[manifest['boundary']['database']]={'expected':manifest['boundary']['sha256'],'actual':digest(manifest['boundary']['database'])}
    intact=all(v['expected']==v['actual'] for v in post.values());save(base/'posthashes.json',post)
    if not intact:code=1;result['status']='INTEGRITY_FAILURE'
    if J:J.store.integrity_check();J.store.close()
    inventory={str(p.relative_to(base)):{'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(base.rglob('*')) if p.is_file()};save(base/'artifact-index.json',inventory)
    result.update(wallSeconds=time.monotonic()-start,completeOperations=J.count if J else 0,posthashesIntact=intact,peakRssKiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,independentMethodClearance=False,physicalLossInterpretation='RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED',automaticRetry=False)
    text=json.dumps(readable(result),indent=2,allow_nan=False)+'\n';(base/'checks.json').write_text(text);sys.stdout.write(text);sys.stdout.flush();return code

if __name__=='__main__':sys.exit(main())
