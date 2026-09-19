#!/usr/bin/env python3
"""Continue complete end invariant pairs in fixed seed coordinates."""
import argparse,contextlib,json,resource,shutil,signal,time
from pathlib import Path
import numpy as np
import scipy.linalg as la
import sympy as sp
import S11c_d_frequency_chart as u
f=u.f;engine=u.engine;boundary=u.boundary;grades=u.grades
PLAN=f.M/'S11c_d_frequency_end_plan.md'
PREFIX='FREQUENCY_END_LAB_HELD_RHO4_CONSTANT'


def load(base):
    data=u.load(base)
    chart,cp,path=f.accepted_packet(f.M/'S11c_d_frequency_chart_checkpoint.json','frequency-chart.pickle')
    uniform,ucp,upath=f.accepted_packet(f.M/'S11c_d_uniform_response_checkpoint.json','uniform-response.pickle')
    for packet in (chart,uniform):
        for name,h in packet['sourceFiles'].items():
            f.require(f.digest(f.ROOT/name)==h,'unchanged consumed end/chart source')
            if name in data['pins']:f.require(data['pins'][name]==h,'shared source join')
            data['pins'][name]=h
    data['operands'].update({str(path):f.digest(path),str(upath):f.digest(upath)})
    for p in (Path(__file__),PLAN,f.M/'S11c_d_frequency_chart_checkpoint.json',f.M/'S11c_d_uniform_response_checkpoint.json'):
        data['pins'][str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for name,h in data['pins'].items():
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,dest)
    data.update(chart=chart,uniform=uniform)
    data['manifest'].update(sourceFiles=data['pins'],inputPackets=data['operands'],scope='Numerical continuation of complete end invariant pairs and analytic trace/forcing/observation maps near the approved seed. No finite interior frequency matrix, inverse or pole result.')
    f.save(base/'inputs.json',data['manifest']);return data


def norm(value):return float(np.max(np.abs(value),initial=0.))


def polynomial(terms,K,Q,dK=None,dQ=None):
    n=K.shape[0];I=np.eye(n,dtype=complex);Z=np.zeros_like(I);value=Z.copy();delta=Z.copy()
    def power(A,dA,p):
        result=np.linalg.matrix_power(A,p)
        derivative=Z if dA is None or p==0 else sum((np.linalg.matrix_power(A,j)@dA@np.linalg.matrix_power(A,p-1-j) for j in range(p)),start=Z.copy())
        return result,derivative
    for (p,t),c in terms:
        kp,dk=power(K,dK,p);qt,dq=power(Q,dQ,t);value+=c*(kp@qt);delta+=c*(dk@qt+kp@dq)
    return value,delta


class Pair:
    """Full R,K,Q chart with a constant left coordinate gauge."""
    def __init__(self,table,seed):
        self.table=table;self.seed=seed;self.n=seed['R'].shape[1];self.I=np.eye(self.n,dtype=complex)
        self.gauge=np.linalg.solve(seed['R'].conj().T@seed['R'],seed['R'].conj().T)
        source=table['source'];self.w=source['frequency'];self.k=source['momentum'];self.q=source['radical']
        self.functions=[]
        for row in table['entries']:
            expressions=[c for _,c in row['numeratorTerms']+row['denominatorTerms']]
            f.require(all(not v.free_symbols-{self.w} for v in expressions),'all end coefficients retain only frequency')
            self.functions.append((row,sp.lambdify(self.w,expressions,'numpy',cse=True),sp.lambdify(self.w,[sp.diff(v,self.w) for v in expressions],'numpy',cse=True)))
        self.wave=sp.Poly(source['wave'],self.k,self.q)
        self.wave_function=sp.lambdify(self.w,[v for _,v in self.wave.terms()],'numpy')
        self.wave_derivative=sp.lambdify(self.w,[sp.diff(v,self.w) for _,v in self.wave.terms()],'numpy')
        self.native=sp.lambdify((self.w,self.k,self.q),source['livePencil'],'numpy',cse=True)
        P=np.asarray(self.native(1.,seed['k'],seed['q']),complex)
        self.row_scale=1+np.max(abs(P),axis=1)
        self.wave_scale=1+abs(seed['q'])**2+100*abs(seed['k'])**2
        self.unknown_scale=np.concatenate((np.full(5*self.n,max(1.,norm(seed['R']))),np.full(self.n**2,max(1.,abs(seed['k']))),np.full(self.n**2,max(1.,abs(seed['q'])))))
        self.coefficient_cache={}

    def pack(self,R,K,Q):return np.concatenate((R.ravel(),K.ravel(),Q.ravel()))
    def unpack(self,x):
        n=self.n;return x[:5*n].reshape(5,n),x[5*n:5*n+n*n].reshape(n,n),x[5*n+n*n:].reshape(n,n)
    def initial(self):return self.pack(self.seed['R'],self.seed['k']*self.I,self.seed['q']*self.I)
    def coefficients(self,w):
        if w not in self.coefficient_cache:
            entries=[]
            for row,fn,dfn in self.functions:
                vals=np.asarray(fn(w),complex).ravel();dvals=np.asarray(dfn(w),complex).ravel();m=len(row['numeratorTerms']);powers=[p for p,_ in row['numeratorTerms']+row['denominatorTerms']]
                entries.append((row['row'],row['column'],list(zip(powers[:m],vals[:m])),list(zip(powers[m:],vals[m:])),list(zip(powers[:m],dvals[:m])),list(zip(powers[m:],dvals[m:]))))
            wt=[p for p,_ in self.wave.terms()]
            self.coefficient_cache[w]=(entries,list(zip(wt,np.asarray(self.wave_function(w),complex).ravel())),list(zip(wt,np.asarray(self.wave_derivative(w),complex).ravel())))
        return self.coefficient_cache[w]

    def equation(self,w,x,dx=None,dw=0.):
        R,K,Q=self.unpack(x);n=self.n;zero=np.zeros_like(K)
        dR,dK,dQ=(np.zeros_like(R),zero,zero) if dx is None else self.unpack(dx)
        values=np.zeros_like(R);delta=np.zeros_like(R);minimum=np.inf
        entries,wave,dwave=self.coefficients(w)
        for i,j,nt,dt,dnt,ddt in entries:
            N,dN=polynomial(nt,K,Q,dK,dQ);D,dD=polynomial(dt,K,Q,dK,dQ)
            if dw:
                dN+=dw*polynomial(dnt,K,Q)[0];dD+=dw*polynomial(ddt,K,Q)[0]
            minimum=min(minimum,float(np.linalg.svd(D,compute_uv=False)[-1]))
            inv=np.linalg.solve(D,self.I);H=N@inv;dH=dN@inv-H@dD@inv
            values[i]+=R[j]@H;delta[i]+=dR[j]@H+R[j]@dH
        W,dW=polynomial(wave,K,Q,dK,dQ)
        if dw:dW+=dw*polynomial(dwave,K,Q)[0]
        residual=self.pack(values/self.row_scale[:,None],self.gauge@R-self.I,W/self.wave_scale)
        derivative=self.pack(delta/self.row_scale[:,None],self.gauge@dR,dW/self.wave_scale)
        return residual,derivative,minimum

    def jacobian(self,w,x):
        return np.column_stack([self.equation(w,x,np.eye(len(x),dtype=complex)[j])[1] for j in range(len(x))])

    def solve(self,w,initial):
        x=initial.copy();history=[]
        for iteration in range(15):
            residual,_,minimum=self.equation(w,x);J=self.jacobian(w,x);scaled=J*self.unknown_scale[None,:];condition=float(np.linalg.cond(scaled));history.append({'iteration':iteration,'residual':norm(residual),'scaledJacobianCondition':condition,'minimumDenominatorSingularValue':minimum})
            if norm(residual)<2e-12:break
            step=np.linalg.solve(scaled,-residual)*self.unknown_scale;accepted=False
            for damping in (1.,.5,.25,.125,.0625):
                trial=x+damping*step
                if norm(self.equation(w,trial)[0])<norm(residual):x=trial;accepted=True;break
            f.require(accepted,'end invariant-pair Newton decrease')
        f.require(norm(residual)<2e-12 and np.isfinite(condition) and minimum>1e-12,'completed regular end invariant pair')
        tangent=np.linalg.solve(J,-self.equation(w,x,dw=1.)[1]);R,K,Q=self.unpack(x)
        commutator=K@Q-Q@K;f.require(norm(commutator)<1e-9*(1+norm(K)*norm(Q)),'commuting lifted end momenta')
        return {'frequency':w,'x':x,'R':R,'K':K,'Q':Q,'tangent':tangent,'residual':residual,'commutator':commutator,'jacobian':J,'history':history,'fixedGauge':self.gauge,'unknownScale':self.unknown_scale,'rowScale':self.row_scale,'waveScale':self.wave_scale}


def seeds(base,data):
    result={};tables=data['chart']['ends']
    for label in ('LEFT','RIGHT'):
        old=data['uniform']['backgrounds'][label];channel=old['response']['channels'][label];modes={v['info']['INDEX']:v for v in old['modes']};selected=channel['outgoing']+channel['incoming'];groups={}
        for item in selected:groups.setdefault(item['RECORD_INDEX'],[]).append(item)
        rows=[]
        for index,items in groups.items():
            original=modes[index];info=original['info'];f.require(len(items)==info['NULLITY'] and sorted(v['BASIS_COLUMN'] for v in items)==list(range(info['NULLITY'])),'whole degenerate cluster retained')
            R=np.column_stack([v['vector'] for v in items]);seed={'index':index,'R':R,'k':complex(info['K']),'q':complex(info['Q']),'kind':items[0]['kind'],'direction':'incoming' if index in {v['RECORD_INDEX'] for v in channel['incoming']} else 'outgoing','items':items,'originalMode':original}
            pair=Pair(tables[label],seed);matrix=np.asarray(pair.native(1.,seed['k'],seed['q']),complex)
            join=matrix-original['pencil'];f.require(norm(join)<1e-10*(1+norm(matrix)),'actual frequency-one pencil equals accepted seed pencil')
            computed=pair.solve(1.+0j,pair.initial());f.require(norm(computed['x']-pair.initial())==0,'reuse accepted modes without a new seed solve')
            seed.update(seedState=computed,pencilReplay=join,seedKernel=matrix@R)
            f.atomic_pickle(base/(label.lower()+f'-seed-{index}.pickle'),seed);rows.append(seed)
        f.require(len(old['modes'])==18 and len(rows)==5 and sum(v['R'].shape[1] for v in rows)==7,'complete end candidate and selected direction census')
        result[label]={'clusters':rows,'allCandidates':old['modes'],'acceptedChannel':channel}
    f.atomic_pickle(base/'end-seeds.pickle',result);return result


def maps(clusters,xposition):
    outgoing=[v for v in clusters if v['seed']['direction']=='outgoing'];incoming=[v for v in clusters if v['seed']['direction']=='incoming']
    R=np.column_stack([v['state']['R'] for v in outgoing]);K=la.block_diag(*[v['state']['K'] for v in outgoing]);Ri=np.column_stack([v['state']['R'] for v in incoming]);Ki=la.block_diag(*[v['state']['K'] for v in incoming])
    inv=np.linalg.solve(R,np.eye(5));D=1j*R@K;trace=D@inv;insert=1j*Ri@Ki-trace@Ri
    # Closed amplitudes remain anchored at their finite boundaries.
    outphase=la.block_diag(*[la.expm(-1j*v['state']['K']*xposition) if v['seed']['kind']=='open' else np.eye(v['state']['K'].shape[0]) for v in outgoing]);inphase=la.expm(1j*Ki*xposition)
    observation=outphase@inv;forcing=insert@inphase;direct=-observation@Ri@inphase
    return {'right':R,'momentum':K,'incomingRight':Ri,'incomingMomentum':Ki,'derivative':D,'trace':trace,'insertion':insert,'extraction':inv,'outgoingPhase':outphase,'incomingPhase':inphase,'forcingAtCommonOrigin':forcing,'observationAtCommonOrigin':observation,'directIncomingSubtraction':direct,'inverseResidual':inv@R-np.eye(5),'traceResidual':trace@R-D,'condition':float(np.linalg.cond(R)),'scope':'Fixed seed current coordinates with analytic complex continuation; no Hermitian flux normalization at complex frequency. Closed directions stay boundary anchored.'}


def continue_pair(base,pair,seedstate,target,step):
    count=max(1,int(np.ceil(abs(target-seedstate['frequency'])/step)));state=seedstate;points=[]
    for j in range(1,count+1):
        w=seedstate['frequency']+(target-seedstate['frequency'])*j/count;initial=state['x']+(w-state['frequency'])*state['tangent'];state=pair.solve(w,initial)
        path=base/f'point-{j}.pickle';f.atomic_pickle(path,state);points.append({'path':str(path),'sha256':f.digest(path)})
    return state,points


def focused(base,data,seeddata):
    records=[]
    for label,end in seeddata.items():
        for seed in end['clusters']:
            pair=Pair(data['chart']['ends'][label],seed);x=seed['seedState']['x'];w=1.+0j
            direction=np.arange(1,len(x)+1,dtype=float)/len(x)+.3j;direction*=pair.unknown_scale
            exact=pair.equation(w,x,direction,dw=.2+.1j)[1];h=1e-6
            finite=(pair.equation(w+h*(.2+.1j),x+h*direction)[0]-pair.equation(w-h*(.2+.1j),x-h*direction)[0])/(2*h)
            derivative_residual=(exact-finite)/(1+norm(exact));f.require(norm(derivative_residual)<2e-7,'independent centered Frechet derivative')
            d=base/(label.lower()+f'-cluster-{seed["index"]}');d.mkdir();fine=d/'fine';fine.mkdir();coarse=d/'coarse';coarse.mkdir()
            a,points=continue_pair(coarse,pair,seed['seedState'],1.-.01j,.01);b,points2=continue_pair(fine,pair,seed['seedState'],1.-.01j,.005)
            difference=(a['x']-b['x'])/pair.unknown_scale;f.require(norm(difference)<1e-9,'independent step refinement for full end cluster')
            mutation=pair.equation(w,pair.pack(seed['R'],seed['k']*pair.I,-seed['q']*pair.I))[0]
            record={'label':label,'index':seed['index'],'nullity':pair.n,'seed':seed,'coarse':a,'fine':b,'derivativeResidual':derivative_residual,'stepDifference':difference,'wrongSheetResidual':mutation,'points':points+points2}
            f.atomic_pickle(base/(label.lower()+f'-focused-{seed["index"]}.pickle'),record);records.append(record)
        actual=maps([{'seed':v,'state':v['seedState']} for v in end['clusters']],-64. if label=='LEFT' else 64.);old=end['acceptedChannel']
        joins={'right':actual['right']-old['right'],'trace':actual['trace']-old['traceMap'],'insertion':actual['insertion']-old['incomingBoundaryData'],'incoming':actual['incomingRight']-old['incomingValues']}
        f.require(max(norm(v) for v in joins.values())<1e-10,'complete accepted seed boundary maps')
        f.atomic_pickle(base/(label.lower()+'-seed-maps.pickle'),{'maps':actual,'joins':joins})
    f.require(any(norm(v['wrongSheetResidual'])>1e-4 for v in records),'actual wrong-sheet control responds')
    f.atomic_pickle(base/'focused.pickle',records);return records


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);ap.add_argument('--focused',action='store_true');args=ap.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('bounded end continuation budget; retain every completed seed and point')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);start=time.monotonic();data=load(base);seeddata=seeds(base,data)
    f.require(args.focused,'production endpoint set is not yet enabled before focused acceptance')
    records=focused(base,data,seeddata)
    f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'source and input post hashes')
    checks={'status':'PASSED_FOCUSED_END_CONTINUATION','sourceFiles':data['pins'],'inputPackets':data['operands'],'clusters':len(records),'directions':sum(v['nullity'] for v in records),'candidates':sum(len(v['allCandidates']) for v in seeddata.values()),'maxDerivativeResidual':max(norm(v['derivativeResidual']) for v in records),'maxStepDifference':max(norm(v['stepDifference']) for v in records),'maxEndResidual':max(norm(v['fine']['residual']) for v in records),'maxJacobianCondition':max(v['fine']['history'][-1]['scaledJacobianCondition'] for v in records),'wrongSheetResponses':sum(norm(v['wrongSheetResidual'])>1e-4 for v in records),'artifacts':{str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts},'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':data['manifest']['scope']}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
