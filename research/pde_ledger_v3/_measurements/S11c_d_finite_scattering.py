#!/usr/bin/env python3
"""Assemble a finite modal-boundary collocation pilot from accepted operands."""
import argparse
import ast
import hashlib
import json
import operator
import pickle
import resource
import shutil
import signal
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import sympy as sp
from sympy.core.function import AppliedUndef

import S11c_d_wide_three_adaptive as prior

ROOT, STORE, M, engine = prior.ROOT, prior.STORE, prior.M, prior.engine
digest, save, atomic_pickle, unpickle = prior.digest, prior.save, prior.atomic_pickle, prior.unpickle
PLAN = M/'S11c_d_finite_scattering_plan.md'
ACCEPTANCE = ROOT/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'
CHECKPOINT = M/'S11c_d_wide_three_adaptive_checkpoint.json'
CHANNELS = M/'S11c_d_matching_channels_checkpoint.json'
NUMERICAL = M/'S11c_d_numerical_action_checkpoint.json'


def require(condition, label):
    if not condition:
        raise ValueError(label)


def accepted_packet(checkpoint, filename):
    record = json.loads(checkpoint.read_text())
    if 'status' in record:
        require(record['status']=='PUBLISHED_ANNEX_VERIFIED', ('accepted source', checkpoint))
    else:
        require(checkpoint in (CHANNELS,NUMERICAL), 'known legacy published checkpoint')
    publication=record['publication']
    require((ROOT/publication['path']).is_symlink() and
            digest(ROOT/publication['path'])==publication['sha256'], 'accepted annex publication')
    path = Path(record['runDirectory'])/filename
    require(digest(path)==record['artifacts'][filename]['sha256'], ('packet hash', path))
    return unpickle(path), record, path


def polynomial_basis(nodes, bound, size, order=0):
    order = operator.index(order)
    bound = float(bound)
    require(order >= 0 and np.isfinite(bound) and bound > 0, 'numeric polynomial derivative domain')
    coefficients = np.polynomial.chebyshev.chebder(np.eye(size), m=order, axis=0)/bound**order
    return np.polynomial.chebyshev.chebval(np.asarray(nodes)/bound, coefficients).T


def boundary_map(packet, orientation):
    outgoing, incoming, excluded = [], [], []
    for candidate in packet['CANDIDATES']:
        info, bases = candidate['INFO'], candidate['BASES']
        k = complex(info['K'])
        if not info['SHEET_MEMBERSHIP'] or not info['BULK_DECAY_DISK_CERTIFIED']:
            excluded.append({'index':info['INDEX'], 'reason':'sheet_or_bulk_decay', 'k':k})
            continue
        if info['EXACT_REAL_NORMAL']:
            channels = [v for v in packet['CHANNELS'] if v['RECORD_INDEX']==info['INDEX']]
            require(len(channels)==info['NULLITY'], 'complete real-root current directions')
            for channel in channels:
                item = dict(channel, k=k, vector=bases['FLUX_RIGHT'][:,channel['BASIS_COLUMN']],
                            kind='open', classifier=info['CLASSIFIER_STATUS'])
                require(channel['DIRECTION'] in ('INCOMING','OUTGOING'), 'resolved open direction')
                (incoming if channel['DIRECTION']=='INCOMING' else outgoing).append(item)
        elif orientation*k.imag>0:
            for j in range(info['NULLITY']):
                outgoing.append({'RECORD_INDEX':info['INDEX'], 'BASIS_COLUMN':j,
                    'k':k, 'vector':bases['RIGHT'][:,j], 'kind':'evanescent',
                    'classifier':info['CLASSIFIER_STATUS']})
        else:
            excluded.append({'index':info['INDEX'], 'reason':'outward_growth', 'k':k})
    right = np.column_stack([v['vector'] for v in outgoing])
    derivative = right@np.diag([1j*v['k'] for v in outgoing])
    require(right.shape==(5,5) and np.linalg.matrix_rank(right)==5, 'five independent outgoing trace directions')
    trace_map = np.linalg.solve(right.T, derivative.T).T
    residual = trace_map@right-derivative
    inc = np.column_stack([v['vector'] for v in incoming])
    inc_derivative = inc@np.diag([1j*v['k'] for v in incoming])
    return {'outgoing':outgoing, 'incoming':incoming, 'excluded':excluded,
            'right':right, 'derivative':derivative, 'traceMap':trace_map,
            'residual':residual, 'condition':float(np.linalg.cond(right)),
            'incomingValues':inc, 'incomingDerivative':inc_derivative,
            'incomingBoundaryData':inc_derivative-trace_map@inc,
            'current':packet['OUTWARD_CURRENT'], 'orientation':orientation}


def source_jets(record, adapter, r):
    expression = adapter.bind(record['symbolicAmplitude'])
    probes = [v for v in expression.atoms(AppliedUndef)
              if v.func.__name__.startswith('s11cdPencilProbe')]
    require(len(probes)==1 and probes[0].args==(r.zp,), 'one actual input-field source')
    probe = probes[0]
    column = int(probe.func.__name__.removeprefix('s11cdPencilProbe'))
    derivatives = expression.atoms(sp.Derivative)
    require(all(v.expr==probe and all(x==r.zp for x,_ in v.variable_count)
                for v in derivatives), 'only source-position probe derivatives')
    degree = max([0]+[sum(n for _,n in v.variable_count) for v in derivatives])
    jets = [sp.diff(probe,r.zp,n) for n in range(degree+1)]
    symbols = sp.symbols('sourceJet0:'+str(degree+1))
    formal = expression.xreplace(dict(zip(jets,symbols)))
    poly = sp.Poly(formal,*symbols)
    require(all(sum(powers)==1 for powers,_ in poly.terms()), 'source amplitude linearity')
    coefficients = [poly.coeff_monomial(v) for v in symbols]
    residual = sp.expand(formal-sum(a*v for a,v in zip(coefficients,symbols)))
    require(residual==0 and all(not a.has(AppliedUndef,sp.Derivative,sp.Subs,sp.Integral)
            and not (a.free_symbols-{r.zp}) for a in coefficients), 'complete source jet extraction')
    return {'column':column, 'coefficients':coefficients, 'residual':residual,
            'originalBoundAmplitude':expression, 'probe':probe,
            'amplitudeUnit':record['amplitudeUnit'], 'integralUnit':record['integralUnit']}


def reused_binding(directory):
    original = json.loads((directory/'inputs.json').read_text())
    name = str(Path(__file__).resolve().relative_to(ROOT))
    for source, sha in original['sourceFiles'].items():
        require(digest(directory/'source'/source)==sha, ('frozen binding source',source))
        if source != name:
            require(digest(ROOT/source)==sha, ('unchanged binding source',source))
    for path, sha in original['operandHashes'].items():
        require(digest(Path(path))==sha, ('unchanged binding operand',path))
    before = ast.parse((directory/'source'/name).read_text())
    after = ast.parse(Path(__file__).read_text())
    functions = ('accepted_packet','boundary_map','source_jets')
    for function in functions:
        old = next(n for n in before.body if isinstance(n,ast.FunctionDef) and n.name==function)
        new = next(n for n in after.body if isinstance(n,ast.FunctionDef) and n.name==function)
        require(ast.dump(old)==ast.dump(new), ('unchanged binding helper',function))
    path = directory/'source-binding.pickle'
    return unpickle(path), original, {'directory':str(directory),
        'inputsSha256':digest(directory/'inputs.json'), 'bindingSha256':digest(path),
        'originalHelperSha256':original['sourceFiles'][name], 'unchangedHelpers':list(functions)}


def load(base, resume_from=None):
    reused, original, reuse = reused_binding(resume_from) if resume_from else (None,None,None)
    packet, accepted, path = accepted_packet(CHECKPOINT,'bound-wide-three-adaptive.pickle')
    for name, sha in accepted['sourceFiles'].items():
        require(digest(ROOT/name)==sha, ('current consumed source', name))
    bound = packet['domains'][3]
    local, local_record, local_path = accepted_packet(NUMERICAL,'bound-operators.pickle')
    channels = {}
    used = {str(CHECKPOINT.relative_to(ROOT)):digest(CHECKPOINT),
            str(NUMERICAL.relative_to(ROOT)):digest(NUMERICAL),
            str(CHANNELS.relative_to(ROOT)):digest(CHANNELS)}
    operands = {str(path):digest(path),str(local_path):digest(local_path)}
    for end, sign in [('LEFT',-1),('RIGHT',1)]:
        data, _, channel_path = accepted_packet(CHANNELS,end.lower()+'.pickle')
        channels[end] = boundary_map(data,sign) if reused is None else reused['channels'][end]
        operands[str(channel_path)] = digest(channel_path)
    source_record = json.loads((M/'S11c_d_reduced_action_source_checkpoint.json').read_text())
    source_path = Path(source_record['runDirectory'])/'reduced-action.pickle'
    require(digest(source_path)==source_record['artifacts']['reduced-action.pickle']['sha256'], 'reduction source hash')
    r, dimensions = prior.domain.momentum.source.native.source.restore_context(unpickle(source_path))
    dimensions.__dict__.update(packet['dimensionState'])
    specification = json.loads((M/'S11c_d_variable_profile_development_input.json').read_text())
    adapter = engine.NumericalReducedAction(SimpleNamespace(r=r),{},specification)
    require(local['input']==specification, 'approved input join')
    if reused is None:
        jets = {i:source_jets(bound['sources'][(0,i)],adapter,r) for i in range(35)}
        checks = []
        sample = np.array([-21.,-7.,-1.,0.,2.,6.,19.])
        for (ti,si),record in sorted(bound['sources'].items()):
            jet = jets[si];width,momentum=bound['tests'][ti]
            field = sp.exp(-(r.zp/width)**2+sp.I*momentum*r.zp)
            actual = sum(a*sp.diff(field,r.zp,n) for n,a in enumerate(jet['coefficients']))
            mutated_index=next(n for n,a in enumerate(jet['coefficients']) if a!=0)
            changed=sum((sp.Rational(101,100)*a if n==mutated_index else a)*sp.diff(field,r.zp,n)
                        for n,a in enumerate(jet['coefficients']))
            a=np.broadcast_to(np.asarray(sp.lambdify(r.zp,actual,'numpy')(sample),complex),sample.shape)
            b=np.broadcast_to(np.asarray(sp.lambdify(r.zp,record['boundAmplitude'],'numpy')(sample),complex),sample.shape)
            mutated=np.broadcast_to(np.asarray(sp.lambdify(r.zp,changed,'numpy')(sample),complex),sample.shape)
            checks.append({'test':ti,'source':si,'actual':a,'accepted':b,'residual':a-b,
                           'scaled':float(np.max(abs(a-b))/(1+np.max(abs(b)))),
                           'coefficientMutationIndex':mutated_index,'coefficientMutationValues':mutated,
                           'coefficientMutationDifference':mutated-a})
        atomic_pickle(base/'source-binding.pickle',{'jets':jets,'gaussianChecks':checks,'channels':channels})
    else:
        require(original['input']==specification and original['operandHashes']==operands, 'binding input/operand joins')
        jets, checks = reused['jets'], reused['gaussianChecks']
        require(set(jets)==set(range(35)) and {(v['test'],v['source']) for v in checks}==set(bound['sources']), 'saved binding census')
        shutil.copyfile(resume_from/'source-binding.pickle',base/'source-binding.pickle')
        require(digest(base/'source-binding.pickle')==reuse['bindingSha256'], 'byte-identical saved binding')
        operands[str(resume_from/'inputs.json')]=reuse['inputsSha256']
        operands[str(resume_from/'source-binding.pickle')]=reuse['bindingSha256']
    require(max(v['scaled'] for v in checks)<1e-11, 'generic source jets reproduce accepted Gaussian bindings')
    require(all(np.max(abs(v['coefficientMutationDifference']))>0 for v in checks), 'source coefficient sensitivity')
    for row in bound['rows']:
        require(row['sourceLimit']==(r.zp,-48,48), 'accepted finite source interval')
        columns={jets[f['sourceIndex']]['column'] for f in row['factors']}
        require(len(columns)==1, 'one input field for each integral row')
    for cell in bound['cells']:
        for row_index,coefficient in cell['terms']:
            require(all(jets[f['sourceIndex']]['column']==cell['column'] for f in bound['rows'][row_index]['factors']), 'native cell/input-field join')
            require(not (coefficient.free_symbols-{r.z,r.regulator}), 'cell coefficient binding')
    pins=dict(accepted['sourceFiles'],**used)
    for p in [PLAN,ACCEPTANCE,M/'S11c_d_focused_completion_plan.md',Path(__file__).resolve()]:
        pins[str(p.relative_to(ROOT))]=digest(p)
    for name in pins:
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((ROOT/name).read_bytes())
    field_units={v['column']:dimensions.measure(v['probe']) for v in jets.values()}
    require(set(field_units)==set(range(5)), 'all field-unit slots')
    save(base/'inputs.json',{'sourceFiles':pins,'operandHashes':operands,'input':specification,'bindingReuse':reuse,
        'boundOrigin':{str(k):str(v) for k,v in adapter.input.origin.items()},
        'sourceJetCounts':{str(i):len(v['coefficients']) for i,v in jets.items()},
        'fieldUnits':{str(i):list(map(str,u)) for i,u in field_units.items()},
        'equationUnits':[list(map(str,u)) for u in bound['equationUnits']],
        'operatorBlockUnits':[[[str(a-b) for a,b in zip(eq,field_units[j])]
                               for j in range(5)] for eq in bound['equationUnits']],
        'numericalContinuumExpansionComputed':False,
        'boundaryCounts':{end:{'incoming':len(v['incoming']),'outgoing':len(v['outgoing']),
            'evanescent':sum(x['kind']=='evanescent' for x in v['outgoing'])} for end,v in channels.items()}})
    return r,bound,local['localMatrices'],jets,channels,pins,operands


class BasisMomentum(engine.BoundedSourceFourierQuadrature.ThreeMomentum):
    def prepare_basis(self,jets,nodes,bound,size,source_order):
        self.jet_data=jets; self.source_nodes,self.source_weights=engine.BoundedSourceFourierQuadrature.rule(
            [-bound,-10.,0.,10.,bound],source_order)
        max_order=max(len(v['coefficients'])-1 for v in jets.values())
        derivatives={n:polynomial_basis(self.source_nodes,bound,size,n) for n in range(max_order+1)}
        self.amplitudes={}
        for si,v in jets.items():
            matrix=np.zeros((len(self.source_nodes),size),complex)
            for n,a in enumerate(v['coefficients']):
                values=np.broadcast_to(np.asarray(sp.lambdify(self.r.zp,a,'numpy')(self.source_nodes),complex),self.source_nodes.shape)
                matrix+=values[:,None]*derivatives[n]
            self.amplitudes[si]=self.source_weights[:,None]*matrix
        self.size=size

    def source_basis(self,si,environment,cache):
        expression=self.sources[(0,si)]['frequency']
        if expression not in self.frequency_functions:
            symbols=tuple(sorted(expression.free_symbols,key=sp.default_sort_key))
            self.frequency_functions[expression]=(symbols,sp.lambdify(symbols,expression,'numpy'))
        symbols,function=self.frequency_functions[expression]
        count=len(next(iter(environment.values())))
        values=np.broadcast_to(np.asarray(function(*(environment[x] for x in symbols)),float),(count,))
        unique,inverse=np.unique(values,return_inverse=True)
        key=(si,unique.tobytes())
        if key not in cache:
            cache[key]=np.exp(-1j*unique[:,None]*self.source_nodes[None,:])@self.amplitudes[si]
        return cache[key][inverse]

    def matrix_group(self,variables,setting,pairs,width,positions,base):
        rows=[v for v in self.rows if tuple(l[0] for l in v['limits'])==variables]
        result=np.zeros((len(rows),len(positions),self.size),complex)
        direct=np.zeros((len(rows),len(positions)),complex)
        mutated_direct=np.zeros_like(direct)
        vector=np.asarray([complex(1/(j+1)**2,(-1)**j/(j+2)**2) for j in range(self.size)])
        mass=0.;nodes=0;count=0;started=time.monotonic()
        for points,weights in self.batches(variables,setting,pairs,width):
            environment={v:points[:,j] for j,v in enumerate(variables)}
            environment[self.r.regulator]=np.full(len(weights),setting['regulator'])
            profiles={};sources={}
            for index,row in enumerate(rows):
                for factor in row['factors']:
                    c=self.coefficient_value(factor['coefficient'],environment,positions,setting,profiles)
                    s=self.source_basis(factor['sourceIndex'],environment,sources)
                    result[index]+=(c*weights[:,None]).T@s
                    direct[index]+=np.sum(c*(s@vector)[:,None]*weights[:,None],axis=0)
                    mutated_weights=weights*1.001
                    mutated_direct[index]+=np.sum(c*(s@vector)[:,None]*mutated_weights[:,None],axis=0)
            mass+=float(np.sum(weights));nodes+=len(weights);count+=1
            if count%64==0:
                atomic_pickle(base/f'partial-layout-{len(variables)}.pickle',{'matrices':result,'direct':direct,'mutatedDirect':mutated_direct,
                    'nodes':nodes,'batches':count,'mass':mass,'setting':setting,'variables':variables})
        residual=result@vector-direct
        packet={'rowIndices':[v['index'] for v in rows],'matrices':result,'direct':direct,
                'actionResidual':residual,'controlVector':vector,'mutatedDirect':mutated_direct,
                'weightMutationDifference':mutated_direct-direct,
                'mass':mass,'massResidual':mass-(2*setting['momentumBound'])**len(variables),
                'nodes':nodes,'batches':count,'variables':variables,'setting':setting,
                'wallSeconds':time.monotonic()-started}
        atomic_pickle(base/f'layout-{len(variables)}.pickle',packet)
        require(np.isfinite(result).all() and np.max(abs(residual))/(1+np.max(abs(direct)))<1e-10, 'finite matrix/action contraction')
        require(not np.any((np.max(abs(direct),axis=1)>1e-12)&
                (np.max(abs(mutated_direct-direct),axis=1)==0)), 'actual momentum measure sensitivity')
        require(abs(packet['massResidual'])<1e-9*(1+abs(mass)), 'finite quadrature mass')
        return packet


def construct(base,data,size,outer,panel,source_order,profile_order):
    r,bound,local,jets,channels,pins,operands=data
    interval=48.;x=np.sort(interval*np.cos(np.arange(size)*np.pi/(size-1)))
    derivative={n:polynomial_basis(x,interval,size,n) for n in local}
    differentiation=[]
    for degree in range(min(size,8)):
        z=sp.Symbol('z');polynomial=sp.chebyshevt(degree,z/interval)
        for n,values in derivative.items():
            expected=np.broadcast_to(np.asarray(sp.lambdify(z,sp.diff(polynomial,z,n),'numpy')(x)),x.shape)
            differentiation.append(float(np.max(abs(values[:,degree]-expected))))
    require(max(differentiation)<1e-9, 'independent polynomial derivative check')
    matrix=np.zeros((5*size,5*size),complex)
    for n,coefficient_matrix in local.items():
        for i in range(5):
            for j in range(5):
                function=sp.lambdify(r.z,coefficient_matrix[i,j],'numpy',cse=True)
                values=np.broadcast_to(np.asarray(function(x),complex),x.shape)
                matrix[i*size:(i+1)*size,j*size:(j+1)*size]+=values[:,None]*derivative[n]
    local_matrix=matrix.copy()
    worker=BasisMomentum(bound['rows'],bound['sources'],r)
    worker.prepare_basis(jets,x,interval,size,source_order)
    setting={'kind':'split','momentumBound':4.,'sourceBound':interval,'profileBound':14.,
        'regulator':0.2,'outerOrder':outer,'panelOrder':panel,'innerOrders':(panel,panel),
        'sourceNodes':source_order,'profileNodes':profile_order}
    width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
    groups=[];row_matrices={}
    for variables in sorted({tuple(l[0] for l in row['limits']) for row in bound['rows']},key=len):
        group=worker.matrix_group(variables,setting,bound['pairs'],width,x,base);groups.append(group)
        row_matrices.update(zip(group['rowIndices'],group['matrices']))
    require(set(row_matrices)==set(range(80)), 'complete nonlocal row matrix census')
    for cell in bound['cells']:
        if cell['test']!=0:continue
        i,j=cell['row'],cell['column']
        for index,c in cell['terms']:
            values=np.broadcast_to(np.asarray(sp.lambdify((r.z,r.regulator),c,'numpy')(x,0.2),complex),x.shape)
            matrix[i*size:(i+1)*size,j*size:(j+1)*size]+=values[:,None]*row_matrices[index]
    original=matrix.copy();incoming_count=sum(len(c['incoming']) for c in channels.values())
    rhs=np.zeros((5*size,incoming_count),complex);offset=0
    for end,index in [('LEFT',0),('RIGHT',size-1)]:
        c=channels[end];count=len(c['incoming'])
        for i in range(5):
            row=i*size+index;matrix[row]=0
            for j in range(5):
                matrix[row,j*size:(j+1)*size]=(derivative[1][index] if i==j else 0)-c['traceMap'][i,j]*derivative[0][index]
            rhs[row,offset:offset+count]=c['incomingBoundaryData'][i]
        offset+=count
    atomic_pickle(base/'finite-system.pickle',{'matrix':matrix,'rhs':rhs,'unreplacedOperator':original,
        'localMatrix':local_matrix,'nodes':x,'derivativeMatrices':derivative,'channels':channels,
        'settings':setting,'sourceFiles':pins,'operandHashes':operands})
    row_scale=np.max(abs(matrix),axis=1);require(np.all(row_scale>0), 'nonempty equation rows')
    scaled=matrix/row_scale[:,None];column_scale=np.linalg.norm(scaled,axis=0)
    require(np.all(column_scale>0), 'nonempty unknown columns')
    balanced=scaled/column_scale[None,:]
    coefficients,_,rank,singular=np.linalg.lstsq(balanced,rhs/row_scale[:,None],rcond=None)
    coefficients/=column_scale[:,None]
    residual=matrix@coefficients-rhs;scaled_residual=residual/row_scale[:,None]
    values=np.stack([derivative[0]@coefficients[i*size:(i+1)*size] for i in range(5)])
    outgoing=[];all_modal={};offset=0
    for end,index in [('LEFT',0),('RIGHT',size-1)]:
        c=channels[end];trace=values[:,index].copy();count=len(c['incoming'])
        trace[:,offset:offset+count]-=c['incomingValues'];offset+=count
        modal=np.linalg.solve(c['right'],trace);all_modal[end]=modal
        for j,item in enumerate(c['outgoing']):
            if item['kind']=='open':outgoing.append((end,item,modal[j]))
    scattering=np.vstack([v[2] for v in outgoing])
    flux=np.zeros((len(outgoing),len(outgoing)),complex)
    for i,(end,a,_) in enumerate(outgoing):
        for j,(other,b,_) in enumerate(outgoing):
            if end==other:flux[i,j]=channels[end]['current'][a['MATRIX_COLUMN'],b['MATRIX_COLUMN']]
    outgoing_flux=np.diag(scattering.conj().T@flux@scattering).real
    incoming_flux=np.asarray([-float(np.real(c['current'][v['MATRIX_COLUMN'],v['MATRIX_COLUMN']]))
        for c in channels.values() for v in c['incoming']])
    require(np.all(incoming_flux>0), 'computed positive incident currents')
    result={'coefficients':coefficients,'fields':values,'rank':int(rank),'singularValues':singular,
        'balancedCondition':float(singular[0]/singular[-1]),'equationResidual':residual,
        'scaledEquationResidual':scaled_residual,'fullOperatorResidual':original@coefficients,
        'rowScale':row_scale,'columnScale':column_scale,'boundaryAnchoredFluxBasisScattering':scattering,
        'modalAmplitudes':all_modal,'outgoingChannelCurrent':flux,'outgoingFlux':outgoing_flux,
        'incomingFlux':incoming_flux,'outgoingFluxRatio':outgoing_flux/incoming_flux,
        'outgoingLabels':[(e,{k:v for k,v in a.items() if k!='vector'}) for e,a,_ in outgoing],
        'polynomialDerivativeResiduals':differentiation,'groups':groups}
    atomic_pickle(base/'finite-solution.pickle',result)
    require(np.isfinite(coefficients).all() and np.isfinite(scattering).all(), 'finite pilot solution')
    return result,setting


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--resume-binding',type=Path)
    parser.add_argument('--size',type=int,default=65);parser.add_argument('--outer',type=int,default=8)
    parser.add_argument('--panel',type=int,default=2);parser.add_argument('--source-order',type=int,default=128)
    parser.add_argument('--profile-order',type=int,default=128);parser.add_argument('--seconds',type=int,default=900)
    args=parser.parse_args();base=args.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('finite scattering pilot time budget; preserve saved operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(args.seconds)
    started=time.monotonic();data=load(base,args.resume_binding.resolve() if args.resume_binding else None)
    result,setting=construct(base,data,args.size,args.outer,args.panel,args.source_order,args.profile_order)
    require(all(digest(ROOT/n)==h for n,h in data[-2].items()), 'post-run sources unchanged')
    require(all(digest(Path(n))==h for n,h in data[-1].items()), 'post-run accepted operands unchanged')
    summary={'runDirectory':str(base),'unknowns':5*args.size,'incidentChannels':result['boundaryAnchoredFluxBasisScattering'].shape[1],
        'outgoingChannels':result['boundaryAnchoredFluxBasisScattering'].shape[0],
        'rank':result['rank'],'balancedCondition':result['balancedCondition'],
        'maximumEquationResidual':float(np.max(abs(result['equationResidual']))),
        'maximumScaledEquationResidual':float(np.max(abs(result['scaledEquationResidual']))),
        'outgoingFlux':result['outgoingFlux'].tolist(),'incomingFlux':result['incomingFlux'].tolist(),
        'outgoingFluxRatio':result['outgoingFluxRatio'].tolist(),'settings':setting,
        'newMomentumNodes':sum(g['nodes'] for g in result['groups']),
        'maximumMatrixActionResidual':max(float(np.max(abs(g['actionResidual']))) for g in result['groups']),
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.glob('*.pickle')},
        'sourceFiles':data[-2],'operandHashes':data[-1],
        'scope':'Finite collocation pilot with approximate modal boundaries, finite source/momentum cutoffs and regulator; continuum expansion and observable convergence pending.'}
    save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
