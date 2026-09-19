#!/usr/bin/env python3
"""Evaluate the finite reduced operator through an affine material chart."""
import argparse,contextlib,copy,json,math,resource,shutil,signal,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import scipy.linalg as la
import sympy as sp
import S11c_d_coordinate_source as source
import S11c_d_profile_form as profile

f=profile.f;engine=f.engine;response=profile.response;boundary=profile.boundary
interior=profile.interior;grades=profile.grades;J=response.J;G=response.G
PLAN=f.M/'S11c_d_coordinate_response_plan.md';PREFIX='COORDINATE_RESPONSE_LAB_HELD_RHO4_CONSTANT'


def load(base):
    data=profile.load(base);r,coefficients,ends,system,reference,baseline,packet,old,adapter,_,pins,operands=data
    coordinate,cp,path=f.accepted_packet(f.M/'S11c_d_coordinate_source_checkpoint.json','coordinate-source.pickle')
    for name,sha in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name)==sha,('accepted coordinate source',name))
        if name in pins:f.require(pins[name]==sha,('common consumed source',name))
        pins[name]=sha
    f.require(set(coordinate['records'])==set(packet['records']),'complete coordinate/native source census')
    for key,item in coordinate['records'].items():
        f.require(item['original']==packet['records'][key]['record']['ORIGINAL'] and item['address']==packet['records'][key]['address'],'actual coordinate source identity')
    operands[str(path)]=f.digest(path)
    for p in (Path(__file__),PLAN,f.M/'S11c_d_coordinate_source_checkpoint.json'):
        pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for n in pins:
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target)
    f.save(base/'inputs.json',{'sourceFiles':pins,'inputPackets':operands,'settings':system['settings'],
        'input':adapter.input.specification,'coordinatePublication':cp['publication'],
        'fieldUnits':[list(map(str,u)) for u in ends['fieldUnits']],'rowUnits':[list(map(str,u)) for u in ends['rowUnits']],
        'scope':'Actual affine transformed integrands and transported trial space, common Eulerian rows and channels inside construction.'})
    return dict(r=r,coefficients=coefficients,ends=ends,system=system,reference=reference,baseline=baseline,packet=packet,
        old=old,adapter=adapter,coordinate=coordinate,pins=pins,operands=operands)


class Chart:
    def __init__(self,data):
        self.r=data['r'];self.saved=data['coordinate'];self.g=self.saved['chart'];self.jets=self.saved['fieldJets'];self.adapter=data['adapter']
        self.d=self.g['normalJacobian'];self.A=self.g['tangentialJacobian'];self.B=self.g['fieldComponents']
        self.kappa=self.adapter.bind(sp.diff(self.g['shear'],self.g['coordinates'][2]));f.require(not self.kappa.free_symbols,'bound shear slope')
        self.rename={new:old for old,(new,value) in self.g['definitions'].items()}
        self.material=engine.NumericalReducedAction(SimpleNamespace(r=self.r),{},self.adapter.input.specification)
        self.material.input.parameters=dict(self.material.input.parameters)
        for j,kt in enumerate(self.g['tangentialCoordinates']):
            self.material.input.parameters[str(kt)]=self.g['F'][j,j]*self.adapter.bind(self.r.tangents[j])
        self.material.input.limits=dict(self.material.input.limits)
        for key,value in self.adapter.input.limits.items():
            image=source.coordinate_change(key,self.g['definitions'],self.g['forward']);self.material.input.limits[image]=value
        self.cuts=data['old']['domain-binding.pickle']['bound']['cutoffBindings']
        self.forward={n:{j:self.adapter.bind(m) for j,m in v.items()} for n,v in self.jets['forwardWithoutPhase'].items()}
        self.inverse={n:{j:self.adapter.bind(m) for j,m in v.items()} for n,v in self.jets['inverseWithoutPhase'].items()}
        self.composition={}
        for n in self.forward:
            for q in range(n+1):
                self.composition[n,q]=sum((self.forward[n][j]*self.inverse[j][q] for j in range(q,n+1)),sp.zeros(5)).applyfunc(sp.expand)
                f.require(self.composition[n,q]==(sp.eye(5) if n==q else sp.zeros(5)),'complete transported field-component/derivative composition')

    def image(self,expr):
        actual=source.coordinate_change(expr,self.g['definitions'],self.g['forward'])
        return self.bind_image(actual)

    def bind_image(self,expr):
        return engine.memo_xreplace(self.material.bind(expr),self.cuts).xreplace(self.rename)

    def physical_bound_image(self,expr):
        return self.bind_image(source.coordinate_change(expr,self.g['definitions'],self.g['forward']))

    def basis(self,positions,bound,size,order,column):
        """Contract the actual material component/derivative basis before a solve."""
        x=np.asarray(positions,float);X=x/float(self.d);phase=np.exp(1j*float(self.kappa)*X)
        derivatives={n:f.polynomial_basis(x,bound,size,n) for n in range(order+1)}
        material={}
        for n in range(order+1):
            values=np.zeros((5,len(x),size),complex)
            for q,matrix in self.inverse[n].items():
                values+=np.asarray(matrix[:,column],complex).ravel()[:,None,None]*derivatives[q][None,:,:]
            material[n]=phase[None,:,None]*values
        pulled=np.zeros((len(x),size),complex)
        for j,matrix in self.forward[order].items():
            pulled+=np.einsum('k,kns->ns',np.asarray(matrix[column,:],complex).ravel(),material[j])/phase[:,None]
        return pulled,material,derivatives[order]

    def maps(self,position):
        X=float(position)/float(self.d);U=np.exp(1j*float(self.kappa)*X)*np.asarray(self.B.inv()/self.A,complex)
        return X,U,np.linalg.inv(U)


def bind(base,data,chart):
    r=data['r'];bound=data['old']['domain-binding.pickle']['bound'];records=data['coordinate']['records'];by_address={tuple(v['address']):v for v in records.values()}
    rows=[];bindings={};limits=[]
    for old in bound['rows']:
        index=old['index'];factors=[];m=len(old['limits'])
        for fi,factor in enumerate(old['factors']):
            rec=by_address['factor',index,fi]
            f.require(rec['original']==factor['symbolicCoefficient'],'actual native factor source')
            image=chart.bind_image(rec['coordinateImage']);coefficient=image/chart.d**m
            factors.append(dict(factor,coefficient=coefficient));bindings['factor',index,fi]={'original':factor['coefficient'],'material':image,'momentumJacobian':chart.d**(-m),'unit':rec['unit']}
        mapped_limits=tuple((v,chart.d*a+chart.kappa,chart.d*b+chart.kappa) for v,a,b in old['limits'])
        mapped_source=(r.zp,old['sourceLimit'][1]/chart.d,old['sourceLimit'][2]/chart.d)
        rows.append(dict(old,factors=factors,limits=mapped_limits,sourceLimit=mapped_source));limits.append({'index':index,'original':old['limits'],'material':mapped_limits,'originalSource':old['sourceLimit'],'materialSource':mapped_source})
    sources={};jets={};pairs=[]
    for si in range(35):
        rec=by_address['source',si];native=rec['sourceJets'];f.require(native['original']==bound['sources'][(0,si)]['symbolicAmplitude'],'actual generic source')
        coefficients=[chart.bind_image(v) for v in rec['materialSourceCoefficients']]
        jets[si]=dict(data['old']['source-binding.pickle']['jets'][si],coefficients=coefficients)
        for ti in (0,1):
            old=bound['sources'][ti,si]
            frequency=chart.d*chart.image(old['symbolicFrequency'])
            sources[ti,si]=dict(old,frequency=frequency)
        pairs.append({'source':si,'column':native['column'],'coefficients':coefficients,'originalJetCount':len(native['coefficients']),'materialSourceAmplitude':rec['materialSourceAmplitude'],'unit':rec['unit']})
    support=sorted({g for v in data['packet']['records'].values() for g in v['record']['COEFFICIENTS']});coefficients={}
    for key,rec in records.items():
        address=tuple(rec['address'])
        if address[0] not in ('local','cell'):continue
        original=data['packet']['records'][key]['record']
        coefficients[address]={g:chart.image(v) for g,v in original['COEFFICIENTS'].items()}
        # The full coefficient is also bound from its actual accepted material image.
        bindings[address]={'original':engine.memo_xreplace(data['adapter'].bind(rec['original']),chart.cuts),'material':chart.bind_image(rec['coordinateImage']),'unit':rec['unit']}
    result={'rows':rows,'sources':sources,'jets':jets,'coefficients':coefficients,'bindings':bindings,'orderedLimits':limits,
        'sourceJoins':pairs,'support':support,'chartComposition':chart.composition,'geometry':chart.g,'settings':data['system']['settings']}
    f.atomic_pickle(base/'material-binding.pickle',result);return result


class MaterialMomentum(f.BasisMomentum):
    def __init__(self,binding,data,chart):
        super().__init__(binding['rows'],binding['sources'],data['r']);self.chart=chart;self.data=data
        self.original=f.BasisMomentum(data['old']['domain-binding.pickle']['bound']['rows'],data['old']['domain-binding.pickle']['bound']['sources'],data['r'])

    def prepare_basis(self,jets,nodes,bound,size,source_order):
        x,w=engine.BoundedSourceFourierQuadrature.rule([-bound,-10.,0.,10.,bound],source_order)
        self.source_nodes=x/float(self.chart.d);self.source_weights=w/float(self.chart.d);self.size=size;self.jet_data=jets;self.amplitudes={};self.basis_checks={}
        max_order=max(len(v['coefficients'])-1 for v in jets.values());cache={}
        for column in range(5):
            for n in range(max_order+1):
                actual,material,expected=self.chart.basis(x,bound,size,n,column);cache[column,n]=actual
                self.basis_checks[column,n]=float(np.max(abs(actual-expected)))
        f.require(max(self.basis_checks.values())<1e-9,'full material source derivative basis')
        for si,jet in jets.items():
            matrix=np.zeros((len(x),size),complex)
            for n,a in enumerate(jet['coefficients']):
                values=np.broadcast_to(np.asarray(sp.lambdify(self.r.zp,a,'numpy')(self.source_nodes),complex),self.source_nodes.shape)
                matrix+=values[:,None]*cache[jet['column'],n]
            self.amplitudes[si]=float(self.chart.d)*self.source_weights[:,None]*matrix
        self.original.prepare_basis(self.data['old']['source-binding.pickle']['jets'],nodes,bound,size,source_order)

    def batches(self,variables,setting,pairs,width):
        # Native panel order is retained; transform each actual point and measure.
        d=float(self.chart.d);s=float(self.chart.kappa)
        for points,weights in self.original.batches(variables,setting,pairs,width):yield d*points+s,d**len(variables)*weights

    def coefficient_value(self,expression,environment,positions,setting,profile_cache):
        material_setting=dict(setting,profileBound=setting['profileBound']/float(self.chart.d))
        return super().coefficient_value(expression,environment,np.asarray(positions)/float(self.chart.d),material_setting,profile_cache)



def local_matrices(data,binding,chart):
    system=data['system'];x=system['nodes'];size=len(x);shape=system['unreplacedOperator'].shape
    result={g:np.zeros(shape,complex) for g in binding['support']};proofs={};derivatives={}
    for n in sorted({a[1] for a in binding['coefficients'] if a[0]=='local'}):
        for j in range(5):
            actual,material,old=chart.basis(x,system['settings']['sourceBound'],size,n,j);derivatives[n,j]=actual;proofs[n,j]=float(np.max(abs(actual-old)))
    for address,series in binding['coefficients'].items():
        if address[0]!='local':continue
        _,n,i,j=address
        for g,a in series.items():
            values=np.broadcast_to(np.asarray(sp.lambdify((data['r'].z,data['r'].regulator),a,'numpy',cse=True)(x/float(chart.d),system['settings']['regulator']),complex),x.shape)
            result[g][i*size:(i+1)*size,j*size:(j+1)*size]+=values[:,None]*derivatives[n,j]
    residuals={g:interior.differences(a,data['coefficients']['matrices']['local'][g],size) for g,a in result.items()}
    f.require(max(v['maximumScaledReferenceFrame'] for v in residuals.values())<1e-10,'actual transformed local coefficient operators')
    return {'matrices':result,'derivativeResiduals':proofs,'comparisons':residuals}


def focused(base,data,binding,chart):
    r=data['r'];system=data['system'];setting=system['settings'];x=system['nodes'][[0,32,64,96,128]];bound=data['old']['domain-binding.pickle']['bound'];worker=MaterialMomentum(binding,data,chart)
    worker.prepare_basis(binding['jets'],x,setting['sourceBound'],17,setting['sourceNodes'])
    width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']));records=[];profile_count=set();sources=set()
    for variables in sorted({tuple(l[0] for l in row['limits']) for row in binding['rows']},key=len):
        ep,ew=next(worker.original.batches(variables,setting,bound['pairs'],width));ep=ep[:16];ew=ew[:16];mp=float(chart.d)*ep+float(chart.kappa);mw=float(chart.d)**len(variables)*ew
        e={v:ep[:,j] for j,v in enumerate(variables)};m={v:mp[:,j] for j,v in enumerate(variables)}
        for env in (e,m):env[r.regulator]=np.full(len(ep),setting['regulator'])
        ec,mc,es,ms={},{},{},{}
        for row in binding['rows']:
            if tuple(l[0] for l in row['limits'])!=variables:continue
            old=bound['rows'][row['index']];material=np.zeros((len(x),17),complex);original=np.zeros_like(material)
            for a,b in zip(row['factors'],old['factors']):
                cv=worker.coefficient_value(a['coefficient'],m,x,setting,mc);ov=worker.original.coefficient_value(b['coefficient'],e,x,setting,ec)
                sv=worker.source_basis(a['sourceIndex'],m,ms);tv=worker.original.source_basis(b['sourceIndex'],e,es)
                material+=(cv*mw[:,None]).T@sv;original+=(ov*ew[:,None]).T@tv;sources.add(a['sourceIndex']);profile_count.update(a['coefficient'].atoms(sp.Integral))
            residual=material-original;scale=1+boundary.norm(original)
            item={'row':row['index'],'variables':variables,'physicalNodes':ep,'materialNodes':mp,'physicalWeights':ew,'materialWeights':mw,'material':material,'original':original,'residual':residual,'scaledResidual':boundary.norm(residual)/scale,'omittedMomentumJacobianDifference':material*float(chart.d)**len(variables)-material}
            f.atomic_pickle(base/f'prefix-row-{row["index"]:02d}.pickle',item);records.append(item)
    f.require({v['row'] for v in records}==set(range(80)) and sources==set(range(35)) and len(profile_count)==6,'complete focused source/profile/integral census')
    f.require(max(v['scaledResidual'] for v in records)<1e-9,'actual mapped momentum/source/profile integrands')
    f.require(any(boundary.norm(v['omittedMomentumJacobianDifference'])>0 for v in records),'actual missing momentum Jacobian sensitivity')
    source_difference={i:worker.amplitudes[i]-worker.original.amplitudes[i] for i in range(35)}
    f.require(boundary.norm(source_difference)<1e-9,'all actual weighted generic source bases')
    local=local_matrices(data,binding,chart);f.atomic_pickle(base/'material-local.pickle',local)
    wrong_phase={}
    for n in (0,1):
        actual,_,expected=chart.basis(system['nodes'],setting['sourceBound'],17,n,0)
        wrong=actual*np.exp(2j*float(chart.kappa)*system['nodes'][:,None]/float(chart.d));wrong_phase[n]=wrong-expected
    f.require(boundary.norm(wrong_phase)>0,'wrong shear-phase sensitivity')
    result={'rows':records,'sourceBasisDifferences':source_difference,'sourceDerivativeChecks':worker.basis_checks,'sourceNodes':worker.source_nodes,'sourceWeights':worker.source_weights,'wrongPhaseDifferences':wrong_phase,'local':local}
    f.atomic_pickle(base/'focused-comparisons.pickle',result);return result


def material_ends(base,data,chart):
    """Derive material boundary and polarized-current routes before response."""
    result={};proofs={};material_ends={};all_tables={};d=float(chart.d);A=float(chart.A);kappa=float(chart.kappa)
    for end,accepted in data['ends']['ends'].items():
        X,U,T=chart.maps(accepted['orientation']*data['system']['settings']['sourceBound']);material=[];common=[]
        for original in accepted['clusters']:
            rm={g:U@v for g,v in original['R'].items()};km={g:d*v for g,v in original['K'].items()};n=len(km[(0,0)]);km[(0,0)]=km[(0,0)]+kappa*np.eye(n)
            m=dict(original,R=rm,K=km,k=d*original['k']+kappa);m['D']={g:1j*v for g,v in J.multiply(rm,km).items()};material.append(m)
            re={g:T@v for g,v in rm.items()};ke={g:v/d for g,v in km.items()};ke[(0,0)]=(km[(0,0)]-kappa*np.eye(n))/d
            de={g:T@(m['D'][g]-1j*kappa*rm[g])/d for g in G};e=dict(original,R=re,K=ke,D=de,k=(m['k']-kappa)/d);common.append(e)
        def boundary_data(clusters):
            outgoing=[v for v in clusters if v['info']['direction']=='outgoing'];incoming=[v for v in clusters if v['info']['direction']=='incoming']
            R=boundary.concatenate(outgoing,'R');D=boundary.concatenate(outgoing,'D');I=boundary.concatenate(incoming,'R');DI=boundary.concatenate(incoming,'D');inv,ir=boundary.inverse(R);trace=J.multiply(D,inv);insert=boundary.subtract(DI,J.multiply(trace,I))
            return dict(outgoing=R,outgoingDerivative=D,incoming=I,incomingDerivative=DI,outgoingInverse=inv,trace=trace,insertion=insert),outgoing,incoming,ir
        md,mout,min_,mi=boundary_data(material);ed,eout,ein,ei=boundary_data(common)
        transported_trace={g:T@v@U/d for g,v in md['trace'].items()};transported_trace[(0,0)]-=1j*kappa/d*np.eye(5)
        transported_insert={g:T@v/d for g,v in md['insertion'].items()}
        phase_in={g:la.block_diag(*(boundary.phase(v,X)[g] for v in min_))*np.exp(-1j*kappa*X) for g in G}
        phase_out={g:la.block_diag(*(boundary.phase(v,-X)[g] for v in mout if v['info']['kind']=='open'))*np.exp(1j*kappa*X) for g in G}
        material_ordered=mout+min_;offsets=np.cumsum([0]+[v['R'][(0,0)].shape[1] for v in material_ordered]);currents={};material_currents={};table_inventory={}
        all_symbols=set().union(*(v.free_symbols for tab in accepted['currentTables'].values() for v in tab.values()))
        names={v.name:v for v in all_symbols};variables=tuple(names[n] for n in ('s11cdCurrentLeftMomentum','s11cdCurrentRightMomentum','s11cdAcousticLeftNormalMomentum','s11cdAcousticRightNormalMomentum'))
        kl,kr,ql,qr=variables;ml=sp.Symbol('s11cdMaterialCurrentLeftMomentum',**kl.assumptions0);mr=sp.Symbol('s11cdMaterialCurrentRightMomentum',**kr.assumptions0)
        dimensions=engine.PHYSICAL_METADATA.dimensions;dimensions.known[ml]=dimensions.measure(kl);dimensions.known[mr]=dimensions.measure(kr)
        pair_arguments=(ml,mr,ql,qr);covectors={kl:(ml-chart.kappa)/chart.d,kr:(mr-chart.kappa)/chart.d}
        for name,table in accepted['currentTables'].items():
            transformed={key:chart.A*(chart.B.T*value.xreplace(covectors)*chart.B)/chart.d**(key[2]+key[3]) for key,value in table.items()}
            f.atomic_pickle(base/(end.lower()+'-'+name+'-material-current-tables.pickle'),{'original':table,'material':transformed,'originalArguments':variables,'materialArguments':pair_arguments,'covectors':covectors,'piolaNormal':chart.A,'tangentialMeasure':chart.A,'fieldComponents':chart.B})
            functions={key:sp.lambdify(pair_arguments,value,'numpy',cse=True) for key,value in transformed.items()}
            values={g:np.zeros((offsets[-1],offsets[-1]),complex) for g in G}
            for i,left in enumerate(material_ordered):
                for j,right in enumerate(material_ordered):
                    f.require(left['q'].imag>0 and right['q'].imag>0,'unchanged actual bulk-depth decay domain')
                    pair=boundary.current_pair(transformed,functions,left,right)
                    for g,v in pair.items():values[g][offsets[i]:offsets[i+1],offsets[j]:offsets[j+1]]=v
            material_currents[name]=values;currents[name]={g:A*v for g,v in values.items()};table_inventory[name]=transformed
        currents['total']={g:currents['slab'][g]+currents['bulk'][g] for g in G};material_currents['total']={g:material_currents['slab'][g]+material_currents['bulk'][g] for g in G}
        differences={name:boundary.subtract(ed[name],accepted[name]) for name in ed}
        differences.update(traceRoute=boundary.subtract(transported_trace,ed['trace']),insertionRoute=boundary.subtract(transported_insert,ed['insertion']),incomingPhase=boundary.subtract(phase_in,accepted['incomingOriginPhase']),outgoingPhase=boundary.subtract(phase_out,accepted['outgoingOriginPhase']),currents={name:boundary.subtract(values,accepted['currents'][name]) for name,values in currents.items()})
        f.require(boundary.norm(differences)<1e-8,'complete common-frame material boundary/phase/current route')
        omitted_measure={name:boundary.subtract(material_currents[name],currents[name]) for name in currents};f.require(boundary.norm(omitted_measure)>0,'actual tangential Fourier measure sensitivity')
        result[end]=dict(accepted,**ed,clusters=common,currents=currents,incomingOriginPhase=phase_in,outgoingOriginPhase=phase_out)
        material_ends[end]={'clusters':material,'boundary':md,'currents':material_currents,'position':X,'eulerianFromMaterial':T,'materialFromEulerian':U,'tangentialMeasure':A}
        proofs[end]={'differences':differences,'inverseResiduals':(mi,ei),'omittedTangentialMeasure':omitted_measure};all_tables[end]=table_inventory
        f.atomic_pickle(base/(end.lower()+'-material-boundary.pickle'),{'material':material_ends[end],'commonEulerian':result[end],'proofs':proofs[end]})
    packet={'ends':dict(data['ends'],ends=result),'materialEnds':material_ends,'proofs':proofs,'currentTableCount':sum(len(tab) for tables in all_tables.values() for tab in tables.values()),'sourceCurrentTableScope':'Actual polarized-current derivative tables, with full field, normal-covector chain and Piola/Fourier measure maps.'}
    f.atomic_pickle(base/'material-boundary.pickle',packet);return packet


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--focused',action='store_true');args=p.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));started=time.monotonic()
    def timeout(*_):raise TimeoutError('affine material route budget; preserve saved operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    data=load(base);chart=Chart(data);binding=bind(base,data,chart)
    f.require(args.focused,'focused binding/quad check required before production implementation')
    result=focused(base,data,binding,chart)
    end_route=material_ends(base,data,chart)
    checks={'status':'FOCUSED_MATERIAL_INTEGRANDS','runDirectory':str(base),'sourceFiles':data['pins'],'inputPackets':data['operands'],'rows':len(result['rows']),'sources':len(result['sourceBasisDifferences']),'prefixResidual':max(v['scaledResidual'] for v in result['rows']),'localResidual':max(v['maximumScaledReferenceFrame'] for v in result['local']['comparisons'].values()),'sourceBasisResidual':boundary.norm(result['sourceBasisDifferences']),'boundaryCurrentResidual':boundary.norm([v['differences'] for v in end_route['proofs'].values()]),'currentTables':end_route['currentTableCount'],'wrongPhaseResponse':boundary.norm(result['wrongPhaseDifferences']),
        'artifacts':{str(p.relative_to(base)):{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts},'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'Bounded actual material integrand/source/local operator and complete boundary/current route check; full material quadrature and scattering solve remain.'}
    f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'unchanged current sources and accepted operands')
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
