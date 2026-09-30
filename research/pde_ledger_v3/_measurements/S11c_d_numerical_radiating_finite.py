"""Finite matrices, transverse-current extraction and the declared controls.

Scientific imports supplied only by the contained pilot worker. No producers.
"""
import gc
import time


def initialize(np_,sp_,integration_):
    global np,sp,I
    np,sp,I=np_,sp_,integration_


def norm(x):return float(np.max(np.abs(x),initial=0.))


def settings(ends):
    ks=[abs(complex(r['k']).real) for e in ends.values() for item in e.values() for r in item['records'] if r['kind']=='transverse']
    k=max(ks);I.require(k+.5<4,'transverse spectrum within base momentum cutoff')
    N=257
    while np.max(np.diff(np.sort(64*np.cos(np.arange(N)*np.pi/(N-1)))))>2*np.pi/k/4:N+=2
    refined=max(385,int(np.ceil(1.5*N)));refined+=int(refined%2==0)
    # A larger domain keeps the same refined N only if it passes the explicit
    # wavelength coverage gate; no unannounced extra basis adjustment.
    base={'kind':'split','momentumBound':4.,'sourceBound':64.,'profileBound':14.,'regulator':.1,'outerOrder':24,'panelOrder':8,'innerOrders':(8,8),'sourceNodes':512,'profileNodes':512,'size':N}
    fine=dict(base,momentumBound=6.,outerOrder=32,panelOrder=12,innerOrders=(12,12),sourceNodes=768,profileNodes=768,size=refined)
    domain=dict(fine,sourceBound=80.);regulator=dict(domain,regulator=.05)
    for x in (base,fine,domain,regulator):
        gap=float(np.max(np.diff(np.sort(x['sourceBound']*np.cos(np.arange(x['size'])*np.pi/(x['size']-1))))));I.require(gap<=2*np.pi/k/4,'declared numerical setting resolves transverse wavelength');x['maximumNodeGap']=gap
    return {'base':base,'refinement':fine,'domain':domain,'regulator':regulator},{'maximumTransverseK':k,'minimumWavelength':2*np.pi/k,'baseSize':N,'refinedSize':refined}


def matrix_group(worker,variables,setting,positions,J,prefix):
    rows=[v for v in worker.rows if tuple(l[0] for l in v['limits'])==variables]
    result=np.zeros((len(rows),len(positions),worker.size),complex);direct=np.zeros((len(rows),len(positions)),complex)
    vector=np.asarray([complex(1/(j+1)**2,(-1)**j/(j+2)**2) for j in range(worker.size)])
    width=float(worker.binding['abel']['width'].subs(worker.r.regulator,setting['regulator']));mass=0.;count=0;nodes=0;started=time.monotonic()
    for points,weights in worker.batches(variables,setting,worker.binding['pairs'],width):
        env={v:points[:,j] for j,v in enumerate(variables)};env[worker.r.regulator]=np.full(len(weights),setting['regulator']);profiles={};sources={}
        for index,row in enumerate(rows):
            for factor in row['factors']:
                coefficient=worker.coefficient_value(factor['coefficient'],env,positions,setting,profiles);source=worker.source_basis(factor['sourceIndex'],env,sources)
                result[index]+=(coefficient*weights[:,None]).T@source
                direct[index]+=np.sum(coefficient*(source@vector)[:,None]*weights[:,None],axis=0)
        mass+=float(weights.sum());nodes+=len(weights);count+=1
        if count%4096==0:
            receipt=J.blob(prefix+f'/partial-{count:09d}.pickle',{'matrices':result,'direct':direct,'mass':mass,'nodes':nodes,'batches':count,'lastPoints':points,'lastWeights':weights,'setting':setting,'variables':variables,'rowIndices':[r['index'] for r in rows]})
            J.report('last-matrix-progress',{'group':prefix,'batches':count,'nodes':nodes,'receipt':receipt,'seconds':time.monotonic()-started})
    residual=result@vector-direct;mass_residual=mass-(2*setting['momentumBound'])**len(variables)
    record={'rowIndices':[r['index'] for r in rows],'matrices':result,'direct':direct,'controlVector':vector,'actionResidual':residual,'mass':mass,'massResidual':mass_residual,'nodes':nodes,'batches':count,'variables':variables,'setting':setting,'seconds':time.monotonic()-started}
    I.require(np.isfinite(result).all() and norm(residual)<1e-10*(1+norm(direct)) and abs(mass_residual)<1e-9*(1+mass),'finite native matrix/action and measure checks',record)
    receipts={row['index']:J.blob(prefix+f'/matrix-row-{row["index"]}.pickle',result[i]) for i,row in enumerate(rows)}
    del record['matrices'];record['matrixReceipts']=receipts
    return record


def matrices(binding,row_matrices,x,setting,r,H):
    N=setting['size'];L=setting['sourceBound'];derivative={n:H.polynomial_basis(x,L,N,n) for n in set(binding['local'])|{0,1}}
    matrix=np.zeros((5*N,5*N),complex)
    for n,coefficients in binding['local'].items():
        for i in range(5):
            for j in range(5):
                c=np.broadcast_to(np.asarray(sp.lambdify(r.z,coefficients[i,j],'numpy',cse=True)(x),complex),x.shape)
                matrix[i*N:(i+1)*N,j*N:(j+1)*N]+=c[:,None]*derivative[n]
    local=matrix.copy();native=matrix.copy()
    def values(c):return np.broadcast_to(np.asarray(sp.lambdify((r.z,r.regulator),c,'numpy',cse=True)(x,setting['regulator']),complex),x.shape)
    for term in binding['termJoins']:
        i,j=term['row'],term['column'];c=binding['actual']['cell',i,j,term['term']]
        matrix[i*N:(i+1)*N,j*N:(j+1)*N]+=values(c)[:,None]*row_matrices[term['integralIndex']]
    for cell in binding['cells']:
        i,j=cell['row'],cell['column']
        for index,c in cell['terms']:native[i*N:(i+1)*N,j*N:(j+1)*N]+=values(c)[:,None]*row_matrices[index]
    residual=matrix-native;I.require(norm(residual)<1e-11*(1+norm(matrix)),'independent native occurrence/row-matrix accumulation',{'residual':residual})
    return {'matrix':matrix,'localMatrix':local,'derivatives':derivative,'occurrenceJoinResidual':residual,'positions':x,'setting':setting,'units':{'fields':binding['fieldUnits'],'rows':binding['rowUnits']}}


def solve(system,ends,J,prefix):
    import scipy.linalg as la
    N=system['setting']['size'];d=system['derivatives'];matrix=system['matrix'].copy();rhs=np.zeros((5*N,4),complex)
    for end,index,offset in (('LEFT',0,0),('RIGHT',N-1,2)):
        c=ends[end]['traceMap']
        for i in range(5):
            row=i*N+index;matrix[row]=0
            for j in range(5):matrix[row,j*N:(j+1)*N]=(d[1][index] if i==j else 0)-c['traceMap'][i,j]*d[0][index]
            rhs[row,offset:offset+2]=c['incomingBoundaryData'][i]
    J.blob(prefix+'/boundary-system.pickle',{'matrix':matrix,'rhs':rhs,'endTraceMaps':{e:v['traceMap'] for e,v in ends.items()},'setting':system['setting']})
    row_scale=np.maximum(np.max(abs(matrix),axis=1),1e-30);balanced=matrix/row_scale[:,None];column_scale=np.maximum(np.linalg.norm(balanced,axis=0),1e-30);balanced/=column_scale[None,:]
    balanced_rhs=rhs/row_scale[:,None]
    lu=la.solve(balanced,balanced_rhs,assume_a='gen',check_finite=True);svd,_,rank,singular=la.lstsq(balanced,balanced_rhs,lapack_driver='gelsd')
    difference=norm(lu-svd)/(1+norm(svd));coefficients=lu/column_scale[:,None];residual=(matrix@coefficients-rhs)/row_scale[:,None]
    checks={'rank':int(rank),'dimension':5*N,'condition':float(singular[0]/singular[-1]),'independentSolveDifference':difference,'maximumScaledResidual':norm(residual),'singularValues':singular}
    J.blob(prefix+'/linear-solve-checks.pickle',checks)
    I.require(rank==5*N and difference<1e-8 and norm(residual)<1e-9 and np.isfinite(coefficients).all(),'finite full-rank response and independent solve',checks)
    fields=np.stack([d[0]@coefficients[i*N:(i+1)*N] for i in range(5)]);incident=np.zeros(4);output=np.zeros((2,4));modal={};cross={}
    for e,(end,index,offset) in enumerate((('LEFT',0,0),('RIGHT',N-1,2))):
        c=ends[end]['traceMap'];trace=fields[:,index].copy();trace[:,offset:offset+2]-=c['incomingValues'];out=la.solve(c['right'],trace)
        full=np.zeros((7,4),complex);full[c['outgoingIndices'],:]=out
        for j,q in enumerate(c['incomingIndices']):full[q,offset+j]=1;incident[offset+j]=-c['current'][q,q].real
        transverse=[i for i in c['outgoingIndices'] if c['currentRecordOrder'][i]['kind']=='transverse']
        rest=[i for i in c['outgoingIndices'] if i not in transverse];tin=c['incomingIndices'];allT=transverse+tin
        output[e]=np.real(np.einsum('ic,ij,jc->c',full[transverse].conj(),c['current'][np.ix_(transverse,transverse)],full[transverse]))
        def mixed(a,b):return 2*np.real(np.einsum('ic,ij,jc->c',full[a].conj(),c['current'][np.ix_(a,b)],full[b]))
        cross[end]={'transverseOther':mixed(allT,rest),'incomingOutgoingTransverse':mixed(tin,transverse)};modal[end]={'outgoing':out,'all':full}
    I.require(np.all(incident>0),'actual positive incident current')
    deficit=1-output.sum(axis=0)/incident
    normalized_cross={e:{name:v/incident for name,v in x.items()} for e,x in cross.items()}
    return {'coefficients':coefficients,'fields':fields,'modalAmplitudes':modal,'incident':incident,'outwardTransverseByEnd':output,'reflected':np.concatenate((output[0,:2],output[1,2:])),'transmitted':np.concatenate((output[1,:2],output[0,2:])),'deficit':deficit,'normalizedCrossCurrents':normalized_cross,'solveChecks':checks,'scaledResidual':residual,'setting':system['setting'],'physicalInterpretation':'RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED'}


def control_report(cases):
    sensitivities={};envelope=np.full(4,1e-6)
    for a in (0.,1.):
        for left,right in (('base','refinement'),('refinement','domain'),('domain','regulator')):
            delta=abs(cases[right,a]['deficit']-cases[left,a]['deficit']);sensitivities[str((left,right,a))]=delta;envelope=np.maximum(envelope,delta)
    result=[]
    for a in (0.,.25,.5,1.):
        case=cases['base',a];floor=abs(cases['base',0.]['deficit']);compare=np.maximum(envelope,floor);D=case['deficit']
        cross=np.max(np.stack([abs(v) for end in case['normalizedCrossCurrents'].values() for v in end.values()]),axis=0)
        partition=cross<=np.maximum(compare,.1*abs(D));status=[]
        for i in range(4):
            if not partition[i]:s='PARTITION_UNRESOLVED'
            elif D[i]<-compare[i]:s='NEGATIVE_DEFICIT_UNRESOLVED'
            elif D[i]>3*compare[i] and envelope[i]<=.2*D[i]:s='FINITE_MODEL_DEFICIT_NUMERICALLY_RESOLVED'
            else:s='NO_DEFICIT_RESOLVED'
            status.append(s)
        result.append({'contrast':a,'deficit':D,'matchedUniformDeficit':cases['base',0.]['deficit'],'contrastDiagnostic':D-cases['base',0.]['deficit'],'numericalEnvelope':envelope,'uniformFloor':floor,'comparisonEnvelope':compare,'crossMaximum':cross,'statuses':status,'negativeBeyondResolution':D < -compare})
    ratios=[]
    for i in range(4):
        positive=all(row['statuses'][i]=='FINITE_MODEL_DEFICIT_NUMERICALLY_RESOLVED' for row in result[1:])
        ratios.append({'incident':i,'status':'RESOLVED_CONTRAST_RATIOS' if positive else 'SCALING_UNRESOLVED','ratios':[float(result[3]['deficit'][i]/result[2]['deficit'][i]),float(result[2]['deficit'][i]/result[1]['deficit'][i])] if positive else None})
    return {'cases':result,'sensitivities':sensitivities,'scaling':ratios,'physicalInterpretation':'RETAINED_ORDER_LOSS_INTERPRETATION_UNRESOLVED','analogLightCalibration':'OPEN'}
