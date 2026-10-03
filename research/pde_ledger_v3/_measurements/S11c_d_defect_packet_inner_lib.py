#!/usr/bin/env python3
"""Bounded H/J/direct inner integrals; scientific calls require the guarded worker.
The old Fourier rules are restored as exact MP tuples, never constructed again.
"""
from fractions import Fraction as F
import heapq


def require(v,message):
    if v is not True:raise ValueError(message)


def point_plan():
    rows=[]
    bases=[('opposite',F(1,3),F(-1,3)),('sum-plus',F(1,3),F(5,3)),
           ('sum-minus',F(-1,3),F(-5,3)),('profile',F(1,3),F(1,3))]
    for label,k,l in bases:
        for shift in [F(0),F(-1,64),F(1,64),F(-1,128),F(1,128)]:
            rows.append({'label':label,'k':str(k),'l':str(l+shift),'shift':str(shift)})
    for axis in ['k','l']:
        for sign in [-1,1]:
            for shift in [F(-1,64),F(1,64),F(-1,128),F(1,128)]:
                r={'label':'grazing-'+axis,'k':'1/5','l':'1/5','shift':str(shift)}
                r[axis]=str(sign+shift);rows.append(r)
    rows += [{'label':'zero','k':'0','l':'0','shift':'0'},
             {'label':'control','k':'1/5','l':'2/5','shift':'0'}]
    require(len(rows)==38 and len({(r['k'],r['l']) for r in rows})==38,'fixed point census')
    return rows


def exact_cuts(k,l,T):
    """a*kappa+b with exact rational a,b. Every root/profile label is retained."""
    k,l=F(k),F(l);cuts={}
    def add(a,b,label):cuts.setdefault((F(a),F(b)),[]).append(label)
    add(0,-T,'left-radius');add(0,T,'right-radius')
    for sign in [-1,1]:
        add(sign-k,0,'height-root-'+str(sign));add(l+sign,0,'reflected-root-'+str(-sign))
    for a,label in [(F(0),'t=0'),(l-k,'t=Q')]:
        add(a,0,label)
        for v in [F(-4,5),F(-2,5),F(-1,5),F(-1,10),F(1,10),F(1,5),F(2,5),F(4,5)]:
            add(a,v,label+'-profile-offset-'+str(v))
    return [{'a':str(a),'b':str(b),'labels':labels} for (a,b),labels in sorted(cuts.items())]


def partition_interval(a,b,n):
    require(type(n)==int and n>=1 and a<b,'valid finite partition')
    # Preserve original endpoint objects. a+(b-a)*n/n may round away from b.
    points=[a]+[a+(b-a)*i/n for i in range(1,n)]+[b]
    require(all(x<y for x,y in zip(points,points[1:])),'distinct open partition')
    return points


def kernel_components(k,l,t,qi,qo,qh,qs,A1,A2,a,mu,beta,W,L,I):
    """Same arithmetic called by the new exact adapter and by both numeric routes."""
    J=a*mu**2*W*L/4*A1*A2*(k+t)*(2*k+t)*qi/(qh*(qo+beta)*(qh+beta)*(qi+beta)*(qh+qi))
    pref=W*L/(4*I)*A1*A2*(-I*mu)/((qi+beta)*(qo+beta))
    reflected=pref*k*(2*l-t)/(qs+qo)
    height=pref*k*(t+2*k)*qi/(qh*(qh+qi))
    quadratic=pref*qi**2/qh
    return [J,reflected,height,quadratic]


class InnerEvaluator:
    def __init__(self,store,oldlib,rules):
        import mpmath
        self.A=mpmath.mp.clone();self.A.dps=30
        self.B=mpmath.mp.clone();self.B.dps=50
        self.store=store;self.encode=oldlib.encode;self.cache={};self.rules={}
        for order in [24,48]:
            r=rules['A-GL'+str(order)];c=self.A
            self.rules[order]=([self.restore(c,v) for v in r['nodes']],[self.restore(c,v) for v in r['weights']])
            require(r['precision']==30,'saved A rule precision')
        r=rules['B-G7-K15'];c=self.B
        x=[self.restore(c,v) for v in r['kronrodNodes']];w=[self.restore(c,v) for v in r['kronrodWeights']]
        gx=[self.restore(c,v) for v in r['gaussNodes']];gw=[self.restore(c,v) for v in r['gaussWeights']]
        self.krule=(x,w,{i:gw[gx.index(t)] for i,t in enumerate(x) if t in gx})
        require(r['precision']==50 and len(self.krule[2])==7,'saved B rule precision and embedded membership')
        for order,(xx,ww) in self.rules.items():require(len(xx)==order==len(ww) and all(-1<t<1 for t in xx) and all(t>0 for t in ww),'saved open A rule')
        require(len(x)==len(w)==15 and all(-1<t<1 for t in x) and all(t>0 for t in w),'saved open B rule')
        self.store.put('restored-rules',{'completeOperands':rules,'constructorsCalled':False,'momentReturnsInherited':True})

    @staticmethod
    def restore(c,v):
        require(set(v)=={'mpf','decimal'} and len(v['mpf'])==4,'exact MPF encoding')
        s,m,e,b=v['mpf'];require(type(s)==int and type(e)==int and type(b)==int and isinstance(m,str),'MP tuple types')
        x=c.make_mpf((s,int(m),e,b));require(c.isfinite(x) and list(x._mpf_)==[s,int(m),e,b],'exact MP tuple restoration');return x

    def emit(self,c,name,value):return self.store.put(name,self.encode(c,value))

    @staticmethod
    def rat(c,x):
        f=F(x);return c.mpf(f.numerator)/f.denominator

    @staticmethod
    def q(c,p):
        d=c.mpf(595)/100-p*p
        return c.sqrt(d) if d>=0 else c.j*c.sqrt(-d)

    @staticmethod
    def profile(c,t):
        z=5*c.pi*t
        # The discarded sinhc tail is <=2 |z|^12/13! for |z|<=10^-6.
        # Even at 30 digits this bound is far below arithmetic roundoff.
        if abs(z)<=c.mpf('1e-6'):
            z2=z*z;s=1+z2/6+z2**2/120+z2**3/5040+z2**4/362880+z2**5/39916800
            return 1/(2*c.pi*s)
        return 5*t/(2*c.sinh(z))

    def middle(self,c,k,l,t,mutate=False):
        qi,qo=self.q(c,k),self.q(c,l);qh=self.q(c,k+t);qs=self.q(c,l-t)
        require(qi!=0 and qo!=0 and qh!=0,'unsampled depth poles')
        used_qs=qh if mutate else qs
        a=1/(1-3*c.j/10);mu=c.mpf(3)/10;beta=a*mu
        aa,bb=self.profile(c,t),self.profile(c,l-k-t)
        values=kernel_components(k,l,t,qi,qo,qh,used_qs,aa,bb,a,mu,beta,c.one,c.mpf(10),c.j)
        require(all(c.isfinite(v) for v in values),'finite middle node')
        return values,{'qi':qi,'qo':qo,'qh':qh,'qs':qs,'used_qs':used_qs,'A1':aa,'A2':bb}

    def panels(self,c,point,T,key):
        kappa=c.sqrt(595)/10;raw=exact_cuts(point['k'],point['l'],T)
        cuts=[(self.rat(c,r['a'])*kappa+self.rat(c,r['b']),r) for r in raw]
        cuts=sorted(((v,r) for v,r in cuts if -T<=v<=T),key=lambda v:v[0])
        self.emit(c,key+'/cut-operands',{'point':point,'T':T,'exactCuts':raw,'resolvedCuts':cuts})
        require(cuts[0][0]==-T and cuts[-1][0]==T,'full radius included')
        require(all(cuts[i][0]<cuts[i+1][0] for i in range(len(cuts)-1)),'distinct exact cuts may not float-alias')
        panels=[]
        for (a,ra),(b,rb) in zip(cuts,cuts[1:]):
            n=int(c.ceil(2*(b-a)));require(n>=1,'nonempty interval')
            points=partition_interval(a,b,n)
            for i in range(n):panels.append((points[i],points[i+1]))
        require(all(a<b for a,b in panels) and all(panels[i][1]==panels[i+1][0] for i in range(len(panels)-1)),'open coverage')
        return panels,{'exactCuts':raw,'resolvedCuts':cuts,'maxPhysicalWidth':'1/2'}

    def gauss(self,key,c,panels,f,order):
        nodes,weights=self.rules[order];sums=[]
        for idx,(a,b) in enumerate(panels):
            mid=(a+b)/2
            for side in [0,1]:
                origin=a if side==0 else b;length=mid-a if side==0 else b-mid
                pts=[];jac=[];vals=[];extras=[]
                for n in nodes:
                    z=(n+1)/2;x=origin+length*z*z if side==0 else origin-length*z*z
                    require(a<x<b and x!=mid,'open squared point')
                    pts.append(x);jac.append(2*length*z)
                    try:v,e=f(x)
                    except BaseException:
                        self.emit(c,key+'/failed-panel',{'bounds':[a,b],'side':side,'failedNode':x,'pointsPrefix':pts,'jacobiansPrefix':jac,'valuesPrefix':vals,'operandsPrefix':extras});raise
                    vals.append(v);extras.append(e)
                value=[c.fsum(w*j*v[col]/2 for w,j,v in zip(weights,jac,vals)) for col in range(len(vals[0]))]
                self.emit(c,key+'/panel/'+str(idx)+'/'+str(side),{'bounds':[a,b],'side':side,'points':pts,'jacobians':jac,'values':vals,'kernelOperands':extras,'return':value})
                sums.append(value)
        total=[c.fsum(v[j] for v in sums) for j in range(len(sums[0]))]
        self.emit(c,key+'/return',{'value':total,'panels':len(panels),'squareSubstitution':True,'order':order});return total

    def adaptive_panel(self,key,c,a,b,f,label):
        nodes,weights,gauss=self.krule;mid=(a+b)/2;half=(b-a)/2
        require(a<mid<b,'physical adaptive precision stagnation')
        points=[mid+half*x for x in nodes]
        require(all(a<x<b for x in points),'physical adaptive open nodes')
        pairs=[]
        for x in points:
            try:pairs.append(f(x))
            except BaseException:
                self.emit(c,key+'/failed-panel',{'bounds':[a,b],'label':label,'failedNode':x,'points':points,'completePrefix':pairs});raise
        v=[p[0] for p in pairs];extra=[p[1] for p in pairs]
        K=[half*c.fsum(w*z[j] for w,z in zip(weights,v)) for j in range(len(v[0]))]
        G=[half*c.fsum(w*v[i][j] for i,w in gauss.items()) for j in range(len(v[0]))]
        error=[abs(x-y) for x,y in zip(K,G)]
        self.emit(c,key+'/panel/'+label,{'bounds':[a,b],'points':points,'values':v,'kernelOperands':extra,'K':K,'G':G,'error':error,'disposition':'ACTIVE_LEAF_UNTIL_REPLACED','localTolerance':None})
        return {'a':a,'b':b,'K':K,'error':error}

    def adaptive(self,key,c,panels,f,budget):
        # One global empirical target per component. A singular endpoint leaf is
        # never assigned a shrinking tol/2. Its error may decrease as sqrt(h).
        require(budget>0 and len(panels)>0,'positive global empirical budget')
        leaves={};heap=[];sequence=0;count=0;refinements=0
        def push(label,record):
            nonlocal sequence
            require(all(e>=0 for e in record['error']),'nonnegative embedded errors')
            leaves[label]=record;heapq.heappush(heap,(-max(record['error']),sequence,label));sequence+=1
        for i,(a,b) in enumerate(panels):
            label='p'+str(i);push(label,self.adaptive_panel(key,c,a,b,f,label));count+=1
        components=len(next(iter(leaves.values()))['error'])
        estimated=[c.fsum(v['error'][j] for v in leaves.values()) for j in range(components)]
        self.emit(c,key+'/global-initial',{'budgetPerComponent':budget,'activeLeaves':list(leaves),'summedEmpiricalErrors':estimated,'localToleranceHalving':False})
        while True:
            # Recompute the sum from all actual leaves before accepting, and
            # periodically to avoid drift in incremental bookkeeping. This is
            # an operation-count check, not a time limit or a refinement cap.
            if refinements%64==0 or any(e<0 for e in estimated) or all(e<=budget for e in estimated):
                error=[c.fsum(v['error'][j] for v in leaves.values()) for j in range(components)]
                self.emit(c,key+'/global-sum/'+str(refinements),{'activeLeaves':list(leaves),'summedEmpiricalErrors':error,'budgetPerComponent':budget,'recomputedFromLeaves':True})
                if all(e<=budget for e in error):break
                estimated=error
            _,_,label=heapq.heappop(heap);parent=leaves[label];a,b=parent['a'],parent['b'];mid=(a+b)/2
            require(a<mid<b,'physical adaptive precision stagnation')
            left=self.adaptive_panel(key,c,a,mid,f,label+'L')
            right=self.adaptive_panel(key,c,mid,b,f,label+'R');count+=2
            require(len(left['error'])==len(right['error'])==components,'consistent vector components')
            del leaves[label];push(label+'L',left);push(label+'R',right);refinements+=1
            estimated=[c.fsum([estimated[j],-parent['error'][j],left['error'][j],right['error'][j]]) for j in range(components)]
            self.emit(c,key+'/refinement/'+str(refinements),{'replacedLeaf':label,'children':[label+'L',label+'R'],'parentError':parent['error'],'childErrors':[left['error'],right['error']],'estimatedGlobalErrors':estimated,'budgetPerComponent':budget,'selection':'largest maximum component error'})
        total=[c.fsum(v['K'][j] for v in leaves.values()) for j in range(components)]
        self.emit(c,key+'/return',{'value':total,'summedEmpiricalErrors':error,'budgetPerComponent':budget,'panelsEvaluated':count,'activeLeaves':list(leaves),'refinements':refinements,'physicalVariable':True,'errorCriterion':'global sum of actual leaf errors; empirical, not a proof'})
        require(all(e<=budget for e in error),'global empirical budget');return total

    def compare(self,key,reference,candidates):
        c=self.B;ref=[c.make_mpc(v._mpc_) if hasattr(v,'_mpc_') else c.make_mpf(v._mpf_) for v in reference]
        comparisons=[]
        for label,vals in candidates.items():
            for i,(a,b) in enumerate(zip(ref,vals)):
                value=c.make_mpc(b._mpc_) if hasattr(b,'_mpc_') else c.make_mpf(b._mpf_)
                tolerance=c.mpf('1e-9')+c.mpf('1e-7')*abs(a);delta=abs(value-a)
                comparisons.append({'route':label,'component':i,'reference':a,'other':value,'difference':delta,'tolerance':tolerance,'passed':bool(delta<=tolerance)})
        receipt=self.emit(c,key,comparisons)
        require(all(x['passed'] for x in comparisons),'unchanged numerical comparison miss; no retry')
        return receipt

    def H(self,Qcoeff):
        Qcoeff=str(F(Qcoeff))
        if Qcoeff in self.cache:
            old=self.cache[Qcoeff];self.store.put('H-reuse/'+str(len(self.cache))+'/'+str(self.reuses),{'exactQCoefficient':Qcoeff,'completedReturn':old[1],'settingsIdentical':True});self.reuses+=1;return old
        key='H/'+Qcoeff.replace('/','_');c=self.A;Q=self.rat(c,Qcoeff)*c.sqrt(595)/10;T=122
        # Exact symmetric subtraction: the j(Q) and chi contacts cancel between +/-t.
        def fa(t):
            am=self.profile(c,Q-t);ap=self.profile(c,Q+t);a=self.profile(c,t)
            return [a*5*(am-ap)/(2*c.j*t)],{'A':a,'A_Q_minus_t':am,'A_Q_plus_t':ap}
        cuts=sorted(set([c.zero,c.mpf(T),abs(Q)]+[c.mpf(i)/10 for i in range(1,9)]))
        panels=[]
        for a,b in zip(cuts,cuts[1:]):
            n=int(c.ceil(2*(b-a)));p=partition_interval(a,b,n);panels.extend(zip(p,p[1:]))
        contact=5*self.profile(c,Q)/4
        self.emit(c,key+'/input',{'QCoefficient':Qcoeff,'Q':Q,'T':T,'panels':panels,'contact':contact})
        a24=[contact+self.gauss(key+'/A24',c,panels,fa,24)[0]]
        a48=[contact+self.gauss(key+'/A48',c,panels,fa,48)[0]]
        c=self.B;Q=self.rat(c,Qcoeff)*c.sqrt(595)/10;trials=[]
        for m in range(1,4097):
            R=c.mpf(10*m);tail=10/(4*c.pi)*c.exp(-2*R/10);trials.append({'R':R,'tail':tail})
            if tail<=c.mpf('1e-14'):break
        else:raise ValueError('physical H radius capacity unavailable')
        def fb(x):
            h=(1+c.tanh(x/10))/4;j=1/(4*c.cosh(x/10)**2)
            return [h*j*c.exp(-c.j*Q*x)/(2*c.pi)],{'h':h,'j':j,'phase':c.exp(-c.j*Q*x)}
        n=int(c.ceil(4*R));p=partition_interval(-R,R,n);panels=list(zip(p,p[1:]))
        self.emit(c,key+'/B/input',{'QCoefficient':Qcoeff,'Q':Q,'radiusTrials':trials,'panels':panels,'profileScale':10,'FourierNormalization':'1/(2pi)'})
        b=self.adaptive(key+'/B',c,panels,fb,c.mpf('1e-11')/(4*38))
        analytic=[5*self.profile(c,Q)*(c.mpf(1)/4-10*c.j*Q/8)]
        compare=self.compare(key+'/comparison',a48,{'A24':a24,'Bphysical':b,'analyticProduct':analytic})
        omission=-5*self.profile(c,Q)/4
        receipt=self.emit(c,key+'/complete',{'A24':a24,'A48':[self.B.make_mpc(a48[0]._mpc_)],'B':b,'analyticProduct':analytic,'comparison':compare,'contactOmissionMovement':omission,'tTailBound':c.mpf(55)/3*c.power(2,-T)})
        self.cache[Qcoeff]=(b[0],receipt);return self.cache[Qcoeff]

    def run(self,plan):
        self.reuses=0;summaries=[]
        for index,point in enumerate(plan):
            key='point/'+str(index);results={}
            self.store.put(key+'/request',{'point':point,'routeSettings':[['A',30,24,122],['A',30,48,122],['A',30,48,124],['B',50,'G7/K15',122]],'noAutomaticRetry':True})
            for name,c,order,T in [('A24',self.A,24,122),('A48',self.A,48,122),('A48T124',self.A,48,124),('B',self.B,None,122)]:
                k=self.rat(c,point['k'])*c.sqrt(595)/10;l=self.rat(c,point['l'])*c.sqrt(595)/10
                panels,cuts=self.panels(c,point,T,key+'/'+name);self.emit(c,key+'/'+name+'/input',{'point':point,'k':k,'l':l,'T':T,'cuts':cuts,'panels':panels,'componentOrder':['J','Dreflected','Dheight','Dquadratic']})
                fun=lambda t:self.middle(c,k,l,t)
                vals=self.adaptive(key+'/'+name,c,panels,fun,c.mpf('1e-11')/(4*38)) if order is None else self.gauss(key+'/'+name,c,panels,fun,order)
                results[name]=vals+[c.fsum(vals[1:])]
            comparison=self.compare(key+'/comparison',results['A48'],{n:v for n,v in results.items() if n!='A48'})
            H,hr=self.H(F(point['l'])-F(point['k']))
            c=self.B;k=self.rat(c,point['k'])*c.sqrt(595)/10;l=self.rat(c,point['l'])*c.sqrt(595)/10;qi=self.q(c,k);qo=self.q(c,l);beta=(30+9*c.j)/109
            mixed=-3*c.j*k*qo*H/(10*(qi+beta)*(qo+beta))+results['B'][0]
            pressure=[mixed,results['B'][-1]]
            normal={'plus':[c.j*qo*v for v in pressure],'minus':[-c.j*qo*v for v in pressure]}
            # These are normalized templates per unit source, not a field/current.
            receipt=self.emit(c,key+'/complete',{'point':point,'routes':results,'comparison':comparison,'H':hr,'pressureMixedAndDirect':pressure,'normalMixedAndDirect':normal})
            summaries.append(receipt)
            if point['label']=='control':
                t=c.mpf(1)/10;base,operands=self.middle(c,k,l,t);wrong,wrongOperands=self.middle(c,k,l,t,True)
                movement=wrong[1]-base[1];normalMove=-2*normal['minus'][1]
                self.emit(c,'controls/reflected-root',{'point':point,'t':t,'baseline':base,'mutated':wrong,'operands':operands,'mutatedOperands':wrongOperands,'movement':movement})
                self.emit(c,'controls/lower-normal',{'point':point,'baseline':normal['minus'][1],'mutated':-normal['minus'][1],'movement':normalMove,'wholeDirectValue':results['B'][-1]})
                require(c.isfinite(movement) and abs(movement)>c.mpf('1e-12'),'responsive reflected-root control')
                require(c.isfinite(normalMove) and abs(normalMove)>c.mpf('1e-12'),'responsive lower-normal control')
                Q=l-k;omission=-5*self.profile(c,Q)/4
                self.emit(c,'controls/H-contact',{'point':point,'baseline':H,'mutated':H+omission,'movement':omission})
                require(c.isfinite(omission) and abs(omission)>c.mpf('1e-12'),'responsive H-contact control')
        self.store.put('complete-point-receipts',summaries)
        return {'points':len(summaries),'distinctHRequests':len(self.cache),'exactHReuses':self.reuses,'controls':3}
