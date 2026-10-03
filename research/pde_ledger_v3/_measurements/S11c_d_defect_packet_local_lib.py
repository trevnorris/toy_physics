#!/usr/bin/env python3
"""New local Gaussian integrals only; run exclusively in the guarded worker."""
from fractions import Fraction as F
import json


def require(v,message):
    if v is not True:raise ValueError(message)


def derivative_coefficients(c,n,p,old):
    q=[c.one]
    for _ in range(n):q=old.poly_add(old.poly_derivative(q),old.poly_mul([c.j*p-c.mpf(5)/128,-c.one/64],q))
    return q


def majorant(n):
    # |p|<=5/2, with derivative of polynomial bounded coefficientwise.
    q=[F(1)]
    for _ in range(n):
        v=[F(0)]*(len(q)+1)
        for i,a in enumerate(q):
            v[i]+=F(325,128)*a;v[i+1]+=a/64
            if i:v[i-1]+=i*a
        q=v
    return q


def profile_derivative(v):
    out=[[F(0),F(0)] for _ in range(len(v)+1)]
    for i,(re,im) in enumerate(v):
        if i:
            for j,a in enumerate((re,im)):
                out[i-1][j]+=i*F(a)/10;out[i+1][j]-=i*F(a)/10
    return [[str(a),str(b)] for a,b in out]


class LocalEvaluator:
    def __init__(self,store,old,inner,rules):
        # Only restore rules and reuse generic quadrature on NEW operands. The
        # old run/H/middle/profile and all rule constructors are never called.
        self.engine=inner.InnerEvaluator(store,old,rules)
        self.A=self.engine.A;self.B=self.engine.B;self.store=store;self.old=old;self.inner=inner

    def emit(self,c,key,value):return self.store.put(key,self.old.encode(c,value))
    def rat(self,c,s):
        f=F(s);return c.mpf(f.numerator)/f.denominator
    def coeff(self,c,pairs):return [self.rat(c,a)+c.j*self.rat(c,b) for a,b in pairs]
    def p(self,c,carrier):return c.zero if carrier=='0' else c.sqrt(595)/10

    def tail(self,c,spec,R,mutant=None):
        def bound(pairs,n):
            C=sum(abs(F(a))+abs(F(b)) for a,b in pairs);q=majorant(n)
            e=c.exp(-R*R/64);mom=[32*e/R,32*e]
            for j in range(2,len(q)):mom.append(32*R**(j-1)*e+32*(j-1)*mom[j-2])
            return 2*self.rat(c,str(C))*sum(self.rat(c,str(a))*v for a,v in zip(q,mom))
        value=bound(spec['coefficientsAscending'],spec['xOrder'])
        if mutant=='Leibniz':value+=bound(profile_derivative(spec['coefficientsAscending']),0)
        return value

    def integrand(self,c,spec,carrier,route,mutant=None):
        p=self.p(c,carrier);n=spec['xOrder'];a=self.coeff(c,spec['coefficientsAscending'])
        q=derivative_coefficients(c,n,p,self.old) if route=='A' else None
        d=self.coeff(c,profile_derivative(spec['coefficientsAscending'])) if mutant=='Leibniz' else None
        require(mutant!='Leibniz' or n==1,'actual first derivative local control')
        def f(x):
            T=c.tanh(x/10);coef=self.old.poly_eval(a,T)
            if route=='A':
                Q=self.old.poly_eval(q,x)
                u=c.exp(-(x+c.mpf(5)/2)**2/128)*c.exp(c.j*p*x)
                v=c.exp(-(x-c.mpf(5)/2)**2/128)*c.exp((1 if mutant=='conjugation' else -1)*c.j*p*x)
                product=u*v
            else:
                z=c.j*p-(x+c.mpf(5)/2)/64
                Q=[c.one,z,z*z-c.one/64,z*z*z-3*z/64][n]
                product=c.exp(-x*x/64-c.mpf(25)/256)
                if mutant=='conjugation':product*=c.exp(2*c.j*p*x)
            correction=c.zero if d is None else self.old.poly_eval(d,T)
            value=product*(coef*Q+correction)
            return [value],{'x':x,'T':T,'coefficient':coef,'waveDerivative':Q,'packetProduct':product,'coefficientDerivativeCorrection':correction,'mutant':mutant}
        return f

    def request(self,spec,carrier,mutant=None):
        c=self.B;key='packet/'+carrier+'/cell/'+str(spec['cellIndex'])+('/mutant-'+mutant if mutant else '')
        require(carrier in ('0','kappa') and 0<=spec['xOrder']<=3,'fixed packet and derivative')
        if spec['zero']:
            require(all(F(a)==F(b)==0 for a,b in spec['coefficientsAscending']),'actual zero polynomial')
            rec=self.emit(c,key+'/complete',{'spec':spec,'carrier':carrier,'value':c.zero,'explicitZero':True,'quadratureCalls':0})
            return {'value':c.zero,'envelope':c.zero,'receipt':rec,'spec':spec}
        trials=[]
        for radius in range(8,4097,8):
            R=c.mpf(radius);tail=self.tail(c,spec,R,mutant);trials.append({'R':R,'tail':tail})
            if tail<=c.mpf('1e-14'):break
        else:raise ValueError('local tail radius capacity unavailable; no retry')
        self.emit(c,key+'/input',{'spec':spec,'carrier':carrier,'mutant':mutant,'radiusTrials':trials,'tailRule':'2 C sum M_j I_j; exp(-25/256) discarded upward','majorant':[[str(v) for v in majorant(n)] for n in range(4)],'R':R,'centers':['-5/2','5/2'],'width':8,'length':10,'nativeTimeTangentsAlreadyAbsorbed':True,'fullPairingUnit':'U_dual_THETA * M_ref * L_ref^-2 * T_ref^-1'})
        values={};tails={}
        for name,ctx,order,rad in [('A24',self.A,24,radius),('A48',self.A,48,radius),('A48Rplus8',self.A,48,radius+8),('B',self.B,None,radius)]:
            # Independent B partition; no A samples or arrays are shared.
            count=4*rad if name!='B' else 3*rad
            points=self.inner.partition_interval(ctx.mpf(-rad),ctx.mpf(rad),count);panels=list(zip(points,points[1:]))
            self.emit(ctx,key+'/'+name+'/input',{'spec':spec,'carrier':carrier,'mutant':mutant,'R':rad,'panels':panels,'order':order,'precision':ctx.dps})
            fun=self.integrand(ctx,spec,carrier,'A' if order else 'B',mutant)
            values[name]=(self.engine.gauss(key+'/'+name,ctx,panels,fun,order) if order else self.engine.adaptive(key+'/'+name,ctx,panels,fun,ctx.mpf('1e-13')/32))[0]
            tails[name]=self.tail(ctx,spec,ctx.mpf(rad),mutant)
        b_record=json.loads(self.store.db.execute('SELECT payload FROM records WHERE name=?',(key+'/B/return',)).fetchone()[0])
        error=self.engine.restore(c,b_record['summedEmpiricalErrors'][0]);v={k:c.make_mpc(x._mpc_) if hasattr(x,'_mpc_') else c.make_mpf(x._mpf_) for k,x in values.items()}
        ref=v['A48'];tol=c.mpf('1e-9')+c.mpf('1e-7')*abs(ref)
        comparisons=[{'route':k,'value':x,'difference':abs(x-ref),'tolerance':tol,'passed':bool(abs(x-ref)<=tol)} for k,x in v.items() if k!='A48']
        analytic=None
        if len(spec['coefficientsAscending'])==1 and mutant is None:
            p=self.p(c,carrier);a=self.coeff(c,spec['coefficientsAscending'])[0];n=spec['xOrder'];z=c.j*p-c.mpf(5)/128
            # Independent closed Gaussian moments, not numerical replay.
            moments=[c.one,z,z*z-c.one/128,z*z*z-3*z/128]
            analytic=a*8*c.sqrt(c.pi)*c.exp(-c.mpf(25)/256)*moments[n]
            comparisons.append({'route':'analyticConstant','value':analytic,'difference':abs(analytic-ref),'tolerance':tol,'passed':bool(abs(analytic-ref)<=tol)})
        envelope=max(x['difference'] for x in comparisons)+max(c.make_mpf(x._mpf_) for x in tails.values())+error
        receipt=self.emit(c,key+'/complete',{'spec':spec,'carrier':carrier,'mutant':mutant,'routes':v,'tails':tails,'BActualEmpiricalError':error,'comparisons':comparisons,'analyticConstant':analytic,'value':ref,'empiricalEnvelope':envelope,'quadratureErrorProof':False})
        require(all(r['passed'] for r in comparisons),'local comparison miss; no retry')
        return {'value':ref,'envelope':envelope,'receipt':receipt,'spec':spec}

    def run(self,specs):
        c=self.B;results={};by_packet={}
        for carrier in ('kappa','0'):
            rows=[self.request(s,carrier) for s in specs];results[carrier]=rows
            grades=[]
            for g in ([0,0],[1,0],[0,1],[1,1]):
                selected=[r for r in rows if r['spec']['grade']==g]
                grades.append({'grade':g,'value':c.fsum(r['value'] for r in selected),'empiricalEnvelopeSum':c.fsum(r['envelope'] for r in selected),'completeCellReceipts':[r['receipt'] for r in selected]})
            by_packet[carrier]=self.emit(c,'packet/'+carrier+'/local-grade-totals',{'grades':grades,'pressureIncluded':False,'fixedUnit':'U_dual_THETA * M_ref * L_ref^-2 * T_ref^-1','finiteGradeWeightsApplied':False})
        # Predeclared actual cells; no search for a favourable control.
        idx=next(i for i,s in enumerate(specs) if s['xOrder']==1 and s['grade']==[0,1])
        controls=[]
        for mutation,cell in [('Leibniz',idx),('conjugation',0)]:
            baseline=results['kappa'][cell];changed=self.request(specs[cell],'kappa',mutation)
            move=changed['value']-baseline['value'];floor=10*(baseline['envelope']+changed['envelope'])
            controls.append({'mutation':mutation,'cell':specs[cell],'baseline':baseline['receipt'],'changed':changed['receipt'],'movement':move,'tenTimesEmpiricalEnvelope':floor,'responsive':bool(abs(move)>floor and abs(move)>c.mpf('1e-12'))})
        control_receipt=self.emit(c,'controls/complete',controls)
        require(all(x['responsive'] for x in controls),'actual local numerical control silent; no retry')
        receipt=self.emit(c,'complete-local-actions',{'packets':by_packet,'controls':control_receipt,'completePacketAction':False,'pressurePending':True})
        return {'packetGradeReceipts':by_packet,'controlReceipt':control_receipt,'completeReceipt':receipt,'localActions':2,'controls':2,'currentOrLoss':None}
