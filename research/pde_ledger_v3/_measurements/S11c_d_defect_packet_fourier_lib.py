#!/usr/bin/env python3
"""Two source-bound Fourier routes. Importing this module performs no science.
All numerical construction/evaluation must be called only by a guarded worker.
"""
import hashlib
import json
import math
import sqlite3
from fractions import Fraction


def require(value, message):
    if value is not True:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


class DurableStore:
    """Immutable full operands/returns, per-record FULL synchronous transactions."""
    def __init__(self,path):
        require(not path.exists(),'refuse existing journal')
        self.db=sqlite3.connect(path)
        self.db.execute('PRAGMA journal_mode=DELETE')
        self.db.execute('PRAGMA synchronous=FULL')
        self.db.execute('CREATE TABLE records (name TEXT PRIMARY KEY, sha256 TEXT NOT NULL, payload BLOB NOT NULL)')
        self.db.commit()

    def put(self,name,value):
        raw=json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
        h=hashlib.sha256(raw).hexdigest()
        with self.db:
            self.db.execute('INSERT INTO records VALUES (?,?,?)',(name,h,raw))
        return {'record':name,'sha256':h,'bytes':len(raw)}

    def close(self):
        self.db.close()


def encode(ctx,value):
    if hasattr(value,'_mpc_'):
        return {'mpc':[encode(ctx,value.real),encode(ctx,value.imag)]}
    if hasattr(value,'_mpf_'):
        sign,mantissa,exponent,bits=value._mpf_
        require(ctx.isfinite(value),'nonfinite numerical evidence')
        return {'mpf':[sign,str(mantissa),exponent,bits],'decimal':ctx.nstr(value,ctx.dps+8)}
    if isinstance(value,dict):return {str(k):encode(ctx,v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(ctx,v) for v in value]
    return value


def poly_eval(c,x):
    out=0
    for a in reversed(c):out=out*x+a
    return out


def poly_add(a,b):
    return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))]


def poly_mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def poly_derivative(a):return [i*a[i] for i in range(1,len(a))] or [0]


def polynomial_request_key(spec):
    # Complete coefficient values and all derivative/normalization/phase operands.
    return digest(spec)


def stieltjes_coefficients():
    """Exact E8 from integral P7 E8 x^j = 0, j=1,3,5,7 on [-1,1]."""
    p=[Fraction(0),Fraction(-35,16),Fraction(0),Fraction(315,16),Fraction(0),Fraction(-693,16),Fraction(0),Fraction(429,16)]
    def moment(j):return sum(a*Fraction(2,i+j+1) for i,a in enumerate(p) if (i+j)%2==0)
    powers=(0,2,4,6)
    rows=[[moment(j+i) for i in powers]+[-moment(j+8)] for j in (1,3,5,7)]
    initial=[row[:] for row in rows]
    for col in range(4):
        piv=next(i for i in range(col,4) if rows[i][col])
        rows[col],rows[piv]=rows[piv],rows[col]
        v=rows[col][col];rows[col]=[a/v for a in rows[col]]
        for i in range(4):
            if i!=col:
                v=rows[i][col];rows[i]=[a-v*b for a,b in zip(rows[i],rows[col])]
    e=[Fraction(0)]*9;e[8]=Fraction(1)
    for i,power in enumerate(powers):e[power]=rows[i][-1]
    require(all(sum(row[i]*e[powers[i]] for i in range(4))==row[-1] for row in initial),'exact Stieltjes moments')
    return e,initial


def kronrod15(ctx,store):
    ec,system=stieltjes_coefficients();e=[ctx.mpf(v.numerator)/v.denominator for v in ec]
    gx,gw=ctx.gauss_quadrature(7,'legendre');gx=list(gx);gw=list(gw)
    brackets=[-ctx.one]+gx+[ctx.one];extra=[];tolerance=ctx.mpf('1e-48')
    for a,b in zip(brackets,brackets[1:]):
        fa=poly_eval(e,a);fb=poly_eval(e,b)
        require(fa*fb<0,'Kronrod interlacing bracket')
        while b-a>tolerance:
            c=(a+b)/2;fc=poly_eval(e,c)
            require(c!=a and c!=b,'Kronrod precision stagnation')
            if fc==0:a=b=c;break
            if fa*fc<0:b=c
            else:a=c;fa=fc
        extra.append((a+b)/2)
    nodes=sorted(gx+extra)
    matrix=ctx.matrix([[x**n for x in nodes] for n in range(15)])
    weights=list(ctx.lu_solve(matrix,ctx.matrix([ctx.mpf(2)/(n+1) if n%2==0 else 0 for n in range(15)])))
    residuals=[sum(w*x**n for x,w in zip(nodes,weights))-(ctx.mpf(2)/(n+1) if n%2==0 else 0) for n in range(24)]
    record={'gaussNodes':gx,'gaussWeights':gw,'kronrodNodes':nodes,'kronrodWeights':weights,'momentResidualsThrough23':residuals,
        'StieltjesCoefficients':[str(v) for v in ec],'StieltjesSystem':[[str(v) for v in row] for row in system],'precision':ctx.dps}
    receipt=store.put('rules/B-G7-K15',encode(ctx,record))
    require(all(-1<x<1 for x in nodes) and all(w>0 for w in weights),'positive open Kronrod rule')
    require(max(map(abs,residuals))<ctx.mpf('1e-42'),'Kronrod rule moments')
    gauss={i:gw[gx.index(x)] for i,x in enumerate(nodes) if x in gx}
    require(len(gauss)==7,'embedded Gauss nodes')
    return nodes,weights,gauss,receipt


def panel_grid(ctx,radius,center,width):
    cuts=sorted(set([-radius,ctx.zero,radius]+([-center] if -radius<-center<radius else [])))
    panels=[]
    for a,b in zip(cuts,cuts[1:]):
        n=int(ctx.ceil((b-a)/width))
        points=[a+(b-a)*i/n for i in range(n+1)]
        panels.extend(zip(points,points[1:]))
    require(panels[0][0]==-radius and panels[-1][1]==radius and all(a<b and b-a<=width*(1+ctx.eps*8) for a,b in panels),'panel coverage and width')
    require(all(panels[i][1]==panels[i+1][0] for i in range(len(panels)-1)),'no panel gaps')
    return panels


class FourierEvaluator:
    """No shared samples, nodes, values, quadrature cache or precision context across routes."""
    def __init__(self,store):
        import mpmath
        self.A=mpmath.mp.clone();self.A.dps=30
        self.B=mpmath.mp.clone();self.B.dps=50
        self.store=store;self.cache={};self.rules={}
        for order in (24,48):
            x,w=self.A.gauss_quadrature(order,'legendre');x=list(x);w=list(w)
            residuals=[sum(a*b**n for b,a in zip(x,w))-(self.A.mpf(2)/(n+1) if n%2==0 else 0) for n in range(2*order)]
            receipt=store.put('rules/A-GL'+str(order),encode(self.A,{'nodes':x,'weights':w,'momentResiduals':residuals,'precision':30}))
            require(all(-1<v<1 for v in x) and all(v>0 for v in w) and max(map(abs,residuals))<self.A.mpf('1e-25'),'Gauss rule certification')
            self.rules[order]=(x,w,receipt)
        self.krule=kronrod15(self.B,store)

    def scalar(self,ctx,text):
        # Only fixed rational numbers and the one saved radical carrier are accepted.
        if text=='sqrt(595)/10':return ctx.sqrt(595)/10
        if text=='-sqrt(595)/10':return -ctx.sqrt(595)/10
        a=Fraction(text);return ctx.mpf(a.numerator)/a.denominator

    def product(self,ctx,spec,argument,shifted):
        coefficients=[ctx.mpc(self.scalar(ctx,v[0]),self.scalar(ctx,v[1])) for v in spec['coefficientsAscending']]
        center=self.scalar(ctx,spec['center']);carrier=self.scalar(ctx,spec['carrier']);s=ctx.mpf(8)
        role=spec['role'];sign=1 if role=='X' else -1
        nu=argument-carrier
        c=ctx.mpf(5) if shifted else ctx.zero
        imaginary=-c*ctx.sign(nu)
        n=spec['spatialOrders'][0] if role=='X' else 0
        nt,n2,n3=spec['timeOrder'],spec['spatialOrders'][1],spec['spatialOrders'][2]
        native=(-3*ctx.j)**nt*(ctx.j/5)**n2*(ctx.j/10)**n3 if role=='X' else ctx.one
        normalization=1/(2*ctx.pi) if role=='X' else ctx.one
        # Q_n is in w=z-center. Derivatives act on u BEFORE multiplying the source coefficient.
        q=[ctx.one]
        for _ in range(n):q=poly_add(poly_derivative(q),poly_mul([ctx.j*carrier,-1/s**2],q))
        # Expand actual Q(y+i*imaginary) and the requested argument derivative powers.
        shifted_q=[ctx.zero]*len(q)
        for j,a in enumerate(q):
            for r in range(j+1):shifted_q[r]+=a*math.comb(j,r)*(ctx.j*imaginary)**(j-r)
        degree=spec['argumentDerivative']
        # d^m/dk^m for X: (-iz)^m; d^m/dl^m for Y(-l): (+iz)^m.
        deriv=[ctx.one]
        for _ in range(degree):deriv=poly_mul(deriv,[-sign*ctx.j*(center+ctx.j*imaginary),-sign*ctx.j])
        all_poly=poly_mul(shifted_q,deriv)
        majorant=[abs(v.real)+abs(v.imag) for v in all_poly]
        norm=sum(abs(v.real)+abs(v.imag) for v in coefficients)
        C0=abs(native)*normalization*norm
        def f(y):
            w=y+ctx.j*imaginary;z=center+w
            phase=ctx.exp(-w*w/(2*s*s)-ctx.j*nu*z)
            return normalization*native*poly_eval(coefficients,ctx.tanh(z/10))*poly_eval(all_poly,y)*phase
        data={'argument':argument,'nu':nu,'center':center,'carrier':carrier,'shift':imaginary,'spatialDerivative':n,
            'native':native,'normalization':normalization,'Q':q,'shiftedQ':shifted_q,'argumentDerivativePolynomial':deriv,
            'fullPolynomial':all_poly,'majorant':majorant,'coefficientL1':norm,'C0':C0,
            'tailExponentialFactor':ctx.exp(c*c/(2*s*s)-c*abs(nu)),'s':s}
        return f,data

    def tail_plan(self,ctx,data):
        s=data['s'];poly=data['majorant'];pref=data['C0']*data['tailExponentialFactor'];trials=[]
        for m in range(1,4097):
            R=m*s;ex=ctx.exp(-R*R/(2*s*s))
            moments=[s*ctx.sqrt(ctx.pi/2)*ctx.erfc(R/(ctx.sqrt(2)*s)),s*s*ex]
            for n in range(2,len(poly)):moments.append(s*s*R**(n-1)*ex+(n-1)*s*s*moments[n-2])
            tail=2*pref*sum(a*b for a,b in zip(poly,moments))
            trials.append({'m':m,'R':R,'oneSidedMoments':moments[:len(poly)],'twoTailBound':tail})
            if tail<=ctx.mpf('1e-14'):
                return R,tail,trials
        raise ValueError('finite transform radius capacity unavailable; no quadrature attempted')

    def route_a(self,key,spec,arg,order):
        c=self.A;f,data=self.product(c,spec,arg,True);R,tail,trials=self.tail_plan(c,data)
        width=min(c.mpf(10)/8,c.mpf(8)/8,c.pi/(4*(1+abs(data['nu']))));panels=panel_grid(c,R,data['center'],width)
        nodes,weights,rule=self.rules[order]
        self.store.put(key+'/A'+str(order)+'/input',encode(c,{'spec':spec,'data':data,'radiusTrials':trials,'panels':panels,'maxWidth':width,'rule':rule}))
        sums=[]
        for i,(a,b) in enumerate(panels):
            mid=(a+b)/2;half=(b-a)/2;points=[mid+half*x for x in nodes];values=[f(x) for x in points]
            value=half*c.fsum(w*v for w,v in zip(weights,values))
            self.store.put(key+'/A'+str(order)+'/panel/'+str(i),encode(c,{'bounds':[a,b],'points':points,'values':values,'value':value,'phaseAdvance':abs(data['nu'])*(b-a)}));sums.append(value)
        total=c.fsum(sums);result={'value':total,'tail':tail,'panels':len(panels),'precision':c.dps,'order':order}
        self.store.put(key+'/A'+str(order)+'/return',encode(c,result));return result

    def route_b(self,key,spec,arg):
        c=self.B;f,data=self.product(c,spec,arg,False);R,tail,trials=self.tail_plan(c,data)
        width=min(c.mpf(10)/5,c.mpf(8)/5,c.pi/(3*(1+abs(data['nu']))));panels=panel_grid(c,R,data['center'],width)
        nodes,weights,gauss,rule=self.krule
        self.store.put(key+'/B/input',encode(c,{'spec':spec,'data':data,'radiusTrials':trials,'initialPanels':panels,'maxWidth':width,'rule':rule,'errorTarget':'1e-13'}))
        stack=[(a,b,c.mpf('1e-13')/len(panels),str(i)) for i,(a,b) in reversed(list(enumerate(panels)))];sums=[];errors=[]
        while stack:
            a,b,target,identifier=stack.pop();mid=(a+b)/2;half=(b-a)/2
            require(a<mid<b,'adaptive precision stagnation')
            points=[mid+half*x for x in nodes];values=[f(x) for x in points]
            fine=half*c.fsum(w*v for w,v in zip(weights,values));coarse=half*c.fsum(w*values[i] for i,w in gauss.items())
            # Embedded difference is an empirical estimate, never a proved error bound.
            err=abs(fine-coarse);accepted=err<=target
            self.store.put(key+'/B/panel/'+identifier,encode(c,{'bounds':[a,b],'points':points,'values':values,'fine':fine,'coarse':coarse,'estimate':err,'allocatedTarget':target,'accepted':accepted}))
            if accepted:sums.append(fine);errors.append(err)
            else:stack.extend([(mid,b,target/2,identifier+'R'),(a,mid,target/2,identifier+'L')])
        total=c.fsum(sums);error=c.fsum(errors)
        result={'value':total,'tail':tail,'empiricalError':error,'acceptedPanels':len(sums),'precision':50}
        self.store.put(key+'/B/return',encode(c,result));require(error<=c.mpf('1e-13'),'summed adaptive error');return result

    def constant_reference(self,spec,arg):
        if len(spec['coefficientsAscending'])!=1:return None
        c=self.B;carrier=self.scalar(c,spec['carrier']);center=self.scalar(c,spec['center']);s=c.mpf(8)
        coeff=spec['coefficientsAscending'][0];coeff=c.mpc(self.scalar(c,coeff[0]),self.scalar(c,coeff[1]))
        native=(-3*c.j)**spec['timeOrder']*(c.j/5)**spec['spatialOrders'][1]*(c.j/10)**spec['spatialOrders'][2] if spec['role']=='X' else c.one
        norm=1/(2*c.pi) if spec['role']=='X' else c.one
        n=spec['spatialOrders'][0] if spec['role']=='X' else 0
        # Polynomial recurrence for derivatives of (i*argument)^n times the analytic Gaussian transform.
        poly=[c.zero]*n+[c.j**n]
        for _ in range(spec['argumentDerivative']):
            poly=poly_add(poly_derivative(poly),poly_mul([s*s*carrier-c.j*center,-s*s],poly))
        sign=1 if spec['role']=='X' else -1
        return coeff*native*norm*sign**spec['argumentDerivative']*poly_eval(poly,arg)*s*c.sqrt(2*c.pi)*c.exp(-s*s*(arg-carrier)**2/2-c.j*(arg-carrier)*center)

    def request(self,spec,argument_text):
        require(spec['role'] in ('X','Y') and spec['argumentDerivative'] in (0,1),'supported product/derivative')
        full={'spec':spec,'argument':argument_text,'Aprecision':30,'Bprecision':50,'orders':[24,48],'comparison':'absolute <1e-12'}
        key=polynomial_request_key(full)
        if key in self.cache:
            require(self.cache[key]['input']==full,'exact cached arguments');return self.cache[key]['return']
        self.store.put(key+'/input',full)
        # Argument descriptors are sum of fixed carrier and a fixed signed offset; no float keys.
        def argument(c):return self.scalar(c,argument_text['base'])+self.scalar(c,argument_text['offset'])*argument_text['sign']
        a=self.route_a(key,spec,argument(self.A),24);b=self.route_a(key,spec,argument(self.A),48);ref=self.route_b(key,spec,argument(self.B))
        c=self.B;av=c.mpc(a['value']);bv=c.mpc(b['value']);rv=ref['value'];analytic=self.constant_reference(spec,argument(c))
        differences={'A24-A48':abs(av-bv),'A48-B':abs(bv-rv),'A24-B':abs(av-rv)}
        if analytic is not None:differences.update({'A24-analytic':abs(av-analytic),'A48-analytic':abs(bv-analytic),'B-analytic':abs(rv-analytic)})
        result={'key':key,'A24':av,'A48':bv,'B':rv,'analyticConstantReference':analytic,'absoluteDifferences':differences,
            'tailBounds':[c.mpf(a['tail']),c.mpf(b['tail']),ref['tail']],'BempiricalError':ref['empiricalError'],
            'pass':all(v<c.mpf('1e-12') for v in differences.values()),'relativeAccuracyClaim':False,'fullActionAccuracyClaim':False}
        receipt=self.store.put(key+'/return',encode(c,result));require(result['pass'],'absolute Fourier comparison failed; stop without replacement orders')
        answer={'receipt':receipt,'values':result};self.cache[key]={'input':full,'return':answer};return answer
