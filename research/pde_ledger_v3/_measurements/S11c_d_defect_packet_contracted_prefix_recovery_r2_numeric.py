"""New contracted numerical evaluator; import is inert, use only in containment.

The old rules are restored, never generated. A24/A48/B50 and all purposes have
separate immutable request namespaces. No old bank or original producer calls.
Error quantities are positive empirical indicators, not interval enclosures.
"""
from S11c_d_defect_packet_contracted_prefix_recovery_r2 import panel_evidence
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
import math


def require(value,message):
    if value is not True:raise ValueError(message)


def encode(value):
    if hasattr(value,'_mpc_'):return {'mpc':[encode_mpf(x) for x in value._mpc_]}
    if hasattr(value,'_mpf_'):return encode_mpf(value._mpf_)
    if isinstance(value,Estimate):return {'value':encode(value.value),'error':encode(value.error)}
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    require(type(value) in (type(None),str,int,bool),'exact codec requires encoded numbers')
    return value


def encode_mpf(t):return {'mpf':[t[0],str(t[1]),t[2],t[3]]}


def decode(c,value):
    if isinstance(value,dict):
        if 'mpf' in value:
            require(set(value)<={'mpf','decimal'} and len(value['mpf'])==4,'original exact MPF')
            s,m,e,b=value['mpf'];require(type(s) is int and s in (0,1) and type(m) is str and type(e) is int and type(b) is int,'MP tuple types')
            original=(s,int(m),e,b);result=c.make_mpf(original)
            require(c.isfinite(result) and result._mpf_==original,'finite identical exact MPF');return result
        if set(value)=={'mpc'}:return c.make_mpc(tuple(decode(c,x)._mpf_ for x in value['mpc']))
        return {k:decode(c,v) for k,v in value.items()}
    if isinstance(value,list):return [decode(c,v) for v in value]
    return value


def rational(c,v):
    z=F(v);return c.mpf(z.numerator)/z.denominator


def complex_pair(c,v):return rational(c,v[0])+c.j*rational(c,v[1])


@dataclass(frozen=True)
class Estimate:
    value:object
    error:object
    def __post_init__(self):require(bool(self.error>=0),'nonnegative indicator')
    def __add__(self,other):
        if not isinstance(other,Estimate):other=Estimate(other,0)
        return Estimate(self.value+other.value,self.error+other.error)
    __radd__=__add__
    def __neg__(self):return Estimate(-self.value,self.error)
    def __sub__(self,other):return self+-other
    def __mul__(self,other):
        if not isinstance(other,Estimate):other=Estimate(other,0)
        return Estimate(self.value*other.value,abs(self.value)*other.error+abs(other.value)*self.error+self.error*other.error)
    __rmul__=__mul__
    def __truediv__(self,constant):
        require(not isinstance(constant,Estimate) and bool(constant!=0),'nonzero fixed numerical divisor')
        return Estimate(self.value/constant,self.error/abs(constant))


def paired_indicator(left,right,own):
    require(own in (24,48),'fixed paired order')
    chosen=left if own==24 else right
    return Estimate(chosen['value'],abs(left['value']-right['value'])+left['error']+right['error'])


def outer_decision(total,nested,eps):
    """Rank leaves by the same guard that is worst in the complete sum."""
    guards={(key,kind):value for key,v in total.items() for kind,value in
            [('total',v.error/eps),('nested',nested[key]/(eps/4))]}
    target=max(guards,key=guards.get);worst,kind=target
    def metric(part,part_nested,key):
        return part[key].error/eps if kind=='total' else part_nested[key]/(eps/4)
    return bool(guards[target]<=1),worst,metric,{'epsilon':eps,'weightedNestedCap':eps/4,'guards':{key+'/'+kind:value for (key,kind),value in guards.items()},'bindingGuard':kind,'noBudgetDonation':True}


def leaf_count(db):
    exists=db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='contracted_leaves'").fetchone()
    return {'initialized':exists is not None,'active':db.execute('SELECT count(*) FROM contracted_leaves WHERE active=1').fetchone()[0] if exists else 0}


def moment_polynomial(c,n,carrier,width):
    coeff=[c.one]
    for _ in range(n):
        new=[c.zero]*(len(coeff)+1)
        for j,a in enumerate(coeff):
            if j:new[j-1]+=j*a
            new[j]+=c.j*carrier*a;new[j+1]-=a/width**2
        coeff=new
    return coeff


def gaussian(c,p,carrier,center,width,n,kind,route,mutant=False):
    """Return full value and Gaussian-stripped amplitude from OWN formula.

    p is physical k/l, not transform argument; Y has argument -p and carrier
    -carrier. n refers only to the physical profile-axis derivative of X.
    Normalized X is k^n*Gu/(2pi); the proven i^n is carried by native alpha.
    """
    require(kind in ('X','Y') and type(mutant) is bool,'declared Fourier selector')
    nu=p-carrier if kind=='X' else carrier-p
    g=c.exp(-width**2*nu**2/2);base=width*c.sqrt(2*c.pi)*c.exp(-c.j*nu*center)
    if kind=='Y':amplitude=base
    elif route=='A':amplitude=base*(carrier if mutant else p)**n/(2*c.pi)
    else:
        # Derivative mutant is the SAME physical replacement on route B:
        # Q_n is replaced by the constant (i*p0)^n, not its original recurrence.
        coeff=[(c.j*carrier)**n] if mutant else moment_polynomial(c,n,carrier,width)
        moments=[c.one]
        if len(coeff)>1:moments.append(-c.j*width**2*nu)
        for r in range(1,len(coeff)-1):moments.append(-c.j*width**2*nu*moments[r]+r*width**2*moments[r-1])
        amplitude=base*sum((a*b for a,b in zip(coeff,moments)),c.zero)/((c.j**n)*2*c.pi)
    return g*amplitude,amplitude,g


def combine(m,qm,beta,Cj,Cd,z):
    """The four separately certified formulas plus the wrong-root mutant."""
    return {'J':Cj*z['Y0']*m/(qm*(qm+beta))*(z['C1T']+m*z['C0T']),
            'Dr':Cd*z['X1']*(z['Y1CT']+m*z['Y0CT']),
            'Dh':Cd*z['Y0']/qm*(z['C2T']+m*z['C1T']),
            'Dq':Cd*z['Y0']/qm*z['X02T'],
            'Dr_wrong_root':Cd*(2*z['X1T']*z['Y1C']+(z['X2T']-m*z['X1T'])*z['Y0C'])}


class Evaluator:
    def __init__(self,store,index,contexts,rules,geometry,math_receipt):
        import mpmath
        self.store=store;self.index=index;self.contexts=contexts;self.geometry=geometry;self.math_receipt=math_receipt
        self.mp=mpmath.mp.clone();self.mp.dps=60
        self.namespaces={};self.serial=0;self.rule_inputs=rules;self.rules={}
        for route,precision in [('A24',30),('A48',30),('B50',50)]:
            c=mpmath.mp.clone();c.dps=precision
            if route=='B50':
                r=rules['B-G7-K15'];xx=[decode(c,x) for x in r['kronrodNodes']];ww=[decode(c,x) for x in r['kronrodWeights']]
                gx=[decode(c,x) for x in r['gaussNodes']];gw=[decode(c,x) for x in r['gaussWeights']]
                embedded={i:gw[gx.index(x)] for i,x in enumerate(xx) if x in gx}
                require(len(embedded)==7 and len(xx)==len(ww)==15,'actual embedded restored rule')
            else:
                r=rules['A-GL'+route[1:]];xx=[decode(c,x) for x in r['nodes']];ww=[decode(c,x) for x in r['weights']];embedded={}
                require(len(xx)==len(ww)==int(route[1:]),'complete restored GL rule')
            require(r['precision']==precision and all(-1<x<1 for x in xx) and all(w>0 for w in ww),'saved rule/precision/open nodes')
            self.rules[route]=(c,xx,ww,embedded)
        with self.store.db:
            self.store.db.execute('''CREATE TABLE contracted_leaves (job TEXT, ordinal INTEGER, active INTEGER, receipt BLOB NOT NULL, PRIMARY KEY(job,ordinal))''')

    def ns(self,route,purpose):
        key=(route,purpose)
        if key not in self.namespaces:
            precision=self.rules[route][0].dps
            self.namespaces[key]=self.store.namespace({'route':route,'settings':{'purpose':purpose,'precision':precision,'privateCheckPrecisions':[30,50,60],'method':'constant-Gaussian-contracted-J-D-v1','mathematicalContextReceipt':self.math_receipt,'fieldsTemplatesUnits':self.contexts,'recordCapBytes':self.index.maximum,'diskReserveBytes':self.index.reserve,'cacheBytes':self.index.cache_limit}})
        return self.namespaces[key]

    def put(self,route,purpose,name,value):
        self.serial+=1
        return self.index._put(self.ns(route,purpose),str(self.serial)+'/'+name,encode(value))

    def request(self,route,purpose,kind,args,settings,callback):
        c=self.rules[route][0];ns=self.ns(route,purpose)
        descriptor={'route':route,'purpose':purpose,'precision':c.dps,'kind':kind,'context':{'mathematicalContextReceipt':self.math_receipt,'exactFamilyContexts':self.contexts},'arguments':encode(args),'settings':settings}
        old=self.index.lookup(ns,descriptor)
        if old is not None:
            self.put(route,purpose,'exact-request-reuse',{'descriptor':descriptor,'inputReceipt':old['inputReceipt'],'returnReceipt':old['receipt'],'fullInputCompared':True})
            return decode(c,old['value'])
        ticket=self.index.begin(ns,descriptor)
        # The full input is committed before callback evaluation; a thrown error
        # keeps this request PENDING. Never automatically recompute that request.
        value=callback();encoded=encode(value);self.index.complete(ticket,encoded);return decode(c,encoded)

    def q(self,c,p):
        rad=c.mpf(595)/100-p*p
        return c.sqrt(rad) if rad>=0 else c.j*c.sqrt(-rad)

    def quad(self,c,v):return rational(c,str(v.a))+rational(c,str(v.b))*c.sqrt(rational(c,str(v.d)))
    def line(self,c,line,m):return rational(c,str(line.slope))*m+self.quad(c,line.intercept)

    def params(self,c):
        a=c.mpf(100)/109+3*c.j/10*(c.mpf(100)/109);mu=c.mpf(3)/10;beta=a*mu
        return beta,a*mu**2*10/4,(-c.j*mu)*10/(4*c.j)

    def transforms(self,route,purpose,n,p,carrier,kind,parent):
        c=self.rules[route][0];center=-c.mpf(5)/2 if kind=='X' else c.mpf(5)/2
        mutant=purpose=='derivative-mutant' and kind=='X';mode='B' if route=='B50' else 'A'
        members=[r for r in self.contexts['entries'] if r['n']==n]
        role='source' if kind=='X' else 'consumer'
        selectors={r['specs'][role]['argumentDerivative'] for r in members}
        require(selectors=={0},'actual original transform selector remains zero')
        settings={'formula':mode,'kind':kind,'n':n,'mutant':mutant,'absoluteGate':'1e-24','strippedGate':'1e-24','argumentDerivative':next(iter(selectors)),'fixedUnit':'X' if kind=='X' else 'Y','parent':parent}
        def work():
            full,amp,g=gaussian(c,p,carrier,center,c.mpf(8),n,kind,mode,mutant)
            other=self.mp.clone();other.dps=30 if mode=='B' else 50
            op=decode(other,encode(p));oc=decode(other,encode(carrier));ot=decode(other,encode(center))
            alt,aamp,ag=gaussian(other,op,oc,ot,other.mpf(8),n,kind,'A' if mode=='B' else 'B',mutant)
            check=self.mp;vf=decode(check,encode(full));va=decode(check,encode(amp));vg=decode(check,encode(g))
            of=decode(check,encode(alt));oa=decode(check,encode(aamp));og=decode(check,encode(ag))
            decisions=[];normalized_error=check.zero
            for item in members:
                # Restore actual per-jet native coefficient, including P_j*i^n.
                coefficient=item['Xmultiplier' if kind=='X' else 'Ymultiplier']
                factor=complex_pair(check,coefficient);own_factor=complex_pair(c,coefficient);other_factor=complex_pair(other,coefficient)
                require(factor!=0,'eligible coefficient factor nonzero')
                # Multiply by the actual native factor in EACH route's own
                # context before conversion to the private comparison context.
                left=decode(check,encode(own_factor*full));right=decode(check,encode(other_factor*alt))
                la=decode(check,encode(own_factor*amp));ra=decode(check,encode(other_factor*aamp))
                difference=abs(left-right);reconstruction=abs(left-vg*la)+abs(right-og*ra)
                stripped=abs(la-ra);bound=check.mpf('1e-24')*max(1,abs(la),abs(ra))
                good=bool(difference<check.mpf('1e-24') and stripped<bound)
                decisions.append({'addressId':item['addressId'],'factor':factor,'actual':left,'oppositeFormula':right,'actualStripped':la,'oppositeStripped':ra,'fullDifference':difference,'strippedDifference':stripped,'strippedBound':bound,'reconstruction':reconstruction,'passed':good})
                normalized_error=max(normalized_error,(difference+reconstruction)/abs(factor))
            # Formula check has its own namespace; none of its values is inserted
            # in a baseline request index. Own-route value is returned unchanged.
            receipt=self.put(route,'formula-check:'+purpose,'Gaussian-formula-check',{'p':p,'carrier':carrier,'center':center,'width':8,'n':n,'kind':kind,'settings':settings,'ownValue':full,'ownAmplitude':amp,'ownGaussian':g,'oppositeValue':alt,'oppositeAmplitude':aamp,'oppositeGaussian':ag,'decisions':decisions})
            require(bool(decisions) and all(v['passed'] for v in decisions),'actual full and stripped Gaussian comparison')
            return {'value':full,'error':c.make_mpf(normalized_error._mpf_),'checkReceipt':receipt}
        result=self.request(route,purpose,'analytic-transform',{'physicalMomentum':p,'carrier':carrier,'center':center,'width':8,'kind':kind,'n':n,'memberIds':[v['addressId'] for v in members]},settings,work)
        return Estimate(result['value'],result['error'])

    def profile(self,c,t):
        z=5*c.pi*t
        if abs(z)<=c.mpf('1e-6'):
            z2=z*z;poly=1+z2/6+z2**2/120+z2**3/5040+z2**4/362880+z2**5/39916800
            value=1/(2*c.pi*poly);err=2*abs(z)**12/math.factorial(13)/(2*c.pi*poly**2)
            return Estimate(value,err),{'z':z,'branch':'sinhc-degree10','polynomial':poly,'remainderUpper':2*abs(z)**12/math.factorial(13)}
        return Estimate(5*t/(2*c.sinh(z)),c.zero),{'z':z,'branch':'native-sinh','empiricalRoundoffNotCertified':True}

    def inner_node(self,route,purpose,n,z,m,carrier,clipped,parent):
        c=self.rules[route][0]
        def work():
            qz,qm=self.q(c,z),self.q(c,m);beta,_,_=self.params(c)
            self.put(route,purpose,'inner-root-operands',{'z':z,'m':m,'qz':qz,'qm':qm,'sum':qz+qm,'beta':beta})
            require(qm!=0 and qz!=0 and qz+qm!=0,'unsampled grazing, no zero-over-zero assignment')
            pr,proof=self.profile(c,m-z)
            xx=self.transforms(route,purpose,n,z,carrier,'X',parent)/(qz+beta)
            yy=self.transforms(route,purpose,n,z,carrier,'Y',parent)/(qz+beta)
            zx=pr*xx;zy=pr*yy;zero=Estimate(c.zero,c.zero);select=lambda x:x if clipped else zero
            result={'Y0':zy,'X1':zx*z,'X1T':select(zx*z),'X2T':select(zx*z*z),'X02T':select(zx*qz*qz),
                    'C0T':select(zx*qz/(qm+qz)),'C1T':select(zx*z*qz/(qm+qz)),'C2T':select(zx*z*z*qz/(qm+qz)),
                    'Y0CT':select(zy/(qm+qz)),'Y1CT':select(zy*z/(qm+qz)),
                    'Y0C':zy/(qm+qz),'Y1C':zy*z/(qm+qz)}
            return {'values':encode(result),'nativeOperands':{'z':z,'m':m,'qz':qz,'qm':qm,'beta':beta,'profile':pr,'profileProof':proof,'XoverE':xx,'YoverE':yy,'clipped':clipped,'nativeNormal':1}}
        result=self.request(route,purpose,'inner-node',{'z':z,'m':m,'carrier':carrier,'n':n,'clipped':clipped}, {'parent':parent,'density':'all certified contractions; Y and X retain their own transform signs'},work)
        values={key:Estimate(v['value'],v['error']) for key,v in result['values'].items()}
        require(all(c.isfinite(v.value) and c.isfinite(v.error) for v in values.values()),'finite complete inner vector')
        return values

    def panel(self,route,purpose,kind,a,b,callback,context,order=None):
        """One bounded open panel at a time. Input/return nodes persist first."""
        c=self.rules[route][0]
        selected=route if route=='B50' else 'A'+str(order)
        _,nodes,weights,embedded=self.rules[selected];mid=(a+b)/2;half=(b-a)/2
        require(a<mid<b,'representable original panel midpoint')
        settings={'kind':kind,'selectedRule':selected,'ownRoute':route,'squareMapped':route!='B50','context':context,'ruleReceipt':self.contexts['ruleReceipts'][selected]}
        def work():
            sums={};errors={};gauss={};gauss_errors={};batch={'points':[],'weights':[],'jacobians':[],'values':[],'operands':[]}
            halves=(0,) if route=='B50' else (0,1)
            for side in halves:
                for index,(node,weight) in enumerate(zip(nodes,weights)):
                    if route=='B50':point=mid+half*node;jac=half;effective=weight*jac
                    else:
                        u=(node+1)/2;jac=2*half*u
                        point=a+half*u*u if side==0 else b-half*u*u;effective=weight*jac/2
                    parent={'panel':settings,'originalA':encode(a),'originalB':encode(b),'side':side,'ruleIndex':index,'ruleNode':encode(node),'representedPoint':encode(point),'jacobian':encode(jac)}
                    self.put(route,purpose,kind+'-node-input',parent)
                    require(a<point<b and jac>0,'strict open positive node and Jacobian')
                    try:value=callback(point,parent)
                    except BaseException:
                        if batch['points']:
                            self.serial+=1;self.index.put_panel(self.ns(route,purpose),str(self.serial)+'/'+kind+'-partial-node-batch',encode(batch))
                        raise
                    require(bool(value),'complete node vector')
                    if sums:require(set(sums)==set(value),'all returned vector components')
                    for key,v in value.items():
                        sums[key]=sums.get(key,c.zero)+effective*v.value
                        errors[key]=errors.get(key,c.zero)+abs(effective)*v.error
                        if index in embedded:
                            ew=embedded[index]*jac
                            gauss[key]=gauss.get(key,c.zero)+ew*v.value
                            gauss_errors[key]=gauss_errors.get(key,c.zero)+abs(ew)*v.error
                    batch['points'].append(point);batch['weights'].append(weight);batch['jacobians'].append(jac);batch['values'].append(value);batch['operands'].append(parent)
                    if len(batch['points'])==48:
                        self.serial+=1;self.index.put_panel(self.ns(route,purpose),str(self.serial)+'/'+kind+'-node-batch',encode(batch));batch={k:[] for k in batch}
            if batch['points']:
                self.serial+=1;self.index.put_panel(self.ns(route,purpose),str(self.serial)+'/'+kind+'-node-batch',encode(batch))
            result={}
            for key in sums:
                embedded_difference=abs(sums[key]-gauss[key]) if route=='B50' else c.zero
                nested=errors[key]+gauss_errors.get(key,c.zero)
                result[key]={'value':sums[key],'error':embedded_difference+nested,'innerError':nested,'embeddedDifference':embedded_difference}
            return {'components':result,'a':a,'b':b,'settings':settings}
        return self.request(route,purpose,kind+'-panel',{'a':a,'b':b},settings,work)

    def register_leaf(self,job,ordinal,payload):
        receipt=payload['receipt'];self.store.db.execute('INSERT INTO contracted_leaves VALUES (?,?,1,?)',(job,ordinal,json.dumps(receipt,sort_keys=True).encode()));self.store.db.commit()

    def leaf_rows(self,job):
        for ordinal,blob in self.store.db.execute('SELECT ordinal,receipt FROM contracted_leaves WHERE job=? AND active=1 ORDER BY ordinal',(job,)):
            receipt=json.loads(blob);value,actual=self.index.read_record(receipt['sequence']);require(actual==receipt,'active leaf complete receipt');yield ordinal,value

    def adaptive(self,route,purpose,kind,initial,callback,context,decision):
        """Global active-leaf sum, largest worst-component contribution; no timer."""
        require(route=='B50','physical adaptive route only');c=self.rules[route][0]
        self.serial+=1;job=str(self.serial)+'/'+kind;ordinal=0
        def add(a,b,meta,parent=None):
            nonlocal ordinal
            panel=self.panel(route,purpose,kind,a,b,lambda x,p:callback(x,p,meta),{**context,'cell':meta})
            record={'a':a,'b':b,'meta':meta,'panel':panel,'parent':parent,'ordinal':ordinal}
            receipt=self.put(route,purpose,'active-leaf',record);self.register_leaf(job,ordinal,{'receipt':receipt});ordinal+=1
        for a,b,meta in initial:add(a,b,meta)
        iteration=0
        while True:
            totals={};nested={}
            for _,raw in self.leaf_rows(job):
                row=decode(c,raw)
                for key,r in row['panel']['components'].items():
                    totals[key]=totals.get(key,Estimate(c.zero,c.zero))+Estimate(r['value'],r['error'])
                    nested[key]=nested.get(key,c.zero)+r['innerError']
            done,target_key,metric,details=decision(totals,nested)
            self.put(route,purpose,'adaptive-global-decision',{'job':job,'iteration':iteration,'totals':totals,'nested':nested,'done':done,'worstComponent':target_key,'details':details})
            if done:return totals,nested
            best=None;score=None
            for leaf_id,raw in self.leaf_rows(job):
                row=decode(c,raw);part={key:Estimate(totals[key].value,r['error']) for key,r in row['panel']['components'].items()}
                part_nested={key:r['innerError'] for key,r in row['panel']['components'].items()}
                candidate=metric(part,part_nested,target_key)
                if score is None or candidate>score:score=candidate;best=(leaf_id,row)
            require(best is not None and score>0,'positive refinable global error; no unknown promotion')
            leaf_id,row=best;a,b=row['a'],row['b'];mid=(a+b)/2
            self.put(route,purpose,'adaptive-refinement-input',{'job':job,'leafId':leaf_id,'a':a,'b':b,'mid':mid,'score':score,'details':details})
            require(a<mid<b,'refinement remains representable; no precision retry')
            add(a,mid,row['meta'],leaf_id);add(mid,b,row['meta'],leaf_id)
            with self.store.db:self.store.db.execute('UPDATE contracted_leaves SET active=0 WHERE job=? AND ordinal=?',(job,leaf_id))
            iteration+=1

    def primitive_vector(self,c,m,totals,selected):
        beta,cj,cd=self.params(c);qm=self.q(c,m);require(qm!=0,'outer open root')
        all_values=combine(m,qm,beta,cj,cd,totals)
        return {key:all_values[key] for key in selected}

    def selected_members(self,n,purpose):
        rows=[r for r in self.contexts['entries'] if r['n']==n]
        if purpose=='wrong-root-mutant':return [r for r in rows if r['addressId']==8347]
        if purpose=='derivative-mutant':return [r for r in rows if r['addressId']==8350]
        return rows

    def addressed(self,c,n,purpose,primitives):
        out={}
        for r in self.selected_members(n,purpose):
            alpha=complex_pair(c,r['alpha'])
            labels=['Dr'] if purpose=='wrong-root-mutant' else r['primitives']
            for p in labels:
                key='Dr_wrong_root' if purpose=='wrong-root-mutant' else p
                if key in primitives:out[str(r['addressId'])+'/'+p]=alpha*primitives[key]
        return out

    def inner(self,route,purpose,n,m,carrier,plan,slab,parent):
        c=self.rules[route][0];K,T,M=plan['K'],plan['T'],plan['M']
        keys=['J','Dr','Dh','Dq'] if purpose=='baseline' else (['J'] if purpose=='derivative-mutant' else ['Dr_wrong_root'])
        context={'n':n,'m':encode(m),'carrier':encode(carrier),'K':K,'T':T,'M':M,'outerParent':parent,'primitiveKeys':keys,'inputClips':'J,Dh,Dq,Dr_wrong_root','outputClips':'Dr','semanticVariables':{'X':'k','Y':'l'},'sharedQuadratureCoordinateDoesNotIdentifyKAndL':True}
        initial=[]
        for i,cell in enumerate(slab['cells']):
            a,b=self.line(c,cell['low'],m),self.line(c,cell['high'],m)
            require(a<b,'no float alias or omitted tiny inner interval')
            initial.append((a,b,{'cellId':i,'clipped':cell['clipped']}))
        require(initial[0][0]==-K and initial[-1][1]==K and all(initial[i][1]==initial[i+1][0] for i in range(len(initial)-1)),'actual full inner coverage')
        callback=lambda z,p,meta:self.inner_node(route,purpose,n,z,m,carrier,meta['clipped'],p)
        if route=='B50':
            eps=c.mpf('1e-11')/80;target=eps*(1+1/abs(self.q(c,m)))/(16*(3*M+4))
            def decision(total,nested):
                primitive=self.primitive_vector(c,m,total,keys);vals=self.addressed(c,n,purpose,primitive)
                ratios={key:v.error/target for key,v in vals.items()};worst=max(ratios,key=ratios.get)
                def metric(part,part_nested,address):return self.addressed(c,n,purpose,self.primitive_vector(c,m,part,keys))[address].error/target
                return bool(ratios[worst]<=1),worst,metric,{'pointTarget':target,'ratios':ratios,'actualPositiveProductPropagation':True}
            result,_=self.adaptive(route,purpose,'inner',initial,callback,context,decision)
        else:
            result={};own=24 if route=='A24' else 48
            for a,b,meta in initial:
                panels={order:self.panel(route,purpose,'inner',a,b,lambda z,p:callback(z,p,meta),{**context,'cell':meta},order) for order in (24,48)}
                require(set(panels[24]['components'])==set(panels[48]['components']),'complete inner paired component set')
                local={}
                for key,chosen in panels[own]['components'].items():
                    left,right=panels[24]['components'][key],panels[48]['components'][key]
                    local[key]=paired_indicator(left,right,own)
                    result[key]=result.get(key,Estimate(c.zero,c.zero))+local[key]
                self.put(route,purpose,'inner-positive-panel-pair',{'a':a,'b':b,'meta':meta,'pairedPanels':panel_evidence(panels),'ownOrder':own,'positiveComponentContributions':local})

        primitive=self.primitive_vector(c,m,result,keys)
        self.put(route,purpose,'completed-inner-contractions',{'context':context,'contractions':result,'primitives':primitive,'addressedDensities':self.addressed(c,n,purpose,primitive)})
        return self.addressed(c,n,purpose,primitive)

    def action(self,route,purpose,n,plan):
        c=self.rules[route][0];carrier=self.quad(c,plan['carrier']);K,T=plan['K'],plan['T'];eps=c.mpf('1e-11')/80
        parent={'n':n,'K':K,'T':T,'carrier':encode(carrier),'purpose':purpose,'epsilon':encode(eps),'geometryReceipt':self.contexts['geometryReceipts'][str(K)+'/'+str(plan['carrier'].b)]}
        initial=[(self.quad(c,s['left']),self.quad(c,s['right']),{'slabId':i}) for i,s in enumerate(plan['slabs'])]
        require(initial[0][0]==-plan['M'] and initial[-1][1]==plan['M'],'full numerical wings')
        require(all(a<b for a,b,_ in initial) and all(initial[i][1]==initial[i+1][0] for i in range(len(initial)-1)),'complete represented outer intervals')
        def callback(m,p,meta):return self.inner(route,purpose,n,m,carrier,plan,plan['slabs'][meta['slabId']],p)
        if route=='B50':
            def decision(total,nested):
                return outer_decision(total,nested,eps)
            totals,nested=self.adaptive(route,purpose,'outer',initial,callback,parent,decision)
        else:
            totals={};nested={}
            for a,b,meta in initial:
                panel=self.panel(route,purpose,'outer',a,b,lambda x,p:callback(x,p,meta),{**parent,'cell':meta},int(route[1:]))
                for key,r in panel['components'].items():
                    totals[key]=totals.get(key,Estimate(c.zero,c.zero))+Estimate(r['value'],r['error']);nested[key]=nested.get(key,c.zero)+r['innerError']
            self.put(route,purpose,'fixed-outer-budget-decision',{'totals':totals,'nested':nested,'epsilon':eps,'cap':eps/4})
            require(all(v<=eps/4 for v in nested.values()),'fixed A weighted nested budget; no refinement retry')
        result={'route':route,'purpose':purpose,'n':n,'K':K,'T':T,'carrier':carrier,'totals':encode(totals),'nested':nested,'empiricalNotRigorous':True}
        self.put(route,purpose,'completed-action-subset',result);return result
