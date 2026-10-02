#!/usr/bin/env python3
"""Saved-source closed direct-kernel limits/bounds; no integral or finite solve."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback
ROOT=Path('/var/projects/toy_physics')
THREADS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')
HELPERS=('require','sha','save','replace_json','posthash_records','containment','Journal','decode','one_symbol','expanded_sinh_arguments')


def require(value,message):
    if value is not True:raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def save(path,value):
    with Path(path).open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())


def helper_ast(source):
    nodes=[n for n in ast.parse(source).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in HELPERS]
    require({n.name for n in nodes}==set(HELPERS),'complete unchanged helper set')
    return ast.Module(body=nodes,type_ignores=[])


def verify_gate(path,manifest_path,manifest):
    g=json.loads(Path(path).read_text())
    require(g['status']=='READY_FOR_ONE_CLOSED_GRAZING_INSTRUMENT','ready gate')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest pins')
    require(g['sourcePins']==manifest['sourcePins'],'source pin census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for key in ('sharedGuard','supervisor','methodRecord','buildReviewRecord','authority'):
        require(sha(g[key])==g[key+'Sha256'],'gate '+key)
    require(g['methodRecord']==manifest['methodRecord'],'method source route')
    method=json.loads(Path(g['methodRecord']).read_text());build=json.loads(Path(g['buildReviewRecord']).read_text())
    authority=json.loads(Path(g['authority']).read_text())
    require(method['jointIndependentMethodClearance'] is True and method['allChecksPassed'] is True,'actual method assessment')
    require(build['independentBuildClearance'] is True and build['allChecksPassed'] is True,'actual substantive build clearance')
    require(all(build['reports'][e]['literalVerdict']=='CLEAR FOR THIS BOUNDED CLOSED-GRAZING BUILD' for e in ('claude','grok')),'literal fresh build verdicts')
    for key in ('workerSha256','manifestSha256','sharedGuardSha256','supervisorSha256'):
        require(build[key]==g[key],'review/gate '+key)
    require(authority['boundedInstrumentAuthorized'] is True and authority['automaticScientificRetry'] is False,'standing bounded authority')
    require(g['scope']==manifest['scope'] and g['scientificRunsAuthorized']==1 and g['durationLimits'] is None,'one bounded no-deadline job')
    return g


def nonnegative_polynomial(J,name,expression,variables):
    """Sufficient exact certificate: nonnegative coefficients and integer powers."""
    J.emit(name+'-input',{'expression':expression,'variables':variables})
    expanded=sp.expand(expression);terms=[]
    for term in sp.Add.make_args(expanded):
        powers=term.as_powers_dict();coefficient=term
        for v in variables:
            power=powers.get(v,sp.S.Zero)
            require(power.is_integer is True and power.is_nonnegative is True,'nonnegative integer power '+name)
            coefficient=coefficient/v**power
        coefficient=sp.cancel(coefficient)
        terms.append({'coefficient':coefficient,'powers':[powers.get(v,sp.S.Zero) for v in variables]})
    J.emit(name+'-certificate',{'expanded':expanded,'terms':terms,'assumptions':'All listed variables are real nonnegative.'})
    require(all(v.is_nonnegative is True for v in variables),'certificate variable domain')
    require(all(x['coefficient'].is_number is True and x['coefficient'].is_nonnegative is True for x in terms),'nonnegative exact coefficients '+name)
    J.zero(name+'-reconstruction',expression,sum(x['coefficient']*sp.prod(v**n for v,n in zip(variables,x['powers'])) for x in terms))
    return terms


def finite_number(J,name,value):
    real,imag=value.as_real_imag();record={'value':value,'real':real,'imaginary':imag,'realFinite':real.is_finite,'imaginaryFinite':imag.is_finite}
    J.emit(name+'-components',record);J.zero(name+'-reconstruction',value,real+sp.I*imag)
    require(not value.free_symbols and real.is_real is True and imag.is_real is True and real.is_finite is True and imag.is_finite is True,name+' exact finite scalar')



def grade_coefficient(modes,coefficients,grade):
    """Join a saved return to its actual mode grade, never a magic index."""
    require(len(modes)==len(coefficients),'mode/return length')
    hits=[value for mode,value in zip(modes,coefficients) if tuple(mode[:2])==tuple(grade)]
    require(len(hits)==1,'unique saved mode grade')
    return hits[0]


def exact_nonzero_number(J,name,value):
    """Exact finite fraction with signed nonzero numerator/denominator components."""
    J.emit(name+'-input',{'value':value})
    num,den=sp.fraction(sp.cancel(value))
    parts=[[sp.cancel(c) for c in z.as_real_imag()] for z in (num,den)]
    J.emit(name+'-fraction',{'numerator':num,'denominator':den,'components':parts,
        'finite':[[c.is_finite for c in row] for row in parts],
        'signedNonzero':[[c.is_positive is True or c.is_negative is True for c in row] for row in parts]})
    require(not value.free_symbols,'constant nonzero certificate '+name)
    for label,z,row in zip(('numerator','denominator'),(num,den),parts):
        J.zero(name+'-'+label+'-components',z,row[0]+sp.I*row[1])
        require(all(c.is_real is True and c.is_finite is True for c in row),'finite components '+name)
        require(any(c.is_positive is True or c.is_negative is True for c in row),'signed nonzero component '+name)
    J.zero(name+'-fraction-reconstruction',value*den,num)
    return {'finiteNonzero':True,'numerator':num,'denominator':den,'components':parts}


def run_science(manifest,J,helpers):
    copies={};cache={};used={}
    for name,r in manifest['savedFiles'].items():
        p=Path(r['path']);require(sha(p)==r['sha256'] and p.stat().st_size==r['bytes'],'saved file '+name)
        dst=J.out/'saved'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dst)
        require(sha(dst)==r['sha256'],'saved copy '+name)
        copies[name]={'source':str(p),'path':str(dst.relative_to(J.out)),'sha256':r['sha256'],'bytes':r['bytes']}
    save(J.out/'saved-copy-index.json',copies)
    def load(name,field=None):
        if name not in cache:cache[name]=json.loads((J.out/'saved'/name).read_text())
        used.setdefault(name,[]).append(field)
        value=cache[name] if field is None else cache[name][field]
        return helpers['decode'](value)
    C=load('physical-factorization-input.json','coefficient');B=load('physical-factorization-return.json','B')
    one=helpers['one_symbol'];source=[C,B]
    k=one(source,'k');H=one(source,'H');Q=one(source,'increment_difference');oldomega=one(source,'omega');rho=one(source,'rho_m')
    oldDepths=[one(source,n) for n in ('q_i','q_h','q_s','q_o')]
    qi,qh,qs,qo=[sp.Symbol('grazing_'+n) for n in ('qi','qh','qs','qo')]
    depthmap=dict(zip(oldDepths,(qi,qh,qs,qo)))
    l,t=sp.symbols('grazing_output grazing_transfer',real=True)
    omega=sp.Symbol('grazing_unrestricted_frequency');delta=sp.Symbol('grazing_delta',nonnegative=True)
    cs=sp.Symbol('grazing_effective_speed',positive=True)
    mass=sp.Rational(1,10);tau=sp.Rational(1,10);lam=sp.Rational(1,100);freq=sp.Integer(3)
    L=sp.Integer(10);W=sp.Integer(1);edge2=sp.Rational(1,20)
    physical=json.loads(Path(manifest['physicalInput']).read_text())
    for key,value in [('rho_m',mass),('tau_A',tau),('Lambda_A_0',lam),('L_W',L),('W_0',W)]:
        require(sp.Rational(physical['parameters'][key])==value,'actual physical binding '+key)
    sheet=load('physical-sheet.json','frequency');edges=load('physical-sheet.json','edge')
    J.emit('binding-frequency-and-domain',{'listedInputFrequency':physical['parameters']['omega'],'actualSavedFrequency':sheet,'savedReuseFrequency':load('effective-speed-reuse-domain.json','newFrequency'),'edges':edges,'newOmega':omega,'omegaAssumptions':omega.assumptions0,'domain':{'cs':[1,2],'k':[-3,3],'l':[-3,3],'delta':[0,sp.Rational(1,10)]},'analyticContinuation':'kernel only; source/consumer fixed at real3'})
    require(sheet==freq==load('effective-speed-reuse-domain.json','newFrequency') and edges==[sp.Rational(1,5),sp.Rational(1,10)],'actual saved frequency/edge')
    require(omega.is_real is None and oldomega.is_real is True,'new unrestricted frequency, not old real assumption')
    mapping={oldomega:omega,rho:mass,H:t,Q:l-k,**depthmap}
    Cu=C.subs(mapping,simultaneous=True);Bu=B.subs(mapping,simultaneous=True)
    J.emit('saved-source-transport',{'coefficient':C,'factor':B,'actualMap':[[a,b] for a,b in mapping.items()],'transportedCoefficient':Cu,'transportedFactor':Bu,'oldBoundaryReplayed':False})
    def algebra():
        difference=sp.together(Cu-t*Bu);num,den=sp.fraction(difference);expanded=sp.expand(num)
        rules={qh**2:qi**2-t*(t+2*k),qs**2:qo**2+t*(2*l-t)}
        J.emit('complex-frequency-factorization-operands',{'difference':difference,'numerator':num,'denominator':den,'rules':[[a,b] for a,b in rules.items()]})
        for term in sp.Add.make_args(expanded):
            require(all(term.as_powers_dict().get(x,0) in (0,1,2) for x in (qh,qs)),'quadratic reduction domain')
        reduced=sp.expand(expanded.subs(rules,simultaneous=True));J.zero('complex-frequency-factorization',reduced,sp.S.Zero)
        lowerMixed=grade_coefficient(load('lower-boundary-operands.json','modes'),load('lower-boundary-return.json'),(1,1))
        lowerTransported=lowerMixed.subs(mapping,simultaneous=True)
        J.zero('minus-unrestricted-source-join',lowerTransported,Cu)
        for face,transported in [('plus',Cu),('minus',lowerTransported)]:
            J.zero(face+'-complex-contact',transported.subs({t:0,qh:qi,qs:qo},simultaneous=True),sp.S.Zero)
        return {'factorizationResidual':reduced,'complexContactZero':True,'scope':'Unrestricted-frequency algebra from saved coefficients, not repeated boundary construction'}
    J.stage('complex-frequency-identities',{'savedC':C,'savedB':B,'frequencyMap':[oldomega,omega]},algebra)
    beta=lam*omega/(mass*(1-sp.I*omega*tau));den=(qi+beta)*(qo+beta)
    Bc=-sp.I*omega*mass/den*(k*(2*l-t)/(qs+qo)+k*(t+2*k)*qi/(qh*(qh+qi))+qi**2/qh)
    R=qi*qo/den;A=lambda x:L*x/(4*sp.sinh(sp.pi*L*x/2))
    pref=W*L/(4*sp.I)*A(t)*A(l-k-t);G=pref*Bc
    J.zero('closed-factor-rearrangement',Bu*R,Bc)
    J.emit('closed-density',{'Bc':Bc,'density':G,'prefactor':pref,'beta':beta,'ordinaryMeasure':'dt','endpointsNotAssignedPointwise':True})
    numericMap={omega:freq};newfromold={H:t,Q:l-k,oldomega:freq,rho:mass,**depthmap}
    rawPre=load('raw-ordered-before-cancel.json');oldTransfer=one([rawPre['heightNumerator']],'increment_transfer');newfromold[oldTransfer]=t
    J.sinh_zero('actual-ordered-prefactor',rawPre['heightScale']*rawPre['heightNumerator'].subs(newfromold,simultaneous=True)*rawPre['jet'].subs(newfromold,simultaneous=True)/rawPre['phaseDivisor'],pref)
    measure=load('edge-delta-reduction.json');require(measure['plainReducedMiddleMeasure'] is True,'plain middle measure')
    lower=load('lower-boundary-operands.json');geometry=load('native-lower-geometry-operands.json')
    require(lower['extensionSign']==lower['labHeightSign']==-1 and lower['reference']==-sp.Rational(1,2),'native lower position/extension')
    normal=geometry['nativeNormal'];sigma=one([normal],'sigma_W');slope=one([normal],'w1_profile_d1')
    J.zero('native-lower-outward-sign',normal[3],-sp.S.One)
    J.zero('native-lower-slope-sign',normal[0],-sigma*slope/2)
    factors={};faceJets={};faceRaw={}
    for face,label in [(1,'plus'),(-1,'minus')]:
        closure=load(label+'-closure-operands.json');tr=load(label+'-trace-domain.json');oldclosed=load(label+'-closed-raw-increment.json');before=load(label+'-closed-before-cancel.json')
        J.zero(label+'-own-saved-mixed-coefficient',before['mixed'].subs(newfromold,simultaneous=True),Cu.subs(omega,freq))
        faceRaw[label]=before['raw'].subs(newfromold,simultaneous=True)
        J.sinh_zero(label+'-actual-raw-kernel',faceRaw[label],(pref*Bu).subs(omega,freq))
        native=closure['nativeDefinition'];z=one([native],'s11cc1_dtn_operator_lab_held_'+label)
        law=sp.expand(native).coeff(z);nativeMap={one([law],'omega'):omega,one([law],'rho_m'):mass,one([law],'Lambda_A_0'):lam,one([law],'tau_A'):tau}
        J.zero(label+'-continued-native-law',law.subs(nativeMap,simultaneous=True),lam/(mass**2*(1-sp.I*omega*tau)))
        J.zero(label+'-saved-real-law',closure['coefficient'],(beta/(mass*omega)).subs(omega,freq))
        J.zero(label+'-saved-external-factor',closure['factor'].xreplace(depthmap),R.subs(omega,freq))
        J.zero(label+'-zero-reference-trace',tr['zeroTrace'],sp.S.One)
        J.zero(label+'-separate-jet-sign',tr['normalJet'].xreplace(depthmap),sp.I*face*qo)
        eta=one([tr['height']],'eta_bg');profile=one([tr['height']],'w1_profile')
        J.zero(label+'-native-height-sign',tr['height'],face*eta*profile/2)
        for key,target in [('physical',G),('reference',G),('normalJet',sp.I*face*qo*G)]:
            J.sinh_zero(label+'-actual-closed-'+key,oldclosed[key].subs(newfromold,simultaneous=True),target.subs(omega,freq))
        J.zero(label+'-reference-factor',before['reference'].xreplace(depthmap),R.subs(omega,freq))
        J.zero(label+'-jet-factor',before['jet'].xreplace(depthmap),sp.I*face*qo*before['reference'].xreplace(depthmap))
        factors[label]=before['reference'].xreplace(depthmap)
        faceJets[label]=sp.cancel(before['jet'].xreplace(depthmap)/factors[label])
        J.zero(label+'-saved-jet-ratio',faceJets[label],tr['normalJet'].xreplace(depthmap))
        J.emit(label+'-face-domain',{'nativeLaw':law,'continuedMap':[[a,b] for a,b in nativeMap.items()],'height':tr['height'],'referenceTrace':tr['zeroTrace'],'jet':tr['normalJet'],'fullLowerNormal':normal,'scope':'Direct (1,1) increment only; full normal/slope provenance saved, first-shape iteration not revalidated'})
    # Join each original physical route before using its analytic outgoing continuation.
    nativeRoots={}
    for label,momentum in [('input',k),('heightRoute',k+t),('slopeRoute',l-t),('output',l)]:
        savedRoot=load('physical-sheet.json',label)
        oldCs=one([savedRoot],'increment_effective_bulk_speed')
        rootMap={**newfromold,oldCs:cs}
        transported=savedRoot.subs(rootMap,simultaneous=True)
        nativeRoots[label]=transported
        J.emit(label+'-native-sheet-operands',{'saved':savedRoot,'actualMap':[[a,b] for a,b in rootMap.items()],'transported':transported,'momentum':momentum})
        require(isinstance(transported,sp.Piecewise) and len(transported.args)==3,'native three-branch outgoing sheet')
        positive,negative,zero=transported.args
        radReal=freq**2/cs**2-edge2-momentum**2
        J.zero(label+'-positive-root',positive.expr,sp.sqrt(radReal))
        J.zero(label+'-decaying-root',negative.expr,sp.I*sp.sqrt(-radReal))
        for branch,sign in [(positive,1),(negative,-1)]:
            condition=branch.cond
            require(condition.rel_op in ('<','>'),'native strict radicand domain')
            signed=(condition.rhs-condition.lhs) if condition.rel_op=='<' else (condition.lhs-condition.rhs)
            J.zero(label+'-domain-'+str(sign),signed,sign*radReal)
        require(zero.expr is sp.nan and zero.cond is sp.true,'saved pointwise zero refusal preserved')
    # New exact inequality certificates supporting the independently assessed analytic proof.
    O=freq+sp.I*delta;betaD=beta.subs(omega,O);D=(1+tau*delta)**2+freq**2*tau**2
    bre=lam/mass*freq/D;bim=lam/mass*(delta*(1+tau*delta)+freq**2*tau)/D
    J.zero('continued-beta-components',betaD,bre+sp.I*bim)
    bm=sp.Rational(3000,11101);km=sp.sqrt(879)/20
    dmax=sp.Rational(1,10);Dmax=sp.Rational(11101,10000)
    J.zero('beta-denominator-gap',Dmax-D,(dmax-delta)*tau*(2+tau*(dmax+delta)))
    J.zero('beta-real-lower-bound-identity',bre-bm,lam/mass*freq*(Dmax-D)/(D*Dmax))
    u,v,x,y=sp.symbols('quadrant_u quadrant_v quadrant_x quadrant_y',nonnegative=True)
    nonnegative_polynomial(J,'quadrant-sum-gain',x*x+y*y+2*u*x+2*v*y,[u,v,x,y])
    J.zero('quadrant-modulus-identity',(u+x)**2+(v+y)**2-(u*u+v*v),x*x+y*y+2*u*x+2*v*y)
    p=sp.Symbol('grazing_real_momentum',real=True);rad=O**2/cs**2-edge2-p*p
    J.zero('outgoing-radicand-components',rad,(freq**2-delta**2)/cs**2-edge2-p*p+sp.I*2*freq*delta/cs**2)
    real,imag=sp.symbols('radicand_real radicand_imaginary',real=True)
    J.zero('radicand-modulus-square-gap',real**2+imag**2-real**2,imag**2)
    aa,b=sp.symbols('distance_absolute_momentum distance_radius',nonnegative=True)
    J.zero('distance-factorization',b*b-aa*aa,(b-aa)*(b+aa))
    kn=(freq**2-delta**2)/cs**2-edge2
    J.zero('kappa-lower-gap',kn-km**2,(freq**2-delta**2)*(1/cs**2-sp.Rational(1,4))+(dmax**2-delta**2)/4)
    J.zero('inverse-cs-bound',1/cs**2-sp.Rational(1,4),(2-cs)*(2+cs)/(4*cs**2))
    J.emit('domain-bound-certificate',{'beta':betaD,'real':bre,'imaginary':bim,'betaMinimum':bm,'kappaMinimum':km,'frequencySquareUpper':freq**2+dmax**2,'qSquareTriangleUpper':freq**2+dmax**2+edge2+9,'qMagnitudeBound':5,'endpointEnclosure':[-6,6],'envelopeCoefficients':[18+3*sp.Abs(t),43+3*sp.Abs(t)],'proof':'delta in [0,1/10], cs in [1,2], |k|,|l|<=3. Positive components give beta >= Re beta >= bm. q^2 has positive imaginary part for delta>0. Principal roots lie first quadrant. The checked quadrant identity gives |u+v|>=max(|u|,|v|). |q|^2>=|Re q^2|; distance product=(kappa+|p|)*dist >= km*dist. kappa_delta<3 and external momenta <=3 enclose endpoints in [-6,6].','endpointMajorant':'For each route: km^(-1/2) SUM over its two endpoints |t-a|^(-1/2).','localSetIntegralBound':'2 sqrt(2 m) per endpoint; multiplying by bounded compact prefactors gives uniform absolute continuity.'})
    require((16-freq**2-dmax**2).is_positive is True and (25-freq**2-dmax**2-edge2-9).is_positive is True,'loose compact magnitude constants')
    # Tail: positive decompositions and an exact elementary primitive, never an integral call.
    z=sp.Symbol('tail_offset',nonnegative=True);s=z+12
    nonnegative_polynomial(J,'tail-polynomial-envelope',12*s-(61+6*s),[z])
    nonnegative_polynomial(J,'tail-distance-gap',s-6-s/2,[z])
    nonnegative_polynomial(J,'tail-profile-argument-ratio',sp.Rational(3,2)*s-(s+6),[z])
    T=sp.Symbol('tail_T',positive=True);rate=10*sp.pi
    Ctail=(W*L/4)*(3*L**2/2)*(4*mass/bm**2)*(24*sp.sqrt(2)/sp.sqrt(km))
    primitive=sp.exp(-rate*T)*(T**3/rate+3*T*T/rate**2+6*T/rate**3+6/rate**4)
    J.zero('tail-antiderivative',sp.diff(primitive,T),-T**3*sp.exp(-rate*T))
    J.emit('uniform-tail-certificate',{'minimumT':12,'constant':Ctail,'rate':rate,'twoSidedBound':2*Ctail*sp.exp(30*sp.pi)*primitive,'ordinaryDensityBound':'C exp(30pi) |t|^3 exp(-10pi|t|), |t|>=12','profileProof':'For |x|>=1: |A(x)|=(L|x|/2)e^(-5pi|x|)/(1-e^(-10pi|x|)) <= L|x|e^(-5pi|x|). pi>3 and exp(x)>=1+x imply 1-e^(-10pi)>1/2. |Q|<=6 and |t|>=12 give |Q-t| between |t|/2 and 3|t|/2.','inverseRootProof':'All endpoints in [-6,6], distance>=|t|/2; each inverse-depth sum <=2sqrt(2)/(sqrt(km)sqrt(|t|)).','primitiveLimit':'The displayed polynomial times exp(-10piT) tends to zero as T tends to infinity.','integralEvaluated':False})
    # Rational a.e. external limits; internal endpoints excluded, as in reviewed proof.
    limit_i=-sp.I*omega*mass*k*(2*l-t)/(beta*(qo+beta)*(qs+qo))
    limit_both=-sp.I*omega*mass*k*(2*l-t)/(beta**2*qs)
    J.zero('input-grazing-limit',Bc.subs(qi,0),limit_i)
    limit_o=-sp.I*omega*mass/(beta*(qi+beta))*(k*(2*l-t)/qs+k*(t+2*k)*qi/(qh*(qh+qi))+qi**2/qh)
    J.zero('output-grazing-limit',Bc.subs(qo,0),limit_o)
    J.zero('simultaneous-grazing-limit',Bc.subs({qi:0,qo:0},simultaneous=True),limit_both)
    J.zero('opposite-momentum-radicand-identity',(l-t)**2-(k+t)**2, (l+k)*(l-k-2*t))
    J.zero('opposite-routes-coincide',((l-t)**2-(k+t)**2).subs(l,-k),sp.S.Zero)
    J.emit('collision-and-L1-conclusion',{'internalEndpoints':[-k+sp.sqrt(kn),-k-sp.sqrt(kn),l+sp.sqrt(kn),l-sp.sqrt(kn)],'sameMomentumAtGrazing':[0,-2*k,2*k,0],'oppositeMomentumAtGrazing':[0,-2*k,0,-2*k],'signedCases':'For k=-kappa the signs of 2k reverse. In l=-k, qs=qh as functions, not merely at endpoints.','inputLimit':limit_i,'outputLimit':limit_o,'bothLimit':limit_both,'reasoning':'Pointwise continuity away from the finite limiting endpoint set plus the certified uniformly integrable moving-root envelope and tight tail gives L1 convergence. At delta>0 contact is exactly zero and C=tB; t PV(1/t)=1. L1 convergence excludes a concentrated Dirac mass. Applies to the kernel only; no differentiation or rate in parameters is asserted.','proofStatus':'Exact source/algebraic certificates plus independently assessed analytic argument; not automated measure-theory proof.','pointwiseEndpointValuesAssigned':False})
    # Formal source-jet scope only. The old k=0 upper Fourier witness does not
    # establish a general or lower-face Fourier binding; no new one is asserted.
    sourceJets={}
    for side in ('Plus','Minus'):
        source00=load('THETA_BALANCE-retained-increment.json','source'+side)
        sourceJets[side.lower()]=source00.free_symbols
        J.emit(side+'-source-jet-domain',{'savedSource':source00,'formalJets':sorted(source00.free_symbols,key=str),
            'scope':'Finite formal source jets; no general/lower-face Fourier convention or incoming eigenmode is certified.'})
    matches=[('LEFT',sp.sqrt(sp.Rational(3,2)),sp.sqrt(595)/10),('RIGHT',sp.sqrt(sp.Rational(150,101)),sp.sqrt(601)/10)]
    for name,speed,normalMomentum in matches:
        J.zero(name+'-modal-match',freq**2/speed**2-edge2,normalMomentum**2)
        require((speed-1).is_positive is True and (2-speed).is_positive is True and (3-normalMomentum).is_positive is True,'match in compact domain')
    for row in manifest['rows']:
        inc=load(row+'-retained-increment.json','mixedCoefficient').xreplace(depthmap);grades=load(row+'-retained-increment.json','grades')
        require(grades==[1,1],'actual retained grades')
        if row.startswith('U'):
            require(inc==0,'saved empty U increment');J.emit(row+'-saved-zero',{'value':inc,'constructionReplayed':False});continue
        eps=load(row+'-retained-increment.json','epsilon');dplus=one([inc],'increment_raw_plus');dminus=one([inc],'increment_raw_minus')
        multipliers={}
        for label,dslot in [('plus',dplus),('minus',dminus)]:
            rawCoefficient=sp.diff(inc,dslot)
            J.emit(row+'-'+label+'-before-slot-division',{'coefficient':rawCoefficient,'reference':factors[label],'epsilon':eps})
            mult=sp.cancel(rawCoefficient/factors[label]/eps)
            J.emit(row+'-'+label+'-slot-operands',{'actualIncrement':inc,'slot':dslot,'epsilon':eps,'rawSlotCoefficient':rawCoefficient,'savedReferenceFactor':factors[label],'multiplier':mult})
            J.zero(row+'-'+label+'-slot-join',rawCoefficient,eps*mult*factors[label])
            require(not any(x in mult.free_symbols for x in (qi,qh,qs,qo,Q,H,oldomega,rho)),'no depth/old input dependence in multiplier')
            require(all(not x.name.startswith('increment_') for x in mult.free_symbols),'no transfer dependence')
            basis=sorted(sourceJets[label],key=str)
            require(mult.free_symbols<=set(basis),'actual source jets only')
            coeffs=[sp.diff(mult,a) for a in basis]
            J.zero(row+'-'+label+'-linear-source-jets',mult,sum(a*c for a,c in zip(basis,coeffs)))
            for j,c in enumerate(coeffs):
                finite_number(J,row+'-'+label+'-formal-jet-coefficient-'+str(j),c)
            J.emit(row+'-'+label+'-finite-formal-multiplier',{'basis':basis,'coefficients':coeffs,
                'multiplier':mult,'tIndependent':True,'omegaHeld':freq,'csIndependent':True,
                'domain':'Fixed finite formal jets. Any later physical Fourier/mode binding needs its own source/sign joins.',
                'generalFourierMapCertified':False,'lowerFaceFourierMapCertified':False})
            multipliers[label]=mult
        J.zero(row+'-both-slot-reconstruction',inc,eps*sum(d*multipliers[label]*factors[label] for label,d in [('plus',dplus),('minus',dminus)]))
        for label,side in [('plus','Plus'),('minus','Minus')]:
            rawKernel=load(row+'-retained-increment.json','rawKernel'+side).subs(newfromold,simultaneous=True)
            J.sinh_zero(row+'-'+label+'-raw-slot-normalization',rawKernel,faceRaw[label])
        rawDensity=load(row+'-retained-increment.json','rawRowDensity').subs(newfromold,simultaneous=True)
        target=eps*sum(multipliers[label] for label in ('plus','minus'))*G.subs(omega,freq)
        J.sinh_zero(row+'-actual-closed-row-density',rawDensity,target)
        J.emit(row+'-row-applicability-scope',{'actualRawRowDensity':rawDensity,'target':target,
            'expandedPhysicalRowRevalidated':False,'rootJoin':'All four native root routes checked separately; no old physical-row construction is replayed.',
            'scope':'L1 closed-kernel factor with finite formal source-jet coefficients, not a physical excitation map.'})

    # New controls with predeclared responses, on the new closed kernel only.
    speed=matches[0][1];kap=matches[0][2];testt=sp.Rational(1,10)
    test_qo=kap;test_qs=sp.sqrt(kap**2-testt**2)
    pole=sp.cancel(qi*Bu).subs(qi,0)
    expectedPole=-sp.I*omega*mass*k*(2*l-t)/(qo*(qs+qo))
    J.zero('missing-factor-pole-residue-identity',pole,expectedPole)
    physicalPole=expectedPole.subs({omega:freq,k:kap,l:0,t:testt,qo:test_qo,qs:test_qs},simultaneous=True)
    residueReal=sp.cancel(physicalPole/sp.I)
    profilePoint=pref.subs({k:kap,l:0,t:testt},simultaneous=True)
    J.zero('control-profile-evenness',A(-kap-testt),A(kap+testt))
    J.emit('control-profile-nonzero',{'prefactor':profilePoint,'positiveProfileValues':[A(testt),A(kap+testt)],'factorDomain':[testt,kap+testt]})
    require(A(testt).is_positive is True and A(kap+testt).is_positive is True,'nonzero full-density control profile')
    J.emit('missing-external-factor-control',{'bareCoefficient':Bu,'poleResidue':pole,'physicalPoint':{'cs':speed,'k':kap,'l':0,'t':testt,'qo':test_qo,'qs':test_qs},'physicalResidue':physicalPole,'residueOverI':residueReal,'positive':residueReal.is_positive,'expectedResponse':'nonzero simple qi pole without closure; not a numerical infinity sample'})
    require(residueReal.is_positive is True,'responsive missing-factor pole')
    controlPoint={cs:speed,k:kap,l:0,t:testt}
    outgoing=nativeRoots['slopeRoute'].subs(controlPoint,simultaneous=True)
    outputDepth=nativeRoots['output'].subs(controlPoint,simultaneous=True)
    J.zero('control-native-slope-depth',outgoing,test_qs)
    J.zero('control-native-output-depth',outputDepth,test_qo)
    wrong=-outgoing
    limitMap={omega:freq,k:kap,l:0,t:testt,qo:outputDepth,qs:outgoing}
    limitPoint=limit_i.subs(limitMap,simultaneous=True)
    wrongMap={**limitMap,qs:wrong}
    wrongLimit=limit_i.subs(wrongMap,simultaneous=True)
    fullPoint=profilePoint*limitPoint;wrongFull=profilePoint*wrongLimit
    J.emit('wrong-sheet-control',{'savedRoute':nativeRoots['slopeRoute'],'point':[[a,b] for a,b in controlPoint.items()],
        'actualMap':[[a,b] for a,b in limitMap.items()],'corruptMap':[[a,b] for a,b in wrongMap.items()],
        'outgoing':outgoing,'wrong':wrong,'actualClosedDensity':fullPoint,'wrongClosedDensity':wrongFull,
        'densityMovement':wrongFull-fullPoint,'movementPerNonzeroProfile':wrongLimit-limitPoint,
        'outgoingPositive':outgoing.is_positive,'wrongNonnegative':wrong.is_nonnegative,
        'expectedResponse':'Nonzero new closed-density movement AND refusal of native first-quadrant predicate.'})
    require(outgoing.is_positive is True and wrong.is_nonnegative is False,'native wrong-quadrant predicate')
    exact_nonzero_number(J,'actual-closed-limit-nonzero',limitPoint)
    exact_nonzero_number(J,'wrong-sheet-movement-nonzero',wrongLimit-limitPoint)
    J.zero('wrong-sheet-full-response-factor',wrongFull-fullPoint,profilePoint*(wrongLimit-limitPoint))
    # Use saved minus jet/reference for actual response; saved plus ratio is the
    # deliberately wrong lower sign. Both ratios joined their own native traces.
    actualRatio=faceJets['minus'].subs(qo,outputDepth)
    wrongRatio=faceJets['plus'].subs(qo,outputDepth)
    actualJet=actualRatio*fullPoint;wrongJet=wrongRatio*fullPoint
    jetMovement=(wrongRatio-actualRatio)*limitPoint
    J.emit('lower-jet-sign-control',{'savedMinusRatio':faceJets['minus'],'savedPlusRatio':faceJets['plus'],
        'outputDepthFromNativeSheet':outputDepth,'inputGrazingClosedFactor':limitPoint,'profileFactor':profilePoint,
        'actualJet':actualJet,'wrongJet':wrongJet,'residual':wrongJet-actualJet,
        'movementPerNonzeroProfile':jetMovement,'expectedResponse':'Saved upper sign substituted for saved lower sign yields nonzero actual jet response.'})
    exact_nonzero_number(J,'lower-jet-movement-nonzero',jetMovement)
    J.zero('lower-jet-full-response-factor',wrongJet-actualJet,profilePoint*jetMovement)
    J.emit('restored-field-index',{'fields':used,'completedFunctionsReplayed':False})
    return {'executionStatus':'COMPLETED_CLOSED_GRAZING_CERTIFICATES','kernelL1ScopeSupported':True,'bothFaces':True,'complexFrequencyContactCertified':True,'sourceRowFrequency':3,'sourceRowScope':'Finite formal jets only; general/lower Fourier and incoming-mode binding not certified','physicalParametersChanged':False,'middleIntegralEvaluated':False,'finiteSolves':0,'productionChanges':False,'physicalLossClaim':False,'completedFunctionsReplayed':False,'analyticArgument':'Reviewed uniform-integrability/tail proof joined to actual expressions; not a formal theorem prover.','exclusions':['kappa=0','beta=0','first-shape iteration grazing','whole operator','untruncated inverse','loss','drain','primitive calibration','pointwise endpoint values','parameter differentiability']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--gate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    manifest=json.loads(args.inputs.read_text());gate=verify_gate(args.gate,args.inputs,manifest)
    pins={**manifest['sourcePins'],str(args.inputs):sha(args.inputs),str(args.gate):sha(args.gate)}
    args.out.resolve().relative_to(ROOT/'_scratch/s11c');args.out.mkdir(exist_ok=False)
    J=None;result={};code=1;started=time.monotonic()
    try:
        source=Path(manifest['helperSource']).read_text();ns={'ast':ast,'hashlib':hashlib,'json':json,'os':os,'Path':Path,'resource':resource,'THREADS':THREADS}
        exec(compile(helper_ast(source),manifest['helperSource']+'#unchanged-helpers','exec'),ns)
        save(args.out/'containment.json',ns['containment']())
        global sp
        import sympy as sp
        from sympy.core.symbol import Str
        ns.update(sp=sp,Str=Str);J=ns['Journal'](args.out)
        result=J.stage('closed-grazing-certificates',{'manifestSha256':gate['manifestSha256'],'methodRecordSha256':gate['methodRecordSha256']},lambda:run_science(manifest,J,ns));code=0
    except BaseException:
        result={'executionStatus':'FAILED_PRESERVED','traceback':traceback.format_exc(),'incompleteOperation':None if J is None else J.active,'automaticRetry':False};save(args.out/'failure.json',result)
    finally:
        if (args.out/'saved-copy-index.json').exists():
            copies=json.loads((args.out/'saved-copy-index.json').read_text());pins.update({str(args.out/v['path']):v['sha256'] for v in copies.values()})
        records={}
        for path,expected in pins.items():
            try:records[path]={'expected':expected,'actual':sha(path),'error':None}
            except OSError as e:records[path]={'expected':expected,'actual':None,'error':str(e)}
        save(args.out/'posthashes.json',records)
        if any(v['expected']!=v['actual'] for v in records.values()):result['integrityFailure']=True;code=1
        result.update(wallSeconds=time.monotonic()-started,scientificAcceptance=False)
        save(args.out/'checks.json',result);sys.stdout.write((args.out/'checks.json').read_text())
    return code

if __name__=='__main__':sys.exit(main())
