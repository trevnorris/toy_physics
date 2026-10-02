#!/usr/bin/env python3
"""Bounded retained reference-response certificates; no integral or finite solve."""
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

HELPERS=('require','sha','save','replace_json','posthash_records','containment','Journal','decode','one_symbol','expanded_sinh_arguments','named','function_source','assignment_source')

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
    require(g['status']=='READY_FOR_ONE_REFERENCE_GRAZING_INSTRUMENT','ready gate')
    require(g['workerSha256']==sha(__file__) and g['manifestSha256']==sha(manifest_path),'worker/manifest pins')
    require(g['sourcePins']==manifest['sourcePins'],'source pin census')
    for p,h in g['sourcePins'].items():require(sha(p)==h,'source pin '+p)
    for key in ('sharedGuard','supervisor','methodRecord','buildReviewRecord','authority'):
        require(sha(g[key])==g[key+'Sha256'],'gate '+key)
    require(g['methodRecord']==manifest['methodRecord'],'method source route')
    method=json.loads(Path(g['methodRecord']).read_text());build=json.loads(Path(g['buildReviewRecord']).read_text())
    authority=json.loads(Path(g['authority']).read_text())
    require(method['status']=='PAIRED_INDEPENDENT_METHOD_CLEARANCE_SOURCE_ONLY' and method['allChecksPassed'] is True,'actual method assessment')
    require(all(method['reviewers'][e]['literalVerdict']=='CLEAR FOR THIS BOUNDED REFERENCE-GRAZING METHOD' for e in ('claude','grok')),'literal paired method verdicts')
    require(method['methodSha256']==sha(manifest['methodPath']),'exact method source')
    require(build['independentBuildClearance'] is True and build['allChecksPassed'] is True,'actual substantive build clearance')
    require(all(build['reports'][e]['literalVerdict']=='CLEAR FOR THIS BOUNDED REFERENCE-GRAZING BUILD' for e in ('claude','grok')),'literal fresh build verdicts')
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

def retained(expression,eta,sigma):
    expanded=sp.expand(expression)
    for term in sp.Add.make_args(expanded):
        powers=term.as_powers_dict()
        for v in (eta,sigma):
            n=powers.get(v,sp.S.Zero)
            require(n.is_integer is True and n.is_nonnegative is True and n in (0,1,2),'native polynomial retained-grade domain')
        coefficient=term/(eta**powers.get(eta,0)*sigma**powers.get(sigma,0))
        require(not coefficient.has(eta,sigma),'grade coefficient independence')
    return sp.Add(*(expanded.coeff(eta,a).coeff(sigma,b)*eta**a*sigma**b
                    for a in (0,1) for b in (0,1)))


def mode_value(modes,values,grade):
    require(len(modes)==len(values),'actual saved mode lengths')
    hits=[v for m,v in zip(modes,values) if tuple(m[:2])==grade]
    require(len(hits)==1,'unique actual saved mode grade')
    return hits[0]


def run_science(manifest,J,helpers):
    copies={};cache={};used={}
    for name,r in manifest['savedFiles'].items():
        p=Path(r['path']);require(sha(p)==r['sha256'] and p.stat().st_size==r['bytes'],'saved '+name)
        dst=J.out/'saved'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dst)
        require(sha(dst)==r['sha256'],'byte copy '+name)
        copies[name]={'source':str(p),'path':str(dst.relative_to(J.out)),'sha256':r['sha256'],'bytes':r['bytes']}
    save(J.out/'saved-copy-index.json',copies)
    def load(name):
        if name not in cache:cache[name]=helpers['decode'](json.loads((J.out/'saved'/name).read_text()))
        used[name]=True
        return cache[name]
    one=helpers['one_symbol'];named=helpers['named']
    native=load('selected-trace.json');routes=load('function-routes.json')['routes']
    for rec in native.values():require(sha(rec['source'])==rec['sourceSha256'],'native source bytes')
    for rec in routes:
        text=Path(ROOT/rec['path']).read_text();fragment=helpers['function_source'](text,rec['function'])
        require(sha(ROOT/rec['path'])==rec['sourceSha256'],'function original source')
        require(fragment.rstrip()==rec['sourceText'].rstrip(),'actual function route '+rec['function'])
    c2=Path(manifest['c2Source']).read_text()
    dtn=sp.sympify(native['dtn']['constructor'],locals={'Str':helpers['Str']})
    first=named(dtn,'FIRST_SHAPE')
    omega=sp.Symbol('reference_unrestricted_frequency');qi,qm,qo=sp.symbols('reference_qi reference_qm reference_qo')
    k,l,m,t=sp.symbols('reference_k reference_l reference_m reference_t',real=True)
    cs=sp.Symbol('reference_cs',positive=True);delta=sp.Symbol('reference_delta',nonnegative=True)
    eta=one([first],'eta_bg');sigma=one([first],'sigma_W')
    h=sp.Function('reference_height_hat');j=sp.Function('reference_slope_hat')
    rho=sp.Rational(1,10);W=sp.S.One;L=sp.Integer(10);tau=sp.Rational(1,10);lam=sp.Rational(1,100);freq=sp.Integer(3);edge2=sp.Rational(1,20)
    mu=rho*omega;a=lam/(rho**2*(1-sp.I*omega*tau));beta=a*mu
    O=freq+sp.I*delta;b0=sp.Rational(3000,11101);kmin=sp.sqrt(879)/20
    phys=json.loads(Path(manifest['physicalInput']).read_text())
    for key,value in [('rho_m',rho),('W_0',W),('L_W',L),('tau_A',tau),('Lambda_A_0',lam)]:require(sp.Rational(phys['parameters'][key])==value,'physical binding '+key)
    sheet=load('physical-sheet.json');require(sheet['frequency']==freq and sheet['edge']==[sp.Rational(1,5),sp.Rational(1,10)],'actual frequency/edge')
    J.emit('physical-domain',{'frequency':freq,'frequencyContinuation':O,'edges':sheet['edge'],'cs':[1,2],'k/l':[-3,3],'delta':[0,sp.Rational(1,10)],'physicalInputFrequencyIsHistorical':phys['parameters']['omega'],'unrestrictedAssumptions':[[x,x.assumptions0] for x in [omega,qi,qm,qo]],'sourceConvention':'per unit native source; eta and sigma independent'})
    require(all(x.is_nonzero is None for x in [qi,qm,qo]) and omega.is_real is None,'unrestricted continuation symbols')
    oldcs=one([first],'c_s0')
    mapping={one([first],'omega'):omega,one([first],'rho_m'):rho,one([first],'W_0'):W,
        one([first],'s11cc1_q_out_input'):qi,one([first],'s11cc1_q_out_output'):qo,
        one([first],'s11cc1_k_input_1'):k,one([first],'s11cc1_k_input_2'):sp.Rational(1,5),one([first],'s11cc1_k_input_3'):sp.Rational(1,10),
        one([first],'s11cc1_w1_profile_hat_transfer'):2*h(l-k)/W,
        one([first],'s11cc1_w1_profile_jet_hat_1'):2*j(l-k),
        one([first],'s11cc1_w1_profile_jet_hat_2'):sp.S.Zero,one([first],'s11cc1_w1_profile_jet_hat_3'):sp.S.Zero}
    transported=first.subs(mapping,simultaneous=True)
    on_shell=transported.subs(omega**2/oldcs**2,qi**2+k**2+edge2)
    Zh=sp.I*mu*(qo-qi)/qo*h(l-k);Zj=mu*k/(qo*qi)*j(l-k)
    J.emit('first-shape-native-transport',{'restored':first,'actualMap':[[x,y] for x,y in mapping.items()],'transported':transported,'dispersionReplacement':[omega**2/oldcs**2,qi**2+k**2+edge2],'onShell':on_shell,'proposed':eta*Zh+sigma*Zj})
    J.zero('native-first-shape',on_shell,eta*Zh+sigma*Zj)
    flat=named(dtn,'FLAT_DIAGONAL');deltas=flat.atoms(sp.DiracDelta)
    J.emit('native-flat-delta-operand',{'flat':flat,'deltas':list(deltas),'coefficientAfterDeltaStrip':flat.xreplace({d:sp.S.One for d in deltas})})
    require(len(deltas)==3,'native flat three-momentum delta census')
    J.zero('native-flat-coefficient',flat.xreplace({d:sp.S.One for d in deltas}).subs(mapping,simultaneous=True),mu/qo)
    lower=load('lower-boundary-operands.json');lower_values=load('lower-boundary-return.json')
    for grade,target in [((1,0),sp.I*mu*(qo-qi)/qo),((0,1),mu*k/(qo*qi))]:
        value=mode_value(lower['modes'],lower_values,grade)
        lm={s:({'k':k,'omega':omega,'rho_m':rho,'q_i':qi,'q_h':qo,'q_s':qo}.get(s.name,s)) for s in value.atoms(sp.Symbol)}
        J.emit('lower-first-'+str(grade[0])+str(grade[1])+'-saved',{'grade':grade,'value':value,'map':[[x,y] for x,y in lm.items()]})
        J.zero('lower-first-'+str(grade[0])+str(grade[1]),value.subs(lm,simultaneous=True),target)
    R=lambda q:q/(q+beta)
    z01,z12,D=sp.symbols('reference_z01 reference_z12 reference_direct_bare_tag')
    Z=sp.Matrix([[mu/qo,z01,D],[0,mu/qm,z12],[0,0,mu/qi]])
    F=sp.Matrix([[mu/(qo+beta),R(qo)*z01*R(qm),R(qo)*D*R(qi)-a*R(qo)*z01*R(qm)*z12*R(qi)],
                 [0,mu/(qm+beta),R(qm)*z12*R(qi)],[0,0,mu/(qi+beta)]])
    J.emit('unrestricted-closure-operands',{'nativeEquation':'(I+a Z) F=Z','Z':Z,'coefficient':a,'proposedF':F,'oldInverseConstructorCalled':False})
    eq=(sp.eye(3)+a*Z)*F-Z
    for row in range(3):
        for col in range(3):J.zero('unrestricted-closure-'+str(row)+str(col),eq[row,col],sp.S.Zero)
    faces={};nativeTraceAssignments={}
    for face,label in [(1,'plus'),(-1,'minus')]:
        old=load(label+'-closure-operands.json');tr=load(label+'-trace-domain.json')
        nativeLaw=sp.expand(old['nativeDefinition']).coeff(one([old['nativeDefinition']],'s11cc1_dtn_operator_lab_held_'+label))
        lawmap={one([nativeLaw],'omega'):omega,one([nativeLaw],'rho_m'):rho,one([nativeLaw],'Lambda_A_0'):lam,one([nativeLaw],'tau_A'):tau}
        J.zero(label+'-unrestricted-law',nativeLaw.subs(lawmap,simultaneous=True),a)
        J.zero(label+'-real-law-anchor',a.subs(omega,freq),old['coefficient'])
        smap={s:({'q_o':qo,'q_h':qm,'q_i':qi,label+'_z01':z01,label+'_z12':z12,label+'_raw_direct':D}.get(s.name,s)) for s in old['matrix'].atoms(sp.Symbol)}
        J.emit(label+'-closure-restoration',{'savedMatrix':old['matrix'],'savedClosed':old['closed'],'symbolMap':[[x,y] for x,y in smap.items()],'newMatrix':Z,'newClosed':F})
        for row in range(3):
            for col in range(3):
                J.zero(label+'-matrix-anchor-'+str(row)+str(col),Z[row,col].subs(omega,freq),old['matrix'][row,col].xreplace(smap))
                J.zero(label+'-closed-anchor-'+str(row)+str(col),F[row,col].subs(omega,freq),old['closed'][row,col].xreplace(smap))
        savedq=one([tr['normalJet']],'q_o');normal_output=tr['normalJet'].xreplace({savedq:qo})
        htmap={one([tr['height']],'eta_bg'):eta,one([tr['height']],'w1_profile'):2*h(l-k)/W}
        labheight=tr['height'].subs(htmap,simultaneous=True)
        J.zero(label+'-height-binding',labheight,face*eta*h(l-k))
        J.zero(label+'-normal-binding',normal_output,sp.I*face*qo)
        n={'sp':sp,'qo':qo,'qi':qi,'MIDDLE_Q':qm,'normal_output':normal_output,'value_coefficient':tr['valueCoefficient'],
            'height_constant':sp.S.Zero,'height_hat':labheight,'height_kernel':labheight,'left':{k:m,qi:qm},'right':{l:m,qo:qm}}
        for target in ['normal_input','normal_middle','trace_two','trace_three']:
            fragment=helpers['assignment_source'](c2,'reference_pressure_kernels',target)
            nativeTraceAssignments[label+'-'+target]=fragment;exec(fragment,n)
        T=n['trace_three'];faces[label]={'T':T,'normal':normal_output,'height':labheight,'trace':tr,'T2':n['trace_two']}
        J.emit(label+'-new-native-trace',{'restoredTrace':tr,'actualHeightMap':[[x,y] for x,y in htmap.items()],'nativeAssignmentSources':{x:nativeTraceAssignments[label+'-'+x] for x in ['normal_input','normal_middle','trace_two','trace_three']},'newTrace':T,'newTwoLegTrace':n['trace_two'],'newDerivationNotRestored':True})
        J.zero(label+'-middle-trace',T[0,1],sp.I*eta*qm*h(l-m))
        J.zero(label+'-input-trace',T[1,2],sp.I*eta*qi*h(m-k))
    # Grade selection of the native physical closure; direct entry is absent in native second.
    left={k:m,qi:qm};right={l:m,qo:qm}
    zleft=(eta*Zh+sigma*Zj).xreplace(left);zright=(eta*Zh+sigma*Zj).xreplace(right)
    physical12=F[1,2].subs(z12,zright);physical01=F[0,1].subs(z01,zleft)
    physical02=F[0,2].subs({D:0,z01:zleft,z12:zright},simultaneous=True)
    P=F.copy();P[0,1]=physical01;P[1,2]=physical12;P[0,2]=retained(physical02,eta,sigma)
    D3=(qo+beta)*(qm+beta)*(qi+beta)
    phs=-sp.I*a*mu**2*k*(qo-qm)/D3
    psh=-sp.I*a*mu**2*m*qi*(qm-qi)/(qm*D3)
    tracehs=-sp.I*mu*k*qm/((qm+beta)*(qi+beta))
    C=-sp.I*mu*k*qo/((qo+beta)*(qi+beta))
    firstref=-sp.I*mu*qi/(qo+beta)*h(l-k)
    sloperef=mu*k/((qo+beta)*(qi+beta))*j(l-k)
    J.zero('physical-mixed-two-assignments',P[0,2].coeff(eta,1).coeff(sigma,1),phs*h(l-m)*j(m-k)+psh*j(l-m)*h(m-k))
    J.zero('height-left-plus-trace',phs+tracehs,C)
    refs={}
    for face,label in [(1,'plus'),(-1,'minus')]:
        T=faces[label]['T'];ref=P.copy()
        ref[0,1]=retained(P[0,1]-T[0,1]*P[1,1],eta,sigma)
        ref[1,2]=retained(P[1,2]-T[1,2]*P[2,2],eta,sigma)
        ref[0,2]=retained(P[0,2]-T[0,1]*ref[1,2],eta,sigma)
        J.emit(label+'-reference-before-guards',{'physical':P,'trace':T,'reference':ref,'normal':faces[label]['normal'],'uncombined':[phs,psh,tracehs],'noDirectInNativeSecond':True})
        remainder=(T*ref-P).applyfunc(lambda v:retained(v,eta,sigma))
        for row in range(3):
            for col in range(3):J.zero(label+'-reference-equation-'+str(row)+str(col),remainder[row,col],sp.S.Zero)
        J.zero(label+'-reference-height',ref[1,2].coeff(eta).subs({m:l,qm:qo},simultaneous=True),firstref)
        J.zero(label+'-reference-slope',ref[1,2].coeff(sigma).subs({m:l,qm:qo},simultaneous=True),sloperef)
        mixed=sp.expand(ref[0,2]).coeff(eta).coeff(sigma)
        J.zero(label+'-reference-mixed',mixed,C*h(l-m)*j(m-k)+psh*j(l-m)*h(m-k))
        # Use the actual affine trace coefficient and native final assignment with
        # a linear formal row action. This is a new kernel routing join, not a source-field contraction.
        pressure_slot,jet_slot,target=sp.symbols(label+'_pressure_slot '+label+'_jet_slot '+label+'_physical_target')
        H=sp.Symbol(label+'_affine_height');v=faces[label]['trace']['valueCoefficient']
        solve_form=(target-H*jet_slot)/v
        J.zero(label+'-affine-solve-equation',v*solve_form+H*jet_slot-target,sp.S.Zero)
        fragment=helpers['assignment_source'](c2,'build_face','reference_pressure')
        # For each input column the height acts on the middle row and outgoing
        # normal continuation. This preserves native ordering before projection.
        slotchecks=[]
        for col in [1,2]:
            action=T[0,1]/(sp.I*face*qm)
            n={'trace_map':{'REFERENCE_VALUE_SOLVE':solve_form,'PHYSICAL_PRESSURE_TARGET':target},'pressure':P[0,col],
               'jet_slot':jet_slot,'normal_jet':sp.I*face*qm*ref[1,col]}
            exec(fragment,n);final=retained(n['reference_pressure'].subs(H,action),eta,sigma)
            J.zero(label+'-native-final-slot-'+str(col),final,ref[0,col]);slotchecks.append(final)
        J.emit(label+'-final-native-slot-routing',{'source':fragment,'affineHeight':H,'savedHeight':faces[label]['height'],'valueCoefficient':v,'equationSolution':solve_form,'joinedColumns':slotchecks,'perUnitSource':True,'sourceCompositionPerformed':False})
        refs[label]=ref
    # Existing certified direct density is restored, not recomputed from its displayed formula.
    direct=load('closed-density.json')
    qr=sp.Symbol('reference_reflected_q')
    directMap={x:({'grazing_unrestricted_frequency':omega,'k':k,'grazing_output':l,'grazing_transfer':t,'grazing_qi':qi,'grazing_qh':qm,'grazing_qs':qr,'grazing_qo':qo}.get(x.name,x)) for x in direct['density'].atoms(sp.Symbol)}
    mappedDirect=direct['density'].subs(directMap,simultaneous=True)
    J.emit('restored-direct-argument-join',{'savedDensity':direct,'map':[[x,y] for x,y in directMap.items()],'transported':mappedDirect,'depthArguments':[[qi,k],[qm,k+t],[qr,l-t],[qo,l]],'measure':'dt; whole-convolution tag below represents this density integrated exactly once; value not evaluated'})
    require(mappedDirect.free_symbols <= {omega,k,l,t,qi,qm,qr,qo},'complete direct parameter map')
    J.zero('restored-direct-beta',direct['beta'].subs(directMap,simultaneous=True),beta)
    directTag=sp.Symbol('certified_closed_direct_whole_convolution')
    wholeH=sp.Symbol('whole_height_slope_convolution');wholeJ=sp.Symbol('whole_iterated_density_integral')
    J.emit('retained-response-census',{'grades':[[0,0],[1,0],[0,1],[1,1]],'discarded':[[2,0],[0,2]],
        'flat':mu/(qi+beta),'heightCoefficient':firstref,'slopeCoefficient':sloperef,
        'mixedIteration':C*wholeH+wholeJ,'taggedTotalMixed':directTag+C*wholeH+wholeJ,
        'certifiedDirectDensityRestored':direct,'directMultiplicity':1,'directNotNativeSecond':True,
        'middleIntegrationOfWholeDirect':False,'oldDirectFunctionsReplayed':False,'normalJet':{label:faces[label]['normal'] for label in faces},'jetKernels':{label:{'flat':faces[label]['normal']*mu/(qo+beta),'height':faces[label]['normal']*firstref,'slope':faces[label]['normal']*sloperef,'mixed':faces[label]['normal']*(directTag+C*wholeH+wholeJ)} for label in faces},'nativeRetentionSource':helpers['function_source'](c2,'retained_shape')})
    # The new contact/PV conversion is separate from the old direct contact certificate.
    A=lambda s:L*s/(4*sp.sinh(sp.pi*L*s/2))
    profiles=load('raw-ordered-before-cancel.json')
    profileMap={x:({'increment_transfer':t,'increment_difference':l-k}.get(x.name,x)) for value in (profiles['heightNumerator'],profiles['jet']) for x in value.atoms(sp.Symbol)}
    J.emit('restored-profile-arguments',{'original':profiles,'map':[[x,y] for x,y in profileMap.items()],'heightNumerator':A(t),'jet':L/2*A(l-k-t),'removableA0':1/(2*sp.pi),'noProfileTransformReplayed':True})
    J.sinh_zero('height-profile-factor-join',profiles['heightNumerator'].subs(profileMap,simultaneous=True),A(t))
    J.sinh_zero('slope-profile-factor-join',profiles['jet'].subs(profileMap,simultaneous=True),L/2*A(l-k-t))
    J.zero('height-scale-join',profiles['heightScale'],W/2)
    J.zero('profile-phase-join',profiles['phaseDivisor'],sp.I)
    Q=l-k;den=(qo+beta)*(qm+beta)*(qi+beta)
    Jdensity=a*mu**2*W*L/4*A(t)*A(Q-t)*(k+t)*(t+2*k)*qi/(qm*den*(qm+qi))
    psh_t=psh.subs(m,k+t)
    originalPV=psh_t*(W/2*A(t)/(sp.I*t))*(L/2*A(Q-t))
    # Replace the root difference at the coefficient level, preserving the pre-cancel operand.
    transformed_psh=-sp.I*a*mu**2*(k+t)*qi*(-t*(t+2*k)/(qm+qi))/(qm*den)
    J.emit('right-height-PV-operands',{'originalCoefficient':psh_t,'originalPV':originalPV,'differenceIdentity':[qm-qi,-t*(t+2*k)/(qm+qi)],'factoredCoefficient':transformed_psh,'density':Jdensity})
    J.zero('root-difference-polynomial',qm**2-qi**2, (qm-qi)*(qm+qi))
    J.zero('shared-dispersion-root-difference',((omega**2/cs**2-edge2-(k+t)**2)-(omega**2/cs**2-edge2-k**2)),-t*(t+2*k))
    difference=sp.together(psh_t-transformed_psh);num,denominator=sp.fraction(difference)
    J.emit('right-height-physical-factorization',{'nativeCoefficient':psh_t,'candidateCoefficient':transformed_psh,'rawDifference':difference,'numerator':num,'denominator':denominator,'physicalRule':[qm**2,qi**2-t*(t+2*k)]})
    expanded=sp.expand(num)
    require(all(term.as_powers_dict().get(qm,0) in (0,1,2) for term in sp.Add.make_args(expanded)),'quadratic middle numerator domain')
    J.zero('right-height-physical-factorization-residual',expanded.subs(qm**2,qi**2-t*(t+2*k)),sp.S.Zero)
    J.zero('right-height-PV-density',transformed_psh*(W/2*A(t)/(sp.I*t))*(L/2*A(Q-t)),Jdensity)
    J.zero('right-height-contact',psh.subs({m:k,qm:qi},simultaneous=True),sp.S.Zero)
    J.zero('physical-left-height-contact',phs.subs({m:l,qm:qo},simultaneous=True),sp.S.Zero)
    contact=C*W/4*(L/2*A(Q))
    J.zero('full-surviving-contact',(phs+tracehs).subs({m:l,qm:qo},simultaneous=True)*W/4*(L/2*A(Q)),contact)
    s=sp.Symbol('reference_left_height_transfer',real=True)
    Hsub=A(s)*(L/2*A(Q-s)-L/2*A(Q))/s
    J.emit('left-height-subtracted-PV',{'contact':W/4*(L/2*A(Q)),'ordinarySubtractedIntegrand':W/(2*sp.I)*Hsub,'variableChange':[m,l-s],'boundsAfterOrientation':'integrate s over the full real line in increasing order','removedTerm':A(s)*(L/2*A(Q))/s,'reason':'A even; symmetric PV of the odd removed term is zero. No integral value evaluated.'})
    J.zero('A-even',A(-s),A(s))
    J.zero('removed-PV-odd',(A(-s)/(-s)),-A(s)/s)
    J.zero('left-height-profile-subtraction',A(s)*(L/2*A(Q-s))/s,Hsub+A(s)*(L/2*A(Q))/s)
    # Native sheet for the actual middle leg; originals retain zero-radicand refusal.
    oldroot=sheet['heightRoute'];rootmap={x:({'k':k,'H':m-k,'increment_effective_bulk_speed':cs}.get(x.name,x)) for x in oldroot.atoms(sp.Symbol)}
    native_middle=oldroot.subs(rootmap,simultaneous=True)
    require(isinstance(native_middle,sp.Piecewise) and len(native_middle.args)==3,'native middle branch census')
    J.emit('native-middle-sheet',{'saved':oldroot,'map':[[x,y] for x,y in rootmap.items()],'transported':native_middle,'complexRadicand':O**2/cs**2-edge2-m*m})
    realrad=freq**2/cs**2-edge2-m*m
    kd=sp.sqrt((freq**2-delta**2)/cs**2-edge2)
    J.zero('middle-complex-radicand',O**2/cs**2-edge2-m*m,kd**2-m*m+sp.I*2*freq*delta/cs**2)
    J.zero('middle-propagating-sheet',native_middle.args[0].expr,sp.sqrt(realrad))
    J.zero('middle-decaying-sheet',native_middle.args[1].expr,sp.I*sp.sqrt(-realrad))
    for n,sign in [(0,1),(1,-1)]:
        cond=native_middle.args[n].cond;signed=cond.lhs-cond.rhs if cond.rel_op=='>' else cond.rhs-cond.lhs
        J.zero('middle-sheet-domain-'+str(n),signed,sign*realrad)
    require(native_middle.args[2].expr is sp.nan,'original endpoint refusal preserved')
    # Reuse certificates, checking actual law/domain operands. New bounds are elementary exact identities.
    domain=load('domain-bound-certificate.json');dbeta=domain['beta'];dmap={x:delta for x in dbeta.atoms(sp.Symbol) if x.name=='grazing_delta'}
    J.zero('reused-beta-actual-law',dbeta.xreplace(dmap),beta.subs(omega,O))
    J.zero('reused-beta-minimum',domain['betaMinimum'],b0);J.zero('reused-radius-minimum',domain['kappaMinimum'],kmin)
    J.emit('reused-domain-certificates',{'domain':domain,'quadrant':load('quadrant-sum-gain-certificate.json'),'priorTail':load('uniform-tail-certificate.json'),'sourceArgumentJoin':'Same cs/delta/edge domain; new middle p=k+t joined explicitly above. Original direct two-route envelope is not reused as the new J envelope.'})
    J.zero('new-a-denominator-modulus',sp.expand((1+tau*delta)**2+(freq*tau)**2-1),2*tau*delta+tau**2*delta**2+(freq*tau)**2)
    require((sp.Rational(4,25)-(rho**2*(freq**2+sp.Rational(1,100)))).is_positive is True,'new mu upper bound')
    J.zero('new-prefactor-constant',sp.Rational(4,25)*W*L/4,sp.Rational(2,5))
    z=sp.Symbol('reference_tail_offset',nonnegative=True);u=z+12
    nonnegative_polynomial(J,'new-tail-root-gap',(u-3)**2-9-u**2/2,[z])
    nonnegative_polynomial(J,'new-tail-polynomial',sp.Rational(15,8)*u**2-(u+3)*(u+6),[z])
    require((sp.Integer(200)**2-(sp.Rational(2,5)*150*sp.Rational(15,8))**2*2).is_positive is True,'new tail constant 200 sufficient')
    xx,yy,uu,vv=sp.symbols('reference_quadrant_x reference_quadrant_y reference_quadrant_u reference_quadrant_v',nonnegative=True)
    J.zero('Holder-quadrant-gap',(xx+uu)**2+(yy+vv)**2-((xx-uu)**2+(yy-vv)**2),4*(xx*uu+yy*vv))
    J.zero('Holder-radicand-difference',(O**2/cs**2-edge2-(k+s)**2)-(O**2/cs**2-edge2-k*k),-s*(2*k+s))
    p=sp.Symbol('reference_profile_nonnegative_argument',nonnegative=True)
    gap=sp.sinh(p)**2-p*sp.cosh(p)+sp.sinh(p)
    J.zero('profile-derivative-gap',sp.diff(gap,p),sp.sinh(p)*(2*sp.cosh(p)-p))
    J.zero('profile-derivative-gap-at-zero',gap.subs(p,0),sp.S.Zero)
    J.zero('cosh-positive-polynomial',p*p-p+2,(p-sp.Rational(1,2))**2+sp.Rational(7,4))
    J.zero('profile-A-derivative-form',sp.diff(A(s),s),L/4*(sp.sinh(sp.pi*L*s/2)-(sp.pi*L*s/2)*sp.cosh(sp.pi*L*s/2))/sp.sinh(sp.pi*L*s/2)**2)
    J.zero('jet-coefficient-variation',qo/(qo+beta)-qi/(qi+beta),beta*(qo-qi)/((qo+beta)*(qi+beta)))
    J.emit('new-bounds',{'Jdensity':Jdensity,'localEnvelope':sp.Rational(2,5)/b0**3*A(t)*A(Q-t)*(sp.Abs(t)+3)*(sp.Abs(t)+6)/sp.Abs(qm),'endpoints':[-k-kd,-k+kd],'endpointMeaning':'Actual moving real comparison endpoints; kmin is only a lower bound on their radius.','localSetBound':'4 sqrt(2 measureE)/sqrt(kmin) times bounded compact prefactors','tailEnvelope':200/b0**3*sp.exp(30*sp.pi)*sp.Abs(t)**3*sp.exp(-10*sp.pi*sp.Abs(t)),'tailThreshold':12,'profileDerivativeBound':L/4,'slopeDerivativeBound':L**2/8,'HLocalIntegrandBoundBeforeWOver2i':L**2/(16*sp.pi),'HOutsideUnitTailBeforeWOver2i':L**2/(2*sp.pi)*sp.exp(-sp.pi*L*sp.Abs(s)/2),'Cbound':sp.Rational(6,5)/b0,
        'heightStrongConstants':{'Lipschitz':sp.Abs(mu*qi)/b0,'Holder':sp.Abs(mu*qi)*sp.sqrt(7)/b0**2},'heightUniformConstants':[2/b0,2*sp.sqrt(7)/b0**2],'jetUniformConstants':[2,4*sp.sqrt(7)/(5*b0**2)],'testSpace':'fixed C_c^1((-3,3)) output test; fixed input k, uniform constants on compact family','proof':'Reviewed moving-endpoint uniform integrability plus explicit tail for J. Dominated convergence for the subtracted first-height distribution; no assertion q(l)*phi(l) is C1.'})
    J.emit('grazing-limit-objects',{'inputZeroJ':Jdensity.subs(qi,0),'outputZeroC':C.subs(qo,0),'inputZeroHeight':firstref.subs(qi,0),'outputZeroJ':Jdensity.subs(qo,0),'bothZeroJ':Jdensity.subs({qi:0,qo:0},simultaneous=True),'sameAndOppositeMatches':['l=k=+/-kappa','l=-k, |k|=kappa'],'contactEndpoint':'t=0 at qi=0; contact already zero for every delta>0, uniform L1 precludes a hidden mass','classes':{'flat':'diagonal delta','firstHeight':'delta plus PV; subtracted action','firstSlope':'ordinary kernel','mixedHeightLeft':'continuous C times H(Q), with explicit surviving contact inside H','mixedHeightRight':'ordinary L1 density J','direct':'restored separately certified L1 addend'},'pointwiseEndpointValuesAssigned':False,'integralsEvaluated':False,'measureTheoryMachineProved':False})
    require(Jdensity.subs(qi,0)==0 and firstref.subs(qi,0)==0 and C.subs(qo,0)==0,'vanishing a.e. limiting formulas')
    kap=sp.sqrt(freq**2/cs**2-edge2);collision_records=[]
    for sign in [-1,1]:
        for relation in [-1,1]:
            ki=sign*kap;lo=relation*ki;tag=str(sign)+'-'+str(relation)
            J.zero('collision-input-'+tag,freq**2/cs**2-edge2-ki**2,sp.S.Zero)
            J.zero('collision-output-'+tag,freq**2/cs**2-edge2-lo**2,sp.S.Zero)
            for where,transfer in [('right-height',sp.S.Zero),('left-height',lo-ki)]:
                J.zero('contact-middle-root-'+where+'-'+tag,freq**2/cs**2-edge2-(ki+transfer)**2,sp.S.Zero)
            collision_records.append({'sign':sign,'outputRelation':relation,'input':ki,'output':lo,'contactTransfers':[sp.S.Zero,lo-ki],'iteratedLimit':'C=0 and J=0 in its L1 class; isolated internal points are not assigned values'})
    J.emit('actual-collision-certificates',collision_records)
    J.zero('H-local-bound-constant',1/(2*sp.pi)*L**2/8,L**2/(16*sp.pi))
    J.emit('profile-bound-reasoning',{'evenProfile':A(s),'removableValue':1/(2*sp.pi),'sinhBound':'sinh(x)>=x for x>=0 gives |A|<=1/(2pi)','derivativeGap':gap,'gapDerivative':sp.sinh(p)*(2*sp.cosh(p)-p),'positiveRemainder':(p-sp.Rational(1,2))**2+sp.Rational(7,4),'coshBound':'cosh(x)>=1+x^2/2; therefore 2cosh(x)-x>=x^2-x+2>0','tailCondition':'|s|>=1 and L=10 imply 1-exp(-pi L |s|)>1/2','proofKind':'Source-bound exact expressions plus the independently assessed elementary inequalities, not an automatic theorem prover'})

    # Native Fourier/units are inherited operands; compare exact source text before dimensional arithmetic.
    dimensions=load('dimension-and-measure-join.json');contract=load('native-fourier-contract.json')
    for fn,field in [('profile_bindings','profiles'),('kernel_apply','kernelApply'),('fourier_profiles','fourierProfiles')]:
        require(helpers['function_source'](c2,fn).rstrip()==dimensions[field].rstrip(),'actual inherited Fourier source '+fn)
    require(contract==dimensions['executableContract'] and contract['invariantEdgeCoordinates']==2,'measure contract identity')
    require(load('edge-delta-reduction.json')['plainReducedMiddleMeasure'] is True,'plain reduced measure')
    vector=lambda xs:sp.Matrix(xs)
    dmu=vector(dimensions['nativeRhoDimension'])+vector(dimensions['nativeOmegaDimension']);dq=vector([-1,0,0]);dh=vector(dimensions['heightHat']);dj=vector(dimensions['jetHat']);dm=vector(dimensions['dt']);expected=vector(dimensions['savedReduced'])
    for name,value in [('height',dmu+dh),('slope',dmu-dq+dj),('mixedHS',dmu+dh+dj+dm),('mixedSH',dmu+dh+dj+dm)]:
        for n in range(3):J.zero('dimension-'+name+'-'+str(n),value[n],expected[n])
    J.emit('new-unit-joins',{'dimensionOrder':'length,time,mass; native rho has length^-4','nativeReduced':expected,'height':dmu+dh,'slope':dmu-dq+dj,'mixed':dmu+dh+dj+dm,'normalJetExtra':dq,'epsilon':'kernel per source amplitude; no source epsilon composition performed','etaSigma':[eta,sigma]})
    # Four controls through actual native closure/trace coefficients at rational-depth physical momenta.
    point={omega:freq,cs:sp.sqrt(sp.Rational(10,7)),k:sp.Rational(3,2),l:sp.Integer(2),m:sp.Rational(30,13),qi:sp.Integer(2),qo:sp.Rational(3,2),qm:sp.Rational(25,26)}
    middleActual=native_middle.subs(point,simultaneous=True)
    J.zero('control-native-middle-root',middleActual,point[qm])
    for name,depth,momentum in [('input',qi,k),('output',qo,l),('middle',qm,m)]:
        J.zero('control-physical-dispersion-'+name,(depth**2-(omega**2/cs**2-edge2-momentum**2)).subs(point,simultaneous=True),sp.S.Zero)
    profile=A((m-k).subs(point))*A((l-m).subs(point))*W*L/(4*sp.I)
    require((m-k).subs(point)!=0 and (l-m).subs(point)!=0,'applicable noncontact momenta')
    require(A((m-k).subs(point)).is_positive is True and A((l-m).subs(point)).is_positive is True,'actual nonzero tanh profile factors')
    reducedKernel=(C/(l-m)+psh/(m-k))
    omittedTrace=phs/(l-m)+psh/(m-k)
    wrongLower=(phs-tracehs)/(l-m)+psh/(m-k)
    omittedMiddle=(phs/R(qm)+tracehs)/(l-m)+(psh/R(qm))/(m-k)
    wrongRoot=reducedKernel.subs(qm,-qm)
    # The trace changes are obtained from the actual restored lower trace T01.
    lowerT=faces['minus']['T'][0,1];lowerNormal=sp.I*(-1)*qm
    # Actual Fj(m,k) uses INPUT k, not the middle output m.
    nativeTraceCoefficient=sp.cancel(-lowerT*(mu*k/((qm+beta)*(qi+beta))*j(m-k))/eta)
    J.zero('native-control-trace-coefficient',nativeTraceCoefficient,tracehs*h(l-m)*j(m-k))
    J.emit('control-actual-native-trace',{'height':faces['minus']['height'],'normal':lowerNormal,'T01':lowerT,'corruptHeight':-faces['minus']['height'],'corruptT01':-lowerT,'physicalMixedCoefficients':[phs,psh],'nativeTraceCoefficient':nativeTraceCoefficient})
    controls={}
    for name,corrupt in [('omit-trace',omittedTrace),('wrong-lower-height',wrongLower),('omit-middle-resolvent',omittedMiddle),('wrong-middle-root',wrongRoot)]:
        baseline=reducedKernel.subs(point,simultaneous=True);changed=corrupt.subs(point,simultaneous=True);move=sp.cancel(changed-baseline)
        J.emit(name+'-control',{'actualPoint':[[x,y] for x,y in point.items()],'baselineKernel':reducedKernel,'corruptKernel':corrupt,'nonzeroProfileFactor':profile,'actualBaseline':profile*baseline,'actualCorrupt':profile*changed,'movementPerNonzeroProfile':move})
        controls[name]=exact_nonzero_number(J,name+'-movement',move)
        J.zero(name+'-actual-response-factor',profile*changed-profile*baseline,profile*move)
    noTraceContact=phs.subs({m:l,qm:qo},simultaneous=True)
    J.zero('omit-trace-contact-zero',noTraceContact,sp.S.Zero)
    exact_nonzero_number(J,'surviving-contact-coefficient',C.subs(point,simultaneous=True))
    physicalFirst=R(qo)*(sp.I*mu*(qo-qi)/qo)*R(qi)
    exact_nonzero_number(J,'omit-trace-first-height-movement',(physicalFirst+sp.I*mu*qi/(qo+beta)).subs(point,simultaneous=True))
    require(point[qm].is_positive is True and (-point[qm]).is_negative is True,'wrong outgoing sheet responds')
    J.emit('restored-field-index',{'files':used,'completedFunctionsReplayed':False,'newDerivations':['unrestricted defining identities','native generic trace and final slot','first/iterated reference kernels','contact/PV transforms','bounds and four new controls']})
    return {'executionStatus':'COMPLETED_REFERENCE_GRAZING_CERTIFICATES','bothFaces':True,'grades':[[0,0],[1,0],[0,1],[1,1]],'controls':list(controls),'integralsEvaluated':False,'finiteSolves':0,'productionChanges':False,'physicalLossClaim':False,'fullSourceConsumerAssembly':False,'analyticArgument':'Independently assessed limiting proof joined to new exact source/inequality certificates; not machine measure-theory proof.','exclusions':['eta^2','sigma^2','kappa=0','beta=0','full slab operator','finite inverse','old cutoff4/6 transfer','loss','drain','primitive calibration','parameter differentiability']}


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
        result=J.stage('reference-grazing-certificates',{'manifestSha256':gate['manifestSha256'],'methodRecordSha256':gate['methodRecordSha256']},lambda:run_science(manifest,J,ns));code=0
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
