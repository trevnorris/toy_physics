#!/usr/bin/env python3
"""Independent S15 construction; Singular alone computes vanishing conditions.

No S15 engine, export, transcript or run record is an input. Run under the
repository guard. Work products and generated Singular inputs go to the caller's
fresh work directory. No elapsed-time or inactivity cutoff is installed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import traceback

import sympy as s

HERE = Path(__file__).resolve().parent
SINGULAR = Path('/var/projects/_tools/singular/bin/Singular')
CORE = HERE / 'differential_ideal.sing'
d, r, R, rin, L = s.symbols('d r R rin L', positive=True)
a1, a2, b1, b2 = s.symbols('a1 a2 b1 b2', real=True)
A1, A2, B1, B2 = s.symbols('A1 A2 B1 B2', real=True)
At1, At2, Bt1, Bt2, Rt, rint = s.symbols('At1 At2 Bt1 Bt2 Rt rint', real=True)
h = s.symbols('h1:4', real=True)
e = s.symbols('e1:4', real=True)
v1 = s.symbols('v11 v12 v13', real=True)
v2 = s.symbols('v21 v22 v23', real=True)
mu = s.Symbol('mu', real=True)
W = s.symbols('W1:4')
varrho, eta, zeta, pbase, pprime, pwhole = s.symbols('varrho eta zeta pbase pprime pwhole')
general = set(W) | {varrho, eta, zeta, pbase, pprime, pwhole}
H = s.Function('H')  # H(p,q) = integral_(d-q)^(d+q) u^p du
K = s.Function('K')  # K(t,p) = integral_rin^R r^t H(p,r) dr
U = s.Function('U')  # U(p,q,sign) = (d+sign*q)^p
T = s.Function('T')  # T(t) = integral_rin^R r^t dr
state_types = (H, K, U)


def emit(tag, value):
    print(f'{tag}: {value}', flush=True)


def step(name):
    emit('WORKING_ON', name)


def text(expr):
    return s.sstr(expr)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dot(x, y):
    return sum((u*v for u, v in zip(x, y)), s.S.Zero)


def mask_name(mask):
    return '{' + ','.join(n for bit, n in ((1, 'drain1'), (2, 'drain2'), (4, 'remainder')) if mask & bit) + '}'


@lru_cache(None)
def angular_tensor(indices):
    """Azimuthal mean of a Cartesian monomial at fixed e.h=mu."""
    n = len(indices)
    if n == 0:
        return s.S.One
    if n == 1:
        return mu*h[indices[0]]
    if n == 2:
        i, j = indices
        return mu**2*h[i]*h[j] + (1-mu**2)*(int(i == j)-h[i]*h[j])/2
    if n == 3:
        i, j, k = indices
        return (mu**3*h[i]*h[j]*h[k] + mu*(1-mu**2)/2 *
                (h[i]*int(j == k)+h[j]*int(i == k)+h[k]*int(i == j)-3*h[i]*h[j]*h[k]))
    raise ValueError('Cartesian angular rank exceeds the constructed moment rules')


@lru_cache(None)
def axial_moment(power, degree):
    # The rule for a literal polynomial is an integration rule, not an exponent
    # stratum. Symbolic power values (including zero) all use H's entire formula.
    if power == 0:
        return s.Rational(1+(-1)**degree, degree+1)
    return sum((s.binomial(degree, j)*(-1)**j*(d*d+r*r)**(degree-j) *
                H(s.expand(power+2*j+1), r)
                for j in range(degree+1)), s.S.Zero) / ((2*d*r)**degree*d*r)


@lru_cache(None)
def angular(expr):
    """Integrate exactly over the sphere; no parameter strata are selected."""
    result = s.S.Zero
    for term in s.Add.make_args(s.expand(expr)):
        powers = term.as_powers_dict()
        lp = powers.get(L, s.S.Zero)
        indices = tuple(i for i in range(3) for _ in range(int(powers.get(e[i], 0))))
        coeff = term / (L**lp * s.prod(e[i] for i in indices))
        for (deg,), c in s.Poly(s.expand(angular_tensor(indices)), mu).terms():
            result += 2*s.pi*coeff*c*axial_moment(lp, deg)
    return s.expand(result)


def radial(expr):
    """Integrate a finite linear combination of r^t H(p,r), and r^t."""
    result = s.S.Zero
    for term in s.Add.make_args(s.expand(expr)):
        hs = list(term.atoms(H))
        if len(hs) > 1:
            raise ValueError('radial moment has more than one angular primitive')
        f = hs[0] if hs else s.S.One
        rest = s.cancel(term/f)
        exponent = rest.as_powers_dict().get(r, s.S.Zero)
        coeff = s.powsimp(rest/r**exponent, force=True)
        if r in coeff.free_symbols:
            raise ValueError('radial coefficient was not separated')
        result += coeff*(K(s.expand(exponent), f.args[0]) if hs else T(s.expand(exponent)))
    return s.expand(result)


def state_derivative(f):
    if f.func == H:
        p, q = f.args
        return U(p, q, 1)-U(p, q, -1)
    if f.func == U:
        p, q, sign = f.args
        return p*f/(d+sign*q)
    if f.func == K:
        t, p = f.args
        return ((t+p+2)*f-R**(t+1)*H(p, R)+rin**(t+1)*H(p, rin))/d
    raise ValueError(f'Unknown state {f}')


def dd(expr):
    ans = s.diff(expr, d)
    for f in sorted(expr.atoms(*state_types), key=s.default_sort_key):
        ans += s.diff(expr, f)*state_derivative(f)
    return s.expand(ans)


def time_rate(momentum, shell):
    w = [v1[j]-v2[j] for j in range(3)]
    hw = dot(h, w)
    ans = dd(momentum)*hw
    for j in range(3):
        ans += s.diff(momentum, h[j])*(w[j]-h[j]*hw)/d
    for amp, rate in zip((A1,A2,B1,B2), (At1,At2,Bt1,Bt2)):
        ans += s.diff(momentum, amp)*rate
    ans += Rt*shell.subs(r, R)-rint*shell.subs(r, rin)
    return s.expand(ans)


def inventory():
    """Construct all label sets before applying the integrand-based scope rule."""
    rho = [(1, B1*L**b1), (2, B2*r**b2), (4, varrho)]
    vel = [(1, [A1*L**(a1-1)*(r*e[j]-d*h[j]) for j in range(3)]),
           (2, [A2*r**a2*e[j] for j in range(3)]), (4, W)]
    records = []
    excluded = []
    def add(name, expr, mode, radius=None):
        # Scope is fixed by the unintegrated labelled field construction.
        if expr.free_symbols & general:
            excluded.append((name, 'general field in integrand'))
        else:
            records.append(dict(name=name, integrand=expr, mode=mode, radius=radius))
    for surface, radius, sign in [('inner', rin, -1), ('outer', R, 1)]:
        for i in range(3):
            pieces = {p: {m: s.S.Zero for m in range(1,8)}
                      for p in ('carried','centre','radius','viscous','pressure_expansion')}
            for lm, den in rho:
                for vm, v in vel:
                    pieces['centre'][lm|vm] += sign*r*r*den*v[i]*dot(v2,e)
                    pieces['radius'][lm|vm] += sign*r*r*den*v[i]*(rint if surface=='inner' else Rt)
                    for wm, z in vel:
                        pieces['carried'][lm|vm|wm] -= sign*r*r*den*v[i]*dot(z,e)
            # P6 is linear in velocity gradients. Gradients are left symbolic
            # here because eta,zeta make every active such group out of scope.
            for vm, _ in vel:
                g = s.Symbol(f'gradient_{vm}_{i}')
                pieces['viscous'][vm] = sign*r*r*(eta+zeta)*g
            pieces['pressure_expansion'][2|4] = -sign*r*r*pbase*e[i]
            pieces['pressure_expansion'][1|2|4] = -sign*r*r*pprime*rho[0][1]*e[i]
            for piece, groups in pieces.items():
                for mask, expr in groups.items():
                    add(f'{surface}/component_{i+1}/{piece}_labels_{mask_name(mask)}', expr, 'surface', radius)
                add(f'{surface}/component_{i+1}/{piece}_whole', sum(groups.values()), 'surface', radius)
            for piece in ('pressure_whole','pressure_remainder','other','rest','transport','crossing'):
                excluded.append((f'{surface}/component_{i+1}/{piece}', 'general field in integrand'))
    for i in range(3):
        groups = {m:s.S.Zero for m in range(1,8)}
        for lm, den in rho:
            for vm, v in vel:
                groups[lm|vm] += r*r*den*v[i]
        for mask, expr in groups.items():
            add(f'region/component_{i+1}/momentum_labels_{mask_name(mask)}', expr, 'region')
            add(f'region/component_{i+1}/rate_labels_{mask_name(mask)}', expr, 'rate')
        for piece in ('momentum_whole','rate_whole','integrated_force_density','rest','transport','crossing'):
            excluded.append((f'region/component_{i+1}/{piece}', 'general field in integrand'))
    return records, excluded


MATH_PREMISES = {
 'angular': 'Azimuthal tensor means on S^2 are: <e_i>=mu*h_i; <e_i e_j>=mu^2*h_i*h_j+(1-mu^2)*(delta_ij-h_i*h_j)/2; <e_i e_j e_k>=mu^3*h_i*h_j*h_k+mu*(1-mu^2)*(h_i*delta_jk+h_j*delta_ik+h_k*delta_ij-3*h_i*h_j*h_k)/2. dOmega=2*pi*dmu after this mean. Substitution u=|r*e-D| gives mu=(d^2+r^2-u^2)/(2*d*r), -dmu=u*du/(d*r).',
 'radial': 'Euler integral for Appell F1: integral_0^1 (1-x*z)^(-beta)*(1-y*z)^(-betaPrime) dz = AppellF1(1;beta,betaPrime;2;x,y), on the real positive-base continuation used here. Integration by parts and d Li2(z)/dz=-log(1-z)/z give the exceptional logarithmic formulas. All bases stay positive for d>R>rin>0.',
 'differentiation': 'On d>R>rin>0, real powers of the positive distances are smooth and real analytic in d; the compact region permits differentiation under the integrals. Amplitudes and radii are C1 at the instant and exponents are constant in time. These are sufficient existence conditions, not claimed necessary.',
 'system': 'H(p,q)=integral_(d-q)^(d+q) u^p du; U(p,q,sign)=(d+sign*q)^p; K(t,p)=integral_rin^R r^t H(p,r) dr; T(t)=integral_rin^R r^t dr. Leibniz differentiation gives H_d=U_plus-U_minus and U_d=p*U/(d+sign*q). Homogeneity of H and integration by parts give d*K_d=(t+p+2)*K-R^(t+1)*H(p,R)+rin^(t+1)*H(p,rin). These identities hold at all real p,t.',
 'scaling': 'For fixed h=D/d, the geometric departure is d*partial_d Q+2*Q. On a connected d-domain its vanishing integrates to Q(d,h)=(2*R/d)^2*Q(2*R,h).'
}


def primitive_definitions():
    emit('FORMING_CASE_SPLIT', 'H(p,q): p+1 != 0 or p+1 = 0; these cases form the primitive, not a vanishing test')
    emit('DEFINITION', 'H(p,q) = Piecewise(((d+q)^(p+1)-(d-q)^(p+1))/(p+1) if p+1 != 0; log((d+q)/(d-q)) if p+1 = 0)')
    emit('FORMING_CASE_SPLIT', 'T(t): t+1 != 0 or t+1 = 0; primitive formation')
    emit('DEFINITION', 'T(t) = Piecewise((R^(t+1)-rin^(t+1))/(t+1) if t+1 != 0; log(R/rin) if t+1 = 0)')
    emit('DEFINITION', 'J(t,q,sign) = (R-rin)*rin^t*(d+sign*rin)^q*AppellF1(1;-t,-q;2;-(R-rin)/rin,-sign*(R-rin)/(d+sign*rin))')
    emit('FORMING_CASE_SPLIT', 'K(t,p): p+1 != 0; p+1 = 0 and t+1 != 0; p+1 = 0 and t+1 = 0; primitive formation')
    emit('DEFINITION', 'K(t,p) = Piecewise((J(t,p+1,1)-J(t,p+1,-1))/(p+1) if p+1 != 0; (R^(t+1)*log((d+R)/(d-R))-rin^(t+1)*log((d+rin)/(d-rin))-J(t+1,-1,1)-J(t+1,-1,-1))/(t+1) if p+1 = 0 and t+1 != 0; Li2(R/d)-Li2(-R/d)-Li2(rin/d)+Li2(-rin/d) if p+1 = 0 and t+1 = 0)')
    emit('DEFINITION', 'U(p,q,sign)=(d+sign*q)^p; h_i=D_i/d; vki=dot_X_k^i; Atk=dot_A_k; Btk=dot_B_k; Rt=dot_R; rint=dot_r_in; rin=r_in')


def explicit(f, at):
    """Render all primitives in named quantities, with their forming strata."""
    def tt(t):
        return f'Piecewise((R**({text(t+1)})-rin**({text(t+1)}))/({text(t+1)}) if ({text(t+1)})!=0; log(R/rin) otherwise)'
    def hh(p,q):
        return f'Piecewise((({at}+{q})**({text(p+1)})-({at}-{q})**({text(p+1)}))/({text(p+1)}) if ({text(p+1)})!=0; log(({at}+{q})/({at}-{q})) otherwise)'
    def jj(t,q,sg):
        return f'((R-rin)*rin**({text(t)})*({at}+({sg})*rin)**({text(q)})*AppellF1(1,-({text(t)}),-({text(q)}),2,-(R-rin)/rin,-({sg})*(R-rin)/({at}+({sg})*rin)))'
    if f.func == H:
        return hh(*f.args)
    if f.func == U:
        p,q,sg=f.args
        return f'({at}+({sg})*({q}))**({text(p)})'
    if f.func == T:
        return tt(f.args[0])
    if f.func == K:
        t,p=f.args
        return (f'Piecewise(({jj(t,p+1,1)}-{jj(t,p+1,-1)})/({text(p+1)}) if ({text(p+1)})!=0; '
                f'(R**({text(t+1)})*log(({at}+R)/({at}-R))-rin**({text(t+1)})*log(({at}+rin)/({at}-rin))-{jj(t+1,-1,1)}-{jj(t+1,-1,-1)})/({text(t+1)}) if ({text(p+1)})=0 and ({text(t+1)})!=0; '
                f'Li2(R/({at}))-Li2(-R/({at}))-Li2(rin/({at}))+Li2(-rin/({at})) otherwise)')
    raise ValueError(f)


def render(expr, at='d'):
    subs={}
    replacements={}
    for i,f in enumerate(sorted(expr.atoms(H,K,U,T),key=s.default_sort_key)):
        key=f'FUNCTIONTOKEN{i}X'
        subs[f]=s.Symbol(key)
        replacements[key]='('+explicit(f,at)+')'
    st=text(expr.xreplace(subs))
    st=re.sub(r'\bd\b', '('+at+')', st)
    for key,val in replacements.items():
        st=st.replace(key,val)
    return st


def polynomial_system(q, departure):
    """Encode only algebraic formation; the vanishing decisions are in Singular."""
    states=set(q.atoms(*state_types)) | set(departure.atoms(*state_types))
    while True:
        more=set().union(*(state_derivative(f).atoms(*state_types) for f in states)) if states else set()
        if more <= states:
            break
        states |= more
    states=sorted(states,key=s.default_sort_key)
    ys={f:s.Symbol(f'y{i}') for i,f in enumerate(states)}
    derivs=[state_derivative(f).xreplace(ys) for f in states]
    f=departure.xreplace(ys)
    # Lift real powers of held radii and the radial elementary primitive to
    # constants with explicit defining expressions. No zero test is performed.
    constants=set()
    for expr in [f,*derivs]:
        constants |= expr.atoms(T)
        for pw in expr.atoms(s.Pow):
            if pw.exp.free_symbols:
                if d in pw.free_symbols:
                    raise ValueError('unencoded separation power')
                constants.add(pw)
    constants=sorted(constants,key=s.default_sort_key)
    cs={v:s.Symbol(f'c{i}') for i,v in enumerate(constants)}
    f=s.cancel(f.xreplace(cs)/s.pi)
    # Every Q has a common pi factor from azimuthal integration. This division
    # is a fixed nonzero scalar, not a parameter division.
    num,den=s.fraction(f)
    ddeltas=[s.cancel(x.xreplace(cs)) for x in derivs]
    P=d*(d-R)*(d+R)*(d-rin)*(d+rin)
    delta=[s.cancel(P*x) for x in ddeltas]
    if any(s.denom(x)!=1 for x in delta):
        raise ValueError('nonpolynomial state derivation')
    values={str(y):explicit(st,'2*R') for st,y in ys.items()}
    for st,c in cs.items():
        values[str(c)]=render(st,'2*R')
    # Denominator clearing is permissible only for a nonzero geometric factor.
    allowed={d,R,rin,d-R,d+R,d-rin,d+rin}
    for base,power in s.factor_list(den)[1]:
        if base not in allowed:
            raise ValueError(f'parameter denominator in formed departure: {base}')
    symbols=set(num.free_symbols) | set(P.free_symbols) | set().union(*(x.free_symbols for x in delta))
    symbols |= set(h) | {d,R,rin} | set(ys.values())
    # h3 is first so the sphere remainder has degree <2 in h3.
    ordered=[h[2],h[1],h[0]]+sorted(symbols-set(h),key=str)
    def singular(x):
        return text(s.expand(x)).replace('**','^')
    source=('ring computation=0,('+','.join(map(str,ordered))+'),dp;\n'
            'option(redSB);\nshort=0;\n'
            f'poly F={singular(num)};\npoly P={singular(P)};\n'
            'ideal stateVars='+(','.join(map(str,ys.values())) if ys else '0')+';\n'
            'ideal stateDelta='+(','.join(singular(x) for x in delta) if delta else '0')+';\n'
            'print("SINGULAR_VERSION: "+string(system("version")));\n'
            +CORE.read_text()+'\nquit;\n')
    return source,values,num,den,P,dict(zip(map(str,ys.values()),map(text,states)))


def condition_render(poly, values):
    return re.sub(r'\b(?:y|c)\d+\b',lambda m:'('+values[m.group(0)]+')',poly)


def compute(record, index, work):
    name=record['name']
    step(name+'/item3/formation')
    emit('DOMAIN', 'd>R>rin>0; a1,a2,b1,b2 real; real amplitudes and velocities; spheres about X2; outward sign inner=-1, outer=+1; both cores excluded')
    emit('EXISTENCE_CONDITION', MATH_PREMISES['differentiation'])
    emit('PREMISE', MATH_PREMISES['angular'])
    emit('FORMING_RULE', 'literal L^0 terms use polynomial mu integration; every symbolic distance exponent uses H with all forming branches retained')
    emit('INTEGRAND', text(record['integrand']))
    shell=angular(record['integrand'])
    emit('ANGULAR_OBJECT', text(shell))
    if record['mode']=='surface':
        q=shell.subs(r,record['radius'])
    else:
        emit('PREMISE', MATH_PREMISES['radial'])
        q=radial(shell)
        if record['mode']=='rate':
            emit('PREMISE', MATH_PREMISES['system'])
            q=time_rate(q,shell)
    emit('PREMISE', MATH_PREMISES['system'])
    emit('PREMISE', MATH_PREMISES['scaling'])
    departure=s.expand(d*dd(q)+2*q)
    primitive_definitions()
    emit('Q',text(q))
    emit('ITEM3',text(departure))
    emit('ITEM3_EXPLICIT',render(departure))
    step(name+'/item4/encoding')
    source,values,num,den,P,states=polynomial_system(q,departure)
    emit('ENCODED_NUMERATOR',text(num))
    emit('CLEARED_DENOMINATOR',text(den))
    emit('DERIVATION_MULTIPLIER',text(P))
    for key,val in states.items():
        emit('STATE',f'{key}={val}')
    for key,val in values.items():
        emit('AT_REFERENCE',f'{key}={val}')
    infile=work/f'object_{index:03d}.sing'
    rawfile=work/f'object_{index:03d}.stdout'
    infile.write_text(source)
    emit('SINGULAR_INPUT_SHA256',sha(infile))
    step(name+'/item4/Singular')
    command=[str(SINGULAR),'--no-stdlib','--no-rc','--no-shell','-q',str(infile)]
    lines=[]
    with rawfile.open('w') as raw:
        child=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in child.stdout:
            raw.write(line); raw.flush()
            print(line,end='',flush=True)
            lines.append(line.rstrip('\n'))
        code=child.wait()
    errors=[line for line in lines if line.lstrip().startswith('?') or 'error occurred' in line]
    if code or errors or 'SINGULAR_COMPLETED: 1' not in lines:
        raise RuntimeError(f'Singular code={code}; diagnostics={errors}')
    conditions=[line.partition(': ')[2] for line in lines if line.startswith('CONDITION_POLYNOMIAL: ')]
    step(name+'/item4/printed_conditions')
    emit('CONDITION_DOMAIN','R>rin>0; a1,a2,b1,b2 real; all amplitudes, their time rates, radii rates and velocities real; sufficient existence condition as printed')
    emit('CONDITION_CONNECTIVE','AND; an empty equation list is True; every Piecewise forming case is included')
    for poly in conditions:
        emit('CONDITION_EQUATION',condition_render(poly,values)+' = 0')
    emit('Q_UNDER_CONDITIONS','(2*R/d)**2*('+render(q,'2*R')+')')
    emit('OBJECT_END',json.dumps({'name':name,'item3':'exact','item4':'conditions_returned','condition_count':len(conditions)},sort_keys=True))
    return {'name':name,'item3':'exact','item4':'conditions_returned','condition_count':len(conditions),
            'input':infile.name,'input_sha256':sha(infile),'raw_output_sha256':sha(rawfile)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-directory',required=True,type=Path)
    args=parser.parse_args()
    if not os.environ.get('S11C_POOLED_GUARD_MANIFEST'):
        raise RuntimeError('A pooled guard invocation is required')
    args.work_directory.mkdir(parents=True,exist_ok=False)
    step('supplied equations and inventory')
    emit('SUPPLIED','rho_br=rho1+rho2+varrho; V=V1+V2+W; G=rho_br*V; carried=-integral rho_br*V_i*(V.n); centre=integral G_i*(dot_X2.n); radius=integral G_i*radius_rate*sign; P6 stress and pressure expansion as specified in the decision list')
    emit('SUPPLIED','rho_k=B_k*r_k^b_k; V_k=A_k*r_k^(a_k-1)*(x-X_k); D=X1-X2; d=|D|; x=X2+r*e; Egeo=d*partial_d Q+2*Q at fixed h=D/d')
    records,excluded=inventory()
    for rec in records:
        emit('IN_SCOPE',rec['name'])
    for name,reason in excluded:
        emit('OUT_OF_SCOPE',name+'; '+reason)
    (args.work_directory/'inventory.json').write_text(json.dumps([r['name'] for r in records],indent=2))
    statuses=[]
    error=None
    try:
        for index,record in enumerate(records):
            statuses.append(compute(record,index,args.work_directory))
    except Exception as exc:
        error=repr(exc)
        emit('ATTEMPT_ERROR',error)
        traceback.print_exc()
        for record in records[len(statuses):]:
            emit('OBJECT_END',json.dumps({'name':record['name'],'item3':'not_obtained_in_this_attempt','item4':'not_obtained','reason':'code_error','attempted': 'formation and differential-ideal construction' if record is records[len(statuses)] else 'not reached after code error'},sort_keys=True))
    finally:
        inputs={str(Path(__file__).resolve()):sha(__file__),str(CORE):sha(CORE),str(SINGULAR):sha(SINGULAR)}
        for module in list(sys.modules.values()):
            file=getattr(module,'__file__',None)
            if file and Path(file).is_file():
                inputs[str(Path(file).resolve())]=sha(file)
        (args.work_directory/'run_data.json').write_text(json.dumps({'statuses':statuses,'error':error,'inputs':inputs,'end_utc':datetime.now(timezone.utc).isoformat()},indent=2))
    emit('ENGINE_END',json.dumps({'objects_in_scope':len(records),'objects_completed':len(statuses),'error':error},sort_keys=True))
    return int(error is not None)


if __name__=='__main__':
    sys.exit(main())
