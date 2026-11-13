from collections import Counter
import re
import sympy as s
from packet_io import *
import S15_reading_math as m
import S15_reading_comparator as c
d=Records()
cat=member(d.tar,'context_checks.json.gz')
args=cat['argument_slots'];time,ld=s.symbols('T L')
wl=(S.parent/'mathematica/S15_separation_departure.wl').read_text()
decl=re.search(r'Do\[register\[centreVelocityChange\[j\],(.+?),\{([a-zA-Z]+)\}\],\{j,(\d+)\}\];',wl)
basisdef=re.search(r'(?m)^.*?(timeDim\s*=\s*\{[\d, ]+\}).*?;',wl)
print('native WL',decl[0],'time',basisdef[1])
ids=[d.n['wl_reading',q][0] for q in ('field_argument_dimension_signatures','family/velocity_changes','family/dimension_signatures_including_changes')]
print('native declaration operands',ids,'before registration time-head presence',('centreVelocityChange' in d.raw(ids[0])),'later recorded RHS',d.raw(ids[1]))
for row in args:
    if row['status']=='no_wl_declared_slot':
        j=int(row['argument'][0].split('_')[2])
        print('supplemental dimension match',j,row['argument'],row['expected'],row['observed'],decl[1],decl[2]);assert 1<=j<=int(decl[3]) and decl[2]=='timeDim' and row['expected']==row['observed']==[[0,0,1]]
# Full generating polynomials and their separate native integrals; cosine moments
# integrate (cos phi)^k on the full circle. Test all actual degrees exactly.
a,b,mu,az=s.symbols('A B mu az',real=True); scope=(('axial',a),('transverse',b),('mu',mu),('azimuth',az))
for i,rid in enumerate(d.n['wl_reading','item3_kernel/generator_operand']):
    tree=m.original.wl_parse(d.raw(rid)); _,args0=m.original.application(tree)
    expr=m.translate(args0[0],'wl',scope);bounds=m.translate(args0[1],'wl',scope)
    degree=s.Poly(expr,a,b).total_degree()
    target=(a*mu+b*s.sqrt(1-mu**2)*s.cos(az))**degree
    raw_integral=0
    for (k,),coef in s.Poly(expr,s.cos(az)).terms():
        mean=s.binomial(k,k//2)/2**k if k%2==0 else s.S.Zero
        raw_integral+=2*s.pi*coef*mean
    rid2=d.n['wl_reading','item3_kernel/generator_integral'][i]
    actual=m.translate(m.original.wl_parse(d.raw(rid2)),'wl',scope)
    result=[s.simplify(expr-target),s.simplify(raw_integral-actual)]
    print('full WL generator',rid,rid2,degree,str(bounds),'residuals',result);assert result==[0,0]
# One-sided source slot and branch guard controls, signed operand, empty sets.
from io import BytesIO
for a,b,aa,bb,item in [('[ (1, 4, -2) ]','{{1, 4, -2}}','[(1,4,-2)]','{{1,-4,-2}}',6),('[]','{}','[]','{{0,0,0}}',6)]:
    for tag,x,y in [('actual',a,b),('changed',aa,bb)]:
        p=x.encode()+b'\n';w=y.encode()+b'\n'
        ids=['py_reading:0','wl_reading:0'];c.LOCATOR={ids[0]:('py_reading',{'offset':0,'length':len(p)},sha(p[:-1])),ids[1]:('wl_reading',{'offset':0,'length':len(w)},sha(w[:-1]))}
        cmp=c.ObjectComparison({'py_reading':BytesIO(p),'wl_reading':BytesIO(w)})
        cmp.compare('E_geo/dimension_operands',item,ids[:1],ids[1:],'cross')
        print('dimension one-sided',tag,x,y,cmp.rows[0]['status'],cmp.rows[0].get('residual'))
