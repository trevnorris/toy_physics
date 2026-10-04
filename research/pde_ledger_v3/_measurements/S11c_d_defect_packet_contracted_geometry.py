"""New one-dimensional contraction mesh. Native calls require containment.

Only the old generic Quad/Line arithmetic is used; no old geometry constructor
is called. Synthetic tests use unrelated fields, windows and carriers.
"""
from fractions import Fraction as F
from itertools import combinations
from S11c_d_defect_packet_geometry_lib import Quad, Line, canonical_lines


def require(v,message):
    if v is not True:raise ValueError(message)


def offsets(scale):return [F(0)]+[F(sign*n,scale) for n in (1,2,4,8) for sign in (-1,1)]


def true_window(K,T,m):return max(-K,m-T),min(K,m+T)


def branches(kappa,K,T,number):
    zero=kappa*0
    low=Line(F(1),zero-T,('clip-low',)) if number==5 else Line(F(0),zero-K,('clip-low',))
    high=Line(F(1),zero+T,('clip-high',)) if number==0 else Line(F(0),zero+K,('clip-high',))
    dg=[(-1,-kappa),(-1,-kappa),(1,kappa),(-1,kappa),(1,-kappa),(1,-kappa)][number]
    return low,high,dg


def graphs(kappa,K,T,carrier,width,length,number):
    zero=kappa*0;lo,hi,(sg,b)=branches(kappa,K,T,number)
    lines=[lo,hi,Line(F(0),zero-K,('box-low',)),Line(F(0),zero+K,('box-high',))]
    for sign in (-1,1):lines.append(Line(F(0),sign*kappa,('root:'+str(sign),)))
    for d in offsets(width):lines.append(Line(F(0),carrier+d,('carrier:'+str(d),)))
    for d in offsets(length):lines.append(Line(F(1),zero+d,('profile:'+str(d),)))
    for root in (-1,1):
        for direction in (-1,1):
            lines.append(Line(F(direction*sg),root*kappa+direction*b,('grazing:'+str(root)+':'+str(direction),)))
    return canonical_lines(lines)


def mesh(kappa,K,T,carrier,width,length):
    require(isinstance(kappa,Quad) and isinstance(carrier,Quad),'exact field')
    require(type(K) is int and type(T) is int and type(width) is int and type(length) is int,'integer radii/scales')
    require(T-K>kappa>0 and width>0 and length>0,'six affinity slabs')
    M=K+T;zero=kappa*0
    initial=[zero-M,zero-(T-K),-kappa,zero,kappa,zero+(T-K),zero+M]
    cuts={v:['affinity:'+str(i)] for i,v in enumerate(initial)};records=[]
    for d in offsets(width):cuts.setdefault(carrier+d,[]).append('outer-carrier:'+str(d))
    for d in offsets(length):cuts.setdefault(carrier+d,[]).append('outer-profile:'+str(d))
    for i,(left,right) in enumerate(zip(initial,initial[1:])):
        lines=graphs(kappa,K,T,carrier,width,length,i);intersections=[]
        for a,b in combinations(range(len(lines)),2):
            if lines[a].slope==lines[b].slope:continue
            x=(lines[b].intercept-lines[a].intercept)/(lines[a].slope-lines[b].slope)
            inside=left<=x<=right
            intersections.append({'lineIds':[a,b],'m':x,'insideAffinitySlab':inside})
            if inside:cuts.setdefault(x,[]).append('intersection:'+str(i)+':'+str(a)+':'+str(b))
        records.append({'index':i,'left':left,'right':right,'graphs':lines,'intersections':intersections})
    ordered=sorted(x for x in cuts if -M<=x<=M);slabs=[]
    for left,right in zip(ordered,ordered[1:]):
        witness=(left+right)/2;i=next(j for j in range(6) if initial[j]<witness<initial[j+1]);lines=records[i]['graphs']
        low,high,_=branches(kappa,K,T,i)
        box=[line for line in lines if -K<=line.at(witness)<=K]
        box.sort(key=lambda line:line.at(witness))
        # Coincident graphs have already been merged with all labels.
        require(len({line.at(witness) for line in box})==len(box),'distinct graphs inside intersection-free slab')
        cells=[]
        for a,b in zip(box,box[1:]):
            midpoint=(a.at(witness)+b.at(witness))/2
            cells.append({'low':a,'high':b,'clipped':low.at(witness)<midpoint<high.at(witness)})
        slabs.append({'left':left,'right':right,'affinity':i,'clipLow':low,'clipHigh':high,'cells':cells})
    result={'K':K,'T':T,'M':M,'kappa':kappa,'carrier':carrier,'width':width,'length':length,'affinity':records,'cuts':[{'m':x,'labels':cuts[x]} for x in ordered],'slabs':slabs}
    audit(result)
    return result


def audit(plan):
    K,T,M=plan['K'],plan['T'],plan['M'];kap=plan['kappa'];slabs=plan['slabs']
    require(len(plan['affinity'])==6 and M==K+T,'six full original affinity slabs')
    require(slabs[0]['left']==-M and slabs[-1]['right']==M,'full wings')
    for i,s in enumerate(slabs):
        left,right=s['left'],s['right'];require(left<right,'strict positive outer slab')
        if i:require(slabs[i-1]['right']==left,'contiguous outer coverage')
        points=[left,(left+right)/2,right];low,high=s['clipLow'],s['clipHigh']
        for m in points:
            actual=true_window(K,T,m)
            require((low.at(m),high.at(m))==actual,'actual max/min clipping window')
            require(-K<=actual[0]<=actual[1]<=K,'window inside original box')
        cells=s['cells'];require(bool(cells),'nonempty full box cells')
        for m in points:
            require(cells[0]['low'].at(m)==-K and cells[-1]['high'].at(m)==K,'complete inner box')
            for j,c in enumerate(cells):
                require(c['low'].at(m)<=c['high'].at(m),'nonnegative affine inner width')
                if j:require(cells[j-1]['high'].at(m)==c['low'].at(m),'contiguous inner coverage')
        mid=points[1];selected=[c for c in cells if c['clipped']]
        require(bool(selected),'positive clipped interval at interior witness')
        require(selected[0]['low'].at(mid)==low.at(mid) and selected[-1]['high'].at(mid)==high.at(mid),'full clipped cell coverage')
        for c in cells:
            val=(c['low'].at(mid)+c['high'].at(mid))/2
            require(c['clipped']==bool(low.at(mid)<val<high.at(mid)),'actual clipped cell membership')
            require(c['low'].at(mid)<c['high'].at(mid),'strict interior width')
        # No omitted crossing is permitted in an open slab. Check all original
        # graphs independently of the planner's selected cells or intersections.
        lines=plan['affinity'][s['affinity']]['graphs']
        for a,b in combinations(lines,2):
            if a.slope!=b.slope:
                cross=(b.intercept-a.intercept)/(a.slope-b.slope)
                require(not left<cross<right,'no missing affine intersection')
        required=graphs(kap,K,T,plan['carrier'],plan['width'],plan['length'],s['affinity'])
        require([v.packed() for v in lines]==[v.packed() for v in required],'complete original graph/label set')
    for d in offsets(plan['width'])+offsets(plan['length']):
        x=plan['carrier']+d
        if -M<=x<=M:require(any(v['m']==x for v in plan['cuts']),'every outer resolution cut')
    require([r['m'] for r in plan['cuts']]==[s['left'] for s in slabs]+[slabs[-1]['right']],'exact cuts/slabs correspondence')
    return True


def wing_control(plan,side):
    from copy import deepcopy
    mutant=deepcopy(plan);index=0 if side<0 else len(plan['slabs'])-1
    s=mutant['slabs'][index];zero=plan['kappa']*0
    s['clipLow']=Line(F(0),zero-plan['K'],('mutant-central-low',))
    s['clipHigh']=Line(F(0),zero+plan['K'],('mutant-central-high',))
    try:audit(mutant)
    except ValueError as error:
        require(str(error)=='actual max/min clipping window','same-path actual window refusal')
        return {'side':side,'original':plan['slabs'][index],'mutant':s,'refusal':str(error)}
    raise ValueError('wing corruption did not refuse')


def packed(value):
    if hasattr(value,'packed'):return value.packed()
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {k:packed(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [packed(v) for v in value]
    return value
