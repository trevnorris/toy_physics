#!/usr/bin/env python3
"""Real-axis radiation rules and source bindings for the authorized finite pilot.

Imported only by a guarded worker. No top-level scientific imports/constructors.
"""
import ast
import hashlib
import math
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace


def stable_panels(points,protected):
    """Coalesce floating aliases, retaining physical endpoints and coverage."""
    result=[];protected=set(protected)
    for value in sorted(set(points)):
        if result and value-result[-1]<=8*max(math.ulp(value),math.ulp(result[-1])):
            if value in protected:
                require(result[-1] not in protected,'distinct protected endpoints must not collapse')
                result[-1]=value
        else:result.append(value)
    require(set(protected)<=set(result),'all exact physical endpoints retained')
    require(result[0]==min(points) and result[-1]==max(points),'unchanged complete physical interval')
    return result


def require(ok,message,evidence=None):
    if not ok:
        error=ValueError(message);error.evidence=evidence;raise error


def native_helpers(engine_path,finite_path,np_,sp_):
    global np,sp,AppliedUndef
    np,sp=np_,sp_
    from sympy.core.function import AppliedUndef
    env={'np':np,'sp':sp,'AppliedUndef':AppliedUndef,'lru_cache':lru_cache}
    proofs={}
    def take(path,names):
        tree=ast.parse(Path(path).read_text());nodes=[]
        for name in names:
            node=next(n for n in tree.body if getattr(n,'name',None)==name)
            proofs[name]=hashlib.sha256(ast.dump(node).encode()).hexdigest();nodes.append(node)
        exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),env)
    take(engine_path,('memo_xreplace','dag_free_symbols','BoundedActionQuadrature','BoundedSourceFourierQuadrature'))
    env['engine']=SimpleNamespace(BoundedSourceFourierQuadrature=env['BoundedSourceFourierQuadrature'])
    env['require']=require
    take(finite_path,('polynomial_basis','BasisMomentum'))
    return SimpleNamespace(**{k:env[k] for k in ('memo_xreplace','dag_free_symbols','BoundedActionQuadrature','BoundedSourceFourierQuadrature','polynomial_basis','BasisMomentum')},proofs=proofs)


def context(native,frequency):
    values=[v['liveFrequencyAndGrades'] for v in frequency['sources']['records'].values()]
    symbols=set().union(*(v.atoms(sp.Symbol) for v in values))
    symbols.update(s for row in native['rows'] for limit in row['limits'] for s in limit[0].atoms(sp.Symbol))
    symbols.update(native['sources'][0,0]['symbolicAmplitude'].atoms(sp.Symbol))
    def find(name):
        found=[s for s in symbols if s.name==name];require(len(found)==1,'unique saved coordinate '+name);return found[0]
    return SimpleNamespace(z=find('s11cdNormalPosition'),zp=find('s11cdSourceNormalPosition'),xi=find('s11cdProfileCoordinate'),regulator=find('s11cdAbelRegulator'))


def source_jets(expression,old,r):
    probes=list(expression.atoms(AppliedUndef));require(len(probes)==1 and probes[0].args==(r.zp,),'single saved field probe')
    probe=probes[0];require(probe.func.__name__.startswith('s11cdPencilProbe'),'source probe identity')
    column=int(probe.func.__name__.removeprefix('s11cdPencilProbe'))
    derivative=expression.atoms(sp.Derivative)
    require(all(v.expr==probe and all(x==r.zp for x,_ in v.variable_count) for v in derivative),'source derivatives only')
    order=max([0]+[sum(n for _,n in v.variable_count) for v in derivative])
    jet_nodes=[sp.Derivative(probe,(r.zp,n)) if n else probe for n in range(order+1)]
    carriers=sp.symbols('numericalSourceJet0:'+str(order+1));formal=expression.xreplace(dict(zip(jet_nodes,carriers)))
    # Coefficients are opaque to the carrier derivative. No general Poly domain.
    coefficients=[sp.diff(formal,c) for c in carriers]
    residual=sp.expand(formal-sum(a*c for a,c in zip(coefficients,carriers)))
    require(residual==0 and all(not a.has(*carriers,AppliedUndef,sp.Derivative,sp.Subs,sp.Integral) and a.free_symbols<={r.zp} for a in coefficients),'complete source-linear reconstruction',{'formal':formal,'coefficients':coefficients,'residual':residual})
    return {'column':column,'coefficients':coefficients,'residual':residual,'originalBoundAmplitude':expression,'probe':probe,'amplitudeUnit':old['amplitudeUnit'],'integralUnit':old['integralUnit']}


def bind_sources(packet,native,grades,contrast,r,H):
    src=packet['sources'];subs={src['frequency']:sp.Integer(3),**{s:sp.Rational(str(contrast))*v for s,v in src['origin'].items()}}
    indexed={tuple(v['address']):v for v in src['records'].values()};grade_index={tuple(v['address']):v['record'] for v in grades['records'].values()}
    require(set(indexed)==set(grade_index) and len(indexed)==375,'complete accepted source-address inventory')
    actual={};joins=[]
    for address,record in indexed.items():
        require(record['original']==grade_index[address]['ORIGINAL'],'original source address join '+str(address))
        value=H.memo_xreplace(record['liveFrequencyAndGrades'],subs)
        require(not value.has(*subs),'complete frequency and contrast binding')
        actual[address]=value;joins.append({'address':address,'source':record['liveFrequencyAndGrades'],'binding':subs,'value':value,'unit':record['unit']})
    jets={};sources={}
    for i in range(35):
        old=native['sources'][0,i];record=indexed['source',i]
        require(record['original']==old['symbolicAmplitude'],'native source amplitude identity')
        jets[i]=source_jets(actual['source',i],old,r)
        frequency=H.memo_xreplace(src['sourceFrequencies'][i]['frequency'],subs)
        require(not src['sourceFrequencies'][i]['dependsOnFrequency'] and frequency==old['frequency'],'saved Fourier source character')
        sources[0,i]=dict(old,frequency=frequency)
    rows=[]
    for old in native['rows']:
        factors=[]
        for j,factor in enumerate(old['factors']):
            address=('factor',old['index'],j);require(indexed[address]['original']==factor['symbolicCoefficient'],'native factor identity')
            factors.append(dict(factor,coefficient=actual[address]))
        require(len({jets[v['sourceIndex']]['column'] for v in factors})==1,'single row field column')
        rows.append(dict(old,factors=factors))
    terms=grades['termJoins'];require(len(terms)==160,'all native cell occurrences')
    for term in terms:
        old=native['rows'][term['integralIndex']]
        require(old['original']==term['originalIntegral'],'native integral occurrence join')
        require(tuple(tuple(v.xreplace(native['cutoffBindings']) for v in lim) for lim in term['remainingLimits'])==old['limits'],'native ordered limits')
        require(indexed['cell',term['row'],term['column'],term['term']]['original']==term['cellCoefficient'],'original cell coefficient identity')
        require(all(jets[f['sourceIndex']]['column']==term['column'] for f in old['factors']),'native cell field column')
    cells=[]
    for cell in native['cells']:
        if cell['test']!=0:continue
        cells.append(dict(cell,terms=[(row,actual['cell',cell['row'],cell['column'],j]) for j,(row,_) in enumerate(cell['terms'])]))
    orders=sorted({a[1] for a in actual if a[0]=='local'})
    local={n:sp.ImmutableMatrix(5,5,[actual['local',n,i,j] for i in range(5) for j in range(5)]) for n in orders}
    profiles=set().union(*(f['coefficient'].atoms(sp.Integral) for row in rows for f in row['factors']))
    require(profiles==set(native['profileUnits']),'same complete finite profile operands')
    return {'contrast':contrast,'rows':rows,'sources':sources,'jets':jets,'cells':cells,'local':local,'actual':actual,'sourceJoins':joins,'termJoins':terms,'pairs':native['pairs'],'abel':native['abel'],'profileUnits':native['profileUnits'],'fieldUnits':packet['fieldUnits'],'rowUnits':packet['rowUnits']}


def radical_inventory(src,J):
    w=src['frequency'];bases=sorted({p.base for v in src['records'].values() for p in v['census']['fractionalPowers']},key=sp.default_sort_key)
    roots={}
    for i,base in enumerate(bases):
        momenta=base.free_symbols-{w};require(len(momenta)==1,'one momentum per actual radical');p=next(iter(momenta))
        bound=sp.expand(base.subs(w,3));poly=sp.Poly(bound,p,domain='EX');a=poly.coeff_monomial(p**2);c=poly.coeff_monomial(1)
        require(sp.expand(bound-a*p*p-c)==0 and a<0 and c>0,'radiating quadratic radical')
        b2=sp.cancel(-c/a);require(b2==sp.Rational(1,25),'actual omega3 threshold joins source b squared')
        roots[p]={'original':base,'bound':bound,'coefficient':a,'constant':c,'bSquared':b2,'b':float(sp.sqrt(b2)),'scale':sp.sqrt(-a),'carrier':sp.Symbol('radiatingQ'+str(i),complex=True)}
    require(len(roots)==3,'three actual momentum legs');J.report('radical-inventory',{str(p):{k:str(v) for k,v in x.items()} for p,x in roots.items()});return roots


def lift_roots(expression,roots):
    bybase={r['bound']:r['carrier'] for r in roots.values()}
    @lru_cache(None)
    def walk(v):
        if v in bybase:return bybase[v]**2
        if v.is_Pow and v.base in bybase and v.exp.is_Rational:
            require((2*v.exp).is_Integer,'half-integer native radical power');return bybase[v.base]**(2*v.exp)
        if not v.args:return v
        args=tuple(walk(a) for a in v.args);return v if args==v.args else v.func(*args)
    return walk(expression)


def endpoint_orders(binding,roots,parameters,J):
    # The positive-real/positive-imaginary q axes cannot meet -C: Re(C),Im(C)>0.
    params={k:sp.Rational(v) for k,v in parameters.items()};scales={v['scale'] for v in roots.values()};require(len(scales)==1,'common actual radical scaling');C=next(iter(scales))*params['Lambda_A_0']*3/(params['rho_m']*(1-sp.I*3*params['tau_A']));require(sp.re(C)>0 and sp.im(C)>0,'physical memory/permeability coupled constant')
    carriers=[r['carrier'] for r in roots.values()];proofs=[];orders=[]
    def valuation(v,q):
        if not v.has(q):return 0
        if v==q:return 1
        if v.is_Add:return min(valuation(a,q) for a in v.args)
        if v.is_Mul:return sum(valuation(a,q) for a in v.args)
        if v.is_Pow and v.exp.is_Integer:
            if v.exp>=0:return int(v.exp)*valuation(v.base,q)
            matches=[]
            for family,degree in ((q,1),(q*q,2),(q+C,0),(1+C/q,-1),((q+C)**2,0),((1+C/q)**2,-2)):
                ratio=sp.cancel(v.base/family)
                if not ratio.free_symbols and ratio!=0:matches.append((family,degree,ratio))
            require(matches,'unsupported q-dependent denominator',{'base':v.base,'q':q})
            family,degree,ratio=matches[0];proofs.append({'base':v.base,'family':family,'constant':ratio,'residual':sp.cancel(v.base-ratio*family),'carrier':q})
            return int(v.exp)*degree
        raise ValueError('unsupported endpoint dependence '+str(v.func))
    for row in binding['rows']:
        for fi,f in enumerate(row['factors']):
            lifted=lift_roots(f['coefficient'],roots)
            values={str(q):valuation(lifted,q) for q in carriers}
            record={'row':row['index'],'factor':fi,'actual':f['coefficient'],'lifted':lifted,'orders':values};J.blob(f'endpoint-orders/row-{row["index"]}-factor-{fi}.pickle',record)
            require(min(values.values())>=-1,'nonintegrable or unsupported uncancelled endpoint order',record);orders.append({'row':row['index'],'factor':fi,'orders':values})
    result={'coupledConstant':C,'positiveReal':sp.re(C),'positiveImaginary':sp.im(C),'denominatorJoins':proofs,'orders':orders,'scope':'Sufficient source-specific endpoint power bounds; smooth finite profile/source factors have no radical-carrier dependence.'}
    J.report('endpoint-order-summary',{'coupledConstant':str(C),'rows':len(orders),'minimumOrder':min(v for r in orders for v in r['orders'].values()),'denominatorJoins':len(proofs)})
    return result


def radiating_class(H):
    class Radiating(H.BasisMomentum):
        def __init__(self,binding,r,roots,setting):
            super().__init__(binding['rows'],binding['sources'],r,workspace_bytes=32*1024**2,batch_nodes=128)
            self.binding=binding;self.roots=roots;self.setting=setting;self.omit_jacobian=False;self.wrong_sheet=None
        def rule(self,lower,upper,order,centers=(),width=None):
            b=next(iter(self.roots.values()))['b'];require(lower==-upper and upper>b,'finite symmetric radiating rule')
            points={lower,upper,-b,b}
            for center in centers:
                require(np.isfinite(center) and width>0,'finite transfer panel center')
                if lower<center<upper:points.add(float(center))
                distance=width
                while distance<2*(upper-lower+abs(center)):
                    points.update(v for v in (center-distance,center+distance) if lower<v<upper);distance*=2
            if order not in self._legendre:self._legendre[order]=np.polynomial.legendre.leggauss(order)
            x,w=self._legendre[order];nodes=[];weights=[];panels=stable_panels(points,(lower,-b,b,upper))
            for lo,hi in zip(panels[:-1],panels[1:]):
                if lo>=-b and hi<=b:
                    a,c=np.arcsin(np.clip([lo/b,hi/b],-1,1));t=(a+c)/2+(c-a)*x/2;k=b*np.sin(t);jac=b*np.cos(t)
                else:
                    sign=1 if lo>=b else -1;a,c=sorted(np.arccosh(np.maximum(1,np.abs([lo,hi])/b)));t=(a+c)/2+(c-a)*x/2;k=sign*b*np.cosh(t);jac=b*np.sinh(t)
                nodes.append(k);weights.append((c-a)*w/2*(1 if self.omit_jacobian else jac))
            nodes=np.concatenate(nodes);weights=np.concatenate(weights)
            require(np.isfinite(nodes).all() and np.all(weights>0),'finite oriented transformed measure')
            return nodes,weights,panels
        def prepare_basis(self,jets,nodes,bound,size,source_order):
            self.jet_data=jets;self.source_nodes,self.source_weights=H.BoundedSourceFourierQuadrature.rule([-bound,-10.,0.,10.,bound],source_order)
            max_order=max(len(v['coefficients'])-1 for v in jets.values())
            derivatives={n:H.polynomial_basis(self.source_nodes,bound,size,n) for n in range(max_order+1)}
            self.amplitudes={};shared={};self.amplitude_aliases={}
            for si,jet in jets.items():
                key=tuple(jet['coefficients'])
                if key in shared:
                    first=shared[key];self.amplitudes[si]=self.amplitudes[first];self.amplitude_aliases[si]=first;continue
                matrix=np.zeros((len(self.source_nodes),size),complex)
                for n,coefficient in enumerate(key):
                    values=np.broadcast_to(np.asarray(sp.lambdify(self.r.zp,coefficient,'numpy')(self.source_nodes),complex),self.source_nodes.shape)
                    matrix+=values[:,None]*derivatives[n]
                self.amplitudes[si]=self.source_weights[:,None]*matrix;shared[key]=si;self.amplitude_aliases[si]=si
            self.size=size
        def prepare_gaussian(self,jets,bound,order,width=8.,momentum=0.):
            self.source_nodes,self.source_weights=H.BoundedSourceFourierQuadrature.rule([-bound,-10.,0.,10.,bound],order)
            x=self.source_nodes;field=np.exp(-(x/width)**2+1j*momentum*x);self.amplitudes={};self.size=1;boundary=[]
            for i,jet in jets.items():
                values=np.zeros(len(x),complex)
                for n,c in enumerate(jet['coefficients']):
                    function=sp.lambdify(self.r.zp,c,'numpy',cse=True,docstring_limit=0)
                    hermite=np.polynomial.hermite.hermval(x/width-0.5j*momentum*width,[0]*n+[1])
                    values+=np.asarray(function(x),complex)*(-1/width)**n*hermite*field
                    edge=np.asarray([-bound,bound]);edge_value=np.asarray(function(edge),complex)*(-1/width)**n*np.polynomial.hermite.hermval(edge/width-0.5j*momentum*width,[0]*n+[1])*np.exp(-(edge/width)**2+1j*momentum*edge)
                    boundary.append(float(np.max(abs(edge_value))))
                self.amplitudes[i]=(self.source_weights*values)[:,None]
            require(max(boundary)<1e-12,'Gaussian finite-source derivative boundary values');return boundary
        def row_batch(self,row,variables,points,positions):
            env={v:points[:,j] for j,v in enumerate(variables)};env[self.r.regulator]=np.full(len(points),self.setting['regulator']);profiles={};sources={};result=np.zeros((len(points),len(positions)),complex)
            for factor in row['factors']:
                expression=factor['coefficient']
                if self.wrong_sheet is not None:
                    lifted=lift_roots(expression,self.roots);replacements={v['carrier']:((-1 if p==self.wrong_sheet else 1)*sp.sqrt(v['bound'])) for p,v in self.roots.items()};expression=lifted.xreplace(replacements)
                c=self.coefficient_value(expression,env,positions,self.setting,profiles);s=self.source_basis(factor['sourceIndex'],env,sources)[:,0];result+=c*s[:,None]
            return result
        def row_integral(self,row,positions,outer_fixed=None):
            variables=tuple(l[0] for l in row['limits']);width=float(self.binding['abel']['width'].subs(self.r.regulator,self.setting['regulator']));points=[];weights=[];result=np.zeros(len(positions),complex);count=0
            def flush():
                nonlocal result,count
                if points:result+=np.asarray(weights)@self.row_batch(row,variables,np.asarray(points),positions);count+=len(points);points.clear();weights.clear()
            def descend(i,env,weight):
                if i<0:
                    points.append(tuple(env[v] for v in variables));weights.append(weight)
                    if len(points)>=self.batch_nodes:flush()
                    return
                v=variables[i]
                if outer_fixed is not None and i==len(variables)-1:
                    env[v]=outer_fixed;descend(i-1,env,weight);env.pop(v);return
                centers=[env[y if x==v else x] for x,y in self.binding['pairs'] if v in (x,y) and (y if x==v else x) in env]
                order=self.setting['innerOrders'][i] if centers else self.setting['outerOrder']
                nodes,masses,_=self.rule(-self.setting['momentumBound'],self.setting['momentumBound'],order,centers,width)
                for node,mass in zip(nodes,masses):env[v]=node;descend(i-1,env,weight*mass)
                env.pop(v)
            descend(len(variables)-1,{},1.);flush();require(np.isfinite(result).all(),'finite complete row integral');return {'value':result,'nodes':count}
    return Radiating


def check_rule(worker,roots,setting):
    b=next(iter(roots.values()))['b'];records=[]
    for centers in ((),(-.31,.11,.41)):
        nodes,weights,panels=worker.rule(-setting['momentumBound'],setting['momentumBound'],setting['outerOrder'],centers,.01)
        mass=float(weights.sum());moment=float(weights@nodes);Q=np.sqrt(np.asarray(b*b-nodes*nodes,complex))
        require(abs(mass-2*setting['momentumBound'])<1e-10 and abs(moment)<1e-9,'transformed positive Fourier measure and negative-leg orientation')
        require(np.isfinite(Q).all() and np.all(Q.real>=0) and np.all(Q.imag>=0) and np.all(abs(nodes)!=b),'physical outgoing branch excludes endpoint nodes')
        records.append({'nodes':nodes,'weights':weights,'panels':panels,'mass':mass,'oddMoment':moment,'physicalQ':Q})
    return records


def endpoint_values(worker,binding,roots,positions,J):
    records=[];b=next(iter(roots.values()))['b']
    for row in binding['rows']:
        variables=tuple(l[0] for l in row['limits'])
        for variable in variables:
            if variable not in roots:continue
            points=[];jac=[];labels=[]
            for end in (-1,1):
                for side in ('interior','exterior'):
                    for delta in (1e-2,1e-3,1e-4,1e-5):
                        k=end*b*(np.cos(delta) if side=='interior' else np.cosh(delta))
                        points.append([k if p==variable else (.43+.17*j) for j,p in enumerate(variables)])
                        jac.append(b*(np.sin(delta) if side=='interior' else np.sinh(delta)));labels.append((end,side,delta))
            points=np.asarray(points);raw=worker.row_batch(row,variables,points,positions);weighted=raw*np.asarray(jac)[:,None]
            record={'row':row['index'],'variable':variable,'points':points,'jacobians':np.asarray(jac),'labels':labels,'completeRowValues':raw,'weightedValues':weighted}
            J.blob(f'endpoint-values/row-{row["index"]}-{variables.index(variable)}.pickle',record)
            require(np.isfinite(weighted).all(),'finite source-bound transformed endpoint sequences',record)
            records.append({'row':row['index'],'variable':str(variable),'weightedMaximum':float(abs(weighted).max())})
    return records


def middle_check(worker,binding,roots,positions,J):
    from scipy.integrate import quad_vec
    rows=[row for row in binding['rows'] if len(row['limits'])==3]
    require(rows,'actual ordered three-momentum source rows')
    layout=[];selected=None
    for row in rows:
        variables=tuple(l[0] for l in row['limits']);middle=[p for p in variables if p.name=='s11cdMiddleNormalMomentum']
        require(len(middle)==1 and variables[-1]==middle[0] and middle[0] in roots,'actual middle leg is outermost')
        layout.append({'row':row['index'],'variables':variables,'middleIndex':2,'limits':row['limits']})
        probe=worker.row_batch(row,variables,np.asarray([[.37,.53,.71],[.47,.67,.83]]),positions)
        if selected is None and np.max(abs(probe))>1e-12:selected=row
    J.blob('middle/layout-and-selection.pickle',{'layout':layout,'selected':selected})
    require(selected is not None,'responsive complete native middle row')
    name='integration/middle/row-'+str(selected['index'])
    main=J.call(name+'/transformed-base',{'row':selected,'positions':positions,'setting':worker.setting},lambda:worker.row_integral(selected,positions))
    # Only the actual middle leg uses adaptive physical-k integration; all inner
    # panels are rebuilt per adaptive node. Shared inner rules are explicit.
    b=next(iter(roots.values()))['b'];K=worker.setting['momentumBound'];parts=[]
    for i,(lo,hi) in enumerate(zip((-K,-b,b),( -b,b,K))):
        def integrate(lo=lo,hi=hi):
            count=0
            def integrand(k):
                nonlocal count
                value=worker.row_integral(selected,positions,float(k));count+=value['nodes'];return value['value']
            value,error,info=quad_vec(integrand,lo,hi,epsabs=1e-9,epsrel=1e-7,quadrature='gk21',norm='max',limit=1000,full_output=True)
            result={'bounds':(lo,hi),'value':value,'error':float(error),'success':bool(info.success),'message':info.message,'evaluations':info.neval,'nestedNodes':count,'intervals':info.intervals,'intervalValues':info.integrals,'intervalErrors':info.errors}
            require(info.success and np.isfinite(value).all(),'independent physical-middle adaptive integration',result);return result
        parts.append(J.call(name+f'/direct-part-{i}',{'row':selected,'setting':worker.setting,'bounds':(lo,hi),'positions':positions,'sharedInnerRules':True},integrate))
    reference=sum(v['value'] for v in parts);error=sum(v['error'] for v in parts);coarse=main['value']
    previous=worker.setting;worker.setting=dict(previous,outerOrder=32)
    fine=J.call(name+'/transformed-refined',{'row':selected,'positions':positions,'setting':worker.setting},lambda:worker.row_integral(selected,positions));worker.setting=previous
    tolerance=1e-8+1e-6*abs(reference);difference=fine['value']-reference
    result={'row':selected['index'],'coarse':coarse,'fine':fine['value'],'reference':reference,'difference':difference,'tolerance':tolerance,'referenceError':error,'sharedInnerRules':'Same native inner rules rebuilt at each physical middle node; not an independent inner-rule accuracy proof.'}
    J.report('middle-integration-comparison',result);J.blob('middle/comparison.pickle',result)
    require(np.all(abs(difference)<=tolerance) and error<=float(np.max(tolerance)),'complete middle integral comparison',result)
    # Source-addressed controls operate on that same nonzero full row.
    worker.wrong_sheet=selected['limits'][-1][0]
    changed=J.call(name+'/wrong-middle-sheet',{'row':selected,'setting':worker.setting,'positions':positions,'negativeLeg':worker.wrong_sheet},lambda:worker.row_integral(selected,positions));worker.wrong_sheet=None
    worker.omit_jacobian=True
    omitted=J.call(name+'/omitted-transformed-jacobians',{'row':selected,'setting':worker.setting,'positions':positions},lambda:worker.row_integral(selected,positions));worker.omit_jacobian=False
    mutation={'wrongSheetMovement':changed['value']-coarse,'omittedJacobianMovement':omitted['value']-coarse}
    J.report('integration-mutations',mutation)
    require(min(float(np.max(abs(v))) for v in mutation.values())>1e-10,'actual integration controls respond',mutation)
    return result,mutation


def gaussian_derivative(x,n,width,momentum):
    return (-1/width)**n*np.polynomial.hermite.hermval(np.asarray(x)/width-.5j*momentum*width,[0]*n+[1])*np.exp(-(np.asarray(x)/width)**2+1j*momentum*np.asarray(x))


def uniform_gaussian(worker,binding,packet,roots,positions,J):
    from scipy.integrate import quad_vec
    require(binding['contrast']==0,'actual bound zero-contrast operator')
    active=[v for v in binding['termJoins'] if binding['actual']['cell',v['row'],v['column'],v['term']]!=0]
    # The required independent reference can be used only if its regulator-free
    # uniform pencil matches the actual regulator-free assembled operator.
    used={v['integralIndex'] for v in active}
    expressions=[binding['actual']['cell',v['row'],v['column'],v['term']] for v in active]
    expressions += [f['coefficient'] for i in used for f in binding['rows'][i]['factors']]
    expressions += list(x for m in binding['local'].values() for x in m)
    require(all(not x.has(worker.r.regulator) for x in expressions),'actual uniform operator has no Abel dependence',{'activeTerms':active,'expressions':expressions})
    end=packet['ends']['REFERENCE']['pencil'];w,k,Q=end['frequency'],end['momentum'],end['radical'];matrix=end['livePencil'].subs(w,3)
    scale=next(iter(roots.values()))['scale'];b=next(iter(roots.values()))['b'];K=worker.setting['momentumBound']
    require(sp.expand(end['wave'].subs(w,3))==sp.expand(Q**2-(scale**2)*(sp.Rational(1,25)-k**2)),'reference original wave/physical depth scaling')
    fn=sp.lambdify((k,Q),matrix,'numpy',cse=True,docstring_limit=0);results=[]
    for momentum in (0.,3*b):
        worker.prepare_gaussian(binding['jets'],worker.setting['sourceBound'],worker.setting['sourceNodes'],8.,momentum)
        row_values={}
        for i in sorted(used):
            row=binding['rows'][i]
            row_values[i]=J.call(f'uniform-Gaussian/{momentum}/row-{i}',{'row':row,'setting':worker.setting,'width':8.,'momentum':momentum,'positions':positions},lambda row=row:worker.row_integral(row,positions))['value']
        # Evaluate every face-driving source column; eW is mandatory, theta is
        # retained too. No transverse zero substitutes for this action test.
        for column in (3,4):
            assembled=np.zeros((5,len(positions)),complex)
            for n,m in binding['local'].items():
                d=gaussian_derivative(positions,n,8.,momentum)
                for i in range(5):assembled[i]+=np.asarray(sp.lambdify(worker.r.z,m[i,column],'numpy',cse=True)(positions),complex)*d
            for term in active:
                if term['column']!=column:continue
                c=binding['actual']['cell',term['row'],column,term['term']]
                assembled[term['row']]+=np.asarray(sp.lambdify(worker.r.z,c,'numpy',cse=True)(positions),complex)*row_values[term['integralIndex']]
            def independent():
                def f(p):
                    q=float(scale)*np.sqrt(complex(b*b-p*p));source=8*np.sqrt(np.pi)*np.exp(-16*(p-momentum)**2)
                    return np.asarray(fn(p,q),complex)[:,column,None]*np.exp(1j*p*positions)[None,:]*source/(2*np.pi)
                parts=[]
                for lo,hi in zip((-K,-b,b),(-b,b,K)):
                    value,error,info=quad_vec(f,lo,hi,epsabs=1e-10,epsrel=1e-8,norm='max',limit=1000,full_output=True,quadrature='gk21')
                    parts.append({'bounds':(lo,hi),'value':value,'error':float(error),'success':bool(info.success),'message':info.message})
                require(all(p['success'] for p in parts),'independent reference Gaussian integral',parts)
                return {'value':sum(v['value'] for v in parts),'error':sum(v['error'] for v in parts),'parts':parts}
            ref=J.call(f'uniform-Gaussian/{momentum}/reference-column-{column}',{'matrix':matrix,'physicalScale':scale,'column':column,'width':8.,'momentum':momentum,'positions':positions,'K':K},independent)
            tolerance=1e-8+1e-6*abs(ref['value']);difference=assembled-ref['value']
            record={'column':column,'momentum':momentum,'assembled':assembled,'reference':ref['value'],'difference':difference,'tolerance':tolerance,'referenceError':ref['error'],'finiteSourceEndpointControl':'All actual source-jet boundary values below 1e-12; width8, finite length64.'}
            J.blob(f'uniform-Gaussian/{momentum}/comparison-column-{column}.pickle',record);J.report(f'uniform-Gaussian-{momentum}-column-{column}',record)
            require(np.max(abs(ref['value']))>1e-10 and np.all(abs(difference)<=tolerance) and ref['error']<=np.max(tolerance),'responsive zero-contrast original-pencil/action comparison',record);results.append(record)
    return results


def source_branch_joins(worker,binding,roots,positions,J):
    records=[]
    for row in binding['rows']:
        variables=tuple(l[0] for l in row['limits']);samples=np.asarray([-.73,-.31,-.13,.07,.29,.61]);env={v:np.roll(samples,j) for j,v in enumerate(variables)};env[worker.r.regulator]=np.full(len(samples),worker.setting['regulator'])
        for p,root in roots.items():
            if p in env:env[root['carrier']]=float(root['scale'])*np.sqrt(np.asarray(root['b']**2-env[p]**2,complex))
        for fi,factor in enumerate(row['factors']):
            original=factor['coefficient'];lifted=lift_roots(original,roots)
            native=worker.coefficient_value(original,env,positions,worker.setting,{})
            mapped=worker.coefficient_value(lifted,env,positions,worker.setting,{})
            residual=(native-mapped)/(1+abs(native));record={'row':row['index'],'factor':fi,'original':original,'lifted':lifted,'environment':env,'positions':positions,'native':native,'physicalRootEvaluation':mapped,'scaledResidual':residual}
            J.blob(f'branch-joins/row-{row["index"]}-factor-{fi}.pickle',record)
            require(np.max(abs(residual))<1e-11,'actual source coefficient/physical-root branch join',record)
            records.append({'row':row['index'],'factor':fi,'maximumScaledResidual':float(abs(residual).max())})
    return records
