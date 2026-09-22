#!/usr/bin/env python3
"""Missing source coefficients by structural linear projection, not old Poly replay."""
import argparse
import ast
import cmath
import inspect
import json
from pathlib import Path
import resource
import signal
import time

import sympy as sp
from sympy.core.function import AppliedUndef
import S11c_d_remaining_case_frequency_end_continuation_inputs as io
import S11c_d_remaining_case_frequency_end_pilot as storage

M, REPO = io.M, io.REPO
F = REPO/'_scratch/s11c/s11c-remaining-case-frequency-20260921'
CP = M/'S11c_d_remaining_case_frequency_scalar_bindings_checkpoint.json'
CP_SHA = 'd1961c25f725959e896e62e49c45d389fdd503575d27803a89055175f6308b50'
PLAN = M/'S11c_d_remaining_case_frequency_source_coefficients_plan.md'
require, same = io.require, io.same


def numeric(expression, coordinate, point, assigned, memo=None):
    """Direct numerical tree interpretation; no projection, substitution or CAS."""
    memo = {} if memo is None else memo
    if expression in memo:
        return memo[expression]
    if expression in assigned:
        value = assigned[expression]
    elif expression == coordinate:
        value = point
    elif expression is sp.I:
        value = 1j
    elif isinstance(expression, sp.Rational):
        value = int(expression.p)/int(expression.q)
    elif isinstance(expression, sp.Float):
        value = float(expression)
    elif isinstance(expression, sp.Add):
        value = sum(numeric(v,coordinate,point,assigned,memo) for v in expression.args)
    elif isinstance(expression, sp.Mul):
        value = 1
        for v in expression.args:
            value *= numeric(v,coordinate,point,assigned,memo)
    elif isinstance(expression, sp.Pow) and isinstance(expression.exp, sp.Integer):
        value = numeric(expression.base,coordinate,point,assigned,memo)**int(expression.exp)
    elif expression.func is sp.tanh:
        value = cmath.tanh(numeric(expression.args[0],coordinate,point,assigned,memo))
    else:
        raise TypeError(('unsupported numerical source node',type(expression),repr(expression)))
    memo[expression] = value
    return value


class Projection:
    """Distribute a linear expression into its actual saved probe nodes.

    Independent subtrees are returned verbatim. Only Add/Mul coefficient
    combinations are new symbolic operations. No derivative node is created.
    """
    def __init__(self, journal):
        self.journal = journal
        self.calls = []
        self.atlas = {}

    def combine(self, function, args, context, consumer):
        if len(args) == 1:
            return args[0], {'kind':'single-saved-operand'}
        requested = {'function':function,'args':tuple(args),'kwargs':{},'unitContext':dict(context,jetOrder=consumer['jetOrder'])}
        key = (function,tuple(args))
        matches = [v for v in self.atlas.get(key,[]) if same(v['input'],requested)]
        number = len(self.calls)
        prefix = 'coefficient-operations/'+str(number)
        inp = self.journal.write(prefix+'/input.pickle',{'call':requested,'consumer':consumer})
        if matches:
            value = matches[0]['value']; owner = matches[0]['owner']; kind = 'REUSED_NEW_COMPLETE_CALL'
        else:
            value = {'Add':sp.Add,'Mul':sp.Mul}[function](*args)
            owner = number; kind = 'NEW_COEFFICIENT_OPERATION'
        output = self.journal.write(prefix+'/value.pickle',value)
        receipt = {'kind':kind,'owner':owner,'input':inp,'value':output,'function':function}
        self.journal.json(prefix+'/completed.json',receipt)
        self.calls.append(receipt)
        if not matches:
            self.atlas.setdefault(key,[]).append({'input':requested,'value':value,'owner':number})
        return value, receipt

    def project(self, expression, node_orders, context, family):
        memo = {}; records = []
        def visit(node, address):
            if node in memo:
                return memo[node]
            op = None
            if node in node_orders:
                result = {node_orders[node]:sp.S.One}
                kind = 'actual-saved-probe-node'
            elif not node.has(AppliedUndef,sp.Derivative):
                result = {None:node}; kind = 'unchanged-independent-subtree'
            elif isinstance(node, sp.Add):
                children = [visit(v,address+(i,)) for i,v in enumerate(node.args)]
                orders = set().union(*(v.keys() for v in children)); result = {}; op = {}
                for order in sorted(orders,key=lambda v:-1 if v is None else v):
                    args = tuple(v[order] for v in children if order in v)
                    result[order],op[order] = self.combine('Add',args,context,{'family':family,'argsAddress':address,'jetOrder':order})
                kind = 'linear-add'
            elif isinstance(node, sp.Mul):
                children = [visit(v,address+(i,)) for i,v in enumerate(node.args)]
                dependent = [i for i,v in enumerate(children) if any(k is not None for k in v)]
                require(len(dependent)==1,'one probe-dependent factor in every source product')
                index = dependent[0]; result = {}; op = {}
                for order,coefficient in children[index].items():
                    args = tuple(coefficient if i==index else children[i][None] for i in range(len(children)))
                    result[order],op[order] = self.combine('Mul',args,context,{'family':family,'argsAddress':address,'jetOrder':order})
                kind = 'linear-product'
            else:
                self.journal.write('families/'+str(family)+'/unsupported-node.pickle',{'node':node,'address':address})
                raise TypeError(('unsupported probe-dependent source node',type(node),address))
            records.append({'node':node,'argsAddress':address,'kind':kind,'result':result,'operations':op})
            memo[node] = result
            return result
        result = visit(expression,())
        self.journal.write('families/'+str(family)+'/projection-tree.pickle',records)
        require(None not in result,'homogeneous linear source; no missing affine contribution')
        return result


def prohibit():
    def forbidden(*args,**kwargs):
        raise RuntimeError('completed native science is outside structural source coefficient projection')
    for name in ('diff','lambdify','cancel','factor','expand','simplify','solve','resultant','gcd','gcdex','integrate','symbols'):
        setattr(sp,name,forbidden)
    # Keep SymPy classes and Basic codec methods intact for faithful restoration.
    for name in ('source_jets','load','construct','main'):
        setattr(io.native.f,name,forbidden)
    for name in ('load','seeds','maps','continue_pair','focused','main'):
        setattr(io.native,name,forbidden)
    io.native.Pair.__init__ = forbidden


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    base=parser.parse_args().run_directory.resolve();base.relative_to(REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    reader,journal=io.Reader(),storage.Journal(base)
    cp=reader.json(CP,CP_SHA);origin=Path(cp['runDirectory']);review=cp['independentSavedReview'];vr=Path(review['runDirectory'])
    require(cp['status']=='ACCEPTED_CASE_FREQUENCY_SCALAR_BINDINGS','accepted actual full scalar inputs')
    reader.retain(origin/'checks.json',cp['checksSha256']);reader.retain(vr/'checks.json',review['checksSha256'])
    name='source-jet-input-catalogue.json';catalogue=reader.json(vr/name,review['artifacts'][name]['sha256'])
    inspection=reader.json(F/'source-jet-operands-recovery-01/complete/checks.json','11fb53943accfa0a500edbb8422e1a1cef11a0d35d14e72bf84878c9e546f8d0')
    require(inspection['newScientificCalls']==0,'completed operand inspection made no coefficients')
    for p in (Path(__file__).resolve(),PLAN,Path(io.__file__),Path(storage.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md'):
        reader.retain(p)
    matrix_route=cp['upstreamCheckpoints']['baselineFrequencyMatrix'];matrix_cp=reader.json(matrix_route['logical'],matrix_route['sha256'])
    nativefile=M/'S11c_d_finite_scattering.py';nativehash=matrix_cp['sourceFiles']['_measurements/S11c_d_finite_scattering.py']
    reader.retain(nativefile,nativehash);reader.retain(Path(matrix_cp['runDirectory'])/'source/_measurements/S11c_d_finite_scattering.py',nativehash)
    tree=ast.parse(nativefile.read_text());bodies={n.name:ast.unparse(n) for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in ('source_jets','BasisMomentum')}
    journal.json('native-caller-and-new-projection.json',{'nativeSource':reader.retain(nativefile),'bodies':bodies,
        'newProjectionSource':inspect.getsource(Projection),'numericCheckSource':inspect.getsource(numeric),
        'nativeSourceJetsCalled':False,'scope':'New structural coefficient projection only for unmatched whole inputs; native source_jets intermediates are not reconstructed.'})
    prohibit();cache={}
    def packet(record):
        path=record.get('logical',record.get('path'));r=reader.retain(path,record['sha256'])
        if r['canonical'] not in cache:cache[r['canonical']]=reader.packet(path,record['sha256'])
        return cache[r['canonical']]
    projected=[];routes={};projection=Projection(journal);maximum=0.;controls=[]
    for request in catalogue['requests']:
        case,si=request['case'],request['sourceIndex'];old=packet(request['nativeSource']['packet']);binding=packet(request['boundAmplitude']['packet'])
        source=old['bound']['sources'][0,si];prior=old['jets'][si];expression=binding['actual']['source',si];probe=prior['probe'];coordinate=probe.args[0]
        full={'expression':expression,'probe':probe,'coordinate':coordinate,'amplitudeUnit':source['amplitudeUnit'],'integralUnit':source['integralUnit']}
        consumer={'case':case,'sourceIndex':si,'fullInputRoutes':request}
        if request['matches']:
            owner=request['matches'][0];value=packet(owner['packet'])[owner['keys'][0]][owner['keys'][1]]
            require(same((expression,probe,source['amplitudeUnit'],source['integralUnit']),
                         (value['originalBoundAmplitude'],value['probe'],value['amplitudeUnit'],value['integralUnit'])),'whole accepted coefficient call inputs')
            kind='SAVED_NATIVE_COMPLETE_JET';target={'packet':owner['packet'],'keys':owner['keys']}
        else:
            matches=[v for v in projected if same(v['input'],full)]
            if matches:
                chosen=matches[0];target=chosen['target'];kind='REUSED_NEW_PROJECTION'
            else:
                family=len(projected);prefix='families/'+str(family)
                journal.write(prefix+'/input.pickle',{'call':full,'consumer':consumer})
                derivatives=expression.atoms(sp.Derivative)
                require(expression.atoms(AppliedUndef)=={probe} and all(v.expr==probe and all(x==coordinate for x,n in v.variable_count) for v in derivatives),'actual one-field derivative nodes')
                node_orders={probe:0}
                for node in derivatives:node_orders[node]=int(sum(n for x,n in node.variable_count))
                degree=max(node_orders.values());context={k:full[k] for k in ('probe','coordinate','amplitudeUnit','integralUnit')}
                result=projection.project(expression,node_orders,context,family)
                require(set(result)==set(range(degree+1)),'all actual derivative orders have structural coefficients')
                coefficients=[result[n] for n in range(degree+1)]
                value={'column':int(probe.func.__name__.removeprefix('s11cdPencilProbe')),'coefficients':coefficients,'probe':probe,
                       'originalBoundAmplitude':expression,'amplitudeUnit':full['amplitudeUnit'],'integralUnit':full['integralUnit'],
                       'algorithm':'saved-expression-linear-tree-projection-v1','structuralProof':{'file':journal.artifacts[prefix+'/projection-tree.pickle']},
                       'scope':'Actual new coefficients. No native Poly/intermediate/residual is claimed or reconstructed.'}
                output=journal.write(prefix+'/value.pickle',value)
                journal.json(prefix+'/completed.json',{'input':journal.artifacts[prefix+'/input.pickle'],'value':output,'operationsCompleted':len(projection.calls)})
                require(all(not a.has(AppliedUndef,sp.Derivative,sp.Subs,sp.Integral) and a.free_symbols <= {coordinate} for a in coefficients),'complete source-only coefficient leaves')
                points=(-41.,-9.,.37,13.,43.,.61+.23j)
                vectors=(tuple(complex(n+1,(-1)**n*.43) for n in range(degree+1)),
                         tuple(complex((-1)**n/(n+2),n+.71) for n in range(degree+1)))
                journal.write(prefix+'/comparison-input.pickle',{'source':full,'coefficients':coefficients,'nodeOrders':node_orders,'points':points,'jetVectors':vectors,'tolerance':2e-12})
                comparisons=[];family_max=0.;coefficient_response=0.;derivative_response=0.
                for point in points:
                    values=[numeric(a,coordinate,point,{}) for a in coefficients]
                    for vector in vectors:
                        assigned={node:vector[n] for node,n in node_orders.items()}
                        actual=numeric(expression,coordinate,point,assigned)
                        reconstructed=sum(a*b for a,b in zip(values,vector))
                        difference=reconstructed-actual;scaled=abs(difference)/(1+abs(actual));family_max=max(family_max,scaled)
                        # Actual numerical changed operands, independent of any CAS recomputation.
                        changed=list(values);changed[0]+=1e-3
                        mutated=sum(a*b for a,b in zip(changed,vector));coefficient_response=max(coefficient_response,abs(mutated-actual))
                        swapped=list(values)
                        if degree:swapped[0],swapped[-1]=swapped[-1],swapped[0]
                        else:swapped[0]=-swapped[0]
                        permuted=sum(a*b for a,b in zip(swapped,vector));derivative_response=max(derivative_response,abs(permuted-actual))
                        comparisons.append({'point':point,'jetVector':vector,'directAmplitude':actual,'coefficientValues':values,
                            'projectedAmplitude':reconstructed,'difference':difference,'scaledDifference':scaled,
                            'coefficientMutationValues':changed,'coefficientMutationAmplitude':mutated,
                            'jetOrderMutationValues':swapped,'jetOrderMutationAmplitude':permuted})
                journal.write(prefix+'/comparison-value.pickle',comparisons)
                summary={'maximumScaledDifference':family_max,'coefficientMutationResponse':coefficient_response,
                         'jetOrderOrSignMutationResponse':derivative_response,'actualComparisonCount':len(comparisons),'tolerance':2e-12}
                journal.json(prefix+'/comparison.json',summary)
                require(family_max<2e-12 and coefficient_response>1e-5 and derivative_response>1e-12,'focused projection/source and changed coefficient/order comparisons')
                maximum=max(maximum,family_max);controls.append(summary)
                target={'packet':output,'keys':[]};projected.append({'input':full,'target':target,'owner':[case,si]});kind='NEW_STRUCTURAL_PROJECTION'
        routes.setdefault(case,[]).append({'consumer':consumer,'kind':kind,'result':target})
    for case,items in routes.items():
        journal.json('cases/'+case+'/source-coefficient-routes.json',items)
        journal.json('cases/'+case+'/summary.json',{'sourceUses':len(items),'savedNative':sum(v['kind']=='SAVED_NATIVE_COMPLETE_JET' for v in items),
            'newOwnerUses':sum(v['kind']=='NEW_STRUCTURAL_PROJECTION' for v in items),'newAliasUses':sum(v['kind']=='REUSED_NEW_PROJECTION' for v in items)})
    journal.json('operation-inventory.json',projection.calls)
    reader.postcheck();journal.json('inputs.json',{'acceptedScalarCheckpoint':reader.retain(CP,CP_SHA),'consumedRoutes':reader.routes})
    counts={'sourceUses':sum(map(len,routes.values())),'savedNativeUses':sum(v['kind']=='SAVED_NATIVE_COMPLETE_JET' for r in routes.values() for v in r),
            'newFullProjections':len(projected),'sharedNewUses':sum(v['kind']=='REUSED_NEW_PROJECTION' for r in routes.values() for v in r),
            'coefficientOperationCalls':len(projection.calls),'newCoefficientOperations':sum(v['kind']=='NEW_COEFFICIENT_OPERATION' for v in projection.calls)}
    require(counts['sourceUses']==120 and counts['savedNativeUses']==100 and counts['newFullProjections']==10,'actual saved catalogue census')
    checks={'status':'COMPLETED_MISSING_FREQUENCY_SOURCE_COEFFICIENTS','counts':counts,'maximumScaledDirectComparison':maximum,
        'focusedControls':controls,'allConsumedHashesUnchanged':True,'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-start,
        'scope':'Source coefficients only at1-0.01i under exploratoryAcceptanceV1. New linear tree projection; no old derivative/Poly/source_jets/seed/end map/science replay. No quadrature or full response acceptance.'}
    journal.json('checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
