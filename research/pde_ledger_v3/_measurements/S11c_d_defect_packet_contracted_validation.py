"""New guarded operand joins; no scientific imports or restored data at import."""
import ast


def require(value,message):
    if value is not True:raise ValueError(message)


def original_selector(unit_return,address,role,coefficients):
    """Select an existing complete interface, including its derivative selector."""
    key='matchingXFamilyInterfaces' if role=='source' else 'matchingYFamilyInterfaces'
    rows=unit_return[key]
    require(len(rows)==1,'one original '+role+' Fourier interface')
    row=rows[0]
    if role=='consumer':
        require(address['addressId'] in row['addresses'],'original Y address membership')
        row=row['spec']
    require(row['role']==('X' if role=='source' else 'Y'),'original Fourier role')
    require(type(row['argumentDerivative']) is int and row['argumentDerivative']==0,'original argumentDerivative zero only')
    require(row['fieldId']==address[role+'Transform']['coefficientId'] and row['coefficientsAscending']==coefficients,'original constant coefficient interface')
    require(row['spatialOrders']==(address['jet']['spatialOrders'] if role=='source' else [0,0,0]),'original spatial derivative orders')
    require(row['timeOrder']==(address['jet']['timeOrder'] if role=='source' else 0),'original time derivative order')
    require(row['center']==('-5/2' if role=='source' else '5/2') and row['width']=='8' and row['profileLength']=='10','original Gaussian and profile scales')
    return row


def whole_densities(raw,manifest,J,sp,C,exact):
    """Attach restored original densities to whole definitions with actual maps.

    These are new cross-stage argument joins. Original density construction,
    saved inner arithmetic and full contraction proofs are not called again.
    """
    R=lambda x:C.restore_scalar(sp,x);old=lambda n:raw['contraction/complete/'+n+'.json']
    source=old('new-original-density-interpretation-input');returned=old('new-original-density-interpretation-return')
    J.emit('restored-original-density-operands',{'input':source,'return':returned,'noConstructorCalled':True})
    require(source['source']['fragmentSha256']==raw['numeric-source-contracts.json']['fragments']['inner-components']['fragmentSha256'],'original kernel AST receipt')
    arguments={n:R(v) for n,v in source['arguments'].items()};densities=[R(v) for v in returned['orderedDensities']]
    require(len(densities)==4,'four original added densities')
    for name,value in zip(('J','Dr','Dh','Dq'),returned['orderedDensities']):
        require(old('new-full-factorization-'+name+'-input')['context']['originalDensity']==value,'actual original factorization density '+name)
    inners={name:raw['saved/inner/new-runtime-'+name+'-input.json'] for name in ('J-arithmetic','D-arithmetic-sum')}
    for name,value in inners.items():
        ret=raw['saved/inner/new-runtime-'+name+'-return.json']
        J.emit('restored-runtime-arithmetic-'+name,{'input':value,'return':ret,'oppositeSidesNotRecomputed':True})
        require(ret['cancelled']=={'text':'0','srepr':'Integer(0)'},'inherited runtime arithmetic zero')
    packet=C.named_symbols(R(inners['J-arithmetic']['left'])+R(inners['D-arithmetic-sum']['left']))
    p={n:packet['packet_'+n] for n in ('k','l','t','qi','qo','qh','qs')}
    aa=R(raw['saved/inner/runtime-a-input.json']['left']);mu=R(raw['saved/inner/runtime-mu-input.json']['left'])
    scale=raw['accepted-units/parameter-quantity-origins.json']['context']['numeric']
    W,L=R(scale['W_0']),R(scale['L_W'])
    A=lambda z:L*z/(4*sp.sinh(sp.pi*L*z/2))
    values={**p,'a':aa,'mu':mu,'beta':aa*mu,'W':W,'L':L,'A1':A(p['t']),'A2':A(p['l']-p['k']-p['t']),'I':sp.I}
    require(set(arguments)==set(values),'complete original density argument map')
    mapping={arguments[n]:v for n,v in values.items() if n!='I'}
    bound=[v.xreplace(mapping) for v in densities]
    exact('new-original-density-to-runtime-J',bound[0],R(inners['J-arithmetic']['left']),{'originalInput':source,'originalReturn':returned,'actualMap':list(mapping.items()),'nativeW':scale['W_0'],'nativeL':scale['L_W']})
    exact('new-original-density-to-runtime-added-D',sum(bound[1:]),R(inners['D-arithmetic-sum']['left']),{'originalInput':source,'actualMap':list(mapping.items()),'threeSeparateAddedPieces':bound[1:]})
    tags=raw['saved/pressure/whole-tags.json'];whole={}
    names={'Jwhole':{'reference_k':'k','reference_l':'l','reference_t':'t','reference_qi':'qi','reference_qo':'qo','reference_qm':'qh'},'Dwhole':{'k':'k','grazing_output':'l','grazing_transfer':'t','grazing_qi':'qi','grazing_qo':'qo','grazing_qh':'qh','grazing_qs':'qs'}}
    for name,trans in names.items():
        alias='saved/whole-origin/'+name+'.json';tag=tags[name];definition=raw[alias]
        J.emit('new-whole-origin-input-'+name,{'tag':tag,'definition':definition,'originalReceipt':manifest['savedInputs'][alias]})
        require(tag['savedDefinition']==definition and tag['sha256']==manifest['savedInputs'][alias]['sha256'],'actual whole definition bytes '+name)
        variable='t' if name=='Jwhole' else 'td'
        require(tag['freeMomenta']==['l','k'] and tag['boundVariable']==variable and tag['middleMomentum']=='k+'+variable,'actual whole variable signature')
        require(tag['reflectedMomentum']==('l-td' if name=='Dwhole' else None),'actual distinct reflected momentum')
        density=R(definition['density']);symbols=C.named_symbols(density)
        freq='reference_unrestricted_frequency' if name=='Jwhole' else 'grazing_unrestricted_frequency'
        require(set(symbols)==set(trans)|{freq},'complete whole density symbol map')
        mp={symbols[n]:p[v] for n,v in trans.items()};mp[symbols[freq]]=sp.Integer(raw['preflight/physical-plan.json']['frequency'])
        mapped=density.xreplace(mp);key='J-arithmetic' if name=='Jwhole' else 'D-arithmetic-sum'
        exact('new-whole-density-to-runtime-'+name,mapped,R(inners[key]['right']),{'tag':tag,'originalWhole':definition,'actualMap':list(mp.items()),'inheritedArithmetic':inners[key],'wholeNotIntegratedAgain':True})
        whole[name]=mapped
    return {'whole':whole,'pieces':dict(zip(('J','Dr','Dh','Dq'),bound)),'packet':p,'originalDensities':densities}


def response_selection(raw,J,sp,C,exact,address,adapter,contraction_return,joined):
    R=lambda x:C.restore_scalar(sp,x);label='-'.join((address['face'],address['slot'],address['component']))
    args=raw['saved/preflight/numeric-factor-'+label+'-arguments.json'];response=R(address['responseCoefficient'])
    J.emit('new-response-selection-input-'+str(address['addressId']),{'actualAddress':address,'adapter':adapter,'originalArguments':args,'contractionReturn':contraction_return,'joinedDensities':joined})
    require(args['actualCompleteFactor']==address['responseCoefficient']==adapter['original'],'complete live native response')
    calls={R(a):R(b) for a,b in args['permittedCalls']};observed=response.atoms(sp.Function)
    require(observed==set(R(v) for v in args['actualCalls']) and observed<=set(calls),'complete original function call set')
    symbols=C.named_symbols(response);p=joined['packet']
    xy={symbols['composition_'+n]:p[n] for n in ('k','l') if 'composition_'+n in symbols}
    mapped=response.xreplace(calls).xreplace(xy)
    exact('new-live-native-template-'+str(address['addressId']),mapped,R(adapter['template']),{'actualAddress':address,'originalCallMap':args,'wholeDefinitions':raw['saved/pressure/whole-tags.json']})
    selected='packet_J' if address['component']=='NATIVE_MIXED_ITERATION' else 'packet_D'
    tagname=('Jwhole_' if selected=='packet_J' else 'Dwhole_')+address['face']
    selected_calls=[v for v in observed if v.func.__name__==tagname]
    require(len(selected_calls)==1,'one actual selected whole call')
    call=selected_calls[0];cs=symbols['composition_cs']
    require(call.args==(symbols['composition_l'],symbols['composition_k'],sp.Integer(3),cs,sp.Rational(1,5),sp.Rational(1,10),sp.Integer(1),sp.Integer(10)),'actual native whole argument order')
    if selected=='packet_J':require(contraction_return['contractedOrdinaryComponent'] is True,'explicit ordinary-only native component selection')
    expression=R(adapter['template']);template_symbols=C.named_symbols(expression)
    chosen=template_symbols[selected];density=joined['whole']['Jwhole' if selected=='packet_J' else 'Dwhole']
    injection={chosen:density}
    if selected=='packet_J':
        require('packet_H' in template_symbols,'surviving H kept separate')
        injection[template_symbols['packet_H']]=sp.S.Zero
    require(set(n for n in template_symbols if n in ('packet_H','packet_J','packet_D'))==set(s.name for s in injection),'only original selected whole and separate H')
    exact('new-whole-selected-density-'+str(address['addressId']),expression.xreplace(injection),density,{'actualTemplate':adapter,'selectedWholeCall':call,'wholeSignature':raw['saved/pressure/whole-tags.json']['Jwhole' if selected=='packet_J' else 'Dwhole'],'injection':list(injection.items()),'HNotZeroPhysically':'H is a separate pending summand; only its formal tag is zeroed to select J','noExtraResolvents':True,'noSecondMiddleIntegral':True})


def source_ast_joins(contracts,numeric_tree,prepare_tree,C,J):
    """Compare new formula bodies, without executing any old callable."""
    J.emit('new-source-ast-join-input',{'originalContracts':contracts,'newNumericSource':ast.unparse(numeric_tree),'newPreparationSource':ast.unparse(prepare_tree)})
    for name in ('exponential_moment','weighted_tail'):
        old=C.fragment(contracts['fragments']['tail-'+name.replace('_','-')])
        nodes=[n for n in ast.walk(prepare_tree) if isinstance(n,ast.FunctionDef) and n.name==name]
        require(len(nodes)==1 and ast.dump(nodes[0])==ast.dump(old),'same original tail function AST '+name)
    contributions=C.fragment(contracts['fragments']['tail-contributions'])
    old_values={ast.dump(n.value) for n in ast.walk(contributions) if isinstance(n,(ast.Assign,ast.AugAssign))}
    new_values=[n.value for n in ast.walk(prepare_tree) if
                isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ('outerbound','middlebound') for t in n.targets)
                and any(isinstance(v,ast.Call) and isinstance(v.func,ast.Name) and v.func.id=='weighted_tail' for v in ast.walk(n.value))]
    new_values += [n.value for n in ast.walk(prepare_tree) if isinstance(n,ast.AugAssign) and isinstance(n.target,ast.Name) and n.target.id=='middlebound']
    require(len(new_values)==4 and all(ast.dump(v) in old_values for v in new_values),'all actual enlarged positive tail contribution ASTs')
    old=C.fragment(contracts['fragments']['inner-profile'])
    new=next(n for n in ast.walk(numeric_tree) if isinstance(n,ast.FunctionDef) and n.name=='profile')
    # The ordinary sinh branch retains the exact original expression; the new
    # small-argument branch separately records its explicit truncation envelope.
    old_returns=[n.value for n in ast.walk(old) if isinstance(n,ast.Return)]
    new_return=new.body[-1].value.elts[0].args[0]
    require(any(ast.dump(new_return)==ast.dump(v) for v in old_returns),'actual native sinh expression AST')
    J.emit('new-profile-tail-source-joins',{'originalProfile':contracts['fragments']['inner-profile'],'newProfile':ast.unparse(new),'weightedTail':contracts['fragments']['tail-weighted-tail'],'moment':contracts['fragments']['tail-exponential-moment'],'contributions':contracts['fragments']['tail-contributions'],'oldFunctionsCalled':False})
