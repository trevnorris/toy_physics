"""Restored units for exact cancellation operands, including computed zeros."""
import re
import sympy as sp


def cancellation_units(engine, proof, strong, dimensions):
    d=dimensions
    fields=[d.known[sp.Function('s11cdReducedField'+v)] for v in ('u1','u2','u3','theta','eW')]
    amplitudes=tuple(next(s for s in d.known if isinstance(s,sp.Symbol) and s.name=='s11cdCurrentPlusAmplitude'+str(i)) for i in range(5))
    def add(*values):return tuple(sum(v) for v in zip(*values))
    def negative(value):return tuple(-v for v in value)
    rows={}
    for i in (3,4):
        j=next(j for j in range(5) if strong[i,j]!=0)
        rows[i]=add(d.measure(strong[i,j]),fields[j])
    def root_unit(label):
        if re.fullmatch(r'face_-?1_mechanical_load',label):return rows[4]
        match=re.fullmatch(r'(mass|mechanical)_(increment|face|join|sum)_([0-4])',label)
        if match:return add(rows[3 if match[1]=='mass' else 4],negative(fields[int(match[3])]))
        raise ValueError(('untyped rational operation',label))
    result={};amplitude_cache={};coefficient_cache={}
    for name,record in proof.items():
        if name.endswith('_reconstruction'):
            unit=root_unit(name.removesuffix('_reconstruction'))
        else:
            label,tail=name.rsplit('_amplitude_',1);ai=int(tail.split('_grade_')[0])
            if label not in amplitude_cache:
                source=proof[label+'_reconstruction']['BEFORE']
                amplitude_cache[label]=engine.polynomial_terms(source,amplitudes)
            if (label,ai) not in coefficient_cache:
                powers,coefficient=amplitude_cache[label][ai]
                _,denominator=sp.together(coefficient).as_numer_denom()
                grades=engine.PHYSICAL_METADATA.generators[1:]
                dynamic=sp.Mul(*(f for f in sp.Mul.make_args(denominator) if f.has(*grades)))
                amplitude_unit=tuple(sum(p*f[j] for p,f in zip(powers,fields)) for j in range(3))
                coefficient_cache[label,ai]=add(root_unit(label),d.measure(dynamic),negative(amplitude_unit))
            unit=coefficient_cache[label,ai]
        before_den=d.measure(record['SOURCE_DENOMINATOR'])
        after_den=d.measure(record['RESULT_DENOMINATOR'])
        result[name]={'BEFORE':unit,'AFTER':unit,
            'SOURCE_NUMERATOR':add(unit,before_den),'SOURCE_DENOMINATOR':before_den,
            'RESULT_NUMERATOR':add(unit,after_den),'RESULT_DENOMINATOR':after_den,
            'CROSS_PRODUCT_RESIDUAL':add(unit,before_den,after_den),
            'SOURCE_DENOMINATOR_DOMAIN':d.zero,'RESULT_DENOMINATOR_DOMAIN':d.zero}
        for kind in ('SOURCE','RESULT'):
            key='ORIGINAL_'+kind+'_DENOMINATOR'
            if key in record:
                result[name][key]=d.measure(record[key])
                result[name][key+'_DOMAIN']=d.zero
    return result


def cancellation_packet(engine, modes, name, key, record, unit):
    """Exact rational grades are represented by both computed polynomial parts."""
    value=record[key];body=engine.cas(value);representation='POLYNOMIAL'
    metadata_body=record[key.removesuffix('_DOMAIN')] if key.endswith('_DOMAIN') else body
    try:
        engine.PHYSICAL_METADATA.coefficients(metadata_body)
        unit_map={path:unit for path,_ in engine.leaves(engine.cas(metadata_body))}
    except NotImplementedError:
        if key.endswith('_DOMAIN'):raise
        numerator,denominator=sp.together(value).as_numer_denom()
        denominator_unit=engine.PHYSICAL_METADATA.dimensions.measure(denominator)
        numerator_unit=tuple(a+b for a,b in zip(unit,denominator_unit))
        metadata_body=engine.cas({'EXACT_NUMERATOR':numerator,'EXACT_DENOMINATOR':denominator})
        unit_map={('EXACT_NUMERATOR',):numerator_unit,('EXACT_DENOMINATOR',):denominator_unit}
        representation='EXACT_RATIONAL_NUMERATOR_DENOMINATOR'
    metadata=modes.numeric_metadata(metadata_body,lambda path:unit_map[path])
    if representation!='POLYNOMIAL':
        metadata=engine.cas([{**{str(k):v for k,v in item},'GRADE_REPRESENTATION':representation,
                              'VALUE_DIMENSION_L_T_M':unit} for item in metadata])
    return {'name':'CANCELLATION_'+name+'_'+key,'body':body,'metadataBody':metadata_body,
            'units':unit_map,'valueUnit':unit,'representation':representation,
            'heavy':key!='CROSS_PRODUCT_RESIDUAL' and not key.endswith('_DOMAIN'),
            'metadata':metadata}
