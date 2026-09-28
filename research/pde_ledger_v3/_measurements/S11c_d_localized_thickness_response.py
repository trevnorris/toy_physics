#!/usr/bin/env python3
"""Fixed-slice thickness forcing transform and outgoing-action candidate.

Prepared for independent method/build review. No launch without a new pinned
approval/gate and the existing whole-job guard plus normalization supervisor.
Restores saved results only; no producer, mode, root, inverse or LU replay.
"""
import argparse
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import pickle
import re
import resource
import signal
import shutil
import time
import traceback

ROOT = Path('/var/projects/toy_physics')
STORE = ROOT / '_scratch/s11c'
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')
GRADES = ((0, 0), (1, 0), (0, 1))


def require(test, message):
    if not test:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def route(path):
    path = Path(path)
    return {'path': str(path), 'canonicalPath': str(path.resolve(strict=True)),
            'bytes': path.stat().st_size, 'sha256': digest(path)}


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


def containment():
    group = next(s[3:] for s in Path('/proc/self/cgroup').read_text().splitlines()
                 if s.startswith('0::'))
    root = Path('/sys/fs/cgroup') / group.lstrip('/')
    result = {k: (root/k).read_text().strip()
              for k in ('memory.max', 'memory.swap.max', 'pids.max')}
    result.update(nice=os.getpriority(os.PRIO_PROCESS, 0),
                  affinity=sorted(os.sched_getaffinity(0)),
                  threads={k: os.environ.get(k) for k in THREADS})
    require(result['memory.max'] == str(2*1024**3) and result['memory.swap.max'] == '0'
            and result['pids.max'] == '32' and result['nice'] >= 15
            and len(result['affinity']) == 1
            and all(v == '1' for v in result['threads'].values()),
            'required whole-job containment absent; no science imported')
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def timeout(*_):
        raise TimeoutError('840-second native limit; preserve incomplete input')
    signal.signal(signal.SIGALRM, timeout); signal.alarm(840)
    return {**result, 'nativeWallSeconds': 840}


class SavedCodec(pickle.Unpickler):
    def find_class(self, module, name):
        if not (module.startswith(('sympy.', 'numpy.'))
                or module in ('sympy', 'numpy', 'builtins', 'collections')):
            raise pickle.UnpicklingError(('unapproved saved class', module, name))
        return super().find_class(module, name)




class Journal:
    """Restore prior returns without calling their functions; journal new work."""
    def __init__(self, out):
        self.out, self.records, self.stack = out, [], []
        self.next_number = 0
        self.restored_parent = self.restored_diagnostic = 0
        self.restored_artifacts = []

    @property
    def active(self):
        return self.stack[-1] if self.stack else None

    def value(self, name, value):
        path = self.out/name; path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            pickle.dump(value, stream, protocol=4)
            stream.flush(); os.fsync(stream.fileno())
        return route(path)

    def copy(self, name, record):
        path = Path(record['path'])
        require(route(path) == record, ('saved bytes changed', path))
        target = self.out/name; target.parent.mkdir(parents=True, exist_ok=True)
        with path.open('rb') as src, target.open('xb') as dest:
            shutil.copyfileobj(src, dest); dest.flush(); os.fsync(dest.fileno())
        result = route(target)
        require(result['sha256'] == record['sha256'], 'copied return identity')
        return result

    def begin(self, name):
        folder = 'operations/%04d-%s' % (self.next_number, name)
        self.next_number += 1
        return folder, {'name': name, 'startedUtc': datetime.now(timezone.utc).isoformat()}

    def restore(self, previous, diagnostic=False):
        name = ('diagnostic-' if diagnostic else '')+previous['name']
        folder, record = self.begin(name)
        record.update(input=self.copy(folder+'/input.pickle', previous['input']),
                      execution=('RESTORED_DIAGNOSTIC_COMPLETE_RETURN' if diagnostic
                                 else 'RESTORED_PRIOR_COMPLETE_RETURN'), priorOperation=previous)
        save(self.out/folder/'started.json', record); self.stack.append(record)
        record['value'] = self.copy(folder+'/value.pickle', previous['value'])
        with Path(record['value']['path']).open('rb') as stream:
            result = SavedCodec(stream).load()
        record['finishedUtc'] = datetime.now(timezone.utc).isoformat()
        save(self.out/folder/'completed.json', record)
        self.records.append(record); self.stack.pop()
        if diagnostic: self.restored_diagnostic += 1
        else: self.restored_parent += 1
        return result

    def artifact(self, name, record, restore=False):
        result = self.copy(name, record)
        self.restored_artifacts.append({'source':record,'copy':result,
                                       'execution':'COPIED_PRIOR_COMPLETE_ARTIFACT'})
        if restore:
            with Path(result['path']).open('rb') as stream:
                return SavedCodec(stream).load()
        return result

    def op(self, name, function, *args):
        folder, record = self.begin(name)
        record.update(input=self.value(folder+'/input.pickle', args),
                      execution='NEW_UNFINISHED_OPERATION')
        save(self.out/folder/'started.json', record); self.stack.append(record)
        result = function(*args)
        record.update(value=self.value(folder+'/value.pickle', result),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json', record)
        self.records.append(record); self.stack.pop()
        return result

    def reuse(self, name, previous, result, args):
        folder, record = self.begin(name)
        record.update(input=self.value(folder+'/input.pickle',args),
                      execution='REUSED_COMPLETE_RETURN_IDENTICAL_OPERANDS',priorOperation=previous)
        save(self.out/folder/'started.json',record);self.stack.append(record)
        record.update(value=self.copy(folder+'/value.pickle',previous['value']),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json',record)
        self.records.append(record);self.stack.pop()
        return result

    def resume_input(self, name, function, previous_input):
        # The failed call never returned. Copy its exact saved input before
        # executing only that unfinished function under the existing guard.
        folder, record = self.begin(name)
        record.update(input=self.copy(folder+'/input.pickle',previous_input),
                      execution='NEW_UNFINISHED_OPERATION_FROM_SAVED_INPUT',
                      priorIncompleteInput=previous_input)
        save(self.out/folder/'started.json',record);self.stack.append(record)
        with Path(record['input']['path']).open('rb') as stream:
            args = SavedCodec(stream).load()
        result = function(*args)
        record.update(value=self.value(folder+'/value.pickle',result),
                      finishedUtc=datetime.now(timezone.utc).isoformat())
        save(self.out/folder/'completed.json',record)
        self.records.append(record);self.stack.pop()
        return result


def construct(spec, out, journal):
    import sympy as sp
    import mpmath as mp
    j = journal
    actual = {name: route(v['path']) for name, v in spec['inputs'].items()}
    save(out/'prehashes.json', actual)
    require(actual == spec['inputs'], 'pinned source bytes changed')
    cache = {}

    def load(name):
        if name not in cache:
            def read(key):
                with Path(actual[key]['path']).open('rb') as stream:
                    return SavedCodec(stream).load()
            cache[name] = j.op('restore-'+name, read, name)
        return cache[name]

    def checked(name, function, *args):
        result = j.op(name, function, *args)
        # Complete scientific return persists before inspecting mathematical flags.
        save(out/(name+'.json'), json_form(result))
        require(all(result['checks'].values()), ('failed saved check', name))
        return result

    def json_form(value):
        if value is None or isinstance(value, (str, bool, int, float)):
            return value
        if isinstance(value, dict):
            return {str(k): json_form(v) for k, v in value.items()}
        if isinstance(value, (tuple, list)):
            return [json_form(v) for v in value]
        if isinstance(value, sp.MatrixBase):
            return {'shape': list(value.shape), 'entries': [
                {'row': i, 'column': h, 'value': sp.sstr(value[i,h])}
                for i in range(value.rows) for h in range(value.cols) if value[i,h] != 0]}
        if value is sp.true or value is sp.false:
            return bool(value)
        return sp.sstr(value)

    candidate, context, units, binding = (load(k) for k in
                                         ('candidate','context','units','binding'))
    kn = context['normalMomentum']; mass = context['fourierMass']
    blocks = {v['index']: v for v in candidate['blocks']}
    require(set(blocks) == {16,17} and mass == 2*sp.pi, 'inherited blocks/Fourier convention')
    require(candidate['fixedPhysicalInput'] == binding['physicalInput'], 'physical input join')
    require(candidate['effectiveReferenceGrades'] == {'eta_bg':0,'sigma_W':0}, 'reference grades')
    ell = sp.Rational(binding['physicalInput']['parameters']['L_W'])
    require(ell.is_positive is True, 'positive saved length')
    r = sp.Dummy('positiveSourceRadical', positive=True)
    T = sp.Dummy('profileT', real=True)
    q = sp.Dummy('dimensionlessTransfer', real=True)
    used_orders = set()
    certificates = {}

    # Fourier transform of sech(x)^2 tanh(x)^n is B0(q)*pn(q).
    # B0 has its actual removable value 2 at q=0. The finite recurrence
    # follows by integration by parts with a vanishing boundary term.
    def basis_polynomial(n, transfer):
        values = [sp.S.One]
        for h in range(n):
            previous = h*values[h-1] if h else sp.S.Zero
            values.append((previous-sp.I*transfer*values[h])/(h+2))
        return values[n]

    def base_transform(transfer):
        return sp.Piecewise((sp.Integer(2), sp.Eq(transfer,0)),
                            (sp.pi*transfer/sp.sinh(sp.pi*transfer/2), True))

    def remove_phase(expression, phase):
        return sp.Add(*(sp.expand_power_exp(sp.powsimp(term*sp.exp(-phase), combine='exp'))
                        for term in sp.Add.make_args(sp.expand_mul(expression))))

    def profile_transform(envelope, coordinate, scale, transfer):
        # No integration engine: only the actually present polynomial profiles.
        converted = envelope.xreplace({sp.tanh(coordinate/scale): T})
        polynomial = sp.Poly(converted, T, domain='EX')
        require(all(not c.has(coordinate, sp.Integral, sp.Limit, sp.Derivative)
                    for c in polynomial.all_coeffs()), 'unrecognized localized source factor')
        quotient, remainder = polynomial.div(sp.Poly(1-T*T,T,domain='EX'))
        rem = [sp.cancel(c) for c in remainder.all_coeffs()]
        endpoint_values = [sp.cancel(polynomial.eval(v)) for v in (-1,1)]
        orders = [powers[0] for powers, c in quotient.terms() if c != 0]
        require(not orders or max(orders) <= 6, 'profile degree outside finite prepared dictionary')
        amplitude = scale*sp.Add(*(c*basis_polynomial(powers[0],transfer)
                                   for powers,c in quotient.terms()))
        return {'checks': {'exactSechSquaredFactor': all(c == 0 for c in rem),
                           'localizedAtBothEnds': all(c == 0 for c in endpoint_values)},
                'polynomial': polynomial.as_expr(), 'quotient': quotient.as_expr(),
                'remainder': rem, 'endpointValues': endpoint_values, 'orders': orders,
                'coordinate': coordinate, 'scale': scale, 'transfer': transfer,
                'amplitudeWithoutB0': amplitude,
                'transform': amplitude*base_transform(transfer),
                'identity': 'FT[sech^2(x)*tanh(x)^n]=B0(q)*p_n(q), exp(-iqx) convention'}

    def transform_local(value, z, k0):
        if value == 0:
            return {'checks': {'zeroSource': True}, 'amplitudeWithoutB0': sp.S.Zero,
                    'orders': [], 'source': value}
        envelope = remove_phase(value, sp.I*k0*z)
        result = profile_transform(envelope,z,ell,ell*(kn-k0))
        result.update(source=value, removedPlanePhase=sp.I*k0*z)
        return result

    def transform_nonlocal(value, z, k0):
        # Two fixed source schemas; every unknown carrier is rejected.
        outer = [v for v in sp.preorder_traversal(value)
                 if isinstance(v,sp.Integral) and len(v.limits) in (2,3)]
        require(len(outer) == 1, 'exactly one native two/three-variable carrier')
        outer = outer[0]
        coefficient = value.xreplace({outer:sp.S.One})
        require(value == coefficient*outer and not coefficient.has(sp.Integral,z),
                'native coefficient extraction')
        require(all(tuple(v[1:]) == (-sp.oo,sp.oo) for v in outer.limits),
                'fixed whole-real native limits')
        variables = tuple(v[0] for v in outer.limits)
        require(not coefficient.has(*variables), 'coefficient independent of every bound variable')
        kout, zs = variables[0], variables[-1]
        require('OutputNormalMomentum' in kout.name and 'SourceNormalPosition' in zs.name,
                'source coordinate roles/order')
        if len(variables) == 2:
            phase = sp.I*kout*z + sp.I*(k0-kout)*zs
            envelope = remove_phase(outer.function,phase).xreplace({kout:kn})
            require(not envelope.has(z), 'output character fully extracted')
            result = profile_transform(envelope,zs,ell,ell*(kn-k0))
            result['amplitudeWithoutB0'] *= mass*coefficient
            schema = 'LOCALIZED_SOURCE_MULTIPLIER'
            delta_factors = 1
        else:
            p = variables[1]
            require('InputNormalMomentum' in p.name, 'input momentum role/order')
            nested = list(outer.function.atoms(sp.Integral))
            require(len(nested) == 1 and len(nested[0].limits) == 1,
                    'one localized profile transform inside the source carrier')
            profile = nested[0]; xi = profile.limits[0][0]
            require(tuple(profile.limits[0][1:]) == (-sp.oo,sp.oo)
                    and 'ProfileCoordinate' in xi.name, 'profile coordinate role')
            phase = sp.I*kout*z + sp.I*(k0-p)*zs
            multiplier = remove_phase(outer.function.xreplace({profile:sp.S.One}),phase)
            require(not multiplier.has(z,zs,xi,sp.Integral), 'plane source leaves a multiplier only')
            profile_envelope = remove_phase(profile.function,-sp.I*ell*(kout-p)*xi)
            require(not profile_envelope.has(kout,p,z,zs), 'localized fixed profile transform')
            input_domain = checked('native-input-multiplier-domain-%04d'%j.next_number,
                                   amplitude_domain,multiplier,p,(kout,))
            result = profile_transform(profile_envelope,xi,sp.S.One,ell*(kn-k0))
            result['inputMultiplierDomain'] = input_domain
            result['amplitudeWithoutB0'] *= mass**2*coefficient*multiplier.xreplace({kout:kn,p:k0})
            schema = 'LOCALIZED_TRANSFER_WITH_SAVED_PLANE_DELTA'
            delta_factors = 2
            result.update(collapsedInputMomentum=p, inheritedPlaneMomentum=k0,
                          sourceMultiplierBeforeCollapse=multiplier)
        result.update(source=value, orderedLimits=outer.limits, schema=schema,
                      fourierDeltaFactors=delta_factors, fourierMass=mass)
        result['transform'] = result['amplitudeWithoutB0']*base_transform(ell*(kn-k0))
        require(not result['amplitudeWithoutB0'].has(z,zs,sp.Integral,sp.Limit,sp.Derivative),
                'no unhandled coordinate or operator in transformed source')
        return result

    def amplitude_domain(amplitude, momentum=None, other_symbols=()):
        momentum = kn if momentum is None else momentum
        # Inspect source multipliers before whole-matrix simplification. This
        # makes the nonzero argument refer to actual factors, not a guess
        # based on the numerator of a combined expression.
        radical_nodes = [v for v in amplitude.atoms(sp.Pow)
                         if v.exp.is_Rational and v.exp.q == 2 and v.has(momentum)]
        radicands = sorted({v.base for v in radical_nodes}, key=sp.default_sort_key)
        require(len(radicands) <= 1, 'multiple momentum radicals outside this fixed native schema')
        substitutions = {}; radical_record = None
        if radicands:
            radicand = radicands[0]; poly = sp.Poly(radicand,momentum)
            a, b = -poly.nth(0), -poly.nth(2)
            require(poly.as_expr() == -a-b*momentum**2 and a.is_positive is True and b.is_positive is True,
                    'saved real-axis evanescent quadratic branch')
            substitutions = {v:(sp.I*r)**(2*v.exp) for v in radical_nodes}
            radical_record = {'sourceRadicand':radicand,'positiveRadicand':a+b*momentum**2,
                              'a':a,'b':b,'sourcePowers':radical_nodes,'branchReplacement':sp.I*r}
        rational = amplitude.xreplace(substitutions)
        require(rational.is_rational_function(momentum,r), 'forcing amplitude rational in momentum and positive radical')
        denominator_records = []
        for node in sorted(rational.atoms(sp.Pow),key=sp.default_sort_key):
            if not node.base.has(momentum,r) or node.exp.is_nonnegative is True:
                continue
            require(node.exp.is_Integer and node.exp.is_negative is True,
                    'only integer inverse powers of smooth source denominators')
            require(not node.base.has(momentum), 'unexpected explicit momentum denominator')
            num, den = sp.fraction(sp.cancel(node.base))
            entries=[]
            for label,part in (('numerator',num),('denominator',den)):
                real = sp.Poly(sp.re(sp.expand_complex(part)),r)
                imag = sp.Poly(sp.im(sp.expand_complex(part)),r)
                certificates_for_part=[]
                for component,poly in (('real',real),('imaginary',imag)):
                    coeffs=[c for c in poly.all_coeffs() if c != 0]
                    positive=bool(coeffs) and all(c.is_positive is True for c in coeffs)
                    negative=bool(coeffs) and all(c.is_negative is True for c in coeffs)
                    certificates_for_part.append({'component':component,'coefficients':coeffs,
                                                   'positive':positive,'negative':negative})
                entries.append({'part':label,'expression':part,'components':certificates_for_part,
                                'nonzeroOnPositiveReal':any(v['positive'] or v['negative'] for v in certificates_for_part)})
            denominator_records.append({'sourceFactor':node,'rationalNumerator':num,'rationalDenominator':den,
                                        'certificates':entries})
        checks={'everyDenominatorNonzero':all(v['nonzeroOnPositiveReal'] for d in denominator_records for v in d['certificates']),
                'noUnboundParameters':not (rational.free_symbols-{momentum,r,*other_symbols}),
                'noRemainingTranscendentals':not rational.has(sp.exp,sp.sinh,sp.cosh,sp.tanh,sp.Integral,sp.Limit)}
        return {'checks':checks,'sourceAmplitude':amplitude,'rationalAmplitude':rational,
                'momentum':momentum,'otherMomentumParameters':other_symbols,
                'radical':radical_record,'denominatorFactors':denominator_records,
                'domainConclusion':'Smooth real amplitude of polynomial growth; with B0(L*(k-k0)), a Schwartz forcing transform.',
                'argument':'Positive-component polynomial denominators have no positive-r zeros and their reciprocal derivatives have polynomial bounds for r>=sqrt(a)>0.'}

    def certificate(amplitude, label):
        if amplitude not in certificates:
            certificates[amplitude] = checked(label,amplitude_domain,amplitude)
        return certificates[amplitude]

    def numerical_matrix(matrix, substitutions):
        return sp.ImmutableMatrix(matrix).applyfunc(lambda v:sp.N(v.subs(substitutions),80))

    def norm(matrix):
        return max([abs(complex(v)) for v in matrix] or [0.0])

    def equation_probe(P, inverse, forcing, point):
        p, inv, f = (numerical_matrix(v,{kn:point}) for v in (P,inverse,forcing))
        response = inv*f; residual = p*response-f
        error = norm(residual)/max(1.0,norm(f))
        return {'checks':{'finiteSourceEquationResidual':error < 1e-10},
                'k':point,'source':p,'savedInverse':inv,'forcing':f,'response':response,
                'residual':residual,'relativeResidual':error,'decimalPrecision':80,
                'scope':'New-source contraction of saved inverse; no solve or inverse reconstruction.'}

    def construct_case(index, f, d, candidate, context, units, binding):
        tag='block%d-grade01'%index
        k0=blocks[index]['k']
        z=next(v for v in f['seed'].free_symbols if v.name=='s11cdNormalPosition')
        source_checks = {'zeroLift':f['lift']==sp.zeros(5),
            'zeroLiftedAction':f['liftedAction']['total']==sp.zeros(5),
            'directSourceSign':f['direct']==-f['forcing']['total'],
            'savedSeedFrame':f['seed']==sp.exp(sp.I*k0*z)*blocks[index]['residue'],
            'domainDirectJoin':d['directForcing']==f['direct'],
            'sourceTermsJoin':d['profileForcingTerms']==f['forcing']['terms'],
            'noProfileAbel':not f['direct'].has(context['profileRegulator']),
            'forcingUnitJoin':f['forcingUnits']==binding['forcingUnits']}
        save(out/(tag+'-source-joins.json'),source_checks)
        require(all(source_checks.values()),'fixed source/end/domain join')
        amplitudes=sp.zeros(5); transformed_terms=[]
        for i in range(5):
            for h in range(5):
                value=-f['forcing']['local'][i,h]
                result=checked(tag+'-local-%d-%d'%(i,h),transform_local,value,z,k0)
                amplitude=result['amplitudeWithoutB0'];used_orders.update(result['orders'])
                certificate(amplitude,tag+'-local-domain-%d-%d'%(i,h))
                amplitudes[i,h]+=amplitude
        for n,term in enumerate(f['forcing']['terms']):
            if term['value']==0:continue
            result=checked(tag+'-native-%03d'%n,transform_nonlocal,-term['value'],z,k0)
            amplitude=result['amplitudeWithoutB0'];used_orders.update(result['orders'])
            certificate(amplitude,tag+'-native-domain-%03d'%n)
            amplitudes[term['row'],term['frameColumn']]+=amplitude
            transformed_terms.append({'address':{k:term[k] for k in ('row','sourceColumn','frameColumn','term')},
                                      'sourceValue':-term['value'],'amplitudeWithoutB0':amplitude})
        amplitude=sp.ImmutableMatrix(amplitudes)
        forcing=amplitude*base_transform(ell*(kn-k0))
        transform={'block':index,'grade':(0,1),'planeMomentum':k0,'source':f['direct'],
                   'amplitudeWithoutB0':amplitude,'forcingTransform':forcing,
                   'transformedNativeTerms':transformed_terms,'sourceChecks':source_checks,
                   'profileZeroTransferValue':2,'fourierMass':mass,
                   'forcingEntryUnits':binding['forcingUnits'],
                   'transformedForcingEntryUnits':[[tuple(a-b for a,b in zip(v,units['spectralMeasureUnit']))
                       for v in row] for row in binding['forcingUnits']]}
        j.value(tag+'-forcing-transform.pickle',transform)
        probes=[checked(tag+'-source-equation-%d'%n,equation_probe,candidate['fixedSymbol'],
                        candidate['fixedInverse'],forcing,point)
                for n,point in enumerate((sp.S.Zero,1/ell,-1/ell))]
        # Source-specific smooth multiplication of the inherited real-axis
        # PV/delta distribution, not evaluation of the pointwise kernel on its diagonal.
        pole_terms=[]
        for pole in candidate['blocks']:
            at_pole=forcing.subs(kn,pole['k'])
            coefficient=pole['deltaCoefficient']*at_pole
            source_at_pole=candidate['fixedSymbol'].subs(kn,pole['k'])
            numeric_residual=numerical_matrix(source_at_pole,{})*numerical_matrix(coefficient,{})
            relative=norm(numeric_residual)/max(1.0,norm(numerical_matrix(coefficient,{})))
            record={'poleIndex':pole['index'],'k':pole['k'],'forcingAtPole':at_pole,
                    'deltaCoefficient':pole['deltaCoefficient'],'responseCoefficient':coefficient,
                    'numericSourceNullResidual':numeric_residual,'relativeResidual':relative}
            j.value(tag+'-pole%d.pickle'%pole['index'],record)
            require(relative<1e-8,'inherited pole term annihilated on this source')
            pole_terms.append(record)
        spectral=candidate['fixedInverse']*forcing
        intervals=candidate['principalValueExclusionIntervals']
        exclusions=set().union(*(a.free_symbols|b.free_symbols for a,b in intervals))
        exclusion=next(v for v in exclusions if v.name=='s11cdOutgoingMomentumExclusion')
        phase=sp.exp(sp.I*kn*z)
        response_pv=sp.ImmutableMatrix(5,5,lambda i,h:sp.Limit(sp.Add(*(
            sp.Integral(phase*spectral[i,h]/mass,(kn,a,b)) for a,b in intervals)),exclusion,0,dir='+'))
        response_delta=sum((sp.exp(sp.I*p['k']*z)*p['responseCoefficient']/mass for p in pole_terms),sp.zeros(5))
        selected=transformed_terms[0]
        changed=forcing.copy().as_mutable()
        changed[selected['address']['row'],selected['address']['frameColumn']] -= selected['amplitudeWithoutB0']*base_transform(ell*(kn-k0))
        movements=[]
        for point in (sp.S.Zero,1/ell,-1/ell):
            difference=numerical_matrix(forcing-sp.ImmutableMatrix(changed),{kn:point})
            response_movement=numerical_matrix(candidate['fixedInverse'],{kn:point})*difference
            movements.append({'k':point,'forcingMovement':difference,'responseMovement':response_movement,
                              'forcingNorm':norm(difference),'responseNorm':norm(response_movement)})
        mutation={'omittedActualNativeTerm':selected,'baselineForcing':forcing,
                  'mutatedForcing':sp.ImmutableMatrix(changed),'probes':movements,
                  'responsive':any(v['forcingNorm']>1e-10 and v['responseNorm']>1e-10 for v in movements)}
        j.value(tag+'-thickness-mutation.pickle',mutation)
        require(mutation['responsive'],'actual thickness-term omission must change transformed response')
        # Contract each field/row/source-column unit using the actual saved tables.
        unit_checks=[]
        for i in range(5):
            for h in range(5):
                for row in range(5):
                    joined=tuple(a+b for a,b in zip(units['inverseEntryUnits'][i][row],binding['forcingUnits'][row][h]))
                    unit_checks.append(joined==tuple(binding['endFieldUnits'][i][h]))
        require(all(unit_checks),'source-specific response unit contraction')
        result={'block':index,'grade':(0,1),'forcingTransform':forcing,'regularSpectralResponse':spectral,
                'principalValueResponse':response_pv,'poleContributions':pole_terms,
                'deltaResponse':response_delta,'response':response_pv+response_delta,
                'fieldEntryUnits':binding['endFieldUnits'],'fourierMass':mass,
                'sourceDomain':'Schwartz forcing from this checked finite localized-profile dictionary and native multipliers.',
                'actionScope':'Candidate smooth-source action of the inherited coupled real-axis PV/delta distribution.',
                'positionIntegralsEvaluated':False,'pointwiseKernelDiagonalValueUsed':False,
                'complexFrequencyRetardedEquivalenceClaimed':False,'currentNormalized':False,'fullFORM':False}
        j.value(tag+'-outgoing-action.pickle',result)
        summary={'block':index,'sourceJoins':source_checks,'nativeTermsTransformed':len(transformed_terms),
                 'profileBasisOrders':sorted(used_orders),'sourceEquationRelativeResiduals':[p['relativeResidual'] for p in probes],
                 'poleRelativeResiduals':[p['relativeResidual'] for p in pole_terms],
                 'actualThicknessMutationResponsive':mutation['responsive'],'allUnitContractionsMatch':all(unit_checks),
                 'forcingSmoothAndRapidlyDecreasingFromSavedCertificates':True,
                 'integralsEvaluated':False,'acceptancePendingSavedOutputInspection':True}
        save(out/(tag+'-summary.json'),summary)
        return summary

    cases=[]
    for index in (17,16):
        tag='block%d-grade01'%index
        cases.append(j.op(tag+'-construct-source-action',construct_case,index,
                          load(tag+'-forcing'),load(tag+'-domain'),candidate,context,units,binding))

    def profile_probe(n,point):
        mp.mp.dps=60
        value=mp.quad(lambda x:mp.sech(x)**2*mp.tanh(x)**n*mp.exp(-1j*point*x),
                      [-mp.inf,-4,0,4,mp.inf])
        exact=basis_polynomial(n,sp.Integer(point))*base_transform(sp.Integer(point))
        expected=complex(sp.N(exact,60));error=abs(complex(value)-expected)
        return {'checks':{'independentProfileTransform':error<1e-10},'order':n,'dimensionlessTransfer':point,
                'quadratureReal':str(mp.re(value)),'quadratureImaginary':str(mp.im(value)),
                'closedForm':exact,'absoluteError':error,'decimalPrecision':60,
                'scope':'Dimensionless dictionary probe, not a new physical input or profile response solve.'}
    for n in sorted(used_orders):
        for point in (0,1,2):
            checked('profile-probe-%d-%d'%(n,point),profile_probe,n,point)
    return {'status':'FIXED_INPUT_THICKNESS_OUTGOING_ACTION_CANDIDATE_BUILT',
            'cases':cases,'profileBasisOrders':sorted(used_orders),'denominatorCertificates':len(certificates),
            'newModeRootOrLUSolves':0,'sourceProducerImports':0,'constructorReplay':False,
            'integratedPositionResponse':False,'newPhysicalInput':False,
            'fullResponseGreenFormA11A12Accepted':False,'newIndependentClearClaimed':False,
            'resultAcceptancePending':True}


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--input-manifest',type=Path,required=True)
    parser.add_argument('--gate-receipt',type=Path,required=True)
    parser.add_argument('--run-directory',type=Path,required=True)
    args=parser.parse_args()
    spec=json.loads(args.input_manifest.read_text());gate=json.loads(args.gate_receipt.read_text())
    require(gate['status']=='READY_FOR_ONE_GUARDED_LOCALIZED_THICKNESS_RESPONSE','method gate incomplete')
    require(gate['workerSha256']==digest(Path(__file__)) and gate['inputManifestSha256']==digest(args.input_manifest),
            'method gate worker/input identity')
    require(gate['scopeExplicitlyApproved'] is True and gate['substantiveReviewFindingsClosed'] is True,
            'approved scope and independent method/build disposition required')
    require(gate['seconds']==900 and gate['nativeSeconds']==840 and gate['automaticRetry'] is False,
            'ordinary bounded stage duration/retry contract')
    for key,path in (('sharedGuardSha256',ROOT/'scripts/s11c_guarded_run.py'),
                     ('supervisorSha256',ROOT/'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py')):
        require(gate[key]==digest(path),'required guard/supervisor identity')
    approval_path=Path(gate['scopeApprovalPath']);disposition_path=Path(gate['reviewDispositionPath'])
    require(gate['scopeApprovalSha256']==digest(approval_path)
            and gate['reviewDispositionSha256']==digest(disposition_path),'approval/disposition identity')
    approval=json.loads(approval_path.read_text())
    require(approval['status']=='AUTHORIZED_ONE_GUARDED_LOCALIZED_THICKNESS_RESPONSE'
            and approval['workerSha256']==digest(Path(__file__))
            and approval['inputManifestSha256']==digest(args.input_manifest),'actual stage approval')
    require(str(args.run_directory.resolve())==gate['resultDirectory'],'new pinned result directory')
    observed=containment()
    out=args.run_directory.resolve();out.relative_to(STORE);out.mkdir(parents=True,exist_ok=False)
    save(out/'native-containment.json',observed);save(out/'input-manifest.json',spec);save(out/'gate-receipt.json',gate)
    journal=Journal(out);started=time.monotonic()
    try:
        result=construct(spec,out,journal)
        post={name:route(v['path']) for name,v in spec['inputs'].items()}
        save(out/'posthashes.json',post);require(post==spec['inputs'],'source posthash mismatch')
        result.update(wallSeconds=time.monotonic()-started,allSourceHashesUnchanged=True,
                      completedOperations=len(journal.records))
        save(out/'operation-index.json',journal.records);save(out/'checks.json',result)
        save(out/'artifact-index.json',{str(p.relative_to(out)):route(p) for p in sorted(out.rglob('*')) if p.is_file()})
        print(json.dumps(result,indent=2,allow_nan=False))
    except BaseException:
        trace=traceback.format_exc();signal.alarm(0)
        post={}
        for name,record in spec['inputs'].items():
            try:post[name]=route(record['path'])
            except OSError as error:post[name]={'path':record['path'],'error':str(error)}
        save(out/'failure-posthashes.json',post)
        save(out/'partial-operation-index.json',journal.records)
        save(out/'failure.json',{'traceback':trace,'completedOperations':len(journal.records),
                                'incompleteOperation':journal.active,'incompleteOperationStack':journal.stack,
                                'wallSeconds':time.monotonic()-started,'automaticRetry':False})
        save(out/'failure-artifact-index.json',{str(p.relative_to(out)):route(p) for p in sorted(out.rglob('*')) if p.is_file()})
        raise


if __name__=='__main__':
    main()
