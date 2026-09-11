#!/usr/bin/env python3
"""Fourier-kernel and path probes of the computed reduced bulk radical.

This instrument consumes a completed, source-pinned one-case engine run.
It emits computed operands/residuals, not a global sheet classification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import resource
import time
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from sympy.core.relational import Relational
from sympy.core.symbol import Str

import S11c_d_mixing_scattering_sympy_audit as engine
from ledger_fold import _restore


class ProbeDimensions(engine.DimensionAnalysis):
    def measure(self, value):
        if isinstance(value, Relational):
            self.equate(self.measure(value.lhs), self.measure(value.rhs))
            return self.zero
        if isinstance(value, sp.logic.boolalg.BooleanFunction):
            for part in value.args:
                self.measure(part)
            return self.zero
        if value.func == sp.meijerg:
            self.equate(self.measure(value.argument), self.zero)
            for part in value.ap + value.bq:
                self.equate(self.measure(part), self.zero)
            return self.zero
        return super().measure(value)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_packet(path, prefix):
    namespace = dict(vars(sp), Str=Str)
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            if line.startswith(prefix):
                return line_number, eval(line.partition(': ')[2], {'__builtins__': {}}, namespace)
    raise ValueError(('missing transcript packet', prefix))


def run():
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-manifest', type=Path, required=True)
    parser.add_argument('--scope-lookup', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--case', default='LAB_HELD_RHO4_CONSTANT')
    args = parser.parse_args()
    root = Path(engine.ROOT)
    manifest = json.loads(args.run_manifest.read_text())
    run_dir = Path(manifest['run_directory'])
    if manifest.get('exit_code') != 0:
        raise ValueError('completed engine run required')
    for key, expected in manifest['source_hashes_before'].items():
        if digest(root/key) != expected or manifest['source_hashes_after'][key] != expected:
            raise ValueError(('run source mismatch', key))
    cache = run_dir/'symbols'/('REFERENCE_'+args.case+'.pickle')
    transcript = run_dir/'preflight.out'
    for path in (cache, transcript):
        if digest(path) != manifest['artifacts'][str(path.relative_to(run_dir))]['sha256']:
            raise ValueError(('run artifact mismatch', str(path)))
    scope = json.loads(args.scope_lookup.read_text())
    if scope['files']['scripts/S11c_c2_exports.py']['sha256'] != digest(root/'scripts/S11c_c2_exports.py'):
        raise ValueError('scope branch bindings do not match the current export')
    input_key = str(args.input.resolve().relative_to(root))
    if digest(args.input) != manifest['source_hashes_before'][input_key]:
        raise ValueError('channel input does not match the engine run')

    full, curl, units, known = pickle.loads(cache.read_bytes())[:4]
    symbols = {s.name: s for s in known if isinstance(s, sp.Symbol)}
    symbols.update({s.name: s for s in full.free_symbols})
    dimensions = ProbeDimensions.__new__(ProbeDimensions)
    dimensions.known = dict(known)
    dimensions.unknown, dimensions.constraints, dimensions.solution = {}, set(), {}
    dimensions.zero = (sp.S.Zero,)*3
    r = SimpleNamespace(symbols=symbols, omega=symbols['omega'], ell=symbols['L_W'])
    engine.PHYSICAL_METADATA = engine.PhysicalMetadata(dimensions, r)
    modes = engine.FullPencilModes(SimpleNamespace(r=r, kn=symbols['s11cdSpectralNormalMomentum']), curl, units)
    _, relation, join = modes.analytic(full)
    specification = json.loads(args.input.read_text())
    parameters = {name: sp.Rational(v) for name, v in specification['parameters'].items()}
    bound = {s: parameters[s.name] for s in relation.free_symbols-{modes.k, modes.q}}
    square = sp.solve(relation, modes.q**2)[0].xreplace(bound)
    complex_k = sp.Dummy('fourierSheetComplexMomentum')
    polynomial = sp.Poly(square.xreplace({modes.k: complex_k}), complex_k)
    points = polynomial.all_roots()
    a_squared = sp.cancel(polynomial.nth(0)/polynomial.nth(2))
    scale = sp.sqrt(a_squared)
    q_unit = dimensions.measure(modes.q)
    k_unit = dimensions.measure(modes.k)

    def numeric(tag, value, unit=lambda path: dimensions.zero):
        engine.emit('FOURIER_SHEET_'+tag, value)
        metadata_body = engine.map_leaves(engine.cas(value), lambda v:
            sp.Piecewise((1,v),(0,True)) if isinstance(v,(Relational,sp.logic.boolalg.BooleanFunction)) else v)
        engine.emit('METADATA_FOURIER_SHEET_'+tag, modes.numeric_metadata(metadata_body, unit))

    numeric('PROVENANCE', {
        'manifestSha256': digest(args.run_manifest), 'cacheSha256': digest(cache),
        'transcriptSha256': digest(transcript), 'inputSha256': digest(args.input),
        'scopeLookupSha256': digest(args.scope_lookup), 'instrumentSha256': digest(Path(__file__)),
        'engineSha256': digest(Path(engine.__file__)), 'case': args.case,
    })
    engine.physical('FOURIER_SHEET_SOURCE_RELATION', relation)
    engine.physical('FOURIER_SHEET_SOURCE_JOIN_RESIDUAL', join,
                    zero_dimensions={(i,): dimensions.zero for i in range(len(join))})
    numeric('BRANCH_POINTS', points, lambda path: k_unit)
    numeric('MOMENTUM_SCALE_SQUARED', a_squared, lambda path: tuple(2*x for x in k_unit))
    numeric('POLYNOMIAL_COEFFICIENTS', polynomial.all_coeffs(),
            lambda path: tuple(2*q-(polynomial.degree()-path[0])*k for q,k in zip(q_unit,k_unit)))
    if polynomial.degree() != 2 or polynomial.nth(1) != 0 or a_squared.is_positive is not True:
        raise NotImplementedError('Fourier-kernel probe requires the even closed-cone chart')

    # Evaluate the actual exported real-axis seed at zero normal momentum.
    seed_rows = []
    for row in scope['closed_rows']:
        for case in row['cases']:
            for source in case['branch_bindings']:
                equation = _restore(source)
                momenta = equation.lhs.args
                mapping = {s: parameters[s.name] for s in equation.rhs.free_symbols-set(momenta)}
                mapping.update(zip(momenta, (parameters['s11cdTangentialMomentum1'],
                                            parameters['s11cdTangentialMomentum2'], sp.S.Zero)))
                seed = equation.rhs.xreplace(mapping)
                seed_rows.append((row['key'], case['case'], str(equation.lhs), seed))
    numeric('EXPORTED_REAL_AXIS_SEEDS', seed_rows,
            lambda path: k_unit if path[-1] == 3 else dimensions.zero)
    # Convert the exported q_out seed to the radical's scale, using the
    # computed coefficient of normal momentum in its dispersion relation.
    q_seed = sp.sqrt(-polynomial.nth(2))*seed_rows[0][3]
    numeric('RADICAL_SEED', q_seed, lambda path: q_unit)
    numeric('SEED_RELATION_RESIDUAL', sp.cancel(q_seed**2-polynomial.nth(0)),
            lambda path: tuple(2*x for x in q_unit))

    p, x = sp.symbols('s11cdFourierSheetP s11cdFourierSheetX', real=True)
    xp, t, z, regulator = sp.symbols('s11cdFourierSheetPositiveX s11cdFourierSheetT '
                                    's11cdFourierSheetZ s11cdFourierSheetRegulator', positive=True)
    beta = sp.Symbol('s11cdFourierSheetBeta', real=True)
    for symbol in (p, x, xp, t, z, regulator, beta):
        dimensions.known[symbol] = dimensions.zero
    normalized_square = sp.cancel(square.subs(modes.k, scale*p)/polynomial.nth(0))
    gamma_integral = sp.Integral(t**sp.Rational(-1, 2)*sp.exp(-z*t), (t, 0, sp.oo))
    gamma_value = gamma_integral.doit()
    gamma_mass = gamma_value.subs(z, 1)
    regulated_ansatz = sp.exp(-regulator*x*x+sp.I*p*x)
    fourier_mass = sp.integrate(sp.integrate(regulated_ansatz, (x, -sp.oo, sp.oo)), (p, -sp.oo, sp.oo))
    heat_factor = gamma_integral.function.subs(z, normalized_square)/gamma_mass
    heat_fourier_integral = sp.Integral(sp.exp(sp.I*p*x)*heat_factor, (p, -sp.oo, sp.oo))/fourier_mass
    heat_fourier = sp.simplify(heat_fourier_integral.doit())
    kernel_integral = sp.Integral(heat_fourier.subs(x, xp), (t, 0, sp.oo))
    kernel = kernel_integral.doit()
    absolute_inner = sp.simplify(sp.integrate(sp.exp(beta*x)*heat_fourier, (x,-sp.oo,sp.oo)))
    decay_coefficient = sp.simplify(sp.diff(sp.log(absolute_inner), t)).limit(t,sp.oo)
    convergence = sp.solve_univariate_inequality(decay_coefficient < 0, beta)
    return_inner = sp.simplify(sp.integrate(sp.exp(-sp.I*p*x)*heat_fourier,(x,-sp.oo,sp.oo)))
    return_transform = sp.integrate(return_inner,(t,0,sp.oo))
    for label,value in (
        ('NORMALIZED_RADICAL_SQUARE', normalized_square), ('GAMMA_OPERANDS',(gamma_integral,gamma_value)),
        ('FOURIER_MASS_OPERANDS',(regulated_ansatz,fourier_mass)),
        ('HEAT_FOURIER_OPERANDS',(heat_fourier_integral,heat_fourier)),
        ('SPATIAL_KERNEL_OPERANDS',(kernel_integral,kernel)),
        ('ABSOLUTE_TRANSFORM_OPERAND',absolute_inner), ('ABSOLUTE_CONVERGENCE',(decay_coefficient,convergence)),
        ('RETURN_TRANSFORM_OPERANDS',(return_inner,return_transform)),
        ('RETURN_RADICAL_RESIDUAL',sp.simplify(return_transform**2*normalized_square-1)),
    ):
        for _, operand in engine.leaves(engine.cas(value)):
            dimensions.equate(dimensions.measure(operand),dimensions.zero)
        numeric(label,value)

    line_number, packet = parse_packet(transcript, 'PY_S11CD_FULL_PENCIL_MODE_INPUT_REFERENCE_')
    candidates = [{str(k):v for k,v in record} for record in packet]
    selected = [(i,c) for i,c in enumerate(candidates)
                if complex(c['K']).real > 0 and 1e-8 < complex(c['K']).imag < float(scale)
                and complex(c['Q']).imag < 0]
    numeric('CANDIDATE_SELECTION', {'sourceLine':line_number,'indices':[i for i,_ in selected]})
    if not selected:
        raise ValueError('no complex candidate in the open Fourier strip matches the recorded diagnostic class')
    selected_index, candidate = selected[0]
    k_end, q_engine = complex(candidate['K']), complex(candidate['Q'])
    square_fn = sp.lambdify(modes.k,square,'numpy')
    relation_bound = relation.xreplace(bound)
    implicit_slope = sp.cancel(-sp.diff(relation_bound,modes.k)/sp.diff(relation_bound,modes.q))
    slope_fn = sp.lambdify((modes.k,modes.q),implicit_slope,'numpy')
    numeric('SELECTED_CANDIDATE', {'index':selected_index,'k':k_end,'q':q_engine,
                                  'enginePhysicalFlag':candidate['PHYSICAL_BULK_SHEET']},
            lambda path: k_unit if path[0]=='k' else q_unit if path[0]=='q' else dimensions.zero)
    numeric('IMPLICIT_PATH_SLOPE', implicit_slope,
            lambda path: tuple(q-k for q,k in zip(q_unit,k_unit)))
    for index,rtol in enumerate((1e-9,1e-12)):
        sol=solve_ivp(lambda u,v:np.asarray([slope_fn(k_end*u,v[0])*k_end]),(0.,1.),
                      np.asarray([complex(q_seed)]),rtol=rtol,atol=rtol*0.01,dense_output=True)
        times=np.linspace(0,1,257)
        values=sol.sol(times)[0]
        result={'rtol':rtol,'status':sol.status,'nfev':sol.nfev,'endQ':values[-1],
                'equationResidual':max(abs(values**2-square_fn(k_end*times))),
                'engineDifference':q_engine-values[-1],'engineSum':q_engine+values[-1]}
        numeric('ODE_TRANSPORT_'+str(index),result,
                lambda path: tuple(2*v for v in q_unit) if path[0]=='equationResidual' else
                q_unit if path[0] in ('endQ','engineDifference','engineSum') else dimensions.zero)

    kernel_fn=sp.lambdify(xp,kernel,'mpmath')
    for index,(digits,tail_exponent) in enumerate(((20,24),(30,34))):
        with mp.workdps(digits):
            point=mp.mpc(str(k_end.real),str(k_end.imag))/mp.mpf(str(scale))
            margin=1-abs(point.imag)
            cutoff=mp.mpf(tail_exponent)/margin
            integrand=lambda v:2*mp.cos(point*v)*kernel_fn(v)
            integral,error=mp.quad(integrand,[0,1,4,12,cutoff],error=True)
            returned_q=mp.mpc(complex(q_seed))/integral
            result={'digits':digits,'cutoff':float(cutoff),'stripMargin':float(margin),
                    'transform':complex(integral),'quadratureError':float(error),
                    'returnedQ':complex(returned_q),'equationResidual':complex(returned_q**2-square_fn(k_end)),
                    'engineDifference':complex(q_engine-returned_q),'engineSum':complex(q_engine+returned_q)}
            numeric('KERNEL_QUADRATURE_'+str(index),result,
                    lambda path: tuple(2*v for v in q_unit) if path[0]=='equationResidual' else
                    q_unit if path[0] in ('returnedQ','engineDifference','engineSum') else dimensions.zero)

    # Two paths through the same local (omega,k) rectangle. On one path the
    # frequency leg starts at real k=0, where S11b supplies its continuation.
    # Neither leg reselects its root by decay or a complex-root half-plane.
    family = relation.xreplace({s:v for s,v in bound.items() if s != r.omega})
    slopes = [sp.cancel(-sp.diff(family,v)/sp.diff(family,modes.q)) for v in (r.omega,modes.k)]
    velocity = sp.lambdify((r.omega,modes.k,modes.q),slopes,'numpy')
    family_fn = sp.lambdify((r.omega,modes.k,modes.q),family,'numpy')
    omega_start = complex(parameters['omega'])
    def segment(start,end,initial):
        delta=[b-a for a,b in zip(start,end)]
        def flow(u,v):
            point=[a+u*b for a,b in zip(start,delta)]
            derivatives=velocity(*point,v[0])
            return np.asarray([sum(a*b for a,b in zip(derivatives,delta))])
        result=solve_ivp(flow,(0,1),np.asarray([initial]),rtol=1e-12,atol=1e-14,dense_output=True)
        grid=np.linspace(0,1,257)
        vals=result.sol(grid)[0]
        residuals=[abs(family_fn(*(a+u*b for a,b in zip(start,delta)),q)) for u,q in zip(grid,vals)]
        return complex(vals[-1]), {'status':result.status,'nfev':result.nfev,
                                  'maximumEquationResidual':max(residuals)}
    for index,direction in enumerate((-1,1)):
        omega_end=omega_start*(1+direction*0.01j)
        base=(omega_start,0j)
        finish=(omega_end,k_end)
        paths=((base,(omega_end,0j),finish),(base,(omega_start,k_end),finish))
        outputs=[]
        for vertices in paths:
            initial=complex(q_seed)
            segments=[]
            for a,b in zip(vertices,vertices[1:]):
                initial,record=segment(a,b,initial)
                segments.append(record)
            outputs.append({'endQ':initial,'segments':segments})
        value={'omegaStart':omega_start,'omegaEnd':omega_end,'endK':k_end,
               'paths':outputs,'pathDifference':outputs[0]['endQ']-outputs[1]['endQ']}
        numeric('FREQUENCY_MOMENTUM_RECTANGLE_'+str(index),value,
                lambda path: k_unit if path[0]=='endK' else
                q_unit if path[0] in ('omegaStart','omegaEnd','pathDifference') or 'endQ' in path else
                tuple(2*v for v in q_unit) if 'maximumEquationResidual' in path else dimensions.zero)
    numeric('DIMENSION_CONSTRAINTS',sorted(dimensions.constraints,key=str))
    numeric('RESOURCE_MEASUREMENTS',{'wallSeconds':time.monotonic()-started,
                                     'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
            lambda path: (0,1,0) if path[0]=='wallSeconds' else dimensions.zero)
    numeric('PROCESS_COMPLETION',True)


if __name__ == '__main__':
    run()
