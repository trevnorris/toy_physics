#!/usr/bin/env python3
"""S11c-d construction workspace; authority: S11c_d_SHARED_PHYSICS.md v10.

Implements a resumable prefix of the reduced scattering construction.
This is not a completed scattering engine. No downstream ledger is published
until the scattering/current/spectral constructions are implemented.

The three mandatory script clauses:
1. The script may PRINT computed objects. It may NOT state conclusions.
2. PRINT the residual; do NOT assert it.
3. Interpretation belongs to the STEP RECORD.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import hashlib
import inspect
import argparse
import json
import pickle
import faulthandler
import os
import re
import resource
import signal
import sys
import time

# Exact, locally computed elimination coefficients can exceed Python's default
# decimal rendering limit. Their transcript representation remains fingerprinted.
if hasattr(sys,'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

import sympy as sp
import numpy as np
from sympy.core.function import AppliedUndef
from sympy.core.symbol import Str

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ledger_fold import load_model, check_consumer, assert_lookups_equal_manifest
from S11c_d_output_codec import PayloadEncoder, emission_index

HERE = Path(__file__).resolve()
ROOT = HERE.parent.parent
CLOSED_KEYS = ('s11cc2ClosedSlabOperator', 's11cc2ClosedCouplingKernel')
IMPORT_KEYS = CLOSED_KEYS + ('L_W', 'omega', 'energy_basis_variable')
# Direct dimension operands; the original root/constant lookups remain intact.
DIMENSION_CARRIERS = (
    's11cc2Coefficientw1Profile', 's11cc2Coefficientm1Profile',
    's11cc2Fieldtheta', 's11cc2FieldeW',
    *(f's11cc2Fieldu{i}' for i in range(1, 4)),
    's11cc2FourierW1ProfileHatTransfer',
    *(f's11cc2FourierW1ProfileJetHat{i}' for i in range(1, 4)),
    *(f's11cc2MiddleMomentum{i}' for i in range(1, 4)),
    's11cc2OutgoingNormalMomentum', 's11cc2Time',
    *(f's11cc2{c}{i}' for c in ('X', 'Y') for i in range(1, 4)),
    *(f's11cc2{c}{f}' for c in ('Trial', 'Test')
      for f in ('A0', 'A1', 'A2', 'Theta', 'E', 'Phi')),
)
IMPORT_KEYS += tuple(k for c in DIMENSION_CARRIERS for k in (c, c+'Dimension'))
BUILD_INPUT_PATHS = (
    HERE, ROOT / 'scripts/S11c_b_exports.py',
    ROOT / 'scripts/S11c_c1_exports.py', ROOT / 'scripts/S11c_c2_exports.py',
    ROOT / 'directives/S11c_d_SHARED_PHYSICS.md', ROOT / 'scripts/ledger_fold.py',
    ROOT / 'directives/S11b_SHARED_PHYSICS.md',
    ROOT / 'scripts/S11c_d_output_codec.py',
)
EMISSION_LINES = {}
PAYLOAD_ENCODER = PayloadEncoder()
PHYSICAL_METADATA = None


def cas(value):
    if isinstance(value, str):
        return Str(value)
    if isinstance(value, dict):
        return sp.Tuple(*(sp.Tuple(cas(k), cas(v)) for k, v in value.items()))
    if isinstance(value, (tuple, list, set, frozenset)):
        return sp.Tuple(*(cas(v) for v in value))
    return sp.sympify(value)


def emit(name, value):
    tag = 'PY_S11CD_' + name
    if tag in EMISSION_LINES:
        raise ValueError(('duplicate tag', tag))
    EMISSION_LINES[tag] = inspect.currentframe().f_back.f_lineno
    print(tag + ': ' + PAYLOAD_ENCODER.encode(sp.srepr(cas(value))), flush=True)


def named(value, key):
    return next(v for k, v in value if str(k) == key)


def leaves(value, path=()):
    """Visit physical leaves without treating association labels as operands."""
    if isinstance(value, (sp.Tuple, tuple, list, sp.MatrixBase)):
        association = all(isinstance(x, (sp.Tuple, tuple)) and len(x) == 2
                          and isinstance(x[0], Str) for x in value)
        for i, item in enumerate(value):
            if association:
                yield from leaves(item[1], path+(str(item[0]),))
            else:
                yield from leaves(item, path+(i,))
    else:
        yield path, value


def map_leaves(value, operation):
    if isinstance(value, sp.Tuple):
        return sp.Tuple(*(map_leaves(x, operation) for x in value))
    if isinstance(value, Str):
        return value
    return operation(value)


def payload_units(value, dimensions, path=()):
    """Preserve supplied tensor units even when a computed map vanishes."""
    if isinstance(dimensions, sp.MatrixBase):
        return {path: tuple(dimensions)}
    out = {}
    association = isinstance(value, sp.Tuple) and all(isinstance(x, sp.Tuple) and len(x)==2
                                                    and isinstance(x[0], Str) for x in value)
    for i, (a,b) in enumerate(zip(value, dimensions)):
        if association:
            out.update(payload_units(a[1], b[1], path+(str(a[0]),)))
        elif not isinstance(a, Str):
            out.update(payload_units(a,b,path+(i,)))
    return out


def memo_xreplace(value, replacements):
    """SymPy xreplace semantics with shared subexpressions visited once."""
    @lru_cache(maxsize=None)
    def visit(node):
        if node in replacements:
            return replacements[node]
        if not node.args:
            return node
        args = tuple(visit(a) for a in node.args)
        return node if args == node.args else node.func(*args)
    return visit(value)


@lru_cache(maxsize=None)
def dag_free_symbols(value):
    """Free symbols with lexical binding, without repeated integral renaming."""
    if isinstance(value, sp.Symbol):
        return frozenset((value,))
    if isinstance(value, sp.Integral):
        result = set(dag_free_symbols(value.function))
        for limit in value.limits:
            if len(limit) == 3:
                result.discard(limit[0])
                result.update(dag_free_symbols(limit[1]))
                result.update(dag_free_symbols(limit[2]))
            else:
                # No inherited integral has an indefinite limit. Keep SymPy's
                # semantics if a future ansatz introduces one.
                return frozenset(value.free_symbols)
        return frozenset(result)
    if isinstance(value, sp.Limit):
        body, variable, endpoint, _ = value.args
        return (dag_free_symbols(body)-{variable}) | dag_free_symbols(endpoint)
    if isinstance(value, sp.Subs):
        return (dag_free_symbols(value.expr)-set(value.variables)) | frozenset().union(
            *(dag_free_symbols(p) for p in value.point))
    return frozenset().union(*(dag_free_symbols(a) for a in value.args))


class DimensionAnalysis:
    """Infer missing input units, then check reduced expressions independently.

    Input carrier dimensions and the input row dimension slots are supplied.
    Additions, phases, derivatives and measures supply linear unit constraints;
    no computed object is assigned its anticipated output dimension.
    """

    def __init__(self, rows, reduction):
        self.known = {}
        self.unknown = {}
        self.constraints = set()
        self.solution = {}
        self.zero = (sp.S.Zero,)*3
        for key in DIMENSION_CARRIERS:
            self.known[rows[key]['value']] = tuple(rows[key+'Dimension']['value'])
        # Fundamental physical-unit declarations in the supplied setup.
        primitives = {'L_W': (1, 0, 0), 'W_0': (1, 0, 0),
                      'omega': (0, -1, 0), 'rho_br': (-3, 0, 1),
                      'rho_m': (-4, 0, 1), 'epsilon_shape': (0, 0, 0),
                      'eta_bg': (0, 0, 0), 'sigma_W': (0, 0, 0)}
        for name, unit in primitives.items():
            self.known[reduction.symbols[name]] = tuple(map(sp.Integer, unit))
        for group in reduction.momentum_groups:
            for atom in group:
                self.known[atom] = tuple(-x for x in self.known[reduction.x[0]])
        for key in CLOSED_KEYS:
            for _, payload in rows[key]['value']:
                # Dimension matrices are leaves of the physical tree, not three
                # unrelated physical quantities.
                def anchor(value, dims):
                    if isinstance(dims, sp.MatrixBase):
                        self.equate(self.measure(value), tuple(dims))
                    elif isinstance(value, sp.Tuple):
                        for a, b in zip(value, dims):
                            if isinstance(a, Str):
                                continue
                            anchor(a, b)
                anchor(named(payload, 'VALUE'), named(payload, 'DIMENSION_L_T_M'))
                for slot in ('COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS'):
                    for eq in named(payload, slot):
                        self.equate(self.measure(eq.lhs), self.measure(eq.rhs))
        equations = sorted(self.constraints, key=sp.default_sort_key)
        variables = sorted(set().union(*(e.free_symbols for e in equations)), key=sp.default_sort_key)
        solved = sp.linsolve(equations, variables)
        emit('INPUT_DIMENSION_CONSTRAINTS', equations)
        emit('INPUT_DIMENSION_SOLVE', solved)
        if solved is sp.S.EmptySet:
            raise ValueError('inconsistent supplied dimension constraints')
        if variables:
            self.solution = dict(zip(variables, next(iter(solved))))
        self.known.update({a: tuple(d.xreplace(self.solution) for d in ds)
                           for a, ds in self.unknown.items()})
        # Units of new coordinates/functions follow their actual ansatz map.
        for a, b in reduction.normal_map.items():
            self.known[b] = self.measure(a)
        for q in reduction.tangents:
            self.known[q] = self.measure(reduction.momentum_groups[0][0])
        for f in reduction.profiles.values():
            self.known[f] = self.measure(rows['s11cc2Coefficientw1Profile']['value'](reduction.x[0]))
        self.known[reduction.xi] = self.zero
        for a in (reduction.regulator, reduction.transfer):
            self.known[a] = self.zero
        for a in reduction.regulated_ansatz.free_symbols:
            self.known[a] = self.zero
        self.known[sp.Function('s11cdSchwartzTest')] = self.zero
        for f in reduction.end_values.values():
            self.known[f] = self.zero
        for old, new in reduction.functions.items():
            self.known[new] = self.known[old]
        # Pre-register field images used by the later pencil construction.
        for key in DIMENSION_CARRIERS:
            old = rows[key]['value']
            if key.startswith(('s11cc2Field', 's11cc2Trial', 's11cc2Test')):
                self.known[sp.Function(key.replace('s11cc2', 's11cdReduced', 1))] = self.known[old]
        self.measure.cache_clear()
        self.constraints.clear()

    def equate(self, a, b):
        if a is not None and b is not None:
            self.constraints.update(sp.expand(x-y) for x, y in zip(a, b) if x != y)

    @lru_cache(maxsize=None)
    def measure(self, value):
        if value == 0:
            return None  # The zero map is polymorphic; never assign units to it.
        if value in self.known:
            return self.known[value]
        if isinstance(value, sp.Symbol) or isinstance(value, AppliedUndef):
            atom = value.func if isinstance(value, AppliedUndef) else value
            if atom in self.known:
                return self.known[atom]
            if isinstance(value, sp.Symbol) and re.fullmatch(r'[wm]1_profile(?:_d[123](?:d[123])*)?', value.name):
                return self.zero
            if atom not in self.unknown:
                i = len(self.unknown)
                self.unknown[atom] = sp.symbols(f's11cdUnit{i}L s11cdUnit{i}T s11cdUnit{i}M')
            return self.unknown[atom]
        if value.is_number or isinstance(value, (Str, sp.logic.boolalg.BooleanAtom)):
            return self.zero
        if isinstance(value, sp.Add):
            ds = [d for a in value.args if (d := self.measure(a)) is not None]
            for d in ds[1:]:
                self.equate(d, ds[0])
            return ds[0] if ds else None
        if isinstance(value, sp.Mul):
            ds = [self.measure(a) for a in value.args]
            return tuple(sum(d[i] for d in ds if d is not None) for i in range(3))
        if isinstance(value, sp.Pow):
            self.equate(self.measure(value.exp), self.zero)
            base = self.measure(value.base)
            if value.exp.is_number:
                return tuple(value.exp*d for d in base)
            self.equate(base, self.zero)
            return self.zero
        if isinstance(value, sp.Derivative):
            d = self.measure(value.expr)
            return tuple(d[i]-sum(n*self.measure(v)[i] for v, n in value.variable_count) for i in range(3))
        if isinstance(value, sp.Subs):
            if isinstance(value.expr, sp.Derivative):
                bound = dict(zip(value.variables, value.point))
                d = self.measure(value.expr.expr)
                return tuple(d[i]-sum(n*self.measure(bound.get(v, v))[i]
                                     for v, n in value.expr.variable_count) for i in range(3))
            for a, b in zip(value.variables, value.point):
                self.equate(self.measure(a), self.measure(b))
            return self.measure(value.expr)
        if isinstance(value, sp.Integral):
            d = self.measure(value.function)
            return tuple(d[i]+sum(self.measure(l[0])[i] for l in value.limits) for i in range(3))
        if isinstance(value, sp.Limit):
            return self.measure(value.args[0])
        if value.func == sp.DiracDelta:
            order = value.args[1] if len(value.args) > 1 else 0
            return tuple(-(order+1)*d for d in self.measure(value.args[0]))
        if value.func in (sp.exp, sp.sin, sp.cos, sp.tanh, sp.log, sp.Heaviside, sp.erf, sp.erfc):
            self.equate(self.measure(value.args[0]), self.zero)
            return self.zero
        if value.func in (sp.conjugate, sp.re, sp.im, sp.Abs):
            return self.measure(value.args[0])
        if value.func == sp.sign:
            return self.zero
        if isinstance(value, sp.Piecewise):
            ds = [self.measure(a.expr) for a in value.args]
            for d in ds[1:]:
                self.equate(d, ds[0])
            return ds[0]
        raise NotImplementedError(('dimension node', value.func))


class PhysicalMetadata:
    def __init__(self, dimensions, reduction):
        self.dimensions, self.reduction = dimensions, reduction
        self.generators = tuple(reduction.symbols[n] for n in ('epsilon_shape', 'eta_bg', 'sigma_W'))

    @lru_cache(maxsize=None)
    def coefficients(self, value):
        """Polynomial arithmetic only in the independent grade bookkeepers."""
        z = (0, 0, 0)
        if not value.has(*self.generators):
            return {z: value} if value != 0 else {}
        if value in self.generators:
            g = [0]*3; g[self.generators.index(value)] = 1
            return {tuple(g): sp.S.One}
        if value.is_Add:
            out = {}
            for a in value.args:
                for g, c in self.coefficients(a).items():
                    out[g] = out.get(g, sp.S.Zero)+c
            return {g: c for g, c in out.items() if c != 0}
        if value.is_Mul or (value.is_Pow and value.exp.is_Integer and value.exp >= 0):
            factors = value.args if value.is_Mul else (value.base,)*int(value.exp)
            out = {z: sp.S.One}
            for a in factors:
                product = {}
                for g, c in out.items():
                    for h, d in self.coefficients(a).items():
                        degree = tuple(x+y for x, y in zip(g, h))
                        product[degree] = product.get(degree, sp.S.Zero)+c*d
                out = product
            return {g: c for g, c in out.items() if c != 0}
        if isinstance(value, (sp.Integral, sp.Derivative, sp.Limit, sp.Subs)):
            operand = value.function if isinstance(value, sp.Integral) else value.args[0]
            out = {}
            for g, c in self.coefficients(operand).items():
                transformed = value.func(c, *value.args[1:])
                if isinstance(value, sp.Derivative):
                    transformed = transformed.doit()
                if transformed != 0:
                    out[g] = transformed
            return out
        if isinstance(value, sp.Piecewise):
            branches = [(self.coefficients(expr), condition) for expr, condition in value.args]
            support = set().union(*(set(c) for c, _ in branches))
            return {g: sp.Piecewise(*((c.get(g, sp.S.Zero), condition)
                                     for c, condition in branches)) for g in support}
        raise NotImplementedError(('nonpolynomial grade dependence', value.func))

    def record(self, value, zero_dimensions=None):
        records = []
        for path, expression in leaves(value):
            if isinstance(expression, Str):
                continue
            ds = self.dimensions.measure(expression)
            if ds is None and zero_dimensions is not None:
                ds = zero_dimensions.get(path)
            coefficients = self.coefficients(expression)
            # Combining eta and sigma orders is performed after the actual
            # homotopy substitution in the coefficient algebra.
            homotopy = {}
            for (e, a, b), c in coefficients.items():
                g = (e, a+b)
                homotopy[g] = homotopy.get(g, sp.S.Zero)+c*(self.reduction.symbols['W_0']/self.reduction.ell)**b
            orders = tuple(g for g, c in sorted(homotopy.items()) if c != 0)
            records.append((path, {'DIMENSION_L_T_M': ds if ds is not None else Str('ZERO_MAP'),
                                   'MULTIGRADE': sorted(coefficients), 'EPSILON_LAMBDA_SUPPORT': orders}))
        return cas(records)


def physical(name, value, *, operands=None, zero_dimensions=None):
    # Keep both operands, but encode large trees with the same carrier-first
    # PIT used for the pencil. This changes serialization, not reduction.
    body = cas(value)
    emit(name, carrier_fingerprint(body) if dag_size(body) > 1200 else body)
    emit('METADATA_'+name, PHYSICAL_METADATA.record(cas(value if operands is None else operands), zero_dimensions))


def fingerprinted(name, value, zero_dimensions=None):
    emit(name, carrier_fingerprint(value))
    emit('METADATA_'+name, PHYSICAL_METADATA.record(value, zero_dimensions))


def dag_size(value):
    seen = set()
    def visit(node):
        if node in seen:
            return
        seen.add(node)
        for child in node.args:
            visit(child)
    visit(value)
    return len(seen)


def dag_substitute(value, mapping):
    """Substitute a field ansatz with derivative/integral linearity retained."""
    @lru_cache(maxsize=None)
    def visit(node):
        if isinstance(node, AppliedUndef) and node.func in mapping:
            return mapping[node.func](*node.args)
        if not node.args or isinstance(node, Str):
            return node
        if isinstance(node, sp.Derivative):
            return sp.diff(visit(node.expr), *node.variable_count)
        if isinstance(node, sp.Integral):
            operand = visit(node.function)
            return sp.Integral(operand, *node.limits) if operand != 0 else operand
        return node.func(*(visit(a) for a in node.args))
    return visit(value)


class ReducedPencil:
    """Strong reduced action and its weak sector extraction.

    Columns are actions on arbitrary normal probe functions, retaining all
    local derivatives and nonlocal normal integrals. They are not a local
    matrix symbol, asymptotic modes, or a solved scattering problem.
    """

    def __init__(self, slab, kernel, reduction):
        self.slab, self.kernel, self.r = slab, kernel, reduction
        self.fields = tuple(sp.Function('s11cdReducedField'+name)
                            for name in ('u1', 'u2', 'u3', 'theta', 'eW'))
        self.strong = sp.Tuple(*named(slab, 'U'), named(slab, 'THETA'), named(slab, 'E_W'))
        self.probes = tuple(sp.Function('s11cdPencilProbe'+str(i)) for i in range(len(self.fields)))
        phase = reduction.wave_phase(reduction.x)
        self.tangent_derivatives = tuple(sp.cancel(sp.diff(phase, x)/phase) for x in reduction.x[:2])
        dual_phase = 1/phase
        self.dual_tangent_derivatives = tuple(sp.cancel(sp.diff(dual_phase, x)/dual_phase)
                                              for x in reduction.x[:2])
        units = PHYSICAL_METADATA.dimensions
        for f, g in zip(self.fields, self.probes):
            units.known[g] = units.known[f]

    def dx(self, expression, i):
        return self.tangent_derivatives[i]*expression if i < 2 else sp.diff(expression, self.r.z)

    def weak_projection(self, vector, sector):
        """Construct test virtual work, then compute its normal Euler derivative.

        The test gradient/curl is an ansatz. Tangential derivatives are obtained
        from the dual harmonic character; normal adjoints follow by varying the
        explicit weak action, rather than inserting projection signs.
        """
        z = self.r.z
        labels = ('Phi',) if sector == 'LONGITUDINAL' else ('A0', 'A1', 'A2')
        tests = [sp.Function('s11cdReducedTest'+name)(z) for name in labels]
        def dual_dx(value, i):
            return self.dual_tangent_derivatives[i]*value if i < 2 else sp.diff(value, z)
        if sector == 'LONGITUDINAL':
            displacement = [dual_dx(tests[0], i) for i in range(3)]
        else:
            displacement = [dual_dx(tests[(i+2)%3], (i+1)%3)
                            -dual_dx(tests[(i+1)%3], (i+2)%3) for i in range(3)]
        action = sum(a*b for a, b in zip(displacement, vector))
        euler = [sp.diff(action, t)-sp.diff(sp.diff(action, sp.diff(t, z)), z) for t in tests]
        return sum(t*e for t, e in zip(tests, euler))

    def trial_ansatz(self, sector):
        def displacement(i, z):
            if sector == 'TRANSVERSE':
                potentials = [sp.Function('s11cdReducedTrialA'+str(j))(z) for j in range(3)]
                def derivative(value, j):
                    return self.tangent_derivatives[j]*value if j < 2 else sp.diff(value, z)
                return derivative(potentials[(i+2)%3], (i+1)%3)-derivative(potentials[(i+1)%3], (i+2)%3)
            if sector == 'LONGITUDINAL':
                phi = sp.Function('s11cdReducedTrialPhi')(z)
                return self.tangent_derivatives[i]*phi if i < 2 else sp.diff(phi, z)
            return sp.S.Zero
        result = {self.fields[i]: lambda z, i=i: displacement(i, z) for i in range(3)}
        for f, active, name in zip(self.fields[3:], ('THETA', 'E_W'), ('Theta', 'E')):
            result[f] = lambda z, active=active, name=name: (
                sp.Function('s11cdReducedTrial'+name)(z) if sector == active else sp.S.Zero)
        return result

    def columns(self):
        epsilon = self.r.symbols['epsilon_shape']
        columns = []
        for j, probe in enumerate(self.probes):
            mapping = {f: (lambda z, i=i: probe(z) if i == j else sp.S.Zero)
                       for i, f in enumerate(self.fields)}
            column = dag_substitute(self.strong, mapping)
            column = map_leaves(column, lambda e: sp.diff(e, epsilon))
            columns.append(column)
        return sp.Tuple(*columns)

    def off_diagonal(self):
        restricted = {sector: dag_substitute(self.slab, self.trial_ansatz(sector))
                      for sector in ('TRANSVERSE', 'THETA', 'E_W', 'LONGITUDINAL')}
        transverse = restricted['TRANSVERSE']
        z = self.r.z
        test = lambda name: sp.Function('s11cdReducedTest'+name)(z)
        forward = {
            'THETA': test('Theta')*named(transverse, 'THETA'),
            'E_W': test('E')*named(transverse, 'E_W'),
            'DIV_U': self.weak_projection(named(transverse, 'U'), 'LONGITUDINAL'),
        }
        reverse = {}
        for sector, label in (('THETA', 'THETA'), ('E_W', 'E_W'), ('LONGITUDINAL', 'DIV_U')):
            vector = named(restricted[sector], 'U')
            reverse[label] = self.weak_projection(vector, 'TRANSVERSE')
        return cas({'TRANSVERSE_TO_THICKNESS': forward, 'THICKNESS_TO_TRANSVERSE': reverse})

    def sector_blocks(self):
        """Weak 2-by-2 sector action, retaining the canonical mixed blocks.

        Each entry is a bilinear differential/integral operator on arbitrary
        trial and test potentials, not a multiplication matrix. The TT and HH
        entries are varied/projected from the reduced slab. Mixed entries are
        the reduced kernel; off_diagonal() remains a separate extraction.
        """
        transverse = dag_substitute(self.slab, self.trial_ansatz('TRANSVERSE'))
        tt = self.weak_projection(named(transverse, 'U'), 'TRANSVERSE')
        hh = {}
        for sector, label in (('THETA', 'THETA'), ('E_W', 'E_W'), ('LONGITUDINAL', 'DIV_U')):
            action = dag_substitute(self.slab, self.trial_ansatz(sector))
            hh[label] = {
                'THETA': sp.Function('s11cdReducedTestTheta')(self.r.z)*named(action, 'THETA'),
                'E_W': sp.Function('s11cdReducedTestE')(self.r.z)*named(action, 'E_W'),
                'DIV_U': self.weak_projection(named(action, 'U'), 'LONGITUDINAL'),
            }
        return cas({'TT': tt, 'TH': named(self.kernel, 'THICKNESS_TO_TRANSVERSE'),
                    'HT': named(self.kernel, 'TRANSVERSE_TO_THICKNESS'), 'HH': hh})


def tree_difference(left, right):
    if isinstance(left, sp.Tuple):
        if not isinstance(right, sp.Tuple) or len(left) != len(right):
            raise ValueError(('operand tree mismatch', left.func, right.func))
        return sp.Tuple(*(tree_difference(a, b) for a, b in zip(left, right)))
    if isinstance(left, Str):
        if left != right:
            raise ValueError(('operand label mismatch', left, right))
        return left
    return left-right


def carrier_fingerprint(value, samples=3):
    """Exact numeric PIT of the scalar carrier algebra, plus symbolic SHA.

    Nonlocal integrals/derivatives/functions are independent shared carriers.
    This does not numerically integrate the operator or decide equality of
    different integral representations. The carrier map is emitted explicitly.
    """
    @lru_cache(maxsize=None)
    def digest(node, bound=()):
        h = hashlib.sha256()
        if isinstance(node, sp.Dummy) and node in bound:
            # Lexical binding identity, not SymPy's process-random dummy_index.
            position = max(i for i, variable in enumerate(bound) if variable == node)
            h.update(b'BOUND_DUMMY')
            h.update(str(len(bound)-1-position).encode())
            h.update(str(sorted(node.assumptions0.items())).encode())
            return h.hexdigest()
        h.update(str(node.func).encode())
        if isinstance(node, sp.Subs):
            inner = bound+tuple(v for v in node.variables if isinstance(v, sp.Dummy))
            for a, scope in ((node.expr, inner), (node.variables, inner), (node.point, bound)):
                h.update(bytes.fromhex(digest(a, scope)))
        elif not node.args:
            h.update(sp.srepr(node).encode())
        else:
            for a in node.args:
                h.update(bytes.fromhex(digest(a, bound)))
        return h.hexdigest()
    carriers = {}
    def atom(node, index):
        key = digest(node)
        carriers[key] = node
        # Deterministic rational points: distinct samples and live inputs.
        seed = hashlib.sha256((key+':'+str(index)).encode()).digest()
        return sp.Rational(int.from_bytes(seed[:2], 'big') % 89+2,
                           int.from_bytes(seed[2:4], 'big') % 83+3)
    @lru_cache(maxsize=None)
    def evaluate(node, index):
        if node.is_number:
            return node
        if node.is_Add:
            return sum(evaluate(a, index) for a in node.args)
        if node.is_Mul:
            return sp.Mul(*(evaluate(a, index) for a in node.args))
        if node.is_Pow and node.exp.is_Integer:
            return evaluate(node.base, index)**node.exp
        return atom(node, index)
    result = tuple((path, digest(e), tuple(sp.N(evaluate(e, i), 24) for i in range(samples)))
                   for path, e in leaves(value))
    # Store the digest of each symbolic carrier, its head, and its free inputs;
    # huge integral operands already occur in the reduced-row output records.
    carrier_record = tuple((k, str(v.func), sorted(dag_free_symbols(v), key=sp.default_sort_key))
                           for k, v in sorted(carriers.items()))
    return cas({'OBJECT_SHA_AND_NUMERIC_PIT': result, 'CARRIER_IDENTITIES': carrier_record})


def bind(fold):
    rows = {key: fold[key] for key in IMPORT_KEYS}
    for key in CLOSED_KEYS:
        provenance = (rows[key]['class'], rows[key]['step'])
        if provenance != ('DERIVED', 'S11c-c2'):
            raise ValueError((key, provenance))
    return rows


def census(value):
    counts = Counter()
    hats = set()
    integrals = set()
    branches = set()
    for node in sp.preorder_traversal(value):
        if isinstance(node, sp.Integral):
            counts['integral_occurrences'] += 1
            integrals.add(node)
        if isinstance(node, sp.DiracDelta):
            counts['delta_occurrences'] += 1
        if isinstance(node, AppliedUndef):
            name = node.func.__name__
            if name.startswith('s11cc2Fourier'):
                counts[name] += 1
                hats.add(node)
            elif name == 's11cc2OutgoingNormalMomentum':
                branches.add(node)
    return {
        'occurrences': dict(sorted(counts.items())),
        'hat_arguments': tuple(sorted(hats, key=sp.default_sort_key)),
        'integral_limits': tuple(sorted({i.limits for i in integrals}, key=str)),
        'distinct_integrals': len(integrals),
        'branch_arguments': tuple(sorted(branches, key=sp.default_sort_key)),
    }


def polynomial_terms(expression, generators):
    """Collect only the named carriers, retaining factored coefficients.

    SymPy's general coefficient-domain discovery otherwise expands the whole
    inherited rational coefficient algebra during this structural operation.
    """
    zero = (0,)*len(generators)
    indices = {g: i for i, g in enumerate(generators)}

    def multiply(a, b):
        result = {}
        for pa, ca in a.items():
            for pb, cb in b.items():
                power = tuple(x+y for x, y in zip(pa, pb))
                result[power] = result.get(power, sp.S.Zero) + ca*cb
        return result

    @lru_cache(maxsize=None)
    def visit(node):
        if node in indices:
            power = list(zero)
            power[indices[node]] = 1
            return {tuple(power): sp.S.One}
        if not node.has(*generators):
            return {zero: node}
        if node.is_Add:
            result = {}
            for arg in node.args:
                for power, coefficient in visit(arg).items():
                    result[power] = result.get(power, sp.S.Zero) + coefficient
            return result
        if node.is_Mul:
            result = {zero: sp.S.One}
            for arg in node.args:
                result = multiply(result, visit(arg))
            return result
        if node.is_Pow and node.exp.is_Integer and node.exp >= 0:
            result = {zero: sp.S.One}
            for _ in range(int(node.exp)):
                result = multiply(result, visit(node.base))
            return result
        raise NotImplementedError(('nonpolynomial profile carrier', node.func))

    return tuple(visit(expression).items())


class EdgeReduction:
    """Partial Fourier transform in an orthonormal chart with n=e_3.

    The profile and field ansatz is the sole manually assembled physical
    expression. All derivative factors and momentum constraints are computed
    from this ansatz and the phases on the actual imported integral nodes.
    FOURIER_PROFILE_BINDINGS is never read to construct a hat replacement.
    """

    def __init__(self, rows):
        self.rows = rows
        self.values = {k: r['value'] for k, r in rows.items()}
        symbols = set().union(*(self.values[k].atoms(sp.Symbol) for k in CLOSED_KEYS))
        self.symbols = {s.name: s for s in symbols}
        self.x = tuple(self.symbols[f's11cc2X{i}'] for i in range(1, 4))
        self.y = tuple(self.symbols[f's11cc2Y{i}'] for i in range(1, 4))
        self.t = self.symbols['s11cc2Time']
        self.ell = self.values['L_W']
        self.omega = self.values['omega']
        self.tangents = sp.symbols('s11cdTangentialMomentum1 s11cdTangentialMomentum2', real=True)
        self.z, self.zp = sp.symbols('s11cdNormalPosition s11cdSourceNormalPosition', real=True)
        self.xi = sp.Symbol('s11cdProfileCoordinate', real=True)
        self.profiles = {p: sp.Function('s11cd' + p.upper() + 'Profile') for p in ('w', 'm')}
        self.functions = {}
        self.normal_transfers = set()
        self.momentum_groups = tuple(
            tuple(self.symbols[pattern.format(i)] for i in range(1, 4))
            for pattern in ('s11cc1_k_output_{}', 's11cc1_k_input_{}', 's11cc2MiddleMomentum{}')
        )
        self.tangent_momenta = tuple(k for group in self.momentum_groups for k in group[:2])
        self.normal_map = {self.x[2]: self.z, self.y[2]: self.zp}
        self.normal_map.update({group[2]: sp.Symbol('s11cd' + name + 'NormalMomentum', real=True)
                                for group, name in zip(self.momentum_groups, ('Output', 'Input', 'Middle'))})
        self.branch_equations = tuple(sorted(set(
            eq for key in CLOSED_KEYS for _, payload in self.values[key]
            for eq in named(payload, 'COMPUTED_BRANCH_BINDINGS')
        ), key=sp.default_sort_key))
        self.branch_map = {}
        for eq in self.branch_equations:
            if eq.lhs in self.branch_map and self.branch_map[eq.lhs] != eq.rhs:
                raise ValueError(('case-dependent branch binding', eq.lhs))
            self.branch_map[eq.lhs] = eq.rhs
        # Distributional normalization computed from the mass of a regulated
        # Fourier kernel. No numerical Fourier weight is supplied here.
        u, q = sp.symbols('s11cdFourierPosition s11cdFourierMomentum', real=True)
        a = sp.Symbol('s11cdFourierRegulator', positive=True)
        self.regulated_ansatz = sp.exp(-a*u**2 + sp.I*q*u)
        self.regulated_transform = sp.integrate(self.regulated_ansatz, (u, -sp.oo, sp.oo))
        self.fourier_mass = sp.integrate(self.regulated_transform, (q, -sp.oo, sp.oo))
        self.normalization_operands = (self.regulated_ansatz, self.regulated_transform, self.fourier_mass)
        self.regulator = sp.Symbol('s11cdAbelRegulator', positive=True)
        self.transfer = sp.Symbol('s11cdDimensionlessTransfer', real=True)
        self.end_values = {(p, end): sp.Limit(f(self.xi), self.xi, end)
                           for p, f in self.profiles.items() for end in (-sp.oo, sp.oo)}
        phase = sp.exp(-sp.I*self.transfer*self.xi)
        positive = sp.integrate(phase*sp.exp(-self.regulator*self.xi),
                                (self.xi, 0, sp.oo), conds='none')
        negative = sp.integrate(phase*sp.exp(self.regulator*self.xi),
                                (self.xi, -sp.oo, 0), conds='none')
        self.abel_halves = (negative, positive)
        self.abel_constant = sp.factor(negative+positive)
        self.abel_even = sp.factor((positive+negative)/2)
        self.abel_odd = sp.factor((positive-negative)/2)
        test = sp.Function('s11cdSchwartzTest')(self.transfer)
        self.distribution_record = {
            'REGULATED_HALF_LINE_OPERANDS': (phase*sp.exp(self.regulator*self.xi),
                                              phase*sp.exp(-self.regulator*self.xi)),
            'REGULATED_HALF_LINE_TRANSFORMS': self.abel_halves,
            'CONSTANT_TRANSFORM': self.abel_constant,
            'HEAVISIDE_EVEN_ODD': (self.abel_even, self.abel_odd),
            'DELTA_MASS': sp.integrate(self.abel_even, (self.transfer, -sp.oo, sp.oo)),
            'PRINCIPAL_VALUE_WEAK_ACTION': sp.Limit(sp.Integral(test*self.abel_odd,
                (self.transfer, -sp.oo, sp.oo)), self.regulator, 0, dir='+'),
        }

    def subtraction(self, profile):
        f = self.profiles[profile](self.xi)
        left, right = (self.end_values[(profile, e)] for e in (-sp.oo, sp.oo))
        jump = right-left
        localized = f-left-jump*sp.Heaviside(self.xi)
        return f, left, right, jump, localized

    def prescribe(self, normal):
        """Abel boundary value, with the weak limit taken AFTER convolution.

        The fixed origin is xi=0. This transformed representation is explicitly
        conditioned on the additional L1 half-line-tail premise in §1c. Jets
        keep their ordinary reduced integrals. No pointwise evaluation of the
        regulator limit is used to erase its delta or principal-value part.
        """
        if not isinstance(normal, sp.Integral):
            return normal
        profile = self.profiles['w'](self.xi)
        if not normal.function.has(profile) or normal.function.has(sp.Derivative):
            return normal
        phases = list(normal.function.atoms(sp.exp))
        if len(phases) != 1:
            raise NotImplementedError(('normal transform phase', normal))
        phase = phases[0]
        scale = sp.cancel(normal.function/(profile*phase))
        if scale.has(self.xi):
            raise NotImplementedError(('normal transform multiplier', scale))
        s = sp.expand(-sp.diff(phase.args[0], self.xi)/sp.I)
        _, left, _, jump, localized = self.subtraction('w')
        return scale*(left*self.abel_constant.subs(self.transfer, s)
                      +jump*self.abel_halves[1].subs(self.transfer, s)
                      +sp.Integral(phase*localized, *normal.limits))

    def weak_limit(self, value):
        return map_leaves(value, lambda v: sp.Limit(v, self.regulator, 0, dir='+')
                          if v.has(self.regulator) else v)

    def wave_phase(self, point):
        return sp.exp(sp.I * (sum(k*x for k, x in zip(self.tangents, point[:2])) - self.omega*self.t))

    def field_phase(self, node):
        phase = self.wave_phase(node.args)
        return sp.Pow(phase, -1) if node.func.__name__.startswith('s11cc2Test') else phase

    @lru_cache(maxsize=None)
    def applied_ansatz(self, node):
        name = node.func.__name__
        if name in ('s11cc2Coefficientw1Profile', 's11cc2Coefficientm1Profile'):
            return self.profiles['w' if 'w1' in name else 'm'](node.args[2]/self.ell)
        if name.startswith(('s11cc2Field', 's11cc2Trial', 's11cc2Test')):
            fn = self.functions.setdefault(node.func, sp.Function(name.replace('s11cc2', 's11cdReduced', 1)))
            return fn(node.args[2])*self.field_phase(node)
        return node

    @lru_cache(maxsize=None)
    def local_atom(self, node):
        if isinstance(node, sp.Symbol):
            match = re.fullmatch(r'([wm])1_profile(?:_((?:d[123])+))?', node.name)
            if match:
                profile, suffix = match.groups()
                value = self.profiles[profile](self.x[2]/self.ell)
                indices = re.findall(r'd([123])', suffix or '')
                differentiated = sp.diff(value, *(self.x[int(i)-1] for i in indices)) if indices else value
                return (self.ell**len(indices)*differentiated).xreplace(self.normal_map)
        if isinstance(node, AppliedUndef):
            result = self.applied_ansatz(node)
            if result != node:
                if node.func in self.functions:
                    result = sp.cancel(result / self.field_phase(node))
                return result.xreplace(self.normal_map)
        if isinstance(node, sp.Derivative):
            funcs = node.expr.atoms(AppliedUndef)
            if len(funcs) == 1:
                f = next(iter(funcs))
                ansatz = self.applied_ansatz(f)
                if ansatz != f:
                    result = sp.diff(node.expr.xreplace({f: ansatz}), *node.variable_count)
                    if f.func in self.functions:
                        result = sp.cancel(result / self.field_phase(f))
                    return result.xreplace(self.normal_map)
        return node

    @lru_cache(maxsize=None)
    def hat(self, node):
        """Fourier coefficient of the profile multiplication ansatz.

        The normalization follows Fourier inversion of the local product in
        the engine's angular phase convention. The imported hat definitions
        remain independent source operands of the subsequent record.
        """
        profile = self.profiles['w'](self.y[2]/self.ell)
        match = re.search(r'JetHat([123])$', node.func.__name__)
        if match:
            profile = self.ell * sp.diff(profile, self.y[int(match[1])-1])
        phase = sp.exp(-sp.I * sum(q*y for q, y in zip(node.args, self.y)))
        phase_normal = phase.xreplace({self.y[0]: 0, self.y[1]: 0})
        if profile == 0:
            normal = profile
        else:
            normal = sp.Integral(phase_normal*profile, (self.y[2], -sp.oo, sp.oo))
        normal = normal.transform(self.y[2], (self.ell*self.xi, self.xi)) if isinstance(normal, sp.Integral) else normal
        self.normal_transfers.add(sp.expand((self.ell*node.args[2]).xreplace(self.normal_map)))
        reduced = (self.prescribe(normal)/self.fourier_mass).xreplace(self.normal_map)
        return reduced, node.args[:2]

    def strip_local(self, value):
        # Replace outer derivative nodes before their applied-function children.
        atoms = value.atoms(sp.Derivative, AppliedUndef, sp.Symbol)
        replacements = {a: self.local_atom(a) for a in atoms}
        return memo_xreplace(value, replacements)

    def branches(self, value):
        return memo_xreplace(value, self.branch_map)

    def profile_definition(self, integral):
        """Reduce an actual source definition as a separate 3-D operand.

        This method is never called by hat() or by the action contraction.
        The return value is the coefficient of the tangential delta factors.
        """
        phases = [e for e in integral.function.atoms(sp.exp)
                  if any(e.has(y) for y in self.y[:2])]
        if len(phases) != 1:
            raise NotImplementedError(('profile phase multiplicity', len(phases)))
        phase = phases[0]
        constraints = tuple(sp.expand(sp.diff(phase.args[0], y)/sp.I) for y in self.y[:2])
        coefficient = self.strip_local(integral.function.xreplace({phase: 1}))
        normal_phase = phase.xreplace({self.y[0]: 0, self.y[1]: 0}).xreplace(self.normal_map)
        normal = sp.Integral(coefficient*normal_phase*self.fourier_mass**len(self.y[:2]),
                             (self.zp, -sp.oo, sp.oo))
        normal = normal.transform(self.zp, (self.ell*self.xi, self.xi))
        if isinstance(normal, sp.Integral) and normal.function == 0:
            normal = normal.doit()
        return self.prescribe(normal), constraints

    @lru_cache(maxsize=None)
    def integral(self, original):
        """Contract every tangential measure using its phase constraints.

        This routine handles action integrals. Profile-definition integrals
        require a different output distribution space and are not passed here.
        """
        variables = tuple(l[0] for l in original.limits)
        if not all(y in variables for y in self.y):
            raise NotImplementedError(('integral domain', original.limits))
        phases = [e for e in original.function.atoms(sp.exp)
                  if any(e.has(y) for y in self.y[:2])]
        if len(phases) != 1:
            raise NotImplementedError(('phase multiplicity', len(phases)))
        phase = phases[0]
        phase_argument = phase.args[0]
        source_phase_argument = phase_argument + self.wave_phase(self.y).args[0]
        source_constraints = tuple(sp.expand(sp.diff(source_phase_argument, y)/sp.I)
                                   for y in self.y[:2])
        integrand = self.branches(self.strip_local(original.function.xreplace({phase: 1})))
        hats = sorted((f for f in integrand.atoms(AppliedUndef)
                       if f.func.__name__.startswith('s11cc2Fourier')), key=sp.default_sort_key)
        hat_data = {h: self.hat(h) for h in hats}
        integrand = integrand.xreplace({h: r for h, (r, _) in hat_data.items() if r == 0})
        hats = [h for h in hats if integrand.has(h)]
        terms = polynomial_terms(integrand, hats) if hats else [((), integrand)]
        integrated_momenta = tuple(k for k in self.tangent_momenta if k in variables)
        normal_limits = tuple((self.normal_map.get(v, v), lo, hi) for v, lo, hi in original.limits
                              if v not in (*self.y[:2], *integrated_momenta))
        terms_out = []
        constraint_records = []
        for powers, coefficient in terms:
            constraints = list(source_constraints)
            hat_product = sp.S.One
            for h, power in zip(hats, powers):
                reduced, tangent_args = hat_data[h]
                for _ in range(power):
                    constraints.extend(tangent_args)
                    hat_product *= reduced
            matrix, rhs = sp.linear_eq_to_matrix(constraints, integrated_momenta)
            if matrix.rows != matrix.cols:
                raise NotImplementedError(('constraint rank', matrix.shape, original.limits))
            determinant = matrix.det()
            solution = sp.linsolve((matrix, rhs), integrated_momenta)
            solutions = list(solution)
            if len(solutions) != 1:
                raise NotImplementedError(('tangential constraint solutions', solution))
            substitution = dict(zip(integrated_momenta, solutions[0]))
            exponent = phase_argument.xreplace(substitution)
            # The common outgoing field character is removed after integrating
            # the conserved tangential measures, before any flux operation.
            exponent += self.wave_phase(self.y).args[0] - self.wave_phase(self.x).args[0]
            exponent = sp.expand(exponent).xreplace(self.normal_map)
            coefficient = coefficient.xreplace(substitution).xreplace(self.normal_map)
            result = coefficient*hat_product*sp.exp(exponent)*self.fourier_mass**len(self.y[:2])/sp.Abs(determinant)
            terms_out.append(sp.Integral(result, *normal_limits) if result != 0 else result)
            constraint_records.append((cas(constraints), matrix, rhs, solution, determinant))
        return sp.Add(*terms_out), cas(constraint_records)

    def value(self, value):
        replacements = {}
        records = []
        for integral in sorted(value.atoms(sp.Integral), key=sp.default_sort_key):
            reduced, constraints = self.integral(integral)
            replacements[integral] = reduced
            records.append((integral, reduced, constraints))
        result = memo_xreplace(self.branches(self.strip_local(memo_xreplace(value, replacements))), self.normal_map)
        return result, records

    def payload(self, source, reduced):
        """Recompute the full five-slot payload in the reduced convention."""
        metadata = PHYSICAL_METADATA.record(reduced)
        grade_support = set()
        for _, row in metadata:
            grade_support.update(tuple(g) for g in named(row, 'MULTIGRADE'))
        def restored_tree(value, source_dimensions):
            if isinstance(source_dimensions, sp.MatrixBase):
                inferred = PHYSICAL_METADATA.dimensions.measure(value)
                return sp.ImmutableMatrix(inferred if inferred is not None else source_dimensions)
            if isinstance(value, Str):
                return value
            return sp.Tuple(*(restored_tree(a,b) for a,b in zip(value,source_dimensions)))
        dimension_tree = restored_tree(reduced,named(source,'DIMENSION_L_T_M'))
        tangent_map = {k: q for group in self.momentum_groups for k, q in zip(group[:2], self.tangents)}
        branch_records = []
        for eq in named(source, 'COMPUTED_BRANCH_BINDINGS'):
            lhs = sp.Function('s11cdOutgoingBulkNormalMomentum')(
                *(a.xreplace(tangent_map).xreplace(self.normal_map) for a in eq.lhs.args))
            rhs = eq.rhs.xreplace(tangent_map).xreplace(self.normal_map)
            PHYSICAL_METADATA.dimensions.known[lhs.func] = PHYSICAL_METADATA.dimensions.measure(rhs)
            branch_records.append(sp.Eq(lhs, rhs, evaluate=False))
        profile_records = []
        for eq in named(source, 'FOURIER_PROFILE_BINDINGS'):
            h = eq.lhs
            lhs = sp.Function(h.func.__name__.replace('s11cc2Fourier', 's11cdFourier'))(
                h.args[2].xreplace(self.normal_map))
            PHYSICAL_METADATA.dimensions.known[lhs.func] = tuple(a+2*b for a,b in zip(
                PHYSICAL_METADATA.dimensions.measure(h),PHYSICAL_METADATA.dimensions.measure(self.tangents[0])))
            # The RHS comes only from hat(), never from the source definition.
            profile_records.append(sp.Eq(lhs, self.hat(h)[0], evaluate=False))
        return cas({'VALUE': reduced, 'MULTIGRADE': sorted(grade_support),
                    'DIMENSION_L_T_M': dimension_tree,
                    'COMPUTED_BRANCH_BINDINGS': branch_records,
                    'FOURIER_PROFILE_BINDINGS': profile_records})


class FourierCarrierReconstruction:
    """Inverse transforms on the supplied smooth interface class.

    Source operands come from the imported three-dimensional definitions;
    image operands come from hat(). Neither route supplies the other's answer.
    A computed Gaussian approximate identity supplies weak delta kernels.
    Abel rational terms are inverted by residues on the two coordinate
    half-lines before removing the regulator. The symmetric value at the
    subtraction origin is retained on both terms of the decomposition.
    """

    def __init__(self, reduction):
        self.r = reduction
        self.q = sp.symbols('s11cdInverseMomentum1:4', real=True)
        self.u = sp.symbols('s11cdInverseSource1:4', real=True)
        self.x = sp.symbols('s11cdInversePosition1:4', real=True)
        self.s = sp.symbols('s11cdInverseTransfer1:4', real=True)
        self.a = sp.Symbol('s11cdInverseGaussianRegulator', positive=True)
        self.p = sp.Symbol('s11cdInversePositiveCoordinate', positive=True)
        self.nu = sp.Symbol('s11cdInverseTestFrequency', real=True)
        self.c = sp.Symbol('s11cdInverseDeltaWeight', real=True)
        self.complex_transfer = sp.Symbol('s11cdInverseComplexTransfer')
        units = PHYSICAL_METADATA.dimensions
        for v in (*self.u, *self.x, *self.s, self.a, self.p, self.nu, self.c, self.complex_transfer):
            units.known[v] = units.zero
        for v in self.q:
            units.known[v] = units.measure(reduction.tangents[0])
        units.measure.cache_clear()
        s, x = self.s[2], self.x[2]
        self.gaussian_ansatz = sp.exp(-self.a*s**2+sp.I*s*x)
        self.gaussian = sp.integrate(self.gaussian_ansatz, (s, -sp.oo, sp.oo))
        self.mass = sp.integrate(self.gaussian, (x, -sp.oo, sp.oo))
        self.delta_ansatz = self.c*sp.DiracDelta(x)
        delta_mass = sp.integrate(self.delta_ansatz, (x, -sp.oo, sp.oo))
        weight = sp.solve(delta_mass-self.mass, self.c)
        if len(weight) != 1:
            raise NotImplementedError(('inverse delta weight', weight))
        self.delta = self.delta_ansatz.subs(self.c, weight[0])
        character = sp.exp(-sp.I*self.nu*x)
        gaussian_character = sp.integrate(self.gaussian*character, (x, -sp.oo, sp.oo))
        delta_character = sp.integrate(self.delta*character, (x, -sp.oo, sp.oo))
        weak_character = sp.limit(gaussian_character, self.a, 0, dir='+')
        moment = sp.integrate(x**2*self.gaussian, (x, -sp.oo, sp.oo))/self.mass
        # The contour orientation is computed from a unit-circle ansatz.
        theta = sp.Symbol('s11cdInverseContourAngle', real=True)
        units.known[theta] = units.zero
        circle = sp.exp(sp.I*theta)
        self.cauchy_weight = sp.integrate(sp.diff(circle, theta)/circle,
                                         (theta, 0, self.mass))
        self.kernel_record = (self.gaussian_ansatz, self.gaussian, self.mass,
                              self.delta_ansatz, delta_mass, cas(weight), self.delta,
                              character, gaussian_character, weak_character, delta_character,
                              moment, sp.limit(moment, self.a, 0, dir='+'),
                              circle, self.cauchy_weight)
        self.kernel_residual = weak_character-delta_character

    def emit_kernel(self):
        self.physical_zeros('INVERSE_FOURIER_WEAK_KERNEL_OPERANDS', self.kernel_record)
        physical('INVERSE_FOURIER_WEAK_KERNEL_CHARACTER_RESIDUAL', self.kernel_residual,
                 zero_dimensions={(): PHYSICAL_METADATA.dimensions.zero})

    def physical_zeros(self, name, value, overrides=None):
        # These inverse-coordinate coefficients are dimensionless. Images and
        # argument-chart residuals have separate units supplied by their maps.
        units = PHYSICAL_METADATA.dimensions
        for substitution in cas(value).atoms(sp.Subs):
            if any(point == 0 for point in substitution.point):
                # Evaluation at the numerical origin preserves the units of
                # the dimensionless inverse-coordinate derivative operand.
                chart_operand = sp.Subs(substitution.expr, substitution.variables,
                    tuple(self.x[2] if point == 0 else point for point in substitution.point))
                units.known[substitution] = units.measure(chart_operand)
        zeros = {path: PHYSICAL_METADATA.dimensions.zero
                 for path, expression in leaves(cas(value)) if expression == 0}
        zeros.update(overrides or {})
        physical(name, value, zero_dimensions=zeros)

    @lru_cache(maxsize=None)
    def template(self, function, definitions):
        equations = [eq for eq in definitions if eq.lhs.func == function]
        if len(equations) != 1:
            raise NotImplementedError(('source Fourier definition count', function, len(equations)))
        equation = equations[0]
        argument_equations = tuple(a-q for a, q in zip(equation.lhs.args, self.q))
        solved = sp.solve(argument_equations, self.r.momentum_groups[0], dict=True)
        if len(solved) != 1:
            raise NotImplementedError(('source Fourier argument chart', solved))
        substitution = solved[0]
        source = equation.rhs.subs(substitution, simultaneous=True)
        argument_residual = tuple(sp.expand(a.subs(substitution)-q)
                                  for a, q in zip(equation.lhs.args, self.q))
        return source, cas(substitution), cas(argument_residual)

    @lru_cache(maxsize=None)
    def source_definition(self, node, definitions):
        source, _, _ = self.template(node.func, definitions)
        return source.xreplace(dict(zip(self.q, node.args)))

    @lru_cache(maxsize=None)
    def source_image(self, node, definitions):
        image, constraints = self.r.profile_definition(self.source_definition(node, definitions))
        if isinstance(image, sp.Integral):
            # Integral linearity separates the source Fourier coefficient
            # while keeping the computed normal measure Jacobian inside.
            jacobian = sp.diff(self.r.ell*self.r.xi, self.r.xi)
            coefficient, dependent = image.function.as_independent(self.r.xi, as_Add=False)
            image = (coefficient/jacobian)*sp.Integral(jacobian*dependent, *image.limits)
        return image, constraints

    @lru_cache(maxsize=None)
    def source_inverse(self, function, definitions):
        r = self.r
        source, substitution, argument_residual = self.template(function, definitions)
        phases = list(source.function.atoms(sp.exp))
        if len(phases) != 1:
            raise NotImplementedError(('source inverse phase multiplicity', len(phases)))
        phase = phases[0]
        profile_map = {f.func: (lambda *args, fn=f.func: r.applied_ansatz(fn(*args)))
                       for f in source.function.atoms(AppliedUndef)}
        coefficient = dag_substitute(source.function.xreplace({phase: 1}), profile_map)
        position_map = dict(zip(r.y, (r.ell*u for u in self.u)))
        momentum_map = dict(zip(self.q, (s/r.ell for s in self.s)))
        position_jacobian = sp.Matrix([position_map[y] for y in r.y]).jacobian(self.u).det()
        momentum_jacobian = sp.Matrix([momentum_map[q] for q in self.q]).jacobian(self.s).det()
        inverse_phase = sp.I*sum(s*x for s, x in zip(self.s, self.x))
        exponent = sp.expand(phase.args[0].xreplace(position_map).xreplace(momentum_map)+inverse_phase)
        constraints = tuple(sp.expand(sp.diff(exponent, s)/sp.I) for s in self.s)
        remainder = sp.simplify(exponent-sp.I*sum(s*c for s, c in zip(self.s, constraints)))
        if remainder != 0:
            raise NotImplementedError(('source inverse phase remainder', remainder))
        kernels = tuple(self.delta.subs(self.x[2], constraint) for constraint in constraints)
        body = coefficient.subs(position_map, simultaneous=True).doit()*position_jacobian*momentum_jacobian
        if body.has(*self.s, *self.q):
            raise NotImplementedError(('momentum-dependent source inverse coefficient', body))
        inverse_integrand = body*sp.prod(kernels)
        value = sp.integrate(inverse_integrand, *((u, -sp.oo, sp.oo) for u in self.u))
        return value, cas((source, substitution, argument_residual, coefficient,
                           position_jacobian, momentum_jacobian, exponent, constraints,
                           inverse_integrand, value))

    @lru_cache(maxsize=None)
    def invert_image(self, image, normal_argument):
        r, s, x = self.r, self.s[2], self.x[2]
        variables = sorted(normal_argument.free_symbols, key=sp.default_sort_key)
        if not variables:
            raise NotImplementedError(('inverse normal argument chart', normal_argument))
        variable = variables[0]
        roots = sp.solve(normal_argument-s/r.ell, variable)
        if len(roots) != 1:
            raise NotImplementedError(('inverse normal argument roots', roots))
        substitution = {variable: roots[0]}
        momentum_jacobian = sp.diff(s/r.ell, s)
        normalized = sp.expand(image.subs(substitution, simultaneous=True)*momentum_jacobian,
                               power_exp=False)
        old_momenta = set(r.normal_map[g[2]] for g in r.momentum_groups)
        if dag_free_symbols(normalized) & old_momenta:
            raise NotImplementedError(('uncontracted inverse normal transfer', normal_argument))
        # Nested Limit objects are independent coefficient carriers here.
        # SymPy otherwise leaves elementary integrals and limits unevaluated.
        constants = {value: sp.Symbol('s11cdInverseEndCoefficient'+str(i))
                     for i, value in enumerate(sorted(normalized.atoms(sp.Limit), key=sp.default_sort_key))}
        for value, symbol in constants.items():
            PHYSICAL_METADATA.dimensions.known[symbol] = PHYSICAL_METADATA.dimensions.measure(value)
        restored_constants = {v: k for k, v in constants.items()}
        coefficient_image = normalized.xreplace(constants)
        integrals = sorted(coefficient_image.atoms(sp.Integral), key=sp.default_sort_key)
        regular, rational, witnesses = sp.S.Zero, sp.S.Zero, []
        for powers, coefficient in polynomial_terms(coefficient_image, integrals):
            if not any(powers):
                rational += coefficient
                continue
            if sum(powers) != 1:
                raise NotImplementedError(('inverse profile integral degree', powers))
            integral = integrals[powers.index(1)]
            phases = list(integral.function.atoms(sp.exp))
            if len(phases) != 1 or coefficient.has(s):
                raise NotImplementedError(('inverse profile convolution form', integral))
            phase = phases[0]
            exponent = sp.expand(phase.args[0]+sp.I*s*x)
            constraint = sp.expand(sp.diff(exponent, s)/sp.I)
            if sp.simplify(exponent-sp.I*s*constraint) != 0:
                raise NotImplementedError(('inverse profile phase remainder', exponent))
            kernel = self.delta.subs(x, constraint)
            body = coefficient*integral.function.xreplace({phase: 1})*kernel
            convolution = sp.integrate(sp.expand(body), *integral.limits)
            regular += convolution
            witnesses.append((integral, coefficient, exponent, kernel, body, convolution))
        raw_rational = rational
        rational = sp.cancel(rational)
        analytic = rational.xreplace({s: self.complex_transfer})
        poles = sp.solve(sp.denom(analytic), self.complex_transfer) if analytic != 0 else []
        upper, lower, residue_rows = sp.S.Zero, sp.S.Zero, []
        for pole in poles:
            imaginary = sp.simplify(sp.im(pole))
            residue = sp.residue(sp.exp(sp.I*self.complex_transfer*x)*analytic, self.complex_transfer, pole)
            residue_rows.append((pole, imaginary, residue))
            if imaginary.is_positive:
                upper += self.cauchy_weight*residue
            elif imaginary.is_negative:
                lower -= self.cauchy_weight*residue
            else:
                raise NotImplementedError(('inverse Abel contour pole', pole))
        if rational != 0 and sp.limit(rational, s, sp.oo) != 0:
            raise NotImplementedError(('inverse Abel polynomial part', rational))
        origin = (upper.subs(x, 0)+lower.subs(x, 0))/2
        coordinates = (self.p, -self.p, sp.S.Zero)
        regulated = tuple(sp.simplify(regular.subs(x, position)+part.subs(x, position))
                          for position, part in zip(coordinates, (upper, lower, origin)))
        weak = tuple(sp.limit(v, r.regulator, 0, dir='+') if v.has(r.regulator) else v
                     for v in regulated)
        tangent_integrals = tuple(sp.Integral(sp.exp(sp.I*s*x)*sp.DiracDelta(s/r.ell)/r.ell,
                                               (s, -sp.oo, sp.oo))
                                 for s, x in zip(self.s[:2], self.x[:2]))
        tangent_inverse = sp.prod(v.doit() for v in tangent_integrals)
        weak = tuple(sp.simplify(tangent_inverse*v).xreplace(restored_constants) for v in weak)
        operands = cas((image, normal_argument, substitution, momentum_jacobian, normalized,
                        witnesses, rational, residue_rows, coordinates, regulated, weak,
                        tangent_integrals, tangent_inverse, raw_rational, restored_constants,
                        analytic, poles))
        return cas(weak), operands

    def carrier(self, node, definitions, suffix):
        r = self.r
        source = self.source_definition(node, definitions)
        image, tangent_arguments = r.hat(node)
        source_image, constraints = self.source_image(node, definitions)
        raw_source_image = r.profile_definition(source)[0]
        source_point, source_operands = self.source_inverse(node.func, definitions)
        reconstructed, inverse_operands = self.invert_image(
            image, node.args[2].xreplace(r.normal_map))
        coordinates = (self.p, -self.p, sp.S.Zero)
        source_values = cas(tuple(source_point.subs(self.x[2], p) for p in coordinates))
        residual = cas(tuple(sp.simplify(a-b) for a, b in zip(reconstructed, source_values)))
        unit = PHYSICAL_METADATA.dimensions.measure(node)
        inverse_unit = tuple(a+3*b for a, b in zip(unit, PHYSICAL_METADATA.dimensions.measure(self.q[0])))
        self.physical_zeros('CARRIER_INVERSE_SOURCE_OPERANDS_'+suffix, source_operands,
                            {(2, i): PHYSICAL_METADATA.dimensions.measure(self.q[i]) for i in range(3)})
        self.physical_zeros('CARRIER_INVERSE_REDUCED_OPERANDS_'+suffix, inverse_operands,
                            {(0,): tuple(a-b for a, b in zip(inverse_unit,
                                PHYSICAL_METADATA.dimensions.measure(self.q[0])))})
        physical('CARRIER_INVERSE_VALUES_'+suffix, (source_values, reconstructed),
                 zero_dimensions={(i, j): inverse_unit for i in range(2) for j in range(3)})
        physical('CARRIER_INVERSE_RESIDUAL_'+suffix, residual,
                 zero_dimensions={(i,): inverse_unit for i in range(3)})
        self.physical_zeros('CARRIER_INVERSE_REMAINDER_CENSUS_'+suffix,
            (len(reconstructed.atoms(sp.Integral)), int(reconstructed.has(r.regulator)),
             len(dag_free_symbols(reconstructed) & set(r.normal_map[g[2]] for g in r.momentum_groups))))
        image_unit = tuple(a-b for a, b in zip(inverse_unit, PHYSICAL_METADATA.dimensions.measure(self.q[0])))
        self.physical_zeros('CARRIER_INVERSE_ARGUMENT_OPERANDS_'+suffix,
                 (node, source, image, tangent_arguments, source_image, constraints),
                 {(2,): image_unit, (4,): image_unit})
        self.physical_zeros('CARRIER_SOURCE_LINEARITY_OPERANDS_'+suffix, (raw_source_image, source_image),
                            {(0,): image_unit, (1,): image_unit})
        delta = sp.prod(sp.DiracDelta(c) for c in constraints)
        physical('CARRIER_SOURCE_IMAGE_RESIDUAL_'+suffix, (image-source_image)*delta,
                 zero_dimensions={(): unit})


class EdgeReconstruction:
    """Reconstruct on the supplied tangentially homogeneous ansatz class.

    The source route eliminates delta constraints one variable at a time.
    It does not consume EdgeReduction's constraint matrix, roots or Jacobian.
    The inverse lift restores the outgoing character and original normal
    coordinate. This is a class-restricted reconstruction, not an inverse on
    arbitrary three-dimensional backgrounds.
    """

    def __init__(self, reduction, carriers):
        self.r = reduction
        self.carriers = carriers

    def lift(self, value):
        r = self.r
        inverse = {v: k for k, v in r.normal_map.items()}
        return map_leaves(value, lambda e: e.xreplace(inverse)*r.wave_phase(r.x))

    @lru_cache(maxsize=None)
    def source_integral(self, original, definitions):
        r = self.r
        variables = tuple(l[0] for l in original.limits)
        phases = [e for e in original.function.atoms(sp.exp)
                  if any(e.has(y) for y in r.y[:2])]
        if len(phases) != 1:
            raise NotImplementedError(('reconstruction character', len(phases)))
        phase = phases[0]
        character = phase.args[0]+r.wave_phase(r.y).args[0]
        constraints_y = [sp.expand(sp.diff(character, y)/sp.I) for y in r.y[:2]]
        source = r.branches(r.strip_local(original.function.xreplace({phase: 1})))
        hats = sorted((f for f in source.atoms(AppliedUndef)
                       if f.func.__name__.startswith('s11cc2Fourier')), key=sp.default_sort_key)
        data = {h: self.carriers.source_image(h, definitions) for h in hats}
        source = source.xreplace({h: image for h, (image, _) in data.items() if image == 0})
        hats = [h for h in hats if source.has(h)]
        momenta = tuple(k for k in r.tangent_momenta if k in variables)
        limits = tuple((r.normal_map.get(v, v), a, b) for v, a, b in original.limits
                       if v not in (*r.y[:2], *momenta))
        terms, witnesses = [], []
        for powers, coefficient in polynomial_terms(source, hats):
            constraints = list(constraints_y)
            multiplier = coefficient
            for h, power in zip(hats, powers):
                image, tangent = data[h]
                multiplier *= image**power
                constraints.extend(q for _ in range(power) for q in tangent)
            original_constraints = tuple(constraints)
            substitution, steps = {}, []
            for variable in reversed(momenta):
                candidates = [(i, q) for i, q in enumerate(constraints) if q.has(variable)]
                if not candidates:
                    raise NotImplementedError(('reconstruction free tangent', variable))
                i, q = candidates[0]
                roots = sp.solve(q, variable)
                if len(roots) != 1:
                    raise NotImplementedError(('reconstruction roots', variable, roots))
                root = roots[0]
                jacobian = sp.diff(q, variable).subs(variable, root)
                multiplier = memo_xreplace(multiplier, {variable: root})/sp.Abs(jacobian)
                constraints.pop(i)
                constraints = [sp.expand(c.subs(variable, root)) for c in constraints]
                substitution = {k: v.subs(variable, root) for k, v in substitution.items()}
                substitution[variable] = root
                steps.append((variable, q, root, jacobian))
            if constraints:
                raise NotImplementedError(('reconstruction remaining deltas', constraints))
            exponent = sp.expand((character-r.wave_phase(r.x).args[0]).subs(substitution))
            body = (multiplier*sp.exp(exponent)*r.fourier_mass**len(r.y[:2])).xreplace(r.normal_map)
            terms.append(sp.Integral(body, *limits) if body != 0 else body)
            witnesses.append((original_constraints, steps))
        return sp.Add(*terms), cas(witnesses)

    def row(self, source, reduced, records, suffix, source_units, definitions):
        r = self.r
        replacements = {}
        for index, (original, image, _) in enumerate(records):
            independent, witness = self.source_integral(original, definitions)
            replacements[original] = independent
            restored, source_on_class = self.lift(image), self.lift(independent)
            residual = tree_difference(restored, source_on_class)
            unit = PHYSICAL_METADATA.dimensions.measure(original)
            name = suffix+'_'+str(index)
            fingerprinted('RECONSTRUCTION_INTEGRAL_OPERANDS_'+name,
                          cas((original, image, restored, source_on_class)),
                          {(i,): unit for i in range(4)})
            fingerprinted('RECONSTRUCTION_INTEGRAL_RESIDUAL_'+name, residual,
                          {(): unit})
            physical('RECONSTRUCTION_DELTA_ELIMINATION_'+name, witness)
        independent = memo_xreplace(r.branches(r.strip_local(memo_xreplace(source, replacements))), r.normal_map)
        restored, source_on_class = self.lift(reduced), self.lift(independent)
        fingerprinted('RECONSTRUCTION_ROW_SOURCE_'+suffix, source, source_units)
        fingerprinted('RECONSTRUCTION_ROW_LIFT_'+suffix, restored, source_units)
        fingerprinted('RECONSTRUCTION_ROW_SOURCE_ON_CLASS_'+suffix, source_on_class, source_units)
        fingerprinted('RECONSTRUCTION_ROW_RESIDUAL_'+suffix,
                      tree_difference(restored, source_on_class), source_units)


class ConstantEndPencil:
    """Constant-end weak limits and normal Fourier symbols of reduced rows.

    The end map acts on *all* profile occurrences in the reduced kernel,
    including its subtracted Fourier representation. It is the constant
    background obtained by translating both kernel coordinates together.
    Short-range jet limits are supplied §1c premises. Fourier masses and all
    delta Jacobians are computed; the Abel limit is taken distributionally.
    No source 3-D row or source Fourier definition is a construction operand.
    """

    def __init__(self, reduction):
        self.r = reduction
        self.kn = sp.Symbol('s11cdSpectralNormalMomentum', real=True)
        self.shift = sp.Symbol('s11cdEndTranslation', real=True)
        units = PHYSICAL_METADATA.dimensions
        units.known[self.kn] = units.measure(reduction.normal_map[reduction.momentum_groups[0][2]])
        units.known[self.shift] = units.measure(reduction.z)
        self.constant_mass = sp.integrate(sp.cancel(reduction.abel_constant),
                                          (reduction.transfer, -sp.oo, sp.oo))
        self.delta_factors = []
        for s in reduction.normal_transfers:
            coefficient, factors = reduction.abel_constant.subs(reduction.transfer, s).as_coeff_Mul()
            self.delta_factors.append((factors.as_powers_dict(),
                                      self.constant_mass*sp.DiracDelta(s)/coefficient))

    def abel_constant_limit(self, value):
        """Replace complete constant Abel kernels, never individual poles."""
        @lru_cache(maxsize=None)
        def visit(node):
            if not node.has(self.r.regulator):
                return node
            if node.is_Mul:
                powers = dict(node.as_powers_dict())
                extra = sp.S.One
                for required, limit in self.delta_factors:
                    while all(powers.get(f, 0)*n > 0 and abs(powers.get(f, 0)) >= abs(n)
                              for f, n in required.items()):
                        for f, n in required.items():
                            powers[f] -= n
                        extra *= limit
                if extra != 1:
                    return extra*visit(sp.Mul(*(f**n for f, n in powers.items())))
            return node.func(*(visit(a) for a in node.args)) if node.args else node
        return visit(value)

    def background(self, value, end=None):
        """End limit on coefficients before weak normal convolution.

        End=None extracts the complete eta^0,sigma^0 coefficient. An end
        replaces the translated background as a whole, including the endpoint
        constants and localized remainder in its canonical subtraction; it
        does not discard a step transform by a pointwise Q-limit.
        """
        r = self.r
        if end is None:
            mapping = {r.symbols[n]: sp.S.Zero for n in ('eta_bg', 'sigma_W')}
        else:
            mapping = {v: r.end_values[(p, end)] for (p, _), v in r.end_values.items()}
        @lru_cache(maxsize=None)
        def visit(node):
            if node in mapping:
                return mapping[node]
            if end is not None and isinstance(node, AppliedUndef) and node.func in r.profiles.values():
                p = next(p for p, f in r.profiles.items() if f == node.func)
                return r.end_values[(p, end)]
            if not node.args or isinstance(node, Str):
                return node
            args = tuple(visit(a) for a in node.args)
            if isinstance(node, sp.Integral):
                return sp.Integral(*args) if args[0] != 0 else sp.S.Zero
            if isinstance(node, sp.Derivative) and args != node.args:
                return sp.diff(args[0], *args[1:])
            if isinstance(node, sp.Subs) and args != node.args:
                return sp.Subs(*args).doit()
            return node if args == node.args else node.func(*args)
        constant = visit(value)
        if end is not None:
            constant = self.abel_constant_limit(constant)
        return constant

    def profile_limit_operands(self, end):
        r = self.r
        return tuple((f((r.z+self.shift)/r.ell),
                      sp.Limit(f((r.z+self.shift)/r.ell), self.shift, end),
                      r.end_values[(p, end)]) for p, f in r.profiles.items())

    def phase_terms(self, value, coordinate):
        """Collect harmonic characters without expanding carrier algebra."""
        @lru_cache(maxsize=None)
        def visit(node):
            if not node.has(coordinate):
                return {sp.S.Zero: node}
            if node.func == sp.exp:
                frequency = sp.expand(sp.diff(node.args[0], coordinate)/sp.I)
                if frequency.has(coordinate):
                    raise NotImplementedError(('nonlinear normal phase', node))
                return {frequency: sp.exp(sp.expand(node.args[0]-sp.I*frequency*coordinate))}
            if node.is_Add:
                result = {}
                for a in node.args:
                    for q, c in visit(a).items():
                        result[q] = result.get(q, sp.S.Zero)+c
                return result
            if node.is_Mul:
                result = {sp.S.Zero: sp.S.One}
                for a in node.args:
                    product = {}
                    for q, c in result.items():
                        for k, b in visit(a).items():
                            qk = sp.expand(q+k)
                            product[qk] = product.get(qk, sp.S.Zero)+c*b
                    result = product
                return result
            if isinstance(node, sp.Piecewise):
                branches = [(visit(expr), condition) for expr, condition in node.args]
                support = set().union(*(set(c) for c, _ in branches))
                return {q: sp.Piecewise(*((c.get(q, sp.S.Zero), condition)
                                         for c, condition in branches)) for q in support}
            raise NotImplementedError(('uncontracted normal coordinate', node.func))
        return visit(value)

    def branches(self, expression, condition=sp.S.true):
        """Finite outgoing-sheet branches; conditions travel with contractions."""
        pieces = expression.atoms(sp.Piecewise)
        if not pieces:
            yield expression, condition
            return
        piece = max(pieces, key=dag_size)
        prior = sp.S.false
        for body, guard in piece.args:
            active = sp.And(condition, guard, sp.Not(prior))
            if active is not sp.S.false:
                yield from self.branches(memo_xreplace(expression, {piece: body}), active)
            prior = sp.Or(prior, guard)

    def integral_symbol(self, integral):
        r = self.r
        variables = tuple(l[0] for l in integral.limits)
        if any(tuple(l[1:]) != (-sp.oo, sp.oo) for l in integral.limits):
            raise NotImplementedError(('non-Fourier end integral', integral.limits))
        momenta = tuple(v for v in variables if v != r.zp)
        if any(v not in r.normal_map.values() for v in momenta):
            raise NotImplementedError(('end integral coordinates', integral.limits))
        evaluated = []
        for body, condition in self.branches(integral.function):
            terms = self.phase_terms(body, r.zp) if r.zp in variables else {None: body}
            pieces = []
            branch_conditions = []
            for frequency, coefficient in terms.items():
                if frequency is not None:
                    coefficient *= r.fourier_mass*sp.DiracDelta(frequency)
                deltas = tuple(sorted(coefficient.atoms(sp.DiracDelta), key=sp.default_sort_key))
                for powers, multiplier in polynomial_terms(coefficient, deltas):
                    constraints = []
                    for power, delta in zip(powers, deltas):
                        if power not in (0, 1) or len(delta.args) != 1:
                            raise NotImplementedError(('distribution product', power, delta))
                        if power:
                            constraints.append(delta.args[0])
                    matrix, rhs = sp.linear_eq_to_matrix(constraints, momenta)
                    if matrix.shape != (len(momenta), len(momenta)):
                        raise NotImplementedError(('normal constraint count', matrix.shape, momenta))
                    determinant = matrix.det()
                    solution = sp.linsolve((matrix, rhs), momenta)
                    roots = list(solution)
                    if len(roots) != 1 or determinant == 0:
                        raise NotImplementedError(('normal constraint rank', determinant, solution))
                    substitution = dict(zip(momenta, roots[0]))
                    pieces.append(multiplier.xreplace(substitution)/sp.Abs(determinant))
                    branch_conditions.append(condition.xreplace(substitution))
            if pieces:
                if len(set(branch_conditions)) != 1:
                    raise NotImplementedError(('different branch supports', branch_conditions))
                evaluated.append((sum(pieces), branch_conditions[0]))
        return sp.Piecewise(*evaluated) if evaluated else sp.S.Zero

    def symbol(self, value, mapping):
        """Evaluate the constant reduced action on a normal harmonic ansatz."""
        plane = dag_substitute(value, mapping)
        @lru_cache(maxsize=None)
        def visit(node):
            if not node.args or isinstance(node, Str):
                return node
            if isinstance(node, sp.Limit):
                return node  # Profile end values are the supplied shape data.
            args = tuple(visit(a) for a in node.args)
            if isinstance(node, sp.Integral):
                return self.integral_symbol(sp.Integral(*args)) if args[0] != 0 else sp.S.Zero
            return node if args == node.args else node.func(*args)
        return visit(plane)

    def strong_matrix(self, pencil, background):
        columns = []
        epsilon = self.r.symbols['epsilon_shape']
        def strip_phase(expression):
            terms = self.phase_terms(sp.diff(expression, epsilon), self.r.z)
            other = {q: c for q, c in terms.items() if q != self.kn and c != 0}
            if other:
                raise NotImplementedError(('end translation character', tuple(other)))
            return terms.get(self.kn, sp.S.Zero)
        for j in range(len(pencil.fields)):
            mapping = {f: (lambda z, i=i: sp.exp(sp.I*self.kn*z) if i == j else sp.S.Zero)
                       for i, f in enumerate(pencil.fields)}
            column = self.symbol(background, mapping)
            column = map_leaves(column, strip_phase)
            columns.append(column)
        return sp.ImmutableMatrix.hstack(*(sp.ImmutableMatrix(c) for c in columns))

    def translate_kernel(self, background):
        """Translate both normal kernel coordinates, keeping the trial fixed."""
        replacements = {self.r.z: self.r.z+self.shift, self.r.zp: self.r.zp+self.shift}
        @lru_cache(maxsize=None)
        def visit(node):
            if isinstance(node, (AppliedUndef, sp.Derivative, sp.Limit)):
                return node
            if node in replacements:
                return replacements[node]
            if isinstance(node, sp.Integral):
                return sp.Integral(visit(node.function), *node.limits)
            args = tuple(visit(a) for a in node.args)
            return node if args == node.args else node.func(*args)
        return visit(background)

    def strong_matrix_units(self, pencil):
        units = PHYSICAL_METADATA.dimensions
        return {(5*i+j,): tuple(a-b for a, b in zip(units.measure(pencil.strong[i]),
                                                   units.known[pencil.fields[j]]))
                for i in range(5) for j in range(5)}

    @staticmethod
    def block_entry(blocks, i, j):
        labels = ('THETA', 'E_W', 'DIV_U')
        if i < 3 and j < 3:
            return named(blocks, 'TT')
        if i < 3:
            return named(named(blocks, 'TH'), labels[j-3])
        if j < 3:
            return named(named(blocks, 'HT'), labels[i-3])
        return named(named(named(blocks, 'HH'), labels[j-3]), labels[i-3])

    def weak_matrix(self, blocks):
        """Full canonical symbol in the inherited curl-potential basis.

        All three transverse potentials are retained. Their curl gauge is
        emitted separately; this redundant representation is not used as a
        determinant-based mode solve.
        """
        labels = ('A0', 'A1', 'A2', 'Theta', 'E', 'Phi')
        trials = tuple(sp.Function('s11cdReducedTrial'+s) for s in labels)
        tests = tuple(sp.Function('s11cdReducedTest'+s) for s in labels)
        amplitudes = sp.symbols('s11cdSpectralTest0:6')
        units = PHYSICAL_METADATA.dimensions
        for f, amplitude in zip(tests, amplitudes):
            units.known[amplitude] = units.known[f]
        epsilon = self.r.symbols['epsilon_shape']
        columns = []
        for j in range(6):
            mapping = {f: (lambda z, i=i: sp.exp(sp.I*self.kn*z) if i == j else sp.S.Zero)
                       for i, f in enumerate(trials)}
            mapping.update({f: (lambda z, i=i: amplitudes[i]*sp.exp(-sp.I*self.kn*z))
                            for i, f in enumerate(tests)})
            evaluated = self.symbol(blocks, mapping)
            column = []
            for i in range(6):
                expression = sp.diff(self.block_entry(evaluated, i, j), amplitudes[i], epsilon)
                terms = self.phase_terms(expression, self.r.z)
                other = {q: c for q, c in terms.items() if q != 0 and c != 0}
                if other:
                    raise NotImplementedError(('weak end translation character', tuple(other)))
                column.append(terms.get(sp.S.Zero, sp.S.Zero))
            columns.append(column)
        return sp.ImmutableMatrix.hstack(*(sp.ImmutableMatrix(c) for c in columns))

    def weak_matrix_units(self, original_blocks):
        units = PHYSICAL_METADATA.dimensions
        labels = ('A0', 'A1', 'A2', 'Theta', 'E', 'Phi')
        result = {}
        for i in range(6):
            for j in range(6):
                source = units.measure(self.block_entry(original_blocks, i, j))
                trial = units.known[sp.Function('s11cdReducedTrial'+labels[j])]
                test = units.known[sp.Function('s11cdReducedTest'+labels[i])]
                result[(6*i+j,)] = tuple(a-b-c for a, b, c in zip(source, trial, test))
        return result

    def curl_gauge_operands(self, pencil):
        """Differentiate the actual curl trial ansatz; expose its null space."""
        amplitudes = sp.symbols('s11cdGaugePotential0:3')
        potentials = {sp.Function('s11cdReducedTrialA'+str(i)):
                      (lambda z, i=i: amplitudes[i]*sp.exp(sp.I*self.kn*z)) for i in range(3)}
        displacement = sp.Tuple(*(pencil.trial_ansatz('TRANSVERSE')[f](self.r.z)
                                  for f in pencil.fields[:3]))
        transformed = dag_substitute(displacement, potentials)
        curl = sp.ImmutableMatrix(3, 3, lambda i, j:
            self.phase_terms(sp.diff(transformed[i], amplitudes[j]), self.r.z).get(self.kn, sp.S.Zero))
        nullspace = curl.nullspace()
        return curl, tuple(nullspace), tuple(curl*v for v in nullspace)

    def field_lift(self, pencil):
        """Differentiate the existing sector ansatz into physical field columns."""
        labels = ('A0', 'A1', 'A2', 'Theta', 'E', 'Phi')
        amplitudes = sp.symbols('s11cdFieldLiftAmplitude0:6')
        trials = {sp.Function('s11cdReducedTrial'+label):
                  (lambda z, j=j: amplitudes[j]*sp.exp(sp.I*self.kn*z))
                  for j, label in enumerate(labels)}
        ansatz = [pencil.trial_ansatz(sector)
                  for sector in ('TRANSVERSE', 'THETA', 'E_W', 'LONGITUDINAL')]
        fields = sp.Tuple(*(sum(a[f](self.r.z) for a in ansatz) for f in pencil.fields))
        wave = dag_substitute(fields, trials)
        return sp.ImmutableMatrix(5, 6, lambda i, j:
            self.phase_terms(sp.diff(wave[i], amplitudes[j]), self.r.z).get(self.kn, sp.S.Zero))


class UniformSlabCurrent:
    """S11b conservative boundary work after the tangential energy reduction.

    The material virtual constraint is applied to the reduced stored-energy
    variation. Normal integration by parts computes its boundary work, then
    the harmonic time translation computes the slab energy-current bilinear.
    Bulk/face responses never enter this varied functional. Their contribution
    is a separate balance-law construction, not part of this slab operand.
    """

    def __init__(self, reduction, energy_row, ends, mass_symbol):
        self.r, self.ends = reduction, ends
        self.mass_symbol = mass_symbol
        self.energy_cases = {str(case): named(payload, 'VALUE')
                             for case, payload in energy_row['value']}
        r = reduction
        self.phase_coordinate = sp.Symbol('s11cdCurrentPhase', real=True)
        self.virtual_parameter = sp.Symbol('s11cdCurrentVariationParameter', real=True)
        self.leg_momenta = sp.symbols('s11cdCurrentLeftMomentum s11cdCurrentRightMomentum', real=True)
        self.fields = tuple(tuple(sp.Function('s11cdCurrent'+side+name)
                                  for name in ('U1', 'U2', 'U3', 'Theta', 'E'))
                            for side in ('Plus', 'Minus'))
        self.variations = tuple(tuple(sp.Function('s11cdCurrentVariation'+side+name)
                                      for name in ('U1', 'U2', 'U3', 'E'))
                                for side in ('Plus', 'Minus'))
        self.amplitudes = tuple(sp.symbols('s11cdCurrent'+side+'Amplitude0:5')
                                for side in ('Plus', 'Minus'))
        units = PHYSICAL_METADATA.dimensions
        units.known[self.phase_coordinate] = units.zero
        units.known[self.virtual_parameter] = units.zero
        for momentum in self.leg_momenta:
            units.known[momentum] = units.measure(ends.kn)
        field_units = [units.known[sp.Function('s11cdReducedField'+name)]
                       for name in ('u1', 'u2', 'u3', 'theta', 'eW')]
        for functions, variations, amplitudes in zip(self.fields, self.variations, self.amplitudes):
            for f, a, unit in zip(functions, amplitudes, field_units):
                units.known[f] = unit
                units.known[a] = unit
            for f, unit in zip(variations, (*field_units[:3], field_units[4])):
                units.known[f] = unit
        self.field_units = field_units
        self.phase = sp.exp(sp.I*(sum(k*x for k, x in zip(r.tangents, r.x[:2]))
                                 -r.omega*r.t+self.phase_coordinate))

    @lru_cache(maxsize=None)
    def retained(self, value):
        """Compute the rectangular first-background-grade Taylor projection."""
        eta, sigma = (self.r.symbols[n] for n in ('eta_bg','sigma_W'))
        return sp.expand(sum(sp.diff(value,eta,a,sigma,b).subs({eta:0,sigma:0})*
                             eta**a*sigma**b for a in range(2) for b in range(2)))

    @lru_cache(maxsize=None)
    def construct(self, anchoring, end):
        r, z = self.r, self.r.z
        body = self.energy_cases[anchoring]
        source = sp.Add(*(row[4] if str(row[0]) in ('W_BG', 'MU_R_BG') else row[3]
                          for row in body[1:]))
        background = {}
        for atom in source.atoms(sp.Symbol):
            match = re.fullmatch(r'([wm])1_profile((?:_d[123](?:d[123])*)?)', atom.name)
            if match:
                constant = r.end_values[(match[1], end)] if end is not None else sp.S.Zero
                indices = re.findall(r'd([123])', match[2])
                background[atom] = sp.diff(constant, r.xi, len(indices)) if indices else constant
        if end is None:
            background.update({r.symbols[n]: sp.S.Zero for n in ('eta_bg', 'sigma_W')})
        uniform_source = source.xreplace(background)
        # Field carriers are an explicit +/- harmonic ansatz. Their spatial
        # factors are differentiated before selecting the zero phase harmonic.
        fields = tuple((plus(z)*self.phase+minus(z)/self.phase)/2
                       for plus, minus in zip(*self.fields))
        source_names = {'u_1': 0, 'u_2': 1, 'u_3': 2, 'theta': 3, 'e_W': 4}
        field_map, parameter_map = {}, {}
        for atom in source.atoms(sp.Symbol):
            name = re.sub(r'^grad_theta_([123])$', r'theta_d\1', atom.name)
            match = re.fullmatch(r'(u_[123]|theta|e_W)((?:_t{1,2})?(?:_?d[123])*)', name)
            if match:
                field, suffix = match.groups()
                value = fields[source_names[field]]
                for index in re.findall(r'd([123])', suffix):
                    value = sp.diff(value, r.x[int(index)-1] if index != '3' else z)
                if '_tt' in suffix:
                    value = sp.diff(value, r.t, 2)
                elif '_t' in suffix:
                    value = sp.diff(value, r.t)
                field_map[atom] = value
                PHYSICAL_METADATA.dimensions.known[atom] = PHYSICAL_METADATA.dimensions.measure(value)
            elif atom.name in r.symbols:
                parameter_map[atom] = r.symbols[atom.name]
                PHYSICAL_METADATA.dimensions.known[atom] = PHYSICAL_METADATA.dimensions.measure(parameter_map[atom])
        harmonic_energy = uniform_source.xreplace(parameter_map).xreplace(field_map)
        energy = sp.expand(self.ends.phase_terms(sp.expand(harmonic_energy),
                                                self.phase_coordinate).get(0, sp.S.Zero))
        # The inherited field is measured against W_0. At a changed end it
        # need not be the local-background thickness fraction. Extract the
        # homogeneous material constraint from the actual mass row with the
        # two supplied mass-transfer drivers disabled, before varying energy.
        no_transfer = {r.symbols[n]:sp.S.Zero for n in ('Lambda_A_0','Lambda_V_0')}
        mass_control = self.mass_symbol.subs(no_transfer).applyfunc(sp.cancel)
        constraint_coefficients = tuple(self.retained(sp.cancel(v/mass_control[3])) for v in mass_control)
        constraint_residual = sp.ImmutableMatrix(1,5,lambda i,j:
            sp.expand(constraint_coefficients[j]*mass_control[3]-mass_control[j]))
        virtual = {}
        for sign, functions, variations in zip((1, -1), self.fields, self.variations):
            displacement = [f(z) for f in variations[:3]]
            thickness = variations[3](z)
            density = sp.S.Zero
            for index, test in zip((0,1,2,4), (*displacement,thickness)):
                coefficient = constraint_coefficients[index].subs(
                    {a:sign*a for a in (*r.tangents,r.omega)}, simultaneous=True)
                # SymPy Limit carries unknown commutativity. It is an end
                # coefficient here, not a normal differential operator.
                constants = {limit:sp.Dummy(real=True) for limit in coefficient.atoms(sp.Limit)}
                restore = {v:k for k,v in constants.items()}
                polynomial = sp.Poly(coefficient.xreplace(constants), self.ends.kn)
                density -= sum(c.xreplace(restore)*sp.diff(test,z,power[0])/sp.I**power[0]
                               for power,c in polynomial.terms())
            virtual.update({f(z): f(z)+self.virtual_parameter*v
                            for f, v in zip(functions, (*displacement, density, thickness))})
        action = -energy
        variation = sp.expand(sp.diff(action.subs(virtual, simultaneous=True).doit(),
                                      self.virtual_parameter).subs(self.virtual_parameter, 0))
        bulk, boundary, reconstruction = sp.S.Zero, sp.S.Zero, sp.S.Zero
        for test in (f(z) for group in self.variations for f in group):
            jets = {test} | {d for d in variation.atoms(sp.Derivative) if d.expr == test}
            for jet in sorted(jets, key=sp.default_sort_key):
                order = sum(n for variable, n in jet.variable_count) if isinstance(jet, sp.Derivative) else 0
                coefficient = sp.diff(variation, jet)
                reconstruction += coefficient*jet
                bulk += (-1)**order*sp.diff(coefficient, z, order)*test
                boundary += sum((-1)**j*sp.diff(coefficient, z, j)*sp.diff(test, z, order-1-j)
                                for j in range(order))
        variation_residual = sp.expand(variation-reconstruction)
        boundary_residual = sp.expand(reconstruction-bulk-sp.diff(boundary, z))
        time_shift = {}
        for sign, functions, variations in zip((1, -1), self.fields, self.variations):
            factor = sp.cancel(sp.diff(self.phase**sign, r.t)/(self.phase**sign))
            time_shift.update({v(z): factor*f(z) for v, f in
                               zip(variations, (*functions[:3], functions[4]))})
        untruncated_current = sp.expand(boundary.subs(time_shift, simultaneous=True).doit())
        current = self.retained(untruncated_current)
        right_momentum, left_momentum = self.leg_momenta[1], self.leg_momenta[0]
        wave = {}
        for sign, momentum, functions, amplitudes in zip(
                (1, -1), (right_momentum, left_momentum), self.fields, self.amplitudes):
            wave.update({f(z): a*sp.exp(sign*sp.I*momentum*z)
                         for f, a in zip(functions, amplitudes)})
        polarized = sp.expand(current.subs(wave, simultaneous=True).doit().subs(z, 0))
        matrix = sp.ImmutableMatrix(5, 5, lambda i, j:
            sp.diff(polarized, self.amplitudes[1][i], self.amplitudes[0][j]))
        return {'SOURCE_ENERGY': source, 'UNIFORM_SOURCE_ENERGY': uniform_source,
                'SOURCE_PARAMETER_ALIGNMENT': tuple(parameter_map.items()),
                'ZERO_TRANSFER_MASS_ROW': mass_control,
                'MATERIAL_CONSTRAINT_COEFFICIENTS': constraint_coefficients,
                'MATERIAL_CONSTRAINT_TRUNCATION_REMAINDER': constraint_residual,
                'MATERIAL_CONSTRAINT_RETAINED_RESIDUAL': constraint_residual.applyfunc(self.retained),
                'TANGENTIAL_ENERGY_REDUCTION': energy, 'VIRTUAL_VARIATION': variation,
                'VARIATION_RECONSTRUCTION_RESIDUAL': variation_residual,
                'NORMAL_BOUNDARY_WORK': boundary, 'NORMAL_BOUNDARY_RESIDUAL': boundary_residual,
                'SLAB_CURRENT_BEFORE_BACKGROUND_PROJECTION': untruncated_current,
                'SLAB_CURRENT_BACKGROUND_REMAINDER': sp.expand(untruncated_current-current),
                'SLAB_CURRENT': current, 'SLAB_CURRENT_MATRIX': matrix}

    def emit(self, anchoring, end, suffix):
        result = self.construct(anchoring, end)
        dimensions = PHYSICAL_METADATA.dimensions
        energy_unit = dimensions.measure(result['TANGENTIAL_ENERGY_REDUCTION'])
        current_unit = dimensions.measure(result['SLAB_CURRENT'])
        for key, value in result.items():
            zero_units = None
            if key in ('VARIATION_RECONSTRUCTION_RESIDUAL', 'NORMAL_BOUNDARY_RESIDUAL'):
                zero_units = {(): energy_unit}
            elif key == 'SLAB_CURRENT_MATRIX':
                zero_units = {(5*i+j,): tuple(a-b-c for a,b,c in
                    zip(current_unit, self.field_units[i], self.field_units[j]))
                    for i in range(5) for j in range(5)}
            elif key in ('MATERIAL_CONSTRAINT_TRUNCATION_REMAINDER', 'MATERIAL_CONSTRAINT_RETAINED_RESIDUAL'):
                zero_units = {(i,):dimensions.measure(result['ZERO_TRANSFER_MASS_ROW'][i]) for i in range(5)}
            elif key == 'SLAB_CURRENT_BACKGROUND_REMAINDER':
                zero_units = {():current_unit}
            physical('CONSERVATIVE_'+key+'_'+suffix, value, zero_dimensions=zero_units)
        physical('CONSERVATIVE_CURRENT_NORMAL_LEGS_'+suffix, self.leg_momenta)
        return result


class SlabEnergyBalance:
    """Actual energy boundary work with the material mass-rate defect retained.

    Only the reduced stored energy is varied. The material constraint acts on
    virtual tests, while the two independent density time rates remain live
    in the energy balance. No response is differentiated with respect to a
    field. The conservative current is an independently reconstructed operand.
    """

    def __init__(self, conservative):
        self.c = conservative
        self.r = conservative.r
        self.tests = tuple(tuple(sp.Function('s11cdBalanceVariation'+side+name)
                                 for name in ('U1', 'U2', 'U3', 'Theta', 'E'))
                           for side in ('Plus', 'Minus'))
        dims = PHYSICAL_METADATA.dimensions
        for group in self.tests:
            for f, unit in zip(group, conservative.field_units):
                dims.known[f] = unit

    def parts(self, expression, tests):
        """Coefficient extraction followed by normal integration by parts."""
        z = self.r.z
        euler, boundary, reconstructed = [], sp.S.Zero, sp.S.Zero
        for test in tests:
            row = sp.S.Zero
            jets = {test} | {d for d in expression.atoms(sp.Derivative) if d.expr == test}
            for jet in sorted(jets, key=sp.default_sort_key):
                order = sum(n for _, n in jet.variable_count) if isinstance(jet, sp.Derivative) else 0
                coefficient = sp.diff(expression, jet)
                reconstructed += coefficient*jet
                row += (-1)**order*sp.diff(coefficient, z, order)
                boundary += sum((-1)**j*sp.diff(coefficient, z, j)*sp.diff(test, z, order-1-j)
                                for j in range(order))
            euler.append(sp.expand(row))
        bulk = sp.Add(*(a*b for a, b in zip(euler, tests)))
        return (sp.ImmutableMatrix(euler), sp.expand(boundary),
                sp.expand(expression-reconstructed),
                sp.expand(reconstructed-bulk-sp.diff(boundary, z)))

    def virtual_map(self, coefficients):
        c, r, z = self.c, self.r, self.r.z
        virtual = {}
        for sign, group in zip((1, -1), self.tests):
            density = sp.S.Zero
            for i in (0, 1, 2, 4):
                coefficient = coefficients[i].subs(
                    {a:sign*a for a in (*r.tangents, r.omega)}, simultaneous=True)
                constants = {v:sp.Dummy(real=True) for v in coefficient.atoms(sp.Limit)}
                restore = {v:k for k, v in constants.items()}
                polynomial = sp.Poly(coefficient.xreplace(constants), c.ends.kn)
                density -= sum(a.xreplace(restore)*sp.diff(group[i](z), z, power[0])/sp.I**power[0]
                               for power, a in polynomial.terms())
            virtual[group[3](z)] = density
        return virtual

    def matrix(self, expression):
        c, z = self.c, self.r.z
        wave = {}
        for sign, momentum, functions, amplitudes in zip(
                (1, -1), (c.leg_momenta[1], c.leg_momenta[0]), c.fields, c.amplitudes):
            wave.update({f(z):a*sp.exp(sign*sp.I*momentum*z) for f, a in zip(functions, amplitudes)})
        polarized = sp.expand(expression.subs(wave, simultaneous=True).doit().subs(z, 0))
        matrix = sp.ImmutableMatrix(5, 5, lambda i, j:
            sp.diff(polarized, c.amplitudes[1][i], c.amplitudes[0][j]))
        reconstructed = (sp.ImmutableMatrix(1, 5, c.amplitudes[1])*matrix*
                         sp.ImmutableMatrix(c.amplitudes[0]))[0]
        return matrix, sp.expand(polarized-reconstructed)

    @lru_cache(maxsize=None)
    def construct(self, anchoring, end):
        c, r, z = self.c, self.r, self.r.z
        source = c.construct(anchoring, end)
        tests = tuple(f(z) for group in self.tests for f in group)
        variation_map = {f(z):f(z)+c.virtual_parameter*v(z)
                         for fields, group in zip(c.fields, self.tests) for f, v in zip(fields, group)}
        variation = sp.expand(sp.diff((-source['TANGENTIAL_ENERGY_REDUCTION']).subs(
            variation_map, simultaneous=True).doit(), c.virtual_parameter).subs(c.virtual_parameter, 0))
        euler, boundary, reconstruction, ibp = self.parts(variation, tests)
        virtual = self.virtual_map(source['MATERIAL_CONSTRAINT_COEFFICIENTS'])
        mechanical_tests = tuple(f(z) for group in self.tests for i, f in enumerate(group) if i != 3)
        virtual_bulk = sp.expand(sp.Add(*(a*b for a, b in zip(euler, tests))).subs(
            virtual, simultaneous=True).doit())
        mechanical, transport, transport_reconstruction, transport_ibp = self.parts(virtual_bulk, mechanical_tests)
        virtual_boundary = sp.expand(boundary.subs(virtual, simultaneous=True).doit()+transport)
        old_tests = {f(z):v(z) for group, old in zip(self.tests, c.variations)
                     for f, v in zip((*group[:3], group[4]), old)}
        virtual_variation_residual = sp.expand(variation.subs(virtual, simultaneous=True).doit().subs(
            old_tests, simultaneous=True).doit()-source['VIRTUAL_VARIATION'])
        virtual_boundary_residual = sp.expand(virtual_boundary.subs(old_tests, simultaneous=True).doit()
                                             -source['NORMAL_BOUNDARY_WORK'])
        rates = {v(z):sp.cancel(sp.diff(c.phase**sign, r.t)/(c.phase**sign))*f(z)
                 for sign, fields, group in zip((1, -1), c.fields, self.tests) for f, v in zip(fields, group)}
        defects = tuple(sp.expand((group[3](z)-virtual[group[3](z)]).subs(rates, simultaneous=True).doit())
                        for group in self.tests)
        actual_boundary = sp.expand(boundary.subs(rates, simultaneous=True).doit())
        chemical_transport = sp.expand(transport.subs(rates, simultaneous=True).doit())
        current = c.retained(actual_boundary+chemical_transport)
        defect_map = {test:sp.S.Zero for test in tests}
        defect_map.update({group[3](z):defect for group, defect in zip(self.tests, defects)})
        correction = c.retained(boundary.subs(defect_map, simultaneous=True).doit())
        current_residual = sp.expand(current-source['SLAB_CURRENT']-correction)
        matrix, polarization_residual = self.matrix(current)
        correction_matrix, correction_polarization_residual = self.matrix(correction)
        chemical_work = sp.expand(sum(-euler[5*i+3]*defect for i, defect in enumerate(defects)))
        rate_variation = sp.expand(variation.subs(rates, simultaneous=True).doit())
        mechanical_power = sp.expand(sum(a*b for a, b in zip(mechanical, mechanical_tests)).subs(
            rates, simultaneous=True).doit())
        balance_residual = c.retained(rate_variation-mechanical_power+chemical_work-sp.diff(current, z))
        # The coefficient connecting an averaged quadratic variation to its
        # positive harmonic is extracted from the same real-field ansatz.
        probe = (c.amplitudes[0][3]*c.phase+c.amplitudes[1][3]/c.phase)/2
        norm = sp.diff(c.ends.phase_terms(sp.expand(probe**2/2), c.phase_coordinate)[0],
                       c.amplitudes[0][3], c.amplitudes[1][3])
        chemical = sp.expand(-euler[8]/norm).coeff(r.symbols['epsilon_shape'], 2)
        chemical_wave = chemical.subs({f(z):a*sp.exp(sp.I*c.ends.kn*z)
                                      for f, a in zip(c.fields[0], c.amplitudes[0])}, simultaneous=True)
        chemical_wave = sp.expand(chemical_wave.doit().subs(z, 0))
        chemical_row = sp.ImmutableMatrix(1, 5, lambda i, j:sp.diff(chemical_wave, c.amplitudes[0][j]))
        chemical_residual = sp.expand(chemical_wave-(chemical_row*sp.ImmutableMatrix(c.amplitudes[0]))[0])
        return {'UNCONSTRAINED_VARIATION':variation, 'ENERGY_EULER_DERIVATIVES':-euler,
                'UNCONSTRAINED_BOUNDARY_WORK':boundary,
                'VARIATION_RECONSTRUCTION_RESIDUAL':reconstruction,
                'ENERGY_IBP_RESIDUAL':ibp, 'MATERIAL_VIRTUAL_DENSITY_MAP':tuple(virtual.items()),
                'CONSTRAINED_MECHANICAL_EULER_DERIVATIVES':mechanical,
                'CHEMICAL_TRANSPORT_BOUNDARY':transport,
                'TRANSPORT_RECONSTRUCTION_RESIDUAL':transport_reconstruction,
                'TRANSPORT_IBP_RESIDUAL':transport_ibp,
                'CONSERVATIVE_VARIATION_RECONSTRUCTION_RESIDUAL':virtual_variation_residual,
                'CONSERVATIVE_BOUNDARY_RECONSTRUCTION_RESIDUAL':virtual_boundary_residual,
                'MATERIAL_DENSITY_RATE_DEFECTS':defects, 'ACTUAL_TIME_BOUNDARY':actual_boundary,
                'CHEMICAL_TRANSPORT_CURRENT':chemical_transport,
                'MASS_RATE_CHEMICAL_WORK':chemical_work, 'MECHANICAL_POWER':mechanical_power,
                'ACTUAL_ENERGY_RATE_VARIATION':rate_variation,
                'MASS_RATE_BOUNDARY_CORRECTION':correction, 'SLAB_CURRENT':current,
                'SLAB_CURRENT_MATRIX':matrix, 'MASS_RATE_CORRECTION_MATRIX':correction_matrix,
                'CURRENT_DECOMPOSITION_RESIDUAL':current_residual,
                'POLARIZATION_RESIDUAL':polarization_residual,
                'CORRECTION_POLARIZATION_RESIDUAL':correction_polarization_residual,
                'ENERGY_BALANCE_RESIDUAL':balance_residual,
                'HARMONIC_VARIATION_NORMALIZATION':norm,
                'CHEMICAL_FUNCTIONAL_DERIVATIVE':chemical,
                'CHEMICAL_FIELD_ROW':chemical_row,
                'CHEMICAL_POLARIZATION_RESIDUAL':chemical_residual}

    def emit(self, anchoring, end, suffix):
        result = self.construct(anchoring, end)
        dims = PHYSICAL_METADATA.dimensions
        energy = dims.measure(self.c.construct(anchoring, end)['TANGENTIAL_ENERGY_REDUCTION'])
        current = dims.measure(result['SLAB_CURRENT'])
        power = tuple(a+b for a, b in zip(energy, dims.measure(self.r.omega)))
        for key, value in result.items():
            units = None
            if key.endswith('_MATRIX'):
                units = {(5*i+j,):tuple(a-b-c for a, b, c in zip(current, self.c.field_units[i], self.c.field_units[j]))
                         for i in range(5) for j in range(5)}
            elif key == 'ENERGY_EULER_DERIVATIVES':
                units = {(i,):tuple(a-b for a, b in zip(energy, self.c.field_units[i % 5])) for i in range(10)}
            elif key == 'CONSTRAINED_MECHANICAL_EULER_DERIVATIVES':
                units = {(i,):tuple(a-b for a, b in zip(energy, self.c.field_units[j]))
                         for i, j in enumerate((0, 1, 2, 4)*2)}
            elif key == 'CHEMICAL_FIELD_ROW':
                units = {(j,):tuple(a-b-c for a, b, c in zip(energy, self.c.field_units[3], self.c.field_units[j]))
                         for j in range(5)}
            elif key == 'ACTUAL_ENERGY_RATE_VARIATION':
                units = {():power}
            elif key.endswith('_RESIDUAL'):
                unit = (power if key == 'ENERGY_BALANCE_RESIDUAL' else current if key in (
                    'CURRENT_DECOMPOSITION_RESIDUAL', 'POLARIZATION_RESIDUAL',
                    'CORRECTION_POLARIZATION_RESIDUAL') else
                    tuple(a+b for a, b in zip(energy, dims.measure(self.r.z))) if
                    key == 'CONSERVATIVE_BOUNDARY_RECONSTRUCTION_RESIDUAL' else
                    tuple(a-b for a, b in zip(energy, self.c.field_units[3])) if
                    key == 'CHEMICAL_POLARIZATION_RESIDUAL' else energy)
                units = {():unit}
            physical('ENERGY_BALANCE_'+key+'_'+suffix, value, zero_dimensions=units)
        return result


class ClosedAcousticEnergy:
    """Acoustic balance and closed face lift on an already reduced end."""

    def __init__(self, balance, modes, strong):
        self.balance, self.modes, self.strong = balance, modes, strong
        self.c, self.r = balance.c, balance.r
        self.depth = sp.Symbol('s11cdAcousticOutwardDepth', nonnegative=True)
        self.height = sp.Symbol('s11cdAcousticDepthCutoff', positive=True)
        self.qlegs = sp.symbols('s11cdAcousticLeftNormalMomentum s11cdAcousticRightNormalMomentum', complex=True)
        self.alegs = sp.symbols('s11cdAcousticLeftAmplitude s11cdAcousticRightAmplitude')
        self.phi = sp.Function('s11cdAcousticPotential')
        self.energy_coefficients = sp.symbols('s11cdAcousticTimeEnergyCoefficient s11cdAcousticGradientEnergyCoefficient')
        dims = PHYSICAL_METADATA.dimensions
        length, frequency = dims.measure(self.r.z), dims.measure(self.r.omega)
        potential = tuple(2*a+b for a, b in zip(length, frequency))
        density = dims.measure(self.r.symbols['rho_m'])
        dims.known[self.depth] = dims.known[self.height] = length
        dims.known[self.phi] = potential
        for q in self.qlegs:
            dims.known[q] = tuple(-a for a in length)
        for a in self.alegs:
            dims.known[a] = potential
        speed = dims.measure(self.r.symbols['c_s0'])
        dims.known[self.energy_coefficients[0]] = tuple(a-2*b for a, b in zip(density, speed))
        dims.known[self.energy_coefficients[1]] = density

    @staticmethod
    @lru_cache(maxsize=1024)
    def cancel_material_coefficient(value):
        numerator, denominator = sp.together(value).as_numer_denom()
        scale, numerator, denominator = sp.cancel((sp.expand(numerator), denominator))
        return scale*numerator/denominator

    def cancel_background_expression(self, value, label):
        """Exact rational reconstruction with all background orders retained.

        Field amplitudes and background powers are collected structurally.
        Background-dependent denominator factors stay explicit; they are not
        Taylor-expanded or treated as independent physical carriers. Only the
        remaining material coefficients enter multivariate cancellation.
        """
        grades = (self.modes.eta, self.modes.sigma)
        amplitudes = self.c.amplitudes[0]
        checkpoint = getattr(self, 'rational_checkpoint', lambda name, source, calculate:calculate())
        def reconstruct():
            amplitude_terms = []
            for ai, (powers, coefficient) in enumerate(polynomial_terms(value, amplitudes)):
                numerator, denominator = sp.together(coefficient).as_numer_denom()
                factors = sp.Mul.make_args(denominator)
                dynamic = sp.Mul(*(f for f in factors if f.has(*grades)))
                static = sp.Mul(*(f for f in factors if not f.has(*grades)))
                background_terms = []
                for exponents, material in polynomial_terms(numerator, grades):
                    operand = material/static
                    name = label+'_amplitude_'+str(ai)+'_grade_'+'_'.join(map(str,exponents))
                    reduced = checkpoint(name, operand,
                        lambda operand=operand:self.cancel_material_coefficient(operand))
                    background_terms.append(reduced*sp.prod(g**p for g,p in zip(grades,exponents)))
                amplitude_terms.append(sp.Add(*background_terms)/dynamic*
                    sp.prod(a**p for a,p in zip(amplitudes,powers)))
            return sp.Add(*amplitude_terms)
        return checkpoint(label+'_reconstruction', value, reconstruct)

    def cancel_background_row(self, row, label):
        return sp.ImmutableMatrix(row.rows, row.cols, [
            self.cancel_background_expression(v,label+'_'+str(i)) for i,v in enumerate(row)])

    @lru_cache(maxsize=None)
    def construct(self, anchoring, end):
        c, r, s = self.c, self.r, self.depth
        rho, speed = (r.symbols[n] for n in ('rho_m', 'c_s0'))
        positions = (*r.x[:2], r.z, s)
        phi = self.phi(r.t, *positions)
        velocity = sp.ImmutableMatrix([sp.diff(phi, x) for x in positions])
        pressure = -rho*sp.diff(phi, r.t)
        wave = sp.diff(phi, r.t, 2)-speed**2*sum(sp.diff(phi, x, 2) for x in positions)
        local_flux = pressure*velocity
        a, b = self.energy_coefficients
        energy_ansatz = a*sp.diff(phi, r.t)**2+b*(velocity.T*velocity)[0]
        balance_ansatz = sp.diff(energy_ansatz, r.t)+sum(sp.diff(local_flux[i], x) for i, x in enumerate(positions))
        wave_evolution = sp.solve(wave, sp.diff(phi, r.t, 2))[0]
        balance_on_wave = sp.expand(balance_ansatz.subs(sp.diff(phi, r.t, 2), wave_evolution))
        jets = sorted(balance_on_wave.atoms(sp.Derivative), key=sp.default_sort_key)
        coefficients = sp.solve(sp.Poly(balance_on_wave, *jets).coeffs(), (a, b), dict=True)[0]
        energy = energy_ansatz.subs(coefficients)
        local_residual = sp.expand(balance_on_wave.subs(coefficients))
        qleft, qright = self.qlegs
        aleft, aright = self.alegs
        kleft, kright = c.leg_momenta
        plus = aright*c.phase*sp.exp(sp.I*(kright*r.z+qright*s))
        minus = aleft/c.phase*sp.exp(-sp.I*(kleft*r.z+qleft*s))
        harmonic = r.symbols['epsilon_shape']*(plus+minus)/2
        def average(expression):
            return sp.expand(c.ends.phase_terms(sp.expand(expression.subs(phi, harmonic).doit()),
                                                c.phase_coordinate).get(0, sp.S.Zero))
        harmonic_energy = average(energy)
        harmonic_flux = local_flux.applyfunc(average)
        harmonic_wave = tuple(sp.cancel(wave.subs(phi, leg).doit()/leg) for leg in (minus, plus))
        depth_factor = sp.exp(sp.I*(qright-qleft)*s)
        current_coefficient = sp.simplify(harmonic_flux[2].subs(r.z, 0)/depth_factor)
        generic_depth = sp.integrate(depth_factor, (s, 0, self.height))
        diagonal_depth = sp.integrate(depth_factor.subs(qleft, qright), (s, 0, self.height))
        finite_depth = sp.Piecewise((generic_depth, sp.Ne(qright, qleft)), (diagonal_depth, True))
        finite_residual = sp.Piecewise((sp.simplify(sp.diff(generic_depth, self.height)-depth_factor.subs(s, self.height)),
                                       sp.Ne(qright, qleft)),
                                      (sp.simplify(sp.diff(diagonal_depth, self.height)-
                                                   depth_factor.subs(qleft, qright).subs(s, self.height)), True))
        decay = sp.Symbol('s11cdAcousticDepthDecay', positive=True)
        oscillation = sp.Symbol('s11cdAcousticDepthOscillation', real=True)
        dims = PHYSICAL_METADATA.dimensions
        dims.known[decay] = dims.known[oscillation] = dims.measure(qright)
        convergent_factor = sp.exp((-decay+sp.I*oscillation)*s)
        # Integrate the exponential using its computed differential
        # coefficient; verify the antiderivative before using its endpoints.
        primitive = convergent_factor/sp.cancel(sp.diff(convergent_factor, s)/convergent_factor)
        convergent_finite = primitive.subs(s, self.height)-primitive.subs(s, 0)
        boundary_modulus = sp.simplify(sp.Abs(primitive.subs(s, self.height)))
        boundary_modulus_limit = sp.limit(boundary_modulus, self.height, sp.oo)
        # The modulus limit bounds both real and imaginary endpoint parts.
        # A nonzero/unresolved bound is retained; it cannot define this limit.
        infinite_depth = (boundary_modulus_limit-primitive.subs(s, 0) if boundary_modulus_limit == 0
                          else sp.Limit(convergent_finite, self.height, sp.oo))
        primitive_residual = sp.simplify(sp.diff(primitive, s)-convergent_factor)
        depth_rate = sp.cancel(sp.diff(depth_factor, s)/depth_factor)
        depth_map = {decay:-sp.re(depth_rate), oscillation:sp.im(depth_rate)}
        depth_join = sp.simplify(sp.expand_complex((-decay+sp.I*oscillation).subs(depth_map)-depth_rate))

        slab = self.balance.construct(anchoring, end)
        chemical = (slab['CHEMICAL_FIELD_ROW']*sp.ImmutableMatrix(c.amplitudes[0]))[0]
        time_rate = sp.cancel(sp.diff(c.phase, r.t)/c.phase)
        surface_density = sp.cancel(c.mass_symbol.subs({r.symbols[n]:0 for n in ('Lambda_A_0', 'Lambda_V_0')})[3]/time_rate)
        mus = sp.cancel(chemical/surface_density)
        amplitude = self.alegs[1]
        outgoing = amplitude*c.phase*sp.exp(sp.I*(self.modes.k*r.z+qright*s))
        face_v = sp.cancel(sp.diff(outgoing, s).subs(s, 0)/(c.phase*sp.exp(sp.I*self.modes.k*r.z)))
        face_p = sp.cancel((-rho*sp.diff(outgoing, r.t)).subs(s, 0)/(c.phase*sp.exp(sp.I*self.modes.k*r.z)))
        kernels = {name:sp.cancel(r.symbols['Lambda_'+name+'_0']/(1-sp.I*r.omega*r.symbols['tau_'+name]))
                   for name in ('A', 'V', 'X')}
        faces = []
        for sign in (-1, 1):
            # Supplied face geometry and real-field harmonic ansatz.
            displacement = sign*r.symbols['W_0']*c.amplitudes[0][4]*c.phase/2
            outward_velocity = sp.cancel(sign*sp.diff(displacement, r.t)/c.phase)
            outward_virtual_lift = sp.diff(sign*displacement/c.phase, c.amplitudes[0][4])
            relative = rho*(face_v-outward_velocity)
            affinity = mus-face_p/rho
            closure = relative-kernels['A']*affinity-kernels['V']*outward_velocity
            coefficient = sp.diff(closure, amplitude)
            constant = closure.subs(amplitude, 0)
            solution = sp.cancel(-constant/coefficient)
            closed = {amplitude:solution}
            faces.append({'ORIENTATION':sign, 'DISPLACEMENT':displacement,
                          'OUTWARD_VELOCITY':outward_velocity, 'VIRTUAL_LIFT':outward_virtual_lift,
                          'AMPLITUDE':solution, 'PRESSURE':sp.cancel(face_p.subs(closed)),
                          'RELATIVE_MASS_FLUX':sp.cancel(relative.subs(closed)),
                          'AFFINITY':sp.cancel(affinity.subs(closed)),
                          'MECHANICAL_LOAD':self.cancel_background_expression(
                              ((face_p+kernels['X']*affinity)*outward_virtual_lift).subs(closed),
                              'face_'+str(sign)+'_mechanical_load'),
                          'AMPLITUDE_EQUATION_COEFFICIENT':coefficient,
                          'AMPLITUDE_EQUATION_RECONSTRUCTION_RESIDUAL':sp.expand(closure-coefficient*amplitude-constant),
                          'CLOSURE_RESIDUAL':sp.cancel(closure.subs(closed))})
        algebraic, relation, joins = self.modes.analytic(self.strong)
        acoustic_row = harmonic_wave[1].subs(kright, self.modes.k)
        radical_scale = sp.sqrt(sp.Poly(relation, self.modes.q).nth(2)/sp.Poly(acoustic_row, qright).nth(2))
        radical_map = {qright:radical_scale*self.modes.q}
        radical_residual = sp.simplify(acoustic_row.subs(radical_map)-relation)
        no_transfer = {r.symbols[n]:0 for n in ('Lambda_A_0', 'Lambda_V_0')}
        no_face = {**no_transfer, r.symbols['Lambda_X_0']:0}
        mass_increment = self.cancel_background_row(
            algebraic[3, :]-algebraic[3, :].subs(no_transfer), 'mass_increment')
        mechanical_increment = self.cancel_background_row(
            algebraic[4, :]-algebraic[4, :].subs(no_face).subs(rho, 0), 'mechanical_increment')
        mass_from_faces = sum(face['RELATIVE_MASS_FLUX'] for face in faces).subs(radical_map)
        mechanical_from_faces = sum(face['MECHANICAL_LOAD'] for face in faces).subs(radical_map)
        mass_face_row = sp.ImmutableMatrix(1, 5, lambda i, j:c.retained(
            self.cancel_background_expression(sp.diff(mass_from_faces,c.amplitudes[0][j]),'mass_face_'+str(j))))
        mechanical_face_row = sp.ImmutableMatrix(1, 5, lambda i, j:c.retained(
            self.cancel_background_expression(sp.diff(mechanical_from_faces,c.amplitudes[0][j]),'mechanical_face_'+str(j))))
        mass_join = self.cancel_background_row(mass_increment-mass_face_row,'mass_join').applyfunc(c.retained)
        mechanical_join = self.cancel_background_row(mechanical_increment-mechanical_face_row,'mechanical_join').applyfunc(c.retained)
        mechanical_sum = self.cancel_background_row(mechanical_increment+mechanical_face_row,'mechanical_sum').applyfunc(c.retained)
        faces = [{key:c.retained(value) for key, value in face.items()} for face in faces]
        return {'ACOUSTIC_WAVE_EQUATION':wave, 'ACOUSTIC_PRESSURE':pressure,
                'ACOUSTIC_VELOCITY':velocity, 'ENERGY_ANSATZ':energy_ansatz,
                'ENERGY_COEFFICIENT_SOLVE':tuple(coefficients.items()), 'ACOUSTIC_ENERGY':energy,
                'ACOUSTIC_LOCAL_CURRENT':local_flux, 'LOCAL_ENERGY_BALANCE_RESIDUAL':local_residual,
                'HARMONIC_FIELD_ANSATZ':harmonic, 'HARMONIC_WAVE_ROWS':harmonic_wave,
                'HARMONIC_ENERGY':harmonic_energy, 'HARMONIC_CURRENT':harmonic_flux,
                'NORMAL_CURRENT_DEPTH_COEFFICIENT':current_coefficient,
                'FINITE_DEPTH_INTEGRAL':finite_depth, 'FINITE_DEPTH_DERIVATIVE_RESIDUAL':finite_residual,
                'DIAGONAL_DEPTH_INTEGRAL':diagonal_depth,
                'CONVERGENT_DEPTH_FACTOR':convergent_factor, 'CONVERGENT_FINITE_DEPTH_INTEGRAL':convergent_finite,
                'DEPTH_PRIMITIVE_RESIDUAL':primitive_residual,
                'DEPTH_BOUNDARY_MODULUS':boundary_modulus, 'DEPTH_BOUNDARY_MODULUS_LIMIT':boundary_modulus_limit,
                'CONVERGENT_INFINITE_DEPTH_INTEGRAL':infinite_depth, 'DEPTH_DECAY_BINDINGS':tuple(depth_map.items()),
                'DEPTH_RATE_JOIN_RESIDUAL':depth_join, 'DEPTH_CONVERGENCE_DOMAIN':sp.Gt(-sp.re(depth_rate), 0),
                'SLAB_SURFACE_DENSITY':surface_density, 'CHEMICAL_AFFINITY_DRIVER':c.retained(mus),
                'MEMORY_KERNELS':kernels, 'FACE_RECORDS':faces,
                'SOURCE_BRANCH_JOINS':joins, 'ACOUSTIC_RADICAL_SCALE':radical_scale,
                'ACOUSTIC_RADICAL_JOIN_RESIDUAL':radical_residual,
                'REDUCED_MASS_CLOSURE_INCREMENT':mass_increment,
                'RECONSTRUCTED_MASS_FACE_ROW':mass_face_row,
                'REDUCED_MECHANICAL_CLOSURE_INCREMENT':mechanical_increment,
                'RECONSTRUCTED_MECHANICAL_FACE_ROW':mechanical_face_row,
                'REDUCED_MASS_FACE_JOIN_RESIDUAL':mass_join,
                'REDUCED_MECHANICAL_FACE_JOIN_RESIDUAL':mechanical_join,
                'REDUCED_MECHANICAL_FACE_SUM':mechanical_sum}

    def emit(self, anchoring, end, suffix):
        result = self.construct(anchoring, end)
        dims = PHYSICAL_METADATA.dimensions
        def shifted(unit, other):
            return tuple(a+b for a, b in zip(unit, other))
        energy_unit = dims.measure(result['ACOUSTIC_ENERGY'])
        for key, value in result.items():
            units = None
            if key == 'LOCAL_ENERGY_BALANCE_RESIDUAL':
                units = {():shifted(energy_unit, dims.measure(self.r.omega))}
            elif key in ('FINITE_DEPTH_DERIVATIVE_RESIDUAL', 'DEPTH_PRIMITIVE_RESIDUAL'):
                units = {():dims.zero}
            elif key == 'DEPTH_BOUNDARY_MODULUS_LIMIT':
                units = {():dims.measure(self.height)}
            elif key == 'DEPTH_RATE_JOIN_RESIDUAL':
                units = {():dims.measure(self.qlegs[0])}
            elif key == 'ACOUSTIC_RADICAL_JOIN_RESIDUAL':
                units = {():tuple(2*a for a in dims.measure(self.r.omega))}
            elif key == 'SOURCE_BRANCH_JOINS':
                units = {(i,):dims.zero for i in range(len(value))}
            elif key == 'FACE_RECORDS':
                unit = dims.measure(result['SLAB_SURFACE_DENSITY'])
                unit = shifted(unit, dims.measure(self.r.omega))
                units = {(i, key):unit for i in range(len(value)) for key in
                         ('CLOSURE_RESIDUAL', 'AMPLITUDE_EQUATION_RECONSTRUCTION_RESIDUAL')}
            elif key in ('REDUCED_MASS_FACE_JOIN_RESIDUAL', 'REDUCED_MECHANICAL_FACE_JOIN_RESIDUAL', 'REDUCED_MECHANICAL_FACE_SUM'):
                row = 3 if key == 'REDUCED_MASS_FACE_JOIN_RESIDUAL' else 4
                field_units = self.c.field_units
                index = next(j for j in range(5) if self.strong[row, j] != 0)
                row_unit = shifted(dims.measure(self.strong[row, index]), field_units[index])
                units = {(j,):tuple(a-b for a, b in zip(row_unit, field_units[j])) for j in range(5)}
            if key == 'DEPTH_CONVERGENCE_DOMAIN':
                emit('CLOSED_ACOUSTIC_'+key+'_'+suffix, value)
                # Predicate grades follow its comparison operands; its
                # output dimension is logical, not the decay-rate dimension.
                emit('METADATA_CLOSED_ACOUSTIC_'+key+'_'+suffix,
                     self.modes.numeric_metadata(value.lhs-value.rhs, lambda p:dims.zero))
            else:
                physical('CLOSED_ACOUSTIC_'+key+'_'+suffix, value, zero_dimensions=units)
        return result


class ClosedCurrentPairing:
    """Two-frequency energy/port balance on a computed constant background.

    The harmonic ansatz polarizes both the stored-energy boundary work and
    the acoustic balance. Reduced row residuals enter a separate source-work
    ansatz. The depth cutoff and interface exchange remain live operands.
    """

    def __init__(self, acoustic, anchoring='LAB_HELD', end=None):
        self.acoustic = acoustic
        self.anchoring, self.end = anchoring, end
        self.balance, self.c, self.r = acoustic.balance, acoustic.c, acoustic.r
        self.modes = acoustic.modes
        self.frequencies = sp.symbols('s11cdPairingLeftFrequency s11cdPairingRightFrequency', real=True)
        dims = PHYSICAL_METADATA.dimensions
        for frequency in self.frequencies:
            dims.known[frequency] = dims.measure(self.r.omega)
        self.residual_amplitudes = tuple(sp.symbols('s11cdPairing'+side+'RowResidual0:5')
                                         for side in ('Plus', 'Minus'))
        for group in self.residual_amplitudes:
            for i, a in enumerate(group):
                j = next(j for j in range(5) if acoustic.strong[i, j] != 0)
                dims.known[a] = tuple(u+v for u, v in zip(
                    dims.measure(acoustic.strong[i, j]), self.c.field_units[j]))

    def polarize(self, expression):
        c = self.c
        amplitudes = (*c.amplitudes[1], *c.amplitudes[0])
        origin = (0,)*len(amplitudes)
        def product(left, right):
            result = {}
            for x, a in left.items():
                for y, b in right.items():
                    degree = tuple(i+j for i, j in zip(x, y))
                    result[degree] = result.get(degree, sp.S.Zero)+a*b
            return result
        @lru_cache(maxsize=None)
        def coefficients(node):
            if not node.has(*amplitudes):
                return {origin:node}
            if node in amplitudes:
                return {tuple(int(a == node) for a in amplitudes):sp.S.One}
            if node.is_Add:
                result = {}
                for child in node.args:
                    for degree, coefficient in coefficients(child).items():
                        result[degree] = result.get(degree, sp.S.Zero)+coefficient
                return result
            if node.is_Mul or (node.is_Pow and node.exp.is_Integer and node.exp >= 0):
                factors = node.args if node.is_Mul else (node.base,)*int(node.exp)
                result = {origin:sp.S.One}
                for factor in factors:
                    result = product(result, coefficients(factor))
                return result
            raise ValueError(('nonpolynomial field amplitude in current', node.func))
        polynomial = coefficients(expression)
        matrix = sp.ImmutableMatrix(5, 5, lambda i, j:polynomial.get(
            tuple(int(n == i)+int(n == 5+j) for n in range(10)), sp.S.Zero))
        reconstructed = (sp.ImmutableMatrix(1, 5, c.amplitudes[1])*matrix*
                         sp.ImmutableMatrix(c.amplitudes[0]))[0]
        reconstructed_coefficients = coefficients(reconstructed)
        residual = sp.Add(*((polynomial.get(degree, sp.S.Zero)-reconstructed_coefficients.get(degree, sp.S.Zero))*
                            sp.Mul(*(a**p for a, p in zip(amplitudes, degree)))
                            for degree in polynomial.keys() | reconstructed_coefficients.keys()))
        return matrix, residual

    @lru_cache(maxsize=None)
    def construct(self, anchoring, end):
        if (anchoring, end) != (self.anchoring, self.end):
            raise ValueError('two-frequency source context differs from requested background')
        c, r, a, b = self.c, self.r, self.acoustic, self.balance
        slab = b.construct(anchoring, end)
        conservative = c.construct(anchoring, end)
        acoustic = a.construct(anchoring, end)
        omega_left, omega_right = self.frequencies
        kleft, kright = c.leg_momenta
        qleft, qright = a.qlegs
        epsilon = r.symbols['epsilon_shape']
        rho = r.symbols['rho_m']
        phases = tuple(sp.exp(sign*sp.I*(sum(k*x for k, x in zip(r.tangents, r.x[:2]))+
                       momentum*r.z-frequency*r.t+c.phase_coordinate))
                       for sign, momentum, frequency in
                       ((1, kright, omega_right), (-1, kleft, omega_left)))
        beat = sp.simplify(phases[0]*phases[1])
        phase_rates = tuple(sp.cancel(sp.diff(beat, x)/beat) for x in (r.t, r.z))
        zero_coordinates = {r.t:0, r.z:0, a.depth:0}

        def average(expression):
            # Shield coefficient carriers before expanding the harmonic
            # polynomial. Deep expansion of a response fraction can move
            # exponential phases into an expanded denominator.
            carriers = {}
            def shield(node):
                if node == c.phase_coordinate:
                    return node
                if not node.has(c.phase_coordinate):
                    if node.is_number:
                        return node
                    return carriers.setdefault(node, sp.Dummy())
                return node.func(*(shield(child) for child in node.args))
            harmonic = shield(expression.doit())
            selected = c.ends.phase_terms(sp.expand(harmonic), c.phase_coordinate).get(0, sp.S.Zero)
            return selected.xreplace({v:k for k, v in carriers.items()}).subs(zero_coordinates)

        def real_field(plus, minus):
            return epsilon*(plus*phases[0]+minus*phases[1])/2

        field_ansatz = tuple(real_field(plus, minus) for plus, minus in zip(*c.amplitudes))
        leg_maps = []
        for sign, frequency, momentum, bulk_momentum, amplitudes in (
                (1, omega_right, kright, qright, c.amplitudes[0]),
                (-1, omega_left, kleft, qleft, c.amplitudes[1])):
            leg_maps.append({r.omega:sign*frequency, self.modes.k:sign*momentum,
                             a.qlegs[1]:sign*bulk_momentum,
                             **{k:sign*k for k in r.tangents},
                             **dict(zip(c.amplitudes[0], amplitudes))})
        time_rates = tuple(sp.cancel(sp.diff(phase, r.t)/phase) for phase in phases)
        test_rates = {test(r.z):rate*field(r.z)
                      for rate, tests, fields in zip(time_rates, b.tests, c.fields)
                      for test, field in zip(tests, fields)}
        boundary = slab['UNCONSTRAINED_BOUNDARY_WORK']+slab['CHEMICAL_TRANSPORT_BOUNDARY']
        slab_current, slab_polarization = b.matrix(boundary.subs(test_rates, simultaneous=True).doit())
        stored_energy, stored_polarization = b.matrix(conservative['TANGENTIAL_ENERGY_REDUCTION'])
        # The supplied S11b kinetic action, with the physical thickness lift
        # taken from the independently computed two-face displacement ansatz.
        face_displacements = [sp.cancel(face['DISPLACEMENT']/c.phase)
                              for face in acoustic['FACE_RECORDS']]
        thickness_plus = face_displacements[1]-face_displacements[0]
        thickness = real_field(thickness_plus.xreplace(leg_maps[0]), thickness_plus.xreplace(leg_maps[1]))
        kinetic_ansatz = (acoustic['SLAB_SURFACE_DENSITY']*sum(sp.diff(f, r.t)**2 for f in field_ansatz[:3])+
                          r.symbols['mu_W']*sp.diff(thickness, r.t)**2)/2
        kinetic_energy, kinetic_polarization = self.polarize(average(kinetic_ansatz))
        slab_energy = stored_energy+kinetic_energy
        mus = tuple(acoustic['CHEMICAL_AFFINITY_DRIVER'].xreplace(mapping) for mapping in leg_maps)
        chemical_field = real_field(*mus)
        residual_fields = tuple(real_field(plus, minus)
                               for plus, minus in zip(*self.residual_amplitudes))
        source_ansatz = sum(sp.diff(field_ansatz[i], r.t)*residual_fields[i] for i in (0, 1, 2, 4))
        source_ansatz += chemical_field*residual_fields[3]
        source_coefficient = average(source_ansatz)
        plus_power_map = sp.ImmutableMatrix(5, 5, lambda i, j:sp.diff(
            source_coefficient, c.amplitudes[1][i], self.residual_amplitudes[0][j]))
        minus_power_map = sp.ImmutableMatrix(5, 5, lambda i, j:sp.diff(
            source_coefficient, self.residual_amplitudes[1][i], c.amplitudes[0][j]))
        algebraic, relation, joins = self.modes.analytic(a.strong)
        physical_q = algebraic.subs(self.modes.q, a.qlegs[1]/acoustic['ACOUSTIC_RADICAL_SCALE'])
        pencils = tuple(physical_q.xreplace(mapping) for mapping in leg_maps)
        source_power = plus_power_map*pencils[0]+pencils[1].T*minus_power_map
        source_reconstruction = sp.expand(source_coefficient-
            (sp.ImmutableMatrix(1, 5, c.amplitudes[1])*plus_power_map*
             sp.ImmutableMatrix(self.residual_amplitudes[0]))[0]-
            (sp.ImmutableMatrix(1, 5, self.residual_amplitudes[1])*minus_power_map*
             sp.ImmutableMatrix(c.amplitudes[0]))[0])

        phi = a.phi(r.t, *r.x[:2], r.z, a.depth)
        bulk_legs = (a.alegs[1]*phases[0]*sp.exp(sp.I*qright*a.depth),
                     a.alegs[0]*phases[1]*sp.exp(-sp.I*qleft*a.depth))
        bulk_ansatz = epsilon*sum(bulk_legs)/2
        wave_rows = tuple(sp.simplify(acoustic['ACOUSTIC_WAVE_EQUATION'].subs(phi, leg).doit()/leg)
                          for leg in bulk_legs)
        pressure_scales = tuple(sp.simplify(acoustic['ACOUSTIC_PRESSURE'].subs(phi, leg).doit()/leg)
                                for leg in bulk_legs)
        velocity_scales = tuple(sp.simplify(acoustic['ACOUSTIC_VELOCITY'][3].subs(phi, leg).doit()/leg)
                                for leg in bulk_legs)
        bulk_energy_coefficient = average(acoustic['ACOUSTIC_ENERGY'].subs(phi, bulk_ansatz))
        bulk_current_coefficients = acoustic['ACOUSTIC_LOCAL_CURRENT'].applyfunc(
            lambda value:average(value.subs(phi, bulk_ansatz)))
        faces = []
        face_leg_objects = []
        for face in acoustic['FACE_RECORDS']:
            legs = tuple({key:value.xreplace(mapping) for key, value in face.items()
                          if key != 'DISPLACEMENT'} for mapping in leg_maps)
            lifts = {a.alegs[1]:legs[0]['AMPLITUDE'], a.alegs[0]:legs[1]['AMPLITUDE']}
            energy_matrix, energy_polarization = self.polarize(bulk_energy_coefficient.xreplace(lifts))
            normal_matrix, normal_polarization = self.polarize(bulk_current_coefficients[2].xreplace(lifts))
            depth_matrix, depth_polarization = self.polarize(bulk_current_coefficients[3].xreplace(lifts))
            pressure = real_field(*(leg['PRESSURE'] for leg in legs))
            outward_velocity = real_field(*(leg['OUTWARD_VELOCITY'] for leg in legs))
            mass_flux = real_field(*(leg['RELATIVE_MASS_FLUX'] for leg in legs))
            affinity = real_field(*(leg['AFFINITY'] for leg in legs))
            mechanical_response = real_field(*(acoustic['MEMORY_KERNELS']['X'].xreplace(mapping)*leg['AFFINITY']
                                               for mapping, leg in zip(leg_maps, legs)))
            velocity_legs = tuple(sp.cancel(sp.diff(leg, a.depth)/leg)*closed['AMPLITUDE']
                                  for leg, closed in zip(bulk_legs, legs))
            face_leg_objects.append(tuple({
                **{key:leg[key] for key in ('AMPLITUDE', 'PRESSURE', 'OUTWARD_VELOCITY',
                                            'RELATIVE_MASS_FLUX', 'AFFINITY')},
                'CHEMICAL_POTENTIAL':mu, 'BULK_VELOCITY':velocity,
                'MECHANICAL_RESPONSE':acoustic['MEMORY_KERNELS']['X'].xreplace(mapping)*leg['AFFINITY'],
                'MASS_RESPONSE_A':acoustic['MEMORY_KERNELS']['A'].xreplace(mapping),
                'MASS_RESPONSE_V':acoustic['MEMORY_KERNELS']['V'].xreplace(mapping)}
                for leg, mu, velocity, mapping in zip(legs, mus, velocity_legs, leg_maps)))
            bulk_velocity = real_field(*velocity_legs)
            # The two supplied port-work forms are separate ansatz operands.
            port_ansatz = (pressure+mechanical_response)*outward_velocity+chemical_field*mass_flux
            split_ansatz = pressure*bulk_velocity+affinity*mass_flux+mechanical_response*outward_velocity
            interface_ansatz = affinity*mass_flux+mechanical_response*outward_velocity
            port_matrix, port_polarization = self.polarize(average(port_ansatz))
            split_matrix, split_polarization = self.polarize(average(split_ansatz))
            interface_matrix, interface_polarization = self.polarize(average(interface_ansatz))
            faces.append({'ENERGY_DENSITY_MATRIX':energy_matrix, 'NORMAL_CURRENT_DENSITY_MATRIX':normal_matrix,
                          'DEPTH_CURRENT_MATRIX':depth_matrix, 'PORT_POWER_MATRIX':port_matrix,
                          'SPLIT_PORT_POWER_MATRIX':split_matrix, 'INTERFACE_POWER_MATRIX':interface_matrix,
                          'PORT_IDENTITY_RESIDUAL':port_matrix-split_matrix,
                          'FACE_CURRENT_JOIN_RESIDUAL':port_matrix-interface_matrix-depth_matrix,
                          'POLARIZATION_RESIDUALS':(energy_polarization, normal_polarization, depth_polarization,
                                                   port_polarization, split_polarization, interface_polarization)})
        def sum_face(key):
            return sum((face[key] for face in faces), sp.zeros(5))
        depth_factor = sp.exp(sp.I*(qright-qleft)*a.depth)
        depth_rate = sp.cancel(sp.diff(depth_factor, a.depth)/depth_factor)
        primitive = depth_factor/depth_rate
        depth_integral = primitive.subs(a.depth, a.height)-primitive.subs(a.depth, 0)
        equal_depth_integral = sp.integrate(depth_factor.subs(qleft, qright), (a.depth, 0, a.height))
        bulk_energy = sum_face('ENERGY_DENSITY_MATRIX')
        bulk_current = sum_face('NORMAL_CURRENT_DENSITY_MATRIX')
        bulk_depth_current = sum_face('DEPTH_CURRENT_MATRIX')
        port_power = sum_face('PORT_POWER_MATRIX')
        interface_power = sum_face('INTERFACE_POWER_MATRIX')
        total_energy = slab_energy+depth_integral*bulk_energy
        total_current = slab_current+depth_integral*bulk_current
        top_power = depth_factor.subs(a.depth, a.height)*bulk_depth_current
        slab_residual = (phase_rates[0]*slab_energy+phase_rates[1]*slab_current+port_power-source_power)
        bulk_density_residual = (phase_rates[0]*bulk_energy+phase_rates[1]*bulk_current+depth_rate*bulk_depth_current)
        finite_residual = phase_rates[0]*total_energy+phase_rates[1]*total_current+interface_power+top_power-source_power
        return {'FREQUENCY_LEGS':self.frequencies, 'NORMAL_LEGS':c.leg_momenta, 'BULK_LEGS':a.qlegs,
                'BEAT_PHASE':beat, 'BEAT_RATES':phase_rates, 'SOURCE_BRANCH_JOINS':joins,
                'ACOUSTIC_WAVE_ROWS':wave_rows, 'KINETIC_ACTION_ANSATZ':kinetic_ansatz,
                'STORED_ENERGY_MATRIX':stored_energy, 'KINETIC_ENERGY_MATRIX':kinetic_energy,
                'SLAB_ENERGY_MATRIX':slab_energy, 'SLAB_CURRENT_MATRIX':slab_current,
                'SLAB_POLARIZATION_RESIDUALS':(slab_polarization, stored_polarization, kinetic_polarization),
                'PLUS_ROW_POWER_MAP':plus_power_map, 'MINUS_ROW_POWER_MAP':minus_power_map,
                'ROW_POWER_RECONSTRUCTION_RESIDUAL':source_reconstruction, 'CLOSED_PENCIL_LEGS':pencils,
                'SOURCE_POWER_MATRIX':source_power, 'FACE_RECORDS':faces,
                'OPEN_BULK_ENERGY_COEFFICIENT':bulk_energy_coefficient,
                'OPEN_BULK_CURRENT_COEFFICIENTS':bulk_current_coefficients,
                'OPEN_BULK_PRESSURE_SCALES':pressure_scales,
                'OPEN_BULK_VELOCITY_SCALES':velocity_scales,
                'FACE_LEG_OBJECTS':face_leg_objects,
                'DEPTH_PHASE_FACTOR':depth_factor, 'DEPTH_PHASE_RATE':depth_rate,
                'DEPTH_PRIMITIVE':primitive,
                'BULK_ENERGY_DENSITY_MATRIX':bulk_energy, 'BULK_NORMAL_CURRENT_DENSITY_MATRIX':bulk_current,
                'BULK_DEPTH_CURRENT_MATRIX':bulk_depth_current, 'PORT_POWER_MATRIX':port_power,
                'INTERFACE_POWER_MATRIX':interface_power, 'GENERIC_DEPTH_INTEGRAL':depth_integral,
                'GENERIC_DEPTH_DOMAIN':sp.Ne(qleft, qright), 'EQUAL_DEPTH_INTEGRAL':equal_depth_integral,
                'EQUAL_DEPTH_DOMAIN':sp.Eq(qleft, qright),
                'FINITE_TOTAL_ENERGY_MATRIX':total_energy, 'FINITE_TOTAL_CURRENT_MATRIX':total_current,
                'TOP_BOUNDARY_POWER_MATRIX':top_power, 'SLAB_BALANCE_RESIDUAL':slab_residual,
                'BULK_DENSITY_BALANCE_RESIDUAL':bulk_density_residual,
                'FINITE_BALANCE_RESIDUAL':finite_residual}


    @staticmethod
    def rational_coefficient(expression):
        """Exact coefficient-field arithmetic after explicit instance binding."""
        from sympy.polys.fields import field
        generators = sorted(expression.free_symbols, key=sp.default_sort_key)
        if not generators:
            return sp.cancel(expression)
        domain = field(generators, sp.QQ_I)[0]
        return domain.from_expr(expression).as_expr()

    def wave_curve_reduction(self, expression, wave_rows):
        """Rational remainder on both wave curves, with denominator retained."""
        carriers = {v:sp.Dummy() for v in sorted(expression.atoms(sp.exp), key=sp.default_sort_key)}
        restore = {v:k for k, v in carriers.items()}
        numerator, denominator = sp.fraction(self.rational_coefficient(expression.xreplace(carriers)))
        right_quotient, right_remainder = sp.div(numerator, wave_rows[0], self.acoustic.qlegs[1])
        left_quotient, remainder = sp.div(right_remainder, wave_rows[1], self.acoustic.qlegs[0])
        reconstruction = sp.expand(numerator-right_quotient*wave_rows[0]-left_quotient*wave_rows[1]-remainder)
        return {'RATIONAL_REMAINDER':self.rational_coefficient(remainder/denominator).xreplace(restore),
                'NUMERATOR':numerator.xreplace(restore), 'DENOMINATOR':denominator.xreplace(restore),
                'WAVE_QUOTIENTS':(right_quotient.xreplace(restore), left_quotient.xreplace(restore)),
                'DIVISION_RECONSTRUCTION_RESIDUAL':reconstruction.xreplace(restore)}

    @staticmethod
    def carrier_expansion(expression):
        """Distribute polynomials while retaining inverse/phase carriers.

        A zero is an exact identity in these actual operands. A nonzero is
        returned unchanged in meaning; no independence or pole extension is
        inferred from treating an operand as a temporary polynomial carrier.
        """
        carriers = {}
        def shield(node):
            if node.func == sp.exp or (node.is_Pow and
                    not (node.exp.is_Integer and node.exp >= 0)):
                return carriers.setdefault(node, sp.Dummy())
            if not node.args:
                return node
            return node.func(*(shield(child) for child in node.args))
        expanded = sp.expand(shield(expression))
        return expanded.xreplace({v:k for k, v in carriers.items()})

    def rational_reconstruction(self, expression, bindings):
        """Clear actual inverse carriers without a multivariate GCD.

        The original denominator is retained with its restored units and
        domain. Binding precedes only the polynomial numerator expansion.
        """
        carriers = tuple(sorted({p for p in expression.atoms(sp.Pow)
            if p.exp.is_Integer and p.exp.is_negative}, key=sp.default_sort_key))
        terms = dict(polynomial_terms(expression, carriers))
        powers = tuple(max(degree[i] for degree in terms) for i in range(len(carriers)))
        denominator = sp.Mul(*(p**(-n) for p, n in zip(carriers, powers)))
        numerator = sp.Add(*(coefficient*sp.Mul(*(p**(-n+d) for p, n, d in zip(carriers, powers, degree)))
                             for degree, coefficient in terms.items()))
        remainder = self.carrier_expansion(numerator.xreplace(bindings))
        return remainder/denominator.xreplace(bindings), denominator, remainder

    def split_balance_checks(self, result, bindings, progress=lambda record:None):
        """Local wave certificate, actual face lifts, then finite-depth sum.

        The acoustic division precedes the response substitution. Large
        closed matrices are compared by polynomial reconstruction, retaining
        their denominator carriers. Only the separate face-row joins and slab
        identity use coefficient-field cancellation at the material binding.
        """
        expand = self.carrier_expansion
        scalar = self.rational_coefficient
        c, a = self.c, self.acoustic
        qleft, qright = a.qlegs
        time_rate, normal_rate = result['BEAT_RATES']
        depth_rate = result['DEPTH_PHASE_RATE']
        waves = result['ACOUSTIC_WAVE_ROWS']
        potential_pair = a.alegs[0]*a.alegs[1]
        open_values = (result['OPEN_BULK_ENERGY_COEFFICIENT'],
                       result['OPEN_BULK_CURRENT_COEFFICIENTS'][2],
                       result['OPEN_BULK_CURRENT_COEFFICIENTS'][3])
        pair_weight = sp.cancel(sp.diff(result['PLUS_ROW_POWER_MAP'][0, 0], self.frequencies[0])/
                                sp.diff(time_rate, self.frequencies[0]))
        coefficients = tuple(dict(polynomial_terms(value, a.alegs)).get((1, 1), sp.S.Zero)
                             for value in open_values)
        checks = {'OPEN_BULK_POLARIZATION_RESIDUAL':tuple(expand(value-coefficient*potential_pair)
                    for value, coefficient in zip(open_values, coefficients)),
                  'OPEN_BULK_PAIR_COEFFICIENTS':coefficients, 'PORT_BILINEAR_WEIGHT':pair_weight}
        open_balance = sp.expand(time_rate*open_values[0]+normal_rate*open_values[1]+depth_rate*open_values[2])
        certificate = self.wave_curve_reduction(open_balance, waves)
        checks.update({'OPEN_BULK_BALANCE_OPERAND':open_balance,
                       'OPEN_BULK_WAVE_QUOTIENTS':certificate['WAVE_QUOTIENTS'],
                       'OPEN_BULK_WAVE_DENOMINATOR':certificate['DENOMINATOR'],
                       'OPEN_BULK_BALANCE_RESIDUAL':certificate['RATIONAL_REMAINDER'],
                       'OPEN_BULK_DIVISION_RESIDUAL':certificate['DIVISION_RECONSTRUCTION_RESIDUAL']})
        # At equal q the two wave rows share a variable. Keep the second
        # row's eliminant; a generic two-independent-q division is insufficient.
        equal_waves = tuple(wave.subs(qleft, qright) for wave in waves)
        eliminant = sp.rem(equal_waves[1], equal_waves[0], qright)
        equal_numerator, equal_denominator = sp.fraction(scalar(open_balance.subs(qleft, qright)))
        quotient_q, remainder_q = sp.div(equal_numerator, equal_waves[0], qright)
        quotient_omega, remainder_equal = sp.div(remainder_q, eliminant, self.frequencies[0])
        checks.update({'EQUAL_DEPTH_WAVE_ROWS':equal_waves, 'EQUAL_DEPTH_WAVE_ELIMINANT':eliminant,
                       'EQUAL_OPEN_WAVE_QUOTIENTS':(quotient_q, quotient_omega),
                       'EQUAL_OPEN_WAVE_DENOMINATOR':equal_denominator,
                       'EQUAL_OPEN_BALANCE_RESIDUAL':scalar(remainder_equal/equal_denominator),
                       'EQUAL_OPEN_DIVISION_RESIDUAL':sp.expand(equal_numerator-
                           quotient_q*equal_waves[0]-quotient_omega*eliminant-remainder_equal)})
        progress({'stage':'open_acoustic_divisions'})
        lift_rows, lift_residuals, reconstruction_residuals, leg_join_residuals = [], [], [], []
        port_residuals, join_residuals, denominators = [], [], []
        port_linear_residuals, port_reconstruction_residuals, depth_coefficient_residuals = [], [], []
        all_port_rows = []
        for face_number, (face, legs) in enumerate(zip(result['FACE_RECORDS'], result['FACE_LEG_OBJECTS'])):
            rows = []
            for leg_number, (leg, amplitudes) in enumerate(zip(legs, c.amplitudes)):
                row = sp.ImmutableMatrix(1, 5, [dict(polynomial_terms(leg['AMPLITUDE'], amplitudes)).get(
                    tuple(int(i == j) for i in range(5)), sp.S.Zero) for j in range(5)])
                rows.append(row)
                lift_residuals.append(expand(leg['AMPLITUDE']-(row*sp.ImmutableMatrix(amplitudes))[0]))
                rho = self.r.symbols['rho_m']
                identities = (
                    leg['PRESSURE']-result['OPEN_BULK_PRESSURE_SCALES'][leg_number]*leg['AMPLITUDE'],
                    leg['BULK_VELOCITY']-result['OPEN_BULK_VELOCITY_SCALES'][leg_number]*leg['AMPLITUDE'],
                    leg['AFFINITY']-leg['CHEMICAL_POTENTIAL']+leg['PRESSURE']/rho,
                    leg['RELATIVE_MASS_FLUX']-rho*(leg['BULK_VELOCITY']-leg['OUTWARD_VELOCITY']),
                    leg['RELATIVE_MASS_FLUX']-leg['MASS_RESPONSE_A']*leg['AFFINITY']-
                        leg['MASS_RESPONSE_V']*leg['OUTWARD_VELOCITY'])
                # Each identity tests the whole five-field linear row.
                for identity in identities:
                    terms = dict(polynomial_terms(identity, amplitudes))
                    leg_join_residuals.append(sp.ImmutableMatrix(1, 5, [scalar(terms.get(
                        tuple(int(i == j) for i in range(5)), sp.S.Zero).xreplace(bindings)) for j in range(5)]))
                denominators.append(tuple(sorted({power.base for value in leg.values()
                    for power in value.atoms(sp.Pow) if power.exp.is_negative}, key=sp.default_sort_key)))
            lift_rows.append(tuple(rows))
            outer = rows[1].T*rows[0]
            for key, coefficient in zip(('ENERGY_DENSITY_MATRIX', 'NORMAL_CURRENT_DENSITY_MATRIX',
                                         'DEPTH_CURRENT_MATRIX'), coefficients):
                reconstruction_residuals.append((face[key]-coefficient*outer).applyfunc(expand))
            # Reconstruct the three supplied port polynomials using the full
            # computed linear rows. Constitutive row joins are reduced before
            # products; no two-frequency response denominator is expanded.
            port_rows = []
            for leg, amplitudes in zip(legs, c.amplitudes):
                linear_rows = {}
                for name in ('PRESSURE', 'OUTWARD_VELOCITY', 'RELATIVE_MASS_FLUX', 'AFFINITY',
                             'CHEMICAL_POTENTIAL', 'BULK_VELOCITY', 'MECHANICAL_RESPONSE'):
                    terms = dict(polynomial_terms(leg[name], amplitudes))
                    row = sp.ImmutableMatrix(1, 5, [terms.get(tuple(int(i == j) for i in range(5)),
                        sp.S.Zero) for j in range(5)])
                    port_linear_residuals.append(expand(leg[name]-(row*sp.ImmutableMatrix(amplitudes))[0]))
                    linear_rows[name] = row
                port_rows.append(linear_rows)
            all_port_rows.append(tuple(port_rows))
            def bilinear(left, right):
                return pair_weight*(left[1].T*right[0]+right[1].T*left[0])
            column = lambda key:tuple(leg[key] for leg in port_rows)
            p, v, velocity, mass, affinity, mu, response = map(column, (
                'PRESSURE', 'BULK_VELOCITY', 'OUTWARD_VELOCITY', 'RELATIVE_MASS_FLUX',
                'AFFINITY', 'CHEMICAL_POTENTIAL', 'MECHANICAL_RESPONSE'))
            reconstructed_interface = bilinear(affinity, mass)+bilinear(response, velocity)
            reconstructed_port = bilinear(p, velocity)+bilinear(response, velocity)+bilinear(mu, mass)
            reconstructed_split = bilinear(p, v)+reconstructed_interface
            port_reconstruction, split_reconstruction, interface_reconstruction = tuple(
                (face[key]-operand).applyfunc(expand) for key, operand in zip(
                    ('PORT_POWER_MATRIX', 'SPLIT_PORT_POWER_MATRIX', 'INTERFACE_POWER_MATRIX'),
                    (reconstructed_port, reconstructed_split, reconstructed_interface)))
            port_reconstruction_residuals.extend((port_reconstruction, split_reconstruction, interface_reconstruction))
            row_checks = leg_join_residuals[-10:]
            e_pressure = (row_checks[0], row_checks[5])
            e_velocity = (row_checks[1], row_checks[6])
            e_affinity = (row_checks[2], row_checks[7])
            e_mass = (row_checks[3], row_checks[8])
            bind_pair = lambda pair:tuple(row.xreplace(bindings) for row in pair)
            port_defect = bilinear(bind_pair(p), e_mass)/rho-bilinear(e_affinity, bind_pair(mass))
            pressure_scales = result['OPEN_BULK_PRESSURE_SCALES']
            velocity_scales = result['OPEN_BULK_VELOCITY_SCALES']
            pressure_from_amplitude = tuple(scale*row for scale, row in zip(pressure_scales, rows))
            depth_coefficient = pair_weight*sum(pressure_scales[i]*velocity_scales[1-i] for i in range(2))
            depth_coefficient_join = expand(depth_coefficient-coefficients[2])
            depth_coefficient_residuals.append(depth_coefficient_join)
            depth_defect = (bilinear(e_pressure, bind_pair(v))+
                           bilinear(bind_pair(pressure_from_amplitude), e_velocity)+
                           depth_coefficient_join.xreplace(bindings)*outer.xreplace(bindings))
            # These residuals compose independent polynomial reconstructions,
            # the five-field constitutive row residuals, and the open current.
            port_residuals.append((port_reconstruction-split_reconstruction).xreplace(bindings)+port_defect)
            join_residuals.append((port_reconstruction-interface_reconstruction-
                reconstruction_residuals[-1]).xreplace(bindings)+port_defect+depth_defect)
            progress({'stage':'closed_face', 'face':face_number})
        checks.update({'FACE_AMPLITUDE_ROWS':tuple(lift_rows), 'FACE_DENOMINATOR_FACTORS':tuple(denominators),
                       'FACE_PORT_ROWS':tuple(all_port_rows),
                       'FACE_LINEAR_LIFT_RESIDUAL':tuple(lift_residuals),
                       'FACE_BULK_RECONSTRUCTION_RESIDUAL':tuple(reconstruction_residuals),
                       'FACE_LEG_JOIN_RESIDUAL':tuple(leg_join_residuals),
                       'FACE_PORT_LINEAR_LIFT_RESIDUAL':tuple(port_linear_residuals),
                       'FACE_PORT_RECONSTRUCTION_RESIDUAL':tuple(port_reconstruction_residuals),
                       'FACE_OPEN_DEPTH_COEFFICIENT_RESIDUAL':tuple(depth_coefficient_residuals),
                       'FACE_PORT_COMPOSED_IDENTITY_RESIDUAL':tuple(port_residuals),
                       'FACE_DEPTH_COMPOSED_JOIN_RESIDUAL':tuple(join_residuals)})
        depth_integral = result['GENERIC_DEPTH_INTEGRAL']
        top_factor = result['DEPTH_PHASE_FACTOR'].subs(a.depth, a.height)
        boundary_residual = sp.simplify(top_factor-1-depth_rate*depth_integral)
        primitive_residual = sp.simplify(sp.diff(result['DEPTH_PRIMITIVE'], a.depth)-result['DEPTH_PHASE_FACTOR'])
        face_join = sum((face['FACE_CURRENT_JOIN_RESIDUAL'] for face in result['FACE_RECORDS']), sp.zeros(5))
        composition_operand = (result['FINITE_BALANCE_RESIDUAL']-result['SLAB_BALANCE_RESIDUAL']-
            depth_integral*result['BULK_DENSITY_BALANCE_RESIDUAL']+face_join-
            boundary_residual*result['BULK_DEPTH_CURRENT_MATRIX'])
        composition_certificates = [self.rational_reconstruction(value, bindings) for value in composition_operand]
        composition = sp.ImmutableMatrix(5, 5, [record[0] for record in composition_certificates])
        composition_denominators = sp.ImmutableMatrix(5, 5, [record[1] for record in composition_certificates])
        composition_numerators = sp.ImmutableMatrix(5, 5, [record[2] for record in composition_certificates])
        equal_energy = (result['SLAB_ENERGY_MATRIX']+result['EQUAL_DEPTH_INTEGRAL']*
                        result['BULK_ENERGY_DENSITY_MATRIX']).subs(qleft, qright)
        equal_current = (result['SLAB_CURRENT_MATRIX']+result['EQUAL_DEPTH_INTEGRAL']*
                         result['BULK_NORMAL_CURRENT_DENSITY_MATRIX']).subs(qleft, qright)
        equal_raw = (time_rate*equal_energy+normal_rate*equal_current+
            (result['INTERFACE_POWER_MATRIX']+result['BULK_DEPTH_CURRENT_MATRIX']-
             result['SOURCE_POWER_MATRIX']).subs(qleft, qright))
        equal_composition = (equal_raw-result['SLAB_BALANCE_RESIDUAL'].subs(qleft, qright)-
            result['EQUAL_DEPTH_INTEGRAL']*result['BULK_DENSITY_BALANCE_RESIDUAL'].subs(qleft, qright)+
            face_join.subs(qleft, qright)).xreplace(bindings).applyfunc(expand)
        checks.update({'DEPTH_PRIMITIVE_RESIDUAL':primitive_residual, 'DEPTH_BOUNDARY_RESIDUAL':boundary_residual,
                       'FINITE_COMPOSITION_RESIDUAL':composition,
                       'FINITE_COMPOSITION_DENOMINATOR_MATRIX':composition_denominators,
                       'FINITE_COMPOSITION_NUMERATOR_RESIDUAL':composition_numerators,
                       'EQUAL_DEPTH_TOTAL_ENERGY_MATRIX':equal_energy,
                       'EQUAL_DEPTH_TOTAL_CURRENT_MATRIX':equal_current,
                       'EQUAL_DEPTH_COMPOSITION_RESIDUAL':equal_composition})
        progress({'stage':'finite_compositions'})
        checks['SLAB_BALANCE_RESIDUAL'] = result['SLAB_BALANCE_RESIDUAL'].xreplace(bindings).applyfunc(scalar)
        open_remainder_coefficient = dict(polynomial_terms(certificate['RATIONAL_REMAINDER'], a.alegs)).get((1, 1), sp.S.Zero)
        equal_remainder_coefficient = dict(polynomial_terms(checks['EQUAL_OPEN_BALANCE_RESIDUAL'], a.alegs)).get((1, 1), sp.S.Zero)
        outer_sum = sum((rows[1].T*rows[0] for rows in lift_rows), sp.zeros(5))
        reconstruction = sum((time_rate*reconstruction_residuals[i]+
            normal_rate*reconstruction_residuals[i+1]+depth_rate*reconstruction_residuals[i+2]
            for i in range(0, len(reconstruction_residuals), 3)), sp.zeros(5))
        bulk_remainder = reconstruction+open_remainder_coefficient*outer_sum
        equal_bulk_remainder = reconstruction.subs(qleft, qright)+equal_remainder_coefficient*outer_sum.subs(qleft, qright)
        closed_join_remainder = sum(join_residuals, sp.zeros(5))
        checks['BULK_COMPOSED_BALANCE_RESIDUAL'] = bulk_remainder.xreplace(bindings)
        checks['FINITE_COMPOSED_BALANCE_RESIDUAL'] = (checks['SLAB_BALANCE_RESIDUAL']+
            (depth_integral*bulk_remainder+composition+boundary_residual*result['BULK_DEPTH_CURRENT_MATRIX']).xreplace(bindings)-
            closed_join_remainder)
        checks['EQUAL_DEPTH_COMPOSED_BALANCE_RESIDUAL'] = (checks['SLAB_BALANCE_RESIDUAL'].subs(qleft, qright)+
            (result['EQUAL_DEPTH_INTEGRAL']*equal_bulk_remainder+equal_composition).xreplace(bindings)-
            closed_join_remainder.subs(qleft, qright))
        progress({'stage':'slab_balance'})
        return checks

    def split_derivative_checks(self, result, checks, bindings, progress=lambda record:None):
        """Differentiate the actual operands on a regular acoustic sheet.

        The differentiated local wave certificate and the product-rule
        reconstructions accompany the finite balance. This retains both row
        power maps, interface exchange, the top boundary and kernel dispersion.
        """
        qleft, qright = self.acoustic.qlegs
        waves = result['ACOUSTIC_WAVE_ROWS']
        time_rate, normal_rate = result['BEAT_RATES']
        values = {key:value.xreplace(bindings) for key, value in result.items()
                  if isinstance(value, sp.MatrixBase)}
        pencils = tuple(value.xreplace(bindings) for value in result['CLOSED_PENCIL_LEGS'])
        plus, minus = values['PLUS_ROW_POWER_MAP'], values['MINUS_ROW_POWER_MAP']
        depth_integral = result['GENERIC_DEPTH_INTEGRAL']
        depth_jet = sp.limit(sp.diff(depth_integral, qright), qright, qleft)
        direct_depth_jet = sp.integrate(sp.diff(result['DEPTH_PHASE_FACTOR'], qright).subs(qright, qleft),
                                       (self.acoustic.depth, 0, self.acoustic.height))
        records = {'EQUAL_DEPTH_RIGHT_PARTIAL_INTEGRAL':depth_jet,
                   'EQUAL_DEPTH_RIGHT_PARTIAL_JOIN_RESIDUAL':sp.simplify(depth_jet-direct_depth_jet)}
        for name, variable in (('NORMAL', self.c.leg_momenta[1]), ('FREQUENCY', self.frequencies[1])):
            prefix = name+'_'
            denominator = sp.diff(waves[0], qright)
            transport = sp.cancel(-sp.diff(waves[0], variable)/denominator)
            def derivative(value):
                return value.diff(variable)+transport.xreplace(bindings)*value.diff(qright)
            def local_derivative_of(value):
                return value.diff(variable)+transport*value.diff(qright)
            operands = tuple(derivative(value) for value in (
                time_rate*values['FINITE_TOTAL_ENERGY_MATRIX'],
                normal_rate*values['FINITE_TOTAL_CURRENT_MATRIX'],
                values['INTERFACE_POWER_MATRIX'], values['TOP_BOUNDARY_POWER_MATRIX']))
            source = derivative(values['SOURCE_POWER_MATRIX'])
            pencil_derivative = derivative(pencils[0])
            source_terms = (derivative(plus)*pencils[0], plus*pencil_derivative,
                            derivative(pencils[1]).T*minus, pencils[1].T*derivative(minus))
            source_reconstruction = (source-sum(source_terms, sp.zeros(5))).applyfunc(self.carrier_expansion)
            # The local operand is differentiated before its wave remainder is
            # taken; this is not a derivative of an already reduced zero.
            local_derivative = local_derivative_of(checks['OPEN_BULK_BALANCE_OPERAND'])
            local_certificate = self.wave_curve_reduction(local_derivative, waves)
            records.update({prefix+'RADICAL_TRANSPORT':transport,
                            prefix+'RADICAL_DERIVATIVE_DENOMINATOR':denominator,
                            prefix+'RADICAL_TANGENCY_RESIDUAL':sp.cancel(local_derivative_of(waves[0])),
                            prefix+'PHASE_DERIVATIVE_WEIGHTS':tuple(derivative(rate) for rate in (time_rate, normal_rate)),
                            prefix+'BALANCE_DERIVATIVE_OPERANDS':operands,
                            prefix+'SOURCE_POWER_DERIVATIVE_MATRIX':source,
                            prefix+'SOURCE_PRODUCT_DERIVATIVE_TERMS':source_terms,
                            prefix+'SOURCE_PRODUCT_DERIVATIVE_RESIDUAL':source_reconstruction,
                            prefix+'CLOSED_PENCIL_DERIVATIVE_MATRIX':pencil_derivative,
                            prefix+'OPEN_BULK_TANGENT_OPERAND':local_derivative,
                            prefix+'OPEN_BULK_TANGENT_WAVE_QUOTIENTS':local_certificate['WAVE_QUOTIENTS'],
                            prefix+'OPEN_BULK_TANGENT_DENOMINATOR':local_certificate['DENOMINATOR'],
                            prefix+'OPEN_BULK_TANGENT_RESIDUAL':local_certificate['RATIONAL_REMAINDER'],
                            prefix+'OPEN_BULK_TANGENT_DIVISION_RESIDUAL':local_certificate['DIVISION_RECONSTRUCTION_RESIDUAL'],
                            prefix+'DEPTH_BOUNDARY_DERIVATIVE_RESIDUAL':sp.simplify(local_derivative_of(
                                result['DEPTH_PHASE_FACTOR'].subs(self.acoustic.depth, self.acoustic.height)-1)-
                                local_derivative_of(result['DEPTH_PHASE_RATE']*depth_integral))})
            progress({'stage':'regular_sheet_derivative', 'variable':name})
        return records

    def check_output_units(self, key, value, result):
        dims = PHYSICAL_METADATA.dimensions
        energy = dims.measure(self.c.construct(self.anchoring, self.end)['TANGENTIAL_ENERGY_REDUCTION'])
        length, frequency = dims.measure(self.r.z), dims.measure(self.r.omega)
        add = lambda u, v:tuple(a+b for a, b in zip(u, v))
        sub = lambda u, v:tuple(a-b for a, b in zip(u, v))
        power = add(energy, frequency)
        bulk_power = sub(power, length)
        potential = dims.measure(self.acoustic.alegs[0])
        field_units = self.c.field_units
        matrix_units = lambda unit:{(5*i+j,):sub(sub(unit, field_units[i]), field_units[j])
                                   for i in range(5) for j in range(5)}
        if key.startswith(('NORMAL_', 'FREQUENCY_')):
            prefix, name = key.split('_', 1)
            variable_unit = dims.measure(self.c.leg_momenta[1] if prefix == 'NORMAL' else self.frequencies[1])
            derivative_power = sub(power, variable_unit)
            if name in ('SOURCE_POWER_DERIVATIVE_MATRIX', 'SOURCE_PRODUCT_DERIVATIVE_RESIDUAL'):
                return matrix_units(derivative_power)
            if name in ('BALANCE_DERIVATIVE_OPERANDS', 'SOURCE_PRODUCT_DERIVATIVE_TERMS'):
                return {(i, *path):unit for i in range(len(value)) for path, unit in matrix_units(derivative_power).items()}
            if name == 'CLOSED_PENCIL_DERIVATIVE_MATRIX':
                row_units = [dims.measure(a) for a in self.residual_amplitudes[0]]
                return {(5*i+j,):sub(sub(row_units[i], field_units[j]), variable_unit)
                        for i in range(5) for j in range(5)}
            if name == 'PHASE_DERIVATIVE_WEIGHTS':
                return {(0,):sub(frequency, variable_unit), (1,):sub(tuple(-v for v in length), variable_unit)}
            if name in ('OPEN_BULK_TANGENT_RESIDUAL', 'OPEN_BULK_TANGENT_OPERAND'):
                return {():sub(bulk_power, variable_unit)}
            if name == 'OPEN_BULK_TANGENT_WAVE_QUOTIENTS':
                denominator_unit = dims.measure(result[prefix+'_OPEN_BULK_TANGENT_DENOMINATOR'])
                quotient_unit = sub(add(sub(bulk_power, variable_unit), denominator_unit), add(frequency, frequency))
                return {(i,):quotient_unit for i in range(len(value))}
            if name == 'OPEN_BULK_TANGENT_DIVISION_RESIDUAL':
                # Bound coefficients still represent the restored units of
                # the physical wave division, including its denominator.
                denominator_unit = dims.measure(result[prefix+'_OPEN_BULK_TANGENT_DENOMINATOR'])
                return {():add(sub(bulk_power, variable_unit), denominator_unit)}
            if name == 'RADICAL_TANGENCY_RESIDUAL':
                return {():sub(add(frequency, frequency), variable_unit)}
            if name == 'DEPTH_BOUNDARY_DERIVATIVE_RESIDUAL':
                return {():tuple(-u for u in variable_unit)}
        if key in ('EQUAL_DEPTH_RIGHT_PARTIAL_INTEGRAL', 'EQUAL_DEPTH_RIGHT_PARTIAL_JOIN_RESIDUAL'):
            return {():add(length, length)}
        if key in ('FINITE_COMPOSITION_RESIDUAL', 'EQUAL_DEPTH_COMPOSITION_RESIDUAL',
                   'FINITE_COMPOSED_BALANCE_RESIDUAL', 'EQUAL_DEPTH_COMPOSED_BALANCE_RESIDUAL'):
            return matrix_units(power)
        if key == 'FINITE_COMPOSITION_NUMERATOR_RESIDUAL':
            return {path:add(unit, dims.measure(result['FINITE_COMPOSITION_DENOMINATOR_MATRIX'][path[0]]))
                    for path, unit in matrix_units(power).items()}
        if key == 'BULK_COMPOSED_BALANCE_RESIDUAL':
            return matrix_units(bulk_power)
        if key in ('OPEN_BULK_BALANCE_RESIDUAL', 'EQUAL_OPEN_BALANCE_RESIDUAL'):
            return {():bulk_power}
        if key in ('OPEN_BULK_DIVISION_RESIDUAL', 'EQUAL_OPEN_DIVISION_RESIDUAL'):
            name = 'OPEN_BULK_WAVE_DENOMINATOR' if key.startswith('OPEN_BULK') else 'EQUAL_OPEN_WAVE_DENOMINATOR'
            return {():add(bulk_power, dims.measure(result[name]))}
        if key in ('DEPTH_PRIMITIVE_RESIDUAL', 'DEPTH_BOUNDARY_RESIDUAL'):
            return {():dims.zero}
        if key == 'OPEN_BULK_POLARIZATION_RESIDUAL':
            return {(i,):unit for i, unit in enumerate((sub(energy, length), power, power))}
        if key == 'FACE_LINEAR_LIFT_RESIDUAL':
            return {(i,):potential for i in range(len(value))}
        if key == 'FACE_AMPLITUDE_ROWS':
            return {(i, leg, j):sub(potential, field_units[j]) for i in range(len(value))
                    for leg in range(2) for j in range(5)}
        if key == 'FACE_BULK_RECONSTRUCTION_RESIDUAL':
            return {(i, *path):unit for i in range(len(value)) for path, unit in
                    matrix_units(sub(energy, length) if i % 3 == 0 else power).items()}
        if key in ('FACE_PORT_COMPOSED_IDENTITY_RESIDUAL', 'FACE_DEPTH_COMPOSED_JOIN_RESIDUAL',
                   'FACE_PORT_RECONSTRUCTION_RESIDUAL'):
            return {(i, *path):unit for i in range(len(value)) for path, unit in matrix_units(power).items()}
        if key == 'FACE_OPEN_DEPTH_COEFFICIENT_RESIDUAL':
            return {(i,):sub(sub(power, potential), potential) for i in range(len(value))}
        if key in ('FACE_PORT_LINEAR_LIFT_RESIDUAL', 'FACE_PORT_ROWS'):
            rho = dims.measure(self.r.symbols['rho_m'])
            velocity = sub(potential, length)
            pressure = add(add(potential, rho), frequency)
            affinity = sub(pressure, rho)
            units = dict(zip(('PRESSURE', 'OUTWARD_VELOCITY', 'RELATIVE_MASS_FLUX', 'AFFINITY',
                             'CHEMICAL_POTENTIAL', 'BULK_VELOCITY', 'MECHANICAL_RESPONSE'),
                            (pressure, velocity, add(rho, velocity), affinity, affinity, velocity, pressure)))
            if key == 'FACE_PORT_LINEAR_LIFT_RESIDUAL':
                return {(i,):tuple(units.values())[i % 7] for i in range(len(value))}
            return {(i, leg, name, j):sub(unit, field_units[j]) for i in range(len(value))
                    for leg in range(2) for name, unit in units.items() for j in range(5)}
        if key == 'FACE_LEG_JOIN_RESIDUAL':
            rho = dims.measure(self.r.symbols['rho_m'])
            velocity = add(potential, tuple(-u for u in length))
            pressure = add(add(potential, rho), frequency)
            affinity = sub(pressure, rho)
            mass_flux = add(rho, velocity)
            units = (pressure, velocity, affinity, mass_flux, mass_flux)
            return {(i, j):sub(units[i % 5], field_units[j]) for i in range(len(value)) for j in range(5)}
        return self.output_units(key, value)

    def output_units(self, key, value):
        """Restore coefficient units from the energy balance and field basis."""
        dims = PHYSICAL_METADATA.dimensions
        energy = dims.measure(self.c.construct(self.anchoring, self.end)['TANGENTIAL_ENERGY_REDUCTION'])
        length, frequency = dims.measure(self.r.z), dims.measure(self.r.omega)
        add = lambda u, v:tuple(a+b for a, b in zip(u, v))
        sub = lambda u, v:tuple(a-b for a, b in zip(u, v))
        power = add(energy, frequency)
        current = add(power, length)
        fields = self.c.field_units
        row_units = [dims.measure(a) for a in self.residual_amplitudes[0]]
        matrix_units = lambda unit:{(5*i+j,):sub(sub(unit, fields[i]), fields[j])
                                   for i in range(5) for j in range(5)}
        if key in ('STORED_ENERGY_MATRIX', 'KINETIC_ENERGY_MATRIX', 'SLAB_ENERGY_MATRIX',
                   'FINITE_TOTAL_ENERGY_MATRIX', 'EQUAL_DEPTH_TOTAL_ENERGY_MATRIX'):
            return matrix_units(energy)
        if key in ('SLAB_CURRENT_MATRIX', 'FINITE_TOTAL_CURRENT_MATRIX', 'EQUAL_DEPTH_TOTAL_CURRENT_MATRIX'):
            return matrix_units(current)
        if key == 'BULK_ENERGY_DENSITY_MATRIX':
            return matrix_units(sub(energy, length))
        if key == 'BULK_DENSITY_BALANCE_RESIDUAL':
            return matrix_units(sub(power, length))
        if key == 'BULK_NORMAL_CURRENT_DENSITY_MATRIX':
            return matrix_units(sub(current, length))
        if key.endswith('_MATRIX') or key in ('SLAB_BALANCE_RESIDUAL', 'FINITE_BALANCE_RESIDUAL',
                                             'EQUAL_DEPTH_BALANCE_RESIDUAL'):
            return matrix_units(power)
        if key == 'PLUS_ROW_POWER_MAP':
            return {(5*i+j,):sub(sub(power, fields[i]), row_units[j]) for i in range(5) for j in range(5)}
        if key == 'MINUS_ROW_POWER_MAP':
            return {(5*i+j,):sub(sub(power, row_units[i]), fields[j]) for i in range(5) for j in range(5)}
        if key == 'CLOSED_PENCIL_LEGS':
            return {(leg, 5*i+j):sub(row_units[i], fields[j]) for leg in range(2)
                    for i in range(5) for j in range(5)}
        if key == 'FACE_RECORDS':
            units = {}
            for i, face in enumerate(value):
                for name in face:
                    if name == 'POLARIZATION_RESIDUALS':
                        units.update({(i, name, j):unit for j, unit in enumerate(
                            (sub(energy, length), sub(current, length), power, power, power, power))})
                    else:
                        unit = (sub(energy, length) if name == 'ENERGY_DENSITY_MATRIX' else
                                sub(current, length) if name == 'NORMAL_CURRENT_DENSITY_MATRIX' else power)
                        units.update({(i, name, *path):unit for path, unit in matrix_units(unit).items()})
            return units
        if key == 'SLAB_POLARIZATION_RESIDUALS':
            return {(i,):unit for i, unit in enumerate((current, energy, energy))}
        if key == 'ROW_POWER_RECONSTRUCTION_RESIDUAL':
            return {():power}
        if key == 'ACOUSTIC_WAVE_ROWS':
            return {(i,):add(frequency, frequency) for i in range(2)}
        if key == 'SOURCE_BRANCH_JOINS':
            return {(i,):dims.zero for i in range(len(value))}
        return None

    def emit(self, anchoring, end, suffix):
        result = self.construct(anchoring, end)
        dims = PHYSICAL_METADATA.dimensions
        for key, value in result.items():
            if key.endswith('_DOMAIN'):
                emit('CLOSED_CURRENT_PAIRING_'+key+'_'+suffix, value)
                emit('METADATA_CLOSED_CURRENT_PAIRING_'+key+'_'+suffix,
                     self.modes.numeric_metadata(value.lhs-value.rhs, lambda path:dims.zero))
            elif isinstance(value, sp.MatrixBase) or key in ('FACE_RECORDS', 'CLOSED_PENCIL_LEGS', 'FACE_LEG_OBJECTS'):
                fingerprinted('CLOSED_CURRENT_PAIRING_'+key+'_'+suffix, cas(value),
                              self.output_units(key, value))
            else:
                physical('CLOSED_CURRENT_PAIRING_'+key+'_'+suffix, value,
                         zero_dimensions=self.output_units(key, value))
        return result



class CurrentSourceControls:
    """Re-enter current construction at reduced rows and the energy input.

    Each instance owns a fresh reduction view, current and closure builder.
    Substitutions are applied before variation, face closure and polarization.
    The original builders and their caches are not changed.
    """

    def __init__(self, baseline, source_energy):
        self.baseline, self.source_energy = baseline, source_energy

    def build(self, replacements):
        from copy import copy
        base = self.baseline
        reduction = copy(base.r)
        reduction.symbols = {name:value.xreplace(replacements)
                             for name,value in base.r.symbols.items()}
        energy_map = {symbol:replacements[base.r.symbols[symbol.name]]
                      for symbol in self.source_energy.free_symbols
                      if symbol.name in base.r.symbols and base.r.symbols[symbol.name] in replacements}
        energy = self.source_energy.xreplace(energy_map)
        strong = base.acoustic.strong.xreplace(replacements)
        ends = copy(base.c.ends)
        ends.r = reduction
        modes = FullPencilModes(ends, base.modes.curl, base.modes.units)
        current = UniformSlabCurrent(reduction, {'value':energy}, ends, strong[3,:])
        balance = SlabEnergyBalance(current)
        acoustic = ClosedAcousticEnergy(balance, modes, strong)
        return ClosedCurrentPairing(acoustic), energy

    @staticmethod
    def spectrum(pairing, bindings):
        """Recompute finite roots and domains of the altered physical pencil."""
        m = pairing.modes
        algebraic, relation, joins = m.analytic(pairing.acoustic.strong)
        physical = algebraic.xreplace(bindings).applyfunc(sp.cancel)
        curve = relation.xreplace(bindings)
        (numerator, denominator), cleared, rows = m.rational_determinant(physical)
        square = sp.solve(curve,m.k**2)[0]
        divisor = sp.Poly(m.k**2-square,m.k)
        remainder = sp.rem(sp.Poly(numerator,m.k),divisor).as_expr()
        if remainder.has(m.k):
            raise NotImplementedError('source-control normal elimination has an odd remainder')
        polynomial = sp.Poly(remainder,m.q)
        if polynomial.is_zero:
            return {'DEFINED':False,'STATUS':'IDENTICALLY_SINGULAR_CONTROL_PENCIL',
                    'PHYSICAL_PENCIL':physical,'POLYNOMIAL':polynomial.as_expr()}
        row_denominator = sp.lcm(rows)
        norm = sp.rem(sp.Poly(sp.expand(row_denominator*row_denominator.xreplace({m.k:-m.k})),m.k),divisor).as_expr()
        exception_polynomials = {
            'DENOMINATOR':sp.gcd(polynomial,sp.Poly(sp.fraction(sp.cancel(norm))[0],m.q)),
            'NORMAL_THRESHOLD':sp.gcd(polynomial,sp.Poly(sp.fraction(sp.cancel(square))[0],m.q)),
            'RADICAL_BRANCH':sp.gcd(polynomial,sp.Poly(m.q,m.q))}
        coverage, roots = EndSpectrumCoverage.isolate(polynomial)
        evaluate = sp.lambdify((m.k,m.q),physical,'numpy',cse=True)
        sheet = BulkSheetPath(curve,m.k,m.q)
        records, paths = [], []
        for index,(q,disk) in enumerate(zip(roots,coverage['ROOT_DISKS'])):
            for sign in (1,-1):
                k = sp.N(sign*sp.sqrt(square.subs(m.q,q)),50)
                matrix = np.asarray(evaluate(complex(k),complex(q)),dtype=complex)
                finite = bool(np.isfinite(matrix).all())
                singular = np.linalg.svd(matrix,compute_uv=False) if finite else None
                threshold = 1e-8*max(1.,singular[0]) if finite else None
                nullity = int(np.sum(singular<threshold)) if finite else None
                membership, path = sheet.classify(complex(k),complex(q))
                paths.append(path)
                record = {'ROOT_DISK_INDEX':index,'NORMAL_LIFT_SIGN':sign,'K':k,'Q':q,
                    'MULTIPLICITY':disk['MULTIPLICITY'],'FINITE_PENCIL':finite,
                    'FIXED_FREQUENCY_SHEET_MEMBERSHIP':membership,
                    'RADICAL_RESIDUAL':sp.N(curve.subs({m.k:k,m.q:q}),25),
                    'ROW_DENOMINATOR_VALUES':tuple(sp.N(v.subs({m.k:k,m.q:q}),25) for v in rows)}
                if finite:
                    record.update({'NULLITY':nullity,'SINGULAR_VALUES':tuple(map(m.number,singular)),
                                   'RANK_THRESHOLD':threshold})
                records.append(record)
        return {'DEFINED':True,'PHYSICAL_PENCIL':physical,'RADICAL_RELATION':curve,'BRANCH_JOIN_RESIDUALS':joins,
            'ELIMINATION_OPERANDS':(cleared,tuple(rows),numerator,denominator,polynomial.as_expr()),
            'EXCEPTION_POLYNOMIALS':{name:p.as_expr() for name,p in exception_polynomials.items()},
            'EXCEPTION_DEGREES':{name:p.degree() for name,p in exception_polynomials.items()},
            'COVERAGE':coverage,'RECORDS':records,'SHEET_PATHS':paths}


class NormalRealityCoverage:
    """Exact axis-root joins in the supplied reference-unit coordinate frame."""

    @staticmethod
    def interval_polynomial(polynomial, interval):
        lo,hi = interval
        result = (sp.S.Zero,sp.S.Zero)
        for coefficient in polynomial.all_coeffs():
            products = [v*w for v in result for w in (lo,hi)]
            result = (min(products)+coefficient,max(products)+coefficient)
        return result

    @classmethod
    def construct(cls, polynomial, curve, denominator, k, q, coverage):
        x,y = sp.symbols('s11cdRealityAxisCoordinate s11cdRealityOtherCoordinate', real=True)
        square = sp.cancel(sp.solve(curve,k**2)[0])
        square_poly = sp.Poly(square,q)
        axes = (('REAL',sp.S.One),('IMAGINARY',sp.I))
        imaginary = sp.expand_complex(square.subs(q,x+sp.I*y)).expand().as_real_imag()[1]
        axis_coefficient = square_poly.nth(2)
        axis_residual = sp.expand(imaginary-2*axis_coefficient*x*y)
        chart = bool(square_poly.degree()==2 and axis_coefficient.is_real and axis_coefficient!=0 and
                     square_poly.nth(1)==0 and square_poly.nth(0).is_real and axis_residual==0)
        factor_scale,factors = sp.sqf_list(polynomial)
        disks = [{str(a):b for a,b in disk} if not isinstance(disk,dict) else disk
                 for disk in coverage['ROOT_DISKS']]
        denominator_norm = sp.Poly(sp.rem(sp.Poly(denominator*denominator.subs(k,-k),k),
            sp.Poly(k**2-square,k)).as_expr(),q)
        denominator_gcd = sp.gcd(polynomial,denominator_norm)
        checks = {'AXIS_DECOMPOSITION_RESIDUAL':axis_residual,
            'WAVE_ELIMINATION_RESIDUAL':sp.cancel(curve.subs(k**2,square)),
            'FACTORIZATION_RESIDUAL':sp.expand(polynomial.as_expr()-factor_scale*
                sp.prod(f.as_expr()**n for f,n in factors)),
            'DEGREE_RESIDUAL':polynomial.degree()-int(coverage['DEGREE']),
            'MULTIPLICITY_RESIDUAL':sum(f.degree()*n for f,n in factors)-int(coverage['COUNT_WITH_MULTIPLICITY'])}
        operands = {'COORDINATE_FRAME':'SUPPLIED_L_T_M_REFERENCE_UNIT_COEFFICIENTS',
            'BOUND_WAVE_COORDINATE':curve,'POLYNOMIAL_COORDINATE':polynomial.as_expr(),
            'NORMAL_SQUARE_COORDINATE':square,'DENOMINATOR_COORDINATE':denominator,
            'DENOMINATOR_NORM_COORDINATE':denominator_norm.as_expr(),
            'DENOMINATOR_GCD_COORDINATE':denominator_gcd.as_expr(),
            'AXIS_IMAGINARY_PART_COORDINATE':imaginary,'AXIS_FACTORIZATION_COORDINATE':sp.factor(imaginary),
            'RADICAL_COORDINATE':x,'OTHER_COORDINATE':y,'REAL_AXIS_CHART_DEFINED':chart}
        axis_records = []
        if chart:
            for factor_index,(factor,multiplicity) in enumerate(factors):
                for axis,multiplier in axes:
                    transformed = sp.expand(factor.monic().as_expr().subs(q,multiplier*x))
                    real,imag = (sp.Poly(v,x,domain=sp.QQ) for v in transformed.as_real_imag())
                    common = sp.gcd(real,imag).monic()
                    intervals = common.intervals(eps=min(d['RADIUS'] for d in disks)**2)
                    entry = {'FACTOR_INDEX':factor_index,'FACTOR_DEGREE':factor.degree(),
                        'MULTIPLICITY':multiplicity,'AXIS':axis,
                        'FACTOR_COORDINATE':factor.monic().as_expr(),
                        'REAL_COORDINATE':real.as_expr(),'IMAGINARY_COORDINATE':imag.as_expr(),
                        'GCD_COORDINATE':common.as_expr(),'GCD_DEGREE':common.degree(),
                        'REAL_ROOT_COUNT':int(common.count_roots(-sp.oo,sp.oo)),
                        'REAL_REMAINDER_COORDINATE':real.rem(common).as_expr(),
                        'IMAGINARY_REMAINDER_COORDINATE':imag.rem(common).as_expr(),
                        'ORIGIN_ROOT':common.eval(0)==0,'ROOTS':[]}
                    axis_square = sp.Poly(square.subs(q,multiplier*x),x,domain=sp.QQ)
                    threshold_gcd = sp.gcd(common,axis_square)
                    entry['NORMAL_SQUARE_COORDINATE'] = axis_square.as_expr()
                    entry['THRESHOLD_GCD_COORDINATE'] = threshold_gcd.as_expr()
                    for interval,root_multiplicity in intervals:
                        left,right = interval
                        bounds = cls.interval_polynomial(axis_square,(left,right))
                        threshold = bool(threshold_gcd.count_roots(left,right))
                        refinements = 0
                        while bounds[0]<=0<=bounds[1] and not threshold and left!=right and refinements<8:
                            left,right = common.refine_root(left,right,eps=(right-left)/100)
                            bounds = cls.interval_polynomial(axis_square,(left,right)); refinements+=1
                        sign = 0 if threshold else 1 if bounds[0]>0 else -1 if bounds[1]<0 else None
                        inside = []
                        for disk_index,disk in enumerate(disks):
                            if int(disk['FACTOR'])!=factor_index:continue
                            margins = tuple(sp.expand(disk['RADIUS']**2-
                                (sp.re(multiplier*t)-sp.re(disk['CENTER']))**2-
                                (sp.im(multiplier*t)-sp.im(disk['CENTER']))**2) for t in (left,right))
                            if all(v>0 for v in margins):inside.append((disk_index,margins))
                        entry['ROOTS'].append({'RADICAL_INTERVAL':(left,right),
                            'INTERVAL_ROOT_COUNT':int(common.count_roots(left,right)),
                            'ROOT_MULTIPLICITY':root_multiplicity,'NORMAL_SQUARE_INTERVAL':bounds,
                            'NORMAL_SQUARE_SIGN':sign if sign is not None else 'UNRESOLVED',
                            'NORMAL_THRESHOLD_ROOT':threshold,'ORIGIN_ROOT':bool(common.eval(0)==0 and left<=0<=right),
                            'REFINEMENTS':refinements,'DISK_JOINS':tuple(inside),
                            'UNIQUE_DISK_JOIN':len(inside)==1})
                    entry['REAL_ROOT_COUNT_RESIDUAL'] = sum(v['ROOT_MULTIPLICITY'] for v in entry['ROOTS'])-entry['REAL_ROOT_COUNT']
                    check_prefix = str(factor_index)+'_'+axis+'_'
                    for name in ('REAL_REMAINDER_COORDINATE','IMAGINARY_REMAINDER_COORDINATE','REAL_ROOT_COUNT_RESIDUAL'):
                        checks[check_prefix+name] = entry[name]
                    for ri,root in enumerate(entry['ROOTS']):
                        checks[check_prefix+'INTERVAL_'+str(ri)+'_ROOT_COUNT_RESIDUAL'] = root['INTERVAL_ROOT_COUNT']-1
                    axis_records.append(entry)
        root_records = []
        for disk_index,disk in enumerate(disks):
            real_clearance = abs(sp.im(disk['CENTER']))-disk['RADIUS']
            imag_clearance = abs(sp.re(disk['CENTER']))-disk['RADIUS']
            matches = [(ai,ri) for ai,axis in enumerate(axis_records) for ri,root in enumerate(axis['ROOTS'])
                       if root['UNIQUE_DISK_JOIN'] and root['DISK_JOINS'][0][0]==disk_index]
            disk_valid = bool(disk['ONE_ROOT_DISK'] and coverage['FINITE_POLYNOMIAL_ROOT_COVERAGE'])
            signs = {axis_records[ai]['ROOTS'][ri]['NORMAL_SQUARE_SIGN'] for ai,ri in matches}
            status,nonzero = 'UNRESOLVED', 'UNRESOLVED'
            if chart and disk_valid:
                if real_clearance>0 and imag_clearance>0:
                    status,nonzero = 'PROVED_NONREAL',True
                elif signs=={1}:
                    status,nonzero = 'PROVED_REAL',True
                elif signs=={-1}:
                    status,nonzero = 'PROVED_NONREAL',True
                elif signs=={0}:
                    status,nonzero = 'PROVED_REAL',False
            root_records.append({'ROOT_DISK_INDEX':disk_index,'FACTOR_INDEX':int(disk['FACTOR']),
                'CENTER':disk['CENTER'],'RADIUS':disk['RADIUS'],
                'REAL_AXIS_CLEARANCE':real_clearance,'IMAGINARY_AXIS_CLEARANCE':imag_clearance,
                'AXIS_ROOT_MATCHES':tuple(matches),'NORMAL_REALITY_STATUS':status,
                'NORMAL_NONZERO':nonzero,'DENOMINATOR_EXCLUDED':denominator_gcd.degree()==0,
                'CERTIFIED_DISK_INPUT':disk_valid})
        return {'OPERANDS':operands,'CHECKS':checks,'AXES':axis_records,'DISKS':root_records}

    @staticmethod
    def unit(path, frequency_unit, length_unit):
        key = next((v for v in reversed(path) if isinstance(v,str)),None)
        if key in ('RADICAL_INTERVAL','CENTER','RADIUS','REAL_AXIS_CLEARANCE','IMAGINARY_AXIS_CLEARANCE'):
            return frequency_unit
        if key=='NORMAL_SQUARE_INTERVAL':return tuple(-2*v for v in length_unit)
        # Disk-join entries are (dimensionless disk index, squared-distance margins).
        if 'DISK_JOINS' in path and len(path)>=2 and path[-2]==1:
            return tuple(2*v for v in frequency_unit)
        return (0,0,0)


class ModalCurrentSubspaces:
    """Full native mode spaces and the computed physical energy forms.

    Numeric arrays are coefficients in the declared reference-unit frame.
    The adjoint row pairing, physical right-field current, and row-power
    bridge are separate objects. No Euclidean overlap is used as a current.
    """

    def __init__(self, pairing, result, material_bindings):
        self.pairing, self.result = pairing, result
        self.modes, self.r = pairing.modes, pairing.r
        self.bindings = dict(material_bindings)
        self.epsilon = self.r.symbols['epsilon_shape']
        self.variables = (*pairing.frequencies, *pairing.c.leg_momenta, *pairing.acoustic.qlegs,
                          pairing.acoustic.height)
        self.field_units = pairing.c.field_units
        d = PHYSICAL_METADATA.dimensions
        self.row_units = [d.measure(v) for v in pairing.residual_amplitudes[0]]
        self.energy_unit = d.measure(pairing.c.construct(pairing.anchoring, pairing.end)['TANGENTIAL_ENERGY_REDUCTION'])
        self.frequency_unit, self.length_unit = d.measure(self.r.omega), d.measure(self.r.z)
        self.power_unit = tuple(a+b for a,b in zip(self.energy_unit, self.frequency_unit))
        self.current_unit = tuple(a+b for a,b in zip(self.power_unit, self.length_unit))

    def prepare(self):
        p, result = self.pairing, self.result
        sources = {'ENERGY_SLAB':'SLAB_ENERGY_MATRIX', 'CURRENT_SLAB':'SLAB_CURRENT_MATRIX',
                   'ENERGY_BULK':'BULK_ENERGY_DENSITY_MATRIX', 'CURRENT_BULK':'BULK_NORMAL_CURRENT_DENSITY_MATRIX',
                   'CURRENT_DEPTH':'BULK_DEPTH_CURRENT_MATRIX', 'POWER_INTERFACE':'INTERFACE_POWER_MATRIX',
                   'POWER_PORT':'PORT_POWER_MATRIX', 'POWER_SOURCE':'SOURCE_POWER_MATRIX',
                   'POWER_MAP_PLUS':'PLUS_ROW_POWER_MAP', 'POWER_MAP_MINUS':'MINUS_ROW_POWER_MAP'}
        expressions, residuals = {}, {}
        for name, key in sources.items():
            def coefficient(value):
                return dict(polynomial_terms(value, (self.epsilon,))).get((2,), sp.S.Zero)
            expressions[name] = result[key].applyfunc(coefficient)
            residuals[key] = (result[key]-self.epsilon**2*expressions[name]).applyfunc(p.carrier_expansion)
        expressions.update({'PENCIL_PLUS':result['CLOSED_PENCIL_LEGS'][0],
                            'PENCIL_MINUS':result['CLOSED_PENCIL_LEGS'][1]})
        qleft, qright = p.acoustic.qlegs
        wave = result['ACOUSTIC_WAVE_ROWS'][0]
        transports = {}
        for name, variable in (('NORMAL', p.c.leg_momenta[1]), ('FREQUENCY', p.frequencies[1])):
            transport = sp.cancel(-sp.diff(wave,variable)/sp.diff(wave,qright))
            transports[name] = transport
            for key, value in tuple(expressions.items()):
                if key.startswith(('NORMAL_', 'FREQUENCY_')):
                    continue
                expressions[name+'_'+key] = value.diff(variable)+transport*value.diff(qright)
        scalar_expressions = {'TIME_RATE':result['BEAT_RATES'][0], 'NORMAL_RATE':result['BEAT_RATES'][1],
            'DEPTH_RATE':result['DEPTH_PHASE_RATE'], 'DEPTH_INTEGRAL':result['GENERIC_DEPTH_INTEGRAL'],
            'DEPTH_INTEGRAL_Q':sp.diff(result['GENERIC_DEPTH_INTEGRAL'],qright),
            'EQUAL_DEPTH_INTEGRAL':result['EQUAL_DEPTH_INTEGRAL'],
            'EQUAL_DEPTH_INTEGRAL_Q':sp.limit(sp.diff(result['GENERIC_DEPTH_INTEGRAL'],qright),qright,qleft),
            'TOP_FACTOR':result['DEPTH_PHASE_FACTOR'].subs(p.acoustic.depth,p.acoustic.height)}
        old_bulk = p.acoustic.construct(p.anchoring,p.end)
        rate_map = dict(old_bulk['DEPTH_DECAY_BINDINGS'])
        decay = next(v for v in rate_map if 'Decay' in v.name)
        oscillation = next(v for v in rate_map if 'Oscillation' in v.name)
        infinite = old_bulk['CONVERGENT_INFINITE_DEPTH_INTEGRAL'].xreplace(
            {-decay+sp.I*oscillation:result['DEPTH_PHASE_RATE']})
        scalar_expressions['INFINITE_DEPTH_INTEGRAL'] = infinite
        scalar_expressions['RADICAL_SCALE'] = old_bulk['ACOUSTIC_RADICAL_SCALE']
        for name, variable in (('NORMAL', p.c.leg_momenta[1]), ('FREQUENCY', p.frequencies[1])):
            scalar_expressions[name+'_TRANSPORT'] = transports[name]
            for label, value in zip(('TIME_WEIGHT','NORMAL_WEIGHT','DEPTH_WEIGHT'),
                    (result['BEAT_RATES'][0],result['BEAT_RATES'][1],result['DEPTH_PHASE_RATE'])):
                scalar_expressions[name+'_'+label] = value.diff(variable)+transports[name]*value.diff(qright)
        bound = {key:value.xreplace(self.bindings) for key,value in {**expressions,**scalar_expressions}.items()}
        unresolved = set().union(*(value.free_symbols for value in bound.values()))-set(self.variables)
        if unresolved:
            raise ValueError(('unbound modal-current inputs',sorted(map(str,unresolved))))
        self.evaluate = {key:sp.lambdify(self.variables,value,'numpy',cse=True) for key,value in bound.items()}
        self.coefficient_residuals = residuals
        self.symbolic_operands = expressions
        self.scalar_operands = scalar_expressions
        return residuals

    def exact_low_degree_lifts(self, records, coverage):
        """Recover small factors from the actual bound pencil, without rerooting the high-degree factor."""
        m = self.modes
        algebraic, relation, _ = m.analytic(self.pairing.acoustic.strong)
        fixed = {**self.bindings,self.r.omega:self.bindings[self.r.omega]}
        physical = algebraic.xreplace(fixed).applyfunc(sp.cancel)
        curve = relation.xreplace(fixed)
        (numerator, denominator), _, _ = m.rational_determinant(physical)
        k_square = sp.solve(curve,m.k**2)[0]
        polynomial = sp.Poly(sp.rem(sp.Poly(numerator,m.k),sp.Poly(m.k**2-k_square,m.k)).as_expr(),m.q)
        factors = sp.sqf_list(polynomial)[1]
        self.normal_reality = NormalRealityCoverage.construct(polynomial,curve,denominator,m.k,m.q,coverage)
        exact = {}
        for index,disk in enumerate(coverage['ROOT_DISKS']):
            disk = {str(k):v for k,v in disk}
            factor,multiplicity = factors[int(disk['FACTOR'])]
            if factor.degree()>2:
                continue
            candidates = sp.solve(factor.as_expr(),m.q)
            selected = [q for q in candidates if abs(complex(sp.N(q-disk['CENTER'],60)))<float(disk['RADIUS'])]
            if len(selected)!=1:
                raise ValueError(('small-factor isolating disk join',index))
            q = selected[0]
            for record_index,record in enumerate(records):
                if int(record['ROOT_DISK_INDEX'])!=index:
                    continue
                k = sp.simplify(record['NORMAL_LIFT_SIGN']*sp.sqrt(k_square.subs(m.q,q)))
                exact[record_index] = {'K':k,'Q':q,'FACTOR':factor.as_expr(),'MULTIPLICITY':multiplicity,
                    'FACTOR_RESIDUAL':sp.simplify(factor.eval(q)),
                    'WAVE_RESIDUAL':sp.simplify(curve.subs({m.k:k,m.q:q})),
                    'PRODUCER_K_DIFFERENCE':sp.N(k-record['K'],30),
                    'PRODUCER_Q_DIFFERENCE':sp.N(q-record['Q'],30),
                    'REAL_K':sp.im(k)==0,'NONZERO_K':sp.Ne(k,0),'POSITIVE_Q_IMAGINARY_PART':sp.im(q)>0}
        return exact, {'POLYNOMIAL_DEGREE_RESIDUAL':polynomial.degree()-int(coverage['DEGREE']),
                       'ROOT_MULTIPLICITY_RESIDUAL':sum(f.degree()*n for f,n in factors)-int(coverage['COUNT_WITH_MULTIPLICITY'])}

    @staticmethod
    def norm(value):
        return float(np.linalg.norm(value))

    def construct(self, records, coverage, frequency, cutoff, progress=lambda record:None):
        self.prepare()
        exact, exact_residuals = self.exact_low_degree_lifts(records,coverage)
        outputs = []
        for index, source in enumerate(records):
            k = complex(exact[index]['K'] if index in exact else source['K'])
            native_q = complex(exact[index]['Q'] if index in exact else source['Q'])
            scale = complex(self.evaluate['RADICAL_SCALE'](frequency,frequency,k.conjugate(),k,0,0,cutoff))
            q = scale*native_q
            point = (frequency,frequency,k.conjugate(),k,q.conjugate(),q,cutoff)
            equal = q.conjugate()==q
            excluded = {'INFINITE_DEPTH_INTEGRAL'} | ({'DEPTH_INTEGRAL','DEPTH_INTEGRAL_Q'} if equal else set())
            values = {key:np.asarray(evaluate(*point),dtype=complex) for key,evaluate in self.evaluate.items() if key not in excluded}
            integral = complex(values['EQUAL_DEPTH_INTEGRAL' if equal else 'DEPTH_INTEGRAL'])
            integral_q = complex(values['EQUAL_DEPTH_INTEGRAL_Q' if equal else 'DEPTH_INTEGRAL_Q'])
            matrix = values['PENCIL_PLUS']
            u,singular,vh = np.linalg.svd(matrix)
            threshold = 1e-8*max(1.,singular[0])
            nullity = int(np.sum(singular<threshold))
            record = {'INDEX':index,'ROOT_DISK_INDEX':int(source['ROOT_DISK_INDEX']),
                'NORMAL_LIFT_SIGN':int(source['NORMAL_LIFT_SIGN']),'K':self.modes.number(k),
                'Q':self.modes.number(native_q),'PHYSICAL_Q':self.modes.number(q),
                'OMEGA':self.modes.number(frequency),'DEPTH_CUTOFF':self.modes.number(cutoff),
                'DEPTH_INTEGRAL':self.modes.number(integral),'EQUAL_DEPTH_BRANCH':equal,
                'SHEET_MEMBERSHIP':source['FIXED_FREQUENCY_SHEET_MEMBERSHIP'],
                'SINGULAR_VALUES':tuple(map(self.modes.number,singular)), 'RANK_THRESHOLD':threshold,
                'NULLITY':nullity,'PRODUCER_NULLITY_RESIDUAL':nullity-int(source['NULLITY']),
                'EXACT_LOW_DEGREE_LIFT':exact.get(index, {'STATUS':'HIGHER_DEGREE_ISOLATED_NUMERIC_ROOT'}),
                'RADICAL_TRANSPORT_DENOMINATOR':self.modes.number(q),
                'BULK_DECAY_NUMERIC':q.imag>0,
                'NORMAL_REALITY_CERTIFICATE':self.normal_reality['DISKS'][int(source['ROOT_DISK_INDEX'])],
                'EXACT_REAL_NORMAL':self.normal_reality['DISKS'][int(source['ROOT_DISK_INDEX'])]['NORMAL_REALITY_STATUS']=='PROVED_REAL',
                'FORMS':{},'RESIDUALS':{},
                'OPERANDS':{key:values[key] for key in ('PENCIL_PLUS','PENCIL_MINUS','NORMAL_PENCIL_PLUS',
                    'FREQUENCY_PENCIL_PLUS','POWER_MAP_PLUS','POWER_MAP_MINUS','ENERGY_SLAB','CURRENT_SLAB',
                    'ENERGY_BULK','CURRENT_BULK','CURRENT_DEPTH','POWER_INTERFACE')}}
            if not nullity:
                outputs.append(record);progress({'mode':index,'nullity':nullity});continue
            right,left = vh.conj().T[:,-nullity:],u[:,-nullity:]
            n_frequency = left.conj().T@values['FREQUENCY_PENCIL_PLUS']@right
            n_normal = left.conj().T@values['NORMAL_PENCIL_PLUS']@right
            n_threshold = 1e-9*max(1.,np.linalg.norm(n_frequency,2))
            normal_threshold = 1e-9*max(1.,np.linalg.norm(n_normal,2))
            record.update({'RIGHT_BASIS_RANK':int(np.linalg.matrix_rank(right,tol=1e-9)),
                           'LEFT_BASIS_RANK':int(np.linalg.matrix_rank(left,tol=1e-9)),
                           'FREQUENCY_PAIRING_RANK':int(np.linalg.matrix_rank(n_frequency,tol=n_threshold)),
                           'FREQUENCY_PAIRING_THRESHOLD':n_threshold,
                           'NORMAL_PAIRING_RANK':int(np.linalg.matrix_rank(n_normal,tol=normal_threshold)),
                           'NORMAL_PAIRING_THRESHOLD':normal_threshold})
            forms,residuals = record['FORMS'],record['RESIDUALS']
            forms.update({'RIGHT':right,'LEFT':left,'N_FREQUENCY':n_frequency,'N_NORMAL':n_normal,
                          'RIGHT_COORDINATE_PROJECTOR':right@right.conj().T,
                          'LEFT_COORDINATE_PROJECTOR':left@left.conj().T})
            residuals.update({'RIGHT_KERNEL':matrix@right,'LEFT_KERNEL':matrix.conj().T@left,
                'RIGHT_BASIS':right.conj().T@right-np.eye(nullity),'LEFT_BASIS':left.conj().T@left-np.eye(nullity),
                'RIGHT_PROJECTOR':forms['RIGHT_COORDINATE_PROJECTOR']@forms['RIGHT_COORDINATE_PROJECTOR']-forms['RIGHT_COORDINATE_PROJECTOR'],
                'LEFT_PROJECTOR':forms['LEFT_COORDINATE_PROJECTOR']@forms['LEFT_COORDINATE_PROJECTOR']-forms['LEFT_COORDINATE_PROJECTOR'],
                'CONJUGATE_PENCIL':values['PENCIL_MINUS']-matrix.conjugate(),
                'CONJUGATE_POWER_MAP':values['POWER_MAP_MINUS']-values['POWER_MAP_PLUS'].conj().T})
            contract = lambda value:right.conj().T@value@right
            energy = values['ENERGY_SLAB']+integral*values['ENERGY_BULK']
            current = values['CURRENT_SLAB']+integral*values['CURRENT_BULK']
            top = complex(values['TOP_FACTOR'])*values['CURRENT_DEPTH']
            forms.update({key:contract(values[key]) for key in ('ENERGY_SLAB','CURRENT_SLAB','ENERGY_BULK',
                'CURRENT_BULK','CURRENT_DEPTH','POWER_INTERFACE','POWER_PORT','POWER_SOURCE')})
            forms.update({'ENERGY_FINITE':contract(energy),'CURRENT_FINITE':contract(current),'POWER_TOP':contract(top)})
            residuals['FINITE_BALANCE'] = contract(complex(values['TIME_RATE'])*energy+
                complex(values['NORMAL_RATE'])*current+values['POWER_INTERFACE']+top-values['POWER_SOURCE'])
            residuals['CURRENT_HERMITIAN'] = forms['CURRENT_FINITE']-forms['CURRENT_FINITE'].conj().T
            residuals['ENERGY_HERMITIAN'] = forms['ENERGY_FINITE']-forms['ENERGY_FINITE'].conj().T
            power_covectors = right.conj().T@values['POWER_MAP_PLUS']
            bridge = power_covectors@left
            defect = power_covectors-bridge@left.conj().T
            forms.update({'POWER_LEFT_BRIDGE':bridge,'POWER_LEFT_DEFECT':defect})
            residuals['POWER_LEFT_SPLIT'] = power_covectors-bridge@left.conj().T-defect
            for label in ('NORMAL','FREQUENCY'):
                transport = complex(values[label+'_TRANSPORT'])
                de = values[label+'_ENERGY_SLAB']+integral*values[label+'_ENERGY_BULK']+integral_q*transport*values['ENERGY_BULK']
                dj = values[label+'_CURRENT_SLAB']+integral*values[label+'_CURRENT_BULK']+integral_q*transport*values['CURRENT_BULK']
                dtop = complex(values['TOP_FACTOR'])*(values[label+'_CURRENT_DEPTH']+
                    cutoff*complex(values[label+'_DEPTH_WEIGHT'])*values['CURRENT_DEPTH'])
                balance_derivative = (complex(values[label+'_TIME_WEIGHT'])*energy+complex(values['TIME_RATE'])*de+
                    complex(values[label+'_NORMAL_WEIGHT'])*current+complex(values['NORMAL_RATE'])*dj+
                    values[label+'_POWER_INTERFACE']+dtop)
                source_derivative = values[label+'_POWER_SOURCE']
                pencil_pair = left.conj().T@values[label+'_PENCIL_PLUS']@right
                weighted_pair = power_covectors@values[label+'_PENCIL_PLUS']@right
                defect_pair = defect@values[label+'_PENCIL_PLUS']@right
                other_source_terms = source_derivative-values['POWER_MAP_PLUS']@values[label+'_PENCIL_PLUS']
                isolated_weight = complex(values[label+('_NORMAL_WEIGHT' if label=='NORMAL' else '_TIME_WEIGHT')])
                isolated_operand = current if label=='NORMAL' else energy
                other_balance_terms = contract(balance_derivative-isolated_weight*isolated_operand)
                reconstructed = (bridge@pencil_pair+defect_pair+contract(other_source_terms)-other_balance_terms)/isolated_weight
                forms.update({label+'_BALANCE_DERIVATIVE':contract(balance_derivative),
                    label+'_SOURCE_DERIVATIVE':contract(source_derivative),label+'_WEIGHTED_PENCIL_PAIRING':weighted_pair,
                    label+'_LEFT_DEFECT_PAIRING':defect_pair,label+'_OTHER_SOURCE_TERMS':contract(other_source_terms),
                    label+'_OTHER_BALANCE_TERMS':other_balance_terms,
                    label+'_RECONSTRUCTED_FORM':reconstructed})
                residuals[label+'_BALANCE_DERIVATIVE'] = contract(balance_derivative-source_derivative)
                residuals[label+'_ROW_PAIRING_BRIDGE'] = weighted_pair-bridge@pencil_pair-defect_pair
                residuals[label+'_CURRENT_ENERGY_RECONSTRUCTION'] = reconstructed-contract(isolated_operand)
            if record['FREQUENCY_PAIRING_RANK']==nullity:
                normalized_left = left@np.linalg.inv(n_frequency).conj().T
                forms['LEFT_FREQUENCY_NORMALIZED'] = normalized_left
                forms['NORMAL_PAIRING_FREQUENCY_NORMALIZED'] = normalized_left.conj().T@values['NORMAL_PENCIL_PLUS']@right
                forms['POWER_LEFT_BRIDGE_FREQUENCY_NORMALIZED'] = bridge@n_frequency
                residuals['NORMALIZED_POWER_LEFT_SPLIT'] = (power_covectors-
                    forms['POWER_LEFT_BRIDGE_FREQUENCY_NORMALIZED']@normalized_left.conj().T-defect)
                residuals['FREQUENCY_NORMALIZATION'] = normalized_left.conj().T@values['FREQUENCY_PENCIL_PLUS']@right-np.eye(nullity)
                residuals['NORMALIZED_LEFT_KERNEL'] = matrix.conj().T@normalized_left
            disk = {str(k):v for k,v in coverage['ROOT_DISKS'][int(source['ROOT_DISK_INDEX'])]}
            decay_certified = bool(scale.real>0 and scale.imag==0 and sp.im(disk['CENTER'])-disk['RADIUS']>0)
            record['BULK_DECAY_DISK_CERTIFIED'] = decay_certified
            certificate = record['NORMAL_REALITY_CERTIFICATE']
            regular_normal = bool(certificate['NORMAL_NONZERO'] is True and certificate['DENOMINATOR_EXCLUDED'] and
                                  record['NORMAL_PAIRING_RANK']==nullity)
            record['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED'] = bool(regular_normal and decay_certified and record['EXACT_REAL_NORMAL'] and
                source['FIXED_FREQUENCY_SHEET_MEMBERSHIP']==sp.true and record['FREQUENCY_PAIRING_RANK']==nullity)
            if decay_certified:
                infinite_integral = complex(self.evaluate['INFINITE_DEPTH_INTEGRAL'](*point))
                record['INFINITE_DEPTH_INTEGRAL'] = self.modes.number(infinite_integral)
                infinite_current = contract(values['CURRENT_SLAB']+infinite_integral*values['CURRENT_BULK'])
                forms['CURRENT_INFINITE'] = infinite_current
                residuals['INFINITE_CURRENT_HERMITIAN'] = infinite_current-infinite_current.conj().T
                hermitian = (infinite_current+infinite_current.conj().T)/2
                eigenvalues,rotation = np.linalg.eigh(hermitian)
                forms.update({'CURRENT_EIGENVALUES':eigenvalues,'CURRENT_ROTATION':rotation})
                residuals['CURRENT_DIAGONALIZATION'] = rotation.conj().T@infinite_current@rotation-np.diag(eigenvalues)
                current_threshold = 1e-9*max(1.,np.linalg.norm(infinite_current,2))
                record['CURRENT_RANK'] = int(np.sum(np.abs(eigenvalues)>current_threshold))
                record['CURRENT_RANK_THRESHOLD'] = current_threshold
                eligible = (regular_normal and record['EXACT_REAL_NORMAL'] and source['FIXED_FREQUENCY_SHEET_MEMBERSHIP']==sp.true and
                    record['CURRENT_RANK']==nullity and record['FREQUENCY_PAIRING_RANK']==nullity and
                    self.norm(residuals['INFINITE_CURRENT_HERMITIAN'])<=1e-9*max(1.,self.norm(infinite_current)))
                record['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED'] = bool(eligible)
                if eligible:
                    transform = rotation@np.diag(1/np.sqrt(np.abs(eigenvalues)))
                    flux_right = right@transform
                    flux_left = normalized_left@np.linalg.inv(transform).conj().T
                    forms.update({'FIELD_TO_FLUX_MAP':transform,'FLUX_RIGHT':flux_right,'FLUX_LEFT':flux_left,
                                  'SIGNED_CURRENT':np.diag(np.sign(eigenvalues))})
                    residuals['SIGNED_CURRENT_NORMALIZATION'] = transform.conj().T@infinite_current@transform-forms['SIGNED_CURRENT']
                    residuals['FLUX_FREQUENCY_NORMALIZATION'] = flux_left.conj().T@values['FREQUENCY_PENCIL_PLUS']@flux_right-np.eye(nullity)
                    residuals['FLUX_RIGHT_KERNEL'] = matrix@flux_right
                    residuals['FLUX_LEFT_KERNEL'] = matrix.conj().T@flux_left
            record['SCALED_RIGHT_KERNEL_NORM'] = self.norm(residuals['RIGHT_KERNEL'])/max(1.,self.norm(matrix)*self.norm(right))
            record['SCALED_LEFT_KERNEL_NORM'] = self.norm(residuals['LEFT_KERNEL'])/max(1.,self.norm(matrix)*self.norm(left))
            record['RESIDUAL_NORMS'] = {key:self.norm(value) for key,value in residuals.items()}
            outputs.append(record)
            progress({'mode':index,'nullity':nullity,'frequencyRank':record['FREQUENCY_PAIRING_RANK'],
                'finiteBalance':record['RESIDUAL_NORMS']['FINITE_BALANCE'],
                'normalDerivative':record['RESIDUAL_NORMS']['NORMAL_BALANCE_DERIVATIVE'],
                'frequencyDerivative':record['RESIDUAL_NORMS']['FREQUENCY_BALANCE_DERIVATIVE'],
                'physicalCurrentNormalization':record.get('PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED',False)})
        return {'RECORDS':outputs,'EXACT_SPECTRUM_RECONSTRUCTION_RESIDUALS':exact_residuals,
                'NORMAL_REALITY_COVERAGE':self.normal_reality}

    def tensor_unit(self, group, key, path, nullity):
        d = PHYSICAL_METADATA.dimensions
        add = lambda u,v:tuple(a+b for a,b in zip(u,v))
        sub = lambda u,v:tuple(a-b for a,b in zip(u,v))
        negative = lambda u:tuple(-a for a in u)
        half_current = tuple(a/2 for a in self.current_unit)
        field = self.field_units
        row = self.row_units
        i = path[0] if path else 0
        if group=='OPERANDS':
            a,b = i//5,i%5
            if key.endswith('PENCIL_PLUS') or key=='PENCIL_MINUS':
                unit = sub(row[a],field[b])
                return sub(unit,negative(self.length_unit) if key.startswith('NORMAL_') else
                           self.frequency_unit) if key.startswith(('NORMAL_','FREQUENCY_')) else unit
            if key=='POWER_MAP_PLUS':return sub(sub(self.power_unit,field[a]),row[b])
            if key=='POWER_MAP_MINUS':return sub(sub(self.power_unit,row[a]),field[b])
            unit = (self.energy_unit if key=='ENERGY_SLAB' else self.current_unit if key=='CURRENT_SLAB' else
                    sub(self.energy_unit,self.length_unit) if key=='ENERGY_BULK' else self.power_unit)
            return sub(sub(unit,field[a]),field[b])
        if key in ('RIGHT','RIGHT_KERNEL','FLUX_RIGHT','FLUX_RIGHT_KERNEL'):
            unit = row[i//nullity] if key.endswith('KERNEL') else field[i//nullity]
            return sub(unit,half_current) if key.startswith('FLUX_') else unit
        if key in ('LEFT','LEFT_KERNEL','LEFT_FREQUENCY_NORMALIZED','NORMALIZED_LEFT_KERNEL','FLUX_LEFT','FLUX_LEFT_KERNEL'):
            unit = negative(field[i//nullity] if key.endswith('KERNEL') else row[i//nullity])
            if key not in ('LEFT','LEFT_KERNEL'):unit=add(unit,self.frequency_unit)
            return add(unit,half_current) if key.startswith('FLUX_') else unit
        if key=='N_FREQUENCY':return negative(self.frequency_unit)
        if key=='N_NORMAL':return self.length_unit
        if key=='NORMAL_PAIRING_FREQUENCY_NORMALIZED':return add(self.length_unit,self.frequency_unit)
        if key=='FIELD_TO_FLUX_MAP':return negative(half_current)
        if key=='CONJUGATE_PENCIL':return sub(row[i//5],field[i%5])
        if key=='CONJUGATE_POWER_MAP':return sub(sub(self.power_unit,row[i//5]),field[i%5])
        if key in ('POWER_LEFT_DEFECT','POWER_LEFT_SPLIT','NORMALIZED_POWER_LEFT_SPLIT'):
            return sub(self.power_unit,row[i%5])
        if key=='POWER_LEFT_BRIDGE_FREQUENCY_NORMALIZED':return self.energy_unit
        if key in ('CURRENT_FINITE','CURRENT_INFINITE','CURRENT_SLAB','CURRENT_EIGENVALUES','CURRENT_HERMITIAN',
                   'INFINITE_CURRENT_HERMITIAN','CURRENT_DIAGONALIZATION') or key.startswith('NORMAL_'):
            return self.current_unit
        if key in ('ENERGY_FINITE','ENERGY_SLAB','ENERGY_HERMITIAN') or key.startswith('FREQUENCY_') and key not in ('FREQUENCY_NORMALIZATION',):
            return self.energy_unit
        if key=='ENERGY_BULK':return sub(self.energy_unit,self.length_unit)
        if key.startswith('POWER_') or key in ('CURRENT_BULK','CURRENT_DEPTH','FINITE_BALANCE'):
            return self.power_unit
        return d.zero

    @staticmethod
    def quadratic_quantity(group,key):
        if group=='OPERANDS':return 'PENCIL' not in key
        if key in ('CURRENT_ROTATION',):return False
        if key=='NORMAL_PAIRING_FREQUENCY_NORMALIZED':return False
        if key in ('FREQUENCY_NORMALIZATION','FLUX_FREQUENCY_NORMALIZATION'):return False
        return key.startswith(('ENERGY_','CURRENT_','POWER_','NORMAL_','FREQUENCY_','INFINITE_CURRENT_')) or key in (
            'FINITE_BALANCE','CONJUGATE_POWER_MAP','NORMALIZED_POWER_LEFT_SPLIT','SIGNED_CURRENT','SIGNED_CURRENT_NORMALIZATION')

    def info_unit(self, path, record):
        d=PHYSICAL_METADATA.dimensions
        key=path[-1]
        if 'NORMAL_REALITY_CERTIFICATE' in path:
            return NormalRealityCoverage.unit(path,self.frequency_unit,self.length_unit)
        if key=='NORMAL_PAIRING_THRESHOLD':return self.length_unit
        if key in ('FACTOR','FACTOR_RESIDUAL'):
            degree=sp.Poly(record['EXACT_LOW_DEGREE_LIFT']['FACTOR'],self.modes.q).degree()
            return tuple(degree*v for v in self.frequency_unit)
        if key in ('K','PRODUCER_K_DIFFERENCE'):return tuple(-v for v in self.length_unit)
        if key in ('Q','OMEGA','PRODUCER_Q_DIFFERENCE'):return self.frequency_unit
        if key in ('PHYSICAL_Q','RADICAL_TRANSPORT_DENOMINATOR'):return tuple(-v for v in self.length_unit)
        if key in ('DEPTH_CUTOFF','DEPTH_INTEGRAL','INFINITE_DEPTH_INTEGRAL'):return self.length_unit
        if key=='WAVE_RESIDUAL':return tuple(2*v for v in self.frequency_unit)
        if key=='FREQUENCY_PAIRING_THRESHOLD':return tuple(-v for v in self.frequency_unit)
        if key=='CURRENT_RANK_THRESHOLD':return self.current_unit
        return d.zero

    def emit(self, computed, provenance, prefix='MODAL_SUBSPACE', context='REFERENCE_LAB_HELD_RHO4_CONSTANT'):
        d=PHYSICAL_METADATA.dimensions
        def convert(value):
            if isinstance(value,np.ndarray):
                if value.ndim==0:return self.modes.number(value.item())
                shape=(len(value),1) if value.ndim==1 else value.shape
                return sp.ImmutableMatrix(*shape,[self.modes.number(v) for v in value.ravel()])
            if isinstance(value,dict):return {key:convert(v) for key,v in value.items()}
            if isinstance(value,(tuple,list)):return tuple(convert(v) for v in value)
            if isinstance(value,(complex,np.complexfloating)):return self.modes.number(value)
            return value
        def output(tag,value,unit,heavy=False):
            body=cas(convert(value))
            emit(tag,carrier_fingerprint(body) if heavy=='carrier' else self.modes.compact_fingerprint(body) if heavy else body)
            emit('METADATA_'+tag,self.modes.numeric_metadata(body,unit))
        output(prefix+'_PROVENANCE',provenance,lambda path:d.zero)
        output(prefix+'_EXACT_SPECTRUM_RECONSTRUCTION_RESIDUALS',
               computed['EXACT_SPECTRUM_RECONSTRUCTION_RESIDUALS'],lambda path:d.zero)
        for key,value in computed['NORMAL_REALITY_COVERAGE'].items():
            output(prefix+'_NORMAL_REALITY_'+key,value,
                   lambda path:NormalRealityCoverage.unit(path,self.frequency_unit,self.length_unit),'carrier' if key in ('OPERANDS','AXES') else False)
        for key,value in self.coefficient_residuals.items():
            units=self.pairing.output_units(key,value)
            output(prefix+'_QUADRATIC_EXTRACTION_'+key,value,lambda path:units[path])
        for record in computed['RECORDS']:
            tag=prefix+'_'+context+'_'+str(record['INDEX'])
            info={key:value for key,value in record.items() if key not in ('FORMS','RESIDUALS','OPERANDS')}
            output(tag+'_RECORD',info,lambda path:self.info_unit(path,info))
            for group in ('OPERANDS','FORMS','RESIDUALS'):
                for key,value in record[group].items():
                    value=convert(value)
                    if self.quadratic_quantity(group,key):value=self.epsilon**2*value
                    output(tag+'_'+group+'_'+key,value,
                           lambda path,g=group,k=key,n=record['NULLITY']:self.tensor_unit(g,k,path,n),group!='RESIDUALS')
        output(prefix+'_DIMENSION_CONSTRAINTS',tuple(d.constraints),lambda path:d.zero)


class AdjointCurrentMap:
    """Power-map row representation and physical/mixed current contractions."""

    def __init__(self, modal):
        self.modal = modal
        self.modes, self.r = modal.modes, modal.r
        self.field, self.row = modal.field_units, modal.row_units
        self.energy, self.power, self.current = modal.energy_unit, modal.power_unit, modal.current_unit
        self.frequency, self.length = modal.frequency_unit, modal.length_unit

    @staticmethod
    def add(*units):return tuple(sum(v) for v in zip(*units))

    @staticmethod
    def negative(unit):return tuple(-v for v in unit)

    def prepare(self):
        m = self.modal
        m.prepare()
        b,p = m.symbolic_operands['POWER_MAP_PLUS'],m.symbolic_operands['PENCIL_PLUS']
        weighted = b*p
        operands = {'POWER_WEIGHTED_PENCIL':weighted}
        residuals = {}
        for label,variable in (('NORMAL',m.pairing.c.leg_momenta[1]),('FREQUENCY',m.pairing.frequencies[1])):
            transport = m.scalar_operands[label+'_TRANSPORT']
            direct = weighted.diff(variable)+transport*weighted.diff(m.pairing.acoustic.qlegs[1])
            product = m.symbolic_operands[label+'_POWER_MAP_PLUS']*p+b*m.symbolic_operands[label+'_PENCIL_PLUS']
            operands[label+'_DIRECT_DERIVATIVE'] = direct
            operands[label+'_PRODUCT_DERIVATIVE'] = product
            residuals[label+'_PRODUCT_RULE'] = (direct-product).applyfunc(m.pairing.carrier_expansion)
        self.symbolic_operands, self.symbolic_residuals = operands,residuals
        self.evaluate = {key:sp.lambdify(m.variables,value.xreplace(m.bindings),'numpy',cse=True)
                         for key,value in operands.items() if key.endswith('DIRECT_DERIVATIVE')}

    def construct(self, source, progress=lambda record:None):
        self.prepare()
        records = []
        add,neg = self.add,self.negative
        zero = (0,0,0)
        velocity = add(self.length,self.frequency)
        for old in source['RECORDS']:
            n = old['NULLITY']
            record = {'INDEX':old['INDEX'],'ROOT_DISK_INDEX':old['ROOT_DISK_INDEX'],
                'NORMAL_LIFT_SIGN':old['NORMAL_LIFT_SIGN'],'NULLITY':n,'K':old['K'],'Q':old['Q'],
                'NORMAL_REALITY_STATUS':old['NORMAL_REALITY_CERTIFICATE']['NORMAL_REALITY_STATUS'],
                'SHEET_MEMBERSHIP':old['SHEET_MEMBERSHIP'],
                'BULK_DECAY_DISK_CERTIFIED':old['BULK_DECAY_DISK_CERTIFIED'],
                'SOURCE_CURRENT_NORMALIZATION_DEFINED':old['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED'],
                'COEFFICIENT_FRAME':'EPSILON_SQUARED_REMOVED_FOR_ROW_MAP_AND_PAIRING',
                'NORM_FRAME':'NUMERIC_L_T_M_REFERENCE_UNIT_COEFFICIENTS','ITEMS':[]}
            items = record['ITEMS']
            def put(group,key,value,units,quadratic=False):
                value = np.asarray(value,dtype=complex)
                if value.ndim==1:value=value.reshape((-1,1))
                dims = [units(i,j) if callable(units) else units for i in range(value.shape[0]) for j in range(value.shape[1])]
                items.append({'GROUP':group,'NAME':key,'VALUE':value,'UNITS':tuple(dims),'EPSILON_POWER':2 if quadratic else 0})
            right = old['FORMS']['RIGHT']
            p,b = old['OPERANDS']['PENCIL_PLUS'],old['OPERANDS']['POWER_MAP_PLUS']
            singular = np.linalg.svd(b,compute_uv=False)
            threshold = 1e-10*max(1.,singular[0])
            rank = int(np.sum(singular>threshold))
            record.update({'POWER_MAP_RANK':rank,'POWER_MAP_RANK_THRESHOLD':threshold,
                'FREQUENCY_PAIRING_RANK':old['FREQUENCY_PAIRING_RANK'],
                'INVERTIBLE_FIELD_MAP_DEFINED':bool(rank==b.shape[0] and old['FREQUENCY_PAIRING_RANK']==n)})
            put('OPERANDS','POWER_MAP',b,lambda i,j:add(self.power,neg(self.field[i]),neg(self.row[j])))
            put('OPERANDS','POWER_MAP_SINGULAR_VALUES',singular,zero)
            put('OPERANDS','PENCIL',p,lambda i,j:add(self.row[i],neg(self.field[j])))
            put('OPERANDS','RIGHT',right,lambda i,j:self.field[i])
            if not record['INVERTIBLE_FIELD_MAP_DEFINED']:
                records.append(record);progress({'mode':old['INDEX'],'powerMapRank':rank,'mapDefined':False});continue
            left = old['FORMS']['LEFT_FREQUENCY_NORMALIZED']
            inverse = np.linalg.solve(b,np.eye(b.shape[0]))
            adjoint = np.linalg.solve(b.conj().T,left)
            weighted = b@p
            bridge = old['FORMS']['POWER_LEFT_BRIDGE_FREQUENCY_NORMALIZED']
            defect = old['FORMS']['POWER_LEFT_DEFECT']
            field_defect = np.linalg.solve(b.conj().T,defect.conj().T)
            row_left_unit = lambda i,j:add(self.frequency,neg(self.row[i]))
            field_left_unit = lambda i,j:add(self.field[i],neg(self.energy))
            q_unit = lambda i,j:add(self.power,neg(self.field[i]),neg(self.field[j]))
            put('OPERANDS','LEFT_FREQUENCY_NORMALIZED',left,row_left_unit)
            put('MAPS','INVERSE_POWER_MAP',inverse,lambda i,j:add(self.row[i],self.field[j],neg(self.power)))
            put('MAPS','ADJOINT_FIELD',adjoint,field_left_unit)
            put('MAPS','POWER_WEIGHTED_PENCIL',weighted,q_unit)
            put('MAPS','PHYSICAL_BRIDGE',bridge,self.energy)
            put('MAPS','POWER_COVECTOR_DEFECT',defect,lambda i,j:add(self.power,neg(self.row[j])))
            put('MAPS','PHYSICAL_FIELD_DEFECT',field_defect,lambda i,j:self.field[i])
            put('RESIDUALS','ROW_TO_FIELD',b.conj().T@adjoint-left,row_left_unit)
            put('RESIDUALS','INVERSE_SOLVE_JOIN',adjoint-inverse.conj().T@left,field_left_unit)
            put('RESIDUALS','ADJOINT_WEIGHTED_KERNEL',weighted.conj().T@adjoint,
                lambda i,j:add(self.frequency,neg(self.field[i])))
            put('RESIDUALS','RIGHT_WEIGHTED_KERNEL',weighted@right,lambda i,j:add(self.power,neg(self.field[i])))
            put('RESIDUALS','PHYSICAL_FIELD_RECONSTRUCTION',adjoint@bridge.conj().T+field_defect-right,
                lambda i,j:self.field[i])
            k,q = complex(old['K']),complex(old['PHYSICAL_Q'])
            omega,height = complex(old['OMEGA']),complex(old['DEPTH_CUTOFF'])
            point = (omega,omega,k.conjugate(),k,q.conjugate(),q,height)
            equal = old['EQUAL_DEPTH_BRANCH']
            excluded = {'INFINITE_DEPTH_INTEGRAL'} | ({'DEPTH_INTEGRAL','DEPTH_INTEGRAL_Q'} if equal else set())
            values = {key:np.asarray(f(*point),dtype=complex) for key,f in self.modal.evaluate.items() if key not in excluded}
            put('RESIDUALS','SOURCE_POWER_MAP_JOIN',values['POWER_MAP_PLUS']-b,
                lambda i,j:add(self.power,neg(self.field[i]),neg(self.row[j])))
            put('RESIDUALS','SOURCE_PENCIL_JOIN',values['PENCIL_PLUS']-p,lambda i,j:add(self.row[i],neg(self.field[j])))
            integral = complex(old['DEPTH_INTEGRAL'])
            integral_q = complex(values['EQUAL_DEPTH_INTEGRAL_Q' if equal else 'DEPTH_INTEGRAL_Q'])
            full_energy = values['ENERGY_SLAB']+integral*values['ENERGY_BULK']
            full_current = values['CURRENT_SLAB']+integral*values['CURRENT_BULK']
            forms = {'ENERGY_FINITE':(full_energy,self.energy,zero),
                     'CURRENT_FINITE':(full_current,self.current,velocity)}
            if old['BULK_DECAY_DISK_CERTIFIED']:
                forms['CURRENT_INFINITE'] = (values['CURRENT_SLAB']+
                    complex(old['INFINITE_DEPTH_INTEGRAL'])*values['CURRENT_BULK'],self.current,velocity)
            for name,(value,unit,mixed_unit) in forms.items():
                mixed = adjoint.conj().T@value@right
                physical = right.conj().T@value@right
                defect_part = field_defect.conj().T@value@right
                put('OPERANDS',name+'_FIELD_MATRIX',value,lambda i,j:add(unit,neg(self.field[i]),neg(self.field[j])),True)
                put('FORMS',name+'_MIXED',mixed,mixed_unit,True)
                put('FORMS',name+'_PHYSICAL',physical,unit,True)
                put('FORMS',name+'_BRIDGE_TERM',bridge@mixed,unit,True)
                put('FORMS',name+'_DEFECT_TERM',defect_part,unit,True)
                put('RESIDUALS',name+'_PHYSICAL_JOIN',physical-old['FORMS'][name],unit,True)
                put('RESIDUALS',name+'_BRIDGE_RECONSTRUCTION',bridge@mixed+defect_part-physical,unit,True)
            for label,var_unit,result_name,result_unit in (
                ('NORMAL',neg(self.length),'CURRENT_FINITE',velocity),
                ('FREQUENCY',self.frequency,'ENERGY_FINITE',zero)):
                db,dp = values[label+'_POWER_MAP_PLUS'],values[label+'_PENCIL_PLUS']
                direct = np.asarray(self.evaluate[label+'_DIRECT_DERIVATIVE'](*point),dtype=complex)
                product = db@p+b@dp
                canonical = left.conj().T@dp@right
                correction = adjoint.conj().T@db@p@right
                pairing = adjoint.conj().T@direct@right
                derivative_unit = lambda i,j:add(q_unit(i,j),neg(var_unit))
                put('OPERANDS',label+'_POWER_MAP_DERIVATIVE',db,
                    lambda i,j:add(self.power,neg(self.field[i]),neg(self.row[j]),neg(var_unit)))
                put('OPERANDS',label+'_WEIGHTED_DIRECT_DERIVATIVE',direct,derivative_unit)
                put('OPERANDS',label+'_WEIGHTED_PRODUCT_DERIVATIVE',product,derivative_unit)
                put('FORMS',label+'_CANONICAL_PAIRING',canonical,result_unit)
                put('FORMS',label+'_WEIGHTED_PAIRING',pairing,result_unit)
                put('FORMS',label+'_WEIGHT_DERIVATIVE_CORRECTION',correction,result_unit)
                put('RESIDUALS',label+'_PRODUCT_RULE',direct-product,derivative_unit)
                put('RESIDUALS',label+'_PAIRING_COVARIANCE',pairing-canonical-correction,result_unit)
                if label=='FREQUENCY':
                    put('RESIDUALS','WEIGHTED_FREQUENCY_NORMALIZATION',pairing-np.eye(n),zero)
                transport = complex(values[label+'_TRANSPORT'])
                de = values[label+'_ENERGY_SLAB']+integral*values[label+'_ENERGY_BULK']+integral_q*transport*values['ENERGY_BULK']
                dj = values[label+'_CURRENT_SLAB']+integral*values[label+'_CURRENT_BULK']+integral_q*transport*values['CURRENT_BULK']
                dtop = complex(values['TOP_FACTOR'])*(values[label+'_CURRENT_DEPTH']+
                    height*complex(values[label+'_DEPTH_WEIGHT'])*values['CURRENT_DEPTH'])
                derivative_balance = (complex(values[label+'_TIME_WEIGHT'])*full_energy+complex(values['TIME_RATE'])*de+
                    complex(values[label+'_NORMAL_WEIGHT'])*full_current+complex(values['NORMAL_RATE'])*dj+
                    values[label+'_POWER_INTERFACE']+dtop)
                source_derivative = values[label+'_POWER_SOURCE']
                isolated = complex(values[label+('_NORMAL_WEIGHT' if label=='NORMAL' else '_TIME_WEIGHT')])
                physical_matrix = full_current if label=='NORMAL' else full_energy
                other_source = adjoint.conj().T@(source_derivative-b@dp)@right
                other_balance = adjoint.conj().T@(derivative_balance-isolated*physical_matrix)@right
                mixed = adjoint.conj().T@physical_matrix@right
                reconstructed = (canonical+other_source-other_balance)/isolated
                put('FORMS',label+'_DERIVATIVE_WORK_PAIRING',adjoint.conj().T@b@dp@right,result_unit,True)
                put('FORMS',label+'_OTHER_SOURCE_WORK',other_source,result_unit,True)
                put('FORMS',label+'_OTHER_BALANCE_WORK',other_balance,result_unit,True)
                put('FORMS',label+'_RECONSTRUCTED_MIXED_FORM',reconstructed,result_unit,True)
                put('RESIDUALS',label+'_MIXED_BALANCE_DERIVATIVE',
                    adjoint.conj().T@(derivative_balance-source_derivative)@right,result_unit,True)
                put('RESIDUALS',label+'_MIXED_CURRENT_ENERGY_RECONSTRUCTION',reconstructed-mixed,result_unit,True)
            transform = np.diag(np.arange(2,n+2,dtype=float)).astype(complex)+np.triu(np.ones((n,n))*(1/3+1j/5),1)
            inverse_t = np.linalg.solve(transform,np.eye(n))
            new_right,new_left = right@transform,left@inverse_t.conj().T
            new_adjoint = np.linalg.solve(b.conj().T,new_left)
            covectors = new_right.conj().T@b
            gram = new_left.conj().T@new_left
            new_bridge = np.linalg.solve(gram.T,(covectors@new_left).T).T
            new_defect = covectors-new_bridge@new_left.conj().T
            new_mixed = new_adjoint.conj().T@full_current@new_right
            new_physical = new_right.conj().T@full_current@new_right
            put('BASIS','CHANGE',transform,zero)
            put('BASIS','RIGHT',new_right,lambda i,j:self.field[i])
            put('BASIS','LEFT_NORMALIZED',new_left,row_left_unit)
            put('BASIS','ADJOINT_FIELD',new_adjoint,field_left_unit)
            put('BASIS','BRIDGE',new_bridge,self.energy)
            put('BASIS','DEFECT',new_defect,lambda i,j:add(self.power,neg(self.row[j])))
            put('BASIS','MIXED_CURRENT',new_mixed,velocity,True)
            put('BASIS','PHYSICAL_CURRENT',new_physical,self.current,True)
            put('RESIDUALS','BASIS_FREQUENCY_PAIRING',new_left.conj().T@values['FREQUENCY_PENCIL_PLUS']@new_right-np.eye(n),zero)
            put('RESIDUALS','BASIS_FIELD_MAP',new_adjoint-adjoint@inverse_t.conj().T,field_left_unit)
            put('RESIDUALS','BASIS_BRIDGE',new_bridge-transform.conj().T@bridge@transform,self.energy)
            put('RESIDUALS','BASIS_DEFECT',new_defect-transform.conj().T@defect,lambda i,j:add(self.power,neg(self.row[j])))
            put('RESIDUALS','BASIS_MIXED_CURRENT',new_mixed-inverse_t@(adjoint.conj().T@full_current@right)@transform,velocity,True)
            put('RESIDUALS','BASIS_PHYSICAL_CURRENT',new_physical-transform.conj().T@old['FORMS']['CURRENT_FINITE']@transform,self.current,True)
            if old['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED']:
                t = old['FORMS']['FIELD_TO_FLUX_MAP']
                flux_current = t.conj().T@(bridge@(adjoint.conj().T@forms['CURRENT_INFINITE'][0]@right)+
                    field_defect.conj().T@forms['CURRENT_INFINITE'][0]@right)@t
                put('FORMS','RECONSTRUCTED_SIGNED_CURRENT',flux_current,zero,True)
                put('RESIDUALS','SIGNED_CURRENT_RECONSTRUCTION',flux_current-old['FORMS']['SIGNED_CURRENT'],zero,True)
            record['RESIDUAL_NORMS'] = {item['NAME']:float(np.linalg.norm(item['VALUE'])) for item in items if item['GROUP']=='RESIDUALS'}
            records.append(record)
            progress({'mode':old['INDEX'],'powerMapRank':rank,'mapDefined':True,
                'maximumResidualNorm':max(record['RESIDUAL_NORMS'].values())})
        return {'RECORDS':records,'SYMBOLIC_OPERANDS':self.symbolic_operands,'SYMBOLIC_RESIDUALS':self.symbolic_residuals}

    def emit(self, result, provenance, prefix='ADJOINT_CURRENT_MAP', context='REFERENCE_LAB_HELD_RHO4_CONSTANT'):
        def output(tag,body,units,heavy=False):
            body = cas(body)
            emit(tag,carrier_fingerprint(body) if heavy=='carrier' else self.modes.compact_fingerprint(body) if heavy else body)
            emit('METADATA_'+tag,self.modes.numeric_metadata(body,units))
        output(prefix+'_PROVENANCE',provenance,lambda path:(0,0,0))
        for group in ('SYMBOLIC_OPERANDS','SYMBOLIC_RESIDUALS'):
            for name,value in result[group].items():
                variable_unit = self.frequency if name.startswith('FREQUENCY_') else self.negative(self.length) if name.startswith('NORMAL_') else (0,0,0)
                output(prefix+'_'+group+'_'+name,value,lambda path:self.add(self.power,
                    self.negative(self.field[path[0]//5]),self.negative(self.field[path[0]%5]),self.negative(variable_unit)),
                    'carrier' if group=='SYMBOLIC_OPERANDS' else False)
        for record in result['RECORDS']:
            tag = prefix+'_'+context+'_'+str(record['INDEX'])
            info = {k:v for k,v in record.items() if k!='ITEMS'}
            output(tag+'_RECORD',info,lambda path:self.negative(self.length) if path[-1]=='K' else self.frequency if path[-1]=='Q' else (0,0,0))
            for item in record['ITEMS']:
                array = item['VALUE']
                body = sp.ImmutableMatrix(*array.shape,[self.modes.number(v) for v in array.ravel()])*self.modal.epsilon**item['EPSILON_POWER']
                output(tag+'_'+item['GROUP']+'_'+item['NAME'],body,lambda path:item['UNITS'][path[0]],item['GROUP']!='RESIDUALS')
        output(prefix+'_DIMENSION_CONSTRAINTS',tuple(PHYSICAL_METADATA.dimensions.constraints),lambda path:(0,0,0))


class TwoEndedMatchingChannels:
    """Source-driven end bases and complete cross-mode current matrices.

    The block assembly below is a coordinate map. Physical currents between
    different roots are evaluated from the polarized source, never supplied
    by a block-diagonal ansatz. No interior matching or S-matrix is solved here.
    """

    def __init__(self, modal_builder, modal, adjoint, end):
        self.builder, self.modal, self.adjoint = modal_builder, modal, adjoint
        self.end = end
        self.orientation = {'LEFT': -1, 'RIGHT': 1}[end]

    def construct(self):
        b = self.builder
        b.prepare()
        candidates, eligible = [], []
        native = self.modal['NATIVE_RECORDS']
        adjoints = {r['INDEX']: r for r in self.adjoint['RECORDS']}
        for record in self.modal['RECORDS']:
            index, n = record['INDEX'], record['NULLITY']
            old, adjoint = native[index], adjoints[index]
            info = {key: record[key] for key in ('INDEX', 'ROOT_DISK_INDEX', 'NORMAL_LIFT_SIGN',
                'K', 'Q', 'PHYSICAL_Q', 'OMEGA', 'NULLITY', 'SHEET_MEMBERSHIP',
                'EXACT_REAL_NORMAL', 'BULK_DECAY_DISK_CERTIFIED',
                'PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED')}
            info.update({key: old[key] for key in
                         ('CLASSIFIER_DEFINED', 'CLASSIFIER_WEIGHTS', 'CLASSIFIER_STATUS')})
            info.update({'OUTWARD_END_ORIENTATION': self.orientation,
                'OUTWARD_IMAGINARY_K_SIGN_NUMERIC': int(np.sign(self.orientation*complex(record['K']).imag)),
                'ADJOINT_FIELD_MAP_DEFINED': adjoint['INVERTIBLE_FIELD_MAP_DEFINED'],
                'RIGHT_RANK_MINUS_NULLITY': int(np.linalg.matrix_rank(record['FORMS']['RIGHT'], tol=1e-9))-n,
                'LEFT_RANK_MINUS_NULLITY': int(np.linalg.matrix_rank(record['FORMS']['LEFT'], tol=1e-9))-n,
                'NATIVE_NULLITY_RESIDUAL': n-int(old['NULLITY'])})
            bases = {key: record['FORMS'][key] for key in
                     ('RIGHT', 'LEFT', 'LEFT_FREQUENCY_NORMALIZED', 'FLUX_RIGHT', 'FLUX_LEFT')
                     if key in record['FORMS']}
            items = [item for item in adjoint['ITEMS'] if
                     (item['GROUP'], item['NAME']) == ('MAPS', 'ADJOINT_FIELD')]
            candidates.append({'INFO': info, 'BASES': bases, 'ADJOINT_FIELD_ITEMS': items})
            if record['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED']:
                eligible.append(record)
        sizes = [r['NULLITY'] for r in eligible]
        offsets = np.cumsum([0]+sizes)
        total = int(offsets[-1])
        field_current = np.zeros((total, total), dtype=complex)
        flux_current = np.zeros_like(field_current)
        coordinate_map = np.zeros_like(field_current)
        source_signed = np.zeros_like(field_current)
        pairs, channels = [], []
        for i, left in enumerate(eligible):
            a = slice(offsets[i], offsets[i+1])
            coordinate_map[a, a] = left['FORMS']['FIELD_TO_FLUX_MAP']
            source_signed[a, a] = left['FORMS']['SIGNED_CURRENT']
            for column in range(left['NULLITY']):
                current = float(np.real(left['FORMS']['SIGNED_CURRENT'][column, column]))
                outward = self.orientation*current
                channels.append({'END': self.end, 'ROOT_DISK_INDEX': left['ROOT_DISK_INDEX'],
                    'NORMAL_LIFT_SIGN': left['NORMAL_LIFT_SIGN'], 'RECORD_INDEX': left['INDEX'],
                    'BASIS_COLUMN': column, 'MATRIX_COLUMN': int(offsets[i]+column),
                    'DIRECTION': 'INCOMING' if outward < 0 else 'OUTGOING' if outward > 0 else 'UNRESOLVED'})
            for j, right in enumerate(eligible):
                c = slice(offsets[j], offsets[j+1])
                qleft, qright = complex(left['PHYSICAL_Q']), complex(right['PHYSICAL_Q'])
                if not (qleft.imag > 0 and qright.imag > 0):
                    raise ValueError('cross-mode infinite-depth current requires the recorded decay domain')
                point = (complex(left['OMEGA']), complex(right['OMEGA']), complex(left['K']).conjugate(),
                         complex(right['K']), qleft.conjugate(), qright, complex(right['DEPTH_CUTOFF']))
                slab = np.asarray(b.evaluate['CURRENT_SLAB'](*point), dtype=complex)
                bulk = np.asarray(b.evaluate['CURRENT_BULK'](*point), dtype=complex)
                depth = complex(b.evaluate['INFINITE_DEPTH_INTEGRAL'](*point))
                depth_rate = complex(b.evaluate['DEPTH_RATE'](*point))
                current = slab+depth*bulk
                if not np.isfinite(current).all() or not np.isfinite(depth):
                    raise ValueError('non-finite cross-mode current operand')
                field_current[a, c] = left['FORMS']['RIGHT'].conj().T@current@right['FORMS']['RIGHT']
                flux_current[a, c] = left['FORMS']['FLUX_RIGHT'].conj().T@current@right['FORMS']['FLUX_RIGHT']
                pairs.append({'LEFT_RECORD_INDEX': left['INDEX'], 'RIGHT_RECORD_INDEX': right['INDEX'],
                    'DEPTH_INTEGRAL': depth, 'DEPTH_RATE': depth_rate, 'CURRENT_SLAB': slab,
                    'CURRENT_BULK': bulk, 'COMPOSED_CURRENT': current})
        transformed = coordinate_map.conj().T@field_current@coordinate_map
        diagonal = np.zeros_like(flux_current)
        for i in range(len(eligible)):
            block = slice(offsets[i], offsets[i+1])
            diagonal[block, block] = flux_current[block, block]
        return {'CANDIDATES': candidates, 'CHANNELS': channels, 'PAIR_OPERANDS': pairs,
            'FIELD_CURRENT': field_current, 'FLUX_CURRENT': flux_current,
            'FIELD_TO_FLUX_MAP': coordinate_map, 'TRANSFORMED_CURRENT': transformed,
            'OUTWARD_CURRENT': self.orientation*flux_current,
            'SOURCE_SIGNED_CURRENT_BLOCKS': source_signed,
            'INTER_ROOT_CURRENT': flux_current-diagonal,
            'RESIDUALS': {'BASIS_CHANGE': flux_current-transformed,
                'SOURCE_DIAGONAL_BLOCKS': diagonal-source_signed,
                'HERMITIAN': flux_current-flux_current.conj().T},
            'COUNTS': {'CANDIDATES': len(candidates), 'BASIS_DIRECTIONS': sum(r['INFO']['NULLITY'] for r in candidates),
                'OPEN_RECORDS': len(eligible), 'OPEN_BASIS_DIRECTIONS': total,
                'INCOMING': sum(r['DIRECTION']=='INCOMING' for r in channels),
                'OUTGOING': sum(r['DIRECTION']=='OUTGOING' for r in channels),
                'EVALUATED_RECORD_PAIRS': len(pairs)}}

    def emit(self, result, provenance, prefix):
        b, modes = self.builder, self.builder.modes
        zero = PHYSICAL_METADATA.dimensions.zero
        def converted(value):
            if isinstance(value, np.ndarray):
                return sp.ImmutableMatrix(*value.shape, [modes.number(v) for v in value.ravel()])
            if isinstance(value, (complex, np.complexfloating)):
                return modes.number(value)
            if isinstance(value, dict):
                return {k: converted(v) for k, v in value.items()}
            if isinstance(value, (tuple, list)):
                return tuple(converted(v) for v in value)
            return value
        def output(name, value, unit=lambda p: zero, quadratic=False, heavy=False):
            body = cas(converted(value))
            if quadratic:
                body = b.epsilon**2*body
            emit(prefix+'_'+name, modes.compact_fingerprint(body) if heavy else body)
            emit('METADATA_'+prefix+'_'+name, modes.numeric_metadata(body, unit))
        output('PROVENANCE', provenance)
        output('COUNTS', result['COUNTS'])
        output('CHANNELS', result['CHANNELS'])
        for candidate in result['CANDIDATES']:
            info = candidate['INFO']; tag = 'CANDIDATE_'+str(info['INDEX'])
            output(tag+'_INFO', info, lambda p: b.info_unit(p, info))
            for key, value in candidate['BASES'].items():
                output(tag+'_'+key, value, lambda p, k=key, n=info['NULLITY']: b.tensor_unit('FORMS', k, p, n), heavy=True)
            for item in candidate['ADJOINT_FIELD_ITEMS']:
                output(tag+'_ADJOINT_FIELD', item['VALUE'], lambda p: item['UNITS'][p[0]], heavy=True)
        for i, pair in enumerate(result['PAIR_OPERANDS']):
            tag = 'PAIR_'+str(i)
            output(tag+'_DOMAIN', {k: v for k, v in pair.items() if not isinstance(v, np.ndarray)},
                lambda p: b.length_unit if p[-1]=='DEPTH_INTEGRAL' else
                tuple(-v for v in b.length_unit) if p[-1]=='DEPTH_RATE' else zero)
            for key in ('CURRENT_SLAB', 'CURRENT_BULK', 'COMPOSED_CURRENT'):
                base_unit = b.power_unit if key=='CURRENT_BULK' else b.current_unit
                output(tag+'_'+key, pair[key],
                    lambda p, u=base_unit: tuple(u[j]-b.field_units[p[0]//5][j]-b.field_units[p[0]%5][j] for j in range(3)),
                    quadratic=True, heavy=True)
        for key in ('FIELD_CURRENT', 'FLUX_CURRENT', 'TRANSFORMED_CURRENT', 'OUTWARD_CURRENT',
                    'SOURCE_SIGNED_CURRENT_BLOCKS', 'INTER_ROOT_CURRENT'):
            output(key, result[key], lambda p, k=key: b.current_unit if k=='FIELD_CURRENT' else zero,
                   quadratic=True, heavy=True)
        output('FIELD_TO_FLUX_MAP', result['FIELD_TO_FLUX_MAP'],
               lambda p: tuple(-v/2 for v in b.current_unit), heavy=True)
        for key, value in result['RESIDUALS'].items():
            output('RESIDUAL_'+key, value, quadratic=True)


class ChannelInput:
    """Explicit profile and numerical parameter input, separate from PIT.

    Parameter values are coefficients in the declared L/T/M reference frame.
    This supplies an instance of the interface class, never a class-wide
    spectral assertion. Every required algebraic carrier must be bound.
    """

    def __init__(self, reduction, specification):
        self.r = reduction
        self.specification = specification
        self.parameters = {name: sp.Rational(value)
                           for name, value in specification['parameters'].items()}
        self.frame = tuple(specification['unit_frame'])
        if len(self.frame) != 3 or len(set(self.frame)) != 3:
            raise ValueError('channel input needs three distinct L/T/M reference-unit names')
        self.profiles = {p: sp.sympify(specification['profiles'][p],
                                      locals={'xi': reduction.xi}) for p in ('w', 'm')}
        if any(f.free_symbols-set((reduction.xi,)) for f in self.profiles.values()):
            raise ValueError('channel profiles contain unbound parameters')
        self.limits = {key: sp.limit(self.profiles[p], reduction.xi, end)
                       for (p, end), key in reduction.end_values.items()}
        if any(value.is_finite is not True or value.is_real is not True
               for value in self.limits.values()):
            raise ValueError('channel input profile limits are unresolved or non-finite')
        eta = self.parameters['eta_bg']
        self.origin = {reduction.symbols['eta_bg']: eta,
                       reduction.symbols['sigma_W']:
                           eta*self.parameters['W_0']/self.parameters['L_W']}
        if self.parameters['L_W'] <= 0 or self.parameters['W_0'] <= 0:
            raise ValueError('channel input length scales must be positive')
        if self.parameters['omega'] <= 0:
            raise ValueError('the implemented channel chart requires positive real frequency')

    def mapping(self, algebraic, relation, live):
        atoms = (algebraic.free_symbols | relation.free_symbols)-set(live)
        missing = sorted(a.name for a in atoms if a.name not in self.parameters)
        if missing:
            raise ValueError(('unbound channel parameters', missing))
        mapping = {a: self.parameters[a.name] for a in atoms}
        for limit in algebraic.atoms(sp.Limit):
            mapping[limit] = self.limits[limit]
        return mapping

    def emit_inputs(self):
        r = self.r
        emit('CHANNEL_INPUT_SPECIFICATION', self.specification)
        emit('CHANNEL_INPUT_SHA256', hashlib.sha256(json.dumps(
            self.specification, sort_keys=True, separators=(',', ':')).encode()).hexdigest())
        emit('CHANNEL_INPUT_UNIT_FRAME', self.frame)
        parameter_symbols = dict(r.symbols)
        parameter_symbols.update({a.name: a for a in r.tangents})
        emit('CHANNEL_INPUT_PARAMETER_VALUES', [
            (name, {'CARRIER': parameter_symbols[name], 'COEFFICIENT': value,
                    'DIMENSION_L_T_M': PHYSICAL_METADATA.dimensions.measure(parameter_symbols[name]),
                    'MULTIGRADE': sorted(PHYSICAL_METADATA.coefficients(value))})
            for name, value in self.parameters.items()])
        physical('CHANNEL_INPUT_BACKGROUND_ORIGIN', tuple(self.origin.values()),
                 zero_dimensions={(i,): (0, 0, 0) for i in range(len(self.origin))})
        for name, expression in self.profiles.items():
            left = self.limits[r.end_values[(name, -sp.oo)]]
            right = self.limits[r.end_values[(name, sp.oo)]]
            derivative = sp.diff(expression, r.xi)
            zero_transfer = sp.integrate(derivative, (r.xi, -sp.oo, sp.oo))
            jump = right-left
            physical('CHANNEL_PROFILE_'+name.upper(), (expression, left, right),
                     zero_dimensions={(i,): (0, 0, 0) for i in range(3)})
            physical('CHANNEL_PROFILE_MOMENT_OPERANDS_'+name.upper(),
                     (derivative, zero_transfer, jump),
                     zero_dimensions={(i,): (0, 0, 0) for i in range(3)})
            physical('CHANNEL_PROFILE_MOMENT_RESIDUAL_'+name.upper(), zero_transfer-jump,
                     zero_dimensions={(): (0, 0, 0)})


class BulkSheetPath:
    """Fixed positive-real-frequency chart attached to the real Fourier line.

    Transport from Re(k) to k along the vertical segment. Branch points and
    segment clearance are computed from the actual reduced radical relation.
    A segment hitting a branch point is unresolved; no root is reselected at
    complex momentum. This chart does not supply complex-frequency pole paths.
    """

    def __init__(self, relation, momentum, radical):
        self.k, self.q = momentum, radical
        self.square = sp.solve(relation, radical**2)[0]
        coordinate = sp.Dummy('s11cdSheetComplexMomentum')
        polynomial = sp.Poly(self.square.xreplace({momentum:coordinate}), coordinate)
        if polynomial.degree() != 2 or any(c.is_real is not True for c in polynomial.all_coeffs()):
            raise NotImplementedError('sheet transport requires the real quadratic bulk radical')
        self.points = tuple(complex(p.evalf(30)) for p in polynomial.all_roots())
        self.evaluate = sp.lambdify(momentum,self.square,'numpy')
        self.cache = {}

    def transport(self, target):
        if target in self.cache:
            return self.cache[target]
        start = complex(target.real)
        delta = target-start
        closest = [start+(min(1.,max(0.,((p-start)*delta.conjugate()).real/abs(delta)**2))*delta
                         if delta else 0) for p in self.points]
        clearance = min(abs(a-b) for a,b in zip(self.points,closest))
        scale = max(abs(target),*(abs(p) for p in self.points))
        result = {'START_K':start,'END_K':target,'BRANCH_POINTS':self.points,
                  'MINIMUM_BRANCH_POINT_DISTANCE':clearance,
                  'GEOMETRIC_RESOLUTION':256*np.finfo(float).eps*scale,
                  'PATH_DEFINED':False}
        # Real-axis value of the already joined positive-frequency operand.
        # This is the only seed selection; subsequent roots are transported.
        seed = complex(self.evaluate(start))**0.5
        result['SEED_Q'] = seed
        if clearance <= result['GEOMETRIC_RESOLUTION'] or abs(seed) == 0:
            result['STATUS'] = 'BRANCH_LOCUS_ON_PATH'
            self.cache[target] = result
            return result
        routes = []
        for fraction in (0.25,0.125):
            u, root, count, residual = 0.,seed,0,abs(seed**2-self.evaluate(start))
            while u < 1 and count < 4096:
                position = start+u*delta
                distance = min(abs(position-p) for p in self.points)
                step = min(1-u,fraction*distance/abs(delta)) if delta else 1.
                if u+step == u:
                    break
                u = min(1.,u+step)
                position = start+u*delta
                candidate = complex(self.evaluate(position))**0.5
                root = min((candidate,-candidate),key=lambda value:abs(value-root))
                residual = max(residual,abs(root**2-self.evaluate(position)))
                count += 1
            routes.append({'STEP_FRACTION':fraction,'STEPS':count,'PARAMETER_REACHED':u,
                           'END_Q':root,'MAXIMUM_RADICAL_RESIDUAL':residual})
        result['REFINEMENTS'] = routes
        result['REFINEMENT_DIFFERENCE'] = routes[0]['END_Q']-routes[1]['END_Q']
        result['PATH_DEFINED'] = all(r['PARAMETER_REACHED'] == 1 for r in routes)
        result['STATUS'] = 'TRANSPORTED' if result['PATH_DEFINED'] else 'PATH_RESOLUTION_LIMIT'
        self.cache[target] = result
        return result

    def classify(self, momentum, radical):
        record = dict(self.transport(momentum))
        if not record['PATH_DEFINED']:
            return Str('UNRESOLVED'),record
        endpoint = record['REFINEMENTS'][-1]['END_Q']
        scale = max(abs(endpoint),abs(radical))
        record['SHEET_DIFFERENCE'] = radical-endpoint
        record['OPPOSITE_SHEET_DIFFERENCE'] = radical+endpoint
        record['RELATIVE_REFINEMENT_DIFFERENCE'] = abs(record['REFINEMENT_DIFFERENCE'])/scale
        record['RELATIVE_SHEET_DIFFERENCE'] = abs(radical-endpoint)/scale
        record['RELATIVE_OPPOSITE_SHEET_DIFFERENCE'] = abs(radical+endpoint)/scale
        tolerance = max(1e-9,256*np.finfo(float).eps*(1+abs(momentum)/record['MINIMUM_BRANCH_POINT_DISTANCE']))
        record['RELATIVE_MATCH_TOLERANCE'] = tolerance
        record['MATCH_RESOLVED'] = (tolerance < 1e-3 and
                                   record['RELATIVE_REFINEMENT_DIFFERENCE'] <= tolerance)
        same = record['RELATIVE_SHEET_DIFFERENCE'] <= tolerance
        opposite = record['RELATIVE_OPPOSITE_SHEET_DIFFERENCE'] <= tolerance
        if record['MATCH_RESOLVED'] and same != opposite:
            return bool(same),record
        record['STATUS'] = 'ROOT_MATCH_UNRESOLVED'
        return Str('UNRESOLVED'),record


class JointBulkSheetPath:
    """Lift explicit paths on the computed quadratic bulk radical curve.

    A real-axis seed is evaluated from a reduced branch operand. Root matching
    and implicit-ODE integration then transport that seed without a decay test.
    Frequency-cut encounters are recorded only in the fixed-real-k chart where
    S11b specifies them; joint paths retain their explicit history instead.
    """

    def __init__(self, relation, frequency, momentum, radical, real_axis_seed):
        self.w,self.k,self.q = frequency,momentum,radical
        self.relation = relation
        self.square = sp.solve(relation,radical**2)[0]
        polynomial = sp.Poly(self.square,frequency,momentum)
        if polynomial.total_degree()!=2 or polynomial.as_expr().free_symbols-{frequency,momentum}:
            raise NotImplementedError('joint transport requires a bound quadratic radical')
        self.coefficients = [(powers,complex(value)) for powers,value in polynomial.terms()]
        self.evaluate = sp.lambdify((frequency,momentum),self.square,'numpy')
        self.seed_expression = real_axis_seed
        self.seed = sp.lambdify((frequency,momentum),real_axis_seed,
                               [{'sqrt':np.lib.scimath.sqrt},'numpy'])
        self.dw,self.dk = sp.symbols('s11cdJointPathFrequencyTangent s11cdJointPathMomentumTangent')
        derivative = sp.Symbol('s11cdJointPathRadicalDerivative')
        self.derivative = sp.solve(sp.diff(relation,frequency)*self.dw+
            sp.diff(relation,momentum)*self.dk+sp.diff(relation,radical)*derivative,derivative)[0]
        self.ode = sp.lambdify((frequency,momentum,radical,self.dw,self.dk),self.derivative,'numpy')

    def segment(self, start, end):
        # Substitute the affine path ansatz into the computed polynomial.
        from numpy.polynomial import Polynomial
        w,k = (Polynomial((a,b-a)) for a,b in zip(start,end))
        polynomial = sum((coefficient*w**i*k**j for (i,j),coefficient in self.coefficients),Polynomial([0j]))
        roots = tuple(complex(v) for v in polynomial.roots())
        clearance = min((abs(z-complex(min(1.,max(0.,z.real)))) for z in roots),default=float('inf'))
        tolerance = 512*np.finfo(float).eps
        hit = any(-tolerance<=z.real<=1+tolerance and abs(z.imag)<=tolerance*max(1.,abs(z)) for z in roots)
        cut_records=[]
        fixed_real_k = start[1]==end[1] and start[1].imag==0
        if fixed_real_k:
            branch_points = np.polynomial.Polynomial([
                sum(c*start[1]**j for (i,j),c in self.coefficients if i==power)
                for power in range(3)]).roots()
            dw=end[0]-start[0]
            for point in branch_points:
                if abs(point.imag)>tolerance*max(1.,abs(point)):continue
                if dw.real:
                    t=(point.real-start[0].real)/dw.real
                    if 0<=t<=1 and (start[0]+t*dw).imag<0:
                        cut_records.append({'PARAMETER':t,'OMEGA':start[0]+t*dw,'BRANCH_OMEGA':complex(point)})
                elif start[0].real==point.real and min(start[0].imag,end[0].imag)<0:
                    cut_records.append({'PARAMETER_INTERVAL':(0.,1.),'BRANCH_OMEGA':complex(point)})
        return polynomial,{'ZERO_PARAMETERS':roots,'PARAMETER_CLEARANCE':clearance,
            'BRANCH_INTERSECTION':hit,'FIXED_REAL_K_CUT_CHART':fixed_real_k,
            'DOWNWARD_FREQUENCY_CUT_ENCOUNTERS':cut_records}

    def trace(self, vertices, seed=None):
        from scipy.integrate import solve_ivp
        vertices=tuple(tuple(complex(x) for x in v) for v in vertices)
        result={'VERTICES':[{'OMEGA':w,'K':k} for w,k in vertices],
                'PATH_DEFINED':False,'SEED_FROM_REAL_AXIS':seed is None,
                'FLOAT_MANTISSA_BITS':np.finfo(float).nmant+1,
                'GEOMETRIC_PARAMETER_TOLERANCE':512*np.finfo(float).eps}
        if len(vertices)<2 or not all(np.isfinite(z) for v in vertices for z in v):
            result['STATUS']='INVALID_PATH_VERTICES'
            return result
        if seed is None and any(z.imag for z in vertices[0]):
            result['STATUS']='NO_REAL_AXIS_SEED'
            return result
        seed=complex(self.seed(*vertices[0])) if seed is None else complex(seed)
        result['SEED_Q']=seed
        result['SEED_RADICAL_RESIDUAL']=seed**2-self.evaluate(*vertices[0])
        segments=[self.segment(a,b) for a,b in zip(vertices,vertices[1:])]
        result['SEGMENTS']=[data for _,data in segments]
        if seed==0 or any(data['BRANCH_INTERSECTION'] for _,data in segments):
            result['STATUS']='BRANCH_LOCUS_ON_PATH'
            return result
        result['SEED_RELATIVE_RESIDUAL']=abs(result['SEED_RADICAL_RESIDUAL'])/max(abs(seed)**2,abs(self.evaluate(*vertices[0])))
        if result['SEED_RELATIVE_RESIDUAL']>result['GEOMETRIC_PARAMETER_TOLERANCE']:
            result['STATUS']='SEED_RELATION_UNRESOLVED'
            return result
        refinements=[]
        for fraction in (0.25,0.125):
            root=seed;count=0;maximum=abs(result['SEED_RADICAL_RESIDUAL']);angle=0.;completed=0
            for polynomial,data in segments:
                t=0.;previous=complex(polynomial(t))
                while t<1 and count<16384:
                    distance=min((abs(t-z) for z in data['ZERO_PARAMETERS']),default=float('inf'))
                    step=min(1-t,fraction*distance)
                    if t+step==t:break
                    t=min(1.,t+step)
                    value=complex(polynomial(t))
                    candidate=value**0.5
                    root=min((candidate,-candidate),key=lambda q:abs(q-root))
                    maximum=max(maximum,abs(root**2-value))
                    angle+=np.angle(value/previous)
                    previous=value;count+=1
                completed+=int(t==1.)
                if t<1:break
            refinements.append({'STEP_FRACTION':fraction,'STEPS':count,'COMPLETED_SEGMENTS':completed,
                'END_Q':root,'MAXIMUM_RADICAL_RESIDUAL':maximum,'RADICAND_ARGUMENT_CHANGE':angle})
        result['REFINEMENTS']=refinements
        result['REFINEMENT_DIFFERENCE']=refinements[-1]['END_Q']-refinements[0]['END_Q']
        if not all(v['COMPLETED_SEGMENTS']==len(segments) for v in refinements):
            result['STATUS']='PATH_RESOLUTION_LIMIT'
            return result
        root=seed;ode_maximum=abs(result['SEED_RADICAL_RESIDUAL']);ode_steps=0;ode_status=[]
        result['ODE_RELATIVE_TOLERANCE']=3e-12
        result['ODE_ABSOLUTE_TOLERANCE']=3e-14*max(abs(seed),np.finfo(float).tiny)
        for start,end in zip(vertices,vertices[1:]):
            dw,dk=(b-a for a,b in zip(start,end))
            solution=solve_ivp(lambda t,q:np.asarray([self.ode(start[0]+t*dw,start[1]+t*dk,q[0],dw,dk)]),
                (0.,1.),np.asarray([root]),rtol=result['ODE_RELATIVE_TOLERANCE'],atol=result['ODE_ABSOLUTE_TOLERANCE'])
            ode_status.append(int(solution.status));ode_steps+=len(solution.t)
            ode_maximum=max(ode_maximum,max(abs(q*q-self.evaluate(start[0]+t*dw,start[1]+t*dk))
                for t,q in zip(solution.t,solution.y[0])))
            root=solution.y[0,-1]
            if solution.status!=0:break
        result.update({'ODE_END_Q':complex(root),'ODE_STEPS':ode_steps,'ODE_STATUSES':ode_status,
            'ODE_MAXIMUM_RADICAL_RESIDUAL':ode_maximum,
            'ODE_DIFFERENCE':complex(root)-refinements[-1]['END_Q'],
            'END_Q':refinements[-1]['END_Q'],
            'CLOSED_COORDINATE_PATH':vertices[0]==vertices[-1],
            'END_TO_SEED_RATIO':refinements[-1]['END_Q']/seed,
            'RADICAND_ARGUMENT_TURNS':refinements[-1]['RADICAND_ARGUMENT_CHANGE']/(2*np.pi)})
        result['PATH_DEFINED']=len(ode_status)==len(segments) and all(v==0 for v in ode_status)
        result['STATUS']='TRANSPORTED' if result['PATH_DEFINED'] else 'ODE_RESOLUTION_LIMIT'
        return result


class RectangularModeJets:
    """Implicit invariant-pair ansatz in the retained two-grade rectangle.

    Matrix normal momenta retain a degenerate mode cluster without choosing
    separate, possibly incompatible eigenvectors for its two perturbations.
    All pencil Taylor coefficients include the implicit bulk-root chain rule.
    """

    grades = ((1, 0), (0, 1), (1, 1))

    @staticmethod
    def polynomial(series, eta, sigma, origin):
        terms = [sp.ImmutableMatrix(value)*(eta-origin[eta])**a*(sigma-origin[sigma])**b
                 for (a,b),value in series.items()]
        return sum(terms,sp.zeros(*terms[0].shape))

    def __init__(self, pencil, relation, k, q, eta, sigma, origin):
        self.k, self.q, self.parameters = k, q, (eta, sigma)
        dummy = sp.Dummy('s11cdImplicitRootDerivative')
        self.root_derivatives = {v: sp.solve(sp.diff(relation, v)
            + sp.diff(relation, q)*dummy, dummy)[0] for v in (k, eta, sigma)}
        def derivative(value, coordinate):
            return value.diff(coordinate)+value.diff(q)*self.root_derivatives[coordinate]
        indices = ((0,0,0), (1,0,0), (0,1,0), (1,1,0),
                   (0,0,1), (1,0,1), (0,1,1), (0,0,2))
        self.pencil_coefficients, self.radical_coefficients = {}, {}
        for index in indices:
            a,b,j = index
            operator, radical = pencil, q
            for coordinate, count in ((eta,a), (sigma,b), (k,j)):
                for _ in range(count):
                    operator = derivative(operator, coordinate)
                    radical = derivative(radical, coordinate)
            self.pencil_coefficients[index] = (operator/sp.factorial(j)).subs(origin).applyfunc(sp.cancel)
            self.radical_coefficients[index] = sp.cancel((radical/sp.factorial(j)).subs(origin))
        self.evaluate = sp.lambdify((k,q), tuple(self.pencil_coefficients.values()), 'numpy', cse=True)
        self.evaluate_radical = sp.lambdify((k,q), tuple(self.radical_coefficients.values()), 'numpy', cse=True)

    @staticmethod
    def multiply(a, b):
        result = {}
        for (i,j), av in a.items():
            for (m,n), bv in b.items():
                index = (i+m,j+n)
                if max(index) <= 1:
                    value = av@bv
                    result[index] = result.get(index, np.zeros_like(value))+value
        return result

    @classmethod
    def equation(cls, coefficients, modes, shift):
        n = modes[(0,0)].shape[1]
        powers = [{(0,0):np.eye(n,dtype=complex)}, shift]
        powers.append(cls.multiply(shift,shift))
        result = {}
        for (a,b,j), matrix in coefficients.items():
            for (c,d), value in cls.multiply(modes,powers[j]).items():
                index = (a+c,b+d)
                if max(index) <= 1:
                    term = matrix@value
                    result[index] = result.get(index,np.zeros_like(term))+term
        return result

    @classmethod
    def pair(cls, coefficients, basis):
        rows, n = basis.shape
        count = (rows+n)*n
        zero_mode = np.zeros_like(basis)
        zero_root = np.zeros((n,n),dtype=complex)
        seed = {(0,0):basis}
        linear_coefficients = {index:coefficients[index] for index in ((0,0,0),(0,0,1))}
        # Differentiate the invariant-pair/gauge ansatz by coefficient probes.
        # Its coefficient at this grade is exactly linear in these unknowns.
        columns = []
        for index in range(count):
            direction = np.zeros(count,dtype=complex)
            direction[index] = 1
            dr, dk = direction[:rows*n].reshape(rows,n), direction[rows*n:].reshape(n,n)
            equation = cls.equation(linear_coefficients, {**seed,(1,0):dr}, {(1,0):dk})[(1,0)]
            columns.append(np.concatenate((equation.ravel(),(basis.conj().T@dr).ravel())))
        jacobian = np.column_stack(columns)
        singular = np.linalg.svd(jacobian,compute_uv=False)
        threshold = np.finfo(float).eps*max(jacobian.shape)*singular[0]
        rank = int(np.count_nonzero(singular>threshold))
        diagnostics = {'JACOBIAN_SINGULAR_VALUES':singular.tolist(), 'JACOBIAN_RANK':rank,
                       'JACOBIAN_RANK_THRESHOLD':threshold, 'UNKNOWN_COUNT':count}
        if rank != count:
            diagnostics['STATUS'] = 'SINGULAR_INVARIANT_PAIR_JACOBIAN'
            return None,None,diagnostics,jacobian
        modes, shift = dict(seed), {}
        residuals = {}
        for grade in cls.grades:
            forcing = cls.equation(coefficients,modes,shift).get(grade,zero_mode)
            rhs = -np.concatenate((forcing.ravel(),zero_root.ravel()))
            solved = np.linalg.solve(jacobian,rhs)
            modes[grade] = solved[:rows*n].reshape(rows,n)
            shift[grade] = solved[rows*n:].reshape(n,n)
            equation = cls.equation(coefficients,modes,shift).get(grade,zero_mode)
            gauge = basis.conj().T@modes[grade]
            residuals[grade] = {'EQUATION':equation, 'GAUGE':gauge,
                               'LINEAR_SYSTEM':jacobian@solved-rhs}
        diagnostics['STATUS'] = 'COMPUTED_INVARIANT_PAIR'
        diagnostics['BASE_EQUATION'] = coefficients[(0,0,0)]@basis
        diagnostics['COEFFICIENT_RESIDUALS'] = residuals
        diagnostics['ROOT_COEFFICIENT_COMMUTATOR'] = shift[(1,0)]@shift[(0,1)]-shift[(0,1)]@shift[(1,0)]
        return modes,shift,diagnostics,jacobian

    def construct(self, k, q, right, left):
        coefficients = {index:np.asarray(value,dtype=complex) for index,value in
                        zip(self.pencil_coefficients,self.evaluate(k,q))}
        r,kr,rd,rj = self.pair(coefficients,right)
        l,kl,ld,lj = self.pair({i:v.conj().T for i,v in coefficients.items()},left)
        result = {'RIGHT_DIAGNOSTICS':rd,'LEFT_DIAGNOSTICS':ld,
                  'RIGHT_JACOBIAN':rj,'LEFT_JACOBIAN':lj,
                  'DEFINED':r is not None and l is not None}
        if not result['DEFINED']:
            return result
        n = right.shape[1]
        identity = np.eye(n,dtype=complex)
        radical = {index:complex(value)*identity for index,value in
                   zip(self.radical_coefficients,self.evaluate_radical(k,q))}
        qseries = self.equation(radical,{(0,0):identity},kr)
        result.update(RIGHT=r,LEFT=l,K={(0,0):k*identity,**kr},Q=qseries)
        overlap = self.multiply({i:v.conj().T for i,v in l.items()},r)
        result['OVERLAP_SINGULAR_VALUES'] = np.linalg.svd(overlap[(0,0)],compute_uv=False).tolist()
        if np.linalg.matrix_rank(overlap[(0,0)],tol=1e-9) == n:
            inverse = {(0,0):np.linalg.inv(overlap[(0,0)])}
            for grade in self.grades:
                known = self.multiply(overlap,inverse).get(grade,np.zeros((n,n),complex))
                inverse[grade] = -inverse[(0,0)]@known
            projector = self.multiply(self.multiply(r,inverse),{i:v.conj().T for i,v in l.items()})
            result['PROJECTOR'] = projector
            square = self.multiply(projector,projector)
            result['PROJECTOR_RESIDUAL'] = {i:square[i]-v for i,v in projector.items()}
            result['INVERSE_OVERLAP_RESIDUAL'] = self.multiply(overlap,inverse)
            result['INVERSE_OVERLAP_RESIDUAL'][(0,0)] -= identity
        return result


class FullPencilModes:
    """Carrier-first spectral PIT on the positive-frequency retarded chart.

    Scalars are sampled before polynomial elimination, while normal momentum
    and the inherited bulk radical remain live algebraic variables. The solve
    is on the complete quotient pencil. A zero coupling block is never an
    input. All candidate roots, including rejected sheet/pole/chart roots, are
    retained with their measured residuals.
    """

    def __init__(self, ends, curl, units):
        self.ends, self.r, self.k = ends, ends.r, ends.kn
        self.curl, self.units = curl, units
        self.q = sp.Symbol('s11cdBulkRadical', complex=True)
        self.eta, self.sigma = (self.r.symbols[n] for n in ('eta_bg', 'sigma_W'))
        self.indices = (1, 2, 3, 4, 5)
        self.embedding = sp.ImmutableMatrix(sp.eye(6)[:, self.indices])
        self.unitless = PHYSICAL_METADATA.dimensions.zero
        self.closed_operands = None
        PHYSICAL_METADATA.dimensions.known[self.q] = PHYSICAL_METADATA.dimensions.measure(self.r.omega)

    def set_closed_operands(self, strong, lift, current, current_builder):
        self.closed_operands = (strong, lift*self.embedding,
                                current['SLAB_CURRENT_MATRIX'], current_builder.leg_momenta,
                                current_builder.field_units,
                                PHYSICAL_METADATA.dimensions.measure(current['SLAB_CURRENT']))

    def numeric_metadata(self, value, unit):
        """Units inherited from the solved equation; grades from coefficients."""
        metadata = {}
        for path, expression in leaves(cas(value)):
            if isinstance(expression, Str):
                continue
            if isinstance(expression, sp.logic.boolalg.BooleanAtom):
                expression = sp.Integer(bool(expression))
            coefficients = PHYSICAL_METADATA.coefficients(expression)
            homotopy = {}
            for (e, a, b), c in coefficients.items():
                homotopy[(e, a+b)] = homotopy.get((e, a+b), sp.S.Zero)+c*(self.r.symbols['W_0']/self.r.ell)**b
            descriptor = (tuple(unit(path)),tuple(sorted(coefficients)),
                          tuple(sorted(g for g,c in homotopy.items() if c != 0)))
            metadata.setdefault(descriptor,[]).append(path)
        return cas([{'PATHS':paths,'DIMENSION_L_T_M':d,'MULTIGRADE':g,'EPSILON_LAMBDA_SUPPORT':h}
                    for (d,g,h),paths in metadata.items()])

    def compact_fingerprint(self, value):
        """Whole-tensor SHA and scalar PIT projections in the numeric unit frame.

        Position-dependent weights fingerprint the coefficient tensor after
        the carrier-first spectral solve. No eigenproblem is replaced by a
        random carrier. Literal equation residuals are emitted separately.
        """
        body = cas(value)
        digest = hashlib.sha256(sp.srepr(body).encode()).hexdigest()
        @lru_cache(maxsize=None)
        def evaluated(expression,index):
            if expression.is_number:
                return complex(expression)
            return complex(expression.subs({self.eta:sp.Rational(index+1,101),
                                            self.sigma:sp.Rational(index+2,103),
                                            self.r.symbols['epsilon_shape']:sp.Rational(index+3,107)}))
        samples = []
        for index in range(3):
            total = 0j
            for path, expression in leaves(body):
                if isinstance(expression,Str):
                    continue
                seed = hashlib.sha256((str(path)+':'+str(index)).encode()).digest()
                weight = (int.from_bytes(seed[:2],'big')%97+1)/101
                total += weight*evaluated(expression,index)
            samples.append(self.number(total))
        return cas({'OBJECT_SHA256':digest,'NUMERIC_UNIT_FRAME_TENSOR_PIT':samples})

    def quotient(self, full, suffix):
        # Axial A0=0 chart, valid at nonzero first tangential momentum.
        # Other charts are emitted; the PIT samples lie in this chart.
        charts = []
        for axis in range(3):
            keep = tuple(i for i in range(6) if i != axis)
            b = sp.ImmutableMatrix(sp.eye(6)[:, keep])
            gauge = self.curl.nullspace()[0]
            gauge = (gauge*sp.lcm([sp.denom(e) for e in gauge])).applyfunc(sp.cancel)
            g = sp.ImmutableMatrix.vstack(gauge, sp.zeros(3, 1))
            denominator = g[axis]
            projection = sp.eye(6)-g*sp.eye(6)[axis, :]/denominator
            projection = projection.applyfunc(sp.cancel)
            restriction = b.T*projection
            charts.append((axis, denominator, b, restriction,
                           (restriction*b-sp.eye(5)).applyfunc(sp.cancel),
                           (b*restriction-projection).applyfunc(sp.cancel)))
        fingerprinted('GAUGE_QUOTIENT_CHARTS_'+suffix, cas(charts),
                      {p: self.unitless for p, _ in leaves(cas(charts))})
        gauge = sp.ImmutableMatrix.vstack(self.curl.nullspace()[0], sp.zeros(3, 1))
        gauge_residual = (full*gauge).applyfunc(sp.expand)
        adjoint_residual = (gauge.T*full).applyfunc(sp.expand)
        fingerprinted('FULL_PENCIL_GAUGE_RIGHT_RESIDUAL_'+suffix, gauge_residual,
                      {(i,): self.units[(6*i,)] for i in range(6)})
        fingerprinted('FULL_PENCIL_GAUGE_LEFT_RESIDUAL_'+suffix, adjoint_residual,
                      {(j,): self.units[(j,)] for j in range(6)})
        quotient = self.embedding.T*full*self.embedding
        quotient_units = {(5*i+j,): self.units[(6*ri+cj,)]
                          for i, ri in enumerate(self.indices) for j, cj in enumerate(self.indices)}
        fingerprinted('GAUGE_QUOTIENT_PENCIL_'+suffix, quotient, quotient_units)
        sector = sp.diag(1, 1, 0, 0, 0)
        commutator = quotient*sector-sector*quotient
        fingerprinted('REFERENCE_SECTOR_COMMUTATOR_'+suffix, commutator, quotient_units)
        return quotient, quotient_units

    def analytic(self, quotient):
        """Continue the actual positive-frequency real-axis branch operands."""
        positive = quotient.xreplace({sp.sign(self.r.omega): sp.S.One})
        piecewise = sorted(positive.atoms(sp.Piecewise), key=sp.default_sort_key)
        differences = tuple(sp.cancel(p.args[0].expr-p.args[-1].expr) for p in piecewise)
        # Both real-axis branches must join this analytic expression; the
        # branch residual is emitted before the structural gate in solve().
        positive = positive.xreplace({p: p.args[0].expr for p in piecewise})
        radical_powers = [p for p in positive.atoms(sp.Pow)
                          if p.exp.is_Rational and p.exp.q == 2 and p.has(self.k)]
        radicals = sorted({sp.sqrt(p.base) for p in radical_powers}, key=sp.default_sort_key)
        if len(radicals) != 1:
            raise NotImplementedError(('spectral radical census', radicals))
        radical = radicals[0]
        relation = self.q**2-radical.base
        replacements = {p: self.q**(2*p.exp) for p in radical_powers}
        algebraic = positive.xreplace(replacements)
        return algebraic, relation, differences

    @staticmethod
    def number(value):
        z = complex(value)
        return sp.Float(z.real, 17)+sp.I*sp.Float(z.imag, 17)

    def sample(self, algebraic, relation, index):
        atoms = sorted((algebraic.free_symbols | relation.free_symbols)-{self.k, self.q, self.eta, self.sigma},
                       key=sp.default_sort_key)
        atoms += sorted(algebraic.atoms(sp.Limit), key=sp.default_sort_key)
        mapping = {}
        for a in atoms:
            seed = hashlib.sha256((sp.srepr(a)+':mode:'+str(index)).encode()).digest()
            mapping[a] = sp.Rational(int.from_bytes(seed[:2], 'big') % 13+2,
                                     int.from_bytes(seed[2:4], 'big') % 11+3)
        # eta/sigma remain differentiation variables; sample points are at
        # their origin for the controlled first-background-grade mode jets.
        return mapping

    @staticmethod
    def rational_determinant(matrix):
        """Clear row denominators before the small exact determinant."""
        denominators, rows = [], []
        for i in range(matrix.rows):
            entries = [sp.cancel(matrix[i, j]) for j in range(matrix.cols)]
            denominator = sp.lcm([sp.denom(e) for e in entries])
            rows.append([sp.cancel(e*denominator) for e in entries])
            denominators.append(denominator)
        polynomial_matrix = sp.ImmutableMatrix(rows)
        determinant = polynomial_matrix.det(method='domain-ge')
        rational = sp.cancel(determinant/sp.prod(denominators))
        return sp.fraction(rational), polynomial_matrix, tuple(denominators)

    def solve(self, full, suffix, end_sign):
        quotient, quotient_units = self.quotient(full, suffix)
        algebraic, relation, differences = self.analytic(quotient)
        physical('SPECTRAL_REAL_AXIS_BRANCH_JOIN_RESIDUAL_'+suffix, differences,
                 zero_dimensions={(i,): self.unitless for i in range(len(differences))})
        fingerprinted('SPECTRAL_ALGEBRAIC_PENCIL_'+suffix, algebraic, quotient_units)
        physical('SPECTRAL_RADICAL_RELATION_'+suffix, relation)
        if any(d != 0 for d in differences):
            raise NotImplementedError('unjoined real-axis bulk sheets')
        for sample_index in range(3):
            self.solve_sample(algebraic, relation, quotient_units, suffix, end_sign, sample_index)

    def solve_input(self, full, suffix, end_sign, channel_input, *, reference=False):
        name = 'INPUT_'+suffix
        quotient, quotient_units = self.quotient(full, name)
        algebraic, relation, differences = self.analytic(quotient)
        physical('SPECTRAL_INPUT_BRANCH_JOIN_RESIDUAL_'+suffix, differences,
                 zero_dimensions={(i,): self.unitless for i in range(len(differences))})
        if any(d != 0 for d in differences):
            raise NotImplementedError('unjoined real-axis bulk sheets')
        mapping = channel_input.mapping(algebraic, relation,
                                         (self.k, self.q, self.eta, self.sigma))
        origin = ({self.eta: sp.S.Zero, self.sigma: sp.S.Zero}
                  if reference else channel_input.origin)
        return self.solve_sample(algebraic, relation, quotient_units, name,
                                 end_sign, 0, carrier_values=mapping,
                                 grade_origin=origin,
                                 record_kind='INPUT_REFERENCE' if reference else 'INPUT_END')

    def solve_sample(self, algebraic, relation, quotient_units, suffix, end_sign, sample_index,
                     *, carrier_values=None, grade_origin=None, record_kind='PIT'):
        name = suffix+'_'+str(sample_index)
        mapping = self.sample(algebraic, relation, sample_index) if carrier_values is None else carrier_values
        origin = {self.eta: sp.S.Zero, self.sigma: sp.S.Zero} if grade_origin is None else grade_origin
        sampled = algebraic.xreplace(mapping)
        relation_sample = relation.xreplace(mapping)
        sheet_path = BulkSheetPath(relation_sample,self.k,self.q)
        k_squared = sp.solve(relation_sample, self.k**2)[0]
        # The determinant is even in k on this input. Polynomial division
        # verifies this algebraic elimination rather than discarding odd terms.
        baseline = sampled.subs(origin).applyfunc(sp.cancel)
        (numerator, denominator), cleared, row_denominators = self.rational_determinant(baseline)
        relation_k = sp.Poly(self.k**2-k_squared, self.k)
        eliminated = sp.rem(sp.Poly(numerator, self.k), relation_k).as_expr()
        if eliminated.has(self.k):
            raise NotImplementedError(('odd spectral determinant', name))
        polynomial = sp.Poly(eliminated, self.q)
        factors = sp.sqf_list(polynomial)[1]
        roots = [(complex(root), int(multiplicity)) for factor, multiplicity in factors
                 for root in sp.nroots(factor, n=30, maxsteps=300)]
        emit('SPECTRAL_'+record_kind+'_CARRIER_VALUES_'+name, mapping)
        physical('SPECTRAL_'+record_kind+'_LOCAL_GRADE_ORIGIN_'+name, tuple(origin.values()),
                 zero_dimensions={(i,): self.unitless for i in range(len(origin))})
        emit('METADATA_SPECTRAL_'+record_kind+'_CARRIER_VALUES_'+name,
             cas([(str(k), {'DIMENSION_L_T_M': PHYSICAL_METADATA.dimensions.measure(k),
                             'MULTIGRADE': tuple(PHYSICAL_METADATA.coefficients(v))}) for k,v in mapping.items()]))
        emit('SPECTRAL_'+record_kind+'_ELIMINATION_'+name, carrier_fingerprint(cas((cleared, row_denominators,
                                                               numerator, denominator, polynomial.as_expr()))))
        # All coefficients in this diagnostic have been expressed in the
        # numerical unit frame. Physical k/q and row/column units are emitted
        # below, separately from the dimensionless polynomial coordinates.
        elimination = cas((cleared, row_denominators,numerator,denominator,polynomial.as_expr()))
        emit('METADATA_SPECTRAL_'+record_kind+'_ELIMINATION_'+name,
             self.numeric_metadata(elimination,lambda p:self.unitless))
        # Numeric arrays use the restored units of each original column/row;
        # they are coefficients, not dimensionless physical modes by fiat.
        labels = ('A0', 'A1', 'A2', 'Theta', 'E', 'Phi')
        dimensions = PHYSICAL_METADATA.dimensions
        column_units = [dimensions.known[sp.Function('s11cdReducedTrial'+labels[j])] for j in self.indices]
        row_units = [tuple(a+b for a, b in zip(quotient_units[(5*i,)], column_units[0])) for i in range(5)]
        emit('SPECTRAL_'+record_kind+'_RESTORED_UNITS_'+name, cas({'K': dimensions.measure(self.k),
             'Q': dimensions.measure(self.q), 'RIGHT_COLUMNS': column_units,
             'LEFT_COLUMNS': [tuple(-x for x in d) for d in row_units],
             'PENCIL': quotient_units}))
        qk = sp.solve(sp.diff(relation_sample, self.k)+sp.diff(relation_sample, self.q)*sp.Symbol('s11cdqk'),
                      sp.Symbol('s11cdqk'))[0]
        lk = baseline.diff(self.k)+baseline.diff(self.q)*qk
        derivatives = [sampled.diff(g).subs(origin) for g in (self.eta, self.sigma)]
        evaluate = sp.lambdify((self.k, self.q), baseline, 'numpy', cse=True)
        evaluate_k = sp.lambdify((self.k, self.q), lk, 'numpy', cse=True)
        evaluate_grades = [sp.lambdify((self.k, self.q), d, 'numpy', cse=True) for d in derivatives]
        qw = sp.solve(sp.diff(relation, self.r.omega)+sp.diff(relation, self.q)*sp.Symbol('s11cdqw'),
                      sp.Symbol('s11cdqw'))[0]
        lw = (algebraic.diff(self.r.omega)+algebraic.diff(self.q)*qw).xreplace(mapping).subs(origin)
        evaluate_w = sp.lambdify((self.k, self.q), lw, 'numpy', cse=True)
        # A second construction evaluates the original pencil at neighboring
        # frequencies. Continue the radical locally from the measured root;
        # no differentiated expression is an operand of this difference check.
        fixed_except_omega = {a: v for a, v in mapping.items() if a != self.r.omega}
        frequency_family = algebraic.xreplace(fixed_except_omega).subs(origin)
        evaluate_frequency = sp.lambdify((self.r.omega, self.k, self.q),
                                         frequency_family, 'numpy', cse=True)
        radical_squared = sp.solve(relation, self.q**2)[0].xreplace(fixed_except_omega)
        evaluate_radical_squared = sp.lambdify((self.r.omega, self.k), radical_squared, 'numpy')
        rectangular_jets = RectangularModeJets(sampled,relation_sample,self.k,self.q,
                                               self.eta,self.sigma,origin)
        derivative_indices = tuple(rectangular_jets.pencil_coefficients)
        derivative_values = cas(tuple(rectangular_jets.pencil_coefficients.values()))
        emit('RECTANGULAR_PENCIL_TAYLOR_INDICES_'+record_kind+'_'+name,derivative_indices)
        emit('METADATA_RECTANGULAR_PENCIL_TAYLOR_INDICES_'+record_kind+'_'+name,
             self.numeric_metadata(cas(derivative_indices),lambda p:self.unitless))
        emit('RECTANGULAR_PENCIL_TAYLOR_OPERANDS_'+record_kind+'_'+name,
             carrier_fingerprint(derivative_values))
        emit('METADATA_RECTANGULAR_PENCIL_TAYLOR_OPERANDS_'+record_kind+'_'+name,
             self.numeric_metadata(derivative_values,lambda p:tuple(a-derivative_indices[p[0]][2]*b
                 for a,b in zip(quotient_units[(p[1],)],dimensions.measure(self.k)))))
        radical_values = cas(tuple(rectangular_jets.radical_coefficients.values()))
        emit('RECTANGULAR_RADICAL_TAYLOR_OPERANDS_'+record_kind+'_'+name,radical_values)
        emit('METADATA_RECTANGULAR_RADICAL_TAYLOR_OPERANDS_'+record_kind+'_'+name,
             self.numeric_metadata(radical_values,lambda p:tuple(a-derivative_indices[p[0]][2]*b
                 for a,b in zip(dimensions.measure(self.q),dimensions.measure(self.k)))))
        field_evaluators = None
        if carrier_values is not None and self.closed_operands is not None:
            strong, lift, slab_current, current_legs, field_units, current_unit = self.closed_operands
            strong_algebraic, strong_relation, strong_join = self.analytic(strong)
            physical('INPUT_CLOSED_FIELD_BRANCH_JOIN_RESIDUAL_'+name, strong_join,
                     zero_dimensions={(i,):self.unitless for i in range(len(strong_join))})
            physical('INPUT_CLOSED_FIELD_RADICAL_RELATION_RESIDUAL_'+name,
                     sp.expand(strong_relation-relation),
                     zero_dimensions={():tuple(2*v for v in dimensions.measure(self.q))})
            evaluate_strong = sp.lambdify((self.k,self.q),
                strong_algebraic.xreplace(mapping).subs(origin), 'numpy', cse=True)
            evaluate_lift = sp.lambdify(self.k, lift.xreplace(mapping), 'numpy', cse=True)
            current_at_input = slab_current.xreplace(mapping).subs(origin)
            field_evaluators = (evaluate_strong, evaluate_lift, current_at_input,
                                current_legs, field_units)
        omega_point = float(mapping[self.r.omega])
        kval = sp.lambdify(self.q, k_squared, 'numpy')
        denominator_value = sp.lambdify((self.k, self.q), denominator, 'numpy')
        outputs = []
        for qroot, multiplicity in roots:
            for kroot in (complex(kval(qroot))**0.5, -complex(kval(qroot))**0.5):
                matrix = np.asarray(evaluate(kroot, qroot), dtype=complex)
                if not np.isfinite(matrix).all():
                    outputs.append(cas({'K': self.number(kroot), 'Q': self.number(qroot),
                                        'FINITE_PENCIL': False}))
                    continue
                left, singular, right_h = np.linalg.svd(matrix)
                tol = 1e-8*max(1., singular[0])
                nullity = int(np.sum(singular < tol))
                physical_sheet, sheet_record = sheet_path.classify(kroot,qroot)
                record = {'K': self.number(kroot), 'Q': self.number(qroot), 'MULTIPLICITY': multiplicity,
                          'SINGULAR_VALUES': list(map(self.number, singular)), 'NULLITY': nullity,
                          'PHYSICAL_BULK_SHEET': physical_sheet,
                          'BULK_SHEET_PATH': sheet_record,
                          'DENOMINATOR': self.number(denominator_value(kroot, qroot)),
                          'RADICAL_RESIDUAL': self.number(complex(relation_sample.subs({self.k:kroot,self.q:qroot}))),
                          'DECAY_AT_END': bool(end_sign*kroot.imag > 1e-10),
                          'REAL_NORMAL_MOMENTUM': bool(abs(kroot.imag) <= 1e-10),
                          'NORMAL_THRESHOLD': bool(abs(kroot) <= 1e-10),
                          'BULK_THRESHOLD': bool(abs(qroot) <= 1e-10)}
                ktotal = complex(kroot**2+sum(mapping[t]**2 for t in self.r.tangents))
                record['HELMHOLTZ_CHART_OPERAND'] = self.number(ktotal)
                record['HELMHOLTZ_CHART_DOMAIN_LIMITATION'] = bool(abs(ktotal) < 1e-9)
                if nullity:
                    right = right_h.conj().T[:, -nullity:]
                    dual = left[:, -nullity:]
                    derivative = np.asarray(evaluate_k(kroot, qroot), dtype=complex)
                    pairing = dual.conj().T@derivative@right
                    omega_derivative = np.asarray(evaluate_w(kroot, qroot), dtype=complex)
                    omega_pairing = dual.conj().T@omega_derivative@right
                    omega_rank = np.linalg.matrix_rank(omega_pairing, tol=1e-9)
                    mode_id = name+'_'+str(len(outputs))
                    if field_evaluators is not None:
                        evaluate_strong, evaluate_lift, current_at_input, current_legs, field_units = field_evaluators
                        lifted = np.asarray(evaluate_lift(kroot), dtype=complex)@right
                        closed_residual = np.asarray(evaluate_strong(kroot,qroot), dtype=complex)@lifted
                        current_matrix = current_at_input.subs({current_legs[0]:kroot.conjugate(),
                                                               current_legs[1]:kroot})
                        physical_fields = sp.ImmutableMatrix(lifted)
                        slab_pairing = physical_fields.conjugate().T*current_matrix*physical_fields
                        slab_pairing = slab_pairing.applyfunc(sp.expand)
                        record['CLOSED_PHYSICAL_FIELD_RESIDUAL'] = sp.ImmutableMatrix(closed_residual)
                        record['CONSERVATIVE_SLAB_CURRENT_PAIRING'] = slab_pairing
                        emit('INPUT_CLOSED_PHYSICAL_FIELDS_'+mode_id, self.compact_fingerprint(physical_fields))
                        emit('METADATA_INPUT_CLOSED_PHYSICAL_FIELDS_'+mode_id,
                             self.numeric_metadata(physical_fields, lambda p:field_units[p[0]//nullity]))
                        emit('INPUT_CONSERVATIVE_SLAB_CURRENT_OPERAND_'+mode_id,
                             self.compact_fingerprint(current_matrix))
                        emit('METADATA_INPUT_CONSERVATIVE_SLAB_CURRENT_OPERAND_'+mode_id,
                             self.numeric_metadata(current_matrix, lambda p:tuple(a-b-c for a,b,c in
                                zip(current_unit,field_units[p[0]//5],field_units[p[0]%5]))))
                    overlap = dual.conj().T@right
                    record.update({'RIGHT_RESIDUAL': sp.ImmutableMatrix(matrix@right),
                                   'LEFT_RESIDUAL': sp.ImmutableMatrix(matrix.conj().T@dual),
                                   'K_DERIVATIVE_PAIRING': sp.ImmutableMatrix(pairing),
                                   'OMEGA_DERIVATIVE_PAIRING': sp.ImmutableMatrix(omega_pairing),
                                   'OMEGA_PAIRING_SINGULAR_VALUES': list(map(self.number,
                                       np.linalg.svd(omega_pairing, compute_uv=False))),
                                   'OMEGA_PAIRING_RANK': int(omega_rank),
                                   'OMEGA_NORMALIZATION_DEFINED': bool(omega_rank == nullity),
                                   'OVERLAP_SINGULAR_VALUES': list(map(self.number, np.linalg.svd(overlap, compute_uv=False)))})
                    difference_pairs, difference_operands = [], []
                    for divisor in (1000, 2000):
                        step = abs(omega_point)/divisor
                        perturbed = []
                        radicals = []
                        for delta in (-step, step):
                            candidate = complex(evaluate_radical_squared(omega_point+delta, kroot))**0.5
                            continued = min((candidate, -candidate), key=lambda value: abs(value-qroot))
                            radicals.append(continued)
                            perturbed.append(np.asarray(evaluate_frequency(
                                omega_point+delta, kroot, continued), dtype=complex))
                        difference_pair = dual.conj().T@((perturbed[1]-perturbed[0])/(2*step))@right
                        difference_pairs.append((self.number(step), sp.ImmutableMatrix(difference_pair),
                                                 sp.ImmutableMatrix(difference_pair-omega_pairing)))
                        difference_operands.append((self.number(step),
                            tuple(map(self.number, radicals)),
                            tuple(sp.ImmutableMatrix(value) for value in perturbed)))
                    record['OMEGA_PAIRING_DIFFERENCE_QUOTIENTS'] = difference_pairs
                    emit('OMEGA_DIFFERENCE_PENCIL_OPERANDS_'+mode_id,
                         self.compact_fingerprint(difference_operands))
                    emit('METADATA_OMEGA_DIFFERENCE_PENCIL_OPERANDS_'+mode_id,
                         self.numeric_metadata(difference_operands, lambda p:
                             dimensions.measure(self.r.omega) if p[1] < 2 else
                             quotient_units[(p[3],)]))
                    if omega_rank == nullity:
                        # Normalize the whole nullspace. The omega derivative
                        # includes the computed bulk-radical chain rule above;
                        # neither a group velocity nor an Euclidean overlap
                        # substitutes for the nonlinear-pencil pairing.
                        normalized_dual = dual@np.linalg.inv(omega_pairing).conj().T
                        normalized_pairing = normalized_dual.conj().T@omega_derivative@right
                        record.update({
                            'OMEGA_NORMALIZED_PAIRING': sp.ImmutableMatrix(normalized_pairing),
                            'OMEGA_NORMALIZATION_RESIDUAL': sp.ImmutableMatrix(
                                normalized_pairing-np.eye(nullity)),
                            'OMEGA_NORMALIZED_LEFT_RESIDUAL': sp.ImmutableMatrix(
                                matrix.conj().T@normalized_dual),
                        })
                        normalized_left_units = [tuple(w-v for w, v in zip(
                            dimensions.measure(self.r.omega), d)) for d in row_units]
                        normalized_left = sp.ImmutableMatrix(normalized_dual)
                        emit('LEFT_MODE_OMEGA_NORMALIZED_'+mode_id,
                             self.compact_fingerprint(normalized_left))
                        emit('METADATA_LEFT_MODE_OMEGA_NORMALIZED_'+mode_id,
                             self.numeric_metadata(normalized_left,
                                 lambda p: normalized_left_units[p[0]//nullity]))
                    if np.linalg.matrix_rank(overlap, tol=1e-9) == nullity:
                        projector = right@np.linalg.solve(overlap, dual.conj().T)
                        # These weights are computed from the full nullspace
                        # projector, not from a guessed channel label.
                        tweight = np.trace(projector[:2,:2])/nullity
                        hweight = np.trace(projector[2:,2:])/nullity
                        mixed_degeneracy = abs(tweight) > 1e-7 and abs(hweight) > 1e-7
                        record.update({'PROJECTOR_RESIDUAL': sp.ImmutableMatrix(projector@projector-projector),
                                       'CLASSIFIER_WEIGHTS': (self.number(tweight), self.number(hweight)),
                                       'CLASSIFIER_DOMAIN_LIMITATION': bool(mixed_degeneracy),
                                       'CLASSIFICATION': 'UNRESOLVED_MIXED_SUBSPACE' if mixed_degeneracy else
                                           ('TRANSVERSE_LIKE' if abs(tweight) > abs(hweight) else 'THICKNESS_LIKE')})
                    else:
                        record['CLASSIFIER_DOMAIN_LIMITATION'] = sp.S.true
                    jets = []
                    jet_defined = np.linalg.matrix_rank(pairing, tol=1e-9) == nullity
                    right_polynomial = sp.ImmutableMatrix(right)
                    left_polynomial = sp.ImmutableMatrix(dual)
                    k_polynomial = sp.ImmutableMatrix(np.eye(nullity)*kroot)
                    projector_polynomial = None
                    if np.linalg.matrix_rank(overlap, tol=1e-9) == nullity:
                        projector_polynomial = sp.ImmutableMatrix(projector)
                    for grade, evaluate_grade in zip((self.eta, self.sigma), evaluate_grades):
                        perturbation = np.asarray(evaluate_grade(kroot, qroot), dtype=complex)
                        if np.linalg.matrix_rank(pairing, tol=1e-9) != nullity:
                            jets.append((grade, Str('SINGULAR_IMPLICIT_MODE_JACOBIAN')))
                            continue
                        shifts, vectors = np.linalg.eig(-np.linalg.solve(pairing, dual.conj().T@perturbation@right))
                        shift_matrix = -np.linalg.solve(pairing, dual.conj().T@perturbation@right)
                        right_forcing = perturbation@right+derivative@right@shift_matrix
                        right_correction = np.linalg.lstsq(np.vstack((matrix, right.conj().T)),
                            np.vstack((-right_forcing, np.zeros((nullity,nullity)))), rcond=None)[0]
                        left_shift = -np.linalg.solve(pairing.conj().T, right.conj().T@perturbation.conj().T@dual)
                        left_forcing = perturbation.conj().T@dual+derivative.conj().T@dual@left_shift
                        left_correction = np.linalg.lstsq(np.vstack((matrix.conj().T, dual.conj().T)),
                            np.vstack((-left_forcing, np.zeros((nullity,nullity)))), rcond=None)[0]
                        local_coordinate = grade-origin[grade]
                        right_polynomial += local_coordinate*sp.ImmutableMatrix(right_correction)
                        left_polynomial += local_coordinate*sp.ImmutableMatrix(left_correction)
                        k_polynomial += local_coordinate*sp.ImmutableMatrix(shift_matrix)
                        if projector_polynomial is not None:
                            inverse_overlap = np.linalg.inv(overlap)
                            overlap_correction = left_correction.conj().T@right+dual.conj().T@right_correction
                            projector_correction = (right_correction@inverse_overlap@dual.conj().T
                                +right@inverse_overlap@left_correction.conj().T
                                -right@inverse_overlap@overlap_correction@inverse_overlap@dual.conj().T)
                            projector_polynomial += local_coordinate*sp.ImmutableMatrix(projector_correction)
                        jet_modes = []
                        for j, shift in enumerate(shifts):
                            vector = right@vectors[:,j]
                            forcing = (perturbation+derivative*shift)@vector
                            augmented = np.vstack((matrix, right.conj().T))
                            correction = np.linalg.lstsq(augmented, np.concatenate((-forcing, np.zeros(nullity))), rcond=None)[0]
                            jet_modes.append((self.number(shift),sp.ImmutableMatrix(matrix@correction+forcing)))
                        jets.append((grade, jet_modes, sp.ImmutableMatrix(matrix.conj().T@left_correction+left_forcing)))
                    record['FIRST_GRADE_IMPLICIT_MODES' if grade_origin is None else
                           'LOCAL_GRADE_IMPLICIT_MODES'] = jets
                    record['IMPLICIT_MODE_JET_DEFINED'] = bool(jet_defined)
                    rectangle = rectangular_jets.construct(kroot,qroot,right,dual)
                    record['RECTANGULAR_MODE_JET_DEFINED'] = rectangle['DEFINED']
                    record['RECTANGULAR_JACOBIAN_COEFFICIENT_DIAGNOSTICS'] = {
                        side:{key:value for key,value in rectangle[side+'_DIAGNOSTICS'].items()
                              if key in ('JACOBIAN_SINGULAR_VALUES','JACOBIAN_RANK',
                                         'JACOBIAN_RANK_THRESHOLD','UNKNOWN_COUNT','STATUS')}
                        for side in ('RIGHT','LEFT')}
                    for side in ('RIGHT','LEFT'):
                        jacobian = sp.ImmutableMatrix(rectangle[side+'_JACOBIAN'])
                        tag = 'RECTANGULAR_'+side+'_JACOBIAN_COEFFICIENTS_'+mode_id
                        emit(tag,self.compact_fingerprint(jacobian))
                        # Coefficients in the declared row/column unit frame;
                        # physical mode/pencil residuals carry restored units below.
                        emit('METADATA_'+tag,self.numeric_metadata(jacobian,lambda p:self.unitless))
                    if rectangle['DEFINED']:
                        old_polynomials = (right_polynomial,left_polynomial,k_polynomial)
                        right_polynomial = rectangular_jets.polynomial(rectangle['RIGHT'],self.eta,self.sigma,origin)
                        left_polynomial = rectangular_jets.polynomial(rectangle['LEFT'],self.eta,self.sigma,origin)
                        k_polynomial = rectangular_jets.polynomial(rectangle['K'],self.eta,self.sigma,origin)
                        q_polynomial = rectangular_jets.polynomial(rectangle['Q'],self.eta,self.sigma,origin)
                        record['RECTANGULAR_EQUATION_COEFFICIENT_RESIDUALS'] = {
                            side:{'G'+''.join(map(str,g)):sp.ImmutableMatrix(value['EQUATION'])
                                  for g,value in rectangle[side+'_DIAGNOSTICS']['COEFFICIENT_RESIDUALS'].items()}
                            for side in ('RIGHT','LEFT')}
                        record['RECTANGULAR_GAUGE_COEFFICIENT_RESIDUALS'] = {
                            side:{'G'+''.join(map(str,g)):sp.ImmutableMatrix(value['GAUGE'])
                                  for g,value in rectangle[side+'_DIAGNOSTICS']['COEFFICIENT_RESIDUALS'].items()}
                            for side in ('RIGHT','LEFT')}
                        record['RECTANGULAR_ROOT_COEFFICIENT_COMMUTATORS'] = tuple(sp.ImmutableMatrix(
                            rectangle[side+'_DIAGNOSTICS']['ROOT_COEFFICIENT_COMMUTATOR']) for side in ('RIGHT','LEFT'))
                        emit('BULK_RADICAL_MODE_JET_'+mode_id,self.compact_fingerprint(q_polynomial))
                        emit('METADATA_BULK_RADICAL_MODE_JET_'+mode_id,
                             self.numeric_metadata(q_polynomial,lambda p:dimensions.measure(self.q)))
                        for label,old,new,unit_fn in (
                            ('RIGHT',old_polynomials[0],right_polynomial,lambda p:column_units[p[0]//nullity]),
                            ('LEFT',old_polynomials[1],left_polynomial,lambda p:tuple(-v for v in row_units[p[0]//nullity])),
                            ('ROOT',old_polynomials[2],k_polynomial,lambda p:dimensions.measure(self.k))):
                            first = new-(self.eta-origin[self.eta])*(self.sigma-origin[self.sigma])*new.diff(
                                self.eta,self.sigma).subs(origin)
                            for kind,value in (('LEGACY_OPERAND',old),('RECTANGULAR_OPERAND',first),('RESIDUAL',first-old)):
                                tag='FIRST_JET_ROUTE_'+label+'_'+kind+'_'+mode_id
                                emit(tag,self.compact_fingerprint(value))
                                emit('METADATA_'+tag,self.numeric_metadata(value,unit_fn))
                        if 'PROJECTOR' in rectangle:
                            projector_polynomial = rectangular_jets.polynomial(rectangle['PROJECTOR'],self.eta,self.sigma,origin)
                            projector_residual = rectangular_jets.polynomial(rectangle['PROJECTOR_RESIDUAL'],self.eta,self.sigma,origin)
                            tag='CLASSIFIER_PROJECTOR_RECTANGLE_RESIDUAL_'+mode_id
                            emit(tag,self.compact_fingerprint(projector_residual))
                            emit('METADATA_'+tag,self.numeric_metadata(projector_residual,lambda p:
                                tuple(a-b for a,b in zip(column_units[p[0]//5],column_units[p[0]%5]))))
                        jet_defined = True
                    record['FIRST_GRADE_PROJECTED_JET_DEFINED'] = record['IMPLICIT_MODE_JET_DEFINED']
                    record['IMPLICIT_MODE_JET_DEFINED'] = bool(jet_defined)
                    if np.linalg.matrix_rank(pairing, tol=1e-9) == nullity:
                        frequency_slopes = np.linalg.eigvals(-np.linalg.solve(pairing, omega_pairing))
                        record['RETARDED_K_FREQUENCY_SLOPES'] = list(map(self.number,frequency_slopes))
                        record['INCOMING'] = [bool(physical_sheet is True and abs(kroot.imag)<1e-10 and end_sign*x.real<0)
                                              for x in frequency_slopes]
                        record['OUTGOING'] = [bool(physical_sheet is True and abs(kroot.imag)<1e-10 and end_sign*x.real>0)
                                              for x in frequency_slopes]
                    coefficient_label = ('JET' if grade_origin is None else 'LOCAL_JET') if jet_defined else 'BASE_COEFFICIENT'
                    for label, value, unit_fn in (
                        ('RIGHT_MODE_'+coefficient_label, right_polynomial, lambda p: column_units[p[0]//nullity]),
                        ('LEFT_MODE_'+coefficient_label, left_polynomial, lambda p: tuple(-v for v in row_units[p[0]//nullity])),
                        ('NORMAL_ROOT_'+coefficient_label, k_polynomial, lambda p: dimensions.measure(self.k))):
                        emit(label+'_'+mode_id, self.compact_fingerprint(value))
                        emit('METADATA_'+label+'_'+mode_id, self.numeric_metadata(value,unit_fn))
                    if projector_polynomial is not None:
                        emit('CLASSIFIER_PROJECTOR_'+coefficient_label+'_'+mode_id, self.compact_fingerprint(projector_polynomial))
                        emit('METADATA_CLASSIFIER_PROJECTOR_'+coefficient_label+'_'+mode_id,
                             self.numeric_metadata(projector_polynomial,
                                 lambda p: tuple(a-b for a,b in zip(column_units[p[0]//5],column_units[p[0]%5]))))
                outputs.append(cas(record))
        # Numeric mode coefficients and literal residuals are intentionally
        # bounded; the symbolic input SHA and exact carrier values accompany them.
        emit('FULL_PENCIL_MODE_'+record_kind+'_'+name, outputs)
        def record_unit(path):
            key = path[1] if len(path)>1 else ''
            if key in ('K',): return dimensions.measure(self.k)
            if key == 'Q': return dimensions.measure(self.q)
            if key == 'BULK_SHEET_PATH':
                item = path[-1]
                if path[2] in ('START_K','END_K','BRANCH_POINTS','MINIMUM_BRANCH_POINT_DISTANCE','GEOMETRIC_RESOLUTION'):
                    return dimensions.measure(self.k)
                if item == 'MAXIMUM_RADICAL_RESIDUAL':
                    return tuple(2*v for v in dimensions.measure(self.q))
                if item in ('SEED_Q','END_Q','REFINEMENT_DIFFERENCE','SHEET_DIFFERENCE','OPPOSITE_SHEET_DIFFERENCE'):
                    return dimensions.measure(self.q)
                return self.unitless
            if key == 'HELMHOLTZ_CHART_OPERAND': return tuple(2*v for v in dimensions.measure(self.k))
            if key == 'RADICAL_RESIDUAL': return tuple(2*v for v in dimensions.measure(self.q))
            if key == 'RECTANGULAR_ROOT_COEFFICIENT_COMMUTATORS':
                return tuple(2*v for v in dimensions.measure(self.k))
            if key == 'RECTANGULAR_EQUATION_COEFFICIENT_RESIDUALS':
                n = int(named(outputs[path[0]],'NULLITY'))
                i = path[-1]//n
                return row_units[i] if path[2]=='RIGHT' else tuple(-v for v in column_units[i])
            if key == 'CONSERVATIVE_SLAB_CURRENT_PAIRING':
                return self.closed_operands[-1]
            if key == 'CLOSED_PHYSICAL_FIELD_RESIDUAL':
                n = int(named(outputs[path[0]], 'NULLITY'))
                i = path[2]//n
                # Strong-equation row units follow the actual field lift.
                strong, lift, _, _, field_units, _ = self.closed_operands
                source = next(j for j in range(5) if strong[i,j] != 0)
                return tuple(a+b for a,b in zip(dimensions.measure(strong[i,source]),field_units[source]))
            if key in ('OMEGA_DERIVATIVE_PAIRING', 'OMEGA_PAIRING_SINGULAR_VALUES'):
                return tuple(-v for v in dimensions.measure(self.r.omega))
            if key == 'OMEGA_PAIRING_DIFFERENCE_QUOTIENTS':
                return tuple((1 if path[3] == 0 else -1)*v
                             for v in dimensions.measure(self.r.omega))
            if key == 'OMEGA_NORMALIZED_LEFT_RESIDUAL':
                n = int(named(outputs[path[0]], 'NULLITY'))
                i = path[2]//n
                return tuple(w-v for w, v in zip(dimensions.measure(self.r.omega), column_units[i]))
            if key in ('RIGHT','LEFT_ADJOINT','RIGHT_RESIDUAL','LEFT_RESIDUAL'):
                n = int(named(outputs[path[0]],'NULLITY'))
                i = path[2]//n
                return {'RIGHT':column_units[i], 'LEFT_ADJOINT':tuple(-v for v in row_units[i]),
                        'RIGHT_RESIDUAL':row_units[i], 'LEFT_RESIDUAL':tuple(-v for v in column_units[i])}[key]
            if key in ('PROJECTOR','PROJECTOR_RESIDUAL'):
                i,j=divmod(path[2],5)
                return tuple(a-b for a,b in zip(column_units[i],column_units[j]))
            # Remaining values are coefficients/diagnostics in the explicitly
            # emitted numerical unit frame, not additional physical fields.
            return self.unitless
        emit('METADATA_FULL_PENCIL_MODE_'+record_kind+'_'+name, self.numeric_metadata(outputs,record_unit))
        return outputs


class EndSpectrumCoverage:
    """Finite algebraic end spectrum in the original five physical fields.

    All carriers are bound before exact rational elimination. Rational Taylor
    bounds isolate the roots of every square-free factor. Their disjoint disks
    and the factor degrees provide finite-root coverage on the algebraic
    radical curve, separately from membership in the Fourier sheet chart.
    This is retained-operator data, not a global profile dispersion relation.
    """

    def __init__(self, modes, strong, full, lift, strong_units):
        self.modes, self.r = modes, modes.r
        self.k, self.q = modes.k, modes.q
        self.strong, self.full, self.lift = strong, full, lift
        self.strong_units = strong_units
        dimensions = PHYSICAL_METADATA.dimensions
        fields = tuple(sp.Function('s11cdReducedField'+name)
                       for name in ('u1','u2','u3','theta','eW'))
        self.field_units = [dimensions.known[f] for f in fields]
        self.row_units = [tuple(a+b for a,b in zip(strong_units[(5*i,)],self.field_units[0]))
                          for i in range(5)]
        self.dimensions = dimensions

    @staticmethod
    def absolute_bounds(value):
        a,b = sp.expand(value).as_real_imag()
        return max(abs(a),abs(b)),abs(a)+abs(b)

    @classmethod
    def isolate(cls, polynomial, digits=50):
        """Exact rational sufficient Rouche inequalities on computed disks."""
        factors = sp.sqf_list(polynomial)[1]
        records, values = [], []
        for factor_index,(factor,multiplicity) in enumerate(factors):
            roots = sp.nroots(factor,n=digits,maxsteps=500)
            refined = sp.nroots(factor,n=digits+30,maxsteps=500)
            for index,root in enumerate(roots):
                real,imag = root.as_real_imag()
                center = sp.Rational(str(real))+sp.I*sp.Rational(str(imag))
                other = min(refined,key=lambda z:abs(complex(z-root)))
                radius = sp.Rational(1,10**(digits//2))
                coefficients = []
                derivative = factor
                for order in range(factor.degree()+1):
                    coefficients.append(sp.expand(derivative.eval(center)/sp.factorial(order)))
                    derivative = derivative.diff()
                lower = cls.absolute_bounds(coefficients[1])[0]*radius
                upper = cls.absolute_bounds(coefficients[0])[1]+sum(
                    cls.absolute_bounds(c)[1]*radius**order for order,c in enumerate(coefficients[2:],2))
                records.append({'FACTOR':factor_index,'FACTOR_DEGREE':factor.degree(),
                    'ROOT_INDEX':index,'MULTIPLICITY':multiplicity,'CENTER':center,'RADIUS':radius,
                    'LINEAR_TERM_LOWER_BOUND':sp.N(lower,20),'REMAINDER_UPPER_BOUND':sp.N(upper,20),
                    'EXACT_BOUND_DIFFERENCE_SIGN':sp.sign(lower-upper),
                    'ONE_ROOT_DISK':bool(lower>upper),
                    'EXACT_TAYLOR_BOUND_SHA256':hashlib.sha256(sp.srepr(cas((coefficients,lower,upper))).encode()).hexdigest(),
                    'PRECISION_REFINEMENT_DIFFERENCE':sp.N(other-root,25),
                    'FACTOR_RESIDUAL':sp.N(factor.eval(other),25)})
                values.append(other)
        separation = []
        for i,a in enumerate(records):
            for j,b in enumerate(records[:i]):
                square = sp.expand((a['CENTER']-b['CENTER'])*sp.conjugate(a['CENTER']-b['CENTER']))
                gap = square-(a['RADIUS']+b['RADIUS'])**2
                separation.append((i,j,sp.sign(gap)))
        covered = sum(r['MULTIPLICITY'] for r in records if r['ONE_ROOT_DISK'])
        result = {'DEGREE':polynomial.degree(), 'DISTINCT_ROOT_COUNT':len(records),
                  'COUNT_WITH_MULTIPLICITY':sum(r['MULTIPLICITY'] for r in records),
                  'ISOLATED_COUNT_WITH_MULTIPLICITY':covered,
                  'DEGREE_COUNT_RESIDUAL':polynomial.degree()-covered,
                  'DISK_SEPARATION_SIGNS':separation,
                  'ALL_DISKS_DISJOINT':all(s>0 for _,_,s in separation),
                  'ROOT_DISKS':records,'DECIMAL_PRECISIONS':(digits,digits+30)}
        result['FINITE_POLYNOMIAL_ROOT_COVERAGE'] = bool(covered==polynomial.degree()
            and result['ALL_DISKS_DISJOINT'] and all(r['ONE_ROOT_DISK'] for r in records))
        return result,values

    def emit(self, tag, value, unit=None, *, heavy=False):
        body = cas(value)
        if heavy:
            numeric = not (body.free_symbols-set(PHYSICAL_METADATA.generators))
            payload = self.modes.compact_fingerprint(body) if numeric else carrier_fingerprint(body)
        else:
            payload = body
        emit(tag,payload)
        emit('METADATA_'+tag,self.modes.numeric_metadata(body,unit or (lambda p:self.dimensions.zero)))

    @staticmethod
    def regularity_criteria(summary, records):
        """Per-bound-input criteria; no generic parameter or sheet coverage."""
        return {'FINITE_ALGEBRAIC_ROOT_COVERAGE':summary['ALGEBRAIC_ROOT_COVERAGE'],
            'DENOMINATOR_DOMAIN':summary['REGULAR_DENOMINATOR_DOMAIN'],
            'NORMAL_LIFT_DOMAIN':summary['REGULAR_NORMAL_LIFT_DOMAIN'],
            'FINITE_NONEMPTY_NULLSPACES':bool(records) and all(
                v['FINITE_PENCIL'] and v.get('NULLITY',0)>0 for v in records),
            'ALGEBRAIC_GEOMETRIC_MULTIPLICITIES':bool(records) and all(
                v.get('ALGEBRAIC_GEOMETRIC_MULTIPLICITY_DIFFERENCE')==0 for v in records),
            'FULL_BASIS_RANKS':bool(records) and all(
                v.get('RIGHT_BASIS_RANK')==v.get('LEFT_BASIS_RANK')==v.get('NULLITY',0)
                and v.get('NULLITY',0)>0 for v in records),
            'NORMAL_DERIVATIVE_PAIRING_RANKS':bool(records) and all(
                v.get('NORMAL_DERIVATIVE_PAIRING_RANK')==v.get('NULLITY',0)
                and v.get('NULLITY',0)>0 for v in records)}

    def construct(self, suffix, end_sign, *, channel_input=None, reference=False, sample_index=0):
        m,d = self.modes,self.dimensions
        join_suffix = ('INPUT_' if channel_input else 'PIT_')+suffix+'_'+str(sample_index)
        algebraic,relation,join = m.analytic(self.strong)
        weak,weak_relation,weak_join = m.analytic(self.full)
        relation_residual = sp.expand(relation-weak_relation)
        self.emit('END_SPECTRUM_BRANCH_JOIN_'+join_suffix,join+weak_join)
        self.emit('END_SPECTRUM_RADICAL_JOIN_'+join_suffix,relation_residual,
                  lambda p:tuple(2*v for v in d.measure(self.q)))
        if any(v!=0 for v in join+weak_join) or relation_residual!=0:
            raise NotImplementedError('end spectrum has unjoined branch operands')
        live = (self.k,self.q,m.eta,m.sigma)
        if channel_input is None:
            mapping = m.sample(weak,relation,sample_index)
            origin = {m.eta:sp.S.Zero,m.sigma:sp.S.Zero}
            kind = 'PIT'
        else:
            mapping = channel_input.mapping(algebraic,relation,live)
            # Include any carrier that is present only in the coordinate lift.
            mapping.update({s:channel_input.parameters[s.name] for s in self.lift.free_symbols-set(live)})
            origin = {m.eta:sp.S.Zero,m.sigma:sp.S.Zero} if reference else channel_input.origin
            kind = 'INPUT'
        name = kind+'_'+suffix+'_'+str(sample_index)
        prefix = 'END_SPECTRUM_'+name
        physical = algebraic.xreplace(mapping).subs(origin).applyfunc(sp.cancel)
        weak_bound = weak.xreplace(mapping).subs(origin).applyfunc(sp.cancel)
        bound_relation = relation.xreplace(mapping)
        unbound_lift = self.lift[:,m.indices]
        lift = unbound_lift.xreplace(mapping)
        dual = unbound_lift.xreplace({s:-s for s in (*self.r.tangents,self.k)}).xreplace(mapping).T
        quotient = weak_bound.extract(m.indices,m.indices)
        pullback = (dual*physical*lift).applyfunc(sp.cancel)
        residual = (quotient-pullback).applyfunc(sp.cancel)
        quotient_units = {(5*i+j,):m.units[(6*ri+cj,)] for i,ri in enumerate(m.indices) for j,cj in enumerate(m.indices)}
        self.emit(prefix+'_BOUND_CARRIERS',tuple((str(k),v) for k,v in mapping.items()),
                  lambda p:d.measure(next(k for k in mapping if str(k)==p[0])))
        self.emit(prefix+'_GRADE_ORIGIN',tuple(origin.items()))
        self.emit(prefix+'_UNIT_FRAME',channel_input.frame if channel_input else ('PIT_L','PIT_T','PIT_M'))
        self.emit(prefix+'_EVALUATION_BINDING',{'BACKGROUND_ORIGIN':tuple(origin.items()),
            'RETAINED_OPERATOR_EVALUATION':True,'CONTINUUM_REEXPANSION_PERFORMED':False,
            'PARAMETER_DOMAIN':{'POSITIVE_REAL_FREQUENCY':mapping[self.r.omega]>0,
                               'FINITE_BOUND_CARRIERS':all(v.is_finite is True for v in mapping.values())}})
        self.emit(prefix+'_SOURCE_GRADE_SUPPORT',tuple((path,tuple(sorted(PHYSICAL_METADATA.coefficients(value))))
            for path,value in leaves(algebraic)))
        self.emit(prefix+'_PHYSICAL_PENCIL',physical,lambda p:self.strong_units[p],heavy=True)
        self.emit(prefix+'_QUOTIENT_PENCIL',quotient,lambda p:quotient_units[p],heavy=True)
        self.emit(prefix+'_PULLBACK_PENCIL',pullback,lambda p:quotient_units[p],heavy=True)
        self.emit(prefix+'_PULLBACK_RESIDUAL',residual,lambda p:quotient_units[p])
        if any(v!=0 for v in residual):
            raise NotImplementedError(('physical/sector end-symbol mismatch',name))
        k_square = sp.solve(bound_relation,self.k**2)[0]
        divisor = sp.Poly(self.k**2-k_square,self.k)
        (numerator,denominator),cleared,row_denominators = m.rational_determinant(physical)
        reduced = sp.rem(sp.Poly(numerator,self.k),divisor).as_expr()
        self.emit(prefix+'_ELIMINATION_K_REMAINDER',reduced.has(self.k))
        if reduced.has(self.k):
            raise NotImplementedError(('non-even physical end determinant',name))
        polynomial = sp.Poly(reduced,self.q)
        if polynomial.is_zero:
            self.emit(prefix+'_DETERMINANT_DOMAIN',{'IDENTICALLY_ZERO':True,'POLYNOMIAL_DEGREE':polynomial.degree()})
            return {'DEFINED':False,'STATUS':'IDENTICALLY_SINGULAR_PHYSICAL_PENCIL','RECORDS':[]}
        # Denominator norm on the two k lifts: exact gcd finds any shared root.
        row_denominator = sp.lcm(row_denominators)
        norm = sp.expand(row_denominator*row_denominator.xreplace({self.k:-self.k}))
        norm = sp.rem(sp.Poly(norm,self.k),divisor).as_expr()
        denominator_polynomial = sp.Poly(sp.fraction(sp.cancel(norm))[0],self.q)
        denominator_common = sp.gcd(polynomial,denominator_polynomial)
        normal_threshold = sp.Poly(sp.fraction(sp.cancel(k_square))[0],self.q)
        normal_common = sp.gcd(polynomial,normal_threshold)
        branch_common = sp.gcd(polynomial,sp.Poly(self.q,self.q))
        self.emit(prefix+'_ELIMINATION_OPERANDS',(cleared,row_denominators,numerator,denominator,polynomial.as_expr()),heavy=True)
        # These are coefficient polynomials in the declared numerical unit
        # frame; physical root and pencil units are emitted independently.
        self.emit(prefix+'_COEFFICIENT_COORDINATE_UNITS',{'K':d.measure(self.k),'Q':d.measure(self.q),
            'PHYSICAL_PENCIL':self.strong_units,'RIGHT_FIELDS':self.field_units,'EQUATION_ROWS':self.row_units})
        self.emit(prefix+'_EXCEPTION_POLYNOMIALS',{'DENOMINATOR_NORM':denominator_polynomial.as_expr(),
            'DENOMINATOR_GCD':denominator_common.as_expr(),'NORMAL_THRESHOLD_GCD':normal_common.as_expr(),
            'RADICAL_BRANCH_GCD':branch_common.as_expr()},heavy=True)
        self.emit(prefix+'_EXCEPTION_DEGREES',{'DENOMINATOR_GCD':denominator_common.degree(),
            'NORMAL_THRESHOLD_GCD':normal_common.degree(),'RADICAL_BRANCH_GCD':branch_common.degree()})
        lift_determinants = (sp.factor(lift.det()),sp.factor(dual.det()))
        lift_determinant_unit = d.measure(unbound_lift.det())
        chart_factor = sp.rem(sp.Poly(sp.expand(sp.prod(lift_determinants)),self.k),divisor).as_expr()
        self.emit(prefix+'_COORDINATE_DETERMINANT_OPERANDS',lift_determinants,
                  lambda p:lift_determinant_unit,heavy=True)
        self.emit(prefix+'_COORDINATE_FACTOR_ON_CURVE',chart_factor,heavy=True)
        self.emit(prefix+'_AXIAL_CHART_IDENTICALLY_SINGULAR',chart_factor==0)
        quotient_fraction,_,_ = m.rational_determinant(quotient)
        quotient_eliminated = sp.rem(sp.Poly(quotient_fraction[0],self.k),divisor).as_expr()
        self.emit(prefix+'_QUOTIENT_ELIMINATION_K_REMAINDER',quotient_eliminated.has(self.k))
        if not quotient_eliminated.has(self.k):
            quotient_polynomial = sp.Poly(quotient_eliminated,self.q)
            self.emit(prefix+'_POLYNOMIAL_DEGREES',{'PHYSICAL':polynomial.degree(),
                'QUOTIENT':quotient_polynomial.degree(),'COORDINATE_FACTOR':sp.Poly(chart_factor,self.q).degree()})
            self.emit(prefix+'_QUOTIENT_ELIMINATION',quotient_polynomial.as_expr(),heavy=True)
        determinant_residual = sp.cancel(quotient_fraction[0]/quotient_fraction[1]
                                        -sp.prod(lift_determinants)*numerator/denominator)
        self.emit(prefix+'_DETERMINANT_PULLBACK_RESIDUAL',determinant_residual)
        certificate,roots = self.isolate(polynomial)
        self.emit(prefix+'_ROOT_COVERAGE',certificate,lambda p:
            d.measure(self.q) if 'ROOT_DISKS' in p and p[-1] in ('CENTER','RADIUS','PRECISION_REFINEMENT_DIFFERENCE') else d.zero)
        # Direct evaluations of the original rational matrix, not the cleared
        # polynomial, supply the physical mode residuals and left/right spaces.
        evaluate = sp.lambdify((self.k,self.q),physical,'numpy',cse=True)
        tangent_q = sp.cancel(-sp.diff(bound_relation,self.k)/sp.diff(bound_relation,self.q))
        normal_derivative = (physical.diff(self.k)+physical.diff(self.q)*tangent_q).applyfunc(sp.cancel)
        evaluate_derivative = sp.lambdify((self.k,self.q),normal_derivative,'numpy',cse=True)
        native_lift = self.lift.xreplace(mapping)
        field_lifts = [native_lift[:,tuple(i for i in range(6) if i!=axis)] for axis in range(3)]
        sheet = BulkSheetPath(bound_relation,self.k,self.q)
        records = []
        for root_index,(qroot,disk) in enumerate(zip(roots,certificate['ROOT_DISKS'])):
            normal_square = k_square.subs(self.q,qroot)
            for direction in (1,-1):
                kroot = sp.N(direction*sp.sqrt(normal_square),50)
                kvalue,qvalue = complex(kroot),complex(qroot)
                matrix = np.asarray(evaluate(kvalue,qvalue),dtype=complex)
                mode_tag = prefix+'_MODE_'+str(len(records))
                record = {'ROOT_DISK_INDEX':root_index,'NORMAL_LIFT_SIGN':direction,
                    'K':kroot,'Q':qroot,'MULTIPLICITY':disk['MULTIPLICITY'],
                    'FINITE_PENCIL':bool(np.isfinite(matrix).all()),
                    'RADICAL_RESIDUAL':sp.N(bound_relation.subs({self.k:kroot,self.q:qroot}),25),
                    'DETERMINANT_NUMERATOR_RESIDUAL':sp.N(numerator.subs({self.k:kroot,self.q:qroot}),25),
                    'ROW_DENOMINATOR_VALUES':tuple(sp.N(v.subs({self.k:kroot,self.q:qroot}),25) for v in row_denominators),
                    'REAL_NORMAL_MOMENTUM':bool(abs(kvalue.imag)<1e-10),
                    'DECAY_AT_END':bool(end_sign*kvalue.imag>1e-10),
                    'NORMAL_THRESHOLD':bool(abs(kvalue)<1e-10), 'RADICAL_BRANCH_POINT':bool(abs(qvalue)<1e-10)}
                membership,path = sheet.classify(kvalue,qvalue)
                record['FIXED_FREQUENCY_SHEET_MEMBERSHIP'] = membership
                self.emit(mode_tag+'_SHEET_PATH',path,lambda p:
                    d.measure(self.k) if p[0] in ('START_K','END_K','BRANCH_POINTS','MINIMUM_BRANCH_POINT_DISTANCE','GEOMETRIC_RESOLUTION') else
                    tuple(2*v for v in d.measure(self.q)) if p[-1]=='MAXIMUM_RADICAL_RESIDUAL' else
                    d.measure(self.q) if p[-1] in ('SEED_Q','END_Q','REFINEMENT_DIFFERENCE','SHEET_DIFFERENCE','OPPOSITE_SHEET_DIFFERENCE') else d.zero)
                if record['FINITE_PENCIL']:
                    left,singular,right_h = np.linalg.svd(matrix)
                    tolerance = 1e-8*max(1.,singular[0])
                    nullity = int(np.sum(singular<tolerance))
                    record.update({'SINGULAR_VALUES':list(map(m.number,singular)),
                                   'RANK_THRESHOLD':tolerance,'NULLITY':nullity,
                                   'ALGEBRAIC_GEOMETRIC_MULTIPLICITY_DIFFERENCE':disk['MULTIPLICITY']-nullity})
                    if nullity:
                        right,dual_mode = right_h.conj().T[:,-nullity:],left[:,-nullity:]
                        record['RIGHT_BASIS_RANK'] = int(np.linalg.matrix_rank(right,tol=1e-9))
                        record['LEFT_BASIS_RANK'] = int(np.linalg.matrix_rank(dual_mode,tol=1e-9))
                        if not record['RADICAL_BRANCH_POINT']:
                            derivative = np.asarray(evaluate_derivative(kvalue,qvalue),dtype=complex)
                            pairing = dual_mode.conj().T@derivative@right
                            record['FINITE_NORMAL_DERIVATIVE_PAIRING'] = bool(np.isfinite(pairing).all())
                            if record['FINITE_NORMAL_DERIVATIVE_PAIRING']:
                                threshold = 1e-9*max(1.,np.linalg.norm(pairing,2))
                                record['NORMAL_DERIVATIVE_PAIRING_THRESHOLD'] = threshold
                                record['NORMAL_DERIVATIVE_PAIRING_RANK'] = int(np.linalg.matrix_rank(pairing,tol=threshold))
                            self.emit(mode_tag+'_NORMAL_DERIVATIVE_PAIRING',sp.ImmutableMatrix(pairing),
                                      lambda p:tuple(-v for v in d.measure(self.k)),heavy=True)
                        self.emit(mode_tag+'_RIGHT',sp.ImmutableMatrix(right),lambda p:self.field_units[p[0]//nullity],heavy=True)
                        self.emit(mode_tag+'_LEFT',sp.ImmutableMatrix(dual_mode),lambda p:tuple(-v for v in self.row_units[p[0]//nullity]),heavy=True)
                        self.emit(mode_tag+'_RIGHT_RESIDUAL',sp.ImmutableMatrix(matrix@right),lambda p:self.row_units[p[0]//nullity])
                        self.emit(mode_tag+'_LEFT_RESIDUAL',sp.ImmutableMatrix(matrix.conj().T@dual_mode),lambda p:tuple(-v for v in self.field_units[p[0]//nullity]))
                        overlap = dual_mode.conj().T@right
                        overlap_rank = np.linalg.matrix_rank(overlap,tol=1e-9)
                        record['LEFT_RIGHT_OVERLAP_RANK'] = int(overlap_rank)
                        record['CLASSIFIER_DEFINED'] = False
                        if overlap_rank==nullity:
                            projector = right@np.linalg.solve(overlap,dual_mode.conj().T)
                            record['PROJECTOR_RANK'] = int(np.linalg.matrix_rank(projector,tol=1e-9))
                            self.emit(mode_tag+'_PROJECTOR',sp.ImmutableMatrix(projector),lambda p:tuple(
                                a-b for a,b in zip(self.field_units[p[0]//5],self.field_units[p[0]%5])),heavy=True)
                            self.emit(mode_tag+'_PROJECTOR_RESIDUAL',sp.ImmutableMatrix(projector@projector-projector),lambda p:tuple(
                                a-b for a,b in zip(self.field_units[p[0]//5],self.field_units[p[0]%5])))
                            charts = [np.array(v.subs(self.k,kroot),dtype=complex) for v in field_lifts]
                            determinants = [np.linalg.det(v) for v in charts]
                            axis = int(np.argmax(np.abs(determinants)))
                            record['FIELD_LIFT_CHART_DETERMINANTS'] = tuple(map(m.number,determinants))
                            record['FIELD_LIFT_CHART_AXIS'] = axis
                            if np.linalg.matrix_rank(charts[axis],tol=1e-9)==5:
                                sector = charts[axis]@np.diag([1,1,0,0,0])@np.linalg.inv(charts[axis])
                                weight = np.trace(sector@projector)/nullity
                                other = np.trace((np.eye(5)-sector)@projector)/nullity
                                record['CLASSIFIER_WEIGHTS'] = (m.number(weight),m.number(other))
                                mixed = abs(weight)>1e-7 and abs(other)>1e-7
                                record['CLASSIFIER_DEFINED'] = not mixed
                                record['CLASSIFIER_STATUS'] = ('MIXED_SUBSPACE_DOMAIN' if mixed else
                                    'TRANSVERSE_LIKE' if abs(weight)>abs(other) else 'THICKNESS_LIKE')
                            else:
                                record['CLASSIFIER_STATUS'] = 'SINGULAR_HELMHOLTZ_FIELD_CHART'
                        else:
                            record['CLASSIFIER_STATUS'] = 'SINGULAR_LEFT_RIGHT_OVERLAP'
                self.emit(mode_tag+'_RECORD',record,lambda p:
                    d.measure(self.k) if p[0]=='K' else d.measure(self.q) if p[0]=='Q' else
                    tuple(-v for v in d.measure(self.k)) if p[0]=='NORMAL_DERIVATIVE_PAIRING_THRESHOLD' else
                    tuple(2*v for v in d.measure(self.q)) if p[0]=='RADICAL_RESIDUAL' else
                    lift_determinant_unit if p[0]=='FIELD_LIFT_CHART_DETERMINANTS' else d.zero)
                records.append(record)
        summary = {'CANDIDATE_COUNT':len(records),
            'FINITE_PENCIL_COUNT':sum(v['FINITE_PENCIL'] for v in records),
            'NULLITY_COUNTS':dict(Counter(v.get('NULLITY',0) for v in records)),
            'SHEET_COUNTS':dict(Counter(str(v['FIXED_FREQUENCY_SHEET_MEMBERSHIP']) for v in records)),
            'CLASSIFIER_COUNTS':dict(Counter(v.get('CLASSIFIER_STATUS','NO_NULLSPACE') for v in records)),
            'ALGEBRAIC_ROOT_COVERAGE':certificate['FINITE_POLYNOMIAL_ROOT_COVERAGE'],
            'REGULAR_DENOMINATOR_DOMAIN':denominator_common.degree()==0,
            'REGULAR_NORMAL_LIFT_DOMAIN':normal_common.degree()==0 and branch_common.degree()==0}
        criteria = self.regularity_criteria(summary,records)
        self.emit(prefix+'_REGULARITY_CRITERIA',criteria)
        summary['REGULAR_ALGEBRAIC_MODE_COVERAGE'] = all(criteria.values())
        summary['PARAMETER_VARIETY_COVERAGE_COMPUTED'] = False
        summary['GLOBAL_SHEET_COVERAGE_COMPUTED'] = False
        self.emit(prefix+'_SUMMARY',summary)
        return {'CERTIFICATE':certificate,'RECORDS':records,'SUMMARY':summary}


class EndModeFrequencyData:
    """Frequency pairing on complete end subspaces, independent of current closure."""

    def __init__(self, modes, strong, units, inputs):
        self.modes,self.strong,self.units,self.inputs = modes,strong,units,inputs
        self.d = PHYSICAL_METADATA.dimensions
        self.field = [self.d.known[sp.Function('s11cdReducedField'+v)] for v in ('u1','u2','u3','theta','eW')]
        self.row = [self.add(units[(5*i,)],self.field[0]) for i in range(5)]
        self.frequency,self.length = self.d.measure(modes.r.omega),self.d.measure(modes.r.ell)
        self.packets = []

    @staticmethod
    def add(*units):
        return tuple(sum(v) for v in zip(*units))

    @staticmethod
    def negative(unit):
        return tuple(-v for v in unit)

    def put(self, name, value, unit=None, heavy=False):
        if isinstance(value,np.ndarray):
            value = sp.ImmutableMatrix(*value.shape,[self.modes.number(v) for v in value.ravel()])
        body = cas(value)
        unit = unit or (lambda path:self.d.zero)
        dimensions = {path:tuple(unit(path)) for path,entry in leaves(body) if not isinstance(entry,Str)}
        self.packets.append({'NAME':name,'VALUE':body,'UNITS':dimensions,'HEAVY':heavy})

    def construct(self, native, coverage, progress=lambda item:None):
        m = self.modes; w,k,q = m.r.omega,m.k,m.q
        algebraic,relation,joins = m.analytic(self.strong)
        mapping = self.inputs.mapping(algebraic,relation,(w,k,q,m.eta,m.sigma))
        graded = algebraic.xreplace(mapping)
        live = graded.subs(self.inputs.origin).applyfunc(sp.cancel)
        wave = relation.xreplace(mapping)
        frequency = self.inputs.parameters[w.name]
        transport = sp.cancel(-wave.diff(w)/wave.diff(q))
        derivative = (live.diff(w)+live.diff(q)*transport).applyfunc(sp.cancel)
        physical = live.subs(w,frequency).applyfunc(sp.cancel)
        bound_wave = wave.subs(w,frequency)
        fixed = algebraic.xreplace({**mapping,w:frequency}).subs(self.inputs.origin).applyfunc(sp.cancel)
        derivative_unit = lambda p:self.add(self.units[p],self.negative(self.frequency))
        self.put('SOURCE_BRANCH_JOIN_RESIDUAL',joins)
        self.put('MATERIAL_PROFILE_BINDINGS',tuple((str(a),b) for a,b in mapping.items()),
                 lambda p:self.d.measure(next(a for a in mapping if str(a)==p[0])))
        self.put('GRADE_ORIGIN',tuple(self.inputs.origin.items()))
        self.put('UNIT_FRAME',self.inputs.frame)
        self.put('SOURCE_GRADE_SUPPORT',tuple((p,tuple(sorted(PHYSICAL_METADATA.coefficients(v)))) for p,v in leaves(algebraic)))
        self.put('LIVE_FREQUENCY_PENCIL',live,lambda p:self.units[p],'carrier')
        self.put('BOUND_PENCIL',physical,lambda p:self.units[p],'carrier')
        self.put('FREQUENCY_DERIVATIVE',derivative,derivative_unit,'carrier')
        self.put('WAVE',wave,lambda p:self.add(self.frequency,self.frequency),'carrier')
        self.put('RADICAL_FREQUENCY_TRANSPORT',transport,None,'carrier')
        self.put('FREQUENCY_WAVE_TANGENCY_RESIDUAL',sp.cancel(wave.diff(w)+wave.diff(q)*transport),lambda p:self.frequency)
        self.put('FREQUENCY_BINDING_RESIDUAL',(physical-fixed).applyfunc(sp.cancel),lambda p:self.units[p])
        self.put('EVALUATION_DOMAIN',{'POSITIVE_FREQUENCY':frequency>0,'FREQUENCY':frequency,
            'RETAINED_OPERATOR_EVALUATION':True,'CONTINUUM_REEXPANSION_PERFORMED':False,
            'GLOBAL_PARAMETER_COVERAGE_COMPUTED':False,'GLOBAL_SHEET_COVERAGE_COMPUTED':False},
            lambda p:self.frequency if p[-1]=='FREQUENCY' else self.d.zero)
        progress({'stage':'frequency_derivative_constructed'})
        (numerator,denominator),cleared,row_denominators = m.rational_determinant(physical)
        square = sp.solve(bound_wave,k**2)[0]
        remainder = sp.rem(sp.Poly(numerator,k),sp.Poly(k**2-square,k)).as_expr()
        polynomial = sp.Poly(remainder,q)
        elimination = (cleared,row_denominators,numerator,denominator,polynomial.as_expr())
        self.put('ELIMINATION_COORDINATES',elimination,None,'carrier')
        certificate = NormalRealityCoverage.construct(polynomial,bound_wave,denominator,k,q,coverage)
        for name,value in certificate.items():
            self.put('NORMAL_REALITY_'+name,value,lambda p:NormalRealityCoverage.unit(p,self.frequency,self.length),
                     'carrier' if name in ('OPERANDS','AXES') else False)
        progress({'stage':'reality_certificates_constructed','disks':len(certificate['DISKS'])})
        evaluate = sp.lambdify((w,k,q),live,'numpy',cse=True)
        diff_evaluate = sp.lambdify((w,k,q),derivative,'numpy',cse=True)
        radical_square = sp.lambdify((w,k),sp.solve(wave,q**2)[0],'numpy',cse=True)
        records = []; wf = float(frequency)
        for index,old in enumerate(native):
            n0 = len(self.packets); prefix = 'MODE_'+str(index)+'_'
            kval,qval = complex(old['K']),complex(old['Q'])
            matrix = np.asarray(evaluate(wf,kval,qval),dtype=complex)
            disk = certificate['DISKS'][int(old['ROOT_DISK_INDEX'])]
            info = {'INDEX':index,'ROOT_DISK_INDEX':int(old['ROOT_DISK_INDEX']),
                'NORMAL_LIFT_SIGN':int(old['NORMAL_LIFT_SIGN']),'K':old['K'],'Q':old['Q'],
                'MULTIPLICITY':int(old['MULTIPLICITY']),'SHEET_MEMBERSHIP':old['FIXED_FREQUENCY_SHEET_MEMBERSHIP'],
                'NORMAL_REALITY_CERTIFICATE':disk,'EXACT_REAL_NORMAL':disk['NORMAL_REALITY_STATUS']=='PROVED_REAL',
                'FINITE_PENCIL':bool(np.isfinite(matrix).all()),'FREQUENCY_NORMALIZATION_DEFINED':False}
            def put(name,value,unit=None,heavy=True):self.put(prefix+name,value,unit,heavy)
            if info['FINITE_PENCIL']:
                left,singular,right_h = np.linalg.svd(matrix)
                threshold = 1e-8*max(1.,singular[0]); n = int(np.sum(singular<threshold))
                info.update({'NULLITY':n,'NATIVE_NULLITY_RESIDUAL':n-int(old['NULLITY']),
                    'ALGEBRAIC_GEOMETRIC_MULTIPLICITY_DIFFERENCE':int(old['MULTIPLICITY'])-n,
                    'COEFFICIENT_FRAME_RANK_THRESHOLD':threshold})
                put('COEFFICIENT_FRAME_SINGULAR_VALUES',singular.reshape(-1,1))
                if n:
                    right,dual = right_h.conj().T[:,-n:],left[:,-n:]
                    info.update({'RIGHT_BASIS_RANK':int(np.linalg.matrix_rank(right,tol=1e-9)),
                                 'LEFT_BASIS_RANK':int(np.linalg.matrix_rank(dual,tol=1e-9))})
                    right_units = lambda p:self.field[p[0]//n]
                    left_units = lambda p:self.negative(self.row[p[0]//n])
                    normalized_left_units = lambda p:self.add(self.frequency,left_units(p))
                    field_map_units = lambda p:self.add(self.field[p[0]//5],self.negative(self.field[p[0]%5]))
                    put('RIGHT_BASIS',right,right_units);put('LEFT_BASIS',dual,left_units)
                    put('RIGHT_KERNEL_RESIDUAL',matrix@right,lambda p:self.row[p[0]//n],False)
                    put('LEFT_KERNEL_RESIDUAL',matrix.conj().T@dual,lambda p:self.negative(self.field[p[0]//n]),False)
                    regular = bool(disk['DENOMINATOR_EXCLUDED'] and not old['RADICAL_BRANCH_POINT'])
                    info['REGULAR_FREQUENCY_CHART'] = regular
                    if regular:
                        dw = np.asarray(diff_evaluate(wf,kval,qval),dtype=complex)
                        pairing = dual.conj().T@dw@right
                        info['FINITE_FREQUENCY_PAIRING'] = bool(np.isfinite(pairing).all())
                        put('FREQUENCY_MATRIX',dw,derivative_unit)
                        put('FREQUENCY_PAIRING',pairing,lambda p:self.negative(self.frequency))
                        for fi,relative in enumerate((1e-4,5e-5)):
                            h = relative*max(1.,abs(wf))
                            qsq = radical_square(wf,kval)
                            shifted = [np.asarray(evaluate(wf+s*h,kval,
                                qval*np.sqrt(complex(radical_square(wf+s*h,kval)/qsq))),dtype=complex) for s in (1,-1)]
                            fd = (shifted[0]-shifted[1])/(2*h)
                            put('FREQUENCY_DIFFERENCE_'+str(fi),fd,derivative_unit)
                            put('FREQUENCY_DIFFERENCE_'+str(fi)+'_RESIDUAL',fd-dw,derivative_unit,False)
                            self.put(prefix+'FREQUENCY_DIFFERENCE_'+str(fi)+'_STEP',h,lambda p:self.frequency)
                        if info['FINITE_FREQUENCY_PAIRING']:
                            ntol = 1e-9*max(1.,np.linalg.norm(pairing,2))
                            rank = int(np.linalg.matrix_rank(pairing,tol=ntol))
                            info.update({'FREQUENCY_PAIRING_RANK':rank,'FREQUENCY_PAIRING_THRESHOLD':ntol,
                                         'FREQUENCY_PAIRING_CONDITION':float(np.linalg.cond(pairing))})
                            if rank==n:
                                normalized = dual@np.linalg.inv(pairing).conj().T
                                projector = right@normalized.conj().T@dw
                                info.update({'FREQUENCY_NORMALIZATION_DEFINED':True,
                                    'FREQUENCY_PROJECTOR_RANK':int(np.linalg.matrix_rank(projector,tol=1e-9))})
                                put('FREQUENCY_NORMALIZED_LEFT',normalized,normalized_left_units)
                                put('FREQUENCY_FIELD_PROJECTOR',projector,field_map_units)
                                put('FREQUENCY_NORMALIZATION_RESIDUAL',normalized.conj().T@dw@right-np.eye(n),None,False)
                                put('NORMALIZED_LEFT_KERNEL_RESIDUAL',matrix.conj().T@normalized,
                                    lambda p:self.add(self.frequency,self.negative(self.field[p[0]//n])),False)
                                put('PROJECTOR_IDEMPOTENCY_RESIDUAL',projector@projector-projector,field_map_units,False)
                                put('PROJECTOR_RANGE_RESIDUAL',projector@right-right,right_units,False)
                                put('PROJECTOR_KERNEL_RESIDUAL',matrix@projector,lambda p:self.units[p],False)
                                # Algorithmic basis changes act on every column, including degenerate spaces.
                                c = np.diag(np.arange(2,n+2).astype(complex))+np.triu(np.full((n,n),.2j),1)
                                b = np.diag(np.arange(3,n+3).astype(complex))+np.tril(np.full((n,n),.3),-1)
                                rt,lt = right@c,dual@b; nt = lt.conj().T@dw@rt
                                ln = lt@np.linalg.inv(nt).conj().T
                                pt = rt@ln.conj().T@dw
                                put('RIGHT_BASIS_CHANGE',c);put('LEFT_BASIS_CHANGE',b)
                                put('CHANGED_FREQUENCY_PAIRING',nt,lambda p:self.negative(self.frequency))
                                put('BASIS_PAIRING_COVARIANCE_RESIDUAL',nt-b.conj().T@pairing@c,lambda p:self.negative(self.frequency),False)
                                put('BASIS_NORMALIZATION_COVARIANCE_RESIDUAL',ln-normalized@np.linalg.inv(c).conj().T,normalized_left_units,False)
                                put('BASIS_PROJECTOR_COVARIANCE_RESIDUAL',pt-projector,field_map_units,False)
            info['COEFFICIENT_FRAME_RESIDUAL_NORMS'] = {p['NAME'][len(prefix):]:float(np.linalg.norm(np.array(p['VALUE'],dtype=complex)))
                for p in self.packets[n0:] if p['NAME'].endswith('_RESIDUAL') and isinstance(p['VALUE'],sp.MatrixBase)}
            def info_units(path):
                if 'NORMAL_REALITY_CERTIFICATE' in path:return NormalRealityCoverage.unit(path,self.frequency,self.length)
                if path[-1]=='K':return self.negative(self.length)
                if path[-1]=='Q':return self.frequency
                if path[-1]=='FREQUENCY_PAIRING_THRESHOLD':return self.negative(self.frequency)
                return self.d.zero
            put('RECORD',info,info_units,False)
            records.append(info);progress({'mode':index,'nullity':info.get('NULLITY'),
                'frequencyNormalizationDefined':info['FREQUENCY_NORMALIZATION_DEFINED']})
        self.put('DIMENSION_CONSTRAINTS',tuple(self.d.constraints))
        return {'RECORDS':records,'NORMAL_REALITY_COVERAGE':certificate,'PACKETS':self.packets,
                'PHYSICAL_PENCIL':physical,'ELIMINATION':elimination,'BOUND_MAPPING':{**mapping,w:frequency},
                'NATIVE_RECORDS':native,'NATIVE_COVERAGE':coverage}

    def emit(self, result, provenance, prefix):
        self.put('PROVENANCE',provenance)
        names = [p['NAME'] for p in result['PACKETS']]+['WRITE_KEYS']
        keys = {name:'s11cd'+''.join(v.title() for v in (prefix+'_'+name).split('_')) for name in names}
        self.put('WRITE_KEYS',keys)
        if len(set(keys.values()))!=len(keys) or set(keys.values())&set(IMPORT_KEYS):
            raise ValueError('end-frequency write-key collision')
        for packet in result['PACKETS']:
            tag = prefix+'_'+packet['NAME'];body = packet['VALUE']
            payload = carrier_fingerprint(body) if packet['HEAVY']=='carrier' else self.modes.compact_fingerprint(body) if packet['HEAVY'] else body
            emit(tag,payload)
            emit('METADATA_'+tag,self.modes.numeric_metadata(body,lambda path:packet['UNITS'][path]))


class EndExceptionalSlice:
    """Operator-derived exceptional conditions on a declared frequency slice.

All other input carriers are bound. Polynomial arithmetic uses numerical
coordinates in the declared unit frame; physical matrix and spectral units
are retained separately. Generic determinant multiplicities and exact minor
identities are distinct from numerical full-basis checks at targeted points.
    """

    def __init__(self, spectrum):
        self.spectrum = spectrum
        self.modes, self.d = spectrum.modes, spectrum.dimensions
        self.w, self.k, self.q = sp.symbols(
            's11cdExceptionFrequencyCoordinate s11cdExceptionNormalCoordinate s11cdExceptionRadicalCoordinate')
        for symbol in (self.w,self.k,self.q):
            self.d.known[symbol] = self.d.zero

    def emit(self, tag, value, unit=None, *, heavy=False, coefficient=False):
        if coefficient:
            # Polynomial coefficients can be exact Gaussian rationals outside
            # binary64 range even when the evaluated physical matrix is modest.
            # Keep this symbolic-coefficient package on the arbitrary-precision
            # carrier path, including when a particular coefficient is constant.
            body=cas(value);name=self.prefix+'_'+tag
            emit(name,carrier_fingerprint(body))
            emit('METADATA_'+name,self.modes.numeric_metadata(body,unit or (lambda p:self.d.zero)))
        else:
            self.spectrum.emit(self.prefix+'_'+tag,value,unit,heavy=heavy)

    @staticmethod
    def real_condition(expression, w):
        polynomial = sp.Poly(sp.fraction(sp.cancel(expression))[0],w,domain=sp.QQ_I)
        real = sp.Poly.from_list([sp.re(c) for c in polynomial.all_coeffs()],w)
        imaginary = sp.Poly.from_list([sp.im(c) for c in polynomial.all_coeffs()],w)
        common = sp.gcd(real,imaginary)
        return real,imaginary,common.monic() if not common.is_zero else common

    @staticmethod
    @lru_cache(maxsize=12)
    def analyze(physical, relation, w, k, q):
        (numerator,denominator),cleared,rows = FullPencilModes.rational_determinant(physical)
        k_square = sp.solve(relation,k**2)[0]
        divisor = sp.Poly(k**2-k_square,k)
        eliminated = sp.cancel(sp.rem(sp.Poly(numerator,k),divisor).as_expr())
        if eliminated.has(k):
            raise NotImplementedError('exceptional-slice normal elimination has an odd remainder')
        polynomial = sp.Poly(sp.fraction(eliminated)[0],q,domain=sp.QQ_I.poly_ring(w))
        if polynomial.is_zero:
            return {'DEFINED':False,'POLYNOMIAL':polynomial,'STATUS':'IDENTICALLY_SINGULAR_SLICE'}
        content,factors = polynomial.sqf_list()
        row_denominator = sp.lcm(rows)
        norm = sp.rem(sp.Poly(sp.expand(row_denominator*row_denominator.xreplace({k:-k})),k),divisor).as_expr()
        tests = {'DENOMINATOR':sp.fraction(sp.cancel(norm))[0],
                 'NORMAL_THRESHOLD':sp.fraction(sp.cancel(k_square))[0],
                 'RADICAL_BRANCH':q}
        loci = {'DETERMINANT_CONTENT':content,'RADICAL_NORMAL_LEADING':sp.Poly(relation,k).LC()}
        rank_records = []
        for i,(factor,multiplicity) in enumerate(factors):
            f = factor.as_expr()
            loci['FACTOR_'+str(i)+'_LEADING'] = factor.LC()
            loci['FACTOR_'+str(i)+'_DISCRIMINANT'] = sp.discriminant(f,q)
            for label,test in tests.items():
                loci['FACTOR_'+str(i)+'_'+label] = sp.resultant(f,test,q)
            for j,(other,_) in enumerate(factors[:i]):
                loci['FACTOR_'+str(j)+'_'+str(i)+'_INTERSECTION'] = sp.resultant(f,other.as_expr(),q)
            # An analytic determinant with local order m has nullity <= m.
            # Vanishing (n-m+1)-minors give the opposite bound on the stated
            # radical/factor domain. No representative point supplies it.
            size = physical.rows-int(multiplicity)+1
            identity_domain = 1<=size<=physical.rows
            basis = sp.groebner((relation,f),k,q,domain=sp.QQ_I.frac_field(w))
            identities = []
            denominators = [sp.denom(sp.cancel(c)) for b in basis.polys for c in sp.Poly(b.as_expr(),k,q).coeffs()]
            for rr in combinations(range(physical.rows),size) if identity_domain else ():
                for cc in combinations(range(physical.cols),size):
                    minor = cleared.extract(rr,cc).det(method='domain-ge')
                    quotients,remainder = basis.reduce(minor)
                    reconstruction = sp.cancel(minor-sum(a*b.as_expr() for a,b in zip(quotients,basis.polys))-remainder)
                    identities.append((rr,cc,minor,tuple(quotients),remainder,reconstruction))
                    denominators.extend(sp.denom(sp.cancel(c)) for a in (*quotients,remainder)
                        for c in sp.Poly(a,k,q).coeffs())
            coefficient_denominator = sp.lcm(denominators)
            loci['FACTOR_'+str(i)+'_RANK_IDENTITY_DENOMINATOR'] = coefficient_denominator
            rank_records.append({'FACTOR_INDEX':i,'GENERIC_DETERMINANT_MULTIPLICITY':multiplicity,
                'RANK_IDENTITY_CONSTRUCTION_DEFINED':identity_domain,
                'TESTED_MINOR_SIZE':size,'MINOR_COUNT':len(identities),
                'GROEBNER_BASIS':tuple(p.as_expr() for p in basis.polys),'MINOR_IDENTITIES':identities,
                'NONZERO_REMAINDER_COUNT':sum(v[4]!=0 for v in identities),
                'IDENTITY_RECONSTRUCTION_RESIDUALS':tuple(v[5] for v in identities),
                'COEFFICIENT_DENOMINATOR':coefficient_denominator})
        real_conditions = {}
        for label,expression in loci.items():
            real,imaginary,common = EndExceptionalSlice.real_condition(expression,w)
            intervals = None if common.is_zero else sp.polys.polytools.intervals(common,eps=sp.Rational(1,10**30))
            real_conditions[label] = {'REAL':real,'IMAGINARY':imaginary,'GCD':common,'INTERVALS':intervals}
        return {'DEFINED':True,'PHYSICAL':physical,'RELATION':relation,'NUMERATOR':numerator,
            'DENOMINATOR':denominator,'CLEARED':cleared,'ROW_DENOMINATORS':rows,
            'POLYNOMIAL':polynomial,'CONTENT':content,'FACTORS':factors,'K_SQUARE':k_square,
            'TESTS':tests,'LOCI':loci,'REAL_CONDITIONS':real_conditions,'RANK_RECORDS':rank_records}

    def construct(self, suffix, channel_input, *, reference=False):
        m,d = self.modes,self.d
        self.prefix = 'END_EXCEPTIONAL_SLICE_INPUT_'+suffix
        algebraic,relation,joins = m.analytic(self.spectrum.strong)
        mapping = channel_input.mapping(algebraic,relation,(m.r.omega,m.k,m.q,m.eta,m.sigma))
        origin = {m.eta:sp.S.Zero,m.sigma:sp.S.Zero} if reference else channel_input.origin
        coordinates = {m.r.omega:self.w,m.k:self.k,m.q:self.q}
        physical = algebraic.xreplace(mapping).subs(origin).xreplace(coordinates).applyfunc(sp.cancel)
        relation = relation.xreplace(mapping).xreplace(coordinates)
        self.emit('UNIT_FRAME',channel_input.frame)
        self.emit('BOUND_CARRIERS',tuple((str(s),v) for s,v in mapping.items()),
                  lambda p:d.measure(next(s for s in mapping if str(s)==p[0])))
        self.emit('GRADE_ORIGIN',tuple((str(s),v) for s,v in origin.items()))
        self.emit('SOURCE_GRADE_SUPPORT',tuple((p,tuple(sorted(PHYSICAL_METADATA.coefficients(v))))
                                              for p,v in leaves(algebraic)))
        self.emit('SPECTRAL_COORDINATE_UNITS',{'FREQUENCY':d.measure(m.r.omega),'NORMAL':d.measure(m.k),
            'RADICAL':d.measure(m.q),'COEFFICIENT_COORDINATES':(self.w,self.k,self.q)})
        self.emit('BRANCH_JOIN_RESIDUALS',joins)
        self.emit('PHYSICAL_MATRIX_COEFFICIENTS',physical,lambda p:self.spectrum.strong_units[p],heavy=True)
        self.emit('RADICAL_RELATION_COEFFICIENTS',relation,heavy=True)
        data = self.analyze(physical,relation,self.w,self.k,self.q)
        self.emit('DETERMINANT_DOMAIN',{'DEFINED':data['DEFINED'],'DEGREE':data['POLYNOMIAL'].degree()})
        if not data['DEFINED']:
            self.emit('STATUS',data['STATUS'])
            return data
        self.emit('ELIMINATION_OPERANDS',(data['CLEARED'],data['ROW_DENOMINATORS'],data['NUMERATOR'],
                  data['DENOMINATOR'],data['POLYNOMIAL'].as_expr()),heavy=True)
        self.emit('SQUARE_FREE_FACTORS',[(f.as_expr(),a) for f,a in data['FACTORS']],heavy=True)
        self.emit('FACTOR_DEGREES_AND_MULTIPLICITIES',[(f.degree(),a) for f,a in data['FACTORS']])
        self.emit('EXCEPTION_TEST_OPERANDS',data['TESTS'],heavy=True)
        self.emit('GENERIC_RANK_IDENTITIES',data['RANK_RECORDS'],heavy=True)
        self.emit('GENERIC_RANK_SUMMARY',[{k:r[k] for k in ('FACTOR_INDEX','GENERIC_DETERMINANT_MULTIPLICITY',
            'RANK_IDENTITY_CONSTRUCTION_DEFINED','TESTED_MINOR_SIZE','MINOR_COUNT','NONZERO_REMAINDER_COUNT','IDENTITY_RECONSTRUCTION_RESIDUALS')}
            for r in data['RANK_RECORDS']])
        target_conditions = {}
        for label,expression in data['LOCI'].items():
            condition = data['REAL_CONDITIONS'][label]
            common = condition['GCD']
            self.emit('LOCUS_'+label,{'ELIMINATION_POLYNOMIAL':expression,'REAL_COEFFICIENT_POLYNOMIAL':condition['REAL'].as_expr(),
                'IMAGINARY_COEFFICIENT_POLYNOMIAL':condition['IMAGINARY'].as_expr(),'REAL_LOCUS_GCD':common.as_expr()},coefficient=True)
            self.emit('LOCUS_'+label+'_REAL_ISOLATION',{'IDENTICALLY_ZERO':common.is_zero,
                'GCD_DEGREE':common.degree(),'REAL_INTERVALS':condition['INTERVALS'] if condition['INTERVALS'] is not None else 'UNRESOLVED_IDENTICAL_LOCUS'})
            if not common.is_zero:
                for root in common.sqf_part().real_roots():
                    if root>=0:
                        target_conditions.setdefault(root,[]).append(label)
        for index,(root,labels) in enumerate(sorted(target_conditions.items(),key=lambda v:float(v[0]))):
            self.target(data,index,root,labels)
        self.emit('COVERAGE_BOUNDARIES',{'BOUND_PARAMETER_FREQUENCY_SLICE':True,
            'REAL_LOCUS_ENUMERATION_DEFINED':all(not c['GCD'].is_zero for c in data['REAL_CONDITIONS'].values()),
            'TARGET_FREQUENCY_COUNT':len(target_conditions),'PARAMETER_VARIETY_ATLAS_COMPUTED':False,
            'COMPLEX_EXCEPTION_ROOT_ISOLATION_COMPUTED':False,'GENERALIZED_DEFECTIVE_MODES_COMPUTED':False,
            'TARGET_PHYSICAL_SHEET_MEMBERSHIP_COMPUTED':False,
            'PROFILE_FREQUENCY_BOUND_POLES_COMPUTED':False})
        return data

    @staticmethod
    @lru_cache(maxsize=24)
    def threshold_data(physical, relation, w, k, q, frequency, normal, radical):
        point = {w:frequency,k:normal,q:radical}
        matrix = physical.subs(point).applyfunc(sp.simplify)
        finite = not any(matrix.has(v) for v in (sp.zoo,sp.nan,sp.oo,-sp.oo))
        result = {'MATRIX':matrix,'FINITE':finite}
        if not finite:
            return result
        right_columns,left_columns = matrix.nullspace(),matrix.conjugate().T.nullspace()
        right = sp.ImmutableMatrix.hstack(*right_columns) if right_columns else sp.ImmutableMatrix.zeros(matrix.cols,0)
        left = sp.ImmutableMatrix.hstack(*left_columns) if left_columns else sp.ImmutableMatrix.zeros(matrix.rows,0)
        q_derivative = sp.cancel(-sp.diff(relation,k)/sp.diff(relation,q))
        derivative = (physical.diff(k)+physical.diff(q)*q_derivative).subs(point).applyfunc(sp.simplify)
        pairing = (left.conjugate().T*derivative*right).applyfunc(sp.simplify)
        result.update({'RIGHT':right,'LEFT':left,'PAIRING':pairing,
            'RIGHT_RESIDUAL':(matrix*right).applyfunc(sp.simplify),
            'LEFT_RESIDUAL':(matrix.conjugate().T*left).applyfunc(sp.simplify),
            'RANKS':{'MATRIX_RANK':matrix.rank(),'RIGHT_NULLITY':right.cols,'LEFT_NULLITY':left.cols,
                'RIGHT_BASIS_RANK':right.rank(),'LEFT_BASIS_RANK':left.rank(),
                'NORMAL_DERIVATIVE_PAIRING_RANK':pairing.rank()}})
        return result

    def target(self, data, index, frequency, labels):
        tag = 'TARGET_'+str(index)
        m,d = self.modes,self.d
        self.emit(tag+'_FREQUENCY',frequency,lambda p:d.measure(m.r.omega))
        self.emit(tag+'_INCIDENT_LOCUS_LABELS',labels)
        self.emit(tag+'_LOCUS_SUBSTITUTION_RESIDUALS',[(label,sp.simplify(data['LOCI'][label].subs(self.w,frequency))) for label in labels])
        specialized = sp.Poly(data['POLYNOMIAL'].as_expr().subs(self.w,frequency),self.q,extension=True)
        self.emit(tag+'_SPECIALIZED_POLYNOMIAL',specialized.as_expr(),heavy=True)
        if frequency==0 or specialized.is_zero:
            self.emit(tag+'_STATUS','UNRESOLVED_ZERO_FREQUENCY_INTERSECTION' if frequency==0 else 'UNRESOLVED_IDENTICAL_DETERMINANT')
            return
        threshold = sp.Poly(data['TESTS']['NORMAL_THRESHOLD'].subs(self.w,frequency),self.q,extension=True)
        common = sp.gcd(specialized,threshold)
        self.emit(tag+'_NORMAL_THRESHOLD_GCD',common.as_expr(),heavy=True)
        self.emit(tag+'_NORMAL_THRESHOLD_GCD_DEGREE',common.degree())
        if common.degree()==0:
            self.emit(tag+'_STATUS','UNRESOLVED_NONTHRESHOLD_EXCEPTION')
            return
        roots = sp.solve(common.as_expr(),self.q)
        self.emit(tag+'_RADICAL_ROOT_COUNT',len(roots))
        self.emit(tag+'_RADICAL_ROOT_DEGREE_RESIDUAL',common.sqf_part().degree()-len(roots))
        for j,qroot in enumerate(roots):
            prefix = tag+'_NORMAL_THRESHOLD_'+str(j)
            normal_roots = sp.solve(data['RELATION'].subs({self.w:frequency,self.q:qroot}),self.k)
            self.emit(prefix+'_NORMAL_LIFT_ROOT_COUNT',len(normal_roots))
            if len(normal_roots)!=1:
                self.emit(prefix+'_STATUS','UNRESOLVED_THRESHOLD_NORMAL_LIFT')
                continue
            point = {self.w:frequency,self.q:qroot,self.k:normal_roots[0]}
            denominator_values = tuple(sp.simplify(v.subs(point)) for v in data['ROW_DENOMINATORS'])
            self.emit(prefix+'_POINT',{'OMEGA':frequency,'K':point[self.k],'Q':qroot},lambda p:
                      d.measure(m.k) if p[0]=='K' else d.measure(m.q) if p[0]=='Q' else d.measure(m.r.omega))
            self.emit(prefix+'_RADICAL_RESIDUAL',sp.simplify(data['RELATION'].subs(point)),
                      lambda p:tuple(2*v for v in d.measure(m.q)))
            self.emit(prefix+'_DENOMINATOR_COEFFICIENT_VALUES',denominator_values)
            if qroot==0:
                self.emit(prefix+'_STATUS','UNRESOLVED_RADICAL_BRANCH_THRESHOLD_INTERSECTION')
                continue
            computed = self.threshold_data(data['PHYSICAL'],data['RELATION'],self.w,self.k,self.q,
                                           frequency,point[self.k],qroot)
            matrix = computed['MATRIX']
            finite = computed['FINITE'] and all(v!=0 for v in denominator_values)
            self.emit(prefix+'_FINITE_PENCIL_DOMAIN',finite)
            if not finite:
                self.emit(prefix+'_STATUS','UNRESOLVED_DENOMINATOR_INTERSECTION')
                continue
            self.emit(prefix+'_PHYSICAL_MATRIX',matrix,lambda p:self.spectrum.strong_units[p],heavy=True)
            # Exact bases use the original rational physical matrix after the
            # computed exceptional substitution, not the cleared determinant.
            right,left = computed['RIGHT'],computed['LEFT']
            nullity = right.cols
            self.emit(prefix+'_RIGHT_BASIS',right,lambda p:self.spectrum.field_units[p[0]//nullity],heavy=True)
            self.emit(prefix+'_LEFT_BASIS',left,lambda p:tuple(-v for v in self.spectrum.row_units[p[0]//left.cols]),heavy=True)
            self.emit(prefix+'_RIGHT_RESIDUAL',computed['RIGHT_RESIDUAL'],lambda p:self.spectrum.row_units[p[0]//nullity])
            self.emit(prefix+'_LEFT_RESIDUAL',computed['LEFT_RESIDUAL'],
                      lambda p:tuple(-v for v in self.spectrum.field_units[p[0]//left.cols]))
            self.emit(prefix+'_NORMAL_DERIVATIVE_PAIRING',computed['PAIRING'],lambda p:tuple(-v for v in d.measure(m.k)),heavy=True)
            self.emit(prefix+'_RANK_RECORD',{**computed['RANKS'],
                'NORMAL_LIFT_COUNT':len(set(normal_roots)),
                'REGULAR_NORMAL_LIFT_DOMAIN':bool(sp.simplify(data['K_SQUARE'].subs(point))!=0 and qroot!=0)})
            self.emit(prefix+'_STATUS','THRESHOLD_SUBSPACE_COMPUTED_GENERALIZED_NORMAL_MODES_UNRESOLVED')


class NormalTaylorChains:
    """Full root-chain spaces of a bound analytic matrix germ.

    Block equations come from substituting a vector Taylor ansatz into the
    computed matrix Taylor series. Exact arithmetic stays in one algebraic
    coefficient field; general-purpose simplify can be costly on these roots.
    A finite cap with no stabilized kernel count is explicitly unresolved.
    """

    @staticmethod
    def coefficient_frame(domain):
        # Convert the expression tree arithmetically. Asking from_sympy to
        # find a fresh minimal polynomial for each whole rational expression
        # needlessly repeats number-field isomorphism and factorization work.
        @lru_cache(maxsize=None)
        def convert(value):
            if value.is_Rational:
                return domain.convert(value)
            if value.is_Add:
                return sum((convert(v) for v in value.args),domain.zero)
            if value.is_Mul:
                result = domain.one
                for factor in value.args:
                    result *= convert(factor)
                return result
            if value.is_Pow and value.exp.is_Integer:
                return convert(value.base)**int(value.exp)
            return domain.from_sympy(value)
        return convert

    @staticmethod
    def construct(jets, domain):
        from sympy.polys.matrices import DomainMatrix
        n = jets[0].rows
        convert = NormalTaylorChains.coefficient_frame(domain)
        matrices = []
        for matrix in jets:
            entries = {}
            for i in range(n):
                row = {j:convert(matrix[i,j]) for j in range(n)}
                row = {j:v for j,v in row.items() if v!=domain.zero}
                if row: entries[i] = row
            matrices.append(DomainMatrix(entries,(n,n),domain))
        blocks, kernels, counts, increments = [], [], [], []
        previous = 0
        for depth in range(1,len(jets)+1):
            rows = [DomainMatrix.hstack(*[matrices[i-j] if i>=j else
                    DomainMatrix.zeros((n,n),domain) for j in range(depth)]) for i in range(depth)]
            block = DomainMatrix.vstack(*rows)
            kernel = block.nullspace(divide_last=True).transpose()
            blocks.append(block); kernels.append(kernel)
            counts.append(kernel.shape[1]); increments.append(counts[-1]-previous)
            previous = counts[-1]
            if increments[-1]==0:
                break
        closed = increments[-1]==0
        chains = []
        chosen = DomainMatrix.zeros((n,0),domain)
        if closed:
            # The first coefficient of ker(T_d) spans chains of length >= d.
            # Choose its complement to already-selected longer root chains.
            for depth in range(len(kernels)-1,0,-1):
                kernel = kernels[depth-1]
                leading = kernel.extract(list(range(n)),list(range(kernel.shape[1])))
                old_count = chosen.shape[1]
                combined = chosen.hstack(leading)
                _,pivots = combined.rref()
                for column in (p-old_count for p in pivots if p>=old_count):
                    vector = kernel.extract(list(range(n*depth)),[column])
                    coefficients = tuple(vector.extract(list(range(n*j,n*(j+1))),[0]) for j in range(depth))
                    residuals = tuple(sum((matrices[a].matmul(coefficients[b-a]) for a in range(b+1)),
                                         DomainMatrix.zeros((n,1),domain)) for b in range(depth))
                    chains.append({'LENGTH':depth,'COEFFICIENTS':tuple(v.to_Matrix() for v in coefficients),
                                   'EQUATION_RESIDUALS':tuple(v.to_Matrix() for v in residuals)})
                    chosen = chosen.hstack(coefficients[0])
        result = {'KERNEL_COUNTS':tuple(counts),'KERNEL_INCREMENTS':tuple(increments),
            'STABILIZED':closed,'BLOCKS':tuple(b.to_Matrix() for b in blocks),
            'KERNELS':tuple(b.to_Matrix() for b in kernels),
            'BLOCK_RANKS':tuple(b.rank() for b in blocks),
            'KERNEL_BASIS_RANKS':tuple(v.rank() for v in kernels),
            'RANK_NULLITY_RESIDUALS':tuple(b.shape[1]-b.rank()-v.rank() for b,v in zip(blocks,kernels)),
            'BLOCK_KERNEL_RESIDUALS':tuple(b.matmul(v).to_Matrix() for b,v in zip(blocks,kernels)),
            'ROOT_SPACE_DIMENSION':counts[0]}
        if closed:
            result.update({'CHAINS':chains,'LEADING_BASIS_RANK':chosen.rank(),
                'CHAIN_LENGTH_SUM':sum(v['LENGTH'] for v in chains),
                'CHAIN_COUNT_RESIDUAL':counts[0]-len(chains),
                'TOTAL_MULTIPLICITY_RESIDUAL':counts[-1]-sum(v['LENGTH'] for v in chains)})
        else:
            result['STATUS'] = 'UNRESOLVED_CHAIN_LENGTH_BEYOND_TAYLOR_CAP'
        return result

    @staticmethod
    @lru_cache(maxsize=24)
    def at_point(physical, relation, w, k, q, frequency, normal, radical, maximum_order=4):
        domain = sp.QQ.algebraic_field(sp.I,frequency,normal,radical)
        convert = NormalTaylorChains.coefficient_frame(domain)
        normal_form = lambda value: domain.to_sympy(convert(value))
        point = {k:normal,q:radical}
        germ = physical.subs(w,frequency)
        curve = relation.subs(w,frequency)
        slope = sp.cancel(-curve.diff(k)/curve.diff(q))
        matrix_derivative, radical_derivative = germ, q
        jets, radical_jets = [], []
        for order in range(maximum_order+1):
            jets.append((matrix_derivative.subs(point)/sp.factorial(order)).applyfunc(normal_form))
            radical_jets.append(normal_form(radical_derivative.subs(point)/sp.factorial(order)))
            right = NormalTaylorChains.construct(jets,domain)
            left = NormalTaylorChains.construct([j.conjugate().T for j in jets],domain)
            if right['STABILIZED'] and left['STABILIZED']:
                break
            matrix_derivative = (matrix_derivative.diff(k)+matrix_derivative.diff(q)*slope).applyfunc(sp.cancel)
            radical_derivative = sp.cancel(radical_derivative.diff(k)+radical_derivative.diff(q)*slope)
        delta = sp.Symbol('s11cdThresholdNormalIncrement')
        radical_ansatz = sum(v*delta**j for j,v in enumerate(radical_jets))
        curve_residual = sp.Poly(sp.expand(curve.subs({k:normal+delta,q:radical_ansatz})),delta)
        return {'JETS':tuple(jets),'RADICAL_JETS':tuple(radical_jets),'RIGHT':right,'LEFT':left,
            'RADICAL_TAYLOR_RESIDUALS':tuple(normal_form(curve_residual.nth(j)) for j in range(len(jets))),
            'COMPUTED_TAYLOR_ORDER':len(jets)-1,'MAXIMUM_TAYLOR_ORDER':maximum_order,
            'COEFFICIENT_FIELD_DEGREE':domain.ext.minpoly.degree()}


class ThresholdModeAudit(EndExceptionalSlice):
    """Generalized normal modes at enumerated finite threshold points."""

    def __init__(self,spectrum,bindings):
        super().__init__(spectrum)
        self.bindings = bindings

    def construct(self,suffix,channel_input,*,reference=False,end_data):
        self.prefix = 'THRESHOLD_MODE_INPUT_'+suffix
        m,d = self.modes,self.d
        self.emit('UNIT_FRAME',channel_input.frame)
        self.emit('SOURCE_SCOPE',{'BOUND_PARAMETER_FREQUENCY_SLICE':True,
            'PROFILE_FREQUENCY_BOUND_POLES_COMPUTED':False,'FLUX_NORMALIZATION_COMPUTED':False})
        data = end_data
        self.emit('DETERMINANT_DOMAIN',data['DEFINED'])
        if not data['DEFINED']:
            return {'DEFINED':False}
        algebraic,relation,joins = m.analytic(self.spectrum.strong)
        mapping = channel_input.mapping(algebraic,relation,(m.r.omega,m.k,m.q,m.eta,m.sigma))
        origin = {m.eta:0,m.sigma:0} if reference else channel_input.origin
        coordinates = {m.r.omega:self.w,m.k:self.k,m.q:self.q}
        self.original = self.spectrum.strong.xreplace(mapping).subs(origin).xreplace(coordinates)
        square = sp.solve(data['RELATION'],self.q**2)[0]
        scale = sp.sqrt(-sp.Poly(square,self.k).nth(2))
        seeds = tuple(scale*rhs.xreplace(dict(zip(lhs.args,(*m.r.tangents,m.k)))).xreplace(mapping).xreplace(coordinates)
                      for lhs,rhs in self.bindings)
        self.seed = seeds[0]
        self.transport = JointBulkSheetPath(data['RELATION'],self.w,self.k,self.q,self.seed)
        self.emit('BOUND_CARRIERS',tuple((str(s),v) for s,v in mapping.items()),
                  lambda p:d.measure(next(s for s in mapping if str(s)==p[0])))
        self.emit('GRADE_ORIGIN',tuple((str(s),v) for s,v in origin.items()))
        self.emit('SOURCE_GRADE_SUPPORT',tuple((p,tuple(sorted(PHYSICAL_METADATA.coefficients(v))))
                                              for p,v in leaves(algebraic)))
        self.emit('SOURCE_BRANCH_SEED',self.seed,lambda p:d.measure(m.q),heavy=True)
        self.emit('SOURCE_BRANCH_JOIN_RESIDUALS',joins)
        self.emit('SOURCE_SEED_JOIN_RESIDUALS',tuple(sp.simplify(v-self.seed) for v in seeds),lambda p:d.measure(m.q))
        targets = {}
        for label,condition in data['REAL_CONDITIONS'].items():
            if not condition['GCD'].is_zero:
                for root in condition['GCD'].sqf_part().real_roots():
                    if root>=0: targets.setdefault(root,[]).append(label)
        summaries = []
        for i,(frequency,labels) in enumerate(sorted(targets.items(),key=lambda item:float(item[0]))):
            tag = 'TARGET_'+str(i)
            self.emit(tag+'_FREQUENCY',frequency,lambda p:d.measure(m.r.omega))
            self.emit(tag+'_LOCUS_LABELS',labels)
            if frequency==0:
                self.emit(tag+'_DOMAIN',{'DEFINED':False,'STATUS':'UNRESOLVED_ZERO_FREQUENCY_INTERSECTION'})
                continue
            specialized = sp.Poly(data['POLYNOMIAL'].as_expr().subs(self.w,frequency),self.q,extension=True)
            threshold = sp.Poly(data['TESTS']['NORMAL_THRESHOLD'].subs(self.w,frequency),self.q,extension=True)
            common = sp.gcd(specialized,threshold)
            roots = sp.solve(common.as_expr(),self.q)
            self.emit(tag+'_ROOT_CENSUS',{'POLYNOMIAL_DEGREE':common.sqf_part().degree(),
                'COMPUTED_ROOT_COUNT':len(roots),'DEGREE_RESIDUAL':common.sqf_part().degree()-len(roots)})
            for j,radical in enumerate(roots):
                name = tag+'_ROOT_'+str(j)
                lifts = sp.solve(data['RELATION'].subs({self.w:frequency,self.q:radical}),self.k)
                self.emit(name+'_NORMAL_LIFTS',lifts,lambda p:d.measure(m.k))
                for h,normal in enumerate(lifts):
                    summaries.append(self.point(data,name+'_LIFT_'+str(h),frequency,normal,radical))
        summary = {'TARGET_FREQUENCY_COUNT':len(targets),'POINT_COUNT':len(summaries),
            'GENERALIZED_POINT_COUNT':sum(v.get('CHAINS_COMPLETE',False) for v in summaries),
            'POINTS':summaries,'GLOBAL_PARAMETER_OR_SHEET_ATLAS_COMPUTED':False,
            'PROFILE_FREQUENCY_BOUND_POLES_COMPUTED':False}
        self.emit('SUMMARY',summary)
        return summary

    def point(self,data,tag,frequency,normal,radical):
        m,d = self.modes,self.d
        point = {self.w:frequency,self.k:normal,self.q:radical}
        self.emit(tag+'_POINT',{'OMEGA':frequency,'K':normal,'Q':radical},lambda p:
                  d.measure(m.r.omega) if p[0]=='OMEGA' else d.measure(m.k) if p[0]=='K' else d.measure(m.q))
        denominators = tuple(sp.cancel(v.subs(point)) for v in data['ROW_DENOMINATORS'])
        self.emit(tag+'_DENOMINATOR_COEFFICIENT_VALUES',denominators)
        finite = all(v!=0 and not v.has(sp.zoo,sp.nan,sp.oo,-sp.oo) for v in denominators)
        domain = {'FINITE_DENOMINATOR_DOMAIN':finite,'ANALYTIC_RADICAL_CHART':radical!=0}
        self.emit(tag+'_DOMAIN',domain)
        if not finite or radical==0:
            return {**domain,'CHAINS_COMPLETE':False}
        computed = NormalTaylorChains.at_point(data['PHYSICAL'],data['RELATION'],self.w,self.k,self.q,
                                               frequency,normal,radical)
        k_unit,q_unit = d.measure(m.k),d.measure(m.q)
        subtract = lambda a,b,j:tuple(x-j*y for x,y in zip(a,b))
        self.emit(tag+'_TAYLOR_DOMAIN',{k:computed[k] for k in
            ('COMPUTED_TAYLOR_ORDER','MAXIMUM_TAYLOR_ORDER','COEFFICIENT_FIELD_DEGREE')})
        for j,matrix in enumerate(computed['JETS']):
            self.emit(tag+'_MATRIX_TAYLOR_'+str(j),matrix,
                      lambda p,j=j:subtract(self.spectrum.strong_units[p],k_unit,j),heavy=True)
            self.emit(tag+'_RADICAL_TAYLOR_'+str(j),computed['RADICAL_JETS'][j],
                      lambda p,j=j:subtract(q_unit,k_unit,j))
            self.emit(tag+'_RADICAL_TAYLOR_RESIDUAL_'+str(j),computed['RADICAL_TAYLOR_RESIDUALS'][j],
                      lambda p,j=j:subtract(tuple(2*x for x in q_unit),k_unit,j))
        for side in ('RIGHT','LEFT'):
            record = computed[side]
            self.emit(tag+'_'+side+'_CENSUS',{k:v for k,v in record.items() if k not in
                ('BLOCKS','KERNELS','BLOCK_KERNEL_RESIDUALS','CHAINS')})
            # Block systems act on Taylor coefficients in the stated unit
            # frame. Their coefficient matrices have no physical dimension;
            # restored coefficient units accompany the physical chains below.
            units = self.spectrum.field_units if side=='RIGHT' else tuple(tuple(-v for v in row) for row in self.spectrum.row_units)
            row_units = self.spectrum.row_units if side=='RIGHT' else tuple(tuple(-v for v in row) for row in self.spectrum.field_units)
            for index,(block,kernel,residual) in enumerate(zip(record['BLOCKS'],record['KERNELS'],record['BLOCK_KERNEL_RESIDUALS'])):
                self.emit(tag+'_'+side+'_UNIT_FRAME_'+str(index)+'_COORDINATE_UNITS',
                    {'ROW_COEFFICIENT_UNITS':tuple(subtract(unit,k_unit,j) for j in range(index+1) for unit in row_units),
                     'COLUMN_COEFFICIENT_UNITS':tuple(subtract(unit,k_unit,j) for j in range(index+1) for unit in units)})
                for label,value in (('BLOCK',block),('KERNEL',kernel),('RESIDUAL',residual)):
                    self.emit(tag+'_'+side+'_UNIT_FRAME_'+str(index)+'_'+label,value,heavy=label!='RESIDUAL')
            jets = computed['JETS'] if side=='RIGHT' else tuple(j.conjugate().T for j in computed['JETS'])
            for index,chain in enumerate(record.get('CHAINS',())):
                name = tag+'_'+side+'_CHAIN_'+str(index)
                self.emit(name+'_LENGTH',chain['LENGTH'])
                for j,(coefficient,residual) in enumerate(zip(chain['COEFFICIENTS'],chain['EQUATION_RESIDUALS'])):
                    self.emit(name+'_COEFFICIENT_'+str(j),coefficient,
                              lambda p,j=j:subtract(units[p[0]],k_unit,j),heavy=True)
                    self.emit(name+'_EQUATION_RESIDUAL_'+str(j),residual,
                              lambda p,j=j:subtract(row_units[p[0]],k_unit,j))
                self.plane_wave(name,chain,jets,normal,units,row_units)
        connection = self.connection(data,tag,frequency,normal,radical,computed)
        summary = {**domain,'CHAINS_COMPLETE':all(computed[s]['STABILIZED'] and
            computed[s]['CHAIN_COUNT_RESIDUAL']==0 and computed[s]['TOTAL_MULTIPLICITY_RESIDUAL']==0
            for s in ('RIGHT','LEFT')),
            'RIGHT_LENGTHS':tuple(c['LENGTH'] for c in computed['RIGHT']['CHAINS']) if computed['RIGHT']['STABILIZED'] else 'UNRESOLVED',
            'LEFT_LENGTHS':tuple(c['LENGTH'] for c in computed['LEFT']['CHAINS']) if computed['LEFT']['STABILIZED'] else 'UNRESOLVED',
            'CONNECTION':connection}
        self.emit(tag+'_SUMMARY',summary)
        return summary

    def plane_wave(self,tag,chain,jets,normal,units,row_units):
        # Polynomial coefficient of the root-function plane-wave ansatz.
        delta,z = sp.symbols('s11cdThresholdAnsatzIncrement s11cdThresholdNormalPosition',real=True)
        length = chain['LENGTH']
        polynomial = sum((c*delta**j for j,c in enumerate(chain['COEFFICIENTS'])),sp.zeros(jets[0].cols,1))
        wave = sp.exp(sp.I*(normal+delta)*z)*polynomial
        mode = (wave.diff(delta,length-1).subs(delta,0)/sp.factorial(length-1))*sp.exp(-sp.I*normal*z)
        mode = mode.applyfunc(sp.expand)
        derivative_symbol = sp.diff(sp.exp(sp.I*delta*z),z)/sp.exp(sp.I*delta*z)/delta
        applied = sum((matrix*mode.diff(z,j)/derivative_symbol**j for j,matrix in enumerate(jets)),
                      sp.zeros(jets[0].rows,1)).applyfunc(sp.expand)
        k_unit = self.d.measure(self.modes.k)
        self.emit(tag+'_FOURIER_DIFFERENTIAL_COEFFICIENT',sp.simplify(derivative_symbol))
        for power in range(length):
            coefficient = mode.applyfunc(lambda v:sp.expand(v).coeff(z,power))
            residual = applied.applyfunc(lambda v:sp.cancel(sp.expand(v).coeff(z,power),extension=True))
            offset = length-1-power
            self.emit(tag+'_PLANE_POLYNOMIAL_'+str(power),coefficient,
                      lambda p:tuple(a-offset*b for a,b in zip(units[p[0]],k_unit)),heavy=True)
            self.emit(tag+'_PLANE_EQUATION_RESIDUAL_'+str(power),residual,
                      lambda p:tuple(a-offset*b for a,b in zip(row_units[p[0]],k_unit)))

    def path_unit(self,path):
        # Reuse the joint-path schema's restored units without substituting
        # the threshold's normal coalescence for a bulk branch point.
        proxy = BulkContinuationAudit.__new__(BulkContinuationAudit)
        proxy.r,proxy.modes,proxy.dimensions = self.modes.r,self.modes,self.d
        return proxy.path_unit(path)

    def connection(self,data,tag,frequency,normal,radical,computed):
        m,d = self.modes,self.d
        point = {self.w:frequency,self.k:normal,self.q:radical}
        source_root = sp.simplify(self.seed.subs(point))
        source = self.original.subs(point).applyfunc(sp.simplify)
        source_residual = (source-computed['JETS'][0]).applyfunc(sp.cancel)
        self.emit(tag+'_SOURCE_ROOT',source_root,lambda p:d.measure(m.q))
        self.emit(tag+'_SOURCE_ROOT_RESIDUAL',sp.simplify(radical-source_root),lambda p:d.measure(m.q))
        self.emit(tag+'_SOURCE_MATRIX',source,lambda p:self.spectrum.strong_units[p],heavy=True)
        self.emit(tag+'_SOURCE_MATRIX_RESIDUAL',source_residual,lambda p:self.spectrum.strong_units[p],heavy=True)
        on_source = bool(sp.simplify(radical-source_root)==0)
        self.emit(tag+'_REAL_AXIS_SHEET_JOIN',{'MATCHES_REDUCED_SEED':on_source,
            'MATRIX_JOIN_ZERO':all(v==0 for v in source_residual),
            'RADICAL_CHART_REGULAR':radical!=0,'FLUX_CHANNEL_LABEL_COMPUTED':False})
        incident = [(i,f,a) for i,(f,a) in enumerate(data['FACTORS'])
                    if sp.simplify(f.as_expr().subs(point))==0]
        charts = []
        valuations = []
        field = sp.QQ.algebraic_field(sp.I,frequency,normal,radical)
        convert = NormalTaylorChains.coefficient_frame(field)
        slope = sp.cancel(-data['RELATION'].diff(self.k)/data['RELATION'].diff(self.q))
        for i,factor,multiplicity in incident:
            derivative = factor.as_expr()
            coefficients = []
            for order in range(9):
                value = field.to_sympy(convert(derivative.subs(point)/sp.factorial(order)))
                coefficients.append(value)
                if value!=0: break
                derivative = sp.cancel(derivative.diff(self.k)+derivative.diff(self.q)*slope)
            resolved = coefficients[-1]!=0
            valuation = len(coefficients)-1 if resolved else None
            self.emit(tag+'_FACTOR_'+str(i)+'_LOCAL_TAYLOR_COEFFICIENTS',coefficients,coefficient=True)
            valuations.append({'FACTOR_INDEX':i,'DETERMINANT_MULTIPLICITY':multiplicity,
                'VALUATION_DEFINED':resolved,'NORMAL_ORDER':valuation if resolved else 'UNRESOLVED_TAYLOR_CAP'})
        self.emit(tag+'_DETERMINANT_LOCAL_VALUATIONS',valuations)
        if all(v['VALUATION_DEFINED'] for v in valuations) and all(computed[s]['STABILIZED'] for s in ('RIGHT','LEFT')):
            total = sum(v['DETERMINANT_MULTIPLICITY']*v['NORMAL_ORDER'] for v in valuations)
            self.emit(tag+'_DETERMINANT_CHAIN_MULTIPLICITY_RESIDUAL',
                {side:total-computed[side]['CHAIN_LENGTH_SUM'] for side in ('RIGHT','LEFT')})
        x = sp.Symbol('s11cdThresholdNormalSquareCoordinate')
        for i,factor,multiplicity in incident:
            polynomial = sp.Poly(sp.resultant(factor.as_expr(),data['RELATION'],self.q),self.k,
                                 domain=sp.QQ_I.poly_ring(self.w)).sqf_part()
            odd = sum(c*self.k**powers[0] for powers,c in polynomial.terms() if powers[0]%2)
            squared = sum(c*x**(powers[0]//2) for powers,c in polynomial.terms() if powers[0]%2==0)
            poly = sp.Poly(squared,x)
            name = tag+'_UNFOLDING_'+str(i)
            self.emit(name+'_FACTOR',factor.as_expr(),coefficient=True)
            self.emit(name+'_NORMAL_PROJECTION',polynomial.as_expr(),coefficient=True)
            self.emit(name+'_ODD_PROJECTION_RESIDUAL',odd,coefficient=True)
            self.emit(name+'_DOMAIN',{'FACTOR_MULTIPLICITY':multiplicity,'NORMAL_SQUARE_DEGREE':poly.degree(),
                'EVEN_PROJECTION':odd==0})
            if odd!=0 or poly.degree()!=1:
                charts.append({'DEFINED':False,'FACTOR_INDEX':i,'STATUS':'UNRESOLVED_NONQUADRATIC_NORMAL_UNFOLDING'})
                continue
            square = sp.cancel(sp.solve(poly.as_expr(),x)[0])
            self.emit(name+'_NORMAL_SQUARE',square,lambda p:tuple(2*v for v in d.measure(m.k)),heavy=True)
            self.emit(name+'_THRESHOLD_RESIDUAL',sp.simplify(square.subs(self.w,frequency)-normal**2),
                      lambda p:tuple(2*v for v in d.measure(m.k)))
            self.emit(name+'_FREQUENCY_SLOPE',sp.simplify(square.diff(self.w).subs(self.w,frequency)),
                      lambda p:tuple(2*a-b for a,b in zip(d.measure(m.k),d.measure(m.r.omega))))
            charts.append(self.local_chart(data,name,frequency,normal,radical,computed,square))
        return {'MATCHES_REDUCED_REAL_AXIS_SEED':on_source,'INCIDENT_FACTOR_COUNT':len(incident),'CHARTS':charts,
                'GLOBAL_SHEET_ATLAS_COMPUTED':False,'FLUX_CHANNEL_LABEL_COMPUTED':False}

    def local_chart(self,data,tag,frequency,normal,radical,computed,square):
        m,d = self.modes,self.d
        disk = self.exception_disk(data,tag,frequency)
        bases = {}
        for side in ('RIGHT','LEFT'):
            chains = computed[side].get('CHAINS')
            if chains is None or not chains or any(c['LENGTH']!=2 for c in chains):
                self.emit(tag+'_LOCAL_MODE_DOMAIN',{'DEFINED':False,'SIDE':side,
                    'STATUS':'UNRESOLVED_NONUNIFORM_LENGTH_TWO_SECANT_CHART'})
                return {'DEFINED':False}
            zeroth = np.asarray(sp.Matrix.hstack(*(c['COEFFICIENTS'][0] for c in chains)).evalf(40),dtype=complex)
            first = np.asarray(sp.Matrix.hstack(*(c['COEFFICIENTS'][1] for c in chains)).evalf(40),dtype=complex)
            gauge = np.linalg.pinv(zeroth)
            bases[side] = (zeroth,gauge,first-zeroth@gauge@first)
        self.emit(tag+'_LOCAL_MODE_DOMAIN',{'DEFINED':True,'SVD_RELATIVE_TOLERANCE':1e-10,
            'MATRIX_DECIMAL_DIGITS':(40,60),'SVD_FLOAT_MANTISSA_BITS':np.finfo(float).nmant+1,
            'COEFFICIENT_GAUGE':'UNIT_FRAME_LEFT_INVERSE_OF_THRESHOLD_BASIS'})
        evaluator = sp.lambdify((self.w,self.k,self.q),data['PHYSICAL'],'numpy',cse=True)
        denominator_evaluator = sp.lambdify((self.w,self.k,self.q),data['ROW_DENOMINATORS'],'numpy',cse=True)
        def node(name,wv,kv,qv,exact=False):
            arguments = {self.w:wv,self.k:kv,self.q:qv}
            if exact:
                operand = data['PHYSICAL'].subs(arguments)
                refined = operand.evalf(60); coarse = operand.evalf(40)
                matrix = np.asarray(refined,dtype=complex)
                self.emit(name+'_MATRIX_REFINEMENT',refined-coarse,lambda p:self.spectrum.strong_units[p])
            else:
                matrix = np.asarray(evaluator(complex(wv),complex(kv),complex(qv)),dtype=complex)
            denominators = np.asarray(denominator_evaluator(complex(wv),complex(kv),complex(qv)),dtype=complex)
            finite = bool(np.all(np.isfinite(matrix)) and np.all(np.isfinite(denominators)) and np.all(denominators!=0))
            self.emit(name+'_POINT',{'OMEGA':wv,'K':kv,'Q':qv},lambda p:
                d.measure(m.r.omega) if p[0]=='OMEGA' else d.measure(m.k) if p[0]=='K' else d.measure(m.q))
            self.emit(name+'_RADICAL_RESIDUAL',sp.N(data['RELATION'].subs(arguments),25),
                      lambda p:tuple(2*v for v in d.measure(m.q)))
            self.emit(name+'_NORMAL_SQUARE_RESIDUAL',sp.N(kv**2-square.subs(self.w,wv),25),
                      lambda p:tuple(2*v for v in d.measure(m.k)))
            self.emit(name+'_DENOMINATOR_COEFFICIENT_VALUES',sp.ImmutableMatrix(denominators))
            self.emit(name+'_FINITE_DOMAIN',finite)
            if not finite: return {'DEFINED':False}
            self.emit(name+'_PHYSICAL_MATRIX',sp.ImmutableMatrix(matrix),lambda p:self.spectrum.strong_units[p],heavy=True)
            u,s,vh = np.linalg.svd(matrix)
            tolerance = 1e-10*max(s)
            rank = int(np.sum(s>tolerance)); nullity = matrix.shape[0]-rank
            result = {'DEFINED':True,'RANK':rank,'NULLITY':nullity}
            self.emit(name+'_UNIT_FRAME_SINGULAR_VALUES',sp.ImmutableMatrix(s))
            self.emit(name+'_RANK_DOMAIN',{'RANK':rank,'NULLITY':nullity,'SVD_THRESHOLD':tolerance})
            for side,raw,operator in (('RIGHT',vh.conj().T[:,rank:],matrix),('LEFT',u[:,rank:],matrix.conj().T)):
                base,gauge,first = bases[side]
                overlap = gauge@raw
                overlap_rank = np.linalg.matrix_rank(overlap,tol=1e-10)
                defined = nullity==base.shape[1] and overlap_rank==base.shape[1]
                self.emit(name+'_'+side+'_GAUGE_DOMAIN',{'DEFINED':bool(defined),'OVERLAP_RANK':int(overlap_rank),
                    'THRESHOLD_SPACE_DIMENSION':base.shape[1]})
                if not defined: result['DEFINED']=False; continue
                basis = raw@np.linalg.inv(overlap)
                residual = operator@basis
                units = self.spectrum.field_units if side=='RIGHT' else tuple(tuple(-v for v in row) for row in self.spectrum.row_units)
                row_units = self.spectrum.row_units if side=='RIGHT' else tuple(tuple(-v for v in row) for row in self.spectrum.field_units)
                self.emit(name+'_'+side+'_BASIS',sp.ImmutableMatrix(basis),lambda p:units[p[0]//nullity],heavy=True)
                self.emit(name+'_'+side+'_EQUATION_RESIDUAL',sp.ImmutableMatrix(residual),lambda p:row_units[p[0]//nullity])
                self.emit(name+'_'+side+'_GAUGE_RESIDUAL',sp.ImmutableMatrix(gauge@basis-np.eye(nullity)))
                self.emit(name+'_'+side+'_THRESHOLD_BASIS_DIFFERENCE',sp.ImmutableMatrix(basis-base),
                          lambda p:units[p[0]//nullity],heavy=True)
                result[side] = basis
            return result
        approach = []
        local_radius = disk['RADIUS']/2 if disk['DEFINED'] else frequency/100
        self.emit(tag+'_LOCAL_RADIUS',local_radius,lambda p:d.measure(m.r.omega))
        for index,divisor in enumerate((1,100)):
            radius = local_radius/divisor
            for direction in (-1,1):
                wv = frequency+direction*radius
                lifts = sp.solve(self.k**2-square.subs(self.w,wv),self.k)
                pair = []
                for index_k,kv in enumerate(lifts):
                    name = tag+'_APPROACH_'+str(index)+'_'+('BELOW' if direction<0 else 'ABOVE')+'_'+str(index_k)
                    qroots = sp.solve(data['RELATION'].subs({self.w:wv,self.k:kv}),self.q)
                    qv = min(qroots,key=lambda root:abs(complex(sp.N(root-radical,30))))
                    trace = self.transport.trace([(complex(frequency),complex(normal)),(complex(wv),complex(kv))],seed=complex(radical))
                    self.emit(name+'_BULK_PATH',trace,self.path_unit)
                    if trace['PATH_DEFINED']:
                        self.emit(name+'_BULK_PATH_ENDPOINT_RESIDUAL',trace['END_Q']-complex(qv),lambda p:d.measure(m.q))
                    pair.append(node(name,wv,kv,qv,exact=True))
                defined = len(pair)==2 and all(v['DEFINED'] for v in pair)
                approach.append(defined)
                self.emit(tag+'_APPROACH_'+str(index)+'_'+str(direction).replace('-','M')+'_DOMAIN',defined)
                if defined:
                    for side in ('RIGHT','LEFT'):
                        denominator = complex(lifts[1]-lifts[0])
                        if side=='LEFT': denominator = denominator.conjugate()
                        secant = (pair[1][side]-pair[0][side])/denominator
                        units = self.spectrum.field_units if side=='RIGHT' else tuple(tuple(-v for v in row) for row in self.spectrum.row_units)
                        count = bases[side][0].shape[1]
                        for label,value in (('SECANT',secant),('CHAIN_OPERAND',bases[side][2]),('RESIDUAL',secant-bases[side][2])):
                            self.emit(tag+'_APPROACH_'+str(index)+'_'+str(direction).replace('-','M')+'_'+side+'_'+label,
                                sp.ImmutableMatrix(value),lambda p:tuple(a-b for a,b in zip(units[p[0]//count],d.measure(m.k))),
                                heavy=label!='RESIDUAL')
        paths = []
        dummy,normal_root = sp.symbols('s11cdThresholdPathFixedCoordinate s11cdThresholdPathNormalRoot')
        normal_transport = JointBulkSheetPath(normal_root**2-square,self.w,dummy,normal_root,sp.sqrt(square))
        radius = float(local_radius)
        for label,angle in (('UPPER',np.pi),('LOWER',-np.pi),('LOOP',2*np.pi)):
            frequencies = [complex(frequency)+radius*np.exp(1j*angle*j/16) for j in range(17)]
            frequencies[0] = complex(float(frequency)+radius)
            frequencies[-1] = frequencies[0] if label=='LOOP' else complex(float(frequency)-radius)
            for direction in (-1,1):
                name = tag+'_'+label+'_'+str(direction).replace('-','M')
                root = direction*complex(sp.N(sp.sqrt(square.subs(self.w,frequencies[0])),30))
                normals = [root]
                for wv in frequencies[1:]:
                    candidate = complex(square.subs(self.w,wv))**0.5
                    root = min((candidate,-candidate),key=lambda value:abs(value-root))
                    normals.append(root)
                normal_path = normal_transport.trace([(wv,0) for wv in frequencies],seed=normals[0])
                # This auxiliary trace transports normal momentum as its root;
                # all momentum/root entries therefore carry the normal unit.
                def normal_path_unit(p):
                    if p[-1] in ('SEED_Q','END_Q','ODE_END_Q','REFINEMENT_DIFFERENCE','ODE_DIFFERENCE','ODE_ABSOLUTE_TOLERANCE'):
                        return d.measure(m.k)
                    if p[-1] in ('SEED_RADICAL_RESIDUAL','MAXIMUM_RADICAL_RESIDUAL','ODE_MAXIMUM_RADICAL_RESIDUAL'):
                        return tuple(2*v for v in d.measure(m.k))
                    return self.path_unit(p)
                self.emit(name+'_NORMAL_PATH',normal_path,normal_path_unit)
                square0 = complex(data['RELATION'].subs({self.w:frequencies[0],self.k:normals[0],self.q:0}))
                qcandidate = (-square0)**0.5
                qstart = min((qcandidate,-qcandidate),key=lambda value:abs(value-complex(radical)))
                initial = self.transport.trace([(complex(frequency),complex(normal)),(frequencies[0],normals[0])],seed=complex(radical))
                self.emit(name+'_BULK_SEED_PATH',initial,self.path_unit)
                bulk = self.transport.trace(list(zip(frequencies,normals)),seed=qstart)
                self.emit(name+'_BULK_PATH',bulk,self.path_unit)
                if initial['PATH_DEFINED']:
                    self.emit(name+'_BULK_SEED_RESIDUAL',initial['END_Q']-qstart,lambda p:d.measure(m.q))
                if normal_path['PATH_DEFINED']:
                    self.emit(name+'_NORMAL_ENDPOINT_RESIDUAL',normal_path['END_Q']-normals[-1],lambda p:d.measure(m.k))
                defined = bool(initial['PATH_DEFINED'] and normal_path['PATH_DEFINED'] and bulk['PATH_DEFINED'])
                if defined:
                    root = qstart
                    node_records = []
                    for node_index,(wv,kv) in enumerate(zip(frequencies,normals)):
                        radicand = -complex(data['RELATION'].subs({self.w:wv,self.k:kv,self.q:0}))
                        candidate = radicand**0.5
                        root = min((candidate,-candidate),key=lambda value:abs(value-root))
                        evaluated = node(name+'_NODE_'+str(node_index),m.number(wv),m.number(kv),m.number(root))
                        node_records.append({key:evaluated[key] for key in ('DEFINED','RANK','NULLITY') if key in evaluated})
                    self.emit(name+'_NODE_DOMAINS',node_records)
                    self.emit(name+'_NODE_BULK_ENDPOINT_RESIDUAL',root-bulk['END_Q'],lambda p:d.measure(m.q))
                    defined = defined and all(v['DEFINED'] for v in node_records)
                paths.append({'LABEL':label,'INITIAL_NORMAL_SIGN':direction,'DEFINED':defined})
        summary = {'DEFINED':all(approach) and all(v['DEFINED'] for v in paths),
            'APPROACH_PAIR_COUNT':len(approach),'APPROACH_DEFINED_COUNT':sum(approach),'PATHS':paths,
            'EXACT_MODE_EXCEPTION_DISK_DEFINED':disk['DEFINED'],
            'AFFINE_BULK_SEGMENTS_BETWEEN_MODE_NODES':True,'GLOBAL_SHEET_ATLAS_COMPUTED':False,
            'PHYSICAL_CURRENT_CHANNEL_LABEL_COMPUTED':False}
        self.emit(tag+'_CONNECTION_SUMMARY',summary)
        return summary

    def exception_disk(self,data,tag,frequency):
        # Exact rational Taylor dominance counts roots of each already-derived
        # exception polynomial in a complex frequency disk. The target's
        # algebraic multiplicity is obtained by polynomial division, not a
        # witness. This keeps additional exceptional subloci out of the local
        # mode neighborhood when the emitted inequalities certify it.
        center = sp.Rational(str(sp.N(frequency,20)))
        initial_radius = sp.Rational(str(sp.N(frequency/50,15)))
        minimal = sp.Poly(sp.minpoly(frequency,self.w),self.w,domain=sp.QQ_I)
        operands = []
        for label,expression in data['LOCI'].items():
            poly = sp.Poly(expression,self.w,domain=sp.QQ_I)
            if poly.is_zero:
                operands.append((label,poly,None,None));continue
            quotient,multiplicity = poly,0
            while quotient.degree()>=minimal.degree():
                divided,remainder = quotient.div(minimal)
                if not remainder.is_zero:break
                quotient=divided;multiplicity+=1
            shifted = poly.shift(center)
            operands.append((label,poly,multiplicity,shifted))
        attempts=[]
        for attempt in range(4):
            radius = initial_radius/10**attempt
            inequalities=[]
            for label,poly,multiplicity,shifted in operands:
                if multiplicity is None:
                    inequalities.append({'LABEL':label,'DEFINED':False,'STATUS':'IDENTICALLY_ZERO_CONDITION'});continue
                lower = EndSpectrumCoverage.absolute_bounds(shifted.nth(multiplicity))[0]*radius**multiplicity
                upper = sum(EndSpectrumCoverage.absolute_bounds(shifted.nth(j))[1]*radius**j
                            for j in range(shifted.degree()+1) if j!=multiplicity)
                sign = sp.sign(lower-upper)
                inequalities.append({'LABEL':label,'DEFINED':bool(sign>0),'TARGET_MULTIPLICITY':multiplicity,
                    'POLYNOMIAL_DEGREE':poly.degree(),'EXACT_DOMINANCE_SIGN':sign,
                    'LOWER_BOUND':sp.N(lower,25),'OTHER_TERMS_UPPER_BOUND':sp.N(upper,25),
                    'BOUND_OPERANDS_SHA256':hashlib.sha256(sp.srepr((lower,upper)).encode()).hexdigest(),
                    'TAYLOR_RECONSTRUCTION_RESIDUAL':(shifted.shift(-center)-poly).as_expr()})
            centered = bool(abs(frequency-center)<radius/2)
            defined = centered and all(v['DEFINED'] for v in inequalities)
            attempts.append({'ATTEMPT':attempt,'DEFINED':defined,'TARGET_IN_INNER_HALF_DISK':centered,
                'RADIUS':radius,'INEQUALITIES':inequalities})
            if defined:break
        self.emit(tag+'_EXCEPTION_DISK_CENTER',center,lambda p:self.d.measure(self.modes.r.omega))
        self.emit(tag+'_EXCEPTION_DISK_ATTEMPTS',attempts,lambda p:
            self.d.measure(self.modes.r.omega) if p[-1]=='RADIUS' else self.d.zero)
        result={'DEFINED':defined,'CENTER':center,'RADIUS':radius,'CONDITION_COUNT':len(operands)}
        self.emit(tag+'_EXCEPTION_DISK_DOMAIN',result,lambda p:
            self.d.measure(self.modes.r.omega) if p[-1] in ('CENTER','RADIUS') else self.d.zero)
        return result


class BulkExceptionalSlice(EndExceptionalSlice):
    """Independent bulk geometry on the bound-carrier frequency slice.

    The physical entry denominators and radical relation supply this family.
    End-mode polynomials enter only the separately emitted intersections.
    Real-axis cells and selected continuation paths do not define a global
    complex-sheet atlas or a physical bound-pole search.
    """

    def __init__(self, spectrum, bindings):
        super().__init__(spectrum)
        self.bindings = bindings
        self.r,self.dimensions=self.modes.r,self.d

    @staticmethod
    def real_normal_projection(relation,denominator,w,k,q):
        projection=sp.Poly(sp.resultant(relation,denominator,q),k,w,domain=sp.QQ_I)
        real=sp.Poly.from_dict({m:sp.re(c) for m,c in projection.terms()},k,w,domain=sp.QQ)
        imaginary=sp.Poly.from_dict({m:sp.im(c) for m,c in projection.terms()},k,w,domain=sp.QQ)
        common=sp.gcd(real,imaginary)
        if common.is_zero:
            return {'DEFINED':False,'PROJECTION':projection.as_expr(),'REAL':real.as_expr(),
                    'IMAGINARY':imaginary.as_expr(),'SHARED':common.as_expr()}
        real_quotient,imaginary_quotient=real.exquo(common),imaginary.exquo(common)
        isolated=(sp.gcd(real_quotient,imaginary_quotient).as_expr()
                  if real_quotient.is_zero or imaginary_quotient.is_zero else
                  sp.resultant(real_quotient.as_expr(),imaginary_quotient.as_expr(),k))
        shared=sp.Poly(common.as_expr(),k,domain=sp.QQ.poly_ring(w))
        content,factors=shared.sqf_list()
        loci={'DENOMINATOR_REAL_NORMAL_ISOLATED_PROJECTION':isolated,
              'DENOMINATOR_REAL_NORMAL_SHARED_CONTENT':content}
        for i,(factor,_) in enumerate(factors):
            name='DENOMINATOR_REAL_NORMAL_SHARED_'+str(i)
            loci[name+'_LEADING']=factor.LC()
            loci[name+'_DISCRIMINANT']=factor.discriminant()
            for j,(other,_) in enumerate(factors[:i]):
                loci['DENOMINATOR_REAL_NORMAL_SHARED_'+str(j)+'_'+str(i)+'_INTERSECTION']=sp.resultant(factor.as_expr(),other.as_expr(),k)
        return {'DEFINED':True,'PROJECTION':projection.as_expr(),'REAL':real.as_expr(),'IMAGINARY':imaginary.as_expr(),
            'SHARED':common.as_expr(),'QUOTIENTS':(real_quotient.as_expr(),imaginary_quotient.as_expr()),
            'RECONSTRUCTION_RESIDUALS':((real-common*real_quotient).as_expr(),(imaginary-common*imaginary_quotient).as_expr()),
            'SHARED_FACTORS':[(f.as_expr(),a) for f,a in factors],'LOCI':loci}

    @staticmethod
    @lru_cache(maxsize=12)
    def analyze(physical, relation, w, k, q):
        denominators = tuple(sp.denom(sp.cancel(v)) for v in physical)
        denominator = sp.lcm(denominators)
        branch_radical=sp.Poly(sp.diff(relation,q),q).monic()
        branch_radicals=sp.solve(branch_radical.as_expr(),q)
        if len(branch_radicals)!=1:
            raise NotImplementedError('bulk slice requires the computed single radical critical point')
        branch = sp.Poly(relation.subs(q,branch_radicals[0]),k)
        normal_zero = relation.subs(k,0)
        norm = sp.Poly(sp.resultant(relation,denominator,k),q,domain=sp.QQ_I.poly_ring(w))
        content,factors = norm.sqf_list() if not norm.is_zero else (sp.S.Zero,[])
        loci = {'BRANCH_NORMAL_LEADING':branch.LC(),
                'BRANCH_NORMAL_DISCRIMINANT':branch.discriminant(),
                'BRANCH_NORMAL_ZERO':branch.eval(0),
                'RADICAL_LEADING':sp.Poly(relation,q).LC(),
                'DENOMINATOR_CONTENT':content}
        for i,(factor,multiplicity) in enumerate(factors):
            name='DENOMINATOR_FACTOR_'+str(i)
            loci[name+'_LEADING']=factor.LC()
            loci[name+'_DISCRIMINANT']=factor.discriminant()
            loci[name+'_BRANCH']=sp.resultant(factor.as_expr(),branch_radical.as_expr(),q)
            loci[name+'_NORMAL_ZERO']=sp.resultant(factor.as_expr(),normal_zero,q)
            for j,(other,_) in enumerate(factors[:i]):
                loci['DENOMINATOR_FACTOR_'+str(j)+'_'+str(i)+'_INTERSECTION']=sp.resultant(factor.as_expr(),other.as_expr(),q)
        real_normal=BulkExceptionalSlice.real_normal_projection(relation,denominator,w,k,q)
        if real_normal['DEFINED']:loci.update(real_normal['LOCI'])
        else:loci['DENOMINATOR_REAL_NORMAL_IDENTICAL_PROJECTION']=real_normal['PROJECTION']
        conditions={}
        for label,expression in loci.items():
            real,imaginary,common=EndExceptionalSlice.real_condition(expression,w)
            conditions[label]={'REAL':real,'IMAGINARY':imaginary,'GCD':common,
                'INTERVALS':None if common.is_zero else sp.polys.polytools.intervals(common,eps=sp.Rational(1,10**30))}
        return {'RELATION':relation,'PHYSICAL':physical,'ENTRY_DENOMINATORS':denominators,
                'DENOMINATOR':denominator,'DENOMINATOR_NORM':norm,'DENOMINATOR_FACTORS':factors,
                'BRANCH':branch,'BRANCH_RADICAL':branch_radical.as_expr(),'NORMAL_ZERO':normal_zero,
                'REAL_NORMAL_PROJECTION':real_normal,'LOCI':loci,'REAL_CONDITIONS':conditions}

    @staticmethod
    def nonnegative_targets(conditions):
        result={}
        for label,condition in conditions.items():
            common=condition['GCD']
            if not common.is_zero:
                for root in common.sqf_part().real_roots():
                    if root>=0:result.setdefault(root,[]).append(label)
        return result

    @staticmethod
    def rational_between(lower,upper):
        if lower is None:return sp.floor(upper)-1
        if upper is None:return sp.ceiling(lower)+1
        scale=sp.S.One
        while True:
            candidate=(sp.floor(lower*scale)+1)/scale
            if lower<candidate<upper:return candidate
            scale*=2

    def point_unit(self,path):
        if path[-1]=='OMEGA':return self.d.measure(self.modes.r.omega)
        if path[-1]=='K':return self.d.measure(self.modes.k)
        if path[-1]=='Q':return self.d.measure(self.modes.q)
        return self.d.zero

    def construct(self,suffix,channel_input,*,reference=False,end_data=None):
        m,d=self.modes,self.d;w,k,q=self.w,self.k,self.q
        self.prefix='BULK_EXCEPTIONAL_SLICE_INPUT_'+suffix
        algebraic,relation,joins=m.analytic(self.spectrum.strong)
        mapping=channel_input.mapping(algebraic,relation,(m.r.omega,m.k,m.q,m.eta,m.sigma))
        origin={m.eta:sp.S.Zero,m.sigma:sp.S.Zero} if reference else channel_input.origin
        coordinates={m.r.omega:w,m.k:k,m.q:q}
        physical=algebraic.xreplace(mapping).subs(origin).xreplace(coordinates).applyfunc(sp.cancel)
        relation=relation.xreplace(mapping).xreplace(coordinates)
        self.emit('UNIT_FRAME',channel_input.frame)
        self.emit('INPUT_BINDING',channel_input.specification)
        self.emit('BOUND_CARRIERS',tuple((str(s),v) for s,v in mapping.items()),
                  lambda p:d.measure(next(s for s in mapping if str(s)==p[0])))
        self.emit('GRADE_ORIGIN',tuple((str(s),v) for s,v in origin.items()))
        self.emit('SOURCE_GRADE_SUPPORT',tuple((p,tuple(sorted(PHYSICAL_METADATA.coefficients(v)))) for p,v in leaves(algebraic)))
        self.emit('SPECTRAL_COORDINATE_UNITS',{'FREQUENCY':d.measure(m.r.omega),'NORMAL':d.measure(m.k),
            'RADICAL':d.measure(m.q),'COEFFICIENT_COORDINATES':(w,k,q)})
        self.emit('BRANCH_JOIN_RESIDUALS',joins)
        self.emit('PHYSICAL_MATRIX_COEFFICIENTS',physical,lambda p:self.spectrum.strong_units[p],heavy=True)
        self.emit('RADICAL_RELATION_COEFFICIENTS',relation,heavy=True)
        data=self.analyze(physical,relation,w,k,q)
        self.emit('GEOMETRY_OPERANDS',{key:data[key] for key in ('ENTRY_DENOMINATORS','DENOMINATOR','NORMAL_ZERO','BRANCH_RADICAL')},heavy=True)
        self.emit('BRANCH_POLYNOMIAL',data['BRANCH'].as_expr(),heavy=True)
        self.emit('DENOMINATOR_NORM',data['DENOMINATOR_NORM'].as_expr(),heavy=True)
        self.emit('DENOMINATOR_FACTORS',[(f.as_expr(),a) for f,a in data['DENOMINATOR_FACTORS']],heavy=True)
        self.emit('DENOMINATOR_FACTOR_DEGREES_AND_MULTIPLICITIES',[(f.degree(),a) for f,a in data['DENOMINATOR_FACTORS']])
        self.emit('REAL_NORMAL_PROJECTION',data['REAL_NORMAL_PROJECTION'],heavy=True)
        self.emit('REAL_NORMAL_PROJECTION_DOMAIN',{'DEFINED':data['REAL_NORMAL_PROJECTION']['DEFINED']})
        if data['REAL_NORMAL_PROJECTION']['DEFINED']:
            self.emit('REAL_NORMAL_PROJECTION_RECONSTRUCTION_RESIDUALS',data['REAL_NORMAL_PROJECTION']['RECONSTRUCTION_RESIDUALS'])
        targets=self.nonnegative_targets(data['REAL_CONDITIONS'])
        for label,condition in data['REAL_CONDITIONS'].items():
            self.emit('LOCUS_'+label,{'ELIMINATION_POLYNOMIAL':data['LOCI'][label],
                'REAL_COEFFICIENT_POLYNOMIAL':condition['REAL'].as_expr(),
                'IMAGINARY_COEFFICIENT_POLYNOMIAL':condition['IMAGINARY'].as_expr(),
                'REAL_LOCUS_GCD':condition['GCD'].as_expr()},coefficient=True)
            self.emit('LOCUS_'+label+'_REAL_ISOLATION',{'IDENTICALLY_ZERO':condition['GCD'].is_zero,
                'GCD_DEGREE':condition['GCD'].degree(),'REAL_INTERVALS':condition['INTERVALS']
                if condition['INTERVALS'] is not None else 'UNRESOLVED_IDENTICAL_LOCUS'})
        end_defined=end_data is not None and end_data['DEFINED']
        end_targets=self.nonnegative_targets(end_data['REAL_CONDITIONS']) if end_defined else {}
        intersections=[]
        if end_defined:
            for bulk_label,bulk_condition in data['REAL_CONDITIONS'].items():
                for end_label,end_condition in end_data['REAL_CONDITIONS'].items():
                    common=sp.gcd(bulk_condition['GCD'],end_condition['GCD'])
                    intersections.append((bulk_label,end_label,common.as_expr(),common.degree()))
        self.emit('END_FAMILY_INTERSECTIONS',intersections,heavy=True)
        self.emit('FAMILY_TARGETS',{'BULK':tuple((root,labels) for root,labels in sorted(targets.items(),key=lambda v:float(v[0]))),
            'END_MODE':tuple((root,labels) for root,labels in sorted(end_targets.items(),key=lambda v:float(v[0])))},
            lambda p:d.measure(m.r.omega) if p[-1]==0 and len(p)==3 else d.zero)
        target_records=[]
        for i,(frequency,labels) in enumerate(sorted(targets.items(),key=lambda v:float(v[0]))):
            tag='TARGET_'+str(i)
            self.emit(tag+'_FREQUENCY',frequency,lambda p:d.measure(m.r.omega))
            self.emit(tag+'_INCIDENT_LOCUS_LABELS',labels)
            self.emit(tag+'_LOCUS_SUBSTITUTION_RESIDUALS',[(label,sp.simplify(data['LOCI'][label].subs(w,frequency))) for label in labels])
            if end_defined:
                self.emit(tag+'_END_LOCUS_INCIDENCE',[(label,sp.simplify(expression.subs(w,frequency))==0)
                    for label,expression in end_data['LOCI'].items()])
            components=[('BRANCH',sp.Poly(data['BRANCH_RADICAL'].subs(w,frequency),q,extension=True))]+[('DENOMINATOR_'+str(j),sp.Poly(f.as_expr().subs(w,frequency),q,extension=True))
                for j,(f,_) in enumerate(data['DENOMINATOR_FACTORS'])]
            for label,polynomial in components:
                name=tag+'_'+label
                self.emit(name+'_SPECIALIZED_POLYNOMIAL',polynomial.as_expr(),heavy=True)
                if polynomial.is_zero:
                    self.emit(name+'_STATUS','UNRESOLVED_IDENTICAL_COMPONENT');continue
                radicals=sp.solve(polynomial.as_expr(),q)
                self.emit(name+'_ROOT_CENSUS',{'DISTINCT_ROOT_COUNT':len(radicals),
                    'DEGREE_COUNT_RESIDUAL':polynomial.sqf_part().degree()-len(radicals)})
                for j,radical in enumerate(radicals):
                    lifts=sp.solve(relation.subs({w:frequency,q:radical}),k)
                    self.emit(name+'_RADICAL_'+str(j)+'_NORMAL_LIFT_COUNT',len(lifts))
                    for h,normal in enumerate(lifts):
                        point={w:frequency,k:normal,q:radical}
                        point_tag=name+'_POINT_'+str(j)+'_'+str(h)
                        values=tuple(sp.simplify(v.subs(point)) for v in data['ENTRY_DENOMINATORS'])
                        regular=all(v!=0 and not v.has(sp.nan,sp.zoo,sp.oo,-sp.oo) for v in values)
                        self.emit(point_tag,{'OMEGA':frequency,'K':normal,'Q':radical},self.point_unit)
                        self.emit(point_tag+'_NORMAL_REALITY_RESIDUAL',sp.simplify(sp.im(normal)),lambda p:d.measure(m.k))
                        self.emit(point_tag+'_RADICAL_RESIDUAL',sp.simplify(relation.subs(point)),lambda p:tuple(2*x for x in d.measure(m.q)))
                        self.emit(point_tag+'_DENOMINATOR_VALUES',values)
                        record={'FINITE_DENOMINATOR_DOMAIN':regular,'PHYSICAL_SHEET_MEMBERSHIP_COMPUTED':False}
                        if end_defined:
                            end_value=sp.simplify(end_data['POLYNOMIAL'].as_expr().subs(point))
                            self.emit(point_tag+'_END_POLYNOMIAL_VALUE',end_value,heavy=True)
                            record['END_POLYNOMIAL_VANISHES']=end_value==0
                        if regular:
                            matrix=physical.subs(point).applyfunc(sp.simplify)
                            self.emit(point_tag+'_PHYSICAL_MATRIX',matrix,lambda p:self.spectrum.strong_units[p],heavy=True)
                            rank=matrix.rank()
                            record.update(MATRIX_RANK=rank,MATRIX_NULLITY=matrix.cols-rank)
                        record['STATUS']='FINITE_BRANCH_POINT_MATRIX_COMPUTED' if regular else 'UNRESOLVED_SINGULAR_ENTRY_DENOMINATOR'
                        self.emit(point_tag+'_DOMAIN',record)
                        target_records.append(record)
        square=sp.solve(relation,q**2)[0]
        scale=sp.sqrt(-sp.Poly(square,k).nth(2))
        seed_operands=tuple(rhs.xreplace(dict(zip(lhs.args,(*m.r.tangents,m.k)))).xreplace(mapping).xreplace(coordinates)
            for lhs,rhs in self.bindings)
        seeds=tuple(scale*v for v in seed_operands)
        self.emit('SOURCE_BRANCH_SEED',seeds[0],lambda p:d.measure(m.q),heavy=True)
        self.emit('SOURCE_SEED_JOIN_RESIDUALS',tuple(sp.simplify(v-seeds[0]) for v in seeds),lambda p:d.measure(m.q))
        transport=JointBulkSheetPath(relation,w,k,q,seeds[0])
        boundaries=sorted(set((sp.S.Zero,*targets,*end_targets)),key=float)
        region_records=[]
        for i,lower in enumerate(boundaries):
            upper=boundaries[i+1] if i+1<len(boundaries) else None
            frequency=self.rational_between(lower,upper)
            self.emit('REGION_'+str(i)+'_FREQUENCY_INTERVAL',{'LOWER':lower,'UPPER':upper if upper is not None else 'UNBOUNDED',
                'WITNESS':frequency,'WITNESS_INSIDE':bool(lower<frequency and (upper is None or frequency<upper))},
                lambda p:d.zero if p[-1]=='WITNESS_INSIDE' else d.measure(m.r.omega))
            region_records.append(self.region(data,'REGION_'+str(i),frequency,transport))
        self.emit('SUMMARY',{'BULK_CRITICAL_FREQUENCY_COUNT':len(targets),'END_CRITICAL_FREQUENCY_COUNT':len(end_targets),
            'TARGET_POINT_COUNT':len(target_records),'FINITE_TARGET_POINT_COUNT':sum(v['FINITE_DENOMINATOR_DOMAIN'] for v in target_records),
            'FREQUENCY_REGION_COUNT':len(region_records),'REGIONS':region_records})
        self.emit('COVERAGE_BOUNDARIES',{'BOUND_PARAMETER_FREQUENCY_SLICE':True,'NONNEGATIVE_FREQUENCY_DOMAIN':True,
            'REAL_LOCUS_ENUMERATION_DEFINED':all(not c['GCD'].is_zero for c in data['REAL_CONDITIONS'].values()),
            'END_FAMILY_INTERSECTIONS_DEFINED':end_defined,'PARAMETER_VARIETY_ATLAS_COMPUTED':False,
            'GLOBAL_COMPLEX_SHEET_ATLAS_COMPUTED':False,'GENERALIZED_DEFECTIVE_MODES_COMPUTED':False,
            'TARGET_PHYSICAL_SHEET_MEMBERSHIP_COMPUTED':False,'PROFILE_FREQUENCY_BOUND_POLES_COMPUTED':False})
        return data

    def numerical_matrix(self,data,tag,point):
        import mpmath as mp
        w,k,q=self.w,self.k,self.q
        substituted=data['PHYSICAL'].subs(point)
        denominators=tuple(sp.N(v.subs(point),40) for v in data['ENTRY_DENOMINATORS'])
        finite=all(not v.has(sp.nan,sp.zoo,sp.oo,-sp.oo) and v!=0 for v in denominators)
        self.emit(tag+'_DENOMINATOR_VALUES',denominators)
        self.emit(tag+'_FINITE_DENOMINATOR_DOMAIN',finite)
        if not finite:
            self.emit(tag+'_STATUS','UNRESOLVED_SINGULAR_ENTRY_DENOMINATOR');return None
        values=[]
        for digits in (40,60):
            matrix=substituted.evalf(digits)
            with mp.workdps(digits):
                numeric=mp.matrix([[mp.mpc(str(sp.re(matrix[i,j])),str(sp.im(matrix[i,j]))) for j in range(matrix.cols)] for i in range(matrix.rows)])
                try:inverse=sp.ImmutableMatrix((numeric**-1).tolist())
                except ZeroDivisionError:inverse=None
            values.append((matrix,inverse))
        def inverse_unit(p):
            i,j=divmod(p[0],data['PHYSICAL'].cols)
            return tuple(a-b for a,b in zip(self.spectrum.field_units[i],self.spectrum.row_units[j]))
        self.emit(tag+'_PHYSICAL_MATRIX',values[1][0],lambda p:self.spectrum.strong_units[p],heavy=True)
        self.emit(tag+'_MATRIX_PRECISION_REFINEMENT',(values[1][0]-values[0][0]).evalf(50),lambda p:self.spectrum.strong_units[p])
        self.emit(tag+'_INVERSE_DOMAIN',{'DECIMAL_DIGITS':(40,60),'INVERSE_COMPUTED':all(v[1] is not None for v in values),
            'ARITHMETIC_SCOPE':'MATRIX_EVALUATION_AT_STORED_POINTS'})
        if any(v[1] is None for v in values):return None
        inverse=values[1][1]
        self.emit(tag+'_INVERSE',inverse,inverse_unit,heavy=True)
        self.emit(tag+'_INVERSE_PRECISION_REFINEMENT',(inverse-values[0][1]).evalf(50),inverse_unit)
        for label,residual,units in (
            ('RIGHT_INVERSE_RESIDUAL',values[1][0]*inverse-sp.eye(inverse.rows),self.spectrum.row_units),
            ('LEFT_INVERSE_RESIDUAL',inverse*values[1][0]-sp.eye(inverse.rows),self.spectrum.field_units)):
            self.emit(tag+'_'+label,residual.evalf(50),lambda p:tuple(a-b for a,b in zip(units[p[0]//inverse.rows],units[p[0]%inverse.rows])))
        return values[1][0],inverse

    def region(self,data,tag,frequency,transport):
        m,d=self.modes,self.d;w,k,q=self.w,self.k,self.q
        relation=data['RELATION'].subs(w,frequency)
        branch=sp.Poly(data['BRANCH'].as_expr().subs(w,frequency),k)
        projected=sp.resultant(relation,data['DENOMINATOR'].subs(w,frequency),q)
        real,imaginary,common=self.real_condition(projected,k)
        self.emit(tag+'_NORMAL_PROJECTION_OPERANDS',(branch.as_expr(),projected,real.as_expr(),imaginary.as_expr(),common.as_expr()),heavy=True)
        if common.is_zero:
            self.emit(tag+'_STATUS','UNRESOLVED_IDENTICAL_NORMAL_PROJECTION')
            return {'NORMAL_CELL_COUNT':0,'NORMAL_CELL_ENUMERATION_DEFINED':False,'BANK_PAIR_COUNT':0}
        real_boundaries=sorted(set((*branch.sqf_part().real_roots(),*common.sqf_part().real_roots())),key=float)
        self.emit(tag+'_REAL_NORMAL_BOUNDARIES',real_boundaries,lambda p:d.measure(m.k))
        endpoints=[None,*real_boundaries,None];cell_count=0
        for i,(lower,upper) in enumerate(zip(endpoints,endpoints[1:])):
            normal=sp.S.Zero if lower is None and upper is None else self.rational_between(lower,upper)
            self.emit(tag+'_CELL_'+str(i)+'_NORMAL_INTERVAL',{'LOWER':lower if lower is not None else 'UNBOUNDED',
                'UPPER':upper if upper is not None else 'UNBOUNDED','WITNESS':normal,
                'WITNESS_INSIDE':bool((lower is None or lower<normal) and (upper is None or normal<upper))},
                lambda p:d.zero if p[-1]=='WITNESS_INSIDE' else d.measure(m.k))
            radicals=sp.solve(relation.subs(k,normal),q)
            self.emit(tag+'_CELL_'+str(i)+'_RADICAL_ROOT_COUNT',len(radicals))
            for j,radical in enumerate(radicals):
                name=tag+'_CELL_'+str(i)+'_LIFT_'+str(j)
                point={w:frequency,k:normal,q:radical}
                self.emit(name+'_POINT',{'OMEGA':frequency,'K':normal,'Q':radical},self.point_unit)
                self.emit(name+'_RADICAL_RESIDUAL',sp.simplify(data['RELATION'].subs(point)),lambda p:tuple(2*v for v in d.measure(m.q)))
                self.numerical_matrix(data,name,point)
            cell_count+=1
        branch_points=sp.solve(branch.as_expr(),k)
        self.emit(tag+'_BRANCH_POINTS',branch_points,lambda p:d.measure(m.k))
        bank_count=0;path_statuses=Counter()
        for i,root in enumerate(branch_points):
            z=complex(root);scale=max(abs(complex(v)) for v in branch_points)
            target=1.5*z if z.imag else z+.5j*scale
            self.emit(tag+'_BANK_'+str(i)+'_TARGET',{'K':target,'OMEGA':frequency},self.point_unit)
            self.emit(tag+'_BANK_'+str(i)+'_OFFSETS',tuple(e*scale for e in (1e-4,1e-6)),lambda p:d.measure(m.k))
            previous=None
            for refinement,epsilon in enumerate((1e-4,1e-6)):
                bank=[]
                for side in (-1,1):
                    endpoint=target+side*epsilon*scale
                    path=transport.trace([(float(frequency),endpoint.real),(float(frequency),endpoint)])
                    path_statuses[path['STATUS']]+=1
                    name=tag+'_BANK_'+str(i)+'_'+str(refinement)+'_'+str(side).replace('-','M')
                    self.emit(name+'_PATH',path,lambda p:BulkContinuationAudit.path_unit(self,p))
                    value=None
                    if path['PATH_DEFINED']:
                        point={w:frequency,k:endpoint,q:path['END_Q']}
                        value=self.numerical_matrix(data,name,point)
                    bank.append((path,value))
                if all(v[1] is not None for v in bank):
                    jump=(bank[1][1][0]-bank[0][1][0]).evalf(50)
                    self.emit(tag+'_BANK_'+str(i)+'_'+str(refinement)+'_MATRIX_JUMP',jump,lambda p:self.spectrum.strong_units[p],heavy=True)
                    if previous is not None:
                        self.emit(tag+'_BANK_'+str(i)+'_'+str(refinement)+'_MATRIX_JUMP_REFINEMENT',jump-previous,
                                  lambda p:self.spectrum.strong_units[p],heavy=True)
                    previous=jump;bank_count+=1
        return {'NORMAL_CELL_COUNT':cell_count,'NORMAL_CELL_ENUMERATION_DEFINED':True,'BANK_PAIR_COUNT':bank_count,
                'BANK_PATH_STATUSES':dict(path_statuses)}


class BulkContinuationAudit:
    """Reduced-source joins, explicit radical paths and rational-matrix banks."""

    emit = EndSpectrumCoverage.emit

    def __init__(self, modes, strong, strong_units, branch_bindings):
        self.modes,self.r = modes,modes.r
        self.dimensions = PHYSICAL_METADATA.dimensions
        self.strong,self.strong_units,self.bindings = strong,strong_units,branch_bindings

    def path_unit(self, path):
        d=self.dimensions;m=self.modes
        if path[-1] in ('OMEGA','BRANCH_OMEGA'):return d.measure(self.r.omega)
        if path[-1]=='K':return d.measure(m.k)
        if path[-1] in ('SEED_Q','END_Q','ODE_END_Q','REFINEMENT_DIFFERENCE','ODE_DIFFERENCE','ODE_ABSOLUTE_TOLERANCE'):
            return d.measure(m.q)
        if path[-1] in ('SEED_RADICAL_RESIDUAL','MAXIMUM_RADICAL_RESIDUAL','ODE_MAXIMUM_RADICAL_RESIDUAL'):
            return tuple(2*v for v in d.measure(m.q))
        return d.zero

    def construct(self, suffix, *, channel_input=None, reference=False, sample_index=0, spectrum=None):
        m,d=self.modes,self.dimensions
        w,k,q=self.r.omega,m.k,m.q
        algebraic,relation,join=m.analytic(self.strong)
        if channel_input is None:
            mapping=m.sample(algebraic,relation,sample_index)
            frequency=mapping.pop(w)
            origin={m.eta:sp.S.Zero,m.sigma:sp.S.Zero}
            frame=('PIT_L','PIT_T','PIT_M')
        else:
            mapping=channel_input.mapping(algebraic,relation,(w,k,q,m.eta,m.sigma))
            frequency=channel_input.parameters['omega']
            origin={m.eta:sp.S.Zero,m.sigma:sp.S.Zero} if reference else channel_input.origin
            frame=channel_input.frame
        prefix='JOINT_SHEET_'+('INPUT_' if channel_input else 'PIT_')+suffix+'_'+str(sample_index)
        def output(label,value,unit=None,heavy=False):
            self.emit(prefix+'_'+label,value,unit,heavy=heavy)
        output('UNIT_FRAME',frame)
        output('BOUND_CARRIERS',[(str(s),v) for s,v in mapping.items()],
               lambda p:d.measure(next(s for s in mapping if str(s)==p[0])))
        output('GRADE_ORIGIN',[(str(s),v) for s,v in origin.items()])
        output('SOURCE_GRADE_SUPPORT',tuple((path,tuple(sorted(PHYSICAL_METADATA.coefficients(value))))
               for path,value in leaves(algebraic)))
        output('INPUT_BINDING',channel_input.specification if channel_input else {'PIT_SAMPLE':sample_index})
        output('OPERATOR_SCOPE',{'RETAINED_OPERATOR_EVALUATION':True,'CONTINUUM_REEXPANSION_PERFORMED':False,
            'FULL_RESOLVENT_CONTOUR_CONSTRUCTED':False,'CONTINUUM_MEASURE_CONSTRUCTED':False})
        original=self.strong.xreplace(mapping).subs(origin)
        bound=algebraic.xreplace(mapping).subs(origin).applyfunc(sp.cancel)
        relation=relation.xreplace(mapping)
        square=sp.solve(relation,q**2)[0]
        scale=sp.sqrt(-sp.Poly(square,k).nth(2))
        seed_operands=tuple(rhs.xreplace(dict(zip(lhs.args,(*self.r.tangents,k)))).xreplace(mapping)
                            for lhs,rhs in self.bindings)
        seeds=tuple(scale*value for value in seed_operands)
        seed_differences=tuple(sp.simplify(value-seeds[0]) for value in seeds)
        output('REDUCED_REAL_AXIS_BRANCH_OPERANDS',seed_operands,lambda p:d.measure(k))
        output('RADICAL_SCALE',scale,lambda p:tuple(a-b for a,b in zip(d.measure(q),d.measure(k))))
        output('SCALED_REAL_AXIS_SEED',seeds[0],lambda p:d.measure(q))
        output('SOURCE_BRANCH_JOIN_RESIDUAL',join)
        output('REDUCED_SEED_JOIN_RESIDUAL',seed_differences,lambda p:d.measure(q))
        output('RELATION',relation,lambda p:tuple(2*v for v in d.measure(q)))
        if any(v!=0 for v in join+seed_differences):
            raise NotImplementedError('joint-path reduced branch operands do not join')
        transport=JointBulkSheetPath(relation,w,k,q,seeds[0])
        d.known[transport.dw]=d.measure(w);d.known[transport.dk]=d.measure(k)
        output('IMPLICIT_DIFFERENTIAL',transport.derivative,lambda p:d.measure(q))
        k0=mapping[self.r.tangents[0]]
        # Polynomial roots retain complex loci even though the source Fourier
        # coordinates carry real-axis assumptions.
        frequency_points=sp.Poly(square.subs(k,k0),w).all_roots()
        momentum_points=sp.Poly(square.subs(w,frequency),k).all_roots()
        output('FREQUENCY_BRANCH_POINTS',frequency_points,lambda p:d.measure(w))
        output('MOMENTUM_BRANCH_POINTS',momentum_points,lambda p:d.measure(k))
        cone=max(frequency_points)
        w0,kbase=float(frequency),float(k0)
        evaluate=sp.lambdify((w,k,q),bound,'numpy',cse=True)
        denominators=tuple(sp.denom(v) for v in bound)
        evaluate_denominators=sp.lambdify((w,k,q),denominators,'numpy',cse=True)
        paths=[];bank_pairs=[];joins=[]
        def trace(label,vertices,seed=None):
            data=transport.trace(vertices,seed=seed)
            output(label,data,self.path_unit)
            paths.append((label,data))
            return data

        if reference:
            for i,factor in enumerate((sp.Rational(1,2),sp.Integer(2),sp.Rational(-1,2),sp.Integer(-2))):
                point={w:factor*cone,k:k0}
                root=sp.simplify(seeds[0].subs(point))
                source_operand=original.subs(point)
                continued_operand=bound.subs({**point,q:root})
                if channel_input is None:
                    # Algebraic PIT coefficients can generate large number
                    # fields. Evaluate both actual operands before comparison.
                    source=source_operand.evalf(40)
                    continued=continued_operand.evalf(40)
                    residual=source-continued
                    refined_source=source_operand.evalf(60)
                    refined_continued=continued_operand.evalf(60)
                    refined_residual=refined_source-refined_continued
                    scales=[max(sp.S.One,abs(a),abs(b)) for a,b in zip(refined_source,refined_continued)]
                    scaled=max(abs(v)/scale for v,scale in zip(refined_residual,scales))
                    output('REAL_AXIS_JOIN_PRECISION_'+str(i),{'DECIMAL_DIGITS':(40,60),
                        'MAXIMUM_SCALED_RESIDUAL':scaled,'ENTRY_COEFFICIENT_SCALE_FLOOR':1,
                        'MATCH_TOLERANCE':sp.Rational(1,10**30)})
                    output('REAL_AXIS_REFINED_MATRIX_RESIDUAL_'+str(i),refined_residual,lambda p:self.strong_units[p])
                    output('REAL_AXIS_SOURCE_REFINEMENT_'+str(i),refined_source-source,lambda p:self.strong_units[p])
                    output('REAL_AXIS_ALGEBRAIC_REFINEMENT_'+str(i),refined_continued-continued,lambda p:self.strong_units[p])
                    joined=bool(scaled<sp.Rational(1,10**30))
                else:
                    source=source_operand.applyfunc(sp.simplify)
                    continued=continued_operand.applyfunc(sp.simplify)
                    residual=(source-continued).applyfunc(sp.simplify)
                    output('REAL_AXIS_JOIN_PRECISION_'+str(i),{'ARITHMETIC':'EXACT'})
                    joined=all(v==0 for v in residual)
                output('REAL_AXIS_POINT_'+str(i),{'OMEGA':point[w],'K':k0,'END_Q':root},self.path_unit)
                output('REAL_AXIS_SOURCE_MATRIX_'+str(i),source,lambda p:self.strong_units[p],True)
                output('REAL_AXIS_ALGEBRAIC_MATRIX_'+str(i),continued,lambda p:self.strong_units[p],True)
                output('REAL_AXIS_MATRIX_RESIDUAL_'+str(i),residual,lambda p:self.strong_units[p])
                output('REAL_AXIS_RADICAL_RESIDUAL_'+str(i),sp.simplify(relation.subs({**point,q:root})),
                       lambda p:tuple(2*v for v in d.measure(q)))
                if not joined:raise NotImplementedError('real-axis source/algebraic matrix join differs')
                for direction in (-1,1):
                    trace('FREQUENCY_RAY_'+str(i)+'_'+str(direction).replace('-','M'),
                          [(complex(point[w]),kbase),(complex(point[w]+direction*sp.I*cone/5),kbase)])
            trace('FREQUENCY_BRANCH_INTERSECTION',[(2*float(cone),kbase),(float(cone),kbase)])
            points=[complex(v) for v in momentum_points]
            dw=1j*min(abs(w0-complex(v)) for v in frequency_points)/10
            dk=1j*min(abs(kbase-v) for v in points)/10
            a=trace('LOCAL_FREQUENCY_THEN_MOMENTUM',[(w0,kbase),(w0+dw,kbase),(w0+dw,kbase+dk)])
            b=trace('LOCAL_MOMENTUM_THEN_FREQUENCY',[(w0,kbase),(w0,kbase+dk),(w0+dw,kbase+dk)])
            if a['PATH_DEFINED'] and b['PATH_DEFINED']:
                output('LOCAL_PATH_ORDER_RESIDUAL',a['END_Q']-b['END_Q'],lambda p:d.measure(q))
            center=complex(cone);radius=float(cone)/3
            loop=[(center+radius*np.exp(2j*np.pi*j/32),kbase) for j in range(32)]
            loop.append(loop[0])
            first=trace('FREQUENCY_WINDING_ONCE',loop)
            if first['PATH_DEFINED']:
                second=trace('FREQUENCY_WINDING_TWICE',loop,seed=first['END_Q'])
                if second['PATH_DEFINED']:
                    output('WINDING_RETURN_OPERANDS',(first['SEED_Q'],first['END_Q'],second['END_Q']),lambda p:d.measure(q))
                    output('WINDING_DOUBLE_RETURN_RESIDUAL',second['END_Q']-first['SEED_Q'],lambda p:d.measure(q))

        def banks(label,coordinate,target,offset_scale):
            previous=None
            for refinement,epsilon in enumerate((1e-3,1e-5,1e-7)):
                bank=[]
                for side in (-1,1):
                    endpoint=target+side*epsilon*offset_scale
                    vertices=([(endpoint.real,kbase),(endpoint,kbase)] if coordinate=='OMEGA' else
                              [(w0,endpoint.real),(w0,endpoint)])
                    data=trace(label+'_'+str(refinement)+'_'+str(side).replace('-','M'),vertices)
                    record={'SIDE':side,'OFFSET':epsilon*offset_scale,'OMEGA':vertices[-1][0],'K':vertices[-1][1],
                            'PATH_DEFINED':data['PATH_DEFINED']}
                    if data['PATH_DEFINED']:
                        point=(*vertices[-1],data['END_Q'])
                        matrix=np.asarray(evaluate(*point),dtype=complex)
                        denominator_values=np.asarray(evaluate_denominators(*point),dtype=complex)
                        record.update({'END_Q':data['END_Q'],'FINITE_MATRIX':bool(np.isfinite(matrix).all()),
                                       'MINIMUM_DENOMINATOR_COEFFICIENT':float(np.min(np.abs(denominator_values)))})
                        output(label+'_MATRIX_'+str(refinement)+'_'+str(side).replace('-','M'),
                               sp.ImmutableMatrix(matrix),lambda p:self.strong_units[p],True)
                        output(label+'_DENOMINATORS_'+str(refinement)+'_'+str(side).replace('-','M'),
                               tuple(denominator_values))
                        bank.append((record,matrix))
                    else:bank.append((record,None))
                output(label+'_BANK_OPERANDS_'+str(refinement),[v[0] for v in bank],lambda p:
                       d.measure(w if coordinate=='OMEGA' else k) if p[-1]=='OFFSET' else self.path_unit(p))
                if all(v[0]['PATH_DEFINED'] and v[0]['FINITE_MATRIX'] for v in bank):
                    gap=bank[1][0]['END_Q']-bank[0][0]['END_Q']
                    jump=bank[1][1]-bank[0][1]
                    output(label+'_RADICAL_JUMP_'+str(refinement),gap,lambda p:d.measure(q))
                    output(label+'_MATRIX_JUMP_'+str(refinement),sp.ImmutableMatrix(jump),lambda p:self.strong_units[p],True)
                    if previous is not None:
                        output(label+'_RADICAL_REFINEMENT_'+str(refinement),
                               tuple(v[0]['END_Q']-old for v,old in zip(bank,previous)),lambda p:d.measure(q))
                    previous=[v[0]['END_Q'] for v in bank]
                    bank_pairs.append({'LABEL':label,'REFINEMENT':refinement,'RADICAL_JUMP':gap})

        if reference:
            for i,point in enumerate(frequency_points):
                banks('FREQUENCY_CUT_BANK_'+str(i),'OMEGA',complex(point-sp.I*cone/3),float(cone))
        unresolved=[] if spectrum is None else [v for v in spectrum['RECORDS']
            if str(v['FIXED_FREQUENCY_SHEET_MEMBERSHIP'])=='UNRESOLVED']
        targets=sorted({complex(v['K']) for v in unresolved},key=lambda z:(z.real,z.imag))
        target_source='UNRESOLVED_NATIVE_CANDIDATES'
        if not targets:
            target_source='COMPUTED_BRANCH_RAY_PROBES'
            targets=[1.5*complex(v) if complex(v).imag else complex(v)+.5j*max(abs(complex(t)) for t in momentum_points)
                     for v in momentum_points]
        output('MOMENTUM_BANK_TARGET_SOURCE',target_source)
        output('MOMENTUM_BANK_TARGETS',targets,lambda p:d.measure(k))
        output('ORIGINAL_UNRESOLVED_CANDIDATES',[(v['K'],v['Q']) for v in unresolved],
               lambda p:d.measure(k if p[-1]==0 else q))
        for i,target in enumerate(targets):
            if reference:trace('MOMENTUM_BRANCH_RAY_INTERSECTION_'+str(i),[(w0,target.real),(w0,target)])
            banks('MOMENTUM_CUT_BANK_'+str(i),'K',target,max(abs(complex(v)) for v in momentum_points))
        resolved=[v for _,v in paths if v['PATH_DEFINED']]
        summary={'PATH_COUNT':len(paths),'DEFINED_PATH_COUNT':len(resolved),
            'PATH_STATUSES':dict(Counter(v['STATUS'] for _,v in paths)),
            'BANK_PAIR_COUNT':len(bank_pairs),'ORIGINAL_UNRESOLVED_CANDIDATE_COUNT':len(unresolved),
            'CUT_ENCOUNTER_COUNT':sum(len(s['DOWNWARD_FREQUENCY_CUT_ENCOUNTERS']) for _,v in paths for s in v.get('SEGMENTS',())),
            'MAXIMUM_REFINEMENT_DIFFERENCE':max((abs(v['REFINEMENT_DIFFERENCE']) for v in resolved),default=sp.S.Zero),
            'MAXIMUM_ODE_DIFFERENCE':max((abs(v['ODE_DIFFERENCE']) for v in resolved),default=sp.S.Zero),
            'MAXIMUM_ODE_RADICAL_RESIDUAL':max((v['ODE_MAXIMUM_RADICAL_RESIDUAL'] for v in resolved),default=sp.S.Zero)}
        output('SUMMARY',summary,lambda p:d.measure(q) if p[-1] in ('MAXIMUM_REFINEMENT_DIFFERENCE','MAXIMUM_ODE_DIFFERENCE') else
               tuple(2*v for v in d.measure(q)) if p[-1]=='MAXIMUM_ODE_RADICAL_RESIDUAL' else d.zero)
        return summary


class EndResolventAudit(BulkContinuationAudit):
    """Inverse end pencils, local Laurent residues and finite-offset banks.

    A normal-momentum Laurent ansatz supplies the derivative pairing. Cauchy
    quadrature is a separate inverse-matrix construction on its local radical
    lift. These are fixed-frequency end data, not profile-frequency poles.
    """

    def __init__(self, modes, strong, strong_units, branch_bindings):
        super().__init__(modes,strong,strong_units,branch_bindings)
        d=self.dimensions
        self.field_units=[d.known[sp.Function('s11cdReducedField'+name)]
                          for name in ('u1','u2','u3','theta','eW')]
        self.row_units=[tuple(a+b for a,b in zip(strong_units[(5*i,)],self.field_units[0])) for i in range(5)]
        self.bank_matrices={};self.bank_operands={};self.poles=[]

    def matrix_unit(self,kind):
        d=self.dimensions;kdim=d.measure(self.modes.k)
        if kind=='inverse':a,b,extra=self.field_units,self.row_units,d.zero
        elif kind=='residue':a,b,extra=self.field_units,self.row_units,kdim
        elif kind=='second_moment':a,b,extra=self.field_units,self.row_units,tuple(2*x for x in kdim)
        elif kind=='field_map':a,b,extra=self.field_units,self.field_units,d.zero
        elif kind=='source_map':a,b,extra=self.row_units,self.row_units,d.zero
        elif kind=='field_residue_map':a,b,extra=self.field_units,self.field_units,kdim
        elif kind=='source_residue_map':a,b,extra=self.row_units,self.row_units,kdim
        elif kind=='pencil_derivative':a,b,extra=self.row_units,self.field_units,tuple(-x for x in kdim)
        else:raise ValueError(('end resolvent matrix unit',kind))
        return lambda p:tuple(x-y+z for x,y,z in zip(a[p[0]//5],b[p[0]%5],extra))

    def output(self,label,value,unit=None,heavy=False):
        if not (heavy and isinstance(value,sp.MatrixBase)):
            return EndSpectrumCoverage.emit(self,self.prefix+'_'+label,value,unit,heavy=heavy)
        # Tensor axes restore every entry's dimension without repeating a
        # leaf-path table for each contour matrix. Grades are computed for the
        # evaluated tensor; the unbound source support is emitted separately.
        unit=unit or (lambda p:self.dimensions.zero)
        body=cas(value);rows,columns=value.shape
        numeric=not (body.free_symbols-set(PHYSICAL_METADATA.generators))
        payload=self.modes.compact_fingerprint(body)if numeric else carrier_fingerprint(body)
        emit(self.prefix+'_'+label,payload)
        row_units=[tuple(unit((i*columns,)))for i in range(rows)]
        column_offsets=[tuple(a-b for a,b in zip(unit((j,)),row_units[0]))for j in range(columns)]
        encoding_residual=sorted({tuple(a-b-c for a,b,c in zip(unit((i*columns+j,)),row_units[i],column_offsets[j]))
                                 for i in range(rows)for j in range(columns)})
        original=self.modes.numeric_metadata(body,unit)
        grades=set();homotopy=set()
        for group in original:
            grades.update(tuple(v)for v in named(group,'MULTIGRADE'))
            homotopy.update(tuple(v)for v in named(group,'EPSILON_LAMBDA_SUPPORT'))
        emit('METADATA_'+self.prefix+'_'+label,{'REPRESENTATION':'MATRIX_AXIS_DIMENSIONS','SHAPE':(rows,columns),
            'ROW_DIMENSIONS_L_T_M':row_units,'COLUMN_DIMENSION_OFFSETS_L_T_M':column_offsets,
            'AXIS_ENCODING_RESIDUAL':encoding_residual,'MULTIGRADE':sorted(grades),'EPSILON_LAMBDA_SUPPORT':sorted(homotopy)})
        if any(any(v)for v in encoding_residual):raise ValueError('nonseparable matrix dimensions')

    def residual(self,label,value,kind):
        value=sp.ImmutableMatrix(value);unit=self.matrix_unit(kind)
        self.output(label,value,unit,True)
        groups={}
        for i,v in enumerate(value):
            dimension=unit((i,))
            groups[dimension]=max(groups.get(dimension,0.),abs(complex(v)))
        records=[{'DIMENSION_L_T_M':key,'MAXIMUM_ABSOLUTE':v} for key,v in sorted(groups.items())]
        self.output(label+'_MAXIMA',records,lambda p:records[p[0]]['DIMENSION_L_T_M']
                    if p[-1]=='MAXIMUM_ABSOLUTE' else self.dimensions.zero)
        return float(np.max(np.abs(np.asarray(value,dtype=complex))))

    def emit(self,tag,value,unit=None,*,heavy=False):
        # The parent emits its existing objects unchanged. Retain the actual
        # computed bank matrices and path endpoints as construction operands.
        super().emit(tag,value,unit,heavy=heavy)
        matrix=re.search(r'_((?:FREQUENCY|MOMENTUM)_CUT_BANK_\d+)_MATRIX_(\d+)_(M1|1)$',tag)
        bank=re.search(r'_((?:FREQUENCY|MOMENTUM)_CUT_BANK_\d+)_BANK_OPERANDS_(\d+)$',tag)
        if matrix:
            label,index,side=matrix.groups()
            self.bank_matrices[(label,int(index),-1 if side=='M1' else 1)]=np.asarray(value,dtype=complex)
        if bank:
            label,index=bank.groups();self.bank_operands[(label,int(index))]=value

    def local_circle(self,center,qcenter,radius,nodes):
        """Cauchy integrals from the circle ansatz and the computed inverse."""
        phase=np.exp(2j*np.pi*np.arange(nodes+1)/nodes)
        offsets=radius*phase;points=center+offsets
        seed=self.transport.trace([(self.frequency,center),(self.frequency,points[0])],seed=qcenter)
        if not seed['PATH_DEFINED']:return {'DEFINED':False,'STATUS':seed['STATUS'],'SEED_PATH':seed}
        root=complex(seed['END_Q']);radicals=[];matrices=[];inverses=[];derivatives=[]
        denominator_minimum=float('inf');condition=0.;inverse_error=0.;radical_error=0.
        for point in points:
            candidate=complex(self.square_evaluate(self.frequency,point))**.5
            root=min((candidate,-candidate),key=lambda v:abs(v-root))
            radicals.append(root)
            matrix=np.asarray(self.evaluate(self.frequency,point,root),dtype=complex)
            derivative=np.asarray(self.derivative_evaluate(self.frequency,point,root),dtype=complex)
            denominator=np.asarray(self.denominator_evaluate(self.frequency,point,root),dtype=complex)
            if not np.isfinite(matrix).all() or not np.isfinite(derivative).all():
                return {'DEFINED':False,'STATUS':'NONFINITE_CONTOUR_OPERATOR'}
            try:inverse=np.linalg.solve(matrix,np.eye(5))
            except np.linalg.LinAlgError:return {'DEFINED':False,'STATUS':'SINGULAR_CONTOUR_OPERATOR'}
            condition=max(condition,float(np.linalg.cond(matrix)))
            inverse_error=max(inverse_error,float(np.max(np.abs(matrix@inverse-np.eye(5)))),
                              float(np.max(np.abs(inverse@matrix-np.eye(5)))))
            radical_error=max(radical_error,abs(root*root-self.square_evaluate(self.frequency,point)))
            denominator_minimum=min(denominator_minimum,float(np.min(np.abs(denominator))))
            matrices.append(matrix);inverses.append(inverse);derivatives.append(derivative)
        inverses=np.asarray(inverses[:-1]);derivatives=np.asarray(derivatives[:-1])
        # The circle derivative is i*(k-center); the 2*pi/N quadrature and
        # 1/(2*pi*i) Cauchy measure produce these computed sample weights.
        angle=sp.Symbol('s11cdEndResolventContourAngle',real=True)
        displacement=sp.Symbol('s11cdEndResolventContourRadius',positive=True)
        circle=displacement*sp.exp(sp.I*angle)
        weight=sp.simplify(sp.diff(circle,angle)*(2*sp.pi/nodes)/(2*sp.pi*sp.I))
        weights=np.asarray(sp.lambdify((displacement,angle),weight,'numpy')(
            radius,2*np.pi*np.arange(nodes)/nodes),dtype=complex)
        residue=np.einsum('n,nij->ij',weights,inverses)
        moment=np.einsum('n,n,nij->ij',weights,offsets[:-1],inverses)
        projector=np.einsum('n,nij,njk->ik',weights,inverses,derivatives)
        determinants=np.asarray([np.linalg.det(v) for v in matrices])
        winding=float(np.sum(np.angle(determinants[1:]/determinants[:-1]))/(2*np.pi))
        return {'DEFINED':True,'STATUS':'CONTOUR_EVALUATED','RESIDUE':residue,
            'SECOND_LAURENT_MOMENT':moment,'DERIVATIVE_INTEGRAL':projector,
            'DERIVATIVE_TRACE':complex(np.trace(projector)),'DETERMINANT_WINDING':winding,
            'ROOT_CLOSURE_RESIDUAL':radicals[-1]-radicals[0],
            'MAXIMUM_RADICAL_RESIDUAL':radical_error,'MINIMUM_DENOMINATOR_COEFFICIENT':denominator_minimum,
            'MAXIMUM_COEFFICIENT_CONDITION':condition,'MAXIMUM_COEFFICIENT_INVERSE_RESIDUAL':inverse_error,
            'SEED_REFINEMENT_DIFFERENCE':seed['REFINEMENT_DIFFERENCE'],
            'SEED_ODE_DIFFERENCE':seed['ODE_DIFFERENCE'],'NODES':nodes,'RADIUS':radius}

    def pole(self,index,record,records,excluded):
        m,d=self.modes,self.dimensions
        center,qcenter=complex(record['K']),complex(record['Q'])
        tag='POLE_'+str(index)
        status={'K':record['K'],'Q':record['Q'],'ROOT_DISK_INDEX':record['ROOT_DISK_INDEX'],
            'ORIGINAL_SHEET_MEMBERSHIP':record['FIXED_FREQUENCY_SHEET_MEMBERSHIP'],
            'ALGEBRAIC_MULTIPLICITY':record['MULTIPLICITY'],'RESIDUE_DEFINED':False,
            'CONTOUR_DEFINED':False,'FLOAT_MANTISSA_BITS':np.finfo(float).nmant+1,
            'COEFFICIENT_SVD_RELATIVE_TOLERANCE':1e-8,'COORDINATE_RELATIVE_TOLERANCE':512*np.finfo(float).eps}
        unit=lambda p:d.measure(m.k) if p[-1] in ('K','RADIUS','MINIMUM_LOCUS_DISTANCE') else (
                      d.measure(m.q) if p[-1]=='Q' else d.zero)
        matrix=np.asarray(self.evaluate(self.frequency,center,qcenter),dtype=complex)
        derivative=np.asarray(self.derivative_evaluate(self.frequency,center,qcenter),dtype=complex)
        self.output(tag+'_ORIGINAL_OPERATOR',sp.ImmutableMatrix(matrix),lambda p:self.strong_units[p],True)
        self.output(tag+'_TOTAL_K_DERIVATIVE',sp.ImmutableMatrix(derivative),self.matrix_unit('pencil_derivative'),True)
        left,singular,right_h=np.linalg.svd(matrix)
        nullity=int(np.sum(singular<1e-8*max(1.,singular[0])))
        status['NULLITY']=nullity
        status['NULLITY_DIFFERENCE_FROM_NATIVE']=nullity-int(record.get('NULLITY',0))
        if nullity and qcenter and center:
            right,dual=right_h.conj().T[:,-nullity:],left[:,-nullity:]
            pairing=dual.conj().T@derivative@right
            rank=int(np.linalg.matrix_rank(pairing,tol=1e-10*max(1.,np.linalg.norm(pairing))))
            status['DERIVATIVE_PAIRING_RANK']=rank
            self.output(tag+'_RIGHT',sp.ImmutableMatrix(right),lambda p:self.field_units[p[0]//nullity],True)
            self.output(tag+'_LEFT',sp.ImmutableMatrix(dual),lambda p:tuple(-v for v in self.row_units[p[0]//nullity]),True)
            self.output(tag+'_DERIVATIVE_PAIRING',sp.ImmutableMatrix(pairing),lambda p:tuple(-v for v in d.measure(m.k)))
            if rank==nullity:
                # R=right*C is the Laurent ansatz. Projecting its constant
                # equation with the left nullspace computes the system for C.
                coefficients=np.linalg.solve(pairing,dual.conj().T)
                residue=right@coefficients
                projector=residue@derivative
                status['RESIDUE_DEFINED']=True
                self.output(tag+'_RESIDUE',sp.ImmutableMatrix(residue),self.matrix_unit('residue'),True)
                self.output(tag+'_DERIVATIVE_PROJECTOR',sp.ImmutableMatrix(projector),self.matrix_unit('field_map'),True)
                self.residual(tag+'_PROJECTOR_RESIDUAL',projector@projector-projector,'field_map')
                self.residual(tag+'_RIGHT_LAURENT_RESIDUAL',matrix@residue,'source_residue_map')
                self.residual(tag+'_LEFT_LAURENT_RESIDUAL',residue@matrix,'field_residue_map')
                self.poles.append({'index':index,'K':center,'Q':qcenter,'RESIDUE':residue,'PROJECTOR':projector})
            else:status['STATUS']='SINGULAR_DERIVATIVE_PAIRING'
        else:status['STATUS']='THRESHOLD_OR_EMPTY_NULLSPACE'
        scale=max(1.,abs(center));tolerance=status['COORDINATE_RELATIVE_TOLERANCE']*scale
        others=[];coincident=[]
        for j,other in enumerate(records):
            if j==index:continue
            distance=abs(complex(other['K'])-center)
            if distance<=tolerance:coincident.append({'INDEX':j,'K_DISTANCE':distance,'Q_DISTANCE':abs(complex(other['Q'])-qcenter)})
            else:others.append(distance)
        distances=[abs(z-center) for z in excluded]+others
        clearance=min(distances,default=0.)
        radius=clearance/5
        status.update({'MINIMUM_LOCUS_DISTANCE':clearance,'RADIUS':radius,'COINCIDENT_K_CANDIDATES':coincident})
        if radius<=tolerance or any(v['Q_DISTANCE']<=512*np.finfo(float).eps*max(1.,abs(qcenter)) for v in coincident):
            status['CONTOUR_STATUS']='UNRESOLVED_LOCAL_SEPARATION'
        else:
            runs=[]
            for radius_index,factor in enumerate((1.,.5)):
                for nodes in (32,64):
                    result=self.local_circle(center,qcenter,factor*radius,nodes)
                    matrices={key:result.pop(key) for key in ('RESIDUE','SECOND_LAURENT_MOMENT','DERIVATIVE_INTEGRAL') if key in result}
                    contour_tag=tag+'_CONTOUR_'+str(radius_index)+'_'+str(nodes)
                    self.output(contour_tag+'_RECORD',result,lambda p:
                        d.measure(m.k) if p[-1]=='RADIUS' else d.measure(m.q) if p[-1] in
                        ('ROOT_CLOSURE_RESIDUAL','SEED_REFINEMENT_DIFFERENCE','SEED_ODE_DIFFERENCE') else
                        tuple(2*v for v in d.measure(m.q)) if p[-1]=='MAXIMUM_RADICAL_RESIDUAL' else d.zero)
                    if result['DEFINED']:
                        self.output(contour_tag+'_RESIDUE',sp.ImmutableMatrix(matrices['RESIDUE']),self.matrix_unit('residue'),True)
                        self.output(contour_tag+'_DERIVATIVE_INTEGRAL',sp.ImmutableMatrix(matrices['DERIVATIVE_INTEGRAL']),self.matrix_unit('field_map'),True)
                        self.residual(contour_tag+'_SECOND_LAURENT_MOMENT',matrices['SECOND_LAURENT_MOMENT'],'second_moment')
                        self.output(contour_tag+'_POLE_COUNT_RESIDUAL',result['DERIVATIVE_TRACE']-int(record['MULTIPLICITY']))
                        self.output(contour_tag+'_WINDING_TRACE_RESIDUAL',result['DERIVATIVE_TRACE']-result['DETERMINANT_WINDING'])
                        if status['RESIDUE_DEFINED']:
                            self.residual(contour_tag+'_MODAL_RESIDUE_RESIDUAL',matrices['RESIDUE']-residue,'residue')
                            self.residual(contour_tag+'_MODAL_PROJECTOR_RESIDUAL',matrices['DERIVATIVE_INTEGRAL']-projector,'field_map')
                        runs.append((radius_index,nodes,result,matrices))
            status['CONTOUR_DEFINED']=len(runs)==4
            status['CONTOUR_STATUS']='FOUR_CONTOURS_EVALUATED' if len(runs)==4 else 'INCOMPLETE_CONTOUR_EVALUATION'
            if len(runs)==4:
                for label,a,b in [('NODE_REFINEMENT_OUTER',runs[0],runs[1]),('NODE_REFINEMENT_INNER',runs[2],runs[3]),
                                  ('RADIUS_REFINEMENT',runs[1],runs[3])]:
                    self.residual(tag+'_'+label,a[3]['RESIDUE']-b[3]['RESIDUE'],'residue')
            if status['RESIDUE_DEFINED']:status['STATUS']='REGULAR_LAURENT_RESIDUE'
        self.output(tag+'_RECORD',status,lambda p:d.measure(m.k) if p[-1]=='K_DISTANCE' else
                    d.measure(m.q) if p[-1]=='Q_DISTANCE' else unit(p))
        return status

    @staticmethod
    def mp_number(value,digits=65):
        import mpmath as mp
        real,imag=sp.N(sp.sympify(value),digits).as_real_imag()
        return mp.mpc(str(real),str(imag))

    def precise_pole(self,index,digits):
        import mpmath as mp
        key=(index,digits)
        if key in self.precise_poles:return self.precise_poles[key]
        native=self.native_records[index]
        point=(self.mp_number(self.exact_frequency),self.mp_number(native['K']),self.mp_number(native['Q']))
        matrix=mp.matrix(self.mp_evaluate(*point));derivative=mp.matrix(self.mp_derivative(*point))
        left,singular,right_h=mp.svd(matrix)
        nullity=sum(v<mp.mpf('1e-8')*max(1,singular[0])for v in singular)
        result={'DEFINED':False,'NULLITY':nullity,'DECIMAL_DIGITS':digits}
        if nullity:
            right,dual=right_h.H[:,matrix.cols-nullity:],left[:,matrix.rows-nullity:]
            pairing=dual.H*derivative*right
            _,s,_=mp.svd(pairing)
            rank=sum(v<mp.mpf('1e-30')*max(1,s[0])for v in s)
            result['PAIRING_NULLITY']=rank
            if rank==0:
                residue=right*(pairing**-1)*dual.H
                result.update({'DEFINED':True,'RESIDUE':residue,'K':point[1],'Q':point[2]})
                tag='POLE_'+str(index)+'_PRECISION_'+str(digits)
                self.output(tag+'_RESIDUE',sp.ImmutableMatrix(residue),self.matrix_unit('residue'),True)
                self.residual(tag+'_RIGHT_LAURENT_RESIDUAL',matrix*residue,'source_residue_map')
                self.residual(tag+'_LEFT_LAURENT_RESIDUAL',residue*matrix,'field_residue_map')
        self.output('POLE_'+str(index)+'_PRECISION_'+str(digits)+'_RECORD',
                    {k:v for k,v in result.items()if k not in ('RESIDUE','K','Q')})
        self.precise_poles[key]=result
        return result

    def inverse_banks(self):
        import mpmath as mp
        with mp.workdps(60):return self.precise_inverse_banks()

    def precise_inverse_banks(self):
        import mpmath as mp
        m,d=self.modes,self.dimensions
        summaries=[];previous={}
        for (label,index),operands in self.bank_operands.items():
            tag=label+'_'+str(index);inverses=[];matrices=[];regular=[];records=[];raw_inverses=[];raw_matrices=[]
            for operand in operands:
                side=int(operand['SIDE']);key=(label,index,side);suffix=tag+'_'+str(side).replace('-','M')
                record=dict(operand);record.update({'INVERSE_DEFINED':False,'SUBTRACTION_DEFINED':False})
                if key not in self.bank_matrices or not operand['PATH_DEFINED']:
                    records.append(record);continue
                raw_matrix=self.bank_matrices[key]
                try:raw_inverse=np.linalg.solve(raw_matrix,np.eye(5))
                except np.linalg.LinAlgError:record['DOUBLE_INVERSE_DEFINED']=False
                else:
                    record['DOUBLE_INVERSE_DEFINED']=bool(np.isfinite(raw_inverse).all())
                    record['DOUBLE_COEFFICIENT_CONDITION']=float(np.linalg.cond(raw_matrix))
                    self.output(suffix+'_DOUBLE_INVERSE',sp.ImmutableMatrix(raw_inverse),self.matrix_unit('inverse'),True)
                    raw_inverses.append(raw_inverse);raw_matrices.append(raw_matrix)
                candidates=[]
                if label.startswith('MOMENTUM'):
                    target=complex(operand['K'])-side*float(operand['OFFSET'])
                    center_root=complex(self.square_evaluate(self.frequency,target))**.5
                    center_root=min((center_root,-center_root),key=lambda z:abs(z-complex(operand['END_Q'])))
                    tolerance=512*np.finfo(float).eps
                    candidates=[(j,v)for j,v in enumerate(self.native_records)
                        if abs(complex(v['K'])-target)<=tolerance*max(1.,abs(target))
                        and abs(complex(v['Q'])-center_root)<=tolerance*max(1.,abs(center_root))]
                evaluated=[]
                for digits in (40,60):
                    with mp.workdps(digits):
                        omega=self.mp_number(self.exact_frequency if label.startswith('MOMENTUM') else operand['OMEGA'])
                        normal=(self.mp_number(candidates[0][1]['K'])+side*self.mp_number(operand['OFFSET'])
                                if candidates else self.mp_number(operand['K']))
                        radical=mp.sqrt(self.mp_square(omega,normal));old_q=self.mp_number(operand['END_Q'])
                        radical=min((radical,-radical),key=lambda v:abs(v-old_q))
                        matrix=mp.matrix(self.mp_evaluate(omega,normal,radical))
                        try:inverse=matrix**-1
                        except ZeroDivisionError:
                            self.output(suffix+'_PRECISION_'+str(digits)+'_DOMAIN',{'INVERSE_DEFINED':False})
                            continue
                        subtraction=mp.zeros(5);matches=[]
                        for candidate_index,_ in candidates:
                            pole=self.precise_pole(candidate_index,digits)
                            if pole['DEFINED']:
                                subtraction+=pole['RESIDUE']/(normal-pole['K']);matches.append(candidate_index)
                        evaluated.append((digits,matrix,inverse,subtraction,matches,omega,normal,radical))
                record['DECIMAL_PRECISIONS']=[v[0]for v in evaluated]
                record['NATIVE_CANDIDATES_AT_TARGET']=[j for j,_ in candidates]
                if len(evaluated)==2:
                    _,matrix,inverse,subtraction,matches,omega,normal,radical=evaluated[-1]
                    record.update({'INVERSE_DEFINED':True,'SUBTRACTION_DEFINED':len(matches)==len(candidates),
                        'SUBTRACTED_NORMAL_POLE_INDICES':matches,
                        'REFINED_COORDINATES':{'OMEGA':omega,'K':normal,'END_Q':radical}})
                    self.output(suffix+'_REFINED_OPERATOR',sp.ImmutableMatrix(matrix),lambda p:self.strong_units[p],True)
                    self.output(suffix+'_OPERATOR_PRECISION_DIFFERENCE',sp.ImmutableMatrix(matrix-mp.matrix(raw_matrix)),
                                lambda p:self.strong_units[p],True)
                    self.output(suffix+'_COORDINATE_PRECISION_DIFFERENCE',{'OMEGA':omega-self.mp_number(operand['OMEGA']),
                        'K':normal-self.mp_number(operand['K']),'END_Q':radical-self.mp_number(operand['END_Q'])},self.path_unit)
                    self.output(suffix+'_INVERSE',sp.ImmutableMatrix(inverse),self.matrix_unit('inverse'),True)
                    self.residual(suffix+'_INVERSE_PRECISION_REFINEMENT',inverse-evaluated[0][2],'inverse')
                    self.residual(suffix+'_LEFT_INVERSE_RESIDUAL',matrix*inverse-mp.eye(5),'source_map')
                    self.residual(suffix+'_RIGHT_INVERSE_RESIDUAL',inverse*matrix-mp.eye(5),'field_map')
                    if label.startswith('MOMENTUM') and record['SUBTRACTION_DEFINED']:
                        self.output(suffix+'_LOCAL_POLE_SUBTRACTION',sp.ImmutableMatrix(subtraction),self.matrix_unit('inverse'),True)
                        regular.append(inverse-subtraction)
                        self.output(suffix+'_REGULAR_PART',sp.ImmutableMatrix(regular[-1]),self.matrix_unit('inverse'),True)
                    else:record['SUBTRACTION_SCOPE']='NOT_APPLICABLE_FREQUENCY_BANK'
                    inverses.append(inverse);matrices.append(matrix)
                records.append(record)
            self.output(tag+'_BANK_RECORDS',records,lambda p:d.measure(self.r.omega if label.startswith('FREQUENCY') else m.k)
                        if p[-1]=='OFFSET' else self.path_unit(p))
            if len(raw_inverses)==2:
                self.residual(tag+'_DOUBLE_INVERSE_JUMP_RESIDUAL',raw_inverses[1]-raw_inverses[0]+
                    raw_inverses[1]@(raw_matrices[1]-raw_matrices[0])@raw_inverses[0],'inverse')
            if len(inverses)==2:
                jump=inverses[1]-inverses[0]
                identity_operand=-inverses[1]*(matrices[1]-matrices[0])*inverses[0]
                self.output(tag+'_INVERSE_JUMP',sp.ImmutableMatrix(jump),self.matrix_unit('inverse'),True)
                self.output(tag+'_INVERSE_JUMP_IDENTITY_OPERAND',sp.ImmutableMatrix(identity_operand),self.matrix_unit('inverse'),True)
                self.residual(tag+'_INVERSE_JUMP_RESIDUAL',jump-identity_operand,'inverse')
                regular_jump=regular[1]-regular[0]if len(regular)==2 else None
                if regular_jump is not None:self.output(tag+'_REGULAR_JUMP',sp.ImmutableMatrix(regular_jump),self.matrix_unit('inverse'),True)
                if label in previous:
                    self.residual(tag+'_INVERSE_JUMP_REFINEMENT',jump-previous[label][0],'inverse')
                    if regular_jump is not None and previous[label][1] is not None:
                        self.residual(tag+'_REGULAR_JUMP_REFINEMENT',regular_jump-previous[label][1],'inverse')
                previous[label]=(jump,regular_jump)
            summaries.append({'LABEL':label,'REFINEMENT':index,'INVERSE_COUNT':len(inverses),
                'SUBTRACTION_DEFINED':all(v['SUBTRACTION_DEFINED']for v in records),
                'SUBTRACTED_POLE_COUNT':sum(len(v.get('SUBTRACTED_NORMAL_POLE_INDICES',()))for v in records)})
        self.output('BANK_SUMMARY',summaries)
        return summaries

    def construct(self,suffix,*,channel_input=None,reference=False,sample_index=0,spectrum=None):
        m,d=self.modes,self.dimensions;w,k,q=self.r.omega,m.k,m.q
        self.prefix='END_RESOLVENT_'+('INPUT_' if channel_input else 'PIT_')+suffix+'_'+str(sample_index)
        if spectrum is not None and spectrum.get('DEFINED') is False:
            self.output('SPECTRUM_DOMAIN',{'DEFINED':False,'STATUS':spectrum['STATUS']})
            return
        self.bank_matrices={};self.bank_operands={};self.poles=[]
        self.native_records=spectrum['RECORDS'];self.precise_poles={}
        algebraic,relation,joins=m.analytic(self.strong)
        if channel_input is None:
            mapping=m.sample(algebraic,relation,sample_index);frequency=mapping.pop(w)
            origin={m.eta:sp.S.Zero,m.sigma:sp.S.Zero};frame=('PIT_L','PIT_T','PIT_M')
        else:
            mapping=channel_input.mapping(algebraic,relation,(w,k,q,m.eta,m.sigma))
            frequency=channel_input.parameters['omega']
            origin={m.eta:sp.S.Zero,m.sigma:sp.S.Zero} if reference else channel_input.origin;frame=channel_input.frame
        bound=algebraic.xreplace(mapping).subs(origin).applyfunc(sp.cancel)
        relation=relation.xreplace(mapping);self.frequency=float(frequency);self.exact_frequency=frequency
        square=sp.solve(relation,q**2)[0]
        scale=sp.sqrt(-sp.Poly(square,k).nth(2))
        seed=scale*self.bindings[0][1].xreplace(dict(zip(self.bindings[0][0].args,(*self.r.tangents,k)))).xreplace(mapping)
        self.transport=JointBulkSheetPath(relation,w,k,q,seed)
        derivative_variable=sp.Symbol('s11cdEndResolventRadicalMomentumDerivative')
        qderivative=sp.solve(sp.diff(relation,k)+sp.diff(relation,q)*derivative_variable,derivative_variable)[0]
        derivative=bound.diff(k)+bound.diff(q)*qderivative
        self.evaluate=sp.lambdify((w,k,q),bound,'numpy',cse=True)
        self.derivative_evaluate=sp.lambdify((w,k,q),derivative,'numpy',cse=True)
        self.square_evaluate=sp.lambdify((w,k),square,'numpy')
        self.mp_evaluate=sp.lambdify((w,k,q),bound,'mpmath',cse=True)
        self.mp_derivative=sp.lambdify((w,k,q),derivative,'mpmath',cse=True)
        self.mp_square=sp.lambdify((w,k),square,'mpmath')
        denominators=tuple(sp.lcm([sp.denom(bound[i,j])for j in range(5)])for i in range(5))
        self.denominator_evaluate=sp.lambdify((w,k,q),denominators,'numpy',cse=True)
        self.output('UNIT_FRAME',frame)
        self.output('BOUND_CARRIERS',[(str(s),v)for s,v in mapping.items()],lambda p:d.measure(next(s for s in mapping if str(s)==p[0])))
        self.output('GRADE_ORIGIN',[(str(s),v)for s,v in origin.items()])
        self.output('SOURCE_GRADE_SUPPORT',tuple((path,tuple(sorted(PHYSICAL_METADATA.coefficients(v))))for path,v in leaves(algebraic)))
        self.output('INPUT_BINDING',channel_input.specification if channel_input else {'PIT_SAMPLE':sample_index})
        self.output('FREQUENCY',frequency,lambda p:d.measure(w))
        self.output('COEFFICIENT_COORDINATE_UNITS',{'FIELDS':self.field_units,'EQUATION_ROWS':self.row_units,
                    'K':d.measure(k),'Q':d.measure(q)})
        self.output('SCOPE',{'CONSTANT_END_RESOLVENT':True,'NORMAL_MOMENTUM_LAURENT_DATA':True,
            'PROFILE_FREQUENCY_POLE_SOLVE':False,'FULL_PROFILE_RESOLVENT':False,'CONTINUUM_MEASURE':False,
            'FLUX_NORMALIZATION':False,'CONTINUUM_REEXPANSION_PERFORMED':False})
        self.output('SOURCE_BRANCH_JOIN_RESIDUAL',joins)
        self.output('OPERATOR',bound,lambda p:self.strong_units[p],True)
        self.output('TOTAL_K_DERIVATIVE',derivative,self.matrix_unit('pencil_derivative'),True)
        self.output('RADICAL_K_DERIVATIVE',qderivative,lambda p:tuple(a-b for a,b in zip(d.measure(q),d.measure(k))),True)
        self.output('RADICAL_DERIVATIVE_RESIDUAL',sp.simplify(sp.diff(relation,k)+sp.diff(relation,q)*qderivative),
                    lambda p:tuple(2*a-b for a,b in zip(d.measure(q),d.measure(k))))
        fixed=relation.subs(w,frequency)
        branch_points=[complex(v)for v in sp.Poly(square.subs(w,frequency),k).all_roots()]
        loci=[];denominator_points=[]
        for denominator in denominators:
            norm=sp.Poly(sp.resultant(denominator.subs(w,frequency),fixed,q),k)
            points=[]
            if norm.degree()>0:
                for factor,_ in sp.sqf_list(norm)[1]:points.extend(sp.nroots(factor,n=30,maxsteps=500))
            loci.append({'COEFFICIENT_POLYNOMIAL':norm.as_expr(),'ROOTS':points})
            denominator_points.extend(complex(v)for v in points)
        self.output('DENOMINATOR_LOCUS_OPERANDS',loci,lambda p:d.measure(k) if 'ROOTS'in p else d.zero,True)
        self.output('BRANCH_POINTS',branch_points,lambda p:d.measure(k))
        self.output('DENOMINATOR_POINTS',denominator_points,lambda p:d.measure(k))
        pole_records=[self.pole(i,v,self.native_records,branch_points+denominator_points)
                      for i,v in enumerate(self.native_records)]
        continuation=super().construct(suffix,channel_input=channel_input,reference=reference,
                                        sample_index=sample_index,spectrum=spectrum)
        banks=self.inverse_banks()
        summary={'NATIVE_CANDIDATE_COUNT':len(self.native_records),'RESIDUE_COUNT':len(self.poles),
            'CONTOUR_CANDIDATE_COUNT':sum(v['CONTOUR_DEFINED']for v in pole_records),
            'NULLITY_DIFFERENCE_COUNT':sum(v['NULLITY_DIFFERENCE_FROM_NATIVE']!=0 for v in pole_records),
            'POLE_STATUSES':dict(Counter(v.get('STATUS','UNRESOLVED')for v in pole_records)),
            'BANK_PAIR_COUNT':len(banks),'BANK_INVERSE_COUNT':sum(v['INVERSE_COUNT']for v in banks),
            'SUBTRACTION_UNRESOLVED_COUNT':sum(not v['SUBTRACTION_DEFINED']for v in banks),
            'ORIGINAL_UNRESOLVED_CANDIDATE_COUNT':continuation['ORIGINAL_UNRESOLVED_CANDIDATE_COUNT']}
        self.output('SUMMARY',summary)
        return summary


def run():
    global PHYSICAL_METADATA
    os.chdir(ROOT)
    faulthandler.register(signal.SIGUSR1)
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=['ALL']+[a+'__'+r for a in
        ('LAB_HELD', 'MATERIAL_ADVECTED') for r in ('RHO4_CONSTANT', 'RHOBR_CONSTANT')], default='ALL')
    parser.add_argument('--dev-symbol-cache', type=Path)
    parser.add_argument('--dev-stop-after-reduction', action='store_true',
                        help='emit a reduction-only development checkpoint')
    parser.add_argument('--dev-reduction-cache', type=Path,
                        help='save imported operands for focused reconstruction checks')
    channel_options = parser.add_mutually_exclusive_group()
    channel_options.add_argument('--channel-input-json',
                                 help='explicit unit_frame, parameters and independent w/m profile input')
    channel_options.add_argument('--channel-input-file', type=Path,
                                 help='JSON file containing the same explicit channel input')
    parser.add_argument('--channel-input-scope', choices=('modes-and-jets','spectrum'), default='modes-and-jets',
                        help='select explicit-input spectrum records with or without the additional legacy mode jets')
    options = parser.parse_args()
    selected = lambda case: options.case == 'ALL' or '__'.join(map(str, case)) == options.case
    fold, audit = load_model('scripts/S11c_b_exports.py', 'scripts/S11c_c1_exports.py', 'scripts/S11c_c2_exports.py')
    closure = check_consumer(fold, IMPORT_KEYS)
    witness = assert_lookups_equal_manifest(bind, fold, IMPORT_KEYS)
    rows = witness['result']
    emit('IMPORT_FOLD', audit)
    emit('IMPORT_LOOKUPS', sorted(witness['lookups']))
    emit('IMPORT_CLOSURE', {k: v for k, v in closure.items() if k != 'resolved_imports'})
    emit('BUILD_INPUT_DIGESTS', {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in BUILD_INPUT_PATHS})
    for row_key in CLOSED_KEYS:
        for case, payload in rows[row_key]['value']:
            if not selected(case):
                continue
            for slot, value in payload:
                name = '_'.join((row_key, *(str(v) for v in case), str(slot)))
                emit('ELEMENT_CENSUS_' + name, census(value))
    reduction = EdgeReduction(rows)
    dimensions = DimensionAnalysis(rows, reduction)
    PHYSICAL_METADATA = PhysicalMetadata(dimensions, reduction)
    if options.dev_reduction_cache:
        options.dev_reduction_cache.parent.mkdir(parents=True, exist_ok=True)
        with options.dev_reduction_cache.open('wb') as stream:
            pickle.dump(({k: dict(v) for k, v in rows.items()}, dimensions.known), stream)
    input_text = (options.channel_input_file.read_text() if options.channel_input_file
                  else options.channel_input_json)
    channel_input = ChannelInput(reduction, json.loads(input_text)) if input_text else None
    if channel_input is not None:
        channel_input.emit_inputs()
    carriers = FourierCarrierReconstruction(reduction)
    reconstruction = EdgeReconstruction(reduction, carriers)
    carriers.emit_kernel()
    physical('FOURIER_NORMALIZATION_COMPUTATION', reduction.normalization_operands)
    emit('INFERRED_INPUT_DIMENSIONS', {str(k): v for k, v in dimensions.known.items()})
    physical('ZERO_JET_DISTRIBUTIONAL_PRESCRIPTION', reduction.distribution_record)
    for profile in reduction.profiles:
        f, left, right, jump, localized = reduction.subtraction(profile)
        physical('ZERO_JET_SUBTRACTION_'+profile.upper(), (f, left, right, jump, localized))
        emit('HALF_LINE_TAIL_PREMISE_'+profile.upper(),
             sp.Lt(sp.Integral(sp.Abs(localized, evaluate=False),
                               (reduction.xi, -sp.oo, sp.oo)), sp.oo, evaluate=False))
    reduced_rows = {}
    reduced_branch_bindings = {}
    for row_key in CLOSED_KEYS:
        for case, payload in rows[row_key]['value']:
            if not selected(case):
                continue
            suffix = '_'.join((row_key, *(str(v) for v in case)))
            value = named(payload, 'VALUE')
            value_units = payload_units(value, named(payload, 'DIMENSION_L_T_M'))
            definitions = tuple(named(payload, 'FOURIER_PROFILE_BINDINGS'))
            for index, hat in enumerate(sorted(payload.atoms(AppliedUndef), key=sp.default_sort_key)):
                if hat.func.__name__.startswith('s11cc2Fourier'):
                    hat_reduction = reduction.hat(hat)
                    hat_unit = dimensions.measure(hat)
                    image_unit = tuple(a+2*b for a,b in zip(hat_unit,dimensions.measure(reduction.tangents[0])))
                    physical('PROFILE_CARRIER_REDUCTION_' + suffix + '_' + str(index),
                             (hat, hat_reduction), operands=(hat, hat_reduction[0]),
                             zero_dimensions={(0,):hat_unit,(1,):image_unit})
                    carriers.carrier(hat, definitions, suffix+'_'+str(index))
            tangent_map = {k: q for group in reduction.momentum_groups
                           for k, q in zip(group[:2], reduction.tangents)}
            for index, equation in enumerate(named(payload, 'COMPUTED_BRANCH_BINDINGS')):
                reduced_branch = equation.rhs.xreplace(tangent_map).xreplace(reduction.normal_map)
                physical('BRANCH_DEFINITION_REDUCTION_' + suffix + '_' + str(index),
                         (equation, reduced_branch), operands=(equation.lhs, equation.rhs, reduced_branch))
                reconstructed_branch = reduced_branch.xreplace({v:k for k,v in reduction.normal_map.items()})
                source_branch = equation.rhs.xreplace(tangent_map)
                physical('BRANCH_RECONSTRUCTION_OPERANDS_'+suffix+'_'+str(index),
                         (equation.rhs, source_branch, reconstructed_branch))
                physical('BRANCH_RECONSTRUCTION_RESIDUAL_'+suffix+'_'+str(index),
                         reconstructed_branch-source_branch,
                         zero_dimensions={(): dimensions.measure(equation.rhs)})
            for index, equation in enumerate(named(payload, 'FOURIER_PROFILE_BINDINGS')):
                reduced_definition, constraints = reduction.profile_definition(equation.rhs)
                source_unit = dimensions.measure(equation.rhs)
                image_unit = tuple(a+2*b for a,b in zip(source_unit,dimensions.measure(reduction.tangents[0])))
                physical('PROFILE_DEFINITION_REDUCTION_' + suffix + '_' + str(index),
                         (equation, reduced_definition, constraints),
                         operands=(equation.lhs, equation.rhs, reduced_definition),
                         zero_dimensions={(0,):source_unit,(1,):source_unit,(2,):image_unit})
                engine_coefficient = reduction.hat(equation.lhs)[0]
                tangent_delta = sp.prod(sp.DiracDelta(c) for c in constraints)
                reconstructed_definition = engine_coefficient*tangent_delta
                source_definition_on_class = reduced_definition*tangent_delta
                physical('PROFILE_DEFINITION_RECONSTRUCTION_OPERANDS_'+suffix+'_'+str(index),
                         (equation.rhs, engine_coefficient, reconstructed_definition, source_definition_on_class),
                         zero_dimensions={(0,):source_unit,(1,):image_unit,(2,):source_unit,(3,):source_unit})
                physical('PROFILE_DEFINITION_RECONSTRUCTION_RESIDUAL_'+suffix+'_'+str(index),
                         reconstructed_definition-source_definition_on_class,
                         zero_dimensions={(): dimensions.measure(equation.rhs)})
            reduced, records = reduction.value(value)
            for index, record in enumerate(records):
                physical('ACTION_INTEGRAL_REDUCTION_' + suffix + '_' + str(index), record,
                         operands=record[:2], zero_dimensions={(i,):dimensions.measure(record[0]) for i in range(2)})
            physical('REDUCED_ACTION_ROWS_' + suffix, reduced, zero_dimensions=value_units)
            reconstruction.row(value, reduced, records, suffix, value_units, definitions)
            reduced_payload = reduction.payload(payload, reduced)
            reduced_branch_bindings[(row_key,case)] = tuple((eq.lhs,eq.rhs)
                for eq in named(reduced_payload,'COMPUTED_BRANCH_BINDINGS'))
            physical('REDUCED_FIVE_SLOT_PAYLOAD_'+suffix, reduced_payload, operands=reduced, zero_dimensions=value_units)
            for slot in ('COMPUTED_BRANCH_BINDINGS','FOURIER_PROFILE_BINDINGS'):
                for index, equation in enumerate(named(reduced_payload,slot)):
                    physical('REDUCED_BINDING_OPERANDS_'+suffix+'_'+slot+'_'+str(index),
                             (equation.lhs,equation.rhs), zero_dimensions={
                                 (1,):dimensions.known[equation.lhs.func]})
            for slot, body in reduced_payload:
                emit('REDUCED_SLOT_CENSUS_'+suffix+'_'+str(slot), census(body))
            # The Abel regulator is removed weakly only after the complete
            # normal convolution. The body remains available to later solvers.
            fingerprinted('REDUCED_ACTION_WEAK_LIMIT_'+suffix, reduction.weak_limit(reduced), value_units)
            reduced_rows[(row_key, case)] = reduced
            emit('REDUCED_ACTION_CENSUS_' + suffix, census(reduced))
            old_coordinates = set((*reduction.x, *reduction.y, reduction.t,
                                   *(k for group in reduction.momentum_groups for k in group)))
            emit('ACTION_COORDINATE_CENSUS_' + suffix,
                 tuple(sorted(dag_free_symbols(reduced) & old_coordinates, key=sp.default_sort_key)))
    emit('REDUCED_DIMENSION_CONSTRAINT_RESIDUALS', sorted(dimensions.constraints, key=sp.default_sort_key))
    unresolved = tuple((str(a), ds) for a, ds in dimensions.unknown.items()
                       if a not in dimensions.known and any(d.free_symbols for d in ds))
    emit('REDUCED_DIMENSION_UNRESOLVED', unresolved)
    if dimensions.constraints or unresolved:
        raise ValueError('reduced dimensional analysis has surfaced unresolved constraints')
    if options.dev_stop_after_reduction:
        emit('DEVELOPMENT_STOP', 'AFTER_REDUCTION')
        emit('RESOURCE_MEASUREMENTS', (time.monotonic()-started,
                                      resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        emit('EMISSION_LINES', emission_index(EMISSION_LINES))
        emit('PROCESS_COMPLETION', sp.Integer(os.getpid()))
        return
    for case, _ in rows[CLOSED_KEYS[0]]['value']:
        if not selected(case):
            continue
        suffix = '_'.join(str(c) for c in case)
        pencil = ReducedPencil(reduced_rows[(CLOSED_KEYS[0], case)],
                               reduced_rows[(CLOSED_KEYS[1], case)], reduction)
        columns = pencil.columns()
        fingerprinted('FULL_REDUCED_PENCIL_ACTION_COLUMNS_'+suffix, columns)
        extracted = pencil.off_diagonal()
        fingerprinted('REDUCED_OPERATOR_OFF_DIAGONAL_OPERAND_'+suffix, extracted)
        fingerprinted('REDUCED_CANONICAL_KERNEL_OPERAND_'+suffix, pencil.kernel)
        residual = tree_difference(extracted, pencil.kernel)
        fingerprinted('REDUCED_OFF_DIAGONAL_EXTRACTION_RESIDUAL_'+suffix, residual)
        blocks = pencil.sector_blocks()
        block_units = {path: dimensions.measure(expr) for path, expr in leaves(blocks)}
        kernel_units = {path: dimensions.measure(expr) for path, expr in leaves(pencil.kernel)}
        strong_units = {path: dimensions.measure(expr) for path, expr in leaves(pencil.strong)}
        fingerprinted('FULL_REDUCED_SECTOR_BLOCK_ACTION_'+suffix, blocks, block_units)
        ends = ConstantEndPencil(reduction)
        physical('CONSTANT_END_FOURIER_MASS_'+suffix, ends.constant_mass)
        strong_symbol_units = ends.strong_matrix_units(pencil)
        weak_symbol_units = ends.weak_matrix_units(blocks)
        curl, gauge, gauge_residual = ends.curl_gauge_operands(pencil)
        curl_unit = dimensions.measure(next(x for x in curl if x != 0))
        fingerprinted('TRANSVERSE_CURL_ANSATZ_SYMBOL_'+suffix, curl,
                      {(i,): curl_unit for i in range(len(curl))})
        physical('TRANSVERSE_POTENTIAL_GAUGE_'+suffix, gauge)
        fingerprinted('TRANSVERSE_CURL_GAUGE_RESIDUAL_'+suffix, cas(gauge_residual),
                      {(i, j): curl_unit for i, v in enumerate(gauge_residual) for j in range(len(v))})
        modes = FullPencilModes(ends, curl, weak_symbol_units)
        lift = ends.field_lift(pencil)
        lift_field_units = [dimensions.known[f] for f in pencil.fields]
        lift_trial_units = [dimensions.known[sp.Function('s11cdReducedTrial'+label)]
                            for label in ('A0','A1','A2','Theta','E','Phi')]
        fingerprinted('CLOSED_PHYSICAL_FIELD_LIFT_'+suffix, lift,
                      {(6*i+j,):tuple(a-b for a,b in zip(lift_field_units[i],lift_trial_units[j]))
                       for i in range(5) for j in range(6)})
        for label, end in (('REFERENCE', None), ('LEFT', -sp.oo), ('RIGHT', sp.oo)):
            if end is not None:
                physical('PROFILE_END_TRANSLATION_LIMIT_OPERANDS_'+label+'_'+suffix,
                         ends.profile_limit_operands(end))
            constant_strong = ends.background(pencil.strong, end)
            constant_kernel = ends.background(pencil.kernel, end)
            constant_blocks = ends.background(blocks, end)
            fingerprinted('CONSTANT_CLOSED_OPERATOR_ACTION_'+label+'_'+suffix,
                          constant_strong, strong_units)
            fingerprinted('CONSTANT_CANONICAL_KERNEL_ACTION_'+label+'_'+suffix,
                          constant_kernel, kernel_units)
            fingerprinted('CONSTANT_FULL_SECTOR_BLOCK_ACTION_'+label+'_'+suffix,
                          constant_blocks, block_units)
            strong_symbol = ends.strong_matrix(pencil, constant_strong)
            full_symbol = ends.weak_matrix(constant_blocks)
            current_builder = UniformSlabCurrent(reduction, rows['energy_basis_variable'], ends,
                                                strong_symbol[3,:])
            current = current_builder.emit(str(case[0]), end, label+'_'+suffix)
            modes.set_closed_operands(strong_symbol, lift, current, current_builder)
            if options.dev_symbol_cache:
                options.dev_symbol_cache.mkdir(parents=True, exist_ok=True)
                with (options.dev_symbol_cache/(label+'_'+suffix+'.pickle')).open('wb') as stream:
                    pickle.dump((full_symbol, curl, weak_symbol_units, dimensions.known,
                                 strong_symbol, rows['energy_basis_variable']['value']), stream)
            fingerprinted('CLOSED_OPERATOR_SYMBOL_'+label+'_'+suffix, strong_symbol, strong_symbol_units)
            fingerprinted('FULL_SECTOR_PENCIL_SYMBOL_'+label+'_'+suffix, full_symbol, weak_symbol_units)
            for direction, rs, cs in (('TH', range(3), range(3, 6)), ('HT', range(3, 6), range(3))):
                coupling = full_symbol.extract(list(rs), list(cs))
                coupling_units = {(3*i+j,): weak_symbol_units[(6*ri+cj,)]
                                  for i, ri in enumerate(rs) for j, cj in enumerate(cs)}
                fingerprinted('K_'+label+'_'+direction+'_'+suffix, coupling, coupling_units)
            translated = ends.translate_kernel(constant_strong)
            translated_symbol = ends.strong_matrix(pencil, translated)
            fingerprinted('SIMULTANEOUS_END_TRANSLATION_OPERAND_'+label+'_'+suffix, translated, strong_units)
            fingerprinted('TRANSLATED_CLOSED_OPERATOR_SYMBOL_'+label+'_'+suffix,
                          translated_symbol, strong_symbol_units)
            fingerprinted('END_TRANSLATION_SYMBOL_RESIDUAL_'+label+'_'+suffix,
                          translated_symbol-strong_symbol, strong_symbol_units)
            unresolved = {
                'NORMAL_INTEGRALS': tuple(sorted((strong_symbol.atoms(sp.Integral)
                                                  | full_symbol.atoms(sp.Integral)), key=sp.default_sort_key)),
                'NORMAL_COORDINATES': tuple(sorted((dag_free_symbols(strong_symbol)
                    | dag_free_symbols(full_symbol)) & {reduction.z, reduction.zp, ends.shift}, key=sp.default_sort_key)),
                'ABEL_REGULATOR_OCCURRENCES': int(strong_symbol.has(reduction.regulator))
                                             + int(full_symbol.has(reduction.regulator)),
            }
            emit('ASYMPTOTIC_SYMBOL_REMAINDER_CENSUS_'+label+'_'+suffix, unresolved)
            if any(unresolved.values()):
                raise ValueError(('uncontracted asymptotic symbol', label, case))
            if label != 'REFERENCE':
                modes.solve(full_symbol, label+'_'+suffix, -1 if label == 'LEFT' else 1)
            if channel_input is not None and options.channel_input_scope=='modes-and-jets':
                modes.solve_input(full_symbol, label+'_'+suffix, -1 if label == 'LEFT' else 1,
                                  channel_input, reference=label == 'REFERENCE')
            spectrum = EndSpectrumCoverage(modes,strong_symbol,full_symbol,lift,strong_symbol_units)
            spectral_result = spectrum.construct(label+'_'+suffix,-1 if label=='LEFT' else 1,
                               reference=label=='REFERENCE',sample_index=0)
            continuation = EndResolventAudit(modes,strong_symbol,strong_symbol_units,
                reduced_branch_bindings[(CLOSED_KEYS[0],case)])
            continuation.construct(label+'_'+suffix,reference=label=='REFERENCE',sample_index=0,spectrum=spectral_result)
            if channel_input is not None:
                spectral_result = spectrum.construct(label+'_'+suffix,-1 if label=='LEFT' else 1,
                                   channel_input=channel_input,reference=label=='REFERENCE')
                continuation.construct(label+'_'+suffix,channel_input=channel_input,
                                       reference=label=='REFERENCE',spectrum=spectral_result)
                exceptional=EndExceptionalSlice(spectrum).construct(label+'_'+suffix,channel_input,reference=label=='REFERENCE')
                BulkExceptionalSlice(spectrum,reduced_branch_bindings[(CLOSED_KEYS[0],case)]).construct(
                    label+'_'+suffix,channel_input,reference=label=='REFERENCE',end_data=exceptional)
                ThresholdModeAudit(spectrum,reduced_branch_bindings[(CLOSED_KEYS[0],case)]).construct(
                    label+'_'+suffix,channel_input,reference=label=='REFERENCE',end_data=exceptional)
    emit('PENCIL_DIMENSION_CONSTRAINT_RESIDUALS', sorted(dimensions.constraints, key=sp.default_sort_key))
    if dimensions.constraints:
        raise ValueError('pencil dimensional analysis has surfaced unresolved constraints')
    emit('IMPLEMENTATION_CHECKPOINT', ('REDUCED_ACTION', 'DIMENSIONS', 'MULTIGRADES',
                                      'ABEL_WEAK_PRESCRIPTION', 'PENCIL_ACTION', 'OFF_DIAGONAL_EXTRACTION',
                                      'FULL_SECTOR_BLOCK_ACTION', 'REFERENCE_AND_END_SYMBOLS',
                                      'REFERENCE_AND_END_CANONICAL_COUPLINGS', 'SIMULTANEOUS_END_TRANSLATION',
                                      'CLASS_RESTRICTED_ACTION_RECONSTRUCTION', 'GAUGE_QUOTIENT_CHARTS',
                                      'ALL_CARRIER_INVERSE_FOURIER_ROUNDTRIPS',
                                      'POSITIVE_FREQUENCY_SPECTRAL_PIT', 'FIRST_GRADE_END_MODE_JETS',
                                      'NULLSPACE_CLASSIFIER_PROJECTORS', 'NONLINEAR_FREQUENCY_PAIRING',
                                      'REGULAR_RECTANGULAR_MODE_JETS',
                                      'PHYSICAL_FIELD_END_SPECTRUM_COVERAGE',
                                      'EXPLICIT_JOINT_BULK_SHEET_PATHS_AND_CUT_BANKS',
                                      'CONSTANT_END_RESOLVENTS_AND_NORMAL_POLE_RESIDUES',
                                      'EXCEPTIONAL_FREQUENCY_SLICE_AND_GENERIC_RANK_IDENTITIES',
                                      'INDEPENDENT_BULK_BRANCH_AND_DENOMINATOR_FREQUENCY_SLICES',
                                      'GENERALIZED_NORMAL_THRESHOLD_CHAINS',
                                      'S11B_CONSERVATIVE_SLAB_CURRENT', 'CLOSED_PHYSICAL_FIELD_LIFT'))
    emit('CHANNEL_INPUT_EXECUTION', channel_input is not None)
    emit('OUTSTANDING_CONSTRUCTIONS', ('FULL_END_SPECTRA_BEYOND_REFERENCE_MODE_JETS',
         'GENERIC_DOMAIN_SHEET_CONTINUATION', 'CLOSED_NONLOCAL_BULK_CURRENT_AND_FLUX_NORMALIZATION',
         'COMPLETE_TWO_ENDED_SCATTERING', 'POLES_RIESZ_OVERLAP', 'SURVIVAL',
         'FLUX_BOOKKEEPING', 'WEAK_COEFFICIENTS', 'SECTION_5_CONTROLS', 'OWN_ROWS_EXPORT'))
    emit('RESOURCE_MEASUREMENTS', (time.monotonic()-started,
                                  resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    emit('EMISSION_LINES', emission_index(EMISSION_LINES))
    emit('PROCESS_COMPLETION', sp.Integer(os.getpid()))


if __name__ == '__main__':
    run()
