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

import sympy as sp
import numpy as np
from sympy.core.function import AppliedUndef
from sympy.core.symbol import Str

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ledger_fold import load_model, check_consumer, assert_lookups_equal_manifest

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
)
EMISSION_LINES = {}
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
    print(tag + ': ' + sp.srepr(cas(value)), flush=True)


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


class EdgeReconstruction:
    """Reconstruct on the supplied tangentially homogeneous ansatz class.

    The source route eliminates delta constraints one variable at a time.
    It does not consume EdgeReduction's constraint matrix, roots or Jacobian.
    The inverse lift restores the outgoing character and original normal
    coordinate. This is a class-restricted reconstruction, not an inverse on
    arbitrary three-dimensional backgrounds.
    """

    def __init__(self, reduction):
        self.r = reduction

    def lift(self, value):
        r = self.r
        inverse = {v: k for k, v in r.normal_map.items()}
        return map_leaves(value, lambda e: e.xreplace(inverse)*r.wave_phase(r.x))

    @lru_cache(maxsize=None)
    def source_integral(self, original):
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
        data = {h: r.hat(h) for h in hats}
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

    def row(self, source, reduced, records, suffix, source_units):
        r = self.r
        replacements = {}
        for index, (original, image, _) in enumerate(records):
            independent, witness = self.source_integral(original)
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


def run():
    global PHYSICAL_METADATA
    os.chdir(ROOT)
    faulthandler.register(signal.SIGUSR1)
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=['ALL']+[a+'__'+r for a in
        ('LAB_HELD', 'MATERIAL_ADVECTED') for r in ('RHO4_CONSTANT', 'RHOBR_CONSTANT')], default='ALL')
    parser.add_argument('--dev-symbol-cache', type=Path)
    channel_options = parser.add_mutually_exclusive_group()
    channel_options.add_argument('--channel-input-json',
                                 help='explicit unit_frame, parameters and independent w/m profile input')
    channel_options.add_argument('--channel-input-file', type=Path,
                                 help='JSON file containing the same explicit channel input')
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
    input_text = (options.channel_input_file.read_text() if options.channel_input_file
                  else options.channel_input_json)
    channel_input = ChannelInput(reduction, json.loads(input_text)) if input_text else None
    if channel_input is not None:
        channel_input.emit_inputs()
    reconstruction = EdgeReconstruction(reduction)
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
    for row_key in CLOSED_KEYS:
        for case, payload in rows[row_key]['value']:
            if not selected(case):
                continue
            suffix = '_'.join((row_key, *(str(v) for v in case)))
            value = named(payload, 'VALUE')
            value_units = payload_units(value, named(payload, 'DIMENSION_L_T_M'))
            for index, hat in enumerate(sorted(payload.atoms(AppliedUndef), key=sp.default_sort_key)):
                if hat.func.__name__.startswith('s11cc2Fourier'):
                    hat_reduction = reduction.hat(hat)
                    hat_unit = dimensions.measure(hat)
                    image_unit = tuple(a+2*b for a,b in zip(hat_unit,dimensions.measure(reduction.tangents[0])))
                    physical('PROFILE_CARRIER_REDUCTION_' + suffix + '_' + str(index),
                             (hat, hat_reduction), operands=(hat, hat_reduction[0]),
                             zero_dimensions={(0,):hat_unit,(1,):image_unit})
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
            reconstruction.row(value, reduced, records, suffix, value_units)
            reduced_payload = reduction.payload(payload, reduced)
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
            if channel_input is not None:
                modes.solve_input(full_symbol, label+'_'+suffix, -1 if label == 'LEFT' else 1,
                                  channel_input, reference=label == 'REFERENCE')
    emit('PENCIL_DIMENSION_CONSTRAINT_RESIDUALS', sorted(dimensions.constraints, key=sp.default_sort_key))
    if dimensions.constraints:
        raise ValueError('pencil dimensional analysis has surfaced unresolved constraints')
    emit('IMPLEMENTATION_CHECKPOINT', ('REDUCED_ACTION', 'DIMENSIONS', 'MULTIGRADES',
                                      'ABEL_WEAK_PRESCRIPTION', 'PENCIL_ACTION', 'OFF_DIAGONAL_EXTRACTION',
                                      'FULL_SECTOR_BLOCK_ACTION', 'REFERENCE_AND_END_SYMBOLS',
                                      'REFERENCE_AND_END_CANONICAL_COUPLINGS', 'SIMULTANEOUS_END_TRANSLATION',
                                      'CLASS_RESTRICTED_ACTION_RECONSTRUCTION', 'GAUGE_QUOTIENT_CHARTS',
                                      'POSITIVE_FREQUENCY_SPECTRAL_PIT', 'FIRST_GRADE_END_MODE_JETS',
                                      'NULLSPACE_CLASSIFIER_PROJECTORS', 'NONLINEAR_FREQUENCY_PAIRING',
                                      'REGULAR_RECTANGULAR_MODE_JETS',
                                      'S11B_CONSERVATIVE_SLAB_CURRENT', 'CLOSED_PHYSICAL_FIELD_LIFT'))
    emit('CHANNEL_INPUT_EXECUTION', channel_input is not None)
    emit('OUTSTANDING_CONSTRUCTIONS', ('ALL_CARRIER_INVERSE_FOURIER_ROUNDTRIPS',
         'FULL_END_SPECTRA_BEYOND_REFERENCE_MODE_JETS',
         'GENERIC_DOMAIN_SHEET_CONTINUATION', 'CLOSED_NONLOCAL_BULK_CURRENT_AND_FLUX_NORMALIZATION',
         'COMPLETE_TWO_ENDED_SCATTERING', 'POLES_RIESZ_OVERLAP', 'SURVIVAL',
         'FLUX_BOOKKEEPING', 'WEAK_COEFFICIENTS', 'SECTION_5_CONTROLS', 'OWN_ROWS_EXPORT'))
    emit('RESOURCE_MEASUREMENTS', (time.monotonic()-started,
                                  resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    emit('EMISSION_LINES', EMISSION_LINES.copy())
    emit('PROCESS_COMPLETION', sp.Integer(os.getpid()))


if __name__ == '__main__':
    run()
