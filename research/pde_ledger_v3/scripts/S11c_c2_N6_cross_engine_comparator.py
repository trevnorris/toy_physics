#!/usr/bin/env python3
"""N6 name/CAS-only cross-engine measurement instrument (frozen T7).

No measured builder is imported. See S11c_c2_SHARED_PHYSICS §§1b,5c,N8.
Schema sources: reconcile coefficient_table/source_table/closed_response and
emit_circuits; covariance source_columns/prolonged_phi; diagnostic
slot_guard/closure_guard/Compiler.scalar; WL carrierMap/sourceMap/weakFromCarrier
and the guard construction at lines 903–931. Representation identifications
(blocks, kernels, waves, components, applied fields) are deliberately absent.

Default per object (one family × anchoring × density): 900 seconds wallclock,
8 GiB RSS, enforced by a parent monitoring an isolated worker. Each rational
residual leaf has a 2 second algebra budget, returning inherited ResidualFailure
on timeout. A killed object is named in deferred accounting and the deferred
ledger, with its input addresses; previously emitted leaves remain visible.
"""
from __future__ import annotations

import argparse
import ast
import functools
import gc
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import re
import resource
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from typing import Iterable

import sympy as sp
from sympy.core.function import AppliedUndef

import S11c_c1_cross_engine_comparator as base
from S11b_cross_engine_comparator import (
    Association, BooleanNotResidualable, Mismatch, ResidualAssociation,
    ResidualFailure, TextAtom, UndecidedResidual, residual,
)

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PY = [ROOT / 'scripts/out' / f'S11c_c2_N6_{s}_sympy.out'
              for s in ('covariance', 'reconcile', 'diagnostic')]
DEFAULT_WL = ROOT / 'mathematica/out/S11c_c2_N6_mathematica_audit.out'
RC_NAMES = ('R_N6 SPLIT_SUM SPLIT_CHECK CARRIER_EULERIAN CARRIER_MATERIAL '
            'CARRIER_CHANNEL CARRIER_BRIDGE_RESIDUAL SOURCE_EULERIAN SOURCE_MATERIAL '
            'SOURCE_CHANNEL SOURCE_BRIDGE_RESIDUAL CROSS_CHANNEL EULERIAN_OPERAND '
            'MATERIAL_OPERAND DIMENSIONS FROZEN_RELATIONS ADVECTION_ABSENCE').split()
COV_NAMES = ('R_COV R_COV_BASELINE R_COV_INCREMENT R_COV_CONTROL_DELTA SOURCE_ACTUAL '
             'SOURCE_PREDICTED SOURCE_BASELINE SOURCE_CONTROL_DELTA FROZEN_PHI '
             'PHI_DOMAIN_CENSUS ACTUAL_CONTROL_PARAMETERS').split()
GUARDS = {f'N6_{kind}_GUARD_{part}' for kind in ('SLOT', 'CLOSURE')
          for part in ('NATIVE', 'CARRIER', 'RESIDUAL')}
SHARED = {f'N6RC_{n}' for n in RC_NAMES} | {f'N6COV_{n}' for n in COV_NAMES} | GUARDS
META = {'N6RC_DIMENSIONS', 'N6RC_FROZEN_RELATIONS', 'N6RC_ADVECTION_ABSENCE',
        'N6COV_FROZEN_PHI', 'N6COV_PHI_DOMAIN_CENSUS', 'N6COV_ACTUAL_CONTROL_PARAMETERS'}
CARRIERS = {f'N6RC_CARRIER_{n}' for n in ('EULERIAN', 'MATERIAL', 'BRIDGE_RESIDUAL')}
RC_SOURCES = {f'N6RC_SOURCE_{n}' for n in ('EULERIAN', 'MATERIAL', 'BRIDGE_RESIDUAL')}
COV_SOURCES = {f'N6COV_{n}' for n in COV_NAMES[:8]} - {'N6COV_R_COV_INCREMENT'}
AXES = ('ANCHORING', 'DENSITY', 'ROUTE', 'ROW', 'ROW_COMPONENT', 'PRESSURE_SLOT', 'PRESSURE_SLOT_2',
        'BLOCK_DIRECTION', 'TEST_FIELD', 'BLOCK_NAME', 'KERNEL_SIGNATURE', 'KERNEL_FAMILY',
        'FACE', 'WAVE', 'ETA', 'SIGMA', 'GUARD_KIND', 'COMPONENT', 'CENSUS_OBJECT', 'MAP_VARIABLE', 'DOMAIN_ATOM', 'FIELD_PATH')
ANCHORINGS = {'LAB_HELD', 'MATERIAL_ADVECTED'}
DENSITIES = {'RHO4_CONSTANT', 'RHOBR_CONSTANT'}
Key = tuple[tuple[str, str], ...]
InputError = base.InputError
AxisError = base.AxisError
split_top = base.split_top
wl_assoc_pairs = base.wl_assoc_pairs
wl_field = base.wl_field


def freeze(v):
    return tuple(map(freeze, v)) if isinstance(v, (list, tuple)) else v


def token(v):
    if isinstance(v, str): return v
    return json.dumps(v, separators=(',', ':'))


def make_key(items: Iterable[tuple[str, object]]) -> Key:
    d = {}
    for axis, value in items:
        if axis not in AXES or axis in d:
            raise AxisError(f'unknown or duplicate axis {axis}')
        if axis == 'ANCHORING' and value not in ANCHORINGS: raise AxisError(f'anchoring {value}')
        if axis == 'DENSITY' and value not in DENSITIES: raise AxisError(f'density {value}')
        if axis == 'ROUTE' and value not in {'EULERIAN', 'MATERIAL'}: raise AxisError(f'route {value}')
        if axis in ('FACE', 'ETA', 'SIGMA', 'ROW_COMPONENT') and isinstance(value, str) and re.fullmatch(r'-?\d+', value):
            value = int(value)
        if axis in ('ETA', 'SIGMA') and (type(value) is not int or value not in (0, 1)):
            raise AxisError(f'grade {axis}={value}')
        if axis == 'ROW_COMPONENT' and (type(value) is not int or value not in (1,2,3)):
            raise AxisError(f'Cartesian row component {value}')
        if axis == 'FACE':
            value = 0 if value == 'SUM' else value
            if type(value) is not int or value not in (-1, 0, 1): raise AxisError(f'face {value}')
        d[axis] = token(value)
    return tuple((a, d[a]) for a in AXES if a in d)


def graded(items, grade, engine):
    if not isinstance(grade, (tuple, list)): raise AxisError('grade is not a list')
    if engine == 'WL':
        if len(grade) != 3 or grade[0] != 1: raise AxisError(f'WL grade {grade}')
        grade = grade[1:]
    if len(grade) != 2: raise AxisError(f'grade {grade}')
    return items + [('ETA', grade[0]), ('SIGMA', grade[1])]


# Only inverse lowerCamel spellings found in actual PY symbols/column labels.
ACTIVE_NAMES: dict[str, str] = {}


def name_token(value, engine):
    return ACTIVE_NAMES.get(value, value) if engine == 'WL' and isinstance(value, str) else value


def row_axes(row, engine):
    """Decode the *strong row container*, never a weak COMPONENT_AXES basis.

    PY diagnostic.flatten zips U[0:3], THETA, E_W with U0,U1,U2,THETA,E_W.
    WL carrierMap labels positions 1:3 U1,U2,U3; guards emit positions 1:5.
    This is an array-index origin, not an equality between representations.
    """
    if engine == 'WL' and type(row) is int:
        if row in (4,5): return [('ROW', {4:'THETA',5:'E_W'}[row])]
        if row not in (1,2,3): raise AxisError(f'WL strong row {row}')
        return [('ROW', 'U'), ('ROW_COMPONENT', row)]
    match = re.fullmatch(r'U([0-3])', row) if isinstance(row,str) else None
    if match:
        component = int(match[1]) + (1 if engine == 'SymPy' else 0)
        if component not in (1,2,3): raise AxisError(f'{engine} strong row {row}')
        return [('ROW', 'U'), ('ROW_COMPONENT', component)]
    return [('ROW', row)]


def decode_py_key(family, column, case=(), probe=None):
    c = column
    if family in CARRIERS:
        row, pressure, face, grade = c
        items = [*row_axes(row, 'SymPy'), ('PRESSURE_SLOT', pressure), ('FACE', face)]
    elif family in RC_SOURCES:
        face, grade = c; items = [('FACE', face)]
    elif family in COV_SOURCES:
        face, wave, grade = c; items = [('FACE', face), ('WAVE', wave)]
    elif family in GUARDS:
        items = [('ROUTE', probe), *row_axes(c[0], 'SymPy')]
        if '_SLOT_' in family:
            kind = c[1]; items += [('GUARD_KIND', kind)]
            if kind == 'linear': grade = c[2]
            elif kind == 'cross':
                items += [('PRESSURE_SLOT', c[2]), ('PRESSURE_SLOT_2', c[3])]; grade = c[4]
            elif kind == 'denominator_dependency':
                return make_key([*case, *items, ('PRESSURE_SLOT', c[2])])
            else: raise AxisError(f'unknown guard kind {kind}')
        else:
            row, face, grade = c; items += [('FACE', face)]
    else:
        block, grade, signature, face = c
        if not isinstance(block, (list, tuple)) or len(block) != 2: raise AxisError('PY block pair')
        items = [('BLOCK_DIRECTION', block[0]), ('TEST_FIELD', block[1]),
                 ('KERNEL_SIGNATURE', signature), ('FACE', face)]
    return make_key([*case, *graded(items, grade, 'PY')])


def wl_literal(raw):
    raw = raw.strip()
    if raw.startswith('{'):
        return [wl_literal(x) for x in wl_list(raw)]
    if raw.startswith('"'): return json.loads(raw)
    if re.fullmatch(r'-?\d+', raw): return int(raw)
    if re.fullmatch(r'[A-Za-z$][A-Za-z0-9$]*', raw): return raw
    raise AxisError(f'not an axis literal: {raw[:120]}')


def decode_wl_key(family, column, case=(), component=None):
    c = wl_literal(column) if isinstance(column, str) else column
    if family in CARRIERS:
        row, face, pressure, grade = c
        items = [*row_axes(row, 'WL'), ('PRESSURE_SLOT', name_token(pressure, 'WL')), ('FACE', face)]
    elif family in RC_SOURCES | COV_SOURCES:
        face, wave, grade = c; items = [('FACE', face), ('WAVE', name_token(wave, 'WL'))]
    elif family in GUARDS:
        if '_SLOT_' in family:
            route, row, face, grade = c
            items = [('ROUTE', route), *row_axes(row, 'WL'), ('FACE', face)]
        else:
            route, row, face, kernel, grade = c
            items = [('ROUTE', route), *row_axes(row, 'WL'), ('FACE', face), ('KERNEL_FAMILY', kernel)]
    else:
        block, kernel, grade, face = c
        items = [('BLOCK_NAME', block), ('KERNEL_FAMILY', kernel), ('FACE', face)]
    if component is not None: items += [('COMPONENT', component)]
    return make_key([*case, *graded(items, grade, 'WL')])


# Fast structural scanner: quoted strings are atomic; delimiters inside them
# cannot alter nesting. Yields spans, not a materialized 400 MB association.
STRUCT = re.compile(r'"(?:\\.|[^"\\])*"|<\||\|>|[\[\]{}(),]|->')


def spans(raw, start, stop, delimiter=','):
    stack = []; last = start
    closing = {']': '[', '}': '{', ')': '(', '|>': '<|'}
    for m in STRUCT.finditer(raw, start, stop):
        t = m.group()
        if t.startswith('"'): continue
        if t in ('[', '{', '(', '<|'): stack.append(t)
        elif t in closing:
            if not stack or stack.pop() != closing[t]: raise InputError('unbalanced WL container')
        elif t == delimiter and not stack:
            yield last, m.start(); last = m.end()
    if stack: raise InputError('unclosed WL container')
    if raw[last:stop].strip(): yield last, stop


def wl_list(raw):
    raw = raw.strip()
    if not (raw.startswith('{') and raw.endswith('}')): raise InputError('WL list required')
    return [raw[a:b].strip() for a, b in spans(raw, 1, len(raw)-1)]


def wl_pairs(raw):
    raw = raw.strip()
    if not (raw.startswith('<|') and raw.endswith('|>')): raise InputError('WL association required')
    for a, b in spans(raw, 2, len(raw)-2):
        parts = list(spans(raw, a, b, '->'))
        if len(parts) != 2: raise InputError('WL association entry needs exactly one rule')
        yield raw[parts[0][0]:parts[0][1]].strip(), raw[parts[1][0]:parts[1][1]].strip()


def wl_dict(raw):
    d = {}
    for k, v in wl_pairs(raw):
        k = wl_literal(k)
        if not isinstance(k, str) or k in d: raise InputError(f'duplicate/nontext field {k}')
        d[k] = v
    return d


@dataclass(frozen=True)
class Address:
    path: str
    offset: int
    length: int
    line: int = 0

    def read(self):
        with open(self.path, 'rb') as f:
            f.seek(self.offset); return f.read(self.length).decode('utf-8').strip()


@dataclass
class PyIndex:
    records: dict = field(default_factory=lambda: defaultdict(list))
    names: set = field(default_factory=set)
    symbols: set = field(default_factory=set)
    local_literals: set = field(default_factory=set)


def metadata_symbols(data):
    """Collect actual Basic atoms from shared metadata's expression strings.

    PROVENANCE source text is not a shared metadata schema and is never parsed.
    This closes the W_bg metadata-vs-DAG vocabulary gap without guessing aliases.
    """
    names = set()
    if isinstance(data, dict): values = [v for k,v in data.items() if not sealed_field(k)]
    elif isinstance(data, list): values = data
    elif isinstance(data, str):
        try:
            value = base.parse_sympy_payload(data)
            if isinstance(value, sp.Basic): names.update(s.name for s in value.atoms(sp.Symbol))
        except (SyntaxError,TypeError,ValueError,NameError): pass
        return names
    else: return names
    for v in values: names.update(metadata_symbols(v))
    return names


# WL jet[field, Sort[indices]] (lines 84–89) and PY symmetric *_d1d2
# records are bare identifier encodings. Applied functions/derivatives keep
# their heads, all arguments, and assumptions. No field-to-bare collapse.
JET_BASES = {'theta':'theta', 'e_W':'eW', 'u_1':'u1', 'u_2':'u2', 'u_3':'u3',
             'zeta_c':'zetaC', 'W_bg':'WBg', 'mu_R_bg':'muRBg',
             'w1_profile':'w1Profile', 'm1_profile':'m1Profile'}
POINT_NAMES = {}


def verified_spelling_maps(index):
    names = base.checked_mechanical_symbol_map(index.symbols)
    aliases = {}
    def register(wl, py):
        if wl in names and names[wl] != py: raise InputError(f'non-injective CAS spelling {wl}')
        if wl in aliases and aliases[wl] != py: raise InputError(f'non-injective jet spelling {wl}')
        aliases[wl] = py
    for py in sorted(index.symbols):
        for prefix, head in JET_BASES.items():
            if py == prefix: register(head, py)
            elif py.startswith(prefix+'_'):
                suffix = py[len(prefix)+1:]
                if re.fullmatch(r'(d[123])+',suffix):
                    indices = re.findall(r'd([123])',suffix)
                    if indices != sorted(indices): continue
                    register(head+'Jet'+''.join(indices),py)
    for wl,py in aliases.items():
        projection=base.mechanical_lower_camel(py)
        if projection!=wl and names.get(projection)==py: del names[projection]
    names.update(aliases)
    if len(set(names.values()))!=len(names): raise InputError('non-injective active CAS spelling map')
    points = {}
    for py,point in index.local_literals:
        # The location is part of the decoded variable, not erased or summed.
        for wl,target in aliases.items():
            if target == py:
                token = wl+'At'+point
                if token in names: raise InputError(f'point/global spelling collision {token}')
                if token in points and points[token] != (py,point): raise InputError('point spelling collision')
                points[token] = (py,point)
    return names,points


def load_py_jsonl(paths):
    index = PyIndex()
    for path in paths:
        with open(path, 'rb') as f:
            line_no = 0
            while True:
                offset = f.tell(); raw = f.readline()
                if not raw: break
                line_no += 1
                if not raw.strip(): continue
                try: r = json.loads(raw)
                except (ValueError, UnicodeError) as e: raise InputError(f'{path}:{line_no}: {e}') from e
                name = r.get('object', '')
                if not name.startswith('S11CC2_') or 'data' not in r: raise InputError('not an N6 JSON row')
                name = name.removeprefix('S11CC2_')
                case = (r.get('anchoring'), r.get('density'))
                if name in SHARED or name.endswith(('_NODES', '_ARITHMETIC_DAG')):
                    make_key(zip(('ANCHORING', 'DENSITY'), case))
                identity = (name, *case, r.get('probe'))
                index.records[identity].append(Address(str(path), offset, len(raw), line_no))
                index.names.add(name)
                if name.endswith('ARITHMETIC_DAG'):
                    for node in r['data']['nodes']:
                        if node['op'] == 'variable':
                            literal = node['args'][0]['literal']
                            if len(literal) == 2 and literal[1] == 'formal_jet':
                                index.symbols.update(base.real_py_symbol_names([literal[0]]))
                            elif len(literal) == 2 and literal[1] in ('global', 'X', 'Y'):
                                index.symbols.add(literal[0])
                                if literal[1] in ('X', 'Y'): index.local_literals.add(tuple(literal))
                if name in META:
                    index.symbols.update(metadata_symbols(r['data']))
                # Pressure/wave labels are actual symbolic names.
                if name in CARRIERS | COV_SOURCES:
                    for col in r['data'].get('columns', []): index.symbols.add(col[1])
    return index


def load_wl(path):
    """Index single-line ` = ` tags; never parse/materialize all tag payloads."""
    out = defaultdict(list)
    with open(path, 'rb') as f:
        line_no = 0
        while True:
            offset = f.tell(); line = f.readline()
            if not line: break
            line_no += 1
            if not line.strip(): continue
            match = re.match(rb'WL_S11CC2_([A-Z0-9_]+) = ', line)
            if not match or not line[match.end():].strip().startswith(b'<|') or not line.rstrip().endswith(b'|>'):
                raise InputError(f'{path}:{line_no}: N6 single-line = grammar required')
            out[match[1].decode()].append(Address(str(path), offset+match.end(), len(line)-match.end(), line_no))
    if not out: raise InputError('empty WL stream')
    return out


def wl_case_addresses(address):
    raw = address.read()
    for a, b in spans(raw, 2, len(raw)-2):
        pair = list(spans(raw, a, b, '->'))
        if len(pair) != 2: raise InputError('WL case rule required')
        ka, kb = pair[0]; va, vb = pair[1]
        case = wl_literal(raw[ka:kb])
        if len(case) != 2: raise AxisError('WL case must be anchoring,density')
        key = make_key(zip(('ANCHORING', 'DENSITY'), case))
        # All emitted WL text is ASCII; validate before using byte offsets.
        if not raw.isascii(): raise InputError('non-ASCII WL stream cannot use byte-address index')
        yield key, Address(address.path, address.offset+va, vb-va, address.line)


# Applied heads and bare identifiers are shielded from inherited geometry
# folds. c1 still supplies its CAS reader, held predicates and BoundIntegral.
# Only Plus/Times/Power are activated: they encode emitted arithmetic circuits.
CAS_WORDS = {'Inactive','Plus','Times','Power','Rational','Integer','Real','Complex',
             'Integrate','Equal','Greater','Less','GreaterEqual','LessEqual','Unequal',
             'FourierTransform','Derivative','HoldForm','I','Pi','Infinity','True','False'}
IDENT = re.compile(r'"(?:\\.|[^"\\])*"|[A-Za-z$][A-Za-z0-9$]*')


def _parse_wl_atomic(raw):
    names = {}; reverse = {}
    def protect(m):
        t = m.group()
        if t.startswith('"') or t in CAS_WORDS: return t
        if t not in names:
            names[t] = 'N6Shield' + str(len(names)); reverse[names[t]] = t
        return names[t]
    protected = IDENT.sub(protect, raw)
    protected = re.sub(r'Inactive\[(Plus|Times|Power)\]\[', r'\1[', protected)
    value = base.parse_wl_value(protected, ())
    def restore(v):
        if isinstance(v, Association): return Association(tuple((k, restore(x)) for k,x in v.entries))
        if isinstance(v, (tuple, list, sp.Tuple)): return tuple(map(restore, v))
        if not isinstance(v, sp.Basic): return v
        def restored_head(name):
            if name in reverse: return reverse[name]
            if name.startswith('HeldInactive') and name[len('HeldInactive'):] in reverse:
                return 'HeldInactive'+reverse[name[len('HeldInactive'):]]
            return name
        v = v.replace(lambda x: isinstance(x, AppliedUndef) and restored_head(x.func.__name__) != x.func.__name__,
                      lambda x: sp.Function(restored_head(x.func.__name__))(*x.args))
        def restored_symbol(name):
            original = reverse[name]
            if original in POINT_NAMES:
                py,point = POINT_NAMES[original]
                return sp.Function('N6LocalPoint')(sp.Symbol(py),sp.Symbol(point))
            return sp.Symbol(name_token(original,'WL'))
        return v.xreplace({x: restored_symbol(x.name) for x in v.atoms(sp.Symbol) if x.name in reverse})
    return restore(value)


def _parse_wl_uncached(raw):
    # The emitted circuit heads are already split by the inherited structural
    # grammar. Reuse c1 for atomic CAS forms; reconstruct arithmetic directly to
    # avoid repeatedly parsing copied subcircuits in a large WL expression.
    match = re.match(r'Inactive\[(Plus|Times|Power)\]\[', raw)
    if match and raw.endswith(']'):
        args = [parse_wl_value(raw[a:b].strip()) for a,b in spans(raw,match.end(),len(raw)-1)]
        if match[1]=='Plus': return sp.Add(*args)
        if match[1]=='Times': return sp.Mul(*args)
        if len(args)!=2: raise InputError('Inactive Power arity')
        return sp.Pow(*args)
    return _parse_wl_atomic(raw)


@functools.lru_cache(maxsize=1024)
def _small_wl_value(raw, spellings, points):
    return _parse_wl_uncached(raw)


def parse_wl_value(raw, key=()):
    raw=raw.strip()
    if len(raw)<=8192:
        return _small_wl_value(raw,tuple(ACTIVE_NAMES.items()),tuple(POINT_NAMES.items()))
    return _parse_wl_uncached(raw)


def parse_py_expr(raw):
    if isinstance(raw, bool): return raw
    if isinstance(raw, (int, float)): return sp.sympify(raw)
    if raw is None: return TextAtom('null')
    if isinstance(raw, list): return tuple(map(parse_py_expr, raw))
    if isinstance(raw, dict): return Association(tuple((k, parse_py_expr(v)) for k,v in raw.items()))
    return base.parse_sympy_payload(raw)


def dag_literal(value):
    if value == ['algebraic_i']: return sp.I
    if len(value) != 2: raise InputError(f'unknown DAG variable {value}')
    name, location = value
    if location == 'formal_jet':
        return base.parse_sympy_payload(name)
    if location == 'global': return sp.pi if name == 'pi' else sp.Symbol(name)
    if location in ('X', 'Y'):
        return sp.Function('N6LocalPoint')(sp.Symbol(name), sp.Symbol(location))
    if name == 'formal_integral_guard':
        return sp.Function('N6FormalIntegralGuard')(sp.Integer(location))
    raise InputError(f'unknown DAG literal variant {value}')


class DAG:
    """Read shared, topologically ordered refs; never treat root copies as DAGs."""
    def __init__(self, data):
        if data is None: raise InputError('missing ARITHMETIC_DAG')
        self.nodes = data['nodes']; self.cache = {}
        for i, n in enumerate(self.nodes):
            if n['op'] not in {'number','variable','add','mul','pow'}: raise InputError(f'DAG op {n["op"]}')
            arity = 1 if n['op'] in ('number','variable') else 2
            if len(n['args']) != arity: raise InputError('DAG op arity')
            if n['op']=='number' and (set(n['args'][0]) != {'literal'} or type(n['args'][0]['literal']) is not int):
                raise InputError('DAG number is an integer literal, never a native boolean')
            if n['op']=='variable' and set(n['args'][0]) != {'literal'}: raise InputError('DAG variable literal')
            if n['op']=='pow' and (set(n['args'][1]) != {'literal'} or type(n['args'][1]['literal']) is not int):
                raise InputError('DAG exponent literal')
            for arg in n['args']:
                if set(arg) not in ({'ref'}, {'literal'}): raise InputError('DAG arg shape')
                if 'ref' in arg and (type(arg['ref']) is not int or not 0 <= arg['ref'] < i):
                    raise InputError(f'DAG nonpreceding ref {arg}')

    def get(self, root):
        if type(root) is not int or not 0 <= root < len(self.nodes): raise InputError('bad DAG root id')
        stack = [root]
        while stack:
            i = stack[-1]
            if i in self.cache: stack.pop(); continue
            node = self.nodes[i]
            missing = [a['ref'] for a in node['args'] if 'ref' in a and a['ref'] not in self.cache]
            if missing: stack.extend(missing); continue
            args = [self.cache[a['ref']] if 'ref' in a else a['literal'] for a in node['args']]
            op = node['op']
            if op == 'number': value = sp.Integer(args[0])
            elif op == 'variable': value = dag_literal(args[0])
            elif op == 'add': value = sp.Add(*args)
            elif op == 'mul': value = sp.Mul(*args)
            else: value = sp.Pow(*args)
            self.cache[i] = value; stack.pop()
        return self.cache[root]


NONZERO_WITNESSED = TextAtom('NONZERO_WITNESSED')
NO_NONZERO_FOUND = TextAtom('NO_NONZERO_FOUND')


def support(value):
    if type(value) is not bool: raise InputError('PY support observation must be boolean')
    return NONZERO_WITNESSED if value else NO_NONZERO_FOUND


def support_residual(a, b):
    if a == NO_NONZERO_FOUND or b == NO_NONZERO_FOUND:
        return UndecidedResidual(a, b, 'one-sided support: NO_NONZERO_FOUND is not a zero proof')
    return residual(a, b)


def wl_support(raw, count):
    """Decode emitted SparseArray CSR, witnessing each component independently.

    No residue is ever passed to residual(). Values are surfaced separately.
    """
    if not raw.startswith('SparseArray['): raise InputError('PIT numerators need SparseArray')
    p = [raw[a:b].strip() for a,b in spans(raw, len('SparseArray['), len(raw)-1)]
    if len(p) != 4: raise InputError('SparseArray arity')
    dims = wl_literal(p[1]); default = int(p[2])
    if len(dims) != 2 or dims[1] != count: raise InputError('PIT component dimensions')
    csr = wl_list(p[3]); index = wl_list(csr[1])
    columns = wl_literal(index[1]); values = wl_literal(csr[2])
    pointers = wl_literal(index[0])
    if len(columns) != len(values) or len(pointers) != dims[0]+1 or pointers[-1] != len(values):
        raise InputError('SparseArray CSR shape')
    seen = [0]*count; positives = [False]*count
    for col, v in zip(columns, values):
        j = col[0]-1
        if not 0 <= j < count: raise InputError('SparseArray column out of bounds')
        seen[j] += 1; positives[j] |= v != 0
    return [support(positives[j] or (default != 0 and seen[j] < dims[0])) for j in range(count)]


@dataclass
class Leaf:
    key: Key
    raw: object
    loader: object
    dimensions: object = None
    observation: object = None
    value: object = None
    error: str | None = None


def build_case(key, raw, loader, **kw): return Leaf(key, raw, loader, **kw)


def materialize(case):
    try: case.value = case.loader()
    except Exception as e: case.error = f'{type(e).__name__}: {e}'
    return case.error is None


def release_case(case):
    if case: case.value = None


@dataclass
class Accounting:
    join: int = 0
    sympy_only: int = 0
    wl_only: int = 0
    duplicate_key: int = 0
    parse_failed: int = 0
    axis_set_mismatch: int = 0
    unmatched_key: int = 0
    deferred_oversize: int = 0
    zero_extract_failures: int = 0
    extracted_sympy: int = 0
    extracted_wl: int = 0
    objects_completed: int = 0


def serialise(v):
    if isinstance(v, dict): return {str(k): serialise(x) for k,x in v.items()}
    if isinstance(v, list): return [serialise(x) for x in v]
    if isinstance(v, Address): return asdict(v)
    if isinstance(v, (str, int, float, bool)) or v is None: return v
    return base.serialise(v)


def emit(kind, **fields):
    print(kind + ' ' + json.dumps({k: serialise(v) for k,v in fields.items()}, sort_keys=True), flush=True)


def seal(engine, family, key, field_name, value, reason='pit_sealed'):
    emit('SURFACED', engine=engine, family=family, key=key, field=field_name,
         operand=value, seal=reason)


def sealed_field(name):
    n = name.lower()
    if any(x in n for x in ('sha256', 'fingerprint', 'digest')): return 'digest_sealed'
    if n in {'primes','probe_numerators','probe_denominators','numerator_denominator',
             'sample_index','sample_points','sampler_entry','sampler_support'}: return 'pit_sealed'
    return None


def extract_py_numeric(family, row, nodes, dag, case):
    d = row['data']; cols = d['columns']; dims = row['dimension']; obs = d['nonzero_modular_numerator']
    if len(dims) != len(cols) or len(obs) != len(cols): raise InputError('PY numeric alignment')
    seal('SymPy', family, case, 'numerator_denominator', d['numerator_denominator'])
    if family not in GUARDS:
        if nodes is None or dag is None: raise InputError('missing _NODES or ARITHMETIC_DAG')
        nd = nodes['data']
        if nd['columns'] != cols or len(nd['root_ids']) != len(cols) or len(nd['root_nodes']) != len(cols):
            raise InputError('PY _NODES column/root alignment')
        for root, copy in zip(nd['root_ids'], nd['root_nodes']):
            if root < 0 or root >= len(dag.nodes) or dag.nodes[root] != copy: raise InputError('DAG root copy mismatch')
    leaves = []
    for i, col in enumerate(cols):
        key = decode_py_key(family, col, case, row.get('probe'))
        if family in GUARDS:
            raw = {'symbolic': 'not emitted by SymPy guard', 'column': col}
            loader = lambda: UndecidedResidual(None, None, 'SymPy guard has no symbolic sibling')
        else:
            root = nd['root_ids'][i]; raw = {'root_id': root, 'column': col}
            loader = lambda root=root: dag.get(root)
        leaves.append(build_case(key, raw, loader, dimensions=parse_py_expr(dims[i]), observation=support(obs[i])))
    return leaves


def extract_wl_numeric(family, raw, case):
    leaves = []
    for col, entry in wl_pairs(raw):
        fields = wl_dict(entry)
        if 'ARITHMETIC' not in fields: raise InputError('nonempty numeric slot missing ARITHMETIC')
        arithmetic = wl_list(fields['ARITHMETIC']); dims = wl_list(fields['DIMENSIONS'])
        components = wl_literal(fields['COMPONENT_AXES'])
        if len(arithmetic) != len(dims) or len(components) != len(dims): raise InputError('WL component alignment')
        observations = wl_support(fields['PROBE_NUMERATORS'], len(components))
        # Numeric payloads have no free metadata exemption: extra fields surfaced.
        for field_name, val in fields.items():
            if field_name not in {'ARITHMETIC','DIMENSIONS','COMPONENT_AXES'}:
                seal('WL', family, [case, col], field_name, val,
                     sealed_field(field_name) or 'surfaced_not_joined')
        scalar_schema = family in CARRIERS | RC_SOURCES | COV_SOURCES
        for i, expr in enumerate(arithmetic):
            # Singleton carrier/source arrays are containers of scalar slots.
            # Weak blocks and guards retain COMPONENT even when it equals 1.
            component = None if scalar_schema and components == [1] else components[i]
            key = decode_wl_key(family, col, case, component)
            leaves.append(build_case(key, expr, lambda expr=expr: parse_wl_value(expr),
                                    dimensions=dims[i], observation=observations[i]))
    return leaves


def field_path(path): return json.dumps(path, separators=(',', ':'))


def field_name(name):
    # Field case/underscore spelling only, no aliases such as kappa_a=ADVECTION.
    return name.lower()


def extract_meta(engine, family, payload, case):
    """Reach every addressable metadata leaf; preserve unlike shapes and names.

    Only MAP/substitution_map is an emitter container adapter mandated by the
    frozen brief. Domain keys themselves remain mechanical names, not jet bridges.
    Lists are structural indexed coordinates; no sum, broadcast, or collapse.
    """
    leaves = []
    def add(path, raw, loader, extra=()):
        leaves.append(build_case(make_key([*case, *extra, ('FIELD_PATH', field_path(path))]), raw, loader))
    def walk_py(v, path=(), extra=()):
        if isinstance(v, dict):
            if not v: add(path, {}, lambda: Association(()), extra)
            for k, item in v.items():
                why = sealed_field(k)
                if why: seal('SymPy', family, case, field_path([*path,k]), item, why); continue
                if family == 'N6COV_FROZEN_PHI' and k == 'substitution_map':
                    for variable, expr in item:
                        walk_py(expr, (*path, 'map'), (*extra, ('MAP_VARIABLE', variable)))
                elif family == 'N6COV_PHI_DOMAIN_CENSUS' and k == 'coverage':
                    for variable, observation, jet_path in item:
                        domain = (*extra, ('DOMAIN_ATOM',variable))
                        walk_py(observation, (*path,'coverage'),domain)
                        walk_py(jet_path, (*path,'coverage_jet_path'),domain)
                elif family == 'N6RC_DIMENSIONS' and path in (('eulerian',), ('material',)):
                    row, pressure = ast.literal_eval(k)
                    walk_py(item, (), (('ROUTE', path[0].upper()), *row_axes(row,'SymPy'), ('PRESSURE_SLOT', pressure)))
                else: walk_py(item, (*path, field_name(k)), extra)
        elif isinstance(v, list):
            if not v: add(path, [], lambda: (), extra)
            for i,item in enumerate(v): walk_py(item, (*path, i), extra)
        else:
            # Expressions arrive as str; plain prose remains text, never eval-to-bool.
            def loader(v=v):
                if isinstance(v, str):
                    if 'map' in path: return parse_py_expr(v)
                    try: return parse_py_expr(v)
                    except (SyntaxError, TypeError, ValueError): return TextAtom(v)
                return parse_py_expr(v)
            add(path, v, loader, extra)
    def walk_wl(raw, path=(), extra=()):
        raw = raw.strip()
        if raw.startswith('<|'):
            pairs = list(wl_pairs(raw))
            if not pairs: add(path, raw, lambda: Association(()), extra)
            for k, v in pairs:
                label = wl_literal(k)
                if isinstance(label, str):
                    why = sealed_field(label)
                    if why: seal('WL', family, case, field_path([*path,label]), v, why); continue
                    if family=='N6COV_PHI_DOMAIN_CENSUS' and path==('coverage',):
                        walk_wl(v,path,(*extra,('DOMAIN_ATOM',name_token(label,'WL'))))
                    else: walk_wl(v, (*path, field_name(label)), extra)
                else:
                    # DIMENSIONS OBJECTS contains the same typed numeric slots.
                    if family == 'N6RC_DIMENSIONS' and len(path) == 2 and path[0] == 'objects':
                        namespace, obj = path[1].split(':', 1)
                        subfamily = ('N6_' if namespace == 'guard' else 'N6'+namespace.upper()+'_') + obj.upper()
                        typed = decode_wl_key(subfamily, label)
                        walk_wl(v, (), (*typed, ('CENSUS_OBJECT', subfamily)))
                    else: walk_wl(v, (*path, token(label)), extra)
        elif raw.startswith('{'):
            parts = wl_list(raw)
            if not parts: add(path, raw, lambda: (), extra)
            for i,v in enumerate(parts):
                pair = list(spans(v, 0, len(v), '->'))
                if len(pair) == 2:
                    variable = v[pair[0][0]:pair[0][1]].strip()
                    expr = v[pair[1][0]:pair[1][1]].strip()
                    label = name_token(wl_literal(variable), 'WL')
                    if family=='N6COV_FROZEN_PHI' and path==('map',):
                        walk_wl(expr,path,(*extra,('MAP_VARIABLE',label)))
                    else: walk_wl(expr, (*path, label), extra)
                else: walk_wl(v, (*path, i), extra)
        else: add(path, raw, lambda raw=raw: parse_wl_value(raw), extra)
    if engine == 'SymPy': walk_py(payload)
    else: walk_wl(payload)
    return leaves


def zero_extract_guard(a, b, a_nonempty, b_nonempty, accounting):
    if a_nonempty and b_nonempty and (not a or not b):
        accounting.zero_extract_failures += 1
        emit('COVERAGE', finding='zero_extract_failures', reason='nonempty shared container extracted no leaves',
             extracted_sympy=len(a), extracted_wl=len(b))


def dimensions_value(engine, value):
    if engine == 'SymPy': return value
    return Association(tuple((k, parse_wl_value(v)) for k,v in wl_dict(value).items()))


def dimension_residual(a, b, budget):
    # computed is a set of vectors, as is WL COEFFICIENT_SUPPORT. Preserve the
    # entire emitted records beside this mechanical vector-support extraction.
    av = a.as_dict().get('computed') if isinstance(a, Association) else a
    bv = b.as_dict().get('COEFFICIENT_SUPPORT') if isinstance(b, Association) else b
    return residual(av, bv, leaf_budget_seconds=budget)


def measured_residual(a, b, budget):
    if isinstance(a, UndecidedResidual) or isinstance(b, UndecidedResidual):
        return UndecidedResidual(a, b, 'symbolic operand unavailable in emitted schema')
    return residual(a, b, leaf_budget_seconds=budget)


def compare_family(family, a, b, *, budget=2., accounting=None, progress=None):
    acc = accounting or Accounting()
    acc.extracted_sympy += len(a); acc.extracted_wl += len(b)
    aa = defaultdict(list); bb = defaultdict(list)
    for leaf in a: aa[leaf.key].append(leaf)
    for leaf in b: bb[leaf.key].append(leaf)
    shapes_a = {frozenset(k for k,v in key) for key in aa}
    shapes_b = {frozenset(k for k,v in key) for key in bb}
    if progress: progress(acc)
    for key in sorted(aa.keys() | bb.keys()):
        lefts, rights = aa.get(key,[]), bb.get(key,[])
        if len(lefts) > 1 or len(rights) > 1:
            acc.duplicate_key += len(lefts)+len(rights)-1
            emit('PAIRING_TABLE', family=family, key=key, pairing='duplicate_key',
                 operand_A=[x.raw for x in lefts], operand_B=[x.raw for x in rights],
                 A_minus_B=UndecidedResidual(None,None,'duplicate typed key; no arbitrary pairing'))
            continue
        left = lefts[0] if lefts else None; right = rights[0] if rights else None
        pairing = 'matched' if left and right else 'sympy_only_column' if left else 'wl_only_slot'
        emit('PAIRING_TABLE', family=family, key=key, pairing=pairing)
        for leaf in (left, right):
            if leaf is not None: materialize(leaf)
        av = left.value if left and not left.error else left.raw if left else None
        bv = right.value if right and not right.error else right.raw if right else None
        errors = [x.error for x in (left,right) if x and x.error]
        if errors:
            acc.parse_failed += len(errors)
            difference = ResidualFailure(av,bv,'parse_failed: '+'; '.join(errors))
        elif left is None or right is None:
            difference = UndecidedResidual(av,bv,'unmatched typed key; no mechanical sibling')
        else: difference = measured_residual(av,bv,budget)
        emit('CASE', family=family, key=key, operand_A=av, operand_B=bv, A_minus_B=difference)
        # All operands/residuals are printed before accounting or coverage guards.
        if left and right: acc.join += 1
        else:
            acc.sympy_only += int(left is not None); acc.wl_only += int(right is not None)
            acc.unmatched_key += 1
            other = bb if left else aa
            shape = {k for k,v in key}
            mismatch = bool(other) and frozenset(shape) not in (shapes_b if left else shapes_a)
            acc.axis_set_mismatch += int(mismatch)
            emit('COVERAGE', family=family, key=key, finding='axis_set_mismatch' if mismatch else 'unmatched_key',
                 reason='no mechanical sibling for full typed key')
        if (left and left.observation) or (right and right.observation):
            ad = dimensions_value('SymPy', left.dimensions) if left else None
            bd = dimensions_value('WL', right.dimensions) if right else None
            sa = left.observation if left else None; sb = right.observation if right else None
            emit('STRUCTURAL', family=family, key=key, operand_A={'dimensions':ad,'support':sa},
                 operand_B={'dimensions':bd,'support':sb},
                 A_minus_B={'dimension_vectors':dimension_residual(ad,bd,budget),
                            'dimension_record':residual(ad,bd,leaf_budget_seconds=budget),
                            'addition_consistency':residual(
                                ad.as_dict().get('consistent') if isinstance(ad,Association) else None,
                                bd.as_dict().get('ADDITION_CONSISTENCY') if isinstance(bd,Association) else None,
                                leaf_budget_seconds=budget),
                            'support':support_residual(sa,sb)} if left and right else
                            UndecidedResidual(ad,bd,'structural comparison requires matched keys'), pit_sealed=True)
        release_case(left); release_case(right)
        if progress: progress(acc)
    return acc


def get_unique(index, identity):
    rows = index.records.get(identity, [])
    if len(rows) > 1: raise InputError(f'duplicate JSON identity {identity}')
    return json.loads(rows[0].read()) if rows else None


def object_work(family, case, pyindex, wladdresses, budget, progress=None):
    acc = Accounting(); a=[]; b=[]; nonempty_a=False; nonempty_b=False
    key = make_key(zip(('ANCHORING','DENSITY'), case))
    rows = [(ident, refs) for ident,refs in pyindex.records.items() if ident[:3] == (family,*case)]
    try:
        dag = None
        if family not in GUARDS | META and rows:
            namespace = family.split('_')[0]
            d = get_unique(pyindex, (namespace+'_ARITHMETIC_DAG', *case, None))
            dag = DAG(d['data'] if d else None)
        for ident, refs in rows:
            if len(refs) != 1:
                acc.duplicate_key += len(refs)-1
                emit('COVERAGE', family=family, key=ident, finding='duplicate_key', operands=refs); continue
            row = json.loads(refs[0].read()); nonempty_a |= bool(row['data'])
            if family in META: a.extend(extract_meta('SymPy',family,row['data'],key))
            else:
                nodes = get_unique(pyindex, (family+'_NODES',*case,None)) if family not in GUARDS else None
                a.extend(extract_py_numeric(family,row,nodes,dag,key))
    except Exception as e:
        acc.parse_failed += 1
        emit('COVERAGE', family=family,key=key,finding='parse_failed',engine='SymPy',reason=f'{type(e).__name__}: {e}',
             input_addresses=[ref for _,refs in rows for ref in refs])
    try:
        if len(wladdresses) > 1:
            acc.duplicate_key += len(wladdresses)-1
            emit('COVERAGE',family=family,key=key,finding='duplicate_key',engine='WL',input_addresses=wladdresses)
        elif wladdresses:
            raw=wladdresses[0].read(); nonempty_b=raw.strip() != '<||>'
            b = extract_meta('WL',family,raw,key) if family in META else extract_wl_numeric(family,raw,key)
    except Exception as e:
        acc.parse_failed += 1
        emit('COVERAGE', family=family,key=key,finding='parse_failed',engine='WL',reason=f'{type(e).__name__}: {e}',
             input_addresses=wladdresses)
    compare_family(family,a,b,budget=budget,accounting=acc,progress=progress)
    zero_extract_guard(a,b,nonempty_a,nonempty_b,acc)
    return acc


def _worker(conn, *args):
    try:
        def progress(acc): conn.send(('progress',asdict(acc)))
        result = object_work(*args, progress=progress)
        conn.send(('done',asdict(result),resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    except BaseException as e:
        conn.send(('error',f'{type(e).__name__}: {e}'))
    finally: conn.close()


def bounded_object(family, case, pyindex, addresses, args):
    ctx=mp.get_context('fork'); parent,child=ctx.Pipe(duplex=False)
    sys.stdout.flush()
    proc=ctx.Process(target=_worker,args=(child,family,case,pyindex,addresses,args.residual_leaf_seconds))
    proc.start(); child.close(); start=time.monotonic(); last=Accounting(); peak=0; done=None; reason=None
    while proc.is_alive() or parent.poll():
        if parent.poll(.05):
            try: message=parent.recv()
            except EOFError: break
            if message[0]=='progress': last=Accounting(**message[1])
            elif message[0]=='done':
                done=Accounting(**message[1]); peak=max(peak,message[2])
                if peak>args.object_rss_mib*1024 or time.monotonic()-start>args.object_seconds:
                    last=done; done=None
                    reason='RSS ceiling exceeded' if peak>args.object_rss_mib*1024 else 'wallclock ceiling exceeded'
                break
            else: reason=message[1]; break
        try:
            with open(f'/proc/{proc.pid}/status') as f:
                memory=re.search(r'VmRSS:\s+(\d+)',f.read())
            rss=int(memory[1]) if memory else 0; peak=max(peak,rss)
        except FileNotFoundError: rss=0
        if rss > args.object_rss_mib*1024: reason='RSS ceiling exceeded'; break
        if time.monotonic()-start > args.object_seconds: reason='wallclock ceiling exceeded'; break
    if proc.is_alive() and done is None: proc.terminate()
    proc.join(2)
    if proc.is_alive(): proc.kill(); proc.join()
    parent.close()
    if done is not None:
        done.objects_completed=1
        return done,peak
    if reason in ('RSS ceiling exceeded','wallclock ceiling exceeded') or proc.exitcode == -9 or (reason and reason.startswith('MemoryError')):
        reason=reason or 'worker killed by memory pressure'; last.deferred_oversize += 1
        elapsed=time.monotonic()-start
        record=dict(family=family,case=case,reason=reason,runtime_seconds=elapsed,peak_rss_kib=peak,
                    rss_ceiling_mib=args.object_rss_mib,wallclock_ceiling_seconds=args.object_seconds,
                    wl_addresses=addresses,py_addresses=[ref for ident,refs in pyindex.records.items()
                    if ident[:3]==(family,*case) for ref in refs])
        emit('COVERAGE',finding='deferred_oversize',operand_A=record['py_addresses'],operand_B=addresses,
             A_minus_B=UndecidedResidual(None,None,reason),**record)
        with open(args.deferred_ledger,'a') as f:
            f.write('\n- N6 comparator deferred object: '+json.dumps(serialise(record),sort_keys=True)+
                    '. No residual target; remaining leaves undecided.\n')
    else:
        last.parse_failed += 1
        emit('COVERAGE',family=family,case=case,finding='parse_failed',reason=reason or f'worker exit {proc.exitcode}')
    return last,peak


def excluded_reason(name):
    if name.startswith('N6_') and name not in GUARDS: return 'no WL sibling'
    if 'REP_INVARIANCE_' in name or 'CONTROL_INDEPENDENCE_' in name or name.endswith('PREMISES'): return 'no WL sibling'
    if name in ('N6RC_PROVENANCE','N6COV_PROVENANCE'): return 'provenance surfaced not joined; digests sealed'
    if name.endswith(('_NODES','_ARITHMETIC_DAG')): return 'symbolic sibling substrate'
    return 'no WL sibling'


def run(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--py',type=Path,nargs='+',default=DEFAULT_PY)
    ap.add_argument('--wl',type=Path,default=DEFAULT_WL)
    ap.add_argument('--family',action='append',choices=sorted(SHARED))
    ap.add_argument('--object-seconds',type=float,default=900.)
    ap.add_argument('--object-rss-mib',type=float,default=8192.)
    ap.add_argument('--residual-leaf-seconds',type=float,default=2.)
    ap.add_argument('--deferred-ledger',type=Path,default=ROOT/'DEFERRED_HEAVY_RUNS.md')
    args=ap.parse_args(argv)
    for n in ('object_seconds','object_rss_mib','residual_leaf_seconds'):
        if not math.isfinite(getattr(args,n)) or getattr(args,n)<=0: ap.error(f'{n} must be finite and positive')
    start=time.monotonic(); total={}; peak=0
    emit('MEASUREMENT_SCOPE',supplied_unfalsifiable='§1–2 and supplied substrate',
         surfaced_not_joined='PIT residues, primes, sample points; digests; SymPy-only families',
         residual_target='none',scope='N6 two-route control; §3d self-energy families absent',
         object_rss_ceiling_mib=args.object_rss_mib,object_wallclock_ceiling_seconds=args.object_seconds,
         residual_leaf_seconds=args.residual_leaf_seconds)
    try:
        pyindex=load_py_jsonl(args.py); wlindex=load_wl(args.wl)
        global ACTIVE_NAMES, POINT_NAMES
        ACTIVE_NAMES,POINT_NAMES=verified_spelling_maps(pyindex)
        emit('SPELLING_MAP',folds=ACTIVE_NAMES,point_literals=POINT_NAMES,injectivity='checked before activation',
             schema_folds=['grade {1,eta,sigma} to (eta,sigma)','SUM to face 0','typed axis order',
                           'MAP/substitution_map container', 'coverage domain atom indexing',
                           'strong-row array index origin (PY zero-based; WL one-based)'],applied_heads='arguments preserved')
        for engine,names in (('SymPy',pyindex.names),('WL',wlindex)):
            emit('LOCAL_INVENTORY',engine=engine,tags=sorted(n for n in names if '_LOCAL_' in n))
        for ident,refs in pyindex.records.items():
            name=ident[0]
            if name in SHARED: continue
            reason=excluded_reason(name)
            for ref in refs:
                emit('SURFACED',engine='SymPy',family=name,identity=ident,reason=reason,
                     operand=json.loads(ref.read()),pit_sealed=any(x in name for x in ('PRIMES','PIT_','SAMPLE_')),
                     digest_sealed='PROVENANCE' in name)
            if reason == 'no WL sibling':
                emit('ACCOUNTING_ROW',family=name,finding='sympy_only',reason=reason,rows=len(refs))
        for name,refs in wlindex.items():
            if name in SHARED: continue
            for ref in refs: emit('SURFACED',engine='WL',family=name,operand=ref.read(),
                                  reason='local or provenance; surfaced not joined',
                                  pit_sealed=name=='N6_LOCAL_PROBE',digest_sealed='PROVENANCE' in name)
        for family in sorted(args.family or SHARED):
            acc=Accounting(); wc=defaultdict(list)
            for address in wlindex.get(family,[]):
                for key,ref in wl_case_addresses(address):
                    kd=dict(key); wc[(kd['ANCHORING'],kd['DENSITY'])].append(ref)
            pc={ident[1:3] for ident in pyindex.records if ident[0]==family}
            cases=pc | wc.keys()
            if not cases:
                emit('COVERAGE',family=family,finding='unmatched_key',reason='declared shared family absent from both inputs')
                acc.unmatched_key+=1
            for case in sorted(cases):
                emit('OBJECT_BEGIN',family=family,case=case)
                found,rss=bounded_object(family,case,pyindex,wc.get(case,[]),args); peak=max(peak,rss)
                for k,v in asdict(found).items(): setattr(acc,k,getattr(acc,k)+v)
            total[family]=acc; emit('ACCOUNTING',family=family,**asdict(acc)); gc.collect()
    except (OSError,ValueError,KeyError,TypeError) as e:
        emit('OPERATIONAL_ERROR',reason=f'{type(e).__name__}: {e}')
        emit('RUN_ACCOUNTING',families=len(total),families_with_join=sum(a.join>0 for a in total.values()),
             families_with_unpaired=sum(a.unmatched_key>0 for a in total.values()),parse_failed=1,
             deferred_oversize=sum(a.deferred_oversize for a in total.values()),
             zero_extract_failures=sum(a.zero_extract_failures for a in total.values()),
             runtime_seconds=time.monotonic()-start)
        return 2
    sums=Counter()
    for a in total.values(): sums.update(asdict(a))
    emit('RUN_ACCOUNTING',families=len(total),families_with_join=sum(a.join>0 for a in total.values()),
         families_with_unpaired=sum(a.unmatched_key>0 for a in total.values()),
         parse_failed=sums['parse_failed'],deferred_oversize=sums['deferred_oversize'],
         zero_extract_failures=sums['zero_extract_failures'],runtime_seconds=time.monotonic()-start,
         peak_rss_kib=max(peak,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    all_deferred=bool(total) and all(a.deferred_oversize and not a.objects_completed for a in total.values())
    return 2 if sums['parse_failed'] or sums['zero_extract_failures'] or all_deferred else 0


if __name__=='__main__': raise SystemExit(run())
