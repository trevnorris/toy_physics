#!/usr/bin/env python3
"""O2 object comparator. Read-only; no engine imports and no result criterion.

The IR is deliberately not a CAS: held heads, constructor options and binders
survive parsing. Only the explicitly supported scalar algebra is sent to SymPy.
An OPEN action is never made into a scalar indeterminate to obtain a residual.
All paths are zero based; a path addresses an emitted occurrence, not a name.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import json
import hashlib
from pathlib import Path
import re
import resource
import sys
import time

import sympy as sp

HERE = Path(__file__).resolve().parent
PY_SOURCE = 'research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py'
WL_SOURCE = 'research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl'
DEFAULT_PY = HERE / 'out/O2_live_balance_sympy_audit.out'
DEFAULT_WL = HERE.parent / 'mathematica/out/O2_live_balance_mathematica_audit.out'


class InputError(ValueError):
    """Missing, ambiguous, or ungrammatical input/configuration."""


@dataclass(frozen=True, slots=True)
class Node:
    head: str
    args: tuple = ()
    value: str = ''
    options: tuple = ()


def atom(head, value):
    return Node(head, value=str(value))


def seq(*items):
    return Node('Sequence', tuple(items))


def record(items):
    keys = [key for key, _ in items]
    if len(keys) != len(set(keys)):
        raise InputError('duplicate record key')
    return Node('Record', tuple(Node('Entry', (value,), key) for key, value in items))


def py_parse(source):
    """Parse lossless constructor grammar without eval or engine execution."""
    def walk(n):
        if isinstance(n, ast.Constant):
            return atom('Boolean' if isinstance(n.value, bool) else
                        'Text' if isinstance(n.value, str) else 'Number', n.value)
        if isinstance(n, ast.Name):
            return atom('Boolean' if n.id in ('true', 'false', 'True', 'False')
                        else 'Name', n.id)
        if isinstance(n, (ast.Tuple, ast.List)):
            return seq(*(walk(a) for a in n.elts))
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            v = walk(n.operand)
            if v.head == 'Number':
                return atom('Number', '-' + v.value)
            raise InputError('non-numeric Python unary expression')
        if not isinstance(n, ast.Call):
            raise InputError('outside lossless constructor grammar: ' + type(n).__name__)
        options = []
        for k in n.keywords:
            if k.arg is None:
                d = ast.literal_eval(k.value)
                if not isinstance(d, dict):
                    raise InputError('constructor options are not a dictionary')
                options.extend((str(key), repr(value)) for key, value in d.items())
            else:
                options.append((k.arg, repr(ast.literal_eval(k.value))))
        options = tuple(sorted(options))
        a = tuple(walk(x) for x in n.args)
        if isinstance(n.func, ast.Call):
            return Node('Apply', (walk(n.func), *a), options=options)
        if not isinstance(n.func, ast.Name):
            raise InputError('non-constructor Python call')
        h = n.func.id
        if h in ('Symbol', 'Dummy', 'Str', 'Function'):
            if len(a) != 1 or a[0].head != 'Text':
                raise InputError('invalid named constructor ' + h)
            return Node({'Symbol': 'Symbol', 'Dummy': 'Dummy', 'Str': 'Text',
                         'Function': 'FunctionName'}[h], value=a[0].value, options=options)
        if h == 'Integer':
            if len(a) != 1 or a[0].head != 'Number':
                raise InputError('invalid Integer')
            return a[0]
        if h == 'Tuple':
            if a and all(v.head == 'Sequence' and len(v.args) == 2
                         and v.args[0].head == 'Text' for v in a):
                return record([(v.args[0].value, v.args[1]) for v in a])
            return Node('Sequence', a)
        if h in ('ImmutableDenseMatrix', 'MutableDenseMatrix'):
            if len(a) != 1 or a[0].head != 'Sequence':
                raise InputError('invalid matrix constructor')
            return Node('Matrix', a[0].args)
        return Node(h, a, options=options)
    try:
        return walk(ast.parse(source, mode='eval').body)
    except (SyntaxError, TypeError, ValueError) as exc:
        raise InputError('SymPy grammar: ' + str(exc)) from exc


# A small InputForm parser: arithmetic, associations, lists and arbitrary held
# application heads. No Mathematica evaluator, automatic head rewriting or alarm.
TOKEN = re.compile(r'\s+|"(?:\\.|[^"\\])*"|<\||\|>|->|==|!=|>=|<=|&&|\|\||'
                   r'[A-Za-z$][A-Za-z0-9$`_]*|(?:\d+(?:\.\d*)?|\.\d+)(?:\*\^[+-]?\d+)?|'
                   r'[\[\]{}(),+*/^<>!\-\']')


class WLParser:
    def __init__(self, source):
        self.tokens = []
        end = 0
        for match in TOKEN.finditer(source):
            if match.start() != end:
                raise InputError('Wolfram grammar near ' + repr(source[end:end+60]))
            end = match.end()
            if not match.group().isspace():
                self.tokens.append(match.group())
        if source[end:].strip():
            raise InputError('Wolfram trailing grammar: ' + repr(source[end:end+60]))
        self.tokens.append('<END>')
        self.i = 0

    def peek(self):
        return self.tokens[self.i]

    def take(self, expected=None):
        v = self.peek()
        if expected is not None and v != expected:
            raise InputError(f'Wolfram expected {expected!r}, found {v!r}')
        self.i += 1
        return v

    def separated(self, closing):
        values = []
        if self.peek() != closing:
            while True:
                values.append(self.expr())
                if self.peek() != ',':
                    break
                self.take(',')
        self.take(closing)
        return tuple(values)

    def expr(self, minimum=0):
        tok = self.take()
        if tok == '(':
            left = self.expr()
            self.take(')')
        elif tok == '{':
            left = Node('Sequence', self.separated('}'))
        elif tok == '<|':
            entries = self.separated('|>')
            if any(e.head != 'Rule' or e.args[0].head != 'Text' for e in entries):
                raise InputError('association requires text-key rules')
            left = record([(e.args[0].value, e.args[1]) for e in entries])
        elif tok.startswith('"'):
            left = atom('Text', ast.literal_eval(tok))
        elif tok in ('+', '-'):
            v = self.expr(35)
            left = v if tok == '+' else Node('Mul', (atom('Number', -1), v))
        elif re.fullmatch(r'\d+(?:\.\d*)?|\.\d+', tok):
            left = atom('Number', tok)
        elif tok in ('True', 'False'):
            left = atom('Boolean', tok)
        elif re.fullmatch(r'[A-Za-z$][A-Za-z0-9$`_]*', tok):
            left = atom('Symbol', tok.removeprefix('O2`'))
        else:
            raise InputError('unexpected Wolfram token ' + repr(tok))
        precedence = {'->': 2, '||': 4, '&&': 5, '==': 10, '!=': 10,
                      '>': 10, '<': 10, '>=': 10, '<=': 10,
                      '+': 20, '-': 20, '*': 30, '/': 30, '^': 40}
        heads = {'->': 'Rule', '==': 'Equality', '!=': 'Unequality', '>': 'StrictGreaterThan',
                 '<': 'StrictLessThan', '>=': 'GreaterThan', '<=': 'LessThan',
                 '&&': 'And', '||': 'Or', '+': 'Add', '*': 'Mul', '^': 'Pow'}
        while True:
            op = self.peek()
            if op == '[' and 50 >= minimum:
                self.take('[')
                left = Node('Apply', (left, *self.separated(']')))
                continue
            if op == "'" and 50 >= minimum:
                self.take()
                left = Node('Prime', (left,))
                continue
            implicit = (op in ('(', '{', '<|') or op.startswith('"') or
                        re.match(r'^[A-Za-z$\d]', op) is not None)
            prec = precedence.get(op, 30 if implicit else -1)
            if prec < minimum:
                break
            if not implicit:
                self.take()
            else:
                op = '*'
            right = self.expr(prec if op in ('^', '->') else prec + 1)
            if op == '-':
                right = Node('Mul', (atom('Number', -1), right))
                op = '+'
            elif op == '/':
                right = Node('Pow', (right, atom('Number', -1)))
                op = '*'
            head = heads[op]
            # Flatten associative syntax, without evaluating operands.
            args = (*left.args, right) if head in ('Add', 'Mul') and left.head == head else (left, right)
            left = Node(head, args)
        return left


def wl_parse(source):
    p = WLParser(source)
    result = p.expr()
    p.take('<END>')
    return result


def read_stream(path, engine):
    """S11c-b's next-tag multiline framing, reimplemented without old residuals.

    Newlines are preserved (including inside quoted strings); every line and
    every local tag is accounted for, and each payload is parsed independently.
    """
    pattern = re.compile(r'^(PY_(?:LOCAL_)?O2_|WL_O2_)([A-Z0-9_]+):\s?(.*)$')
    result, chunks, current = {}, [], None
    def finish():
        if current is None:
            return
        if current in result:
            raise InputError('duplicate tag ' + current)
        payload = '\n'.join(chunks)
        result[current] = py_parse(payload) if engine == 'py' else wl_parse(payload)
    try:
        with Path(path).open() as stream:
            for line in stream:
                line = line.rstrip('\r\n')
                m = pattern.fullmatch(line)
                if m:
                    if not m[1].startswith(engine.upper() + '_'):
                        raise InputError('wrong engine prefix')
                    finish()
                    current = ('LOCAL/' if 'LOCAL_' in m[1] else '') + m[2]
                    chunks = [m[3]]
                elif engine == 'wl' and current is not None:
                    chunks.append(line)
                elif line.strip():
                    raise InputError('unframed input: ' + line[:80])
        finish()
    except OSError as exc:
        raise InputError(str(exc)) from exc
    if not result:
        raise InputError('empty stream ' + str(path))
    return result


@dataclass(frozen=True)
class NameRow:
    py: str
    wl: str
    spec: str
    py_line: int
    wl_line: int
    kind: str = 'operand'


# Binding table, not an assertion that the objects named have equal values.
NAME_TABLE = tuple(NameRow(*r) for r in [
    ('x1','x1','§1 Cartesian coordinate 1',145,27,'coordinate'),
    ('x2','x2','§1 Cartesian coordinate 2',145,27,'coordinate'),
    ('x3','x3','§1 Cartesian coordinate 3',145,27,'coordinate'),
    ('t','t','§1 lab time',146,63,'coordinate'),
    ('w','w','§1 bulk coordinate',146,114,'coordinate'),
    ('ell','ell','§3.1 fixed reduction scale',147,33,'parameter'),
    ('c0','c0','§3.1 asymptotic speed',147,35,'parameter'),
    ('GM','GM','§7 independent orbital parameter',148,319,'parameter'),
    ('V_r','VR','§1 radial material velocity',152,29,'profile'),
    ('o2_rho_br_live','RhoBr','§1 live brane density',152,29,'profile'),
    ('mu_perp','MuPerp','§3.1 optical stiffness',152,30,'profile'),
    ('xi_w','XiW','§3.1 graph displacement',152,30,'profile'),
    ('h','H','§3.1 reduced embedding field',152,30,'profile'),
    ('delta','DeltaOpt','§3.1 optical speed shift',152,31,'profile'),
    ('j_n','Jn','§5 outward mass loss',152,31,'profile'),
    ('f','Fbulk','§3.1 bulk density contrast',152,31,'profile'),
    ('o2_c_gamma_live','CGamma','§3.1 live optical speed',309,34,'profile'),
    ('B_A13','BA13','§3.2 branch',155,45),
    ('I_br_live','IBrLive','§3.2 inertia',155,46),
    ('P_br_cons','PBrCons','§3.2 conservative momentum',155,47),
    ('T_br_live','TBrLive','§3.2 live stress',155,48),
    ('T_br_cons','TBrCons','§3.2 conservative stress',155,49),
    ('N_br_live','NBrLive','§3.2 normal response',155,50),
    ('A_rot_live','ARotLive','§3.2 rotational response',156,51),
    ('R_ref_strain_live','RRefStrainLive','§3.2 reference evolution',156,41),
    ('M_perp','MPerp','§3.2 O1 stiffness response',156,54),
    ('R_br','RBr','§3.2 O7 density response',156,55),
    ('E_h_live','EhLive','§3.2 O4 embedding relation',156,56),
    ('J_map','JMap','§3.3 O6 map',156,52),
    ('Pi_n','PiN','§3.3 O3 exchange momentum',157,70),
    ('H_core','HCore','§3.3 O5 holder',157,53),
    ('E_br_live','EBrLive','§6 energy density operand',157,72),
    ('J_E_live','JELive','§6 energy current operand',157,73),
    ('P_ref_relax_live','PRefRelaxLive','§6 relaxation power',157,75),
    ('P_convert_exchange_live','PConvertExchangeLive','§6 conversion power',158,76),
    ('P_boundary_live','PBoundaryLive','§6 boundary power',158,77),
    ('S_E_net','SENet','§6 supplier',158,78),
    ('P_E_supply','PESupply','§6 supply budget',158,79),
    ('C_ref','CRef','§8.6 energy reference',159,80),
    ('S12_local_source_controller','S12LocalConversionReturnControllers','§3.3 local source inventory',159,43),
    ('S12_boundary_domain','S12MouthCollarReturnIRBulkBoundary','§3.3 boundary inventory',159,44),
    ('T_hold_s','THold','§3.3 complete native support operand',231,66),
    ('o2_s_face','s','§5 native face label',217,61,'binder'),
    ('material_action_compatibility','UnresolvedStressInertiaNormalIdentifications',
     '§3.2,4 unresolved identification among material stress, inertia and normal response',162,133),
    ('energy_accounting_overlap','UnresolvedEnergyOccurrenceIdentifications',
     '§6 unresolved overlaps among material, boundary, conversion and supply energy occurrences',162,289),
    ('S12_reaction_system','OPENReactionSystem','§3.3,5 S12 additional momentum reaction system',160,209),
    ('face_support_partition','UnresolvedSupportPartition','§3.3 unresolved partition of complete face/support loading',163,177),
])


@dataclass(frozen=True)
class JoinRow:
    label: str
    py: tuple
    wl: tuple
    py_line: int
    wl_line: int
    spec: str
    layout: str = ''


def J(label, py, wl, pl, wl_line, spec, layout=''):
    return JoinRow(label, tuple(py.split('/')), tuple(wl.split('/')), pl, wl_line, spec, layout)


# Each selection is an actual emitted occurrence. No expression is manufactured
# from another emitted object to fill a missing counterpart. Matrix transpose is
# only the declared tangent storage convention (PY ambient rows; WL tangent rows).
JOIN_TABLE = [
    J('coordinates','BASIS/coordinates','BASIS_MEASURES_GEOMETRY/Coordinates',340,114,'§1'),
    J('basis','BASIS/ambient_basis','BASIS_MEASURES_GEOMETRY/Basis',341,113,'§1'),
    J('tangents','BASIS/graph_tangent','BASIS_MEASURES_GEOMETRY/Geometry/tangent_vectors',122,109,'§1,3.1','transpose'),
    J('graph_normal','BASIS/graph_normal','BASIS_MEASURES_GEOMETRY/Geometry/graph_normal',127,110,'§1,3.1'),
    J('density_measure','MEASURES/density_measure','BASIS_MEASURES_GEOMETRY/DensityMeasure',344,115,'§1'),
    J('native_measure','MEASURES/native_face_area','BASIS_MEASURES_GEOMETRY/NativeFaceMeasure',225,115,'§5'),
    J('native_map','MEASURES/map','BASIS_MEASURES_GEOMETRY/NativeFaceMap',346,116,'§5'),
    J('native_context','MEASURES/face_context','BASIS_MEASURES_GEOMETRY/NativeDependence',220,116,'§5'),
    J('native_face_set','MEASURES/face_set','MECHANICAL_LOAD/NativeFaceSet',219,179,'§5'),
    J('metric','GEOMETRY/metric','BASIS_MEASURES_GEOMETRY/Geometry/g_ij',123,97,'§3.1'),
    J('inverse_metric','GEOMETRY/inverse','BASIS_MEASURES_GEOMETRY/Geometry/g_inverse',124,98,'§3.1'),
    J('metric_determinant','GEOMETRY/determinant','BASIS_MEASURES_GEOMETRY/Geometry/det_g',125,99,'§3.1'),
    J('field_identity','GEOMETRY/identity','BASIS_MEASURES_GEOMETRY/FieldIdentity',354,117,'§3.1'),
    J('native_tangents','GEOMETRY/native_tangent','BASIS_MEASURES_GEOMETRY/Geometry/native_face_tangents',223,111,'§5'),
    J('native_area','GEOMETRY/native_area','BASIS_MEASURES_GEOMETRY/Geometry/native_face_measure',225,112,'§5'),
    J('native_normal','GEOMETRY/native_normal','BASIS_MEASURES_GEOMETRY/Geometry/native_face_normal',226,110,'§5'),
    J('material_velocity','MATERIAL_VELOCITY','BASIS_MEASURES_GEOMETRY/MaterialVelocity',129,118,'§1'),
    J('momentum_density','MOMENTUM_DENSITY','MATERIAL_MOMENTUM/DensityAction',194,146,'§4'),
    J('momentum_current','MOMENTUM_FLUX','MATERIAL_MOMENTUM/CurrentAction',195,149,'§4'),
    J('momentum_storage','MOMENTUM_STORAGE','MATERIAL_MOMENTUM/Storage',197,152,'§4'),
    J('momentum_transport','MOMENTUM_TRANSPORT','MATERIAL_MOMENTUM/SpatialTransport',198,153,'§4'),
    J('material_compatibility','MATERIAL_COMPATIBILITY','MATERIAL_MOMENTUM/NormalAccounting',212,162,'§3.2,4'),
    J('internal_force','INTERNAL_FORCE','INTERNAL_MATERIAL_FORCE/Components',211,168,'§4'),
    J('bulk_load','NATIVE_BULK_LOAD','MECHANICAL_LOAD/NativeBulkTraction',230,107,'§3.3,5'),
    J('native_hold','NATIVE_HOLD_LOAD','MECHANICAL_LOAD/NativeCompleteLoad',234,175,'§3.3,5'),
    J('mechanical_load','MECHANICAL_LOAD','MECHANICAL_LOAD/CoordinateComponents',247,194,'§4,5'),
    J('carried_momentum','CARRIED_MOMENTUM','EXCHANGE_MOMENTUM/Carried',253,206,'§5'),
    J('source_partners','SOURCE_PARTNERS','EXCHANGE_MOMENTUM/AdditionalPartners',254,207,'§5'),
    J('drive_body_entries','DRIVE/body_entries','DRIVE_PROVENANCE/SeparateBodyForceEntries',257,215,'§2'),
    J('drive_local_inventory','DRIVE/local_source','DRIVE_PROVENANCE/LocalSourceInventory',167,217,'§3.3,5'),
    J('drive_boundary_inventory','DRIVE/boundary','DRIVE_PROVENANCE/BoundaryInventory',168,217,'§3.3,5'),
    J('drive_occurrences','DRIVE/occurrences','DRIVE_PROVENANCE/EntryDependencies',373,218,'§2,4,5'),
    J('drive_representation','DRIVE/representation','DRIVE_PROVENANCE/Drive',374,216,'§2'),
    J('hold_inplane','HOLD_INPLANE','B_HOLD_LIVE/InPlane',261,237,'§4,5'),
    J('hold_bulk','HOLD_W','B_HOLD_LIVE/BulkCoordinate',261,237,'§4,5'),
    J('hold_normal','HOLD_GRAPH_NORMAL','B_HOLD_LIVE/GraphNormalProjection',376,238,'§4,5'),
    J('mass_current','MASS_INPUT/current','MASS_INPUT/Flux',264,242,'§3.1'),
    J('mass_divergence','MASS_INPUT/divergence','MASS_INPUT/Divergence',265,243,'§3.1'),
    J('mass_equation','MASS_INPUT/supplied_equation','MASS_INPUT/Equation',266,247,'§3.1'),
    J('mass_measure','MASS_INPUT/measure','MASS_INPUT/Measure',379,248,'§1'),
    J('mass_density','MASS_INPUT/density_input/density/0','MASS_INPUT/DensityOperand',172,245,'§3.1'),
    J('mass_residual','MASS_RESIDUAL','MASS_INPUT/Residual',267,247,'§3.1'),
    J('power_graph_velocity','MATERIAL_POWER_PAIRING/velocity','FORCE_POWER_PAIRINGS/GraphVelocity',381,301,'§6'),
    J('material_work','MATERIAL_POWER_PAIRING/work','FORCE_POWER_PAIRINGS/MaterialWork',282,263,'§6'),
    J('face_velocity','FACE_POWER_PAIRING/velocity','FORCE_POWER_PAIRINGS/NativeApplicationVelocity',270,256,'§6'),
    J('face_native_work','FACE_POWER_PAIRING/native_power','FORCE_POWER_PAIRINGS/NativeFaceWork',273,261,'§6'),
    J('face_coordinate_work','FACE_POWER_PAIRING/reduced_power','FORCE_POWER_PAIRINGS/CoordinateFaceWork',275,262,'§6'),
    J('face_shared_map','FACE_POWER_PAIRING/map','FORCE_POWER_PAIRINGS/SharedMap',389,304,'§5,6'),
    J('energy_storage','ENERGY_STORAGE','B_E_STEADY/Storage',294,280,'§6'),
    J('energy_transport','ENERGY_TRANSPORT','B_E_STEADY/Transport',295,281,'§6'),
    J('energy_power','ENERGY_POWER','B_E_STEADY/PowerOccurrence',298,288,'§6'),
    J('energy_balance','ENERGY_STEADY','B_E_STEADY/Object',307,300,'§6'),
    J('optical_inputs','COUPLED_INPUTS/optical','COUPLED_INPUTS_MODEL_POINT/OpticalIdentifications',310,34,'§3.1'),
    J('coupled_embedding','COUPLED_INPUTS/embedding','B_HOLD_LIVE/O4Identity',395,239,'§3.2'),
    J('missing_grades','COUPLED_INPUTS/unknown_grades_and_derivative_scales','COUPLED_INPUTS_MODEL_POINT/MissingGradesScales',312,324,'§7'),
    J('supplier','COUPLED_INPUTS/supplier_obligation/0','B_E_STEADY/Supplier',402,309,'§6'),
    J('supply_budget','COUPLED_INPUTS/supplier_obligation/1','B_E_STEADY/Budget',402,309,'§6'),
    J('interfaces','COUPLED_INPUTS/ownership','COUPLED_INPUTS_MODEL_POINT/LaterInterfaces',403,342,'§10'),
    J('restrictions','MODEL_POINT/restrictions','COUPLED_INPUTS_MODEL_POINT/Restrictions',316,326,'§7,8'),
    J('premises','MODEL_POINT/premises','COUPLED_INPUTS_MODEL_POINT/Premises',409,316,'§2'),
    J('epsilon','MODEL_POINT/epsilon','COUPLED_INPUTS_MODEL_POINT/Counting/epsilon',150,319,'§7'),
    J('optical_box','MODEL_POINT/optical_monomial_indices','COUPLED_INPUTS_MODEL_POINT/Counting/OpticalMonomialBox',414,320,'§7'),
    J('recorded_grades','MODEL_POINT/recorded_grades','COUPLED_INPUTS_MODEL_POINT/Counting/Grades',416,321,'§7'),
    J('local_name_register','LOCAL_NAMES','LOCAL_NAMES',586,344,'serialization inventory'),
]
for i, key in enumerate(('V_r','rho_br','mu_perp','xi_w','h','delta','j_n','f')):
    JOIN_TABLE.append(J('profile_' + key, f'PROFILES/{i}',
                        'COUPLED_INPUTS_MODEL_POINT/Profiles/' + key,152,
                        (29,29,30,30,30,31,31,31)[i],'§1,3.1'))
# Six explicitly paired coupled-input occurrences, not a join by name.
for i, j in ((0,0),(8,1),(10,2),(13,3),(11,4),(9,5)):
    JOIN_TABLE.append(J(f'coupled_operand_{i}',f'COUPLED_INPUTS/operands/{i}',
                        f'COUPLED_INPUTS_MODEL_POINT/CoupledInputs/{j}',155,315,'§3.2,3.3'))
# The repeated section/entry emissions provide independent occurrences for the
# register and derivative objects. @n selects the nth constructor argument,
# including @0 for an applied head; it never evaluates that constructor.
JOIN_TABLE.extend([
    J('lab_time','BASIS/time','MATERIAL_MOMENTUM/DifferentiatedSection/OtherDependence/@3',146,131,'§1'),
    J('outward_mass_loss','MASS_INPUT/outward_loss','EXCHANGE_MOMENTUM/MaterialIdentification/@1',378,212,'§3.1,5 outward material loss'),
    J('differential_measure','MATERIAL_INPUT_DIFFERENTIALS/measure','FORCE_POWER_PAIRINGS/Measure',207,304,'§1'),
    J('reference_operand','COUPLED_INPUTS/operands/7','MATERIAL_MOMENTUM/DifferentiatedSection/MaterialReference',156,144,'§3.2'),
    J('exchange_operand','COUPLED_INPUTS/operands/12','EXCHANGE_MOMENTUM/MaterialIdentification/@2',157,212,'§5'),
    J('energy_reference','COUPLED_INPUTS/operands/21','B_E_STEADY/EnergyReference',159,309,'§8.6'),
    J('nonpassive_obligation','COUPLED_INPUTS/operands/28','B_E_STEADY/NonPassiveObligation',161,310,'§6'),
    J('material_identifications','COUPLED_INPUTS/operands/29','MATERIAL_MOMENTUM/DifferentiatedSection/OtherDependence/@6',162,133,'§3.2'),
    J('energy_overlap','COUPLED_INPUTS/operands/30','B_E_STEADY/Entries/2/Object/@2/1',162,289,'§6'),
    J('open_calculus','COUPLED_INPUTS/open_action_semantics','MATERIAL_MOMENTUM/Calculus',400,163,'§3.2,4'),
    J('density_response_in_mass','MASS_INPUT/density_input/density/1','MATERIAL_MOMENTUM/DifferentiatedSection/OtherDependence/@1/8',172,128,'§3.2'),
    J('stiffness_response_in_mass','MASS_INPUT/density_input/optical_stiffness/1','MATERIAL_MOMENTUM/DifferentiatedSection/OtherDependence/@1/9',173,128,'§3.2'),
    J('stiffness_profile_in_mass','MASS_INPUT/density_input/optical_stiffness/0','MATERIAL_MOMENTUM/DifferentiatedSection/Profiles/mu_perp',173,130,'§1,3.1'),
    J('paired_force','MATERIAL_POWER_PAIRING/force','B_E_STEADY/Entries/2/Object/@3/MaterialWork/@3/ForceAction',381,265,'§6'),
    J('generalized_rates','MATERIAL_POWER_PAIRING/unfixed_normal_generalized_rates','B_E_STEADY/Entries/2/Object/@3/MaterialWork/@3/GeneralizedRates',277,266,'§6'),
])
for label,py_path,wl_path,pl,wl_line in (
    ('density_profile','density/0','Profiles/rho_br',172,130),
    ('stiffness_profile','optical_stiffness/0','Profiles/mu_perp',173,130),
    ('density_response','density/1','OtherDependence/@1/8',172,128),
    ('stiffness_response','optical_stiffness/1','OtherDependence/@1/9',173,128)):
    JOIN_TABLE.append(J('constitutive_'+label,'COUPLED_INPUTS/constitutive_inputs/'+py_path,
                        'B_E_STEADY/Entries/0/Object/@1/@3/'+wl_path,pl,wl_line,'§3.2,6'))
_WORK_REDUCTION = 'B_E_STEADY/Entries/2/Object/@3/MechanicalBoundaryWork/@1/@1/@2/@3'
for i in range(4):
    JOIN_TABLE.append(J(f'paired_native_traction_{i}',f'FACE_POWER_PAIRING/traction/{i}/0',
                        f'{_WORK_REDUCTION}/NativeIntegrand/@2/@{i}/@1',234,261,'§5,6'))
JOIN_TABLE.append(J('paired_native_measure','FACE_POWER_PAIRING/native_area',
                    _WORK_REDUCTION+'/NativeMeasure',225,184,'§5,6'))
for py_index,wl_index in ((1,0),(2,1),(3,2),(4,3),(5,4),(6,5),(22,10),(23,11)):
    JOIN_TABLE.append(J(f'register_material_{py_index}',f'COUPLED_INPUTS/operands/{py_index}',
                        f'MATERIAL_MOMENTUM/DifferentiatedSection/OtherDependence/@1/{wl_index}',155,127,'§3.2,3.3'))
for py_index,wl_index in ((14,0),(15,14),(16,1),(17,2),(18,3),(19,4),(20,5)):
    JOIN_TABLE.append(J(f'register_energy_{py_index}',f'COUPLED_INPUTS/operands/{py_index}',
                        f'B_E_STEADY/Entries/2/Object/@2/0/{wl_index}',157,271,'§6'))
for i,key in enumerate(('V_r','rho_br','xi_w')):
    JOIN_TABLE.append(J(f'differentiated_profile_{i}',f'MATERIAL_INPUT_DIFFERENTIALS/profiles/{i}',
                        'MATERIAL_MOMENTUM/DifferentiatedSection/Profiles/'+key,184,145,'§1,4'))
    for j in range(3):
        JOIN_TABLE.append(J(f'profile_gradient_{i}_{j}',f'MATERIAL_INPUT_DIFFERENTIALS/gradients/{i}/{j}',
                            f'B_HOLD_LIVE/Entries/1/Object/0/@{j}/@2/Profiles/{key}',204,135,'§1,4'))
for i in range(4):
    for j in range(3):
        JOIN_TABLE.append(J(f'velocity_gradient_{i}_{j}',f'MATERIAL_INPUT_DIFFERENTIALS/velocity_gradient/{i}/{j}',
                            f'B_HOLD_LIVE/Entries/1/Object/0/@{j}/@2/Velocity/{i}',205,136,'§1,4'))

# Provenance has different grain: one trace per SymPy publication versus one
# Origin per Wolfram association, with further entry origins for its balances.
for py_tag,wl_tag in (
    ('BASIS','BASIS_MEASURES_GEOMETRY'),('MOMENTUM_DENSITY','MATERIAL_MOMENTUM'),
    ('INTERNAL_FORCE','INTERNAL_MATERIAL_FORCE'),('MECHANICAL_LOAD','MECHANICAL_LOAD'),
    ('CARRIED_MOMENTUM','EXCHANGE_MOMENTUM'),('DRIVE','DRIVE_PROVENANCE'),
    ('HOLD_INPLANE','B_HOLD_LIVE'),('MASS_INPUT','MASS_INPUT'),
    ('FACE_POWER_PAIRING','FORCE_POWER_PAIRINGS'),('ENERGY_STEADY','B_E_STEADY'),
    ('MODEL_POINT','COUPLED_INPUTS_MODEL_POINT')):
    JOIN_TABLE.append(J('trace_'+py_tag,'TRACE/'+py_tag,wl_tag+'/Origin',421,
                        {'BASIS_MEASURES_GEOMETRY':118,'MATERIAL_MOMENTUM':164,
                         'INTERNAL_MATERIAL_FORCE':171,'MECHANICAL_LOAD':200,
                         'EXCHANGE_MOMENTUM':213,'DRIVE_PROVENANCE':219,'B_HOLD_LIVE':240,
                         'MASS_INPUT':250,'FORCE_POWER_PAIRINGS':305,'B_E_STEADY':312,
                         'COUPLED_INPUTS_MODEL_POINT':343}[wl_tag],'§9 provenance'))
for tag,index in (('MOMENTUM_STORAGE',0),('MOMENTUM_TRANSPORT',1),('NATIVE_HOLD_LOAD',3),('SOURCE_PARTNERS',5)):
    JOIN_TABLE.append(J('entry_trace_'+tag,'TRACE/'+tag,f'B_HOLD_LIVE/Entries/{index}/Origin',421,228+index,'§9 provenance'))
for tag,index in (('ENERGY_STORAGE',0),('ENERGY_TRANSPORT',1),('ENERGY_POWER',2)):
    JOIN_TABLE.append(J('entry_trace_'+tag,'TRACE/'+tag,f'B_E_STEADY/Entries/{index}/Origin',430,297+index,'§9 provenance'))
JOIN_TABLE = tuple(JOIN_TABLE)


def validate_tables(joins, names):
    for side in ('py', 'wl'):
        paths = [getattr(row, side) for row in joins]
        if len(paths) != len(set(paths)):
            raise InputError('duplicate join row: ' + side)
        for i, p in enumerate(paths):
            for q in paths[i+1:]:
                if p[:len(q)] == q or q[:len(p)] == p:
                    raise InputError('overlapping join rows: ' + side + repr((p,q)))
        values = [getattr(row, side) for row in names]
        if len(values) != len(set(values)):
            raise InputError('duplicate name binding: ' + side)
    if len({r.label for r in joins}) != len(joins):
        raise InputError('duplicate join label')


# Role bindings identify actions, never their values or constitutive arguments.
# Component offsets follow the emitted ambient/tangent basis, not tree position.
@dataclass(frozen=True)
class ActionRow:
    py: str
    wl: str
    py_line: int
    wl_line: int
    spec: str


ACTION_TABLE = []
for py,wl,pl,ww,spec in (
    ('MomentumDensity','MomentumDensity',194,146,'§4 momentum density'),
    ('InternalForce','InternalForce',211,168,'§4 internal material force'),
    ('FullFaceSupportLoad','NativeCompleteHold',234,175,'§5 complete native load'),
    ('OutwardSourcePartner','OutwardAdditionalMomentumPartner',254,207,'§5 additional exchange'),
    ('FaceApplicationVelocity','NativeApplicationVelocity',270,256,'§6 native application velocity')):
    for i in range(4):
        ACTION_TABLE.append(ActionRow('OPEN_'+py+'_'+str(i),wl+'['+str(i+1)+']',pl,ww,spec))
for a in range(4):
    for i in range(3):
        ACTION_TABLE.append(ActionRow(f'OPEN_MomentumFlux_{a}_{i}',f'MomentumCurrent[{a+1},{i+1}]',195,149,'§4 momentum current'))
for i in range(3):
    ACTION_TABLE.append(ActionRow(f'OPEN_MaterialEnergyFlux_{i}',f'MaterialEnergyCurrent[{i+1}]',291,278,'§6 energy current'))
ACTION_TABLE.extend(ActionRow(*r) for r in (
    ('OPEN_NativeFaceSet','NativeBoundingFaceSet',219,179,'§5 complete native face set'),
    ('OPEN_ApplyNativeFaceReduction','NativeToCoordinateDensity',243,181,'§5,6 native reduction'),
    ('OPEN_CarriedMomentumW','OutwardCarriedBulkMomentum',249,202,'§5 bulk carried exchange'),
    ('OPEN_MaterialEnergyDensity','MaterialEnergyDensity',290,276,'§6 energy density'),
    ('OPEN_MaterialStressNormalWork','MaterialStressNormalRotationalWork',282,263,'§6 material work'),
    ('OPEN_JointPowerAccounting','JointNetPowerOccurrence',298,288,'§6 joint power'),
    ('OPEN_SumOverAllNativeFaces','held:Total',140,192,'§5,6 complete face aggregation'),
    ('OPEN_T_bulk_n_s_live','native:TBulkNormalLive',229,68,'§3.3,5 native bulk amplitude'),
    ('OPEN_MaterialCompatibility','identification:material',212,162,'§3.2,4 material compatibility'),
    ('OPEN_UnfixedNormalGeneralizedRates','operand:UnspecifiedRotationalNormalRates',277,266,'§6 generalized rates')))

ACTION_TABLE = tuple(ACTION_TABLE)
# These are secondary analyses of the declared balance objects, not new stream
# joins or new physical equations. Their citations are the same construction cuts.
BALANCE_TABLE = tuple(r for r in JOIN_TABLE if r.label in
                      ('hold_inplane','hold_bulk','hold_normal','energy_balance'))


def validate_secondary(actions,balances):
    for side in ('py','wl'):
        for table in (actions,balances):
            values = [getattr(r,side) for r in table]
            if len(values) != len(set(values)):
                raise InputError('duplicate secondary table binding: '+side)
    if len({r.label for r in balances}) != len(balances):
        raise InputError('duplicate balance label')


def children(n):
    if n.head == 'Record':
        return {e.value: e.args[0] for e in n.args}
    if n.head in ('Sequence', 'Matrix'):
        return {str(i): a for i, a in enumerate(n.args)}
    return {}


def path_children(n):
    parts = children(n)
    if parts or n.head in ('Record','Sequence','Matrix'):
        return parts
    return {'@'+str(i):v for i,v in enumerate(n.args)}


def extract(stream, path):
    value = stream.get(path[0])
    for step in path[1:]:
        if value is None:
            return None
        value = path_children(value).get(str(step))
    return value


def parsed_leaves(n):
    """Independent parse census: terminal IR nodes, including heads/options.

    Empty containers count as one; options are counted individually. Applied
    names and derivative/binder indices are part of this count, never erased.
    """
    return (sum(parsed_leaves(a) for a in n.args) if n.args else 1) + len(n.options)


def algebra_leaves(n):
    # Constructor options are retained in structure, but not consumed by scalar
    # subtraction. Do not inflate arithmetic coverage by counting those options.
    return sum(algebra_leaves(a) for a in n.args) if n.args else 1


def data(n):
    return [n.head, n.value, list(n.options), [data(a) for a in n.args]]


class Bindings:
    def __init__(self, rows):
        self.maps = {'py': {r.py: r.py for r in rows}, 'wl': {r.wl: r.py for r in rows}}
        self.kinds = {r.py: r.kind for r in rows}

    def name(self, value, engine):
        return self.maps[engine].get(value, engine + '::' + value)

    def apply(self, n, engine):
        return Node(n.head, tuple(self.apply(a, engine) for a in n.args),
                    self.name(n.value, engine) if n.head in ('Symbol','FunctionName','Dummy') else n.value,
                    n.options)


class Unsupported(Exception):
    pass


def algebra(n, bindings, engine):
    """Exact supported algebra only; no unknown heads cast to scalar symbols."""
    h, a = n.head, n.args
    if h == 'Number':
        return sp.Rational(n.value)
    if h == 'Symbol':
        key = bindings.name(n.value, engine)
        if bindings.kinds.get(key) not in ('coordinate','parameter'):
            raise Unsupported('named operand has no supplied scalar value')
        return sp.Symbol(key, real=True)
    if h in ('Add','Mul','Pow','Rational'):
        values = [algebra(x, bindings, engine) for x in a]
        return {'Add': sp.Add, 'Mul': sp.Mul, 'Pow': sp.Pow, 'Rational': sp.Rational}[h](*values)
    if h == 'Apply':
        head = a[0]
        if head.head in ('FunctionName','Symbol'):
            key = bindings.name(head.value, engine)
            if engine == 'wl' and head.value == 'Sqrt' and len(a) == 2:
                return sp.sqrt(algebra(a[1], bindings, engine))
            if bindings.kinds.get(key) == 'profile':
                if len(a) < 2:
                    raise Unsupported('unapplied profile')
                return sp.Function(key)(*(algebra(x,bindings,engine) for x in a[1:]))
        # InputForm Derivative[n][F][argument] keeps the applied argument live.
        if (engine == 'wl' and head.head == 'Apply' and head.args[0].head == 'Apply'
                and head.args[0].args[0] == atom('Symbol','Derivative')):
            orders = head.args[0].args[1:]
            function = head.args[1]
            key = bindings.name(function.value,engine)
            if len(orders) == 1 and len(a) == 2 and bindings.kinds.get(key) == 'profile':
                var = sp.Symbol('_profile_argument', real=True)
                derivative = sp.Derivative(sp.Function(key)(var),(var,int(orders[0].value)),evaluate=False)
                return sp.Subs(derivative,var,evaluation_point(a[1],bindings,engine)).doit()
        raise Unsupported('held, OPEN or unsupported applied head')
    if h == 'Derivative':
        v = algebra(a[0],bindings,engine)
        specs = []
        for x in a[1:]:
            if x.head == 'Sequence':
                specs.append((algebra(x.args[0],bindings,engine),int(x.args[1].value)))
            else:
                specs.append(algebra(x,bindings,engine))
        return sp.Derivative(v,*specs,evaluate=False)
    if h == 'Subs' and len(a) == 3:
        def members(v):
            return v.args if v.head == 'Sequence' else (v,)
        if len(members(a[1])) != len(members(a[2])) or not members(a[2]):
            raise Unsupported('substitution evaluation point absent or arity differs')
        # Dummy differentiation coordinates are local algebraic binders only.
        local = Bindings(())
        local.maps = {k:dict(v) for k,v in bindings.maps.items()}
        local.kinds = dict(bindings.kinds)
        for v in members(a[1]):
            if v.head not in ('Symbol','Dummy'):
                raise Unsupported('non-symbol substitution binder')
            local.maps[engine][v.value] = '_profile_argument'
            local.kinds['_profile_argument'] = 'coordinate'
        def convert(v):
            if v.head == 'Dummy':
                v = Node('Symbol', value=v.value)
            return algebra(v,local,engine)
        # Rewrite Dummy nodes for this bound scope, preserving all orders/points.
        def dummy(v):
            return Node('Symbol' if v.head == 'Dummy' else v.head,
                        tuple(dummy(x) for x in v.args),v.value,v.options)
        return sp.Subs(algebra(dummy(a[0]),local,engine),
                       tuple(convert(x) for x in members(a[1])),
                       tuple(evaluation_point(x,bindings,engine) for x in members(a[2]))).doit()
    raise Unsupported('non-algebraic ' + h)


def evaluation_point(n,bindings,engine):
    """Evaluation points are operands, including when they are composite."""
    return algebra(n,bindings,engine)


def layout(n):
    if n.head == 'Matrix':
        if all(row.head == 'Sequence' and len(row.args) == 1 for row in n.args):
            return seq(*(row.args[0] for row in n.args))
        return seq(*n.args)
    return n


def fingerprint(n):
    # Merkle representation retains argument order inside each engine; it never
    # pairs one engine's child with the other engine's child by that position.
    return hashlib.sha256(json.dumps([n.head,n.value,n.options,
        [fingerprint(x) for x in n.args]],separators=(',',':')).encode()).hexdigest()


def bag_delta(a,b):
    return [[key,a.get(key,0)-b.get(key,0)] for key in sorted(a.keys() | b.keys())
            if a.get(key,0) != b.get(key,0)]


def differences(a,b):
    """Unpaired subtree fingerprints, not a positional cross-engine alignment."""
    aa,bb = Counter({fingerprint(a):1}),Counter({fingerprint(b):1})
    for digest,count in bag_delta(aa,bb):
        yield {'subtree_sha256':digest,'py_minus_wl_count':count,
               'correspondence':'unpaired subtree; complete mapped operand printed above'}


def wl_structural_role(n):
    """Specific constructor forms cited by ACTION_TABLE, including non-OpenAction
    representations. Others remain unpaired; no general head alias is inferred.
    """
    def name(v):
        return v.value.removeprefix('wl::')
    def operand_name(v):
        if v.head=='Apply' and name(v.args[0])=='OPEN' and len(v.args)>1:
            return name(v.args[1])
        return ''
    if n.head!='Apply':
        return None
    head=n.args[0]
    if head.head=='Apply':
        if name(head.args[0])=='Inactive' and len(head.args)==2 and name(head.args[1])=='Total':
            return 'held:Total'
        if name(head.args[0])=='OpenNativeField' and len(head.args)==2:
            operand=operand_name(head.args[1])
            if operand:
                return 'native:'+operand
    if name(head)=='UnresolvedIdentification' and len(n.args)==5 and operand_name(n.args[1]) in ('NBrLive','N_br_live'):
        return 'identification:material'
    if name(head)=='OPEN' and len(n.args)>1 and name(n.args[1])=='UnspecifiedRotationalNormalRates':
        return 'operand:UnspecifiedRotationalNormalRates'
    return None


def action_key(n,engine):
    if n.head != 'Apply':
        return None
    head = n.args[0]
    if engine == 'py' and head.value.startswith('OPEN_'):
        return head.value
    if engine == 'wl' and head == atom('Symbol','OpenAction') and len(n.args)>1:
        role = n.args[1]
        if role.head == 'Symbol':
            return role.value
        if role.head == 'Apply' and all(v.head=='Number' for v in role.args[1:]):
            return role.args[0].value+'['+','.join(v.value for v in role.args[1:])+']'
        return json.dumps(data(role),separators=(',',':'))
    return wl_structural_role(n) if engine=='wl' else None


def action_occurrences(n,include_structural=False):
    """Yield mapped OPEN nodes with balance orientation or argument status."""
    occurrences = []
    def signed_actions(v,sign=1,dependency=False):
        # Linear syntax carries the enclosing expression's orientation. This is
        # the engines' held calculus/aggregation syntax, not a constitutive law:
        # PY audit 134-140, 242-244; WL audit 142-143, 181-193.
        if v.head == 'Mul':
            for child in v.args:
                if child.head == 'Number' and child.value.startswith('-'):
                    sign *= -1
            for child in v.args:
                if child.head != 'Number':
                    signed_actions(child,sign,dependency)
        elif v.head == 'Apply':
            head = v.args[0]
            linear_arguments = ()
            if (head.head == 'Apply' and len(head.args) == 2
                    and head.args[0].value in ('Inactive','wl::Inactive')):
                operator = head.args[1].value.removeprefix('wl::')
                if operator in ('Total','Map','D'):
                    linear_arguments = (1,)
            elif head.value.removeprefix('wl::') in ('Function','OpenFirstVariation'):
                linear_arguments = (2,) if head.value.removeprefix('wl::') == 'Function' else (1,)
            elif head.value.removeprefix('py::') == 'OPEN_SumOverAllNativeFaces':
                linear_arguments = (2,)
            role = None
            if head.value.removeprefix('py::').startswith('OPEN_'):
                role = head.value
            elif head.value.removeprefix('wl::') == 'OpenAction' and len(v.args) > 1:
                role = json.dumps(data(v.args[1]),separators=(',',':'))
            if role is None and include_structural:
                role = wl_structural_role(v)
            if role is not None:
                occurrences.append((v,role,sign,dependency))
            # Unknown/OPEN action arguments are dependencies. Only explicitly
            # identified linear bodies inherit the sign; Map's set and binders
            # do not. Dependency mode remains distinct even inside a wrapper.
            for i,child in enumerate(v.args[1:],1):
                if i in linear_arguments:
                    signed_actions(child,sign,dependency)
                else:
                    signed_actions(child,1,True)
        elif v.head in ('Derivative','Lambda'):
            body = 0 if v.head == 'Derivative' else 1
            for i,child in enumerate(v.args):
                signed_actions(child,sign if i == body else 1,dependency or i != body)
        else:
            for child in v.args:
                signed_actions(child,sign if v.head in ('Add','Sequence','Matrix') else 1,dependency)
    signed_actions(n)
    return occurrences


def structure(n):
    names, heads, options, orientations = Counter(), Counter(), Counter(), Counter()
    argument_actions = Counter()
    def visit(v):
        heads[v.head] += 1
        if v.head in ('Symbol','FunctionName','Dummy'):
            names[v.value] += 1
        options.update(v.options)
        for child in v.args:
            visit(child)
    visit(n)
    for v,role,sign,dependency in action_occurrences(n):
        if dependency:
            argument_actions[role] += 1
        else:
            orientations[(role,sign)] += 1
    terms = n.args if n.head == 'Add' else (n,)
    factors_by_term = []
    for term in terms:
        factors = term.args if term.head == 'Mul' else (term,)
        factors_by_term.append([v.value for v in factors if v.head == 'Number'] or ['implicit +1'])
    return {'named_operands_and_heads':dict(sorted(names.items())),
            'node_heads':dict(sorted(heads.items())),
            'constructor_options':[[list(k),v] for k,v in sorted(options.items())],
            'term_numeric_factors':factors_by_term,
            'open_action_orientation_counts':[[role,sign,count] for (role,sign),count in sorted(orientations.items())],
            'open_argument_occurrence_counts':[[role,count] for role,count in sorted(argument_actions.items())]}


def role_maps(actions):
    return {'py':{r.py:r.py for r in actions},'wl':{r.wl:r.py for r in actions}}


def object_key(n):
    """Readable semantic key, separate from the lossless serialization census."""
    return json.dumps(data(n),separators=(',',':'))


def feature_tree(n,bindings,engine,scope=None):
    """Canonical syntax for object inventories, never a scalar residual.

    Only declared names are shared. Normalize constructor spelling (PY audit
    60-74, 152; WL audit 22-31), rational arithmetic, and profile derivatives.
    Unknown heads stay engine-qualified. No OPEN action or binder is evaluated.
    Original metadata/structure remain in the lossless trees and their deltas.
    """
    scope = {} if scope is None else scope
    def convert(v,local=scope):
        return feature_tree(v,bindings,engine,local)
    def members(v):
        return v.args if v.head=='Sequence' else (v,)
    h,a=n.head,n.args
    if h in ('Symbol','FunctionName','Dummy'):
        return scope.get(n.value,atom('Name',bindings.name(n.value,engine)))
    if h=='Number':
        return atom('Number',Fraction(n.value))
    if h=='Apply':
        head=a[0]
        if engine=='wl' and head.head=='Prime' and len(a)==2:
            order=0
            function=head
            while function.head=='Prime':
                order+=1
                function=function.args[0]
            key=bindings.name(function.value,engine)
            if bindings.kinds.get(key)=='profile':
                return Node('ProfileDerivative',(seq(atom('Number',order)),seq(convert(a[1]))),key)
        # Derivative[n1,...][F][points...] is one derivative application,
        # including its head: never visit F or its placeholders as operands.
        if (engine=='wl' and head.head=='Apply' and head.args[0].head=='Apply'
                and head.args[0].args[0]==atom('Symbol','Derivative')
                and len(head.args)==2):
            key=bindings.name(head.args[1].value,engine)
            orders=head.args[0].args[1:]
            if bindings.kinds.get(key)=='profile' and len(orders)==len(a)-1:
                return Node('ProfileDerivative',(seq(*(convert(v) for v in orders)),
                            seq(*(convert(v) for v in a[1:]))),key)
        if head.head in ('Symbol','FunctionName'):
            key=bindings.name(head.value,engine)
            if bindings.kinds.get(key)=='profile':
                return Node('LiveProfile',tuple(convert(v) for v in a[1:]),key)
            if engine=='wl' and head.value=='Sqrt' and len(a)==2:
                return convert(Node('Pow',(a[1],atom('Number','1/2'))))
            if engine=='wl' and head.value=='Function' and len(a)==3:
                return convert(Node('Lambda',a[1:]))
            if engine=='wl' and head.value=='OpenAction':
                # The role has its own declared table; it is not an operand.
                return Node('ActionArguments',tuple(convert(v) for v in a[2:]))
            if engine=='wl' and head.value=='OpenFirstVariation':
                return Node('FirstVariation',tuple(convert(v) for v in a[1:]))
        if (engine=='wl' and head.head=='Apply' and len(head.args)==2
                and head.args[0]==atom('Symbol','Inactive') and head.args[1]==atom('Symbol','D')):
            return Node('HeldDerivative',tuple(convert(v) for v in a[1:]))
        # Heads may themselves contain derivatives, native fields, or binders.
        return Node('Apply',tuple(convert(v) for v in a))
    if h=='Derivative' and a:
        body=convert(a[0])
        if body.head=='LiveProfile':
            orders=[0]*len(body.args)
            for spec in a[1:]:
                parts=members(spec)
                variable=convert(parts[0])
                if len(parts)>2 or body.args.count(variable)!=1:
                    break
                order=int(parts[1].value) if len(parts)==2 and parts[1].head=='Number' else 1
                orders[body.args.index(variable)] += order
            else:
                return Node('ProfileDerivative',(seq(*(atom('Number',v) for v in orders)),
                            seq(*body.args)),body.value)
        return Node('Derivative',tuple(convert(v) for v in a))
    if h in ('Subs','Lambda'):
        variables=members(a[1] if h=='Subs' else a[0])
        local=dict(scope)
        placeholders=[]
        for i,v in enumerate(variables):
            placeholder=atom('Bound',f'{len(scope)}:{i}')
            local[v.value]=placeholder
            placeholders.append(placeholder)
        body=convert(a[0] if h=='Subs' else a[1],local)
        if h=='Subs':
            points=tuple(convert(v) for v in members(a[2]))
            if body.head=='ProfileDerivative' and len(points)==len(placeholders) and points:
                replacements=dict(zip(placeholders,points))
                def substitute(v):
                    return replacements.get(v,Node(v.head,tuple(substitute(x) for x in v.args),v.value))
                return substitute(body)
            return Node('Subs',(seq(*placeholders),body,seq(*points)))
        return Node('Lambda',(seq(*placeholders),body))
    args=tuple(convert(v) for v in a)
    if h=='Rational':
        return atom('Number',Fraction(args[0].value)/Fraction(args[1].value))
    if h in ('Add','Mul'):
        flat=[x for v in args for x in (v.args if v.head==h else (v,))]
        numbers=[Fraction(v.value) for v in flat if v.head=='Number']
        rest=[v for v in flat if v.head!='Number']
        number=Fraction(0 if h=='Add' else 1)
        for v in numbers:
            number=number+v if h=='Add' else number*v
        if number!=(0 if h=='Add' else 1) or not rest:
            rest.append(atom('Number',number))
        args=tuple(sorted(rest,key=object_key))
        if len(args)==1:
            return args[0]
    if h=='Pow' and len(args)==2 and all(v.head=='Number' for v in args):
        power=Fraction(args[1].value)
        if power.denominator==1 and (Fraction(args[0].value) or power>=0):
            return atom('Number',Fraction(args[0].value)**int(power))
    return Node(h,args,n.value)


OBJECT_FIELDS = frozenset(('named_OPEN_operands','live_arguments','binders'))


def features(n,bindings,engine):
    """Inventories of objects, not counts of engine-specific syntax nodes."""
    mapped = bindings.apply(n,engine)
    names,live,binders = Counter(),Counter(),Counter()
    def visit(v,head=False):
        if v.head=='Name' and not head and bindings.kinds.get(v.value) not in ('coordinate','parameter','profile','binder'):
            names[v.value] = 1
        if v.head in ('LiveProfile','ProfileDerivative'):
            live[object_key(v)] = 1
        if v.head in ('ProfileDerivative','Derivative','Subs','Lambda','FirstVariation','HeldDerivative'):
            binders[object_key(v)] = 1
        for i,x in enumerate(v.args):
            visit(x,head=v.head=='Apply' and i==0)
    visit(feature_tree(n,bindings,engine))
    return {'head':Counter({json.dumps([mapped.head,mapped.value,
                data(mapped.args[0]) if mapped.head=='Apply' else None],separators=(',',':')):1}),
            'named_OPEN_operands':names,'live_arguments':live,'binders':binders,
            'argument_trees':Counter({fingerprint(mapped):1})}


def merge_fields(current,fields):
    for field,values in fields.items():
        inventory=current.setdefault(field,Counter())
        if field in OBJECT_FIELDS:
            inventory.update({key:1 for key in values if key not in inventory})
        else:
            inventory.update(values)


def field_deltas(left,right):
    return {field:bag_delta(left.get(field,{}),right.get(field,{}))
            for field in sorted(left.keys() | right.keys())}


def action_comparison(a,b,bindings,actions):
    maps = role_maps(actions)
    grouped = {}
    for engine,node in (('py',a),('wl',b)):
        group = {}
        for original,_,sign,dependency in action_occurrences(node,include_structural=True):
            key = action_key(original,engine)
            role = maps[engine].get(key,engine+'::'+key)
            fields = features(original,bindings,engine)
            fields['head'] = Counter({key:1})
            fields['orientation'] = Counter({('argument' if dependency else str(sign)):1})
            fields['role'] = Counter({role:1})
            current = group.setdefault(role,{})
            merge_fields(current,fields)
        grouped[engine] = group
    result=[]
    for role in sorted(grouped['py'].keys() | grouped['wl'].keys()):
        left,right = grouped['py'].get(role,{}),grouped['wl'].get(role,{})
        result.append({'role':role,'py':left,'wl':right,
                       'unpaired_reason':None if left and right else
                           ('no cross-engine role binding; representation-specific action' if '::' in role else
                            'no occurrence of this declared role in the other operand'),
                       'differences':field_deltas(left,right)})
    return result


def additive_terms(n):
    """Distribute only syntactic sums/products; do not evaluate held bodies."""
    if n.head=='Add':
        return [term for v in n.args for term in additive_terms(v)]
    if n.head=='Mul':
        terms=[()]
        for v in n.args:
            terms=[(*prefix,w) for prefix in terms for w in additive_terms(v)]
        return [Node('Mul',t) for t in terms]
    return [n]


def contains_open(n,engine):
    if action_key(n,engine) is not None:
        return True
    return any(contains_open(v,engine) for v in n.args)


def term_orientation(n):
    sign=1
    if n.head=='Number' and n.value.startswith('-'):
        return -1
    if n.head=='Mul':
        for v in n.args:
            sign *= term_orientation(v)
    elif n.head=='Derivative':
        sign *= term_orientation(n.args[0])
    elif n.head=='Lambda':
        sign *= term_orientation(n.args[1])
    elif n.head=='Apply':
        head=n.args[0]
        if head.value in ('OpenFirstVariation','Function','OPEN_SumOverAllNativeFaces'):
            body=1 if head.value=='OpenFirstVariation' else 2
            sign *= term_orientation(n.args[body])
        elif (head.head=='Apply' and head.args[0].value=='Inactive' and
              head.args[1].value in ('Total','Map','D')):
            sign *= term_orientation(n.args[1])
    return sign


def balance_comparison(a,b,bindings,actions):
    maps=role_maps(actions)
    sides,closed,groups = {},{},{}
    for engine,node in (('py',a),('wl',b)):
        entries=[]
        closed_terms=[]
        group={}
        for term in additive_terms(node):
            is_open=contains_open(term,engine)
            if not is_open:
                closed_terms.append(term)
            mapped=bindings.apply(term,engine)
            # Roles are a multiset of the declared actions in this term. Unknown
            # actions remain engine-qualified; the OPEN-free group is the closed
            # portion of this declared balance, not an assumed physical closure.
            def roles(v):
                key=action_key(v,engine)
                if key is not None:
                    return [maps[engine].get(key,engine+'::'+key)]
                return [r for child in v.args for r in roles(child)]
            role=json.dumps(sorted(roles(term))) if is_open else 'OPEN_free'
            fields=features(term,bindings,engine)
            fields['orientation']=Counter({str(term_orientation(term)):1})
            fields['role']=Counter({role:1})
            entries.append({'role':role,'orientation':term_orientation(term),
                            'open_free':not is_open,'operand':data(mapped)})
            current=group.setdefault(role,{})
            merge_fields(current,fields)
        sides[engine]=entries
        groups[engine]=group
        closed[engine]=Node('Add',tuple(closed_terms)) if closed_terms else atom('Number',0)
    result={'outcome':'balance','entries':sides,'entry_differences':{
        role:field_deltas(groups['py'].get(role,{}),groups['wl'].get(role,{}))
        for role in sorted(groups['py'].keys() | groups['wl'].keys())},
        'closed_operands':{e:data(n) for e,n in closed.items()}}
    # The caller emits these operands before forming this supplemental residual.
    return result,closed


def emit_balance(output,label,a,b,bindings,actions):
    def emit_component(a,b,component):
        a,b=layout(a),layout(b)
        if a.head=='Sequence' or b.head=='Sequence':
            if a.head!=b.head or len(a.args)!=len(b.args):
                emit(output,{'kind':'balance_comparison','row':label,'component':component,
                             'comparison':{'outcome':'not_formed','reason':'component layout differs'}})
                return
            for i in range(len(a.args)):
                emit_component(a.args[i],b.args[i],(*component,i))
            return
        result,closed=balance_comparison(a,b,bindings,actions)
        emit(output,{'kind':'balance_entries','row':label,'component':component,**result})
        comparison,counts=compare(closed['py'],closed['wl'],bindings)
        emit(output,{'kind':'balance_comparison','row':label,'component':component,
                     'comparison':{'outcome':'balance','entry_differences':result['entry_differences'],
                                   'closed_residual':comparison,'closed_compared_leaves':counts}})
    emit_component(a,b,())


def compare(a,b,bindings):
    result,counts = _compare(a,b,bindings)
    if result['outcome'] == 'not_formed':
        result['structural_delta'] = list(differences(bindings.apply(a,'py'),bindings.apply(b,'wl')))
    return result,counts


def _compare(a,b,bindings):
    """Three outcomes: exact residual; not formed; absent counterpart.

    No text/boolean/name equality is promoted to algebraic zero. Container
    siblings are independent, so a boolean cannot suppress an algebraic sibling.
    """
    a,b = layout(a),layout(b)
    ca,cb = children(a),children(b)
    if ca or cb or a.head in ('Record','Sequence') or b.head in ('Record','Sequence'):
        if a.head != b.head:
            return {'outcome':'not_formed','reason':'container structure differs'},(0,0)
        result, counts = {}, [0,0]
        for key in sorted(ca.keys() | cb.keys()):
            if key not in ca or key not in cb:
                result[key] = {'outcome':'not_formed','reason':'nested sibling absent',
                               'py_present':key in ca,'wl_present':key in cb}
            else:
                result[key],count = compare(ca[key],cb[key],bindings)
                counts[0] += count[0]
                counts[1] += count[1]
        return {'outcome':'container','children':result},tuple(counts)
    if a.head in ('Text','Boolean','Name','FunctionName') or b.head in ('Text','Boolean','Name','FunctionName'):
        return {'outcome':'not_formed','reason':'text, native boolean or name is not a subtractable operand'},(0,0)
    # Relational operands retain the relation head, with each side compared.
    relations = {'Equality','Unequality','StrictGreaterThan','StrictLessThan','GreaterThan','LessThan'}
    if a.head in relations or b.head in relations:
        if a.head != b.head or len(a.args) != len(b.args):
            return {'outcome':'not_formed','reason':'relational structure differs'},(0,0)
        result, counts = [], [0,0]
        for x,y in zip(a.args,b.args):
            r,c = compare(x,y,bindings)
            result.append(r)
            counts = [counts[i]+c[i] for i in range(2)]
        return {'outcome':'relation','head':a.head,'operands':result},tuple(counts)
    if a.head == b.head == 'Symbol':
        return {'outcome':'not_formed','reason':'bindings alone are not value evidence'},(0,0)
    try:
        left,right = algebra(a,bindings,'py'),algebra(b,bindings,'wl')
    except Unsupported as exc:
        return {'outcome':'not_formed','reason':str(exc)},(0,0)
    # cancel treats intact applied functions and derivative applications as
    # algebraic generators internally, without replacing/printing bare symbols.
    delta = sp.cancel(left-right)
    return {'outcome':'exact','py_value':str(left),'wl_value':str(right),
            'residual':str(delta),'points':[],
            'point_policy':'exact symbolic subtraction; no numeric samples'},(algebra_leaves(a),algebra_leaves(b))


def remainder(stream, paths):
    """Partition every parsed stream independently around the declared join cuts."""
    def walk(n,path):
        if path in paths:
            return
        if any(p[:len(path)] == path for p in paths):
            parts = path_children(n)
            if not parts:
                yield path,n
            else:
                for key,value in parts.items():
                    yield from walk(value,(*path,key))
        else:
            yield path,n
    for tag,value in stream.items():
        yield from walk(value,(tag,))


# Explicit explanations for engine-only emission families. More specific
# explanations are installed below for the remaining physical occurrences.
UNJOINED_REASONS = {
    'py': {
        'BASIS/graph_dual':'WL emits metric inverse and tangents, but no dual-tangent matrix; constructing their product is outside comparison (§1).',
        'MEASURES/graph_area_factor':'WL emits det(g), but no separate graph-area-factor object; taking its square root would construct a counterpart (§1).',
        'MEASURES/native_chart_domain':'PY has an OPEN chart-domain action; WL uses an abstract native point and unrestricted section, with no chart-domain action (§5).',
        'MEASURES/native_reduction_map':'PY emits a reduction-map instance separately from application; WL emits reduction applications, with no separate map-instance object (§5).',
        'MEASURES/chart_coordinates':'WL uses the abstract point pNative, with no three chart-coordinate objects (§5).',
        'MEASURES/orientation':'PY declares an outward regular native chart; WL has no native chart-orientation declaration (§5).',
        'GEOMETRY/embedding':'WL constructs graphEmbedding internally but does not emit the embedding vector; emitted tangent/identity objects are joined (§1).',
        'GEOMETRY/native_embedding':'WL native geometry is an OPEN map action, with no immersion-coordinate vector (§5).',
        'GEOMETRY/native_cofactor':'WL has no native cofactor vector; its native normal is an OPEN map action (§5).',
        'GEOMETRY/native_regular_domain':'WL has no immersion-chart regularity inequality (§5).',
        'GEOMETRY/native_normal_norm':'WL has no separately computed native-normal norm (§5).',
        'GEOMETRY/native_normal_tangent_pairing':'WL has no separately computed native normal/tangent contractions (§5).',
        'MATERIAL_INPUT_DIFFERENTIALS/material_velocity_rate':'WL emits profile material derivatives and section tangents, but no contracted material acceleration vector (§1,4).',
        'MATERIAL_POWER_PAIRING/graph_velocity_contraction':'WL retains force and velocity as arguments of joint work; it emits no standalone graph-velocity contraction (§6).',
        'MATERIAL_POWER_PAIRING/contraction_scope':'WL has no separate qualification for a standalone graph-velocity contraction (§6).',
        'MATERIAL_POWER_PAIRING/rotational':'WL retains rotational work within its joint material-work action, with no separately attributable rotational-work action (§6).',
        'FACE_POWER_PAIRING/native_reduction_map':'WL has no separately emitted reduction-map instance for work, only its application (§5,6).',
        'FACE_POWER_PAIRING/native_area':'Additional PY native-area occurrence; WL keeps native measure within work-reduction applications, with no additional separate area object (§5,6).',
        'COUPLED_INPUTS/live_density_name':'PY density-name collision explanation; WL has no equivalent collision-audit object.',
        'COUPLED_INPUTS/native_reduction':'PY standalone no-sheet/slab-choice declaration; WL encodes the reduction as an OPEN application without a separate declaration (§5).',
        'MODEL_POINT/premise_status':'PY separately repeats the adopted-premise label/date; WL includes its date only inside the joined Premises object (§2).',
        'COUPLED_INPUTS/operands/24':'PY separately registers the S12 reaction system; WL carries it within additional momentum-partner actions, with no second register occurrence (§5).',
        'COUPLED_INPUTS/operands/25':'PY separately registers S12 energy reaction/supply; WL has no separately named energy-reaction/supply operand (§6).',
        'COUPLED_INPUTS/operands/26':'PY separately registers mouth/core data; WL carries boundary/core operands but no separate mouth/core-data object (§3.3).',
        'COUPLED_INPUTS/operands/27':'PY separately registers bulk_state; WL emits EOS inputs and unrestricted native state, no separate general bulk-state operand (§1).',
        'COUPLED_INPUTS/operands/31':'PY separately registers face-support partition; WL carries unresolved support partition within complete-hold actions, no separate register occurrence (§3.3).',
    },
    'wl': {
        'BASIS_MEASURES_GEOMETRY/Domain':'WL explicitly emits the r>0 domain; PY emits no corresponding inequality object (§1).',
        'INTERNAL_MATERIAL_FORCE/GraphNormalProjection':'PY emits internal-force components, but no separate normal projection of internal force (§4).',
        'MECHANICAL_LOAD/NativeBulkAmplitude':'WL separately emits its native bulk-amplitude field; PY emits amplitude only inside joined bulk-traction components (§5).',
        'MECHANICAL_LOAD/NativeHoldOperand':'WL separately emits its native hold field; PY emits the descriptor only within complete-hold actions (§5).',
        'MECHANICAL_LOAD/GraphNormalProjection':'PY emits mechanical-load components, but no separate mechanical-load normal projection (§4,5).',
        'EXCHANGE_MOMENTUM/Orientation':'WL separately declares outward loss; PY carries it in the joined carried action/trace, no standalone declaration (§5).',
        'DRIVE_PROVENANCE/GMInterface':'WL separately emits S16 matching dependency; PY includes it only in the joined interface list (§10).',
        'B_HOLD_LIVE/Status':'WL separately labels its named balance; PY has no separate conditional-balance label object (§9).',
        'MASS_INPUT/VelocityOperand':'WL separately repeats the in-plane mass velocity; PY emits it only as part of mass current and graph velocity (§3.1).',
        'MASS_INPUT/RHS':'Mass-law right-hand side, not outward loss: PY has no separately emitted same-role RHS occurrence outside the already joined mass equation (§3.1,5; PY audit 266,378; WL audit 244,247).',
        'MASS_INPUT/Qualification':'WL separately declares the mass-law qualification; PY includes it in the joined model restriction inventory (§7).',
        'MASS_INPUT/NativeIdentification':'WL separately repeats the O6 operand for mass; PY has no additional mass-native-identification entry (§5).',
        'COUPLED_INPUTS_MODEL_POINT/BulkInputs':'WL emits the three bulk EOS/sound/contrast relations; PY emits f and bulk_state but no bulk EOS relation objects (§3.1).',
        'COUPLED_INPUTS_MODEL_POINT/Anchoring':'WL emits a structured LABHELD object; PY includes its declaration within the joined restriction inventory (§3.1,7).',
        'COUPLED_INPUTS_MODEL_POINT/Counting/Truncation':'WL explicitly emits a no-truncation token; PY includes its declaration in the joined restriction inventory (§7).',
        'COUPLED_INPUTS_MODEL_POINT/TransferLimits':'Additional structured WL transfer-limit list; PY has one combined restriction list, already joined, with no second transfer-limit object (§7).',
        'COUPLED_INPUTS_MODEL_POINT/HistoricalDomains':'Additional structured WL historical-domain record; PY carries its historical declarations in the joined restriction list, with no second domain record (§8).',
        'COUPLED_INPUTS_MODEL_POINT/O4EquationIdentityCount':'Additional WL unsettled-count token; PY carries it within the joined coupled-embedding object, with no separate count token (§3.2).',
        'COUPLED_INPUTS_MODEL_POINT/RelaxationOwner':'Additional WL unassigned-owner token; PY carries ownership within the joined interface list, with no separate owner token (§10).',
    },
}


def unjoined_reason(engine,path):
    exact = UNJOINED_REASONS[engine].get('/'.join(path))
    if exact:
        return exact
    if path[0].startswith('LOCAL/'):
        return 'SymPy fold/publication diagnostic; Wolfram has no fold/publication object (§9).'
    if path[-1] == 'Origin' or path[0] == 'TRACE':
        return 'Additional emission-level provenance occurrence; the other engine has no additional independently paired provenance record at this granularity (§9).'
    if engine == 'wl' and path[:2] in (('B_HOLD_LIVE','Entries'),('B_E_STEADY','Entries')):
        return 'Additional WL accounting-entry occurrence or remaining section/entry metadata; PY emits assembled components and actions, not a separate entry-list copy (§4,6).'
    if engine == 'wl' and path[:2] == ('MATERIAL_MOMENTUM','DifferentiatedSection'):
        return 'Remaining WL unrestricted-section occurrence; PY has no additional separately emitted section copy beyond the joined profiles, gradients and OPEN-action arguments (§4).'
    if engine == 'wl' and path[:2] == ('MATERIAL_MOMENTUM','MaterialProfileDerivatives'):
        return 'WL emits contracted material derivatives of scalar profiles; PY emits partial gradients and material acceleration, not these contracted scalar objects (§1,4).'
    if engine == 'wl' and path[:2] == ('EXCHANGE_MOMENTUM','MaterialIdentification'):
        return 'Remaining WL same-exchanged-material descriptor; PY has no separate descriptor beyond the joined exchange actions (§5).'
    return 'No declared object counterpart; additional occurrence requires source-level coverage review: ' + '/'.join(path)


def emit(stream, value):
    stream.write(json.dumps(value,ensure_ascii=False,separators=(',',':')) + '\n')
    stream.flush()


def run(py_path,wl_path,output,accounting,joins=JOIN_TABLE,names=NAME_TABLE,
        actions=ACTION_TABLE,balances=BALANCE_TABLE):
    validate_tables(joins,names)
    validate_secondary(actions,balances)
    bindings = Bindings(names)
    py,wl = read_stream(py_path,'py'),read_stream(wl_path,'wl')
    partition = {'py':Counter(),'wl':Counter()}
    compared = {'py':Counter(),'wl':Counter()}
    for row in joins:
        a,b = extract(py,row.py),extract(wl,row.wl)
        info = {'row':row.label,'py_path':row.py,'wl_path':row.wl,
                'accounting':'joined' if a is not None and b is not None else 'counterpart_absent',
                'parsed_leaves':[parsed_leaves(a) if a is not None else 0,
                                 parsed_leaves(b) if b is not None else 0]}
        # Print complete operands before computing any residual or count guard.
        emit(output,{'kind':'operands',**info,'py':data(a) if a else None,'wl':data(b) if b else None})
        if a is None or b is None:
            result,counts = {'outcome':'not_formed','reason':'declared counterpart absent'},(0,0)
        else:
            aa,bb = bindings.apply(a,'py'),bindings.apply(b,'wl')
            emit(output,{'kind':'structure','row':row.label,'py':structure(aa),'wl':structure(bb),
                         'mapped_py':data(aa),'mapped_wl':data(bb),
                         'differences':list(differences(aa,bb)),
                         'action_comparison':action_comparison(a,b,bindings,actions)})
            if row.layout == 'transpose' and a.head == 'Matrix':
                a = seq(*(seq(*(r.args[i] for r in a.args)) for i in range(len(a.args[0].args))))
            result,counts = compare(a,b,bindings)
        emit(output,{'kind':'residual','row':row.label,'comparison':result})
        info['compared_leaves'] = counts
        for i,engine in enumerate(('py','wl')):
            tag = getattr(row,engine)[0]
            partition[engine][tag] += info['parsed_leaves'][i]
            compared[engine][tag] += counts[i]
        emit(output,{'kind':'accounting',**info})
        emit(accounting,info)
    for row in balances:
        # Secondary analyses use their declared paths, and emit their own operands.
        a,b=extract(py,row.py),extract(wl,row.wl)
        if a is not None or b is not None:
            emit(output,{'kind':'balance_operands','row':row.label,
                         'py_path':row.py,'wl_path':row.wl,
                         'py':data(a) if a else None,'wl':data(b) if b else None})
            if a is not None and b is not None:
                emit_balance(output,row.label,a,b,bindings,actions)
            else:
                emit(output,{'kind':'balance_comparison','row':row.label,
                             'comparison':{'outcome':'not_formed','reason':'balance counterpart absent'}})
    for engine,stream in (('py',py),('wl',wl)):
        for path,value in remainder(stream,{getattr(row,engine) for row in joins}):
            info = {'accounting':'unjoined','engine':engine,'path':path,
                    'reason':unjoined_reason(engine,path),'parsed_leaves':parsed_leaves(value),
                    'compared_leaves':0}
            emit(output,{'kind':'unjoined',**info,'operand':data(value)})
            emit(accounting,info)
            partition[engine][path[0]] += info['parsed_leaves']
        for tag,value in stream.items():
            summary = {'kind':'stream_object_count','engine':engine,'tag':tag,
                       'parsed_leaves':parsed_leaves(value),'accounted_leaves':partition[engine][tag],
                       'compared_leaves':compared[engine][tag],
                       'unaccounted_leaves':parsed_leaves(value)-partition[engine][tag]}
            emit(output,summary)
            emit(accounting,summary)
    return py,wl


def catalog(path,joins=JOIN_TABLE,names=NAME_TABLE,actions=ACTION_TABLE,balances=BALANCE_TABLE):
    with Path(path).open('w') as out:
        emit(out,{'path_convention':'zero-based sequence indices; record keys; @n constructor argument (@0 applied head); tag first',
                  'leaf_convention':'terminal parsed IR nodes plus constructor options; empty container one',
                  'feature_convention':'named_OPEN_operands, live_arguments and binders are sets of canonical objects (1=present); raw serialization counts/trees are separate; argument order retained within each live object; unknown heads stay engine-qualified',
                  'joins':[dict(label=r.label,py=r.py,wl=r.wl,spec=r.spec,layout=r.layout,
                                py_citation=f'{PY_SOURCE}:{r.py_line}',wl_citation=f'{WL_SOURCE}:{r.wl_line}')
                           for r in joins],
                  'actions':[dict(py=r.py,wl=r.wl,spec=r.spec,py_citation=f'{PY_SOURCE}:{r.py_line}',wl_citation=f'{WL_SOURCE}:{r.wl_line}') for r in actions],
                  'balances':[dict(label=r.label,py=r.py,wl=r.wl,spec=r.spec,py_citation=f'{PY_SOURCE}:{r.py_line}',wl_citation=f'{WL_SOURCE}:{r.wl_line}') for r in balances],
                  'names':[dict(py=r.py,wl=r.wl,spec=r.spec,kind=r.kind,
                                py_citation=f'{PY_SOURCE}:{r.py_line}',wl_citation=f'{WL_SOURCE}:{r.wl_line}')
                           for r in names]})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--py',type=Path,default=DEFAULT_PY)
    p.add_argument('--wl',type=Path,default=DEFAULT_WL)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--accounting',type=Path,required=True)
    p.add_argument('--catalog',type=Path,required=True)
    args = p.parse_args()
    start = time.monotonic()
    try:
        catalog(args.catalog)
        with args.output.open('x') as out,args.accounting.open('x') as accounting:
            run(args.py,args.wl,out,accounting)
        print(json.dumps({'runtime_seconds':time.monotonic()-start,
                          'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))
        return 0
    except (InputError,OSError) as exc:
        print('operational_error: ' + str(exc),file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
