#!/usr/bin/env python3
"""Directive-owned O2 K1--K13 (K6a--c) single-site source ablations.

Run ONLY through scripts/s11c_guarded_run.py. No time limit, retries or physical
expected values. The baseline invokes the live engine unchanged and publishes
its delta. Each copy changes exactly one construction statement. All section-9
payloads are reparsed from literal emission; every baseline/corrupted/difference
triple is printed, including zeros. Scratch contains sources and native streams.
"""
import argparse
import hashlib
import importlib.util
from pathlib import Path
import sys
import contextlib

sys.dont_write_bytecode = True
import sympy as sp
from sympy.core.symbol import Str
from sympy.functions.elementary.piecewise import ExprCondPair

HERE = Path(__file__).resolve()
LIVE = HERE.with_name('O2_live_balance_sympy_audit.py')

# Full statement substitutions; no result/tag/expected-value mutation.
KNIVES = {
 'K1': ('    mass_density = rho  # KNIFE_K1',
        "    mass_density = sp.Symbol('rho_mass_constant')  # KNIFE_K1"),
 'K2': ('    graph_velocity = tangent * v  # KNIFE_K2',
        '    graph_velocity = vector((*v, sp.S.Zero))  # KNIFE_K2'),
 'K3': ('    carried_v = v  # KNIFE_K3',
        "    carried_v = vector(sp.Function('independent_carried_' + str(i))(*x, t) for i in range(3))  # KNIFE_K3"),
 'K4': ('    metric = sp.ImmutableMatrix(tangent.T * tangent)  # KNIFE_K4',
        '    metric = sp.ImmutableMatrix(sp.eye(3))  # KNIFE_K4'),
 'K5': ("    ref_stress = op['R_ref_strain_live']  # KNIFE_K5",
        "    ref_stress = Str('reference_history_removed')  # KNIFE_K5"),
 'K6a': ('    momentum_profiles = sp.Tuple(vr, rho, xi)  # KNIFE_K6',
         "    momentum_profiles = sp.Tuple(sp.Symbol('V_r_constant'), rho, xi)  # KNIFE_K6"),
 'K6b': ('    momentum_profiles = sp.Tuple(vr, rho, xi)  # KNIFE_K6',
         "    momentum_profiles = sp.Tuple(vr, sp.Symbol('rho_momentum_constant'), xi)  # KNIFE_K6"),
 'K6c': ('    momentum_profiles = sp.Tuple(vr, rho, xi)  # KNIFE_K6',
         "    momentum_profiles = sp.Tuple(vr, rho, sp.Symbol('xi_w_constant'))  # KNIFE_K6"),
 'K7': ('    native_normal = native_cofactor / native_area  # KNIFE_K7',
        '    native_normal = vector((0, 0, 0, 1))  # KNIFE_K7'),
 'K8': ('    transport = vector(sum(derivative(flux[a, i], x[i]) for i in range(3))\n                       for a in range(4))  # KNIFE_K8',
        '    transport = sp.ImmutableMatrix(sp.zeros(4, 1))  # KNIFE_K8'),
 'K9': ('    paired_face_velocity = face_velocity  # KNIFE_K9',
        "    paired_face_velocity = vector(sp.Function('independent_power_' + str(a))(*x, t) for a in range(4))  # KNIFE_K9"),
 'K10': ("    energy_transport_channels = (op['J_E_live'],)  # KNIFE_K10",
         '    energy_transport_channels = ()  # KNIFE_K10'),
 'K11': ("    relaxation = op['P_ref_relax_live']  # KNIFE_K11",
         '    relaxation = sp.S.Zero  # KNIFE_K11'),
 'K12': ('    body_entries = ()  # KNIFE_K12',
         "    body_entries = ((-1, vector(action('SeparateBodyForce_' + str(a), state) for a in range(4))),)  # KNIFE_K12"),
 'K13': ("    carried_w = action('CarriedMomentumW', op['Pi_n'], jn, op['J_map'],\n                       op['N_br_live'], graph_velocity, geom, state, source,\n                       boundary, Str('outward_native_relative_mass_current'),\n                       Str('premise_3_local_material_velocity_at_each_transfer'))  # KNIFE_K13",
         '    carried_w = jn * graph_velocity[3]  # KNIFE_K13'),
}


def scalar_difference(baseline, corrupted):
    """Exact algebraic residual, with unresolved actions as formal operands.

    Protect calculus/functional atoms during rational normalization; no OPEN
    constitutive identity is invented. Restore them before printing. This also
    avoids recursively expanding their large operand registers as polynomials.
    """
    if baseline == corrupted:
        return sp.S.Zero
    atoms = set()
    def collect(value):
        if isinstance(value, (sp.Function, sp.Derivative, sp.Subs, sp.Lambda)):
            atoms.add(value)
        else:
            for arg in value.args:
                collect(arg)
    collect(baseline)
    collect(corrupted)
    protect = {atom: sp.Dummy('o2_difference_atom') for atom in sorted(atoms, key=sp.default_sort_key)}
    residual = corrupted.xreplace(protect) - baseline.xreplace(protect)
    residual = sp.cancel(sp.expand(residual))
    return residual.xreplace({value: key for key, value in protect.items()})


def difference(baseline, corrupted):
    """Normalize every scalar entry; unchanged matrices/containers emit zero.

    A nonzero formal residual is not a verdict about an OPEN constitutive law.
    Metadata uses a structural edit pair rather than subtraction of text.
    """
    if baseline == corrupted:
        return sp.S.Zero
    if isinstance(baseline, sp.MatrixBase) and isinstance(corrupted, sp.MatrixBase):
        if baseline.shape == corrupted.shape:
            entries = [difference(a,b) for a,b in zip(baseline,corrupted)]
            return (sp.S.Zero if all(item == 0 for item in entries) else
                    sp.ImmutableMatrix(baseline.rows, baseline.cols, entries))
    if isinstance(baseline, sp.Tuple) and isinstance(corrupted, sp.Tuple):
        if len(baseline) == len(corrupted):
            entries = [difference(a,b) for a,b in zip(baseline,corrupted)]
            return sp.S.Zero if all(item == 0 for item in entries) else sp.Tuple(*entries)
    if isinstance(baseline, sp.Equality) and isinstance(corrupted, sp.Equality):
        return difference(sp.Tuple(baseline.lhs, baseline.rhs),
                          sp.Tuple(corrupted.lhs, corrupted.rhs))
    if isinstance(baseline, sp.Expr) and isinstance(corrupted, sp.Expr):
        return scalar_difference(baseline, corrupted)
    return sp.Tuple(Str('STRUCTURAL_EDIT'), baseline, corrupted)


def repair_evidence(objects, engine):
    """Computed observations for the six repair findings; no expected values."""
    native_bound = set(sp.symbols('o2_s_face o2_q_face_1 o2_q_face_2 o2_q_face_3'))
    # Compare names because native chart coordinates carry real assumptions.
    native_names = {atom.name for atom in native_bound}
    assembled = ('MECHANICAL_LOAD', 'HOLD_INPLANE', 'HOLD_W',
                 'HOLD_GRAPH_NORMAL', 'ENERGY_STEADY')
    free_native = sp.Tuple(*(sp.Tuple(Str(tag), engine.text_tuple(*sorted(
        atom.name for atom in objects[tag].free_symbols if atom.name in native_names)))
        for tag in assembled))
    material_tags = ('MOMENTUM_DENSITY', 'MOMENTUM_FLUX', 'MOMENTUM_STORAGE',
                     'MOMENTUM_TRANSPORT', 'INTERNAL_FORCE', 'ENERGY_STORAGE',
                     'ENERGY_TRANSPORT', 'ENERGY_POWER', 'ENERGY_STEADY')
    operands = sp.symbols('R_br M_perp')
    dependencies = sp.Tuple(*(sp.Tuple(Str(tag), sp.Tuple(*(
        sp.sympify(objects[tag].has(operand)) for operand in operands))) for tag in material_tags))
    momentum_tags = ('MOMENTUM_DENSITY', 'MOMENTUM_FLUX', 'MOMENTUM_STORAGE',
                     'MOMENTUM_TRANSPORT', 'MATERIAL_INPUT_DIFFERENTIALS')
    xi_derivatives = sp.Tuple(*(sp.Tuple(Str(tag), sp.Tuple(*sorted(
        (item for item in objects[tag].atoms(sp.Derivative)
         if any(fn.func.__name__ == 'xi_w' for fn in item.expr.atoms(sp.Function))),
        key=sp.default_sort_key))) for tag in momentum_tags))
    return engine.record(free_native_names=free_native,
                         material_O7_O1=dependencies,
                         momentum_xi_derivatives=xi_derivatives)


def parse_stream(path, engine):
    values = {}
    namespace = {'__builtins__': {}, **vars(sp), 'Str': Str,
                 'ExprCondPair': ExprCondPair, **engine._RELATIONALS}
    seen = set()
    with path.open() as stream:
        for line in stream:
            tag, payload = line.rstrip('\n').split(': ', 1)
            if tag in seen:
                raise ValueError('duplicate tag ' + tag)
            seen.add(tag)
            if tag.startswith('PY_O2_') and tag[6:] in engine.SECTION9:
                values[tag[6:]] = eval(payload, namespace)
    if set(values) != set(engine.SECTION9):
        raise ValueError('section 9 stream coverage ' + str(path))
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch', required=True, type=Path)
    args = parser.parse_args()
    args.scratch.resolve().relative_to(LIVE.parents[3] / '_scratch')
    args.scratch.mkdir(parents=True, exist_ok=False)
    source = LIVE.read_text()
    spec = importlib.util.spec_from_file_location('o2_live', LIVE)
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    fold, _ = engine.load_model(*(str(p) for p in engine.FOLD_PATHS))
    baseline_path = args.scratch / 'baseline.stdout'
    with baseline_path.open('w') as stream, contextlib.redirect_stdout(stream):
        engine.run(fold=fold)
    baseline = parse_stream(baseline_path, engine)
    engine.emit('REPAIR_BASELINE', repair_evidence(baseline, engine), local=True)
    (args.scratch / 'baseline.py').write_text(source)
    engine.emit('ABLATION_SOURCE', sp.Tuple(Str(str(LIVE)),
                 Str(hashlib.sha256(source.encode()).hexdigest())), local=True)
    engine.emit('ABLATION_KNIVES', engine.text_tuple(*KNIVES), local=True)
    for knife, (old, new) in KNIVES.items():
        matches = source.count(old)
        engine.emit('ABLATION_' + knife + '_SITE', sp.Tuple(Str(old), Str(new),
                    sp.Integer(matches)), local=True)
        if matches == 0:
            engine.emit('ABLATION_' + knife + '_CONSTRUCTION',
                        engine.text_tuple('CONSTRUCTION_NOT_PRESENT', old), local=True)
            continue
        if matches != 1:
            raise ValueError('ambiguous knife site ' + knife)
        mutated = source.replace(old, new, 1)
        copy_path = args.scratch / (knife + '.py')
        copy_path.write_text(mutated)
        # Same file-location context as the production module; only the named
        # source statement differs. No import/export paths or emissions change.
        module = {'__file__': str(LIVE), '__name__': 'o2_copy_' + knife}
        exec(compile(mutated, str(copy_path), 'exec'), module)
        capture = args.scratch / (knife + '.stdout')
        with capture.open('w') as stream, contextlib.redirect_stdout(stream):
            module['run'](publish_delta=False, fold=fold)
        corrupted = parse_stream(capture, engine)
        if knife == 'K6c':
            engine.emit('REPAIR_K6C', repair_evidence(corrupted, engine), local=True)
        for name in engine.SECTION9:
            engine.emit('ABLATION_' + knife + '_' + name,
                        engine.record(baseline=baseline[name], corrupted=corrupted[name],
                                      difference=difference(baseline[name], corrupted[name])),
                        local=True)
    engine.emit('ABLATION_LIVE_SOURCE', sp.sympify(LIVE.read_text() == source), local=True)


if __name__ == '__main__':
    main()
