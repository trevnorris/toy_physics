"""Exact last resort for a mixed material prism and a nonnegative spatial Q1 field.

NLSAT decides a rational polynomial feasibility problem on the closed reference
prism. At a mesh knot only the basis functions with strictly positive weight
contribute; an endpoint basis includes its own exterior domain face.
The returned numbers are display approximations. The exact SMT model and its
algebraic signs are retained for witnesses on a tangency.
"""
import ctypes
from fractions import Fraction
from itertools import product
import re

import numpy as np

from mouth_active_prism_orientation import _coefficients


def _number(x):
    v = Fraction(x)
    p, q = abs(v.numerator), v.denominator
    term = str(p) if q == 1 else f"(/ {p} {q})"
    return term if v >= 0 else f"(- {term})"


def _define(name, body, sort="Real"):
    return f"(define-fun {name} () {sort} {body})"


def _combination(weights, terms):
    return "(+ " + " ".join(f"(* {w} {t})" for w, t in zip(weights, terms)) + ")"


class _Session:
    """Isolated libz3 SMT-LIB session. The distro runtime supplies libz3.so.4."""
    def __enter__(self):
        lib = ctypes.CDLL("libz3.so.4")
        ptr, string = ctypes.c_void_p, ctypes.c_char_p
        for name, args, result in (
            ("Z3_mk_config", [], ptr), ("Z3_del_config", [ptr], None),
            ("Z3_mk_context", [ptr], ptr), ("Z3_del_context", [ptr], None),
            ("Z3_set_error_handler", [ptr,ptr], None),
            ("Z3_eval_smtlib2_string", [ptr,string], string),
            ("Z3_get_error_code", [ptr], ctypes.c_int),
            ("Z3_get_error_msg", [ptr,ctypes.c_int], string),
        ):
            fn = getattr(lib,name)
            fn.argtypes, fn.restype = args, result
        config = lib.Z3_mk_config()
        self.ctx = lib.Z3_mk_context(config)
        lib.Z3_del_config(config)
        if not self.ctx:
            raise RuntimeError("Cannot initialize exact real-algebraic decision")
        lib.Z3_set_error_handler(self.ctx,None)
        self.lib = lib
        return self

    def __call__(self, program):
        result = self.lib.Z3_eval_smtlib2_string(self.ctx,program.encode("ascii"))
        error = self.lib.Z3_get_error_code(self.ctx)
        if error:
            raise ArithmeticError(self.lib.Z3_get_error_msg(self.ctx,error).decode())
        text = result.decode().strip()
        if "(error " in text:
            raise ArithmeticError(text)
        return text

    def __exit__(self, *unused):
        self.lib.Z3_del_context(self.ctx)


def _hats(support, points):
    axes = support.axes
    region = support.box(points)
    choices = [range(s.start, s.stop) for s in region]
    names = ("px","pt","pa")
    definitions, boxes, order_terms = [], [], []
    for node in product(*choices):
        value = support.values[node]
        if value == 0:
            continue
        tests, hat_names = [], []
        for d, (axis, center, coord) in enumerate(zip(axes,node,names)):
            mid = _number(axis[center])
            hat = f"h{d}_{center}"
            left = center > 0
            right = center < len(axis)-1
            halves = []
            if left:
                low = _number(axis[center-1])
                tests.append(f"(> {coord} {low})")
                halves.append(f"(/ (- {coord} {low}) (- {mid} {low}))")
            if right:
                high = _number(axis[center+1])
                tests.append(f"(< {coord} {high})")
                halves.append(f"(/ (- {high} {coord}) (- {high} {mid}))")
            # The validated prism lies inside the closed carrier domain. Outer
            # endpoints need only the one existing side of the hat.
            tent = f"(ite (<= {coord} {mid}) {halves[0]} {halves[1]})" if left and right else halves[0]
            item = _define(hat, f"(let ((h {tent})) (ite (> h 0) h 0))")
            if item not in definitions:
                definitions.append(item)
            hat_names.append(hat)
        boxes.append("(and " + " ".join(tests) + ")")
        order_terms.append(f"(* {_number(value)} {' '.join(hat_names)})")
    return definitions,boxes,order_terms


def _decimal_witness(text,names):
    # Only simple symbols are requested. Z3 prints a negative decimal as
    # (- a), including for algebraic results rounded for display.
    values={}
    for name in names:
        val=re.search(r"\("+name+r"\s+(\(-\s+[\d.]+\??\)|[\d.]+\??)\)",text)
        if not val:
            raise ArithmeticError(f"Unexpected value of {name} in {text}")
        s=val.group(1).replace("?","")
        values[name]=-float(s[2:-1]) if s.startswith("(-") else float(s)
    return values


def decide(points, scale, support):
    """Return (certified, None) or (active_fold, model), on this entire prism.

    points are exact grid-plus-offset rational positions [plane,corner,axis];
    scale is the positive reference parameter-cell volume. support has the
    _Support contract (absolute nodal values in [0,1] and tensor-product axes).
    A library/solver failure raises a numerical issue, never a physical fold.
    """
    if scale<=0 or np.shape(points)!=(2,3,3):
        raise ValueError("Expected a positive reference cell and six chart points")
    for d,axis in enumerate(support.axes):
        if min(points[:,:,d].flat)<axis[0] or max(points[:,:,d].flat)>axis[-1]:
            raise ValueError("Material prism exceeds the physical order chart")
    hats,boxes,terms = _hats(support,points)
    if not boxes:
        return "certified",None

    program = ["(set-logic QF_NRA)",
               "(declare-const b1 Real)", "(declare-const b2 Real)", "(declare-const y Real)",
               _define("b0","(- 1 b1 b2)"),_define("l","(- 1 y)"),
               "(assert (and (>= b0 0) (>= b1 0) (>= b2 0) (>= y 0) (<= y 1)))"]
    bary = ("b0","b1","b2")
    for d,coord in enumerate(("px","pt","pa")):
        terms_by_corner = [_combination(("l","y"),map(_number,points[:,j,d])) for j in range(3)]
        program.append(_define(coord,_combination(bary,terms_by_corner)))
    coefficients = _coefficients(points)
    weights = ("(* l l)","(* 2 l y)","(* y y)")
    det_by_corner = [_combination(weights,(_number(x/scale) for x in row)) for row in coefficients]
    program.extend([_define("detF",_combination(bary,det_by_corner)),
                    "(assert (<= detF 0))",
                    "(assert (or " + " ".join(boxes) + "))",
                    "(check-sat-using qfnra-nlsat)"])
    with _Session() as session:
        status = session("\n".join(program))
        if status == "unsat":
            return "certified",None
        if status != "sat":
            raise ArithmeticError(f"Undecided exact order-support query: {status}")
        names = ("b0","b1","b2","y","px","pt","pa","detF","ord")
        readback = "(get-value ("+" ".join(names)+"))"
        exact = session("\n".join(hats + [_define("ord","(+ "+" ".join(terms)+")"),
                                         "(set-option :pp.decimal false)",readback]))
        flags=session("(get-value ((> ord 0) (< detF 0) (= detF 0)))")
        truth=[v=="true" for v in re.findall(r"\b(true|false)\b",flags)]
        if len(truth)!=3 or not truth[0] or truth[1]==truth[2]:
            raise ArithmeticError(f"Exact witness fails physical signs: {flags}")
        raw=session("(set-option :pp.decimal true)\n(set-option :pp.decimal_precision 28)\n"+readback)
        v=_decimal_witness(raw,names)
        if v["ord"]==0. or (truth[1] and v["detF"]==0.):
            v=_decimal_witness(session("(set-option :pp.decimal_precision 350)\n"+readback),names)
    return "active_fold",dict(barycentric=[v[k] for k in bary],angular_fraction=v["y"],
                              placement_x_t_angle=[v[k] for k in ("px","pt","pa")],
                              chart_determinant=v["detF"],exact_determinant_sign=-int(truth[1]),
                              order=v["ord"],exact_order_positive=True,exact_model=exact)
