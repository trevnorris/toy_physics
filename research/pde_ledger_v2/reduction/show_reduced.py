#!/usr/bin/env python3
"""Readable view of the registry: every quantity, its defining equation, and that
equation fully reduced down to the things you must actually supply.

Answers two questions directly:
  1. which quantities must I define for a simulation to run?
  2. for the ones I don't, what do they reduce to?
"""
import sys, pathlib, yaml, sympy as sp

HERE = pathlib.Path(__file__).parent
OPS = {"Add": lambda *a: sp.Add(*a), "Sub": lambda a, b: a - b,
       "Mul": lambda *a: sp.Mul(*a), "Div": lambda a, b: a / b,
       "Pow": lambda a, b: a ** b, "Sqrt": lambda a: sp.sqrt(a),
       "Neg": lambda a: -a}


def build(node, sym):
    """prefix-v1 -> sympy."""
    if isinstance(node, (int, float)):
        return sp.Integer(node) if isinstance(node, int) else sp.Float(node)
    if isinstance(node, list):
        head = node[0]
        if head == "Q":
            return sym.setdefault(node[1], sp.Symbol(node[1].split(".")[-1], positive=True))
        if head in OPS:
            return OPS[head](*[build(c, sym) for c in node[1:]])
        raise ValueError(f"unknown operator {head!r}")
    raise ValueError(f"unexpected node {node!r}")


def main():
    quantities = yaml.safe_load((HERE / "quantities.yaml").read_text())["quantities"]
    relations = yaml.safe_load((HERE / "relations.yaml").read_text())["relations"]
    sym, defs, status = {}, {}, {}

    for rel in relations:
        out, res = rel.get("designated_output"), rel.get("residual")
        if not out or not res:
            continue
        expr = build(res, sym)                      # residual == 0
        target = sym.setdefault(out, sp.Symbol(out.split(".")[-1], positive=True))
        sol = sp.solve(sp.Eq(expr, 0), target, dict=True)
        if sol:
            defs[out] = sp.simplify(sol[0][target])
            status[out] = rel["provenance_status"]

    def reduce_fully(expr, seen=None):
        """Substitute defined quantities until only undefined ones remain."""
        seen = seen or set()
        for _ in range(40):
            subs = {sym[q]: defs[q] for q in defs
                    if q not in seen and sym.get(q) is not None and sym[q] in expr.free_symbols}
            if not subs:
                break
            expr = sp.simplify(expr.subs(subs))
        return expr

    print("=" * 78)
    print("MUST BE SUPPLIED  (no defining equation -> a simulation input)")
    print("=" * 78)
    supplied = []
    for q in quantities:
        if q["qid"] not in defs:
            supplied.append(q)
            print(f"  {q['symbol']:<14} dim={q['dimension']['exponents']}  [{q['kind']}]")

    print()
    print("=" * 78)
    print("DERIVED  (has a defining equation -> computed, not supplied)")
    print("=" * 78)
    for q in quantities:
        qid = q["qid"]
        if qid in defs:
            print(f"\n  {q['symbol']}   [{status[qid]}]")
            print(f"      as written : {q['symbol']} = {defs[qid]}")
            red = reduce_fully(defs[qid])
            if red != defs[qid]:
                print(f"      reduced    : {q['symbol']} = {red}")
            leftovers = sorted(str(s) for s in red.free_symbols)
            print(f"      rests on   : {', '.join(leftovers) if leftovers else '(pure number)'}")

    print()
    print("=" * 78)
    print(f"SUMMARY: {len(supplied)} must be supplied, {len(defs)} are derived, "
          f"{len(quantities)} tracked total")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
