"""Compare the exact algebraic outputs from the two rim calculations."""

import argparse
from pathlib import Path

import sympy as sp


def payloads(path):
    rows = {}
    for line in path.read_text().splitlines():
        name, value = line.split(" = ", 1)
        rows[name.strip('"')] = value.strip('"')
    return rows


def algebra(text, symbols):
    clean = text.replace("{", "[").replace("}", "]").replace("^", "**")
    value = sp.sympify(clean, locals=symbols)
    return value if isinstance(value,list) else [value]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sympy_stdout",type=Path)
    parser.add_argument("wolfram_stdout",type=Path)
    args = parser.parse_args()
    py, wl = payloads(args.sympy_stdout), payloads(args.wolfram_stdout)
    a, mu, B, Zn, Zt, Kl, Ke, c, Y = sp.symbols("a mu B Zn Zt Kl Ke c Y")
    sy = {str(s):s for s in (a,mu,B,Zn,Zt,Kl,Ke,c,Y)}
    sy.update(dict(b=B,m=mu,zn=Zn,zt=Zt,kl=Kl,ke=Ke,contact=c,y=Y,Pi=sp.pi))
    pairs = [
        ("normal and tangential virtual-work coefficients","boundary normal and tangential loads"),
        ("integrated exterior load on any sphere","exterior resultant"),
        ("exterior stiffness for a carried full trace","full-follow stiffness"),
        ("exterior energy","on-shell volume energy"),
        ("incoming, exterior, external hold","incoming, exterior, external hold"),
        ("local force residual H+P-S","net local residual"),
    ]
    differences = []
    for left,right in pairs:
        v, w = algebra(py[left],sy), algebra(wl[right],sy)
        if len(v) != len(w):
            raise ValueError(f"Output dimensions differ for {left}")
        delta = [sp.factor(i-j) for i,j in zip(v,w)]
        print(f"{left}: cross-engine residual = {delta}",flush=True)
        differences.extend(delta)
    lam=sp.symbols("lambda",positive=True)
    poisson=lam/(2*(lam+mu))
    classical=24*sp.pi*a*mu*(1-poisson)/(5-6*poisson)
    derived, = algebra(py["exterior stiffness for a carried full trace"],sy)
    residual=sp.factor(derived.subs(B,lam+2*mu)-classical)
    print(f"full-follow comparison to adhering sphere = {residual}",flush=True)
    if any(differences) or residual != 0:
        raise ValueError("Nonzero residual in the shared readouts or adhering-sphere check")
