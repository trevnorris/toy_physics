#!/usr/bin/env python3
"""Self-contained reproduction of the S11c-c1 cross-engine reconcile bridges.

Reads ONLY the two committed engine transcripts (git-annex; `datalad get` first) — no scratchpad
dependency. Extracts the joined operands the comparator printed and applies the physically-justified
representational identifications (reconcile record §2), then reports the collapse. ⛔ The identifications
are read from the engine SOURCES, not tuned to force zero; the adversarial corruptions (run with --adversarial)
show the collapse is load-bearing.

Usage:
  python3 S11c_c1_reconcile_reproduce.py \
    --py  research/pde_ledger_v3/scripts/out/S11c_c1_bulk_closure_sympy_audit.out \
    --wl  research/pde_ledger_v3/mathematica/out/S11c_c1_bulk_closure_mathematica_audit.out
Requires a comparator run first is NOT needed — it re-runs the comparator family-scoped to get the operands,
OR (default) parses a pre-saved run log via --runlog. To regenerate the run log:
  python3 research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py \
     --family DTN_KERNEL --family FACE_RESPONSE_COEFFS > /tmp/reconcile_run.out
  python3 S11c_c1_reconcile_reproduce.py --runlog /tmp/reconcile_run.out
"""
import argparse, random, re, subprocess, sys, tempfile, os
import sympy as sp
from sympy import Add, Mul, Pow, Integer, Rational, Symbol, I, Function, Tuple

HFT = Function('HeldInactiveFourierTransform'); qOut = Function('qOut')
w1Profile = Function('w1Profile'); rhoField = Function('rhoBrBgRho4Constant')
for _n in ('chiOne', 'chiTwo', 'chiThree'):
    globals()[_n] = Function(_n)
NS = dict(Add=Add, Mul=Mul, Pow=Pow, Integer=Integer, Rational=Rational, Symbol=Symbol, I=I,
          Function=Function, Tuple=Tuple, DiracDelta=sp.DiracDelta, HeldInactiveFourierTransform=HFT,
          qOut=qOut, w1Profile=w1Profile, rhoBrBgRho4Constant=rhoField,
          chiOne=chiOne, chiTwo=chiTwo, chiThree=chiThree, oo=sp.oo, conjugate=sp.conjugate)

kO = [Symbol('kOne'), Symbol('kTwo'), Symbol('kThree')]
kI = [Symbol('kPrimeOne'), Symbol('kPrimeTwo'), Symbol('kPrimeThree')]
QOUT, QIN, HAT, RHO, MUTH, EPS = (Symbol('QOUT'), Symbol('QIN'), Symbol('HAT'),
                                  Symbol('RHO'), Symbol('MUTH'), Symbol('epsilon_shape'))


def plon(e):
    return e.xreplace({s: Symbol(s.name) for s in e.free_symbols})


def canon(expr, jet_sign=+1, freeze_two_legs=False):
    """Map opaque<->live under the reconcile identifications. jet_sign/freeze are the adversarial knobs."""
    e = plon(expr)
    if freeze_two_legs:  # adversarial: collapse WL's two legs into one (breaks the two-momentum structure)
        e = e.replace(lambda x: x.func == qOut, lambda x: Symbol('QSAME'))
    else:
        e = e.replace(lambda x: x.func == qOut and tuple(x.args[1]) == tuple(kO), lambda x: QOUT)
        e = e.replace(lambda x: x.func == qOut and tuple(x.args[1]) == tuple(kI), lambda x: QIN)
    e = e.replace(lambda x: x.func == HFT, lambda x: HAT)
    e = e.replace(lambda x: x.func == rhoField, lambda x: RHO)
    r = {Symbol('s11cc1_q_out_output'): (Symbol('QSAME') if freeze_two_legs else QOUT),
         Symbol('s11cc1_q_out_input'): (Symbol('QSAME') if freeze_two_legs else QIN),
         Symbol('s11cc1_w1_profile_hat_transfer'): HAT,
         Symbol('rho_br_bg_rho4_constant'): RHO}
    for i in range(3):
        r[Symbol(f's11cc1_k_input_{i+1}')] = kI[i]
        r[Symbol(f's11cc1_w1_profile_jet_hat_{i+1}')] = jet_sign * I * Symbol('L_W') * (kO[i] - kI[i]) * HAT
    for s in e.free_symbols:
        if 'mu_theta' in s.name:
            r[s] = MUTH
    return e.xreplace(r)


def numeval(A, B, trials, diagonal=False, seed=0):
    random.seed(seed)
    worst = 0.0
    for _ in range(trials):
        R = lambda: Rational(random.randint(3, 40), random.randint(2, 15))
        env = {Symbol(k): R() for k in ['L_W', 'W_0', 'c_s0', 'rho_m', 'eta_bg', 'Lambda_A_0',
                'Lambda_V_0', 'Lambda_X_0', 'tau_A', 'tau_V', 'tau_X', 'epsilon_shape']}
        env[Symbol('sigma_W')] = env[Symbol('eta_bg')] * env[Symbol('W_0')] / env[Symbol('L_W')]
        env[Symbol('omega')] = Integer(80)
        kv = [R() / 6 for _ in range(3)]
        for i in range(3):
            env[kO[i]] = kv[i]
            env[kI[i]] = kv[i] if diagonal else R() / 7
        cs0, om = env[Symbol('c_s0')], env[Symbol('omega')]
        env[QOUT] = sp.sqrt(om**2 / cs0**2 - sum(env[s]**2 for s in kO))
        env[QIN] = sp.sqrt(om**2 / cs0**2 - sum(env[s]**2 for s in kI))
        env[Symbol('QSAME')] = env[QIN]
        env[HAT] = R() + I * R(); env[RHO] = R(); env[MUTH] = R() + I * R()
        miss = (A.free_symbols | B.free_symbols) - set(env)
        if miss:
            return None, sorted(str(m) for m in miss)
        a = complex(sp.N(A.xreplace(env))); b = complex(sp.N(B.xreplace(env)))
        worst = max(worst, abs(a - b) / (abs(a) + abs(b) + 1e-30))
    return worst, None


def extract(runlog, family, key_substr):
    lines = open(runlog).read().split("\n")
    for i, l in enumerate(lines):
        if l.startswith(f"CASE family={family} ") and key_substr in l:
            a = lines[i+1][12:] if lines[i+1].startswith("operand_A = ") else None
            b = lines[i+2][12:] if lines[i+2].startswith("operand_B = ") else None
            if a and b and "<MISSING>" not in a and "<MISSING>" not in b:
                return a, b
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runlog", help="a comparator run log (DTN_KERNEL + FACE_RESPONSE_COEFFS)")
    ap.add_argument("--py"); ap.add_argument("--wl")
    ap.add_argument("--comparator", default="research/pde_ledger_v3/scripts/S11c_c1_cross_engine_comparator.py")
    ap.add_argument("--adversarial", action="store_true")
    args = ap.parse_args()
    runlog = args.runlog
    if runlog is None:
        if not (args.py and args.wl):
            sys.exit("give --runlog OR --py and --wl (to regenerate the run)")
        fd, runlog = tempfile.mkstemp(suffix=".out"); os.close(fd)
        with open(runlog, "w") as f:
            subprocess.run([sys.executable, args.comparator, "--py", args.py, "--wl", args.wl,
                            "--family", "DTN_KERNEL", "--family", "FACE_RESPONSE_COEFFS"],
                           stdout=f, check=True)
    # Seal 3 (two-momentum kernel) — all 4 cases, off-diagonal
    print("== Seal 3: DTN_KERNEL (two-momentum), off-diagonal ==")
    for anc in ("LAB_HELD", "MATERIAL_ADVECTED"):
        for face in ("1", "-1"):
            key = f"(BRANCH={anc}, FACE={face}, LEAF=KERNEL_EXPRESSION)"
            a, b = extract(runlog, "DTN_KERNEL", key)
            if a is None:
                print(f"  {anc} FACE={face}: no joined case"); continue
            w, _ = numeval(canon(eval(a, NS)), canon(eval(b, NS)), 8)
            print(f"  {anc} FACE={face}: worst_rel={w:.1e} {'AGREE' if w < 1e-9 else 'NONZERO'}")
    # Seal 1 (response coefficients) — LAB_HELD off-diagonal; MATERIAL on-diagonal (flat-leg convention)
    print("== Seal 1: FACE_RESPONSE_COEFFS (ε·A vs B) ==")
    for anc, diag in (("LAB_HELD", False), ("MATERIAL_ADVECTED", True)):
        key = f"(BRANCH={anc}, FACE=-1, DENSITY=RHO4_CONSTANT, LEAF=PRESSURE.FIRST_SHAPE_KERNEL.MU_THETA_COEFFICIENT)"
        a, b = extract(runlog, "FACE_RESPONSE_COEFFS", key)
        if a is None:
            print(f"  {anc}: no joined case"); continue
        w, miss = numeval(EPS * canon(eval(a, NS)), canon(eval(b, NS)), 8, diagonal=diag)
        tag = "on-diagonal" if diag else "off-diagonal"
        print(f"  {anc} μθ-coeff ({tag}): worst_rel={w if w is not None else miss} "
              f"{'AGREE' if (w is not None and w < 1e-9) else 'NONZERO/incomplete'}")
    if args.adversarial:
        print("== Adversarial (must be NONZERO) ==")
        a, b = extract(runlog, "DTN_KERNEL", "(BRANCH=LAB_HELD, FACE=1, LEAF=KERNEL_EXPRESSION)")
        A, B = eval(a, NS), eval(b, NS)
        for label, kw in (("wrong jet sign", dict(jet_sign=-1)), ("freeze two legs", dict(freeze_two_legs=True))):
            w, _ = numeval(canon(A, **kw), canon(B), 8)
            print(f"  {label}: worst_rel={w:.1e} {'BITES (good)' if w > 1e-9 else 'VACUOUS (bad)'}")


if __name__ == "__main__":
    main()
