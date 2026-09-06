"""Orchestrator's INDEPENDENT adjudication of section F (uniform limit).
Question: in the uniform limit, does the genuine closure-induced coupling
[= extract(close)|_uniform] decouple to zero, leaving only the -extract(open)
bare-slot bookkeeping? (fresh-Claude leg: YES; Grok: NO, a Z0.mu_theta Integral survives.)

Method (independent of the leg's 'Trial'-label heuristic): the increment =
extract(close) - extract(open). extract(close) substitutes ALL delta_p slots, so
it carries NO bare delta_p_{plus,minus}/d_w_delta_p symbols; every bare-delta_p term
therefore comes only from -extract(open). Zero those bare slots -> what remains is
extract(close)|_uniform. .doit() the Integrals and test == 0.
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from S11c_c2_stdout_loader import value_of, to_dict  # committed loader (needs the .out; see loader's OUT path)
import sympy as sp

BARE = {'delta_p_minus','delta_p_plus','d_w_delta_p_minus','d_w_delta_p_plus'}

def bare_zero_subs(e):
    sub={}
    # zero any Symbol or AppliedUndef whose name is a bare open slot
    for a in e.atoms(sp.Symbol):
        if a.name in BARE: sub[a]=0
    for a in e.atoms(sp.core.function.AppliedUndef):
        if a.func.__name__ in BARE: sub[a]=0
    # also symbols/functions carrying the bare name as prefix (jets like d_w_delta_p_plus)
    for a in e.atoms(sp.Symbol):
        if any(a.name==b for b in BARE): sub[a]=0
    return e.subs(sub) if sub else e

for tag in ['PY_S11CC2_UNIFORM_LIMIT_OPERAND_LAB_HELD_RHO4_CONSTANT',
            'PY_S11CC2_UNIFORM_LIMIT_OPERAND_MATERIAL_ADVECTED_RHO4_CONSTANT',
            'PY_S11CC2_UNIFORM_LIMIT_OPERAND_LAB_HELD_RHOBR_CONSTANT',
            'PY_S11CC2_UNIFORM_LIMIT_OPERAND_MATERIAL_ADVECTED_RHOBR_CONSTANT']:
    d=to_dict(value_of(tag)); print('====',tag)
    for outer in d:
        if not isinstance(d[outer],dict): continue
        for inner in d[outer]:
            e=sp.expand(d[outer][inner])
            if e==0:
                print('   [%s/%s] block ZERO'%(outer,inner)); continue
            # names present, and whether Integrals appear
            names=sorted({a.name for a in e.atoms(sp.Symbol) if a.name in BARE})
            n_int=len(e.atoms(sp.Integral))
            closed_part=bare_zero_subs(e)
            closed_part=sp.expand(closed_part)
            # try to evaluate any surviving integrals
            try:
                closed_doit=sp.expand(closed_part.doit())
            except Exception as ex:
                closed_doit=('doit-failed:%s'%ex)
            is_zero = (closed_part==0)
            is_zero_doit = (closed_doit==0) if not isinstance(closed_doit,str) else 'n/a'
            print('   [%s/%s] full_block_zero=%s bare_slots=%s n_Integral=%d  '
                  'closed_part_zero(expand)=%s  closed_part_zero(doit)=%s'
                  %(outer,inner, e==0, names, n_int, is_zero, is_zero_doit))
            if not is_zero:
                # show a short head of what survives, to characterize it
                surv = closed_part if not isinstance(closed_doit,str) else closed_part
                head = sp.srepr(surv)[:180]
                print('        SURVIVING(closed_part) head:', str(surv)[:200])
