"""Orchestrator's independent checks for E (N6) and G (adjointness/directionality)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from S11c_c2_stdout_loader import value_of, to_dict  # committed loader (needs the .out; see loader's OUT path)
import sympy as sp

def blocks(tag):
    """Return {outer/inner: expr} for a nested block object."""
    d=to_dict(value_of(tag)); out={}
    if isinstance(d,dict):
        for o in d:
            if isinstance(d[o],dict):
                for i in d[o]: out['%s/%s'%(o,i)]=d[o][i]
            else: out[o]=d[o]
    return out

# ---------- G: block directionality of the self-energy (full object, not uniform limit) ----------
print('==== G: SELF_ENERGY_INCREMENT block directionality (LAB_HELD_RHO4) ====')
inc=blocks('PY_S11CC2_SELF_ENERGY_INCREMENT_LAB_HELD_RHO4_CONSTANT')
for k in sorted(inc):
    e=inc[k]
    z=(sp.expand(e)==0)
    print('   %-32s zero=%s  bytes~%d'%(k, z, len(sp.srepr(e))))

print('==== G: CLOSED_COUPLING_KERNEL block directionality (LAB_HELD_RHO4) ====')
ker=blocks('PY_S11CC2_CLOSED_COUPLING_KERNEL_LAB_HELD_RHO4_CONSTANT')
for k in sorted(ker):
    e=ker[k]
    z=(sp.expand(e)==0)
    print('   %-32s zero=%s  bytes~%d'%(k, z, len(sp.srepr(e))))

# ---------- E: does the N6 rep-invariance residual carry ONLY sigma_W-sector terms? ----------
# If residual.subs(sigma_W->0) == 0, the entire remnant carries sigma_W => leading O(eps),O(eps.eta)
# rep-invariance HOLDS and the non-invariance is confined to the sigma_W drain sector.
print('==== E: REP_INVARIANCE_RESIDUAL sigma_W-sector confinement (RHO4) ====')
for tag in ['PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHO4_CONSTANT',
            'PY_S11CC2_REP_INVARIANCE_RESIDUAL_RHOBR_CONSTANT']:
    try:
        b=blocks(tag)
    except KeyError:
        print('   MISSING tag',tag); continue
    # find the sigma_W symbol among atoms
    allexpr=sp.Add(*[v for v in b.values()]) if b else sp.S.Zero
    sig=[a for a in allexpr.atoms(sp.Symbol) if a.name in ('sigma_W','sigmaW','sigma_w')]
    print('   ',tag,'blocks=%d sigma_W_symbol=%s'%(len(b), sig))
    for k in sorted(b):
        e=sp.expand(b[k])
        if e==0:
            print('      %-32s zero'%k); continue
        e0 = e.subs({s:0 for s in sig}) if sig else e
        e0 = sp.expand(e0)
        print('      %-32s full_zero=%s  sigmaW->0_zero=%s'%(k, e==0, e0==0))
