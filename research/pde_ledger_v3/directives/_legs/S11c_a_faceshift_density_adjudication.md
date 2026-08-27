# Independent from-spec adjudication — S11c-a T-e shifted density trace

## The object to adjudicate
S11c-a has two independent CAS engines that each build the same physics blind, in their own basis:
- PY (SymPy): `research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py`
- WL (Wolfram): `research/pde_ledger_v3/mathematica/S11c_a_interface_geometry_mathematica_audit.wl`

They emit a `FACE_SHIFT` family (spec object **T-e**, §430): the shifted-face trace operator obtained
from **§3c** for every bulk field consumed downstream. For the **bulk density** field, at
`(BRANCH=LAB_HELD, FACE=MINUS, DOF=DELTA_W)`, the two engines emit (verbatim, from the committed
transcripts):

```
PY :  -W_0*delta_rho_4D_face_minus_dw*eta_bg*w1_profile/2  +  delta_rho_4D_face_minus

WL :   W_0^2*e_W*eta_bg*(d^2/dw^2 rhoBulkBackground)*w1_profile/4
       - W_0*e_W*(d/dw rhoBulkBackground)/2
       - W_0*eta_bg*(d/dw rhoBulkPerturbation)*w1_profile/2
       + rhoBulkPerturbation(x1,x2,x3,(-W_0/2,time))
```

where in WL, `rhoBulkBackground` and `rhoBulkPerturbation` are abstract functions of
`(x1,x2,x3,{w,time})`; `e_W ≡ δW/W_0` is the reduced shape coordinate (spec §45). In the WL source,
`rhoBulkZero[coords_,normal_] := rhoBulkBackground @@ Append[coords,{normal,time}]` (line 435), and
`rhoBulkBackground` has **no definition** anywhere in the file (grep it). The PY side builds its density
background from `density_pair(...)` and evaluates the trace via `compose_physical_trace(...)`.

## What to decide (the whole question)
The two operands share the density **perturbation** content. They differ in that **WL additionally
carries terms proportional to `d/dw rhoBulkBackground` and `d^2/dw^2 rhoBulkBackground`** (the
background normal derivatives), which PY does not carry.

**From the spec alone, derive the correct §3c shifted trace of the bulk density and decide which engine
is faithful.** Specifically:
1. Per §3c's trace-linearisation law `δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s·∂_w f⁰(x,h_s⁰)`, what is the
   correct expression for the shifted density trace, and what is the value of the background normal
   derivative `∂_w ρ⁰` that appears in the shift term? Ground your answer in §2b (the two density
   representatives), §2d (the supplied background state `𝔅⁰`), and §3c's own sentence about what the
   supplied density background depends on. Quote the governing text.
2. Is WL's inclusion of a symbolic, undefined `rhoBulkBackground` whose `∂_w`/`∂²_w` are carried as
   live nonzero quantities consistent with §3c's requirement that *"every background face value or
   normal derivative appearing in this law is obtained by differentiating a member of the supplied
   background state 𝔅⁰; none may be introduced as a free premise"*? Or is PY's operand (which carries
   no such background term) the faithful one?
3. State plainly: does the WL density trace differ from PY only by a quantity the spec forces to a
   definite value, or is there a genuine physical term one engine is missing? If the former, what is
   that value and why?

## What you are handed
- The spec: `research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md` — read §1a, §1b, §2a–§2d,
  §3a–§3c, and the T-e line at §430. This is the source of truth. Form your own view of the correct
  §3c density trace from it BEFORE looking at either engine's operand above.
- The two engine source files named above (read the density-trace construction in each: PY
  `physical_trace_fields` / `compose_physical_trace`; WL `rhoBulkZero` / `traceSource` / the FACE_SHIFT
  emission).
- The committed transcripts (`scripts/out/S11c_a_interface_geometry_sympy_audit.out`,
  `mathematica/out/S11c_a_interface_geometry_mathematica_audit.out`) if you want the raw tags.

## Required method
- Derive the correct §3c density trace yourself from the supplied background state. If you assert a
  value for `∂_w ρ⁰`, show the member of `𝔅⁰` (§2b/§2d) you differentiated and the derivative — a prose
  claim ("I computed it and got zero") is discarded; write a runnable SymPy (or Wolfram) snippet and
  paste its literal stdout, and give its absolute path.
- Do NOT assume either engine is correct; the disagreement between them is the measurement. Do NOT try
  to make them agree — report what the spec forces.
- Physics filter: report a conclusion only if it bears on whether the emitted T-e density trace is
  physically right or §3c-compliant. "It would be wrong for a w-dependent background" is only relevant
  if the spec permits a w-dependent density background — check whether it does.
- If you find the spec is ambiguous about whether T-e should ground the background normal derivative or
  leave it symbolic for a downstream consumer, say so explicitly and quote the ambiguous text — a spec
  ambiguity is itself a finding.

## Return
A short verdict: (a) the correct §3c shifted density trace; (b) the value of `∂_w ρ⁰` with the `𝔅⁰`
member you differentiated; (c) which engine is §3c-faithful and why; (d) whether the difference is a
spec-forced value or a genuine missing/extra physical term; (e) any spec ambiguity. Cite line numbers
and paste your script's literal stdout.
