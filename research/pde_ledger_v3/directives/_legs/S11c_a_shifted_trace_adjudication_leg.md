# Independent physics derivation — S11c-a shifted trace of a bulk field at the moving interface

## Your task
Derive, FROM THE SPEC ALONE and with a CAS script, the first-order linearised trace of a generic bulk
field at the S11c-a moving interface, and answer two specific yes/no questions about which terms survive.
Do NOT read either engine's implementation; do NOT look for an "expected" answer. Derive it yourself.

## Source of truth (read this, nothing else about the answer)
`research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md`
Read especially: §1 (lines ~44-56, the bulk fields and the steady bulk-normal drain `v_bulk_normal_0`),
§2a (lines ~163-200, the background profiles: `W_bg(y) = W0[1 + η w1(ξ)]`, `σ_W = η W0/L_W`, and the
`(ε,η,σ_W)` multigrade with the truncation "first order in wave amplitude and first shape order in each
background bookkeeper"), §2d (lines ~258-275, the supplied background state with `V_s^0 = 0`, `J_s^0 = 0`,
`𝒜_s^0 = 0`), and §3c (lines ~373-383, the trace-linearisation law and the dynamic window).

## The setup (all from the spec; restate in your own symbols)
A bulk field is `f(x,w,t) = f0(x,w) + ε·δf(x,w,t)`, where `f0` is the steady background and `δf` the wave
perturbation (bookkeeper `ε`). The interface for face `s∈{+1,−1}` sits at the normal coordinate
`w = h_s(x,t)`, whose background position is `h_s^0(x) = s·W_bg(x)/2` with `W_bg(x) = W0[1 + η·w1(ξ)]`,
and whose perturbation is `δh_s(x,t)` (the wave face displacement, O(ε)). The spec's §3c trace law is
`δ[f(x,h_s)] = δf(x,h_s^0) + δh_s·∂_w f0(x,h_s^0)`, to be applied to every traced bulk field, and §3c
explicitly says "do not evaluate first at `w = sW0/2` and then discard the face shift."

## Method — REQUIRED (a prose derivation is discarded)
Write a symbolic script (sympy or wolframscript) that:
1. Builds `f(x, h_s(x,t), t)` as an explicit composition, with `h_s = s·W0(1+η·w1(x))/2 + ε·δh(x,t)`,
   `f0(x,w)` and `δf(x,w,t)` generic symbolic functions, and takes the trace by substituting `w → h_s`.
2. Expands to FIRST order in `ε` (wave amplitude) AND to FIRST order in `η` (background shape order),
   keeping the mixed grade `ε·η`; drops `ε²`, `η²`, and higher.
3. Prints the resulting linearised trace, term by term, with its `(ε,η)` grade.
Save the script and its LITERAL stdout to named absolute paths and report both. Do NOT assert a conclusion
the script did not print.

## The two questions to answer from YOUR script's output
Q1 — background normal-jet term: Does the linearised trace contain a term proportional to
     `δh_s · ∂_w f0` (the wave face-shift times the NORMAL DERIVATIVE OF THE BACKGROUND field)? Then judge
     its status FOR THE BULK VELOCITY specifically, given the spec's scope: the background face-normal bulk
     velocity/flux obey `V_s^0 = 0`, `J_s^0 = 0`, and the drain `v_bulk_normal_0` is described (§1) as an
     inherited rest-frame scope limit with "the convective bulk problem is not reopened." Under that scope,
     is `∂_w f0` for the bulk velocity nonzero-and-contributing, or zero-and-dropped? Cite the spec lines
     that decide it. (§3c line ~375 fires the law for a "nonzero background face value OR derivative" —
     weigh that against §1/§2d.)
Q2 — shape-order correction of the perturbation term: Does the perturbation term `δf(x, h_s^0)` acquire a
     first-shape-order (`η`) correction because `h_s^0 = s·W_bg/2` itself depends on `η` through `W_bg`
     (i.e. a term `∝ η·(∂h_s^0/∂η)·∂_w δf = η·(s·W0·w1/2)·∂_w δf`)? Does §3c's "do not evaluate first at
     `w=sW0/2`" mandate keeping it, or is it out of the requested truncation?

Report: your script paths + literal stdout; Q1 answer with the deciding spec lines; Q2 answer with the
deciding spec lines. Physics filter — report only what bears on which terms physically survive.
