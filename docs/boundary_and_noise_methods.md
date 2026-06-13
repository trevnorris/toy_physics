# Boundary-Absorption & Noise-Control Methods for the Stage-1 Solver — Literature Synthesis

**Status:** reference note (methodology only; target-blind). Compiled 2026-06-12 from a three-agent
literature sweep (boundary-treatment method catalog; noise failure-mode catalog; closest published
analogues). Informs the §J.3 boundary-control gate (step 5) and the WP3-tangent wave-reflection work
(step 8). **No §H target / §F extraction content** — this is purely about *how* to truncate cleanly
and suppress numerical noise.

**Motivating problem.** The brief (§5, §9) records that the prior branch-realization attempt came out
"stable and open but **noise-dominated**" and could not be certified. Repeated fluid/superfluid
attempts in this program have failed the same way. This note collects the published state of the art
so we don't reinvent boundary absorption, and so we attack the *actual* dominant noise sources rather
than the assumed one.

---

## 0. Decision of record (2026-06-12, user)

**Defer the complex absorbing potential (CAP) / PML to step 8.** The Stage-1 WP1 branch is *stationary
+ isotropic* and therefore has **no propagating waves to absorb**, so step 5's job is to demonstrate
the finite-domain truncation does not *contaminate* the interior — a real soft-confinement "sponge"
(or simply far-field placement) suffices there. The genuine **outgoing-wave reflection coefficient**
is a property of the wave-like WP3 linearized tangent (step 8); the proper absorber (CAP / PML) is
designed and characterized **there**, not now. See §3.

---

## 1. The reframe: for a STATIONARY Newton solve, the boundary is not the #1 noise source

Almost all "absorbing boundary" literature is written for **time-domain** propagation, where a
reflected wave returns and corrupts the interior over later timesteps. Our solve is a **stationary
Newton–Krylov** solve — no time axis — so the mechanism by which a boundary corrupts the interior is
different (it enters through the elliptic/algebraic coupling of the discrete operator, as a
*converged-but-wrong* field, not as a returning wave). The ranked culprits for our class of solve,
with how the Stage-1 solver currently fares against each:

| # | Noise source | Literature remedy | Stage-1 status (as of step 4 / step-3b) |
|---|---|---|---|
| 1 | **Under-converged inner GMRES** (stiff Jacobian + gauge near-null space → loose inner tol leaves residual structure that *looks like* signal) | physics-based / Laplacian–Thomas–Fermi preconditioner; deflate gauge mode; Eisenstat–Walker forcing | **Not our problem.** We converge to ~1e-8–1e-12, far below discretization (~1e-3). Fixed *tight* GMRES tol (no Eisenstat–Walker) — slightly inefficient, but *safer* re: noise. step-3b preconditioner is the recommended lever. |
| 2 | **Gauge / constraint drift** in the EM sector (gauge null space → spurious field content) | gauge fixing; elliptic/projection divergence cleaning; deflation | **Largely de-risked.** H=Z gauge-fixing term `(1/ξ)∇(Z∇·A)` + A₀ Dirichlet anchor remove the gauge null-space. (Worth a Jacobian-spectrum spot-check at step 7.) |
| 3 | **Boundary closure** (reflecting / too-near truncation) | absorbing layer (sponge/CAP/PML), far-field placement, higher-order ABC | The step-5 question. For the *stationary* branch = contamination, not reflection (step 5 measures it). Wave reflection deferred to step 8. |
| 4 | **Spurious / near-zero discrete modes** (non-structure-preserving stencil) | structure-preserving / compatible FV; OSS-SUPG stabilization vanishing on stationary states | Conservative telescoping FV with exact r² measure (good); the compatible-diffusion point (Barsukow et al. 2025) is worth knowing. |
| 5 | **Aliasing of the stiff `|ψ|⁸ψ` term** near the throat | de-aliasing (padding rule); conservative flux form; refine | FV (not spectral) → classical aliasing less acute; step-4 order-2 convergence indicates current scales are resolved. |
| 6 | **Continuation near a fold** | pseudo-arclength + deflated continuation; fold detection | Natural-parameter continuation-in-K (not yet pseudo-arclength). Upgrade if a fold is hit. |

**Conclusion:** our stationary + tightly-converged + conservative-FV + gauge-fixed + preconditioned
construction *structurally avoids* the top suspects that sink time-domain attempts. The historical
"noise always gets in" is most consistent with time-domain reflection and/or under-convergence —
neither of which is our current regime.

---

## 2. Certification backbone (what makes "signal above noise" meaningful)

Both surveys converge on the same recipe — and it is already our step-4/5/7 sequence:

1. **Grid-convergence study** (Richardson / GCI) — order of accuracy + discretization error (step 4).
2. **Boundary-distance study** — vary the truncation radius; interior signal must be invariant as the
   boundary moves out (step 5).
3. **Measured (not modeled) noise floor** — drive the inner Newton/GMRES residual *below* the
   discretization error so iterative error cannot masquerade as signal; read the round-off-limited
   floor at the finest grid (step 7 error budget). *Caution (Roy review): estimating round-off by
   adding white noise is unreliable — measure it empirically.*
4. **Live constraint-functional monitor** — adopt the draining-throat precedent (Federici et al.
   2005): monitor a constraint/identity functional over the whole domain (they held it < 1e-6) as a
   running noise certificate alongside the primary residual.

---

## 3. Boundary-treatment families (for step 8's absorber, and the stationary-branch choice)

Constraints that decide it for us: **stationary** (not time-marching), **stiff polytropic** (non-cubic)
nonlinearity, **gauge-covariant + coupled Maxwell**, **conservative FV** on an r²-measure radial × axial
mesh. These rule out the elegant-but-inapplicable options:

- **Exact transparent BC / discrete TBC / DtN** — exact & nonreflecting but **nonlocal in time** (and
  in space in >1D); exact nonlinear TBC exists **only for the cubic/integrable NLS**. Not applicable
  to a stationary polytropic GPE.
- **Real "sponge" / Rayleigh damping `+σ(x)·field`** — composes trivially with Newton (static,
  additive, nonlinearity-agnostic, per-field) — **but a real potential is not a true absorber for a
  Schrödinger/GPE field; it reflects.** Fine as *soft confinement / interior protection* for a
  no-wave stationary branch; wrong tool for outgoing waves. **This is what step 5 currently uses
  (off by default), correctly scoped.**
- **Complex absorbing potential (CAP) `−iσ(x)`** — the right additive absorber; static, composes with
  Newton, per-field. Use a **transmission-free profile**: Manolopoulos (2002) — single width
  parameter (set to the de Broglie wavelength at min absorbed energy), provable reflection
  <1% at E_min, <0.01% for E≥3E_min. Equivalent optimum = smooth-exterior-scaling-derived CAP
  (Shemer–Brisker–Moiseyev 2005). **Recommended absorber for step 8.** Caveat: the wavelength sizing
  is heuristic for a stiff nonlinear field — tune empirically.
- **PML for matter waves** — impedance-matched (reflectionless at all angles in principle), demonstrated
  for the GPE (Antoine–Geuzaine–Tang 2020) and coupled NLS systems (Zheng 2007) — but the mature,
  well-posed Schrödinger PMLs are **time-dependent** (auxiliary ODEs; "quiescent-state" instability);
  a *stationary nonlinear-GPE* PML is essentially unproven. Higher payoff, higher risk than CAP.
- **Exterior complex scaling (ECS)** — the literature's "perfect absorber" (Scrinzi 2010, 14-digit) and
  natively *stationary* (resonance eigenproblem heritage) — but reflection-free only on discretizations
  that resolve the discontinuity at the scaling radius (FE/DVR), **not a plain conservative-FV stencil
  across R₀**; nonlinear/gauge-coupled ECS is unproven. Longer-horizon option.
- **Higher-order / self-adaptive local ABC** — our **Robin exit IS the lowest-order local ABC** (order-0/1
  truncation of the exact DtN; matches one wavenumber at normal incidence). Its single-wavenumber
  limitation is a credible reflection mechanism for a stiff nonlinear field's spread of effective
  wavenumbers → upgrade path is a self-adaptive local ABC (Trofimov et al. 2023) tracking the local k.

**Step-8 plan:** transmission-free CAP (Manolopoulos) paired with the Robin exit; escalate to a
self-adaptive higher-order ABC or a (stationary-adapted) matter-wave PML if residual low-momentum
reflection persists.

---

## 4. Closest published analogues (no single template exists)

No published work combines a gauge-covariant stiff GPE + localized Maxwell + (r,w) FV + stationary
Newton-through-an-open-throat. The methodology is assembled from three precedents:

- **Coupled gauge-field BVP + boundary handling** — gauged Q-balls / charged boson stars: Loiko &
  Shnir, arXiv:2207.02646; "Mapping Gauged Q-Balls," arXiv:2103.06905 (direct-BVP-over-shooting via a
  **compactified radial coordinate** + Neumann-at-origin / vacuum-Dirichlet-at-infinity); spinning
  gauged Q-balls on 2D (r,θ) via Newton (FIDISOL/CADSOL).
- **Stationary GPE solver core** — ATUS-PRO (arXiv:1506.07710): constrained Newton + FEM,
  cylinder-symmetry (r,z) reduction, μ-continuation, residual-of-stationary-GPE convergence;
  deflated continuation for localized GPE branches (Charalampidis/Farrell/Kevrekidis,
  arXiv:1612.08145 / 2004.10446); JFNK + preconditioning (Knoll & Keyes, JCP 193:357, 2004).
- **Open outgoing exit + noise certification** — draining-bathtub (Federici et al., gr-qc/0503089):
  horizon-excision / ingoing-radiation exit + **constraint-violation monitoring < 1e-6** as the
  quality gate. Cautionary realized case: the 2024 "giant quantum vortex" (Nature 630, 451) — a
  **hard outer wall produces standing-wave / bound-state contamination** (argues against reflecting
  truncation).

---

## 5. Curated references (verify the pre-2010 classics against originals before citing in a paper)

**Must-read reviews:**
- Antoine, Arnold, Besse, Ehrhardt, Schädle (2008), *CiCP* 4(4):729–796 — TBC/ABC taxonomy; the
  time(/space)-nonlocality tradeoff; the cubic-NLS nonlinear case.
- Antoine, Bao, Besse (2013), *Comput. Phys. Commun.* 184:2621 — GPE computational-methods review.
- Knoll & Keyes (2004), *J. Comput. Phys.* 193:357 — JFNK survey; preconditioning is the lever.

**Absorbers / boundaries:**
- Manolopoulos (2002), *J. Chem. Phys.* 117:9552 — transmission-free CAP (single param, provable bounds). **← step-8 default.**
- Shemer, Brisker, Moiseyev (2005), *Phys. Rev. A* 71:032716 — optimal reflection-free CAP = discretized smooth ECS.
- Klein, Antoine, Besse, Ehrhardt (2011), *CiCP* 10(5):1280 — ABCs for **stationary** N-D Schrödinger with nonlinearities (rare stationary+nonlinear match).
- Antoine, Geuzaine, Tang (2020), *CNSNS* 90:105406 — PML for the **GPE** (nonlinear + rotation).
- Zheng (2007), *J. Comput. Phys.* 227:537 — PML for NLS / coupled Schrödinger systems.
- Scrinzi (2010), *Phys. Rev. A* 81:053845 — infinite-range ECS perfect absorber (benchmark).
- Scrinzi, Stimming, Mauser (2014), *J. Comput. Phys.* 269:372 — PML ≠ ECS for Schrödinger.
- Modave, Delhez, Geuzaine (2014), *IJNME* 99(6):410 — discrete optimization of absorbing-layer profiles ("what σ₀ to pick").
- Trofimov et al. (2023), *Math. Methods Appl. Sci.* 46(15):16006 — self-adaptive local ABC tracking local k (Robin upgrade path).
- Farrell & Leonhardt (2004), arXiv:cond-mat/0404020 — PML for matter waves / BEC.

**Noise / conditioning / structure / certification:**
- Antoine, Levitt, Tang (2017), *J. Comput. Phys.* 343:92 — preconditioned stationary GPE (Laplacian/TF preconditioner).
- Gaul et al. (2013), arXiv:1303.5692 — deflation of near-null / gauge modes in Krylov.
- Munz et al. (2000), *JCP* 161:484 + Dedner et al. (2002), *JCP* 175:645 — divergence/constraint cleaning (use elliptic/projection form for a stationary solve).
- Barsukow, Ricchiuto, Torlo (2025), *NMPDE* 41:e23167 — structure-preserving FV; ordinary diffusion corrupts the set of discrete stationary states.
- Roy (2005) verification review + Eça & Hoekstra (2014) GCI — noise-floor / order-of-accuracy certification.

**Provenance / caveats:** compiled from arXiv/journal web search + a Wiley-weighted academic corpus
(which under-covered the atomic-physics CAP/ECS/PML-Schrödinger canon — those came from web/abstract
search). Verify the pre-2010 classics (Knoll–Keyes 2004; Munz 2000; Dedner 2002; Manolopoulos 2002)
against the primary sources before citing in any published paper.
