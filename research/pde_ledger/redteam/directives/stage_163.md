---
unit_id: 163
batch: IV.6
created_at: 2026-05-27T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 163

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes outside the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl:28-79`

**Issue:**

The Mathematica script currently mirrors the SymPy script line-for-line: identical variable order, identical residual choices, identical hand-written form of `deltaG = gMinus*(dlnGq - dlnGs + 1/2*dlnKs - 1/2*dlnKq)`, identical `dF = (D[fComp, g] /. g -> gMinus)*dg + (D[fComp, r] /. g -> gMinus)*dr` construction. Cross-engine agreement is therefore structurally guaranteed rather than an independent check. Two of the load-bearing identities must be re-derived in Mathematica by a route distinct from the SymPy choreography.

**Required change:**

Edit only the Mathematica script. Make the following structural changes; keep the existing assertions in place (they may stay and continue to pass) and add the additional independent-derivation assertions alongside them.

### Change 1: re-derive `gPrime` via the implicit-function theorem from `fComp = 0`

Currently (line 33): `gPrime = D[gMinus, r];`

After the existing `gPrime` definition, add an independent derivation using implicit differentiation of `fComp == 0`. Append after current line 33 a block like:

```mathematica
(* Independent route: implicit-function derivative of g_-(r) from F(g,r)=0 *)
gPrimeImplicit = -((D[fComp, r])/(D[fComp, g])) /. g -> gMinus;
expectZero["gPrime matches implicit-function route", gPrime - gPrimeImplicit];
```

This route does not use `D[gMinus, r]` — it differentiates `F` partially in `g` and `r` separately and forms the implicit-function ratio. If the SymPy `gp = sp.diff(gminus, r)` had a sign error, this independent check would catch it.

### Change 2: re-derive `delta_perp` microscopic form via a Series expansion route

Currently (lines 47-54) the script hand-writes `deltaR = r*(dlnLam - 1/2*dlnKs - 1/2*dlnKq)` and `deltaG = gMinus*(dlnGq - dlnGs + 1/2*dlnKs - 1/2*dlnKq)`, then computes `deltaPerpMicro = FullSimplify[deltaG - gPrime*deltaR]`.

Add an independent block that builds `delta r` and `delta g` from the parent expressions
```
r = lam/Sqrt[Ks*Kq],   g = gq*Sqrt[Ks]/(gs*Sqrt[Kq])
```
by the chain rule on `Log[...]`, NOT by hand:

```mathematica
(* Independent route: derive delta r, delta g via chain-rule on Log *)
Clear[Ks, Kq, lam, gs, gq, eps];
rExpr  = lam/Sqrt[Ks*Kq];
gExpr  = gq*Sqrt[Ks]/(gs*Sqrt[Kq]);
(* perturb each parameter as p -> p (1 + eps dlnP), expand to O(eps) *)
pertRule = {
  Ks -> Ks (1 + eps dlnKs),
  Kq -> Kq (1 + eps dlnKq),
  lam -> lam (1 + eps dlnLam),
  gs -> gs (1 + eps dlnGs),
  gq -> gq (1 + eps dlnGq)
};
deltaRSeries = Coefficient[Series[rExpr /. pertRule, {eps, 0, 1}] // Normal, eps];
deltaGSeries = Coefficient[Series[gExpr /. pertRule, {eps, 0, 1}] // Normal, eps];
(* Substitute the parent-ratio identifications r = lam/Sqrt[Ks Kq], g = gq Sqrt[Ks]/(gs Sqrt[Kq]) *)
deltaRSubst = FullSimplify[deltaRSeries] /. {lam -> r*Sqrt[Ks*Kq]};
deltaGSubst = FullSimplify[deltaGSeries] /. {gq -> gMinus*gs*Sqrt[Kq]/Sqrt[Ks]};
expectZero["delta r series matches hand form", deltaRSubst - r*(dlnLam - dlnKs/2 - dlnKq/2)];
expectZero["delta g series matches hand form", deltaGSubst - gMinus*(dlnGq - dlnGs + dlnKs/2 - dlnKq/2)];
(* Now build delta_perp via the series-derived deltas *)
deltaPerpSeries = FullSimplify[deltaGSubst - gPrime*deltaRSubst];
expectZero["microscopic delta_perp via series route",
  deltaPerpSeries -
    (gMinus*(dlnGq - dlnGs - dlnLam + dlnKs) + (dlnKs + dlnKq - 2*dlnLam)/(4*s))];
```

This route uses `Series[..., {eps, 0, 1}]` and `Coefficient` — Mathematica machinery that has no SymPy counterpart in the existing script. The substitutions back to `r` and `gMinus` are deterministic identifications of the parent ratios with their abstract symbols.

After these blocks, leave the rest of the script (lines 57-104 in the current file) unchanged. Do not delete the existing `deltaPerpMicro` check — both routes coexisting is the point.

### Restore $Assumptions properly

After Change 2, the `Clear` of `Ks, Kq, lam, gs, gq` must be reset before the script continues, since the original `$Assumptions = Element[{r, sigmaStar, ...}, Reals]` block at original line 58 does not include them. Append after Change 2:

```mathematica
Clear[Ks, Kq, lam, gs, gq, eps, pertRule, rExpr, gExpr, deltaRSeries, deltaGSeries, deltaRSubst, deltaGSubst, deltaPerpSeries];
```

so the later sections (deltaC, dPi) work in a clean scope.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 163` and confirm:

1. The script still exits 0.
2. The two new `PASS` lines appear: `gPrime matches implicit-function route` and `microscopic delta_perp via series route`.
3. The two intermediate `PASS` lines also appear: `delta r series matches hand form`, `delta g series matches hand form`.
4. The pre-existing `PASS` lines (delta F, delta R, microscopic delta_perp identity, delta C, delta Pi split) all still appear.

The captured `.txt` output file is regenerated by the same exec step.
