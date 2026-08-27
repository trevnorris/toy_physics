# Measurements — Grok review of S11c_a_py_projection_current_fix_directive.md
Date: 2026-08-26

## Claim: `dw_delta_j_bulk` does not enter `projection_terms`; enters only the trace at L651

```bash
python3 - <<'PY'
import pathlib, re
src = pathlib.Path('research/pde_ledger_v3/scripts/S11c_a_interface_geometry_sympy_audit.py').read_text()
# projection_terms body
m = re.search(r'^def projection_terms\b.*?(?=\n\n_PROJECTION_CACHE)', src, re.M|re.S)
body = m.group(0)
print('dw_delta_j_bulk in projection_terms:', 'dw_delta_j_bulk' in body)
print('affine_bulk_perturbation in projection_terms:', 'affine_bulk_perturbation' in body)
for i,l in enumerate(src.splitlines(),1):
    if 'dw_delta_j_bulk' in l or 'affine_bulk_perturbation(j_bulk' in l:
        print(f'{i}: {l.rstrip()}')
PY
```

Literal stdout (abridged to load-bearing lines):
```
dw_delta_j_bulk in projection_terms: False
affine_bulk_perturbation in projection_terms: False
177: dw_delta_j_bulk = {
257: SYMBOL_DIMENSIONS.update(... dw_delta_j_bulk ...)
651: affine_bulk_perturbation(j_bulk[i], dw_delta_j_bulk[face][i], face)
```

`projection_terms` current lines (quoted):
- 1152: `current_divergence = sp.Add(*(grad_j_bulk[i][i] for i in range(3)))` — in-plane only
- 1157/1159: `j_bulk[i]` under window / window_gradient
- 1160/1168: `WINDOW_NORMAL_CURRENT` uses bare `j_bulk[3]`

## Claim: `affine_bulk_perturbation` is face-keyed; `projection_terms` has no face argument

```text
affine_bulk_perturbation(reference_value, reference_normal_jet, face)  # L578-585
  reference_height = face/2 * W0
  return reference_value + (w - reference_height) * reference_normal_jet

dw_delta_j_bulk = {s: tuple(...) for s in FACES}   # L177-180, face-keyed
j_bulk = tuple(...)                                 # L176, NOT face-keyed

def projection_terms(branch, dof, representative, *, dynamic, ablate_direction=None)
  # no face parameter; cases keyed (branch, dof, representative) at L1194
```

## Claim: existing projection is post-IBP; adding pre-IBP `Ω ∂_w j` on top of affine `j` in `WINDOW_NORMAL` double-counts

Probe (`O=1+2w+3w²`, `jw=j0+(w-h)dj` on `[w1,w2]`):
```
IBP residual pre - ibp = 0
double-count sum (affine j in -∫j∂wO + explicit ∫O∂wj) = nonzero polynomial in (w1,w2,j0,dj,h)
```

`window_normal = plus_bg - minus_bg` matches `∂_w Ω` for `Ω=O(w-W0/2, -w-W0/2)` (L1095-1144, L1160).

## Claim: §1b / §3c scoping (spec quotes)

§1b L73-81: `j=ρ_4D v_bulk`, `∂_tρ_4D+∇₄·j=0`; projection = integrate against `Ω` and IBP in `w`.
§3c L375-392: shifted-trace law + zero-background language for **traced** bulk fields at `h_s⁰`; does not set `∂_w δj=0` under the slab integral.

## Claim: parallel WL engine uses a bulk `w`-function, not face-affine jets

WL L442-448, L805-818: `currentWPerturbation @@ Append[coordinates,{normal,time}]` inside
`Inactive[Integrate][D[window,normalCoordinate] currentW, ...]`.
