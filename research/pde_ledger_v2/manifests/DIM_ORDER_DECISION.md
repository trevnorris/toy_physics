# Manifest dimension interchange — DECISION (2026-07-24)

## Retraction

The v1 statement that `(L,M,T)` was the audit scripts' common order and
therefore required zero transposition is retracted. The script corpus is
mixed. In particular, the four v2 pilots use `[L,T,M]`, while stage027 uses
`(L,M,T)`.

No positional tuple order is canonical across stage boundaries.

## Canonical interchange form

The v2 canonical form is the named dimension map:

```json
{"L":"-2","M":"1","T":"-2"}
```

Keys may be omitted only when their exponent is zero. Exponents are exact
rational strings. The named map is authoritative for interchange, C1, and C4.

Every symbol also records:

- `dim_source_order`: the positional order used by that stage's audit script;
  and
- `dim_source`: the script path and stable locus where the tuple/dimension is
  defined.

Extraction reads the source tuple using `dim_source_order` and transposes it to
the named map. C4 repeats that operation as an independent source certificate.

## Observed source orders

| Stage | Script source order | Evidence |
|---|---|---|
| 030 | `LTM` | SymPy `Dim` docstring and fields at lines 172–198: `[L,T,M]` |
| 031 | `LTM` | SymPy `Dim` docstring and fields at lines 140–166: `[L,T,M]` |
| 032 | `LTM` | SymPy dimension block: `dim_L=(1,0,0)`, `dim_E=(2,-2,1)` |
| 044 | `LTM` | SymPy `dimensional_tooth`, line 416: `[L,T,M]` |
| 027 | `LMT` | SymPy audit line 793 explicitly prints `(L,M,T)` |

Thus the required pilot table is 030/031/032/044 = `LTM`; stage027 is a
confirmed example that the wider corpus is mixed.

## Non-symmetric sentinel

Dimension-order tests must use distinct M and T exponents. The standard
sentinel is

`a_B = L^-2 M^1 T^-2`.

Its named form is `{"L":"-2","M":"1","T":"-2"}`. In an `LTM` source script,
the positional tuple is `(-2,-2,1)`. A checker that reads it as `LMT` produces
the wrong named map and must fail the source-certificate comparison.

The register header or any manifest-v1 tuple must not be used to infer a
stage's source order. Read the stage's own audit script and record its locus.
