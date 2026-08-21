# Measurements for S11b_wl_engine_review_disposition.md (rule 2)

## F-WL-1 — the correct `ZPERM_SLICE_MAP` coefficient sign, derived from the spec §4 affinity alone

```
$ python3 research/pde_ledger_v3/steps/_measurements/s11b_zperm_slicemap_sign.py
J on mu_s=0 slice            : -Lambda_A*dp/rho_m + Lambda_V*V
Lambda_p = coeff of dp in J  : -Lambda_A/rho_m
coeff of V in J on slice     : Lambda_V  (must equal Lambda_V)
Lambda_p + Lambda_A/rho_m    : 0  (0 => Lambda_p = -Lambda_A/rho_m)
Lambda_p - Lambda_A/rho_m    : -2*Lambda_A/rho_m  (0 => Lambda_p = +Lambda_A/rho_m)
```

Reads: `Lambda_p = coeff of dp in J on the mu_s=0 slice = -Lambda_A/rho_m`, velocity channel
unchanged (`coeff of V = Lambda_V`), and `Lambda_p + Lambda_A/rho_m = 0`. The WL engine emits
`+Lambda_A/rho_m` (residual-coefficient extraction, .wl L858-859 + equationZeroForm L367); SymPy
G13 (sympy_audit.py:1329-1334) and both SymPy legs give `-Lambda_A/rho_m`. F-WL-1 confirmed.

## F-WL-2 / F-WL-3 — tautological/decoupled checks: backed by the committed .wl code +
the two legs' FORM-ablation transcripts (byte-identical output under a traction-sign structural
change). PRESSURE_WORK_SIGN_CHECK: .wl L1168-1181 (slabPressurePower = -1/2 Re[X],
outgoingBulkPower = +1/2 Re[X], same X => residual identically 0). energyEquationRules hand-typed
at .wl L1116-1135 (not the assembled fullSystem).
