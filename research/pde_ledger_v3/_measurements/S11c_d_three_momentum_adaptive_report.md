# S11c-d independent adaptive outer quadrature

Source/profile refinement is published and annex-verified at 9ecb6ed8/16838693.
Its final raw source/profile changes are below 6.14e-19 / 3.64e-17. The next
calculation independently replaces only the outer Gauss rule with adaptive GK21;
inner momentum, source and profile rules remain at their accepted resolutions.
The physical engine and all prior constructors remain unchanged.

The new helper derives its conditional batch recursion and group evaluator from
the native method bodies. Exact AST reversal removes only the fixed-outer
recursion entry and the omitted measure in the conditional volume. Conditional
point units are the original row units minus one momentum unit. Four isolated
workers cover two fields and two disjoint outer half-intervals, with one native
thread each. Every conditional point and partial sum is saved before guards.

Preflight is prepared: compare production-order conditioned prefixes with saved
native operands, compare complete underresolved Gauss grouping and full actions,
then run an explicitly underresolved four-worker adaptive/emission smoke test.
No new adaptive physical result has yet been accepted. Production tolerance is
5e-11 absolute per half in the declared unit frame, zero relative tolerance;
the smoke test has its separate loose numerical tolerance and limited scope.
Physical domain/tail and Abel checks remain work before scattering and poles.
