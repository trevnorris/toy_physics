#!/usr/bin/env python3
"""
Independent SymPy verification of the 2PN self+static test-mass reduction.

What it checks
--------------
1) Build the same two-body candidate used in the Wolfram harness.
2) Extract the linear-in-m_A worldline coefficient.
3) Reduce to the test-mass limit and convert to (epsU, epsV2).
4) Build the direct self+static one-body candidate from the Bernoulli worldline.
5) Verify the residual is exactly zero.

This is a diagnostic to separate algebra/physics from any WL-specific extraction issue.
"""

import sympy as sp

# Symbols
mA, mB, Mcentral = sp.symbols("mA mB Mcentral", positive=True)
Gconst, cLight, rAB = sp.symbols("Gconst cLight rAB", positive=True)
lambdaRho = sp.symbols("lambdaRho")
vA2, vB2, vAB, vAn, vBn = sp.symbols("vA2 vB2 vAB vAn vBn")
epsU, epsV2, epsPhi = sp.symbols("epsU epsV2 epsPhi", positive=True)
lam = sp.symbols("lam")

# Two-body candidate through 2PN
L1PNTwoBodyFrozen = (
    sp.Rational(1, 2) * mA * vA2
    + sp.Rational(1, 2) * mB * vB2
    + (mA * vA2**2 + mB * vB2**2) / (8 * cLight**2)
    + Gconst * mA * mB / rAB
    + (Gconst * mA * mB) / (cLight**2 * rAB)
    * (
        sp.Rational(3, 2) * (vA2 + vB2)
        - sp.Rational(7, 2) * vAB
        - sp.Rational(1, 2) * vAn * vBn
    )
    - (Gconst**2 * mA * mB * (mA + mB)) / (2 * cLight**2 * rAB**2)
)

L2PNSelfStaticCandidate = (
    (mA * vA2**3 + mB * vB2**3) / (16 * cLight**4)
    + (7 * Gconst * mA * mB) / (8 * cLight**4 * rAB) * (vA2**2 + vB2**2)
    + (6 * Gconst**2 * mA * mB) / (cLight**4 * rAB**2) * (mB * vA2 + mA * vB2)
    + (lambdaRho * Gconst**3 * mA * mB * (mA**2 + mB**2)) / (2 * cLight**4 * rAB**3)
)

L2BodyCandidateThrough2PN = sp.expand(L1PNTwoBodyFrozen + L2PNSelfStaticCandidate)

# Robust linear-in-mA extraction
L2BodyLinearInmA = sp.expand(sp.diff(L2BodyCandidateThrough2PN, mA).subs(mA, 0))

# Test-mass reduction
L2BodyCandidateScaledForTestMassRaw = sp.expand(
    (L2BodyLinearInmA / cLight**2).subs(
        {
            mB: Mcentral,
            vA2: cLight**2 * epsV2,
            vB2: 0,
            vAB: 0,
            vAn: 0,
            vBn: 0,
        }
    )
)

LtestMassSubRules = {
    (Gconst * Mcentral) / (cLight**2 * rAB): epsU,
    (Gconst**2 * Mcentral**2) / (cLight**4 * rAB**2): epsU**2,
    (Gconst**3 * Mcentral**3) / (cLight**6 * rAB**3): epsU**3,
}

LtestMassCandidate2PNScaledFromTwoBody = sp.expand(
    sp.simplify(L2BodyCandidateScaledForTestMassRaw.subs(LtestMassSubRules))
)

# Direct self-sector Bernoulli worldline + static test-mass block
LselfExactMinimalScaled = -(1 + epsPhi) * sp.sqrt(1 - epsV2 / (1 + 4 * epsPhi))
Lself2PNScaled = sp.expand(
    sp.series(
        LselfExactMinimalScaled.subs({epsPhi: lam * epsPhi, epsV2: lam * epsV2}),
        lam,
        0,
        4,
    )
    .removeO()
    .subs(lam, 1)
)
Lself2PNScaledU = sp.expand(Lself2PNScaled.subs(epsPhi, -epsU))
LstaticTestMass2PNScaled = sp.expand(-epsU**2 / 2 + (lambdaRho / 2) * epsU**3)
LtestMassCandidate2PNScaledDyn = sp.expand(
    Lself2PNScaledU + LstaticTestMass2PNScaled + 1
)

residual = sp.expand(
    sp.simplify(LtestMassCandidate2PNScaledFromTwoBody - LtestMassCandidate2PNScaledDyn)
)

print("Linear-in-mA two-body coefficient:")
print(L2BodyLinearInmA)
print()

print("Test-mass candidate from two-body reduction:")
print(LtestMassCandidate2PNScaledFromTwoBody)
print()

print("Direct self+static test-mass construction:")
print(LtestMassCandidate2PNScaledDyn)
print()

print("Residual:")
print(residual)
print()

assert residual == 0
print("PASS: SymPy confirms the two-body reduction matches the direct self+static test-mass construction exactly.")
