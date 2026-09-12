#!/usr/bin/env python3
"""Dimensionless analytic-matrix controls for generalized-chain coverage."""
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import S11c_d_mixing_scattering_sympy_audit as engine


def run():
    sp=engine.sp
    x=sp.Symbol('syntheticNormalIncrement',real=True)
    ansatzes={'simple':sp.diag(x,1),'jordan':sp.Matrix([[x,1],[0,x]]),
              'mixed':sp.diag(x**2,x**3)}
    records=[]
    for name,matrix in ansatzes.items():
        determinant=sp.Poly(matrix.det(),x)
        valuation=min(powers[0] for powers,coefficient in determinant.terms() if coefficient!=0)
        for order in ((2,3) if name=='mixed' else (3,)):
            jets=[matrix.diff(x,j).subs(x,0)/sp.factorial(j) for j in range(order+1)]
            data=engine.NormalTaylorChains.construct(jets,sp.QQ_I)
            records.append({'ansatz':name,'matrix':str(matrix),'maximumTaylorOrder':order,
                'determinant':str(determinant.as_expr()),'determinantValuation':valuation,
                'kernelCounts':data['KERNEL_COUNTS'],'kernelIncrements':data['KERNEL_INCREMENTS'],
                'stabilized':data['STABILIZED'],'lengths':[v['LENGTH'] for v in data['CHAINS']] if data['STABILIZED'] else 'UNRESOLVED',
                'leadingBasisRank':data.get('LEADING_BASIS_RANK','UNRESOLVED'),
                'determinantMultiplicityResidual':valuation-data['CHAIN_LENGTH_SUM'] if data['STABILIZED'] else 'UNRESOLVED',
                'chainCountResidual':data.get('CHAIN_COUNT_RESIDUAL','UNRESOLVED'),
                'blockKernelResiduals':[str(v) for v in data['BLOCK_KERNEL_RESIDUALS']],
                'chainResiduals':[[str(v) for v in c['EQUATION_RESIDUALS']] for c in data['CHAINS']] if data['STABILIZED'] else 'UNRESOLVED'})
    print(json.dumps({'instrumentSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'engineSha256':hashlib.sha256(Path(engine.__file__).read_bytes()).hexdigest(),
        'dimensionLTM':[0,0,0],'multigrade':[[0,0,0]],'epsilonLambdaSupport':[[0,0]],
        'syntheticMatrixControls':records},indent=2))


if __name__=='__main__':run()
