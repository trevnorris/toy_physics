#!/usr/bin/env python3
"""P1–P4 compact native identification. Run only under the host resource guard.
Selected original AST helpers on small exact fixtures; no production imports or operands.
"""
import ast
import hashlib
import json
import platform
from pathlib import Path
from types import SimpleNamespace
import numpy as np

BASE=Path(__file__).resolve().parents[1]
M=BASE/'_measurements'
PATHS={'current':M/'S11c_d_continuum_currents.py','response':M/'S11c_d_continuum_response.py',
       'spec':BASE/'directives/S11c_d_SHARED_PHYSICS.md'}
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    before={str(p.relative_to(BASE)):sha(p) for p in PATHS.values()}
    ast_hashes={}
    def selected(key,names):
        nodes=[n for n in ast.parse(PATHS[key].read_text()).body if isinstance(n,ast.FunctionDef) and n.name in names]
        assert {n.name for n in nodes}==set(names)
        for node in nodes:node.decorator_list=[]
        tree=ast.Module(body=nodes,type_ignores=[])
        ast_hashes[key+':'+','.join(names)]=hashlib.sha256(ast.dump(tree).encode()).hexdigest()
        return tree
    res={'np':np};exec(compile(selected('response',['adjoint','lambda_series']),'<selected response>', 'exec'),res)
    grades=((0,0),(1,0),(0,1),(1,1))
    cur={'np':np,'G':grades,'response':SimpleNamespace(adjoint=res['adjoint'])}
    exec(compile(selected('current',['multiply','quadratic','subtract','amplitude_bookkeeping','quotient']),'<selected current>','exec'),cur)
    checks=[]
    def check(name,actual,expected,kind='identity'):
        assert np.array_equal(actual,expected),name
        checks.append({'name':name,'kind':kind,'passed':True})
    def reject(name,actual,wrong):
        assert not np.array_equal(actual,wrong),'Undetected control: '+name
        checks.append({'name':name,'kind':'wrong-formula control','passed':True})
    arr=lambda x:np.array([[x]],complex)
    rectangle={g:arr(v) for g,v in zip(grades,[1,2,3,4])}
    parts=cur['amplitude_bookkeeping'](rectangle)
    for g in grades:
        check('component_reconstruction_'+str(g),parts['reconstructionResidual'][g],arr(0))
        check('delta_includes_zero_and_first_jets_'+str(g),parts['delta'][g],parts['zeroJetContrast'][g]+parts['firstJet'][g])
    path=res['lambda_series'](rectangle,2)
    check('real_path_a0_a1_a2',np.array([path[d][0,0] for d in range(3)]),np.array([1,8,8]))
    reject('wrong_mixed_ratio_power',path[2],arr(16))
    reject('delta_omits_zero_jet',parts['delta'][(1,0)],parts['firstJet'][(1,0)])
    # Full non-diagonal two-channel complex forms; columns are separate supplied incidents.
    amp={g:np.array(x,complex) for g,x in zip(grades,[ [[1,1j],[1j,2]],[[2,0],[1,1j]],[[0,1],[2j,1]],[[1,2],[0,1j]] ])}
    metric={g:np.array(x,complex) for g,x in zip(grades,[ [[2,1j],[-1j,3]],[[1,2],[2,-1]],[[0,1j],[-1j,2]],[[3,1],[1,2]] ])}
    full=cur['quadratic'](amp,metric);a=res['lambda_series'](amp,2);j=res['lambda_series'](metric,2);q=res['lambda_series'](full,2)
    independent={d:sum((a[i].conj().T@j[k]@a[l] for i in range(3) for k in range(3) for l in range(3) if i+k+l==d),np.zeros((2,2),complex)) for d in range(7)}
    for d in range(7):check('all_full_current_terms_degree_'+str(d),q[d],independent[d])
    t=.5
    eval_series=lambda s:sum((t**d*x for d,x in s.items()),np.zeros_like(s[0]))
    check('actual_polynomial_contraction',eval_series(q),eval_series(a).conj().T@eval_series(j)@eval_series(a))
    check('explicit_higher_remainder',eval_series(q)-sum((t**d*q[d] for d in range(3)),np.zeros_like(q[0])),sum((t**d*q[d] for d in range(3,7)),np.zeros_like(q[0])))
    reject('omit_offdiagonal_current',q[2],sum((a[i].conj().T@np.diag(np.diag(j[k]))@a[l] for i in range(3) for k in range(3) for l in range(3) if i+k+l==2),np.zeros((2,2),complex)))
    # Independent scalar convention witnesses corresponding to the formal controls.
    scalar_amp={g:arr(v) for g,v in zip(grades,[1,2,0,3])}
    scalar_metric={g:arr(v) for g,v in zip(grades,[1,4,0,5])}
    sq=res['lambda_series'](cur['quadratic'](scalar_amp,scalar_metric),1)
    check('scalar_degrees_0_to_6',np.array([sq[d][0,0] for d in range(7)]),np.array([1,8,31,72,107,96,45]))
    reject('omit_current_variation',sq[2],arr(10))
    reject('omit_a2',sq[2],arr(25))
    reject('omit_baseline_interference',sq[2],arr(9))
    reject('discard_higher_terms',sq[3],arr(0))
    sp=cur['amplitude_bookkeeping'](scalar_amp)
    induced=res['lambda_series'](cur['quadratic'](sp['delta'],scalar_metric),1)
    baseline=res['lambda_series'](cur['quadratic'](sp['baseline'],scalar_metric),1)
    check('induced_second_coefficient',induced[2],arr(4))
    check('total_minus_baseline_second',sq[2]-baseline[2],arr(26))
    reject('conflate_subtracted_and_induced',sq[2]-baseline[2],induced[2])
    no_baseline={**scalar_amp,(0,0):arr(0)}
    nb=res['lambda_series'](cur['quadratic'](no_baseline,scalar_metric),1)
    check('baseline_free_second',nb[2],arr(4))
    # Complex diagonal quotient is tested per incident column, not as a full matrix inverse.
    den={0:np.diag([2,2j]),1:np.diag([1,1j]),2:np.diag([2,-2j])}
    num={0:np.array([[2,7],[8,2j]],complex),1:np.array([[3,9],[10,3j]],complex),2:np.array([[5,11],[12,5j]],complex)}
    coeff,residual=cur['quotient'](num,den)
    expected=[np.array([1,1],complex),np.array([1,1],complex),np.array([1,3],complex)]
    for d in range(3):
        check('quotient_coefficient_'+str(d),coeff[d],expected[d]);check('quotient_equation_'+str(d),residual[d],np.zeros(2,complex))
    reject('omit_incident_j1',coeff[1],np.diag(num[1])/np.diag(den[0]))
    reject('omit_incident_j2',coeff[2],(np.diag(num[2])-np.diag(den[1])*coeff[1])/np.diag(den[0]))
    complex_result,_=cur['quotient']({0:arr(1j)},{0:arr(2)},degree=0)
    check('native_complex_quotient_retained',complex_result[0],np.array([.5j]))
    # Observe the old helper's invalid-domain behavior in isolation; do not change it.
    with np.errstate(divide='ignore',invalid='ignore'):
        invalid,_=cur['quotient']({0:arr(1)},{0:arr(0)},degree=0)
    assert not np.isfinite(invalid[0]).all()
    checks.append({'name':'zero_leading_denominator_requires_external_guard','kind':'domain witness','passed':True,'observed':'nonfinite original-helper output, outside correspondence domain'})
    eps=2
    scaled={g:eps*x for g,x in scalar_amp.items()}
    check('epsilon_squared_flux',cur['quadratic'](scaled,scalar_metric)[(0,0)],eps**2*cur['quadratic'](scalar_amp,scalar_metric)[(0,0)])
    reject('epsilon_linear_flux',cur['quadratic'](scaled,scalar_metric)[(0,0)],eps*cur['quadratic'](scalar_amp,scalar_metric)[(0,0)])
    assert all(sha(BASE/p)==h for p,h in before.items())
    report={'status':'PASS','source_sha256':before,'instrument_sha256':sha(Path(__file__)),'selected_ast_sha256':ast_hashes,'runtime':{'python':platform.python_version(),'numpy':np.__version__},'checks':checks,
     'scope':'Supplied finite rectangle, unrestricted current convolution, real path, degree-six quadratic and nonzero-leading complex per-column quotient.',
     'limits':['No production imports or saved scientific operands.','Small binary-exact synthetic fixtures; translation outside Lean kernel.','Native quotient retains complex values and has no zero guard; physical reality and admissible denominators are application obligations.','No parent-theory Taylor claim, numerical solve or strong-edge extrapolation.']}
    (M/'S11_lean_bookkeeping_source_checks.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
