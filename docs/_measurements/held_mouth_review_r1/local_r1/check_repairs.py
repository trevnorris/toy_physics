"""Review the identified held-mouth repairs on frozen source; start with s11c_guarded_run.py."""
from contextlib import redirect_stdout
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import sys
import traceback
import numpy as np
base = Path(__file__).resolve().parent
frozen = base.parent/'r1/candidate'
tool = frozen/'software/stage1_solver/tools'
sys.dont_write_bytecode = True
sys.path[:0] = [str(tool), str(Path('software/stage1_solver/tools').resolve())]
for line in (frozen.parent/'candidate.sha256').read_text().splitlines():
    h,p=line.split(maxsplit=1)
    assert sha256((frozen/p).read_bytes()).hexdigest()==h
print('frozen 10 files',len((frozen.parent/'candidate.sha256').read_text().splitlines()))
for p in sorted((frozen/'docs/_measurements').glob("*.txt")):
    s=p.read_text(); path,size,h=re.search(r'Stdout (\S+) \((\d+) bytes, SHA256 (\w+)',s).groups()
    raw=Path(path).read_bytes()
    assert raw == s.split("Literal stdout\n")[1].encode() and len(raw)==int(size) and sha256(raw).hexdigest()==h
    outcome=json.loads((Path(path).parent/'outcome.json').read_text())
    assert outcome['limitsVerified'] and outcome['exitCode']==outcome['childOutcome']['exitCode']==0
    hlist=re.findall(r'^([0-9a-f]{64})  (.+)$',s,re.M)
    for k,v in hlist:
        x=frozen/v if (frozen/v).is_file() else Path(v)
        assert sha256(x.read_bytes()).hexdigest()==k
    print('verified',p.stem,'stdout bytes',size,'listed hashes',len(hlist))
u='unequal_mouth_boundary_work'; b='held_crossflow_budget'; t='mouth_transmission_work'; c='mouth_confinement_projection'
cases=[
 (u,'freeze_prepared_normalizer','0.25 / np.sum(ring_volume * result)','0.25 / np.sum(ring_volume * isolated_current * (pair/alone)**2)','assert_allclose(uniform, uniform_expected'),
 (u,'swap_fixed_law','fixed_rate=(i == 0)','fixed_rate=(i == 1)','assert_allclose(current, prepared_current'),
 (u,'corrupt_paired_j','current = currents[i][indices]','current = 1.01*currents[i][indices]','assert_allclose(current, prepared_current'),
 (u,'drop_reference','((pair - 1)**2 - (alone - 1)**2) / 2','(pair-1)**2 / 2','assert_allclose(f, fd'),
 (u,'wrong_annular_measure','dv = 2 * np.pi * r * dx**2','dv = 2 * np.pi * dx**2','assert_allclose('),
 (b,'flip_feed','feed = incident * rate','feed = -incident * rate','assert_allclose(odd[0],'),
 (b,'flip_boundary_stress','opening = carry + stress','opening = carry - stress','assert_allclose(odd[1],'),
 (b,'flip_layer_stress','direct = direct_carry + direct_stress','direct = direct_carry - direct_stress','assert_allclose(direct - direct_carry,'),
 (b,'erase_layer_stress','direct = direct_carry + direct_stress','direct = direct_carry','assert_allclose(direct - direct_carry,'),
 (b,'omit_mask_rim','dx4(abs(field[-1])**2 * mouth, dx)','dx4(abs(field[-1])**2, dx) * mouth','assert_allclose(direct - direct_carry,'),
 (b,'alter_layer_carry','direct = direct_carry + direct_stress','direct_carry = 1.01*direct_carry\n        direct = direct_carry + direct_stress','assert_allclose(odd[2],'),
 (b,'halve_secant','odd = (pos - neg) / (2 * speed)','odd = (pos - neg) / (4 * speed)','assert_allclose(odd[0],'),
 (t,'nonlinear_trace_chain','return sp.diff(value, y) + nx * sp.diff(value, n)','return sp.diff(value, y)','Interface momentum balance failed'),
 (t,'one_side_quantum','return carry + c*','return carry - c*','Interface momentum balance failed'),
 (t,'change_general_lambda','receiving = sp.factor(-g * (amp-1) + amp_jump)','receiving = sp.factor(-g * (amp-1) + amp_jump) + g*lam*(lam-1)','Interface momentum balance failed'),
 (t,'zero_law','receiving = sp.factor(-g * (amp-1) + amp_jump)','receiving = lam*sp.factor(-g * (amp-1) + amp_jump)','Interface momentum balance failed'),
 (t,'double_analytic_odd','/ (4*speed)','/ (2*speed)','assert_allclose(slope, differenced'),
 (t,'zero_analytic_odd','density_work = (','density_work = 0*(','assert_allclose(slope, differenced'),
 (t,'omit_negative_energy','signed_force(amplitude[::-1])','signed_force(amplitude)','assert_allclose(slope, differenced'),
 (t,'off_aperture_current','if (drive != speed','current[0,0,0] += 1e-4; current[0,0,1] -= 1e-4\n    if (drive != speed','Unexpected mouth geometry'),
 (c,'do_not_evaluate_zero','sp.diff(radial_force, speed).subs(speed, 0).subs(vr, vr_symbol)','sp.diff(radial_force, speed).subs(vr, vr_symbol)','Inconsistent angular projection'),
 (c,'incorrect_cubic_scale','sp.diff(radial_force, speed, 3)/6','sp.diff(radial_force, speed, 3)/2','Inconsistent angular projection'),
 (c,'non_radial_baseline','density = n0 + speed','density = n0*(1+mu) + speed','Inconsistent angular projection'),
 (c,'endpoint_sign','near_force = sp.diff(elastic, slope)','near_force = -sp.diff(elastic, slope)','Inconsistent angular projection'),
]
def run_mutant(mod,key,old,new,expects):
    src=(tool/(mod+'.py')).read_text()
    assert src.count(old)==1,(key,old)
    dest=base/'mutants';dest.mkdir(exist_ok=True)
    p=dest/(mod+'_'+key+'.py');p.write_text(src.replace(old,new))
    log=io.StringIO()
    try:
        with redirect_stdout(log):
            exec(compile(p.read_text(),str(p),'exec'),{'__name__':'__main__','__file__':str(p)})
    except (AssertionError,ArithmeticError,ValueError) as err:
        tr=traceback.format_exc()
        (p.with_suffix('.trace')).write_text(tr)
        assert expects in tr,(key,expects,tr)
        print('REJECT',key,'at',expects,flush=True)
    else:
        raise AssertionError(('unexpected survival',key))
    finally:p.with_suffix('.stdout').write_text(log.getvalue())
for case in cases:run_mutant(*case)
print('verified discriminating mutations',len(cases))
# Exercise the general variance expression on a noncentered admissible trace; the
# archived drifts have a nearly zero dipole mean and the uniform current witnesses
# check a different variation.
import unequal_mouth_boundary_work as uw
x=(np.arange(18)-8.5)/7
r=(np.arange(9)+.5)/7
X=x[:,None];cut=np.maximum(1-X**2-r*r,0);j=(1+.76*X)*cut**3
dv=np.pi*((r+1/14)**2-(r-1/14)**2)/7
j*=.25/np.sum(j*dv)
amp=1+.05*X
_,_,var,fd,uniform=uw.probe(x,r,1,amp,amp,j,j,True)
w=X*cut**3; wrong=2*np.sum(dv*j*w*w)
print('skew test correct variance / probed / raw uncentered / uniform',var,fd,wrong,uniform)
assert abs(var-fd)<1e-12 and abs(var-wrong)>1e-5 and abs(uniform)<1e-10
