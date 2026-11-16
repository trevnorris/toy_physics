"""One-sided source/field perturbations on private, retained copies."""
from pathlib import Path
import contextlib, io, sys, json, traceback
import numpy as np
base=Path(__file__).resolve().parent
cand=base.parent/'r0/candidate/software/stage1_solver/tools'
sys.path.insert(0,str(cand)); sys.path.insert(1,str(Path('software/stage1_solver/tools').resolve()))
U='unequal_mouth_boundary_work'; B='held_crossflow_budget'; T='mouth_transmission_work'; C='mouth_confinement_projection'
tests=[]
def mut(s,key,old,new,shouldfail,reason):
    tests.append((s,key,old,new,shouldfail,reason))

mut(U,"sign", "force = -np.sum","force = np.sum",True,"force/difference")
mut(U,"mixed","mixed = 2 *","mixed = 3 *",True,"mixed/difference")
mut(U,"zero_variation","eta = axial / radius * base(0)**3","eta = axial * 0",True,"degenerate variation")
mut(U,"calibration","if right[\"reference_fields\"] != pair[\"reference_fields\"]:","right = {**right, 'reference_fields':{}}\n        if right[\"reference_fields\"] != pair[\"reference_fields\"]:",True,"matched reference")
mut(U,"pair_geometry",'if int(data["d"]) != distance','if int(data["d"]) + 1 != distance',True,"d")
mut(U,"pair_radius",'float(data["target_radius"]) != 1.5','float(data["target_radius"])+1 != 1.5',True,"pair radius")
mut(U,"left_center",'control["centers"].tolist() != [center]','control["centers"].tolist() != [center+1]',True,"off-center control")
mut(U,"left_rate",'control["rates"].tolist() != [0.25]','control["rates"].tolist() != [0.125]',True,"left rate")
mut(U,"right_radius",'float(control["radius"]) != radius','float(control["radius"])+1 != radius',True,"right radius")
mut(U,"x_match",'np.testing.assert_array_equal(x, xc)','xc = xc + 1\n                    np.testing.assert_array_equal(x, xc)',True,"local x coordinates")
mut(U,"r_match",'np.testing.assert_array_equal(r, rc)','rc = rc + 1\n                    np.testing.assert_array_equal(r, rc)',True,"local r coordinates")
mut(U,"wrong_rate","rate = np.sum(dv * current)","current = np.array(current)*2\n    rate = np.sum(dv * current)",False,"current array doubled after checksum")
mut(U,"wrong_law","fixed_rate=(i == 0)","fixed_rate=(i == 1)",False,"left/right phase laws exchanged")
mut(U,"wrong_measure","dv = 2 * np.pi * r * dx**2","dv = np.pi * dx**2 * np.ones_like(r)",False,"annular Jacobian removed both sides")
mut(U,"fixed_rate_fd_omitted","return result * (rate / np.sum(dv * result)) if fixed_rate else result","return result",None,"one-sided fixed-rate normalization deletion")
mut(U,"flat_baseline","force = -np.sum","force = -np.sum",False,"placeholder") # replace with actual shared reference perturbation
tests[-1]=(U,"unit_baseline","difference = ((pair - 1)**2 - (alone - 1)**2) / 2","difference = (pair - 1)**2 / 2",False,"substitutes uniform isolated amplitude instead of matched state")
mut(B,"feed_scale","feed = incident * rate","feed = -incident * rate",False,"feed/aperture sign not checked")
mut(B,"aperture_stress","opening = carry + stress","opening = carry - stress",False,"one-sided stress sign")
mut(B,"direct_flux","direct = direct_carry + direct_stress","direct = direct_carry - direct_stress",False,"one-sided layer stress sign")
mut(B,"odd_normalizer","(2 * speed)","(4 * speed)",False,"difference denominator")
mut(B,"parity","neg, -pos, rtol=", "neg + 0.1, -pos, rtol=",True,"antisymmetry")
mut(B,"guard","not guard[\"limitsVerified\"]","not False",True,"guard")
mut(B,"checksum","!= digest:","!= digest[::-1]:",True,"digest")
mut(B,"calibration_speed",'float(data["speed"]) != 0','float(data["speed"]) + 1 != 0',True,"calibration metadata")
mut(B,"calibration_parity","np.testing.assert_allclose(zero_density,","zero_density[0] *= 1.01\n    np.testing.assert_allclose(zero_density,",True,"zero-drive parity")
mut(B,"speed","float(data[\"speed\"]) != speed","float(data[\"speed\"])+1 != speed",True,"speed metadata")
mut(B,"field_shape","field.shape != (65, 32, 32, 32)","field.shape != (64, 32, 32, 32)",True,"depth")
mut(B,"current_shape","raw = weight *","current = np.roll(current,1,0)\n        raw = weight *",True,"nonuniform current shape")
mut(B,"support","if not np.all(current[mouth] > 0)","current[~mouth] += 1e-16\n        if not np.all(current[mouth] > 0)",True,"strict off-aperture support")
mut(B,"inward","if not np.all(current[mouth] > 0)","current[mouth] = -current[mouth]\n        if not np.all(current[mouth] > 0)",True,"inward")
mut(B,"reflect_validation","field = field[:, ::-1]","field = field[:, ::-1].copy(); field[2,3,4,5] += 0.01",True,"signed stationary PDE")
mut(B,"feed_phase","field = field[:, ::-1]","field = field[:, ::-1].copy(); field[0,0,0,0] *= 1.01",True,"feed envelope")
mut(B,"mouth_mask","field = field[:, ::-1]","field = field[:, ::-1]; mouth = np.roll(mouth,1,0)",True,"mask")
mut(B,"rate_invariant","field = field[:, ::-1]","field = field[:, ::-1]; dz = np.nan",None,"nonfinite normal spacing in signed reflected partner")
mut(T,"imaginary_jump","condition = c * derivative_jump -", "condition = -sp.I*amp + c * derivative_jump -",True,"current continuity")
mut(T,"nonlinear","return sp.diff(value, y) + nx * sp.diff(value, n)","return sp.diff(value, y)",False,"varying density suppressed in both chain rules")
mut(T,"log_density","drive = sp.simplify(c * n * dx(log_jump) / 2)","drive = sp.simplify(c * n * dx(log_jump) / 3)",True,"surface stress coefficient")
mut(T,"edge_sign","translated = sp.simplify(drive - dx(edge_flux))","translated = sp.simplify(drive + dx(edge_flux))",True,"one-sided divergence")
mut(T,"lambda_one","receiving = sp.factor(-g * (amp-1) + amp_jump)","receiving = sp.factor(-g * (amp-1) - amp_jump)",True,"receiving at 1")
mut(T,"lambda_zero","receiving = sp.factor(-g * (amp-1) + amp_jump)","receiving = sp.factor(-g * (amp-1) + amp_jump)*lam",False,"lambda=0 and finite lambda")
mut(T,"one_side_force","slope = dx**3 * np.sum","slope = -dx**3 * np.sum",True,"energy difference")
mut(T,"odd_norm","/ (4*speed)","/ (2*speed)",False,"odd sign/speed normalization both branches")
mut(T,"profile_null","density_work = (","density_work = 0 * (",False,"erase numerical force and energy")
mut(T,"rate","rate = dx**3 * np.sum(current)","rate = dx**3 * np.sum(current)*0.5",True,"inward rate")
mut(T,"speed","drive != speed","drive != speed+1",True,"drive")
mut(T,"radius","mouth, base(0) > 0","np.roll(mouth,1,0), base(0) > 0",True,"support geometry")
mut(T,"inward","np.all(current[mouth] > 0)","np.all(-current[mouth] > 0)",True,"inward")
mut(T,"support_leak","if drive != speed","current[0,0,0] += 0.0001; current[0,0,1] -= 0.0001\n    if drive != speed",False,"nonzero signed flux outside aperture, total J fixed")
mut(C,"rim","rim = x * mu +","rim = 2*x * mu +",True,"rim and body")
mut(C,"volume","volume_shift = 2 * np.pi" if False else "volume_shift = 2 * sp.pi","volume_shift = 1 + 2 * sp.pi",True,"volume")
mut(C,"reaction","reaction = sp.simplify(work -","reaction = sp.simplify(work +",True,"force on geometry sign")
mut(C,"odd_measure","2 * sp.pi * sp.integrate(work,","3 * sp.pi * sp.integrate(work,",True,"radial measure")
mut(C,"trace","second = first + s**2 * (1 - s)**2","second = first + s * (1 - s)**2",True,"equal endpoint slopes")
mut(C,"interior","second = first + s**2 * (1 - s)**2","second = first + 2*s**2 * (1 - s)**2",True,"interior shape")
mut(C,"collar","collar_work = sp.integrate(","collar_work = 1+sp.integrate(",True,"shrinking regular load")
mut(C,"endpoint","contact_weight = collar.subs(w, 0)","contact_weight = collar.subs(w, width)",True,"endpoint evaluation")
mut(C,"elastic","near_force = sp.diff(elastic, slope)","near_force = -sp.diff(elastic, slope)",True,"lower-end sign")
mut(C,"undriven","density = n0 + speed","density = n0 * (1+mu) + speed",False,"undriven radial-force statement, anisotropic input")
mut(C,"affine_only","density = n0 + speed * (n1 * mu + n2 * sp.legendre(2, mu))","density = n0 + speed * (n1 * mu + n2 * sp.legendre(2, mu)) + speed**3*mu*n1",True,"nonlinear odd U contamination through missing evaluation at U=0")
def run():
    dest=base/'mutated'; dest.mkdir(exist_ok=True)
    bad=[]
    for mod,key,old,new,expect,why in tests:
        src=(cand/(mod+'.py')).read_text()
        assert src.count(old)==1,(mod,key,old)
        path=dest/(mod+'_'+key+'.py'); path.write_text(src.replace(old,new))
        ns={"__name__":"__main__","__file__":str(path)}; buf=io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                exec(compile(path.read_text(),str(path),'exec'),ns)
            died=False; detail='survived'
        except Exception as ex:
            died=True
            detail=type(ex).__name__+' '+str(ex).replace('\n',' ')[:180]
        path.with_suffix('.stdout').write_text(buf.getvalue())
        print(f"{mod} {key}: {detail} | {why}",flush=True)
        if expect is not None and expect!=died:
            bad.append((key,died,expect))
    print("unexpected change in outcome",bad)
    if bad: raise ValueError(bad)
if __name__=='__main__':run()
