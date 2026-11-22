from pathlib import Path
import math, json, itertools, hashlib
import numpy as np
def load(path): return json.loads(Path(path).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
base=Path('_scratch/s11c/selected_shape_sensitivity01')
measure=Path('research/pde_ledger_v3/_measurements')
m=load(measure/'mouth_selected_shape_sensitivity.json')
def verify_guard(d,code=0):
 out=load(d/'outcome.json')
 assert out['exitCode']==code and out['childOutcome']['exitCode']==code
 assert load(d/'child-outcome.json')['exitCode']==code
 lim=load(d/'limit-validation.json')
 assert lim['verified'] and int(lim['actual']['memory.swap.max'])==0
 duration=load(d/'duration-validation.json')
 assert duration['verified'] and duration['wallDeadlineSeconds'] is None
for name,d in m['runs'].items():
 assert d['stdout_sha256']==sha(base/name/'stdout')
 assert d['as_run_main_sha256']==sha('software/stage1_solver/tools/mouth_selected_shape_sensitivity.py' if name.endswith('_corner') else base/'first_solver.py')
 assert d['literal_stdout']==(base/name/'stdout').read_text().splitlines()
for k,h in m['sha256'].items(): assert sha(k)==h
def geometric_area(R,h,b):
 a=np.linspace(0,np.pi/2,3001); t=np.linspace(0,1,3001)
 rc=1+.5*np.sin(np.arctan(.25/(4*np.pi*.4)))
 neck=rc-.5; rise=1-neck
 r1=rc-.5*np.sin(a); z1=.5*(1-np.cos(a))
 r2=neck+rise*(10*t**3-15*t**4+6*t**5)+.49*t*t*(1-t)**3
 r3=(1+.5*(3*t**2-2*t**3))*(1+b*64*(t*(1-t))**3)
 radius=R*np.concatenate((r1,r2[1:],r3[1:]))
 depth=h*np.concatenate((z1,.5+.7*t[1:],1.2+1.85*t[1:]))
 return 4*np.pi*np.sum(np.hypot(np.diff(depth),np.diff(radius))*(radius[:-1]**2+radius[:-1]*radius[1:]+radius[1:]**2)/3)
obs={}
max_area=0
for d in sorted(k.parent for k in base.glob('r*/outcome.json') if k.parent.name!='readout'):
 name=d.name
 if name in m['excluded']:
  verify_guard(d,1);assert (d/'stdout').stat().st_size==0; assert 'self.bulk' in (d/'stderr').read_text()
  continue
 verify_guard(d)
 raw=[json.loads(i) for i in (d/'stdout').read_text().splitlines()]
 shape=next(i for i in raw if i['test']=='shape');st=raw[0]
 assert st['test']=='stationary_action' and abs(st['feed']-st['bulk'])<1e-8 and st['bias']==.4510084282857874
 modes=[i for i in raw if i['test']=='refined_dipole_residual']
 assert len(modes)==4 and all(0<=e['errors'][-1]<3e-8 for e in modes)
 area=geometric_area(shape['radius'],shape['depth_scale'],shape['bulge'])
 max_area=max(max_area,abs(area/shape['material_wall_area']-1))
 assert max_area<1e-6
 windows=[i for i in raw if i['test']=='shape_response'];assert len(windows)==2
 for w in windows: assert abs(w['C']+w['B']+w['T'])<5e-8 and w['D']>0 and w['M']>0
 obs[shape['radius'],shape['mode'],round(shape['step'],6)]=(shape,windows)
assert len(m['excluded'])==2 and len(obs)==16
for radius,mode in itertools.product((1,2),('depth','bulge')):
 for h in (.025,.05):
  minus,plus=(obs[radius,mode,z] for z in (-h,h))
  logstep=math.log(plus[0]['depth_scale']/minus[0]['depth_scale']) if mode=='depth' else math.log((1+h)/(1-h))
  e=(plus[0]['finite_phase_energy']-minus[0]['finite_phase_energy'])/(2*h)
  a=(plus[0]['material_wall_area']-minus[0]['material_wall_area'])/(2*h)
  s=[math.log((wp['D']/wp['M'])/(wm['D']/wm['M']))/logstep for wm,wp in zip(minus[1],plus[1])]
  row=next(rr for rr in m['signed_pairs'] if rr['radius']==radius and rr['mode']==mode and rr['amplitude']==h)
  assert all(abs(v-w['D_over_M'])<2e-8 for v,w in zip(s,row['log_slopes']))
  print(radius,mode,h,s,e,a)
print('maximum independent wall-area discrepancy',max_area)
pbase=Path('_scratch/s11c/free_shape_radius01/energy')
verify_guard(pbase)
vals=load(pbase/'stdout')
for l in vals['rigid_reference_ports']:
 assert math.isclose(l['feed'],l['bulk'],abs_tol=2e-10) and l['native_net_port_power']==0
 assert abs(l['minimum_power_times_conversion_capacity_at_this_rate']-l['feed']*l['bulk'])<1e-14
for name,key in vals['input_sha256'].items(): assert key==sha(name)
c=vals['known_order_field']
assert math.isclose(c['intake']**2/c['capacity'],c['lower_bound']) and 0<c['lower_bound']<c['measured_loss']
for c in vals['powered_annulus_controls']: assert abs(c['power_per_inner_area']/c['loss_per_inner_area']-1)<4e-10
print('shape and energy checks passed')
