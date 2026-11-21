"""Test one-way current through a moving-coordinate aperture and its held-out slope."""
import argparse
import json
from pathlib import Path
import sympy as sp


def patch():
    # A genuinely varying, divergence-free local flow. The x,y coordinates
    # are in the brane; an unused third in-brane coordinate is homogeneous.
    x,y,z,w,a,b=sp.symbols('x y z w a b',real=True)
    h=a*(1-x*x-y*y)*(1+b*x)
    jx=x+z
    jy=-2*y+z
    jz=1+z+x*y
    assert sp.diff(jx,x)+sp.diff(jy,y)+sp.diff(jz,z)==0
    # Active transport of the oriented number-current element by w=z+h(x,y).
    lift=jz+jx*h.diff(x)+jy*h.diff(y)
    physical=[f.subs(z,w-h) for f in (jx,jy,lift)]
    div=sp.simplify(sum(f.diff(v) for f,v in zip(physical,(x,y,w),strict=True)))
    crossing=sp.factor((physical[2]-physical[0]*h.diff(x)-physical[1]*h.diff(y)).subs(w,h))
    print('reference crossing, graph crossing, lab divergence =',jz.subs(z,0),crossing,div)
    assert crossing==jz.subs(z,0) and div==0
    # Test a property this construction need not preserve. Initially potential
    # flow v=e_x, unit density, h=a(1-r²) near the mid-plane y=0.
    curl_y=sp.simplify(-sp.diff(sp.diff(a*(1-x*x-y*y),x),x))
    print('resting-coordinate uniform stream curl_y, lifted Euclidean curl_y =',0,curl_y)
    assert curl_y==2*a


def guarded_lines(directory):
    outcome=json.loads((directory/'outcome.json').read_text())
    assert outcome['exitCode']==outcome['childOutcome']['exitCode']==0 and outcome['limitsVerified']
    return [json.loads(l) for l in (directory/'stdout').read_text().splitlines() if l.startswith('{')]


def predict(front_dir,bulk_dir):
    rows=guarded_lines(front_dir)
    front={row['d']:row for row in rows if row['grid']==[96,6,24] and row['elements']==32 and not row['transport']}
    bulk={row['d']:row for row in guarded_lines(bulk_dir) if 'jump_moment' in row}
    assert set(front)==set(bulk)=={8,16}
    # In the lifted current closure the pressure jump is sigma/K times the
    # measured current moment. Use 8R only to set this single ratio and keep
    # 16R entirely out of that fit. The receiver here has radial discharges.
    measured={d:r['pair']['flat_current_moment']-r['alone']['flat_current_moment'] for d,r in front.items()}
    ratio=measured[16]/measured[8]
    scale=bulk[8]['jump_moment']/measured[8]
    prediction=scale*measured[16]
    print(json.dumps(dict(material_current_moments=measured,independent_receiver_control_sigma_over_K_fit8=scale,
        current16_over_8=ratio,control_jump16_over_8=bulk[16]['jump_moment']/bulk[8]['jump_moment'],
        predicted_jump16=prediction,control_jump16=bulk[16]['jump_moment'],
        prediction_over_control=prediction/bulk[16]['jump_moment']),allow_nan=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('front_run',type=Path)
    parser.add_argument('bulk_run',type=Path)
    args=parser.parse_args()
    patch()
    predict(args.front_run,args.bulk_run)
