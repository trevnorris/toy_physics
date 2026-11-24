"""Six frozen total shape derivatives with the pointwise smooth gauge at the R4 stop."""
import argparse
import copy
import json
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from mouth_r4_smooth_static import Gauge, LD, frozen, EXPECTED, digest


def run(area, output):
    with TemporaryDirectory() as directory:
        tmp=Path(directory)
        driver,t,fields,labels,order=frozen(area,tmp)
        from mouth_cp10_nested_forces import integral,kernel
        from mouth_cp10_field_shape import Bath
        from mouth_cp10_spatial_crystal import signed_geometry, MODES
        from mouth_r4_depth_precision import extended
        from mouth_r4_eulerian_solid import SpatialCrystal
        body=t.b
        solid=SpatialCrystal(body)
        gauge=Gauge(t,driver.reference_depth,order)
        ref=kernel(tmp)
        axes=gauge.axes
        length=t.shape.length
        depth=axes[1]
        dz=np.where(depth*length<length-.12,(depth*length/(length-.12))**25,1.)-depth
        energy,grad,_,_=solid.evaluate(labels,order)
        args=(solid.points,solid.weights)
        native=integral(ref,body,axes,dz,labels,order,args,True)
        if abs(energy-native['energy'])>1e-12:
            raise ArithmeticError('Frozen crystal accumulation differs')
        elastic=np.array(native['analytic'])
        elastic[0]-=np.sum(grad[...,1]*dz[None,:,None],dtype=LD)
        bath=Bath(t,driver.base,fields,driver.reference_depth)
        ctx=(driver.prep,driver.reference_depth,body,t.shape.params,t.shape.extra,t.polar_shift)
        def field_parts(dx):
            parts=bath.energy(dx)
            # All shape changes perturb the extensive reservoir integrals.
            # Reevaluate those two terms in promoted arithmetic before differencing.
            shape=copy.copy(t.shape)
            dd=np.asarray(dx,dtype=np.result_type(dx,LD))
            shape.params=t.shape.params+dd[:4]*[1,.25,.25,.25]
            shape.extra=t.shape.extra+.25*dd[4]
            extended_bath=copy.copy(bath)
            extended_bath.model=copy.copy(bath.model)
            extended_bath.model.polar_shift=t.polar_shift+dd[5]
            parts[3:5]=extended(extended_bath,shape,LD(length)*np.exp(dd[0]),driver.reference_depth)
            return parts
        rows=[]
        for i,unit in enumerate(np.eye(6)):
            tangent=unit*1e-24j
            bath_force=-field_parts(tangent).imag.sum()/LD(1e-24)
            gauge_force=-gauge.evaluate(labels,tangent).sum().imag/LD(1e-24)
            analytic=elastic[i]+bath_force
            row=dict(mode=MODES[i],elastic=float(elastic[i]),bath=float(bath_force),unit_gauge=float(gauge_force),banks=[])
            for h in (2e-5,1e-5):
                field,cryst,regularizer=[],[],[]
                for sign in (-1,1):
                    dx=sign*h*unit
                    field.append(field_parts(dx))
                    body2=signed_geometry(ctx,i,sign*h)
                    next_axes=tuple(body2.buffers['field_'+v] for v in 'xta')
                    shifted=labels.copy()
                    shifted[...,1]+=(next_axes[1]-axes[1])[None,:,None]
                    cryst.append(integral(ref,body2,next_axes,dz,shifted,order,args,False)['energy'])
                    regularizer.append(gauge.evaluate(labels,dx))
                fc=(LD(cryst[0])-LD(cryst[1]))/(2*LD(h))
                fb=np.sum(field[0]-field[1])/(2*LD(h))
                fg=np.sum(regularizer[0]-regularizer[1])/(2*LD(h))
                row['banks'].append(dict(h=h,crystal=float(fc),bath=float(fb),unit_gauge=float(fg),
                    epsilon=[dict(epsilon=e,analytic=float(analytic+e*gauge_force),central=float(fc+fb+e*fg),
                                  signed_error=float((fc-elastic[i])+(fb-bath_force)+e*(fg-gauge_force)))
                             for e in (1e-8,1e-10)]))
            rows.append(row)
        data=dict(command=sys.orig_argv,source_sha256=digest(Path(__file__)),state_sha256=EXPECTED['held_stop'],
                  held_fields=True, stationary=False, shapes=rows,
                  precision='Promoted energy sums, extensive reservoir order terms, gauge shape and signed differences. Native two-metric crystal kernel and tabulated vortex use float64 inputs.',
                  max_abs_at_each_width=[max(abs(probe['signed_error']) for row in rows for probe in row['banks'][j]['epsilon']) for j in (0,1)])
        output.write_text(json.dumps(data,indent=2)+'\n')
        print(data['max_abs_at_each_width'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('area',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    run(args.area.resolve(),args.output.resolve())
