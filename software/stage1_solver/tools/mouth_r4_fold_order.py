"""Read crystalline order on saved radius-four placements and the identified near-singular prism."""
import argparse
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import RegularGridInterpolator

FILES = {}


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def input_file(path, expected=None):
    path = path.resolve()
    value = digest(path)
    if expected is not None and value != expected:
        raise ValueError(f'Changed input {path}')
    FILES[str(path)] = value
    return path


def read(path, expected=None):
    return json.loads(input_file(path, expected).read_bytes())


def span(values):
    return [float(np.min(values)), float(np.max(values))]


def axes_model(shape, reference=None):
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_last_rings_fluid import GRIDS
    from mouth_supercore_joined_front import Joined

    class Nodes(Joined):
        def matrix(self, order, elastic_density):
            # Joined has finished numbering all six blocks. No field equation is needed.
            pass

    m = Nodes(shape, None, 1, GRIDS['double'], precondition=False)
    if reference is not None:
        rr, z = axes_for_depth(shape.length, reference)
        if any(len(a) != len(b) or a[0] != b[0] or a[-1] != b[-1]
               for a, b in zip((*rr,*z), (*m.rr,*m.z))):
            raise ValueError('Changed spatial topology')
        m.rr, m.z = rr, z
    return m


def body_axes(solid, registration):
    from mouth_free_shape_solid import Solid
    nt,nx,na,_ = solid.shape
    b = SimpleNamespace(cut=np.exp(-.5*np.exp(registration)))
    if (nx-1)%3:
        raise ValueError('Unexpected inner mesh')
    b.x = Solid.radial_mesh(b, (nx-1)//3)
    b.t = np.linspace(0, 1, nt)
    b.a = np.linspace(0, np.pi/2, na)
    return b


def chart_data(m, absolute, base=None):
    axes = (m.rr[0], m.a, m.z[1]/m.shape.length)
    values = absolute[m.ids[1,0]]
    if base is not None:
        # The earlier solid sampled the sum of two independently gridded fields.
        # Overlay the meridional basis knots, on which that sum is again Q1.
        bt,bx = (np.linspace(0, 1, n) for n in base.shape)
        xx,aa,tt = np.meshgrid(np.union1d(axes[0],bx),axes[1],np.union1d(axes[2],bt),indexing='ij')
        values = (RegularGridInterpolator(axes,values)((xx,aa,tt)) +
                  RegularGridInterpolator((bt,bx),base)((tt,xx)))
        axes = xx[:,0,0], aa[0,:,0], tt[0,0,:]
    return axes, values


def all_low_after(t, values, threshold):
    """Endpoint of last above-threshold interval, all x and angle, then <= threshold."""
    if np.any(values[...,-1]>threshold):
        return None
    ids = np.nonzero(np.any((values>threshold).reshape(-1,len(t)),axis=0))[0]
    if not len(ids):
        return 0.
    j = ids[-1]
    a,b = values[...,j],values[...,j+1]
    keep = a>threshold
    return float(t[j]+(t[j+1]-t[j])*np.max((a[keep]-threshold)/(a[keep]-b[keep])))


def radial_sections(axes, values, depths):
    x,a,t = axes
    result = []
    for h in depths:
        j = np.clip(np.searchsorted(t,h,side='right')-1,0,len(t)-2)
        w = (h-t[j])/(t[j+1]-t[j])
        c = (1-w)*values[:,:,j]+w*values[:,:,j+1]
        result.append(dict(t=h, disk=span(c), axis=float(c[0,0]), wall=span(c[-1])))
    return result


def describe(axes, values):
    x,a,t = axes
    out = {}
    sets = dict(disk=values, axis=values[:1,:1], wall=values[-1:])
    for name, arr in sets.items():
        out[name] = {str(h): all_low_after(t, arr, h) for h in (.5,.1,.01,0.)}
    # One physical center per z rather than duplicated polar coordinate entries.
    unique = np.vstack((values[0,:1],values[1:].reshape(-1,len(t))))
    fractions = {str(h): int(np.sum(unique<=h)) for h in (.5,.1,.01,0.)}
    return dict(terminal_from_fraction=out,shaft_distinct_nodes=unique.size,shaft_nodes_at_or_below=fractions,
                sections=radial_sections(axes,values,[0.,.1,.25,.5,.75,.97765,1.]))


def physical(shape, shift, chart):
    from mouth_joined_angular_work import bump
    x,t,a = np.moveaxis(np.asarray(chart),-1,0)
    z = t*shape.length
    r = x*shape.radius(z)[0]+shift*x**3*bump(z,shape.length)[0]*(3*np.cos(a)**2-1)/2
    return np.stack((r,z,a),axis=-1)


def mesh_ring(b, solid, axes, values, m, shift):
    from mouth_joint_endpoint_audit import material_map
    x,t,a = np.broadcast_arrays(b.x[-2],b.t[-2],b.a)
    chart = np.stack((x,t,a),axis=-1) + solid[-2,-2]
    mapped = physical(m.shape,shift,chart)
    check = material_map(m.shape,shift,b,solid,t,x,a)
    if np.max(abs(check-mapped))>2.e-12:
        raise ArithmeticError('Ring differs from native physical map')
    chi = RegularGridInterpolator(axes,values)(chart[:,[0,2,1]])
    return dict(reference_t_x=[float(b.t[-2]),float(b.x[-2])],physical_r=span(mapped[:,0]),
                landed_t=span(chart[:,1]),chi_at_nodes=chi.tolist(),
                enclosing_all_segments=box_bound(axes, values, chart[:,[0,2,1]]))


def exact(value):
    return Fraction.from_float(float(value))


def outward(y, lower):
    v = float(y)
    if (lower and exact(v)>y) or (not lower and exact(v)<y):
        v = float(np.nextafter(v, -np.inf if lower else np.inf))
    return v


def exact_value(axes, values, point):
    rules = []
    for axis,q in zip(axes,point):
        i = int(np.clip(np.searchsorted(axis,float(q),side='right')-1,0,len(axis)-2))
        # A rational coordinate may round onto the adjacent double knot.
        if i>0 and q<exact(axis[i]):
            i -= 1
        w = (q-exact(axis[i]))/(exact(axis[i+1])-exact(axis[i]))
        if not 0<=w<=1:
            raise ValueError(f'Outside interpolation grid: {q}')
        rules.append(((i,1-w),(i+1,w)))
    total = Fraction()
    for (i,v),(k,w),(j,s) in product(*rules):
        if v and w and s:
            total += v*w*s*exact(values[i,k,j])
    return total


def box_bound(axes, values, corners):
    """Exact Q1 extrema on box enclosing all mapped points (P1 triangle times Q1 angle)."""
    samples = []
    for d,ax in enumerate(axes):
        lo,hi = min(corners[:,d]),max(corners[:,d])
        lo,hi = Fraction(lo),Fraction(hi)
        samples.append(sorted({lo,hi,*(exact(v) for v in ax if lo<exact(v)<hi)}))
    vals = [exact_value(axes,values,p) for p in product(*samples)]
    return [outward(min(vals),True),outward(max(vals),False)]


def prism(b, solid, axes, values, m, shift, where, detailed=False):
    from mouth_joint_endpoint_audit import material_map
    j,i,k = where['cell_t_x_angle']
    chosen = ((1,1),(0,1),(1,0)) if where['triangle'] else ((0,0),(1,0),(0,1))
    labels = np.array([[[b.x[i+di],b.t[j+dj],b.a[k+dk]] for dk in (0,1)] for di,dj in chosen])
    offsets = np.array([[solid[j+dj,i+di,k+dk] for dk in (0,1)] for di,dj in chosen])
    corners = labels+offsets
    f = where['angular_fraction']
    idx = where['vertex']
    q = (1-f)*labels[idx,0]+f*labels[idx,1]
    p = (1-f)*corners[idx,0]+f*corners[idx,1]
    native = material_map(m.shape,shift,b,solid,np.array([q[1]]),np.array([q[0]]),np.array([q[2]]))[0]
    rz = physical(m.shape,shift,p)
    if np.max(abs(native-rz))>2.e-12:
        raise ArithmeticError('Prism point differs from native physical map')
    # Exact operations on the saved binary nodal labels and offsets. Include
    # both angular faces; an interior parent point is a convex combination.
    rational = np.vectorize(exact,otypes=[object])
    verts = rational(labels)+rational(offsets)
    center = (1-exact(f))*verts[idx,0]+exact(f)*verts[idx,1]
    v = exact_value(axes,values,center[[0,2,1]])
    bound = box_bound(axes,values,verts.reshape(-1,3)[:,[0,2,1]])
    row = dict(reference_x_t_angle=q.tolist(),landed_x_t_angle=p.tolist(),physical_r_z_angle=rz.tolist(),
               depth_fraction=float(p[1]),axis_distance=float(rz[0]),wall_distance=float(
                   physical(m.shape,shift,[1.,p[1],p[2]])[0]-rz[0]),
               point_chi=[outward(v,True),outward(v,False)],whole_prism_chi_enclosure=bound)
    if detailed:
        fe = RegularGridInterpolator(axes,values)
        # Actual samples accompany the enclosing box, which can contain points
        # outside the deformed triangular prism.
        bary = np.array([(l/16,h/16,1-(l+h)/16) for l in range(17) for h in range(17-l)])
        lines = np.einsum('bi,ifd->bfd',bary,corners)
        u = np.linspace(0,1,17)[:,None]
        sample = (1-u)*lines[:,:1,:]+u*lines[:,1:,:]
        chi = fe(sample[:,:,[0,2,1]])
        row.update(vertices_landed_x_t_angle=corners.tolist(),
                   vertices_chi=fe(corners[:,:,[0,2,1]]).tolist(),sample_count=chi.size,
                   sampled_chi_range=span(chi),point_chi_rational=str(v),
                   physical_vertex_r_z_angle=physical(m.shape,shift,corners).tolist())
        for lo,hi in (row['sampled_chi_range'],row['point_chi']):
            if lo < bound[0]-1.e-14 or hi > bound[1]+1.e-14:
                raise ArithmeticError('Order samples outside prism bound')
    return row


def collar(m,absolute):
    t = .97765
    z = t*m.shape.length
    radius = float(m.shape.radius(z)[0])
    if z/m.shape.length<=.62:
        raise ValueError('Active polar wall perturbation at section')
    profiles={}
    for side,name in ((0,'inside'),(1,'outside')):
        x = m.rr[side]
        axis = x,m.a,m.z[1]/m.shape.length
        v = absolute[m.ids[1,side]]
        xx,aa = np.meshgrid(x,m.a,indexing='ij')
        fe = RegularGridInterpolator(axis,v)
        c = fe((xx,aa,np.full_like(xx,t)))
        r = radius*x if not side else radius+(28-radius)*x
        # At fixed depth, Q1 extrema along radius and angle occur on these knots.
        profiles[name] = dict(physical_radial_knots=r.tolist(),angular_min=c.min(axis=1).tolist(),
                              angular_max=c.max(axis=1).tolist(),whole_section_range=span(c))
    return dict(depth_fraction=t,physical_z=z,wall_radius=radius,**profiles)


def plot_profile(axes, values):
    xg = np.linspace(0,1,41)
    tg = np.unique(np.r_[np.linspace(0,1,51),.97765])
    x,a,t = np.meshgrid(xg,axes[1],tg,indexing='ij')
    c = RegularGridInterpolator(axes,values)((x,a,t))
    return dict(radial_chart=xg.tolist(),depth_fraction=tg.tolist(),angle_knots=axes[1].tolist(),
                angle_min_t_by_x=c.min(axis=1).T.tolist(),angle_max_t_by_x=c.max(axis=1).T.tolist())


def shape_from(archive):
    from mouth_sealed_shape_relax import Meridian
    if int(archive['radius']) != 4 or int(archive['winding']) != 1:
        raise ValueError('Unexpected preparation')
    return Meridian(4., archive['shape_params'],float(archive['shape_extra']))


def held_history(args, initial, registration):
    from mouth_continuous_order import seed_order, topology
    fine = args.original/'fine'
    held = args.held/'saved'
    record = read(held/'record.json')
    fine_record = read(fine/'endpoint.json')
    with np.load(input_file(fine/fine_record['file'],record['fine'])) as f:
        shape = shape_from(f)
        solid, shift, narrow_wound = f['solid'],f['x'][5],f['field']
    with np.load(input_file(fine/'endpoint_controls.npz',record['controls'])) as f:
        narrow_rest = f['rest']
    m = axes_model(shape)
    number,_ = topology(m)
    b = body_axes(solid,registration)
    ref = {}
    for tag,field in [('narrow_rest',narrow_rest),('narrow_wound',narrow_wound)]:
        axes,c = chart_data(m,field[2],initial)
        ref[tag] = dict(readout=describe(axes,c),ring=mesh_ring(b,solid,axes,c,m,shift))
    for tag in ('rest','wound'):
        with np.load(input_file(args.wider/(tag+'.npz'))) as f:
            axes,c = chart_data(m,f['field'][2],initial)
            ref['wider_'+tag] = dict(readout=describe(axes,c),ring=mesh_ring(b,solid,axes,c,m,shift))
    # Only the side's wall and receiver's lower disk share a spatial node
    # with the shaft; both traces of the old receiver are pinned at 1/2.
    # Other receiver coefficients do not enter a shaft value. Supply a
    # constant extension of those traces for this shaft-only seed reconstruction.
    def trace(x,z):
        return np.full(np.broadcast_shapes(np.shape(x),np.shape(z)),.5)
    trace.outside = lambda r,z,above: trace(r,z)
    seed = seed_order(m,initial,trace,narrow_rest)
    axes,c = chart_data(m,seed[2])
    ref['continuous_rest_shaft_seed'] = dict(readout=describe(axes,c),ring=mesh_ring(b,solid,axes,c,m,shift))
    for name in ('rest','wound'):
        with np.load(input_file(held/('intermediate_'+name+'.npz'))) as f:
            absolute = f['absolute_order'][number]
            iter_num = int(f['iteration'])
        axes,c = chart_data(m,absolute)
        ref['continuous_'+name+'_iteration_'+str(iter_num)] = dict(readout=describe(axes,c),ring=mesh_ring(b,solid,axes,c,m,shift))
        with np.load(input_file(held/(name+'.npz'),record['solutions'][name]['sha256'])) as f:
            absolute = f['field'][2]
        axes,c = chart_data(m,absolute)
        ref['continuous_'+name] = dict(readout=describe(axes,c),ring=mesh_ring(b,solid,axes,c,m,shift),
                                       fixed_depth_collar=collar(m,absolute))
    return ref


def run(args):
    from mouth_joint_endpoint_audit import orientation
    decision = read(args.gate)
    checkpoint = Path(decision['checkpoint']).resolve()
    run_folder = args.run/'coarse'
    if checkpoint != (run_folder/decision['latest']['file']).resolve():
        raise ValueError('Selected input differs from gate')
    inv = read(args.original/'coarse/invocation.json')['command']
    registered = Path(inv[inv.index('-u')+2])
    rows = [json.loads(s) for s in input_file(run_folder/'stdout').read_text().splitlines()]
    start = next(row for row in rows if row['kind']=='input')
    accepted = {r['iteration']:r for r in rows if r['kind']=='accepted'}
    with np.load(input_file(registered,start['state_sha'])) as f:
        initial, registration = f['order'], float(f['registration'].item())
    final = decision['latest']
    target = final['orientation']['extremizers']['min']
    result = dict(reader_hash=digest(Path(__file__)),command=sys.argv,selected_iteration=decision['iteration'],
                  selected_sha256=decision['sha256'],selected_extremizer=target,held={},coarse=[])
    result['held'] = held_history(args,initial,registration)
    for it in range(decision['iteration']+1):
        row = accepted[it]
        path = run_folder/row['file']
        if it==decision['iteration'] and row['sha256']!=decision['sha256']:
            raise ValueError('Accepted state differs from gate')
        with np.load(input_file(path,row['sha256'])) as f:
            m = axes_model(shape_from(f),start['reference_depth'])
            solid, absolute, shift = f['solid'],f['field'][2],f['x'][5]
        axes, c = chart_data(m,absolute)
        b = body_axes(solid,registration)
        geom = orientation(m.shape,shift,b,solid)
        if geom['extremizers']['min'] != row['orientation']['extremizers']['min'] or (
            geom['prism_chart_determinant_range']!=row['orientation']['prism_chart_determinant_range']):
            raise ValueError(f'Orientation not reproduced for {it}')
        detail = (it==decision['iteration'])
        one = dict(iteration=it,chart_min=geom['prism_chart_determinant_range'][0],depth=m.shape.length,
                   selected_prism=prism(b,solid,axes,c,m,shift,target,detailed=detail),
                   worst_here=geom['extremizers']['min'],readout=describe(axes,c),
                   last_free_ring=mesh_ring(b,solid,axes,c,m,shift))
        if detail:
            one['fixed_depth_collar'] = collar(m,absolute)
            one['profile'] = plot_profile(axes,c)
        result['coarse'].append(one)
    result['input_hashes'] = FILES
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for key,row in result['held'].items():
        print('held',key,'terminal_start',row['readout']['terminal_from_fraction'],
              'collar',row['readout']['sections'][-2],'ring_chi',span(row['ring']['chi_at_nodes']))
    for row in result['coarse']:
        fold = row['selected_prism']
        print('coarse',row['iteration'],'det',row['chart_min'],'selected_chi',fold['point_chi'],
              'prism',fold['whole_prism_chi_enclosure'],'disk_terminal',row['readout']['terminal_from_fraction']['disk'])
    print('selected_fold',json.dumps(result['coarse'][-1]['selected_prism']))
    print('output_sha256',digest(args.output),flush=True)


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('run','gate','original','held','wider','output'):
        p.add_argument(key,type=Path)
    args = p.parse_args()
    pin = args.run/'pin'
    for path,h in read(args.run/'sources.json').items():
        input_file(pin/path,h)
    sys.path.insert(0,str(pin/'software/stage1_solver/tools'))
    run(args)
