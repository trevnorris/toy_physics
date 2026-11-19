"""Locate compensating winding across a feeding collar with a local energy sink.

Use gauge-invariant loop angles on all sides of triangular prisms. The
incoming material is uniform on the outer cut, and six encounters prepare
it along the flow toward the inner cut. Tilt the mesh independently of the
uniform incoming direction to avoid placing the critical zero on a radial
mesh edge. All spatial normals point away from the center at both cuts.
"""
import numpy as np
from scipy.spatial.transform import Rotation

from mouth_material_links import section
from mouth_collar_order_preparation import calibration, predicted


def states(vectors):
    n = vectors/np.linalg.norm(vectors,axis=1)[:,None]
    x,y,z = n.T
    north = np.column_stack((1+z,x+1j*y))
    south = np.column_stack((x-1j*y,1-z))
    rays = np.where((z>0)[:,None],north,south)
    return rays/np.linalg.norm(rays,axis=1)[:,None]


def face_angle(rays):
    v = np.sum(rays.conj()*np.roll(rays,-1,axis=1),axis=2)
    assert np.min(abs(v))>1e-6
    return np.angle(np.prod(v/abs(v),axis=1))


def scan(level,eta,rotation,survival):
    normals,triangles = section(level)
    normals = Rotation.from_rotvec(rotation).apply(normals)
    previous = None
    total = 0
    locations = []
    for layer in range(7):
        rays = states(predicted(normals,eta,survival,layer))
        # Arbitrary local site phases changed on every shell.
        j = np.arange(len(rays))
        rays *= np.exp(1j*(.791*j+.31*layer*np.sin(j)))[:,None]
        faces = rays[triangles]
        current = face_angle(faces)
        if previous is None:
            outer = sum(current)/(2*np.pi)
        else:
            sides = np.zeros_like(current)
            for i in range(3):
                j=(i+1)%3
                sides += face_angle(np.stack([faces[:,i],faces[:,j],previous[:,j],previous[:,i]],axis=1))
            prism = (last-current+sides)/(2*np.pi)
            integers = np.rint(prism).astype(int)
            assert np.max(abs(prism-integers))<1e-9
            for idx in np.flatnonzero(integers):
                directions = normals[triangles[idx]]
                dot = np.clip(-eta*directions[:,2],-1.,1.)
                locations.append((layer,integers[idx],np.max(np.arccos(dot))))
            total += sum(integers)
        last,previous = current,faces
    inner = sum(last)/(2*np.pi)
    print(f"prisms level={level} eta={eta:+d} rot_x={rotation[0]:g} exterior_degree={outer:.8g} "
          f"interior_cut_degree={inner:.8g} net_collar_index={total} occupied_prisms={len(locations)} "
          f"event_intervals={[row[0] for row in locations]} "
          f"max_vertex_distance_from_opposed_pole_rad={max(row[2] for row in locations):.8g}")
    assert abs(outer)<1e-9 and abs(inner-eta)<1e-9 and total==-eta
    assert len(locations)==1 and locations[0][:2]==(4,-eta)


def main():
    survival=calibration()
    for level in (2,3,4):
        for sign in (-1,1):
            for rotation in ([.26,-.19,.17],[-.13,.24,.31]):
                scan(level,sign,rotation,survival)


if __name__=="__main__":
    main()
