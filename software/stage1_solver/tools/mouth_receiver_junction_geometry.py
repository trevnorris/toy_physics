"""Build common-face 4D material coupons from the measured undriven transit endpoints."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


PAIRINGS = {"training": (40, 0), "heldout": (42, 4)}
SHELL = 1.30
MOBILITY = 0.00022797992459479204


def read(path):
    return json.loads(path.read_bytes())


def digest(raw):
    return sha256(raw).hexdigest()


def source(base, history):
    name = f"id{history}_zero"
    raw = (base / "runs" / (name + ".json")).read_bytes()
    checksum = digest(raw)
    assert checksum == read(base / "ZERO_RUNS.json")["results"][name]
    run = json.loads(raw)
    assert (run["history"], run["strain"], run["pressure_external"], run["temperature"]) == (history, 0, 0, 1.088)
    return run, checksum


def face_frame(box):
    """Minimal proper rotation on column vectors that brings the face normal to +w."""
    w = np.eye(4)[3]
    n = np.linalg.solve(box, w)
    n /= np.linalg.norm(n)
    if n @ w < 0:
        raise ValueError("Native cell points away from the specified interface normal")
    skew = np.outer(w, n) - np.outer(n, w)
    rotation = np.eye(4) + skew + (skew @ skew) / (1. + n @ w)
    oriented = box @ rotation.T
    np.testing.assert_allclose(oriented[:3, 3], 0, atol=2e-12)
    np.testing.assert_allclose(rotation @ rotation.T, np.eye(4), atol=2e-12)
    np.testing.assert_allclose(np.linalg.det(rotation), 1.)
    oriented[:3, 3] = 0.
    return n, oriented


def shell_vectors(state):
    """Yield site index, periodic neighbor IDs and displacement vectors inside the shell."""
    box = np.asarray(state["box"])
    pos = np.asarray(state["pos"])
    inv = np.linalg.inv(box)
    x = ((pos @ inv) % 1.) @ box
    # Central cell to all images within the cutoff. The face spacing bounds each fractional offset.
    extents = np.ceil(1. + SHELL * np.linalg.norm(inv, axis=0)).astype(int) - 1
    shifts = np.array(list(product(*(range(-k, k+1) for k in extents)))) @ box
    cloud = (x[None, :, :] + shifts[:, None, :]).reshape(-1, 4)
    neighbors = cKDTree(cloud).query_ball_point(x, SHELL)
    for i, indices in enumerate(neighbors):
        bonds = cloud[indices] - x[i]
        radii = np.linalg.norm(bonds, axis=1)
        mask = (radii > 1e-10) & (radii < SHELL)
        yield i, np.array(indices)[mask] % len(x), bonds[mask]


def local_order(state):
    """Return each site's mean pairwise U6(cos(angle))/7 over its distinct shell bonds."""
    out = np.full(len(state["pos"]), np.nan)
    for i, _, bonds in shell_vectors(state):
        unit = bonds / np.linalg.norm(bonds, axis=1)[:,None]
        if len(unit) >= 2:
            cos2 = np.clip(unit @ unit.T, -1., 1.) ** 2
            kernel = ((64.*cos2 - 80.)*cos2 + 24.)*cos2 - 1.
            np.fill_diagonal(kernel, 0.)
            out[i] = kernel.sum() / (7. * len(unit) * (len(unit)-1))
    return out


def reference(base, output):
    result = {}
    for i in (40, 0, 1, 2, 42, 4):
        run, checksum = source(base, i)
        rows = []
        for block in run["blocks"]:
            score = local_order(block["state"])
            if not np.isfinite(score).all():
                raise ValueError(f"Unclassifiable shell in homogeneous preparation {i}")
            rows.append({"mean": float(score.mean()), "site_quantiles": np.quantile(score, [.1, .5, .9]).tolist()})
        box = np.asarray(run["final"]["box"])
        normal, oriented = face_frame(box)
        tensor = np.average([b["mean"]["stress_tensile"] for b in run["blocks"][-4:]], axis=0,
                            weights=[b["steps"] for b in run["blocks"][-4:]])
        result[i] = dict(sha256=checksum, face_normal=normal.tolist(), face_height=float(oriented[3, 3]),
                         pressure_on_native_face=float(-normal @ tensor @ normal),
                         tangent_traction_norm=float(np.linalg.norm(tensor @ normal-(normal @ tensor @ normal)*normal)),
                         face_traction_vector=(tensor @ normal).tolist(), mark6=rows)
    def means(i):
        return [r["mean"] for r in result[i]["mark6"]]
    crystal, receiver = np.mean(means(40)), np.mean(means(0))
    threshold = float((crystal+receiver)/2.)
    direction = np.sign(crystal-receiver)
    passed = all(min(direction*(np.array(means(h))-threshold)) > 0 for h in (40,42))
    passed &= all(max(direction*(np.array(means(h))-threshold)) < 0 for h in (0,1,2,4))
    summary = {"training_midpoint": threshold, "ordered_increasing": bool(direction>0),
               "separates_all_homogeneous_snapshot_means": bool(passed), "case": result}
    output.write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
    for h in result:
        print(f"history {h}, 16 snapshot shell means {min(means(h)):.6f} to {max(means(h)):.6f}, "
              f"Pnn {result[h]['pressure_on_native_face']:.6f}, tangent traction norm {result[h]['tangent_traction_norm']:.6f}")
    print(f"training midpoint {threshold:.7f}, all snapshots separated {passed}")


def write_input(path, box, pos, label, info):
    payload = json.dumps(info, sort_keys=True, allow_nan=False).encode()
    family = digest(payload)
    assert len(label) == len(pos) and np.all((label == 0) | (label == 1))
    with path.open("x") as f:
        f.write(f"{len(pos)} 1.088 {MOBILITY:.17g} {family}\n")
        np.savetxt(f, np.vstack((box, np.linalg.inv(box), pos)), fmt="%.17g")
        np.savetxt(f, label, fmt="%d")
    path.with_suffix(path.suffix+".json").write_bytes(payload + b"\n")


def prepare(base, output, case, material):
    crystal, receiver = PAIRINGS[case]
    h = crystal if material == "ordered" else receiver
    ref, cref = source(base, crystal)
    run, parent = source(base, h)
    native = np.asarray(run["final"]["box"])
    normal, oriented = face_frame(native)
    area = face_frame(np.asarray(ref["final"]["box"]))[1][:3,:3]
    common = oriented.copy()
    common[:3,:3] = area
    common[3,3] = np.linalg.det(native) / np.linalg.det(area)
    np.testing.assert_allclose(np.linalg.det(common), np.linalg.det(native))
    # Logarithmic fit strain after a proper rotation of the source, at fixed constituent density.
    singular = np.linalg.svd(np.linalg.solve(oriented, common), compute_uv=False)
    pos = ((np.asarray(run["final"]["pos"]) @ np.linalg.inv(native)) % 1.) @ common
    info = dict(case=case, material=material, source=parent, history=h, reference=cref,
                tangent_vectors=area.tolist(), original_face_normal=normal.tolist(),
                input_face_height=float(common[3,3]), native_face_height=float(oriented[3,3]),
                fit_principal_stretches=singular.tolist(), native_to_common_volume_ratio=float(np.linalg.det(common)/np.linalg.det(native)))
    write_input(output, common, pos, np.full(len(pos), int(material == "receiver")), info)
    print(json.dumps(info, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="task", required=True)
    for task in ("parents", "prepare"):
        item = sub.add_parser(task)
        item.add_argument("base", type=Path)
        item.add_argument("output", type=Path)
        if task == "prepare":
            item.add_argument("--case", required=True, choices=PAIRINGS)
            item.add_argument("--material", required=True, choices=("ordered", "receiver"))
    args = parser.parse_args()
    if args.task == "parents":
        reference(args.base, args.output)
    else:
        prepare(args.base, args.output, args.case, args.material)
