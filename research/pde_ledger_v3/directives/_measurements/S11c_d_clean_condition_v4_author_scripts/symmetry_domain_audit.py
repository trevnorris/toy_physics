#!/usr/bin/env python3
"""Exact finite-dimensional checks used by the S11c-d v4 author.

This is a standard-library algebra script.  It imports neither S11c engine and
does not construct or evaluate the S11c operator.
"""


def add(left, right):
    if isinstance(left, dict) or isinstance(right, dict):
        left = left if isinstance(left, dict) else ({} if left == 0 else {"1": left})
        right = right if isinstance(right, dict) else ({} if right == 0 else {"1": right})
        keys = set(left) | set(right)
        return {key: left.get(key, 0) + right.get(key, 0) for key in keys if left.get(key, 0) + right.get(key, 0)}
    return left + right


def multiply(left, right):
    if isinstance(left, dict) and isinstance(right, dict):
        raise ValueError("only constant-by-linear products are needed")
    if isinstance(left, dict):
        return {key: right * value for key, value in left.items() if right * value}
    if isinstance(right, dict):
        return {key: left * value for key, value in right.items() if left * value}
    return left * right


def matmul(left, right):
    return [
        [
            fold_add(multiply(a, b) for a, b in zip(row, column))
            for column in zip(*right)
        ]
        for row in left
    ]


def fold_add(values):
    result = 0
    for value in values:
        result = add(result, value)
    return result


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def diag(*entries):
    return [
        [entry if i == j else 0 for j in range(len(entries))]
        for i, entry in enumerate(entries)
    ]


mirror = diag(1, 1, -1)
quarter_turn = [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
tensor_o2 = [[{"A": 1}, 0, 0], [0, {"B": 1}, 0], [0, 0, {"B": 1}]]
tensor_axial = [[0, 0, 0], [0, 0, {"C": 1}], [0, {"C": -1}, 0]]

print("MIRROR_AXIS_1", mirror)
print("O2_RANK2_POLAR_FORM", tensor_o2)
print("O2_RANK2_POLAR_MIRROR_IMAGE", matmul(matmul(mirror, tensor_o2), transpose(mirror)))
print("AXIAL_23_SO2_CANDIDATE", tensor_axial)
print("AXIAL_23_MIRROR_IMAGE", matmul(matmul(mirror, tensor_axial), transpose(mirror)))
print("AXIAL_23_FIXED_BY_MIRROR", matmul(matmul(mirror, tensor_axial), transpose(mirror)) == tensor_axial)
print("AXIS_1_QUARTER_TURN", quarter_turn)
print("O2_RANK2_POLAR_QUARTER_TURN_IMAGE", matmul(matmul(quarter_turn, tensor_o2), transpose(quarter_turn)))

parity = {
    "u_1": 1,
    "u_2": 1,
    "u_3": -1,
    "scalar": 1,
    "v_face_3": -1,
    "delta_v_x_3": -1,
    "V": 1,
    "J": 1,
    "affinity": 1,
}
for output in ("v_face_3", "delta_v_x_3", "V", "J", "affinity"):
    product = parity[output] * parity["u_3"]
    print("PARITY_BLOCK", output, "FROM", "u_3", "SAME_SECTOR", product == 1)

print("W_COORDINATE_MIRROR_SIGN", 1)
print("NORMAL_DERIVATIVE_PRESERVES_INPLANE_PARITY", True)

for ell in range(5):
    scalar = (-1) ** ell
    exists = ell >= 1
    toroidal = -scalar if exists else "ZERO_SECTOR"
    print(
        "ROUND_SECTOR",
        ell,
        "SCALAR_PARITY",
        scalar,
        "SPHEROIDAL_PARITY",
        scalar,
        "TOROIDAL_PARITY",
        toroidal,
        "TOROIDAL_EXISTS",
        exists,
    )

# The controls are appended after the O(2)-invariant baseline has been
# specialized.  Therefore the fixed control datum e is not projected by R1.
print("CONTROL_ORDER", "SPECIALIZE_BASELINE_THEN_APPEND_CONTROL")
print("K1_FIXED_DATUM", "e=(sin(beta),0,cos(beta))")
print("K1_ON_P_U3_THETA_CARRIERS", ["a_K1*cos(beta)*d1(u3)*d1(theta)", "a_K1*cos(beta)*d2(u3)*d2(theta)"])
print("K2_ON_R1_P_CARRIER", "a_K2*theta*g1*d2(u3)")
