"""Local reconnection of directed material bonds on two finite cubic banks.

Sites (x,y,z,w) occupy w=0 and w=s at unit lattice spacing. Each interior
site initially has one outgoing and one incoming bond per spatial axis.
The open x boundaries can feed a material current independently of the arrows.
"""

from fractions import Fraction
from itertools import product


L = 8
RADII = (1, 2, 4)


def intact_bonds(s):
    edges = set()
    for xyz in product(range(-L, L + 1), repeat=3):
        for axis in range(3):
            if xyz[axis] < L:
                neighbor = list(xyz)
                neighbor[axis] += 1
                for w in (0, s):
                    edges.add(((*xyz, w), (*neighbor, w)))
    return edges


def section_flux(edges, radius, bank):
    """Outward counts on one bank's box, resolved by same-bank and seam bonds."""
    def inside(v):
        x, y, z, w = v
        return w == bank and -radius <= x <= radius + 1 and max(abs(y), abs(z)) <= radius

    surface = exchange = 0
    for u, v in edges:
        crossing = int(inside(u)) - int(inside(v))
        if u[3] == v[3]:
            surface += crossing
        else:
            exchange += crossing
    return surface, exchange


def divergence(edges):
    result = {}
    for tail, head in edges:
        result[tail] = result.get(tail, 0) + 1
        result[head] = result.get(head, 0) - 1
    return result


def feed(edges, lanes, s):
    """One unit of imposed material per tick distributed among donor inlet lanes.

    A unit packet starts at each left donor inlet, follows +x bonds (including
    a possible change of bank) and contributes only at the right receiver.
    """
    successor = {u: v for u, v in edges if v[0] == u[0] + 1}
    arrived = 0
    for y, z in lanes:
        v = (-L, y, z, 0)
        while v in successor:
            v = successor[v]
        arrived += v == (L, y, z, s)
    return Fraction(arrived, len(lanes))


def run(s, radius, initial):
    lanes = [(y, z) for y, z in product(range(-radius, radius + 1), repeat=2)
             if y*y + z*z <= radius*radius]
    removed, bridges = set(), set()
    for y, z in lanes:
        for here, other in ((0, s), (s, 0)):
            tail = (0, y, z, here)
            removed.add((tail, (1, y, z, here)))
            bridges.add((tail, (1, y, z, other)))
    cut = initial - removed
    repaired = cut | bridges
    print(f"s={s:+d} R={radius} cut_each={len(initial - cut)//2} "
          f"bridges_from_B_R={sum(u[3] == 0 for u, v in bridges)},"
          f"{sum(u[3] == s for u, v in bridges)} "
          f"received_intake_intact_cut_repaired={feed(initial, lanes, s)},"
          f"{feed(cut, lanes, s)},{feed(repaired, lanes, s)}")

    # A missing R->B partner leaves endpoints on different banks. A separate
    # missing chain in the donor reaches past the enclosing surface.
    missing = repaired - {((0, 0, 0, s), (1, 0, 0, 0))}
    chain_start = -L + 1
    chain = {((x, 0, 0, 0), (x + 1, 0, 0, 0)) for x in range(chain_start, 0)}
    baseline = divergence(initial)
    for name, graph in (("cut", cut), ("repair", repaired), ("missing", missing),
                        ("remote", repaired - chain)):
        degree = divergence(graph)
        delta = {v: degree.get(v, 0) - d for v, d in baseline.items()}
        tips = tuple(sum(delta[(x, y, z, w)] for y, z in lanes)
                     for w in (0, s) for x in (0, 1))
        bank_flux = [section_flux(graph, radius + 1, w) for w in (0, s)]
        # Independently sum vertex divergence over each section.
        interior = [sum(value for (x, y, z, w), value in degree.items()
                        if w == bank and -radius - 1 <= x <= radius + 2
                        and max(abs(y), abs(z)) <= radius + 1) for bank in (0, s)]
        print(f"  {name} tip_delta_B0_B1_R0_R1={tips} "
              f"outward_surface_exchange_B_R={bank_flux} enclosed_div_B_R={interior} "
              f"remote_tip={delta[(chain_start, 0, 0, 0)]} "
              f"changed_sites={sum(v != 0 for v in delta.values())}")
        assert [sum(flux) for flux in bank_flux] == interior


if __name__ == "__main__":
    for orientation in (1, -1):
        intact = intact_bonds(orientation)
        for physical_radius in RADII:
            run(orientation, physical_radius, intact)
