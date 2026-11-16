"""Compare archived fine-grid pair readings and change only the target's outlet law."""

import re
import sys
from decimal import Decimal, getcontext
from pathlib import Path


def numpy_rows(path):
    rows = {}
    mobility = {}
    for section in path.read_text().split("separation=")[1:]:
        header, *lines = section.splitlines()
        case, degree = header.split(" degree=")
        if degree != "14":
            continue
        distance = 0 if case == "None" else int(Decimal(case))
        for line in lines:
            if line.startswith("k="):
                fields = dict(item.split("=") for item in line.split())
                rows[distance, int(fields["k"])] = [
                    Decimal(fields[key]) for key in ("H_target", "a_minus", "a_plus", "a_free")
                ]
            if "mobility_target_row=" in line:
                mobility[distance] = [Decimal(x) for x in line.split("mobility_target_row=")[1].split()]
    return rows, mobility


def wl_vector(text):
    values = re.sub(r"`[0-9.]+", "", text).replace("*^", "e")
    return [Decimal(x.strip()) for x in values.split(",")]


def wolfram_rows(path):
    rows = {}
    mobility = {}
    fine = False
    for line in path.read_text().splitlines():
        if line.startswith('"centers="'):
            centers, degree = re.fullmatch(r'"centers=""\{([^}]+)\}"" degree=(\d+)"', line).groups()
            distance = int(centers.split(",")[-1])
            fine = degree == "9"
        elif fine and line.startswith('"k="'):
            k, data = re.fullmatch(
                r'"k="(\d)" \{H_target,a_minus,a_plus,a_free\}=""\{([^}]+)\}"" controls=""\{True, True, True, True\}"',
                line,
            ).groups()
            rows[distance, int(k)] = wl_vector(data)
        elif fine and line.startswith('"mobility_target_row="'):
            data = re.fullmatch(r'"mobility_target_row=""\{([^}]+)\}"', line).group(1)
            mobility[distance] = wl_vector(data)
    return rows, mobility


def main(numpy_stdout, wolfram_stdout):
    getcontext().prec = 45
    nrows, nmob = numpy_rows(numpy_stdout)
    wrows, wmob = wolfram_rows(wolfram_stdout)
    distances = {0, 10, 20}
    expected = {(d, k) for d in distances for k in (0, 1)}
    assert set(nrows) == set(wrows) == expected and set(nmob) == set(wmob) == distances
    row_errors = [abs(x-y) for key in expected for x, y in zip(nrows[key], wrows[key], strict=True)]
    mob_errors = []
    for d in distances:
        assert len(wmob[d]) == (2 if d else 1)
        mob_errors.extend(abs(x-y) for x, y in zip(nmob[d], wmob[d], strict=True))
    print(f"{len(row_errors)} force/acceleration readings, max absolute difference {max(row_errors):.4E}")
    print(f"{len(mob_errors)} target mobility entries, max absolute difference {max(mob_errors):.4E}")
    if max(row_errors + mob_errors) > Decimal("1e-12"):
        raise ArithmeticError("engine outputs disagree")
    for d in (10, 20):
        # k affects only that outlet's load. Shapes, prescribed intakes and mobility
        # are common. Positive lab x runs from the left target to its right companion.
        def response(rows, mob):
            B0, B1 = [-rows[d, k][0] for k in (0, 1)]
            mu_tt, mu_ts = mob[d]
            return (B0-B1, (mu_tt-mu_ts)*B0, mu_tt*B1-mu_ts*B0)
        first, second = response(nrows, nmob), response(wrows, wmob)
        if max(abs(x-y) for x, y in zip(first, second)) > Decimal("1e-12"):
            raise ArithmeticError("fixed-source comparison disagrees")
        export, a00, a10 = second
        print(f"d={d}, target export change={export:.12g}, source k=0: a(target k=0)={a00:.12g}, a(target k=1)={a10:.12g}")
        print(f"  a/U_monopole for target k=0,1: {a00*d**2:.12g}, {a10*d**2:.12g}")


if __name__ == "__main__":
    main(*(Path(arg) for arg in sys.argv[1:]))
