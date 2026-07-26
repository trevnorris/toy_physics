class DimBad:
    l = 0


def bad_scope(DimBad):
    bad_value = DimBad()
    return bad_value


class DimGood:
    l = 0


good_value = DimGood()
