class DimBad:
    l = 0


def helper():
    return None


if True:
    DimBad = helper

bad_value = DimBad()


class DimGood:
    l = 0


good_value = DimGood()
