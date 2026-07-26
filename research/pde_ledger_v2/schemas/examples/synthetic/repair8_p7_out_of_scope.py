def holder():
    class DimBad:
        l = 0

    return DimBad


bad_value = DimBad()


class DimGood:
    l = 0


good_value = DimGood()
