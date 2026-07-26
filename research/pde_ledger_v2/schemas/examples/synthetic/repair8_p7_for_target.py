class DimBad:
    l = 0


for DimBad in [lambda: None]:
    pass


bad_value = DimBad()


class DimGood:
    l = 0


good_value = DimGood()
