class DimBad:
    l = 0


class mod:
    Dim = DimBad


bad_module_value = mod.Dim()


class Holder:
    Dim = DimBad


obj = Holder()
bad_object_value = obj.Dim()


class DimGood:
    l = 0


good_value = DimGood()
