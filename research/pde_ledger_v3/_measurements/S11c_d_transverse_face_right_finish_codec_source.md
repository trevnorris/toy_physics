# Saved matrix storage: source-only codec evidence

The static opcode census found DomainMatrix and IntegerRing only in the two
saved REFERENCE/LEFT face-premise returns, within mutable matrix `_rep` state.
Each shown storage argument is `{3: {0: 1}}`, shape `(5, 1)`, over IntegerRing.
The worker admits these exact two classes in addition to its existing codec.
No general function/eval permission is added. Runtime restoration remains guarded.
This inspection did not import SymPy, call a reducer or restore any payload.

Source: /home/trevnorris/.local/lib/python3.10/site-packages/sympy/polys/matrices/domainmatrix.py
SHA256: 32f8d647343894acd7dab0a5f2646da6221749eeb8d7fc7a38cad3afc8ac70f0

```python
    def __new__(cls, rows, shape, domain, *, fmt=None):
        """
        Creates a :py:class:`~.DomainMatrix`.

        Parameters
        ==========

        rows : Represents elements of DomainMatrix as list of lists
        shape : Represents dimension of DomainMatrix
        domain : Represents :py:class:`~.Domain` of DomainMatrix

        Raises
        ======

        TypeError
            If any of rows, shape and domain are not provided

        """
        if isinstance(rows, (DDM, SDM, DFM)):
            raise TypeError("Use from_rep to initialise from SDM/DDM")
        elif isinstance(rows, list):
            rep = DDM(rows, shape, domain)
        elif isinstance(rows, dict):
            rep = SDM(rows, shape, domain)
        else:
            msg = "Input should be list-of-lists or dict-of-dicts"
            raise TypeError(msg)

        if fmt is not None:
            if fmt == 'sparse':
                rep = rep.to_sdm()
            elif fmt == 'dense':
                rep = rep.to_ddm()
            else:
                raise ValueError("fmt should be 'sparse' or 'dense'")

        # Use python-flint for dense matrices if possible
        if rep.fmt == 'dense' and DFM._supports_domain(domain):
            rep = rep.to_dfm()

        return cls.from_rep(rep)

    def __reduce__(self):
        rep = self.rep
        if rep.fmt == 'dense':
            arg = self.to_list()
        elif rep.fmt == 'sparse':
            arg = dict(rep)
        else:
            raise RuntimeError # pragma: no cover
        args = (arg, rep.shape, rep.domain)
        return (self.__class__, args)

```

Source: /home/trevnorris/.local/lib/python3.10/site-packages/sympy/polys/domains/integerring.py
SHA256: e28cb8f714e2f2157aaa1f0259402805c827d9a26a82ac3287b6d906c8c67f02

```python
class IntegerRing(Ring, CharacteristicZero, SimpleDomain):
    r"""The domain ``ZZ`` representing the integers `\mathbb{Z}`.

    The :py:class:`IntegerRing` class represents the ring of integers as a
    :py:class:`~.Domain` in the domain system. :py:class:`IntegerRing` is a
    super class of :py:class:`PythonIntegerRing` and
    :py:class:`GMPYIntegerRing` one of which will be the implementation for
    :ref:`ZZ` depending on whether or not ``gmpy`` or ``gmpy2`` is installed.

    See also
    ========

    Domain
    """

    rep = 'ZZ'
    alias = 'ZZ'
    dtype = MPZ
    zero = dtype(0)
    one = dtype(1)
    tp = type(one)


    is_IntegerRing = is_ZZ = True
    is_Numerical = True
    is_PID = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def __init__(self):
        """Allow instantiation of this domain. """

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        if isinstance(other, IntegerRing):
            return True
        else:
            return NotImplemented
```

