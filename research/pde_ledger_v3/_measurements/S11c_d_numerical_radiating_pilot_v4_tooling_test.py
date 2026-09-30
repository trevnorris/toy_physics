"""Exact restored-array comparison regression; synthetic stdlib objects only."""
import pickle
import unittest
from S11c_d_numerical_radiating_saved_operations import exact_structure
from S11c_d_numerical_radiating_pilot_v3_tooling_test import Tests as PriorRegressionTests


class ArrayStandIn:
    """Multi-element array's ambiguous scalar equality, without NumPy import."""
    def __init__(self, data=b'12345678', dtype='f8', shape=(1,)):
        self.data, self.dtype, self.shape = data, dtype, shape
    def tobytes(self):
        return self.data
    def __eq__(self, other):
        raise ValueError('The truth value of an array is ambiguous')


class Tests(unittest.TestCase):
    def test_shared_identity_masked_error_until_restoration(self):
        array = ArrayStandIn()
        old = {'source': [{'values': array, 'index': 7}]}
        shared = {'source': [{'values': array, 'index': 7}]}
        self.assertTrue(old == shared)
        restored = pickle.loads(pickle.dumps(shared))
        with self.assertRaisesRegex(ValueError, 'truth value'):
            bool(old == restored)
        self.assertTrue(exact_structure(old, restored))

    def test_changed_array_bytes_dtype_shape_rejected(self):
        source = {'values': ArrayStandIn()}
        for changed in (ArrayStandIn(b'12345679'), ArrayStandIn(dtype='i8'),
                        ArrayStandIn(shape=(1, 1))):
            self.assertFalse(exact_structure(source, {'values': changed}))
            self.assertFalse(exact_structure({'values': changed}, source))

    def test_changed_scalar_container_keys_and_length_rejected(self):
        source = {'items': [ArrayStandIn(), 3]}
        for changed in ({'items': [ArrayStandIn(), 4]},
                        {'different': [ArrayStandIn(), 3]},
                        {'items': (ArrayStandIn(), 3)},
                        {'items': [ArrayStandIn()]}):
            self.assertFalse(exact_structure(source, changed))
        self.assertTrue(exact_structure(source, pickle.loads(pickle.dumps(source))))


if __name__ == '__main__':
    unittest.main()
