"""Select serialized b provenance components without restoring giant trees."""
import ast


def origin_components(serialized):
    """Return every case/slot/origin value with its exact serialized operand.

    Only the Tuple/Str association envelope is interpreted. The caller restores
    the small KINETIC operands; all other components remain exact source slices.
    Reject duplicate keys or unexpected envelope syntax rather than dropping it.
    """
    encoded = serialized.encode('utf-8')
    if b'\n' in encoded:
        raise ValueError('single-line srepr required')

    def arguments(node, name):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == name and not node.keywords):
            raise ValueError(('association syntax', name))
        return node.args

    def key(node):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'Str':
            values = arguments(node, 'Str')
            if len(values) != 1 or not isinstance(values[0], ast.Constant) or not isinstance(values[0].value, str):
                raise ValueError('string key syntax')
            return values[0].value
        return tuple(key(value) for value in arguments(node, 'Tuple'))

    def association(node):
        result = {}
        for entry in arguments(node, 'Tuple'):
            pair = arguments(entry, 'Tuple')
            if len(pair) != 2:
                raise ValueError('association pair arity')
            label = key(pair[0])
            if label in result:
                raise ValueError(('duplicate key', label))
            result[label] = pair[1]
        return result

    result = {}
    for case, payload in association(ast.parse(serialized, mode='eval').body).items():
        for slot, value in association(payload).items():
            entries = association(value).items() if slot == 'VALUE' else [(None, value)]
            for origin, operand in entries:
                path = (case, slot) + ((origin,) if origin is not None else ())
                result[path] = encoded[operand.col_offset:operand.end_col_offset].decode('utf-8')
    return result
