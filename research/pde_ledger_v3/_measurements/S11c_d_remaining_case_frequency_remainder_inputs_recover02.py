#!/usr/bin/env python3
"""Resume unfinished metadata, recognizing actual annotated namespace source."""
import ast
from pathlib import Path
import S11c_d_remaining_case_frequency_remainder_inputs_recover as prior

FIRST_SHA = '63fbcfed6328fba72aef58108352aca0a0d29087f06c41183b13e234a8f1a987'


def printer_metadata(reader, printer_file, printer):
    namespace = Path(prior.old.inspect.getsourcefile(prior.old.sp.lambdify))
    module = ast.parse(namespace.read_text())
    def names(node):
        return node.targets if isinstance(node, ast.Assign) else [node.target] if isinstance(node, ast.AnnAssign) else []
    assignment = next(n for n in module.body if any(isinstance(t, ast.Name) and t.id == 'NUMPY_TRANSLATIONS' for t in names(n)))
    translation = ast.literal_eval(assignment.value)['Heaviside']
    prior.old.require(translation == 'heaviside', 'actual native numpy namespace translation')
    method = printer._print_Heaviside
    return {'file': reader.retain(printer_file), 'dedicatedMethodPresent': callable(method),
            'body': prior.old.inspect.getsource(method) if callable(method) else None,
            'namespaceFile': reader.retain(namespace), 'namespaceAssignment': ast.unparse(assignment),
            'actualTranslation': translation, 'lambdifyInvoked': False}


def main():
    helper = Path(prior.__file__).resolve()
    prior.old.require(prior.old.saved.digest(helper) == FIRST_SHA, 'immutable first metadata continuation')
    module = ast.parse(helper.read_text()); fn = next(n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    changed = 0
    for n in ast.walk(fn):
        if isinstance(n, ast.Constant) and n.value == 'S11c_d_remaining_case_frequency_remainder_inputs_recover':
            n.value = 'S11c_d_remaining_case_frequency_remainder_inputs_recover02'; changed += 1
    prior.old.require(changed == 1, 'only new continuation metadata name')
    env = dict(vars(prior), __file__=str(Path(__file__).resolve()), printer_metadata=printer_metadata,
               FAILED=prior.old.F/'remainder-inputs-recovery-02/failed-file-inventory.json', FAILED_SHA='28b398dfda1d84c3b0ef799c929311f87ae38ec38e93aa3ed9e9b278173acfa3')
    exec(compile(ast.Module(body=[fn], type_ignores=[]), '<annotated-source-metadata-continuation>', 'exec'), env); env['main']()


if __name__ == '__main__': main()
