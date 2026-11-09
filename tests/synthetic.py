"""Backward-compatible test/CLI import of the shared public synthetic example."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location('capzyme_synthetic_example',
    Path(__file__).resolve().parents[1]/'examples/synthetic/generate.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
ADAPTER, sequence, generate = module.ADAPTER, module.sequence, module.generate


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    print(generate(parser.parse_args().directory))
