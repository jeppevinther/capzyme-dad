"""Installed dispatcher; standalone tools never import Snakemake."""
from pathlib import Path
import runpy
import sys


def main():
    checkout = Path(__file__).resolve().parent.parent
    if (checkout / 'toolbox.py').is_file():
        root = checkout
    else:
        from importlib.util import find_spec
        root = Path(next(iter(find_spec('capzyme_resources').submodule_search_locations)))
    if len(sys.argv) > 1 and sys.argv[1] == 'workflow':
        sys.argv = [str(root / 'capzyme.py'), *sys.argv[2:]]
        script = root / 'capzyme.py'
    elif len(sys.argv) == 1 or sys.argv[1] in ('-h', '--help'):
        print('CapZyme-DAD: Data Analysis and Display for CapZyme-seq\n\n'
              'Commands: count, bigwig, trim, rna-align, genome-align, collapse, index-record\n'
              'Optional coordinator: workflow {doctor,example,check,preview,run}\n'
              'Use capzyme COMMAND --help for inputs, options and outputs.')
        return
    else:
        script = root / 'toolbox.py'
    runpy.run_path(str(script), run_name='__main__')
