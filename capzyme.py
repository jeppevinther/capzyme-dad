#!/usr/bin/env python3
"""Beginner launcher for existing CapZyme-seq Snakemake workflows.

This adapter checks configuration and invokes Snakemake; scientific operations
remain in the existing standalone tools/shared functions. See docs/getting_started.md.
"""
import argparse
from contextlib import contextmanager
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'workflow/scripts'))
ROUTES = {'rna': 'Snakefile', 'sequential': 'Snakefile',
          'genome': 'genome.smk', 'bam': 'existing_bam.smk'}


@contextmanager
def working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def external(path, label):
    path = Path(path).expanduser().resolve()
    if path == ROOT or ROOT in path.parents:
        raise ValueError(f'{label} must be outside the source checkout')
    return path


def doctor(args):
    ok = sys.version_info >= (3, 12)
    print(f'{"OK" if ok else "MISSING"} Python {sys.version.split()[0]}: {sys.executable}')
    yaml_ok = importlib.util.find_spec('yaml') is not None
    print(f'{"OK" if yaml_ok else "MISSING"} PyYAML (configuration reader)')
    ok &= yaml_ok
    tools = ['snakemake', 'samtools', 'bedtools']
    if args.route in ('rna', 'sequential'):
        tools += ['cutadapt', 'bowtie2', 'bowtie2-build', 'bowtie2-inspect']
    if args.route in ('genome', 'sequential'):
        inspector = args.aligner+'-inspect'
        if args.aligner == 'hisat2' and shutil.which('hisat2-inspect-s'):
            inspector = 'hisat2-inspect-s'  # Same vendor-wrapper workaround as the existing tool.
        tools += [args.aligner, inspector]
    for tool in dict.fromkeys(tools):
        executable = shutil.which(tool)
        if not executable:
            print(f'MISSING {tool}: activate the environment described in docs/getting_started.md')
            ok = False
            continue
        result = subprocess.run([executable, '--version'], text=True, capture_output=True,
                                errors='replace', timeout=30)
        lines = (result.stdout or result.stderr).strip().splitlines()
        valid = result.returncode == 0
        if tool == 'snakemake':
            version = re.search(r'(\d+)\.(\d+)\.(\d+)', result.stdout)
            valid = valid and bool(version) and tuple(map(int, version.groups())) >= (9, 9, 0)
        print(f'{"OK" if valid else "FAILED"} {tool}: {lines[0] if lines else "no version output"} [{executable}]')
        ok &= valid
    print('Doctor checks executable availability/version calls, not reference compatibility or a fresh installation.')
    return 0 if ok else 1


def checked_config(args):
    try:
        import yaml
    except ImportError as exc:
        raise ValueError('PyYAML is missing; activate the documented environment and run doctor') from exc
    config = Path(args.config).expanduser().resolve()
    work = external(args.workdir or config.parent, 'Working directory')
    if not work.is_dir():
        raise ValueError(f'Working directory does not exist: {work}')
    try:
        data = yaml.safe_load(config.read_text())
    except yaml.YAMLError as exc:
        raise ValueError(f'Invalid YAML configuration: {exc}') from exc
    if not isinstance(data, dict):
        raise ValueError('Configuration must be a YAML/JSON mapping')
    with working_directory(work):
        if args.route in ('rna', 'sequential'):
            from configuration import load_settings
            settings = load_settings(data, ROOT)
            if bool(settings['genome']) != (args.route == 'sequential'):
                raise ValueError('Use --route sequential with a genome block; use --route rna without one')
        elif args.route == 'genome':
            from genome_configuration import genome_only_settings
            settings = genome_only_settings(data, ROOT)
        else:
            from bam_configuration import load_bam_settings
            settings = load_bam_settings(data, ROOT)
    with working_directory(work):
        from workflow_bigwig import plans
        plans(data.get('bigwigs'), settings, ROOT, bam=args.route=='bam')
    print(f'Configuration valid: {args.route}; samples in order: {", ".join(settings["samples"])}')
    print(f'Working directory: {work}\nOutputs: {settings["output"]}')
    return config, work, settings


def execute(args):
    config, work, settings = checked_config(args)
    executable = shutil.which('snakemake')
    if not executable:
        raise ValueError('Snakemake is missing; activate the documented environment and run doctor')
    command = [executable, '--snakefile', str(ROOT/'workflow'/ROUTES[args.route]),
               '--configfile', str(config), '--directory', str(work), '--executor', 'local',
               '--cores', str(args.cores), '--resources', f'mem_mb={args.memory_mb}', '--printshellcmds']
    if args.operation == 'preview':
        command.append('--dry-run')
    elif args.rerun_incomplete:
        command.append('--rerun-incomplete')
    print('Command: '+shlex.join(command), flush=True)
    return subprocess.run(command, cwd=work).returncode


def example(args):
    destination = external(args.directory, 'Example directory')
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise ValueError('Example directory must be new or empty; existing files are never overwritten')
    spec = importlib.util.spec_from_file_location('capzyme_synthetic_example', ROOT/'examples/synthetic/generate.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    config = module.generate(destination)
    if args.with_bigwigs:
        data=json.loads(config.read_text())
        data['bigwigs']={'rna':{'cohort':['B','A'],'strands':'plus','methods':{'raw':{},'rpm':{}}}}
        config.write_text(json.dumps(data,indent=2)+'\n')
    print(f'Synthetic example created: {destination}\nConfiguration: {config}')
    print('Next: capzyme workflow check --config '+shlex.quote(str(config)))
    return 0


def positive(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError('Use a positive integer')
    return value


def main():
    parser = argparse.ArgumentParser(prog='capzyme workflow', description=__doc__)
    sub = parser.add_subparsers(dest='operation', required=True)
    p = sub.add_parser('doctor', help='Check required executables; does not install software')
    p.add_argument('--route', choices=ROUTES, default='rna')
    p.add_argument('--aligner', choices=['hisat2','bowtie2'], default='hisat2', help='Genome aligner for doctor only')
    p = sub.add_parser('example', help='Create a small artificial RNA-only example outside the checkout')
    p.add_argument('--directory', required=True)
    p.add_argument('--with-bigwigs', action='store_true', help='Enable raw/RPM tracks for the nonempty B,A samples')
    for name, help_text in [('check','Validate configuration without executing jobs'),
                            ('preview','Show the planned jobs; no analysis runs'),
                            ('run','Execute or resume the configured workflow')]:
        p = sub.add_parser(name, help=help_text)
        p.add_argument('--config', required=True, help='YAML/JSON configuration path')
        p.add_argument('--route', choices=ROUTES, default='rna')
        p.add_argument('--workdir', help='External working directory; default: directory containing config')
        if name != 'check':
            p.add_argument('--cores', type=positive, default=4, help='Total CPU scheduling budget (default 4)')
            p.add_argument('--memory-mb', type=positive, default=16384, help='Total memory scheduling budget, not an enforced cap')
        if name == 'run':
            p.add_argument('--rerun-incomplete', action='store_true', help='Retry incomplete jobs after diagnosing interruption/failure')
    args = parser.parse_args()
    try:
        if args.operation == 'doctor': return doctor(args)
        if args.operation == 'example': return example(args)
        if args.operation == 'check': checked_config(args); return 0
        return execute(args)
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        parser.exit(1, f'capzyme: {exc}\n')
    except KeyboardInterrupt:
        parser.exit(130, '\nInterrupted. Keep inputs/configuration and outputs; see docs/results_and_troubleshooting.md to resume.\n')


if __name__ == '__main__':
    raise SystemExit(main())
