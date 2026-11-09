#!/usr/bin/env python3
"""Choose an independent CapZyme operation; Snakemake is not required."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / 'workflow/scripts'))
from tool_settings import main

if __name__ == '__main__':
    main()
