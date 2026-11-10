"""Record the actual portable resources shipped in each wheel."""
import hashlib
import json
from pathlib import Path
from setuptools import setup
from setuptools.command.build_py import build_py

class BuildResources(build_py):
    def run(self):
        super().run()
        root=Path(self.build_lib)/'capzyme_resources'
        manifest={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in sorted(root.rglob('*')) if p.is_file() and p.name!='source_snapshot.json'}
        (root/'source_snapshot.json').write_text(json.dumps({'kind':'installed-wheel',
            'distribution':'capzyme-dad','version':'0.1.0.dev3','files_sha256':manifest},indent=2)+'\n')

setup(cmdclass={'build_py':BuildResources})
