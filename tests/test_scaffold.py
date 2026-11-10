"""Repository contract checks; no sequencing files or workflow execution."""
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]

class ScaffoldTests(unittest.TestCase):
    def test_required_layout(self):
        for name in ("workflow/Snakefile", "workflow/rules", "workflow/scripts",
                     "config", "tests", "docs", "examples", "legacy_scripts",
                     "R", "README.md", "AGENTS.md", ".gitignore", "LICENSE"):
            with self.subTest(path=name):
                self.assertTrue((ROOT / name).exists())

    def test_data_and_private_files_are_ignored(self):
        paths = ["nested/reads.fastq.gz", "nested/reads.fq.gz", "nested/a.bam",
                 "nested/a.bam.bai", "nested/ref.fa.gz", "nested/ref.fasta",
                 "nested/ref.1.bt2", "nested/ref.1.ht2", "nested/a.bigWig",
                 "results/counts.tsv", "config/local/experiment.yaml", ".env"]
        proc = subprocess.run(["git", "check-ignore", "--no-index", "--stdin"],
                              cwd=ROOT, input="\n".join(paths) + "\n",
                              text=True, capture_output=True, check=True)
        self.assertEqual(set(proc.stdout.splitlines()), set(paths))

    def test_portable_runtime_templates(self):
        for folder in ("workflow", "config", "profiles", "R"):
            for path in (ROOT / folder).rglob("*"):
                if path.is_file() and "__pycache__" not in path.parts and path.suffix not in (".pyc", ".pyo"):
                    text = path.read_text().lower()
                    for forbidden in ("/projects/", "/users/", "unicph.domain", "scarball"):
                        with self.subTest(path=path, forbidden=forbidden):
                            self.assertNotIn(forbidden, text)

if __name__ == "__main__":
    unittest.main()
