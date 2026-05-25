import json
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = REPO_ROOT / "e156-submission" / "config.json"


class PortabilityContracts(unittest.TestCase):
    def test_submission_config_uses_repo_relative_root(self) -> None:
        payload = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

        self.assertEqual(payload["path"], "..")
        resolved_root = (CONFIG_PATH.parent / payload["path"]).resolve()
        self.assertEqual(resolved_root, REPO_ROOT.resolve())

    def test_r_parity_gate_uses_environment_or_path_discovery(self) -> None:
        text = (REPO_ROOT / "tests" / "test_r_parity.py").read_text(encoding="utf-8")

        self.assertIn('os.environ.get("RSCRIPT")', text)
        self.assertIn('shutil.which("Rscript")', text)
        self.assertNotIn(r'C:\Program Files\R\R-4.5.2\bin\Rscript.exe', text)


if __name__ == "__main__":
    unittest.main()
