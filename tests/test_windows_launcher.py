import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
POWERSHELL = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / \
    "WindowsPowerShell" / "v1.0" / "powershell.exe"


@unittest.skipUnless(os.name == "nt" and POWERSHELL.is_file(), "Windows PowerShell required")
class LauncherTests(unittest.TestCase):
    def run_launcher(self, existing=False, fail_create=False, fail_run=""):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fake = root / "conda.cmd"
            log = root / "calls.txt"
            fake.write_text(
                '@echo off\n'
                'if "%1"=="env" (\n'
                '  if "%FAKE_EXISTING%"=="1" (echo {"envs":["C:\\\\conda\\\\envs\\\\basin-project"]}) '
                'else (echo {"envs":[]})\n'
                '  exit /b 0\n'
                ')\n'
                'if "%1"=="create" (\n'
                '  echo create>>"%FAKE_LOG%"\n'
                '  if "%FAKE_FAIL_CREATE%"=="1" goto failed\n'
                '  exit /b 0\n'
                ')\n'
                'if "%1"=="run" (\n'
                '  echo %7>>"%FAKE_LOG%"\n'
                '  if "%7"=="%FAKE_FAIL_RUN%" goto failed\n'
                '  exit /b 0\n'
                ')\n'
                'exit /b 42\n'
                ':failed\n'
                'exit /b 31\n',
                encoding="ascii")
            env = os.environ.copy()
            env.update({"CONDA_EXE": str(fake), "FAKE_LOG": str(log),
                        "FAKE_EXISTING": "1" if existing else "0",
                        "FAKE_FAIL_CREATE": "1" if fail_create else "0",
                        "FAKE_FAIL_RUN": fail_run,
                        "PATH": str(POWERSHELL.parent.parent.parent)})
            result = subprocess.run(
                [str(POWERSHELL), "-NoProfile", "-ExecutionPolicy", "Bypass",
                 "-File", str(ROOT / "run_all.ps1")],
                cwd=ROOT, env=env, text=True, capture_output=True, timeout=30)
            calls = log.read_text(encoding="ascii").splitlines() if log.exists() else []
            return result, calls

    def test_creates_missing_environment_before_running(self):
        result, calls = self.run_launcher()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(calls, ["create", "basinlab.experiments.check_environment",
                                 "unittest", "basinlab.run_all"])

    def test_existing_environment_is_not_modified(self):
        result, calls = self.run_launcher(existing=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("create", calls)
        self.assertEqual(len(calls), 3)

    def test_failed_creation_stops_before_running(self):
        result, calls = self.run_launcher(fail_create=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(calls, ["create"])

    def test_failed_tests_stop_before_experiments(self):
        result, calls = self.run_launcher(existing=True, fail_run="unittest")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(calls, ["basinlab.experiments.check_environment", "unittest"])


if __name__ == "__main__":
    unittest.main()
